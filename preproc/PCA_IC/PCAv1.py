 
from mpi4py import MPI
import numpy as np
from netCDF4 import Dataset
from datetime import datetime

from commons.mask import Mask

from matplotlib import pyplot as plt

class mpi:
    status=MPI.Status()
    comm = MPI.COMM_WORLD.Clone()
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    @staticmethod
    def send(message, dest, tag=0, counters=None):
        t0=time.time()
        mpi.comm.send(message, dest=dest, tag=tag)
        if counters is not None:
            counters.add(mpi=time.time()-t0)
    
    @staticmethod
    def recv(buf=None, source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, counters=None, sleep=0.05):
        t0=time.time()
        
        message=mpi.comm.recv(buf=buf, source=source, tag=tag, status=mpi.status)
    
        #request=mpi.comm.irecv(buf=buf, source=source, tag=tag)
        #recieved=False
        #while not recieved:
            #recieved, message=request.test(status=mpi.status)
            #time.sleep(sleep)
            
        if counters is not None:
            counters.add(mpi=time.time()-t0)
        return message
    
    @staticmethod
    def Barrier():
        log('Barrier')
        mpi.comm.Barrier()
        
class hMask:
    def __init__(self, mask, local_r=30):
        self.mask=mask
        self.local_r=local_r
        
class Cell:
    def __init__(self, cell_line, start=0, width=1):
        self.cell_line=cell_line
        self.width=width
        self.x_start=start
        self._vol=None
        
    @property
    def hmask(self):
        return self.cell_line.hmask
        
    @property
    def x_next(self):
        return self.x_start+self.width
    
    @property
    def y_start(self):
        return self.cell_line.y_start
    
    @property
    def y_next(self):
        return self.cell_line.y_next
    
    @property
    def x_start_wb(self):
        temp=self.x_start - self.hmask.local_r
        return 0 if temp <0 else temp
    
    @property
    def x_next_wb(self):
        temp=self.x_next+self.hmask.local_r
        return self.hmask.mask.shape[1] if temp > self.hmask.mask.shape[1] else temp
    
    @property
    def y_start_wb(self):
        return self.cell_line.y_start_wb
    
    @property
    def y_next_wb(self):
        return self.cell_line.y_next_wb
    
    @property
    def vol(self):
        if self._vol is None:
            self._vol=self.hmask.mask[self.y_start_wb:self.y_next_wb, self.x_start_wb:self.x_next_wb].sum()
        return _vol
    
    def inside(self):
        return self.hmask.mask[self.y_start:self.y_next,self.x_start:self.x_next].sum()
    
    def outside(self):
        return self.hmask.mask[self.y_start_wb:self.y_next_wb,self.x_start_wb:self.x_next_wb].sum()
    
class CellLineException(Exception):
    pass
    
class CellLine:
    def __init__(self,lines, start=0, height=1, cells=None, ncells=None):
        self.lines= lines
        self.y_start=start
        self.height=height
        self._maxvol=None
        if cells is None:
            if ncells is None:
                raise Exception("cells and ncells cannot be both None")
            cells,self._maxvol=self._init_cells(ncells)
        self.cells=cells
        
    def _init_cells(self, ncells):
        l=self.hmask.mask.shape[1]//ncells
        delta=self.hmask.mask.shape[1]-l*ncells
        if delta>0:
            c=1
            delta-=1
        else:
            c=0
        cells=[Cell(self, start=0, width=l+c)]
        for i in range(ncells-1):
            if delta>0:
                c=1
                delta-=1
            else:
                c=0
            cells.append(Cell(self, start=cells[-1].x_next, width=l+c))
        return cells, None
    
    def _init_cells(self, ncells):
        array=self.hmask.mask[self.y_start_wb:self.y_next_wb].sum(0)
        array_inside=self.hmask.mask[self.y_start:self.y_next].sum(0)
        j=0
        start=np.empty([ncells], dtype=np.int64)
        stop=np.empty([ncells], dtype=np.int64)
        for i in range(ncells):
            while array_inside[j]==0:
                j+=1
            start[i]=j
            j+=1
            stop[i]=j
        j=self.hmask.mask.shape[1]
        while array_inside[j-1]==0:
            j-=1
        stop[-1]=j
        #print(start)
        #print(stop)
        imax_old=ncells
        saves=[np.zeros([ncells]), [None]*ncells]
        while True:
            expanded_i=start-self.hmask.local_r
            expanded_i[expanded_i<0]=0
            expanded_f=stop+self.hmask.local_r
            expanded_f[expanded_f>self.hmask.mask.shape[1]]=self.hmask.mask.shape[1]
            expanded=np.array([array[i:j].sum() for i,j in zip(expanded_i, expanded_f)])
            imax=np.argmax(expanded)
            if imax_old<ncells:
                if imax<imax_old:
                    saves[0][imax_old]=1/ expanded_max
                    saves[1][imax_old]=(start_old, stop_old)
                elif imax>imax_old:
                    saves[0][:imax]=0
                    saves[1][:imax]=[None]*imax
            expanded_max=expanded[imax]
            start_old=start.copy()
            stop_old=stop.copy()
            if expanded[0]==expanded_max:
                break
            start[imax]+=1
            stop[imax-1]=start[imax]            
            while array_inside[start[imax]] == 0:
                start[imax]+=1
            #print(imax)
            #print(start)
            #print(stop)
            if start[imax]>=stop[imax]:
                break
        saves[0][0]=1/ expanded_max
        saves[1][0]=(start_old, stop_old)
        imax=np.argmax(saves[0])
        start, stop=saves[1][imax]
        cells=[]
        for i in range(ncells):
            cells.append(Cell(self, start=start[i], width=stop[i]-start[i]))
        return cells, int(1/saves[0][imax]+0.5)
    
    @property
    def maxvol(self):
        if self._maxvol is None:
            self._maxvol = max([cell.vol for cell in self.cells])
        return self._maxvol
            
            
    @property
    def ncells(self):
        return len(self.cells)
    
    @property
    def hmask(self):
        return self.lines.hmask
    
    @property
    def y_next(self):
        return self.y_start+self.height
    
    @property
    def y_start_wb(self):
        temp=self.y_start - self.hmask.local_r
        return 0 if temp <0 else temp
    
    @property
    def y_next_wb(self):
        temp=self.y_next+self.hmask.local_r
        return self.hmask.mask.shape[0] if temp > self.hmask.mask.shape[0] else temp
        
        
class Decomposition:
    def __init__(self, hmask, hlines=None, nodes=64, nlines=None):
        self.hmask=hmask #mask deve essere piu' larga che alta
        self._maxvol=None
        if hlines is None:
            hlines, self._maxvol=self._init_hlines(nodes, nlines)
        self.hlines=hlines
        
    def _init_hlines(self, nodes, nlines):
        if nlines is None:
            l=np.sqrt(self.hmask.mask.size/nodes)
            nlines=(self.hmask.mask.shape[0]/l+0.5).astype(np.int0)
            
        ncells=nodes//nlines
        reminder=nodes-nlines*ncells
        if reminder>0:
            deltax=1
            reminder-=1
        else:
            deltax=0
        l=self.hmask.mask.shape[0]//nlines
        reminder_y=self.hmask.mask.shape[0]-nlines*l
        if reminder_y>0:
            c=1
            reminder_y-=1
        else:
            c=0
        hlines=[CellLine(self, start=0, height=l+c, cells=None, ncells=ncells+deltax)]
        for i in range(nlines-1):
            if reminder>0:
                deltax=1
                reminder-=1
            else:
                deltax=0
            if reminder_y>0:
                c=1
                reminder_y-=1
            else:
                c=0
            hlines.append(CellLine(self, start=hlines[-1].y_next, height=l+c, cells=None, ncells=ncells+deltax))
        return hlines, None
    
    def _init_hlines(self, nodes, nlines):
        if nlines is None:
            l=np.sqrt(self.hmask.mask.size/nodes)
            nlines=(self.hmask.mask.shape[0]/l+0.5).astype(np.int0)
            
        ncells=nodes//nlines
        reminder=nodes-nlines*ncells
        cells_per_line=np.ones([nlines], dtype=np.int64)*ncells
        cells_per_line[:reminder]+=1
        hlines, maxvol=self._hlines_by_cell_n(cells_per_line)
        return hlines, maxvol
    
    def _hlines_by_cell_n(self, cells_per_line):
        nlines=len(cells_per_line)
        #array=self.hmask.mask[self.y_start_wb:self.y_next_wb].sum(0)
        array_inside=self.hmask.mask.sum(1)
        j=0
        start=np.empty([nlines], dtype=np.int64)
        stop=np.empty([nlines], dtype=np.int64)
        for i in range(nlines):
            while array_inside[j]==0:
                j+=1
            start[i]=j
            j+=1
            stop[i]=j
        j=self.hmask.mask.shape[0]
        while array_inside[j-1]==0:
            j-=1
        stop[-1]=j
        imax_old=nlines
        imax=nlines-1
        saves=[np.zeros([nlines]), [None]*nlines]
        while True:  
            hlines=[CellLine(self, start=start[i], height=stop[i]-start[i], cells=None, ncells=cells_per_line[i]) for i in range(nlines)]
            expanded=np.array([line.maxvol for line in hlines])
            imax=np.argmax(expanded)
            log(f'init_hlines: {start}')
            if imax_old<nlines:
                if imax<imax_old:
                    saves[0][imax_old]=1/ expanded_max
                    saves[1][imax_old]=hlines_old
                elif imax>imax_old:
                    saves[0][:imax]=0
                    saves[1][:imax]=[None]*imax
            expanded_max=expanded[imax]
            hlines_old=hlines
            if expanded[0]==expanded_max:
                break
            start[imax]+=1
            stop[imax-1]=start[imax]            
            while array_inside[start[imax]] == 0:
                start[imax]+=1
            if start[imax]>=stop[imax]:
                break
        saves[0][0]=1/ expanded_max
        saves[1][0]=hlines_old
        imax=np.argmax(saves[0])
        hlines=saves[1][imax]
        return hlines, int(1/saves[0][imax]+0.5)
    
    @property
    def maxvol(self):
        if self._maxvol is None:
            self._maxvol = max([line.maxvol for line in self.hlines])
        return self._maxvol
    
    @property
    def nodes(self):
        temp=0
        for line in self.hlines:
            temp+=line.ncells
        return temp
    
    def __str__(self):
        s=''
        for i, line in enumerate(self.hlines):
            s+=f'line {i} ({line.ncells} cells, inside_y=[{line.y_start} : {line.y_next}], border_y=[{line.y_start_wb} : {line.y_next_wb}]):\n'
            for j, cell in enumerate(line.cells):
                s+=f'    cell {j}: inside: x=[{cell.x_start} : {cell.x_next}], y=[{cell.y_start} : {cell.y_next}], n={cell.inside()}; border: x=[{cell.x_start_wb} : {cell.x_next_wb}], y=[{cell.y_start_wb} : {cell.y_next_wb}], n={cell.outside()}\n'
        return s
    
    def plot(self):
        fig,ax=plt.subplots()
        ax.imshow(self.hmask.mask.astype(np.int64), origin='lower')
        ylim=ax.get_ylim()
        xlim=ax.get_xlim()
        for line in self.hlines:
            ax.plot(xlim, (line.y_start_wb,)*2,'r')
            ax.plot(xlim, (line.y_next_wb,)*2, 'r')
            for cell in line.cells:
                ax.plot((cell.x_start_wb,)*2, (line.y_start_wb,line.y_next_wb),'r')
                ax.plot((cell.x_next_wb,)*2, (line.y_start_wb,line.y_next_wb),'r')
                ax.plot((cell.x_start,)*2, (line.y_start,line.y_next),'w')
                ax.plot((cell.x_next,)*2, (line.y_start,line.y_next),'w')
            ax.plot(xlim, (line.y_start,)*2,'w')
            ax.plot(xlim, (line.y_next,)*2, 'w')
            
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        plt.show()
            
def log(msg):
    print(str(datetime.now())+': '+str(msg))
    
    
if __name__=='__main__':
    if mpi.rank==0:
        maskfile="/g100_work/OGS_devC/SSPADA/GHOSH/T1/wrkdir/MODEL/meshmask.nc"
        log('start')
        mask=Mask(filename=maskfile, maskvarname="tmask", zlevelsvar="nav_lev", ylevelsmatvar="nav_lat", xlevelsmatvar="nav_lon", dzvarname="e3t", loadtmask=True)
        log('mask loaded')
        log(mask.mask.shape)
        hmask=hMask(mask.mask.sum(0))
        dec=Decomposition(hmask, nodes=33)
        log(dec)
        dec.plot()
        #facciamo una funzione per dividere in n parti unguali una riga. l'obiettivo e' minimizzare il massimo numero di celle per nodo. Poi si prova spostando le linee orizzontali e cambiando il numero di parti per linea. e si ricomincia.
