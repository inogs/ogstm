 
from mpi4py import MPI
from mpi4py.util import dtlib
import numpy as np
from scipy.linalg import eigh
from netCDF4 import Dataset
from datetime import datetime
import time

from commons.mask import Mask

from matplotlib import pyplot as plt

import glob
import sys
import argparse

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
        
    @staticmethod
    def scatter(obj, root=0, counters=None):
        t0=time.time()
        message=mpi.comm.scatter(obj, root)
        if counters is not None:
            counters.add(mpi=time.time()-t0)
        return message
    
    @staticmethod
    def bcast(obj, root=0, counters=None):
        t0=time.time()
        message=mpi.comm.bcast(obj, root)
        if counters is not None:
            counters.add(mpi=time.time()-t0)
        return message
    
    @staticmethod
    def gather(obj, root=0, counters=None):
        t0=time.time()
        message=mpi.comm.gather(obj, root)
        if counters is not None:
            counters.add(mpi=time.time()-t0)
        return message
    
    @staticmethod
    def Gatherv(sendbuf, shapes=None, root=0, counters=None):
        t0=time.time()
        
        if shapes is None:
            shape=np.array(sendbuf.shape, dtype=np.int64)
            dt=dtlib.from_numpy_dtype(shape.dtype)
            if mpi.rank==root:
                shapes=np.empty([mpi.size, len(shape)], dtype=shape.dtype)
                mpi.comm.Gather([shape, len(shape), dt], [shapes, len(shape), dt], root)
            else:
                mpi.comm.Gather([shape, len(shape), dt], None, root)
        
        dt=dtlib.from_numpy_dtype(sendbuf.dtype)
        if mpi.rank==root:
            counts=np.prod(shapes,-1)
            counts_sum=np.cumsum(counts)
            recvbuf=np.empty([counts_sum[-1]], dtype=sendbuf.dtype)
            mpi.comm.Gatherv([sendbuf, sendbuf.size, dt], [recvbuf,counts, dt], root)
            lens=[0]+list(counts_sum)
            message=[recvbuf[start: stop].reshape(shapes[i]) for i, (start, stop) in enumerate(zip(lens[:-1],lens[1:]))]
        else:
            mpi.comm.Gatherv([sendbuf, sendbuf.size, dt], None, root)
            message=None
            
        if counters is not None:
            counters.add(mpi=time.time()-t0)
        return message
    
    @staticmethod
    def Igatherv(sendbuf, shapes=None, root=0, counters=None):
        t0=time.time()
        
        if shapes is None:
            shape=np.array(sendbuf.shape, dtype=np.int64)
            dt=dtlib.from_numpy_dtype(shape.dtype)
            if mpi.rank==root:
                shapes=np.empty([mpi.size, len(shape)], dtype=shape.dtype)
                request=mpi.comm.Igather([shape, len(shape), dt], [shapes, len(shape), dt], root)
                request.Wait()
            else:
                mpi.comm.Igather([shape, len(shape), dt], None, root)
        
        dt=dtlib.from_numpy_dtype(sendbuf.dtype)
        if mpi.rank==root:
            counts=np.prod(shapes,-1)
            counts_sum=np.cumsum(counts)
            recvbuf=np.empty([counts_sum[-1]], dtype=sendbuf.dtype)
            request=mpi.comm.Igatherv([sendbuf, sendbuf.size, dt], [recvbuf,counts, dt], root)
            request.Wait()
            lens=[0]+list(counts_sum)
            message=[recvbuf[start: stop].reshape(shapes[i]) for i, (start, stop) in enumerate(zip(lens[:-1],lens[1:]))]
        else:
            mpi.comm.Igatherv([sendbuf, sendbuf.size, dt], None, root)
            message=None
            
        if counters is not None:
            counters.add(mpi=time.time()-t0)
        return message
    
    @staticmethod
    def allgather(obj, counters=None):
        t0=time.time()
        message=mpi.comm.allgather(obj)
        if counters is not None:
            counters.add(mpi=time.time()-t0)
        return message
    
    @staticmethod
    def Allgatherv(sendbuf, shapes=None, counters=None):
        t0=time.time()
        
        if shapes is None:
            shapes=np.empty([mpi.size, len(sendbuf.shape)], dtype=np.int64)
            dt=dtlib.from_numpy_dtype(shapes.dtype)
            mpi.comm.Allgather([np.array(sendbuf.shape, dtype=shapes.dtype), len(sendbuf.shape), dt], [shapes, len(sendbuf.shape), dt])
        
        counts=np.prod(shapes,-1)
        counts_sum=np.cumsum(counts)
        recvbuf=np.empty([counts_sum[-1]], dtype=sendbuf.dtype)
        dt=dtlib.from_numpy_dtype(sendbuf.dtype)
        assert sendbuf.size == counts[mpi.rank]
        mpi.comm.Allgatherv([sendbuf, sendbuf.size, dt], [recvbuf, counts, dt])
        lens=[0]+list(counts_sum)
        message=[recvbuf[start: stop].reshape(shapes[i]) for i, (start, stop) in enumerate(zip(lens[:-1],lens[1:]))]
        
        if counters is not None:
            counters.add(mpi=time.time()-t0)
        return message

    @staticmethod
    def stats():
        log(f'rank = {mpi.rank}, size = {mpi.size}.')
        
        
        
class hMask:
    def __init__(self, mask, local_r=30):
        self.mask=mask
        self.local_r=local_r
        
class Cell:
    def __init__(self, hmask, x_start, x_next, y_start, y_next):
        self.x_start=x_start
        self.x_next=x_next
        self.y_start=y_start
        self.y_next=y_next
        self._set_properties(hmask)
        
    def _set_properties(self,hmask):
        temp=self.x_start - hmask.local_r
        self.x_start_wb=0 if temp <0 else temp
        
        temp=self.x_next+hmask.local_r
        self.x_next_wb= hmask.mask.shape[1] if temp > hmask.mask.shape[1] else temp
        
        temp=self.y_start - hmask.local_r
        self.y_start_wb= 0 if temp <0 else temp
        
        temp=self.y_next+hmask.local_r
        self.y_next_wb= hmask.mask.shape[0] if temp > hmask.mask.shape[0] else temp
        
        self._vol=None
        self._perim=None
        
    def vol(self, hmask):
        if self._vol is None:
            self._vol=hmask.mask[self.y_start_wb:self.y_next_wb, self.x_start_wb:self.x_next_wb].sum()
        return self._vol
    
    def perim(self, hmask):
        if self._perim is None:
            self._perim=self.vol(hmask)-hmask.mask[self.y_start:self.y_next, self.x_start:self.x_next].sum()
        return self._perim

    def trim(self, hmask, ax=None):
        rect=hmask.mask[self.y_start : self.y_next, self.x_start:self.x_next]
        if ax is None or ax==1:
            ysum=rect.sum(0)
            i=0
            while ysum[i]==0:
                i+=1
            self.x_start+=i
            i=-1
            while ysum[i]==0:
                i-=1
            self.x_next+=i+1
        if ax is None or ax==0:
            xsum=rect.sum(1)
            i=0
            while xsum[i]==0:
                i+=1
            self.y_start+=i
            i=-1
            while xsum[i]==0:
                i-=1
            self.y_next+=i+1
            
        self._set_properties(hmask)
        
        rect=hmask.mask[self.y_start_wb : self.y_next_wb, self.x_start_wb:self.x_next_wb]
        if ax is None or ax==1:
            ysum=rect.sum(0)
            i=0
            while ysum[i]==0:
                i+=1
            self.x_start_wb+=i
            i=-1
            while ysum[i]==0:
                i-=1
            self.x_next_wb+=i+1
        if ax is None or ax==0:
            xsum=rect.sum(1)
            i=0
            while xsum[i]==0:
                i+=1
            self.y_start_wb+=i
            i=-1
            while xsum[i]==0:
                i-=1
            self.y_next_wb+=i+1
        
        return self
    
    def halve(self, hmask):
        epsilon=0.01
        rect=hmask.mask[self.y_start : self.y_next, self.x_start:self.x_next]
        ysum=rect.sum(0)
        low=self.x_start+1
        high=self.x_next-1
        value=0
        perim=0
        while low <= high:
            start = low + (high - low)//2
            stop=start
            while ysum[start-self.x_start]==0:
                start+=1
                if start>=self.x_next:
                    right_cell=None
                    right=0
                    r_perim=0
                    break
            else:
                right_cell=Cell(hmask, start, self.x_next, self.y_start, self.y_next).trim(hmask, 0)
                right=right_cell.vol(hmask)
                r_perim=right_cell.perim(hmask)
            while ysum[stop-self.x_start-1]==0:
                stop-=1
                if stop<=self.x_start:
                    left_cell=None
                    left=0
                    l_perim=0
                    break
            else:
                left_cell=Cell(hmask, self.x_start, stop, self.y_start, self.y_next).trim(hmask, 0)
                left=left_cell.vol(hmask)
                l_perim=left_cell.perim(hmask)
            
            temp=left-right
            maxvol=1/max([left,right])
            tot_perim=1/(l_perim+r_perim)
            if np.abs(value/maxvol-1)<=epsilon:
                if tot_perim>perim:
                    couple=left_cell,right_cell
                    value=maxvol
                    perim=tot_perim
            elif maxvol>value:
                couple=left_cell,right_cell
                value=maxvol
                perim=tot_perim
            if temp==0:
                break
            elif temp < 0:
                low = start + 1
            else:
                high = stop - 1
                
        xsum=rect.sum(1)
        low=self.y_start+1
        high=self.y_next-1
        while low <= high:
            start = low + (high - low)//2
            stop=start
            while xsum[start-self.y_start]==0:
                start+=1
                if start>=self.y_next:
                    right_cell=None
                    right=0
                    r_perim=0
                    break
            else:
                right_cell=Cell(hmask, self.x_start, self.x_next, start , self.y_next).trim(hmask, 1)
                right=right_cell.vol(hmask)     
                r_perim=right_cell.perim(hmask)     
            while xsum[stop-self.y_start-1]==0:
                stop-=1
                if stop<=self.y_start:
                    left_cell=None
                    left=0
                    l_perim=0
                    break
            else:
                left_cell=Cell(hmask, self.x_start, self.x_next, self.y_start, stop ).trim(hmask, 1)
                left=left_cell.vol(hmask)
                l_perim=left_cell.perim(hmask)
                
            temp=left-right
            maxvol=1/max([left,right])
            tot_perim=1/(l_perim+r_perim)
            if np.abs(value/maxvol-1)<=epsilon:
                if tot_perim>perim:
                    couple=left_cell,right_cell
                    value=maxvol
                    perim=tot_perim
            elif maxvol>value:
                couple=left_cell,right_cell
                value=maxvol
                perim=tot_perim
            if temp==0:
                break
            elif temp < 0:
                low = start + 1
            else:
                high = stop - 1
        try:
            return couple
        except:
            raise Exception('cella troppo piccola')
        
    def __str__(self):
        return f"Cell(x = [{self.x_start} : {self.x_next}], y = [{self.y_start} : {self.y_next}],\n  x_wb = [{self.x_start_wb} : {self.x_next_wb}], y_wb = [{self.y_start_wb} : {self.y_next_wb}])"

        
class Decomposition:
    def __init__(self, hmask, cells=None, nodes=64):
        self.hmask=hmask
        self._maxvol=None
        if cells is None:
            cells=[Cell(hmask, 0, hmask.mask.shape[1] , 0, hmask.mask.shape[0]).trim(hmask)]
            while len(cells)<nodes:
                cells.sort(key=lambda cell:cell.vol(hmask))
                cell=cells.pop()
                cells.extend(list(cell.halve(hmask)))
        self.cells=cells
    
    @property
    def maxvol(self):
        if self._maxvol is None:
            self._maxvol = max([cell.vol(self.hmask) for cell in self.cells])
        return self._maxvol
    
    @property
    def nodes(self):
        return len(self.cells)
    
    def __str__(self):
        s=f'Nodes={self.nodes}, maxvol={self.maxvol}. Cells:\n'
        for i, cell in enumerate(self.cells):
            s+=f' cell {i}: inside: x=[{cell.x_start} : {cell.x_next}], y=[{cell.y_start} : {cell.y_next}]; border: x=[{cell.x_start_wb} : {cell.x_next_wb}], y=[{cell.y_start_wb} : {cell.y_next_wb}]; vol={cell.vol(self.hmask)}\n'
        return s
    
    def plot(self):
        fig,ax=plt.subplots()
        ax.imshow(self.hmask.mask.astype(np.int64), origin='lower')
        ylim=ax.get_ylim()
        xlim=ax.get_xlim()
        for cell in self.cells:
            ax.plot((cell.x_start_wb,cell.x_next_wb), (cell.y_start_wb,)*2,'r')
            ax.plot((cell.x_start_wb,cell.x_next_wb), (cell.y_next_wb,)*2, 'r')
            ax.plot((cell.x_start_wb,)*2, (cell.y_start_wb,cell.y_next_wb),'r')
            ax.plot((cell.x_next_wb,)*2, (cell.y_start_wb,cell.y_next_wb),'r')
        for cell in self.cells:
            ax.plot((cell.x_start,)*2, (cell.y_start,cell.y_next),'w')
            ax.plot((cell.x_next,)*2, (cell.y_start,cell.y_next),'w')
            ax.plot((cell.x_start,cell.x_next), (cell.y_start,)*2,'w')
            ax.plot((cell.x_start,cell.x_next), (cell.y_next,)*2, 'w')
            
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        plt.show()
            
def log(msg, node=None):
    if mpi.rank==node or node is None:
        print(f'rank = {mpi.rank} - '+str(datetime.now())+': '+str(msg))
        sys.stdout.flush()
    
def readnc(filename, varname, cell):
    #log('reading '+ filename)
    #with Dataset(filename=filename, mode='r', parallel=True, comm=mpi.comm) as ds:
    with Dataset(filename=filename, mode='r') as ds:
        m=ds.variables[varname]
        if len(m.shape) == 4:
            mask=m[0,:,cell.y_start_wb:cell.y_next_wb, cell.x_start_wb:cell.x_next_wb]
        elif len(m.shape) == 3:
            mask=m[:,cell.y_start_wb:cell.y_next_wb, cell.x_start_wb:cell.x_next_wb]
        else:
            raise Exception('strano file nc')
        
    log(filename + ' read',0)
        
    return np.transpose(np.array(mask), axes=(2,1,0)).copy()

def writenc(filename, varname, cell, bmask, mask_shape, data, fill_value=1.e+20):
    #log('writing '+filename)
    with Dataset(filename=filename, mode='w', parallel=True, comm=mpi.comm, compression='zlib') as nc:
        x = nc.createDimension('x', mask_shape[0])
        y = nc.createDimension('y', mask_shape[1])
        dims=('time',)
        if len(bmask.shape)>2:
            z = nc.createDimension('z', mask_shape[2])
            dims=('time','z','y','x')
        else:
            dims=('time','y','x')
        time = nc.createDimension('time', 1)
        v=nc.createVariable(varname, np.float64, dims, fill_value=fill_value)
        v.set_collective(True)
        data[np.logical_not(bmask)]=fill_value
        if len(bmask.shape)>2:
            data=data.transpose((2,1,0))
            #v[0, :,
              #cell.y_start:cell.y_next,
              #cell.x_start:cell.x_next] = data[:,
                                               #cell.y_start-cell.y_start_wb:cell.y_next-cell.y_start_wb,
                                               #cell.x_start-cell.x_start_wb:cell.x_next-cell.x_start_wb]
        
        else:
            data=data.transpose()
            #log(data.shape)
            #log(cell)
            #log(data[cell.y_start-cell.y_start_wb:cell.y_next-cell.y_start_wb,cell.x_start-cell.x_start_wb:cell.x_next-cell.x_start_wb].shape)
            #log(v)
            #v[0, 
              #cell.y_start:cell.y_next,
              #cell.x_start:cell.x_next] = data[cell.y_start-cell.y_start_wb:cell.y_next-cell.y_start_wb,
                                               #cell.x_start-cell.x_start_wb:cell.x_next-cell.x_start_wb]
        v[0, ...,
          cell.y_start:cell.y_next,
          cell.x_start:cell.x_next] = data[...,
                                           cell.y_start-cell.y_start_wb:cell.y_next-cell.y_start_wb,
                                           cell.x_start-cell.x_start_wb:cell.x_next-cell.x_start_wb]
    log(filename + ' written',0)
    
class IO:
    def __init__(self, cells, mask_shape, fill_value=1.e+20, compression='zlib', complevel=4):
        self.filename=None
        self.varname=None
        self.tracers=None
        self.bmasks=None
        self.dim2flats=None
        self.index=0
        self.cells=cells
        self.mask_shape=mask_shape
        self.fill_value=fill_value
        self.compression=compression
        self.complevel=complevel
        
    def write(self, filename, varname, flatten_tracer, bmasks, dim2flats):   
        if False:
            tracers=mpi.gather(flatten_tracer, self.index)
            
            if mpi.rank==self.index:
                self.filename=filename
                self.varname=varname
                self.tracers=tracers
                self.bmasks=bmasks
                self.dim2flats=dim2flats
                log(filename+" ready for writing")
        else:
            if mpi.rank==self.index:
                lens=[[dim2flat.reshape([-1])[-1]+1] for dim2flat in dim2flats]
                
                #lens_sum=np.array(lens).cumsum()
                #tracers=np.empty(lens_sum[-1])
                #request=mpi.comm.Igatherv([flatten_tracer, flatten_tracer.size, MPI.DOUBLE], [tracers, lens, MPI.DOUBLE], self.index)
                #request.Wait()
                #lens=[0]+list(lens_sum)
                #self.tracers=[tracers[start: stop] for start, stop in zip(lens[:-1],lens[1:])]
                
                self.tracers=mpi.Gatherv(flatten_tracer, shapes=lens, root=self.index)
                self.filename=filename
                self.varname=varname
                self.bmasks=bmasks
                self.dim2flats=dim2flats
                log(filename+" ready for writing")
            else:
                #mpi.comm.Igatherv([flatten_tracer, flatten_tracer.size, MPI.DOUBLE], None, self.index)
                mpi.Gatherv(flatten_tracer, shapes=0, root=self.index)
            
        self.index+=1
        if self.index==mpi.size:
            self.write_all()        
    
    def write_all(self):
        if mpi.rank>=self.index:
            self.index=0
            return
        with Dataset(filename=self.filename, mode='w', compression=self.compression, complevel=self.complevel) as nc:
            x = nc.createDimension('x', self.mask_shape[0])
            y = nc.createDimension('y', self.mask_shape[1])
            dims=('time',)
            if len(self.bmasks[0].shape)>2:
                z = nc.createDimension('z', self.mask_shape[2])
                dims=('time','z','y','x')
            else:
                dims=('time','y','x')
            time = nc.createDimension('time', 1)
            v=nc.createVariable(self.varname, np.float64, dims, fill_value=self.fill_value)
            
            for cell, tracer, bmask, dim2flat in zip(self.cells, self.tracers, self.bmasks, self.dim2flats):
                data=tracer[dim2flat]
                data[np.logical_not(bmask)]=self.fill_value
                data=data.transpose(range(len(bmask.shape)-1,-1,-1))
                v[0, ...,
                    cell.y_start:cell.y_next,
                    cell.x_start:cell.x_next] = data[...,
                                                    cell.y_start-cell.y_start_wb:cell.y_next-cell.y_start_wb,
                                                    cell.x_start-cell.x_start_wb:cell.x_next-cell.x_start_wb]
        log(self.filename + ' written')
        self.index=0
    
def main():
    maskfile="/g100_work/OGS_devC/SSPADA/GHOSH/T1/wrkdir/MODEL/meshmask.nc"
    local_r=30
    edge=local_r
    n_eigh=100
    EnsSize=24
    #prefix='/g100_scratch/userexternal/sspada00/restart_PCA/test/input/RST'
    #prefix="/g100_work/OGS_devC/SSPADA/GHOSH/T1/wrkdir/MODEL/RESTART/"
    template='/g100_scratch/userexternal/sspada00/restart_PCA/test/input/RST{member}.{var}.nc'
    output='/g100_scratch/userexternal/sspada00/restart_PCA/test/output/'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nodes", type=int, default=64)
    args = parser.parse_args()
    nodes=args.nodes
    
    mpi.stats()
    
    if mpi.rank==0:
        log('start')
        mask=Mask(filename=maskfile, maskvarname="tmask", zlevelsvar="nav_lev", ylevelsmatvar="nav_lat", xlevelsmatvar="nav_lon", dzvarname="e3t", loadtmask=True)
        log('mask loaded')
        mask_shape=mask.mask.shape[::-1]
        log(mask_shape)
        hmask=hMask(mask.mask.sum(0), local_r=local_r*2)
        #hmask=hMask(mask.mask[0,:,:], local_r=local_r)
        dec=Decomposition(hmask, nodes=nodes)
        log(dec)
        if False:
            dec.plot()
            mpi.comm.Abort()
        
        mpi.bcast(mask_shape)
        #cell=mpi.scatter(dec.cells[:mpi.size])
        cells=mpi.bcast(dec.cells[:mpi.size])
    else:
        mask_shape=mpi.bcast(None)
        #cell=mpi.scatter(None)
        cells=mpi.bcast(None)
    cell=cells[mpi.rank]
    
    log(f'mask_shape = {mask_shape}:',0)
    mask=readnc(maskfile, 'tmask', cell)
    #with Dataset(filename=maskfile, mode='r', parallel=True, comm=mpi.comm) as ds:
        #m=ds.variables['tmask']
        #if len(m.shape) == 4:
            #mask=m[0,:,cell.y_start_wb:cell.y_next_wb, cell.x_start_wb:cell.x_next_wb]
        #elif len(m.shape) == 3:
            #mask=m[:,cell.y_start_wb:cell.y_next_wb, cell.x_start_wb:cell.x_next_wb]
        #else:
            #raise Exception('strana mask')
    #mask=np.transpose(mask.data, axes=(2,1,0))
    
    bmask=mask.astype(bool)
    three2flat=np.cumsum(bmask)
    n_points=three2flat[-1]
    three2flat=three2flat.reshape(mask.shape)-1
    #xyz=np.empty([n_points,3])
    #for x in mask.shape[0]:
        #for y in mask.shape[1]:
            #for z in mask.shape[2]:
                #temp=three2flat[x,y,z]
                #if temp!=-1:
                    #xyz[temp]=[x,y,z]
                    
    two2flat=np.cumsum(bmask[...,0])
    n_points2d=two2flat[-1]
    two2flat=two2flat.reshape(mask.shape[:2])-1
    #xy=np.empty([n_points2d,2])
    #for x in mask.shape[0]:
        #for y in mask.shape[1]:
            #temp=two2flat[x,y]
            #if temp!=-1:
                #xy[temp]=[x,y]
                
    #bmask=mask.astype(bool)
    
    filelist=glob.glob(template.format(member='*', var='P1l'))
    if [member.count('P1l') for member in filelist]!=[1]*len(filelist):
        raise Exception("Bad template")
    members=[member.replace('P1l', '{var}', 1) for member in filelist]
    members.sort()
    varlist=glob.glob(members[0].format(var='*'))
    i=members[0].index('{var}')
    j=len(members[0])-i-5
    varlist=[member[i:-j] for member in varlist]
    varlist.sort()
    nvars=len(varlist)
    log(f'nvars={nvars}',0)
    
    size=len(members)
    if size*nvars != len(glob.glob(template.format(member='*', var='*'))):
        raise Exception('qualcosa non torna nel numero dei files')
    
    #filelist=glob.glob(prefix+'000.*.nc')
    #varlist=[f[len(prefix)+4:-3] for f in filelist]
    #varlist.sort()
    #nvars=len(varlist)
    #log(f'nvars={nvars}',0)
    
    #size=len(glob.glob(prefix+'*.nc'))
    #if size % nvars != 0:
        #raise Exception('qualcosa non torna nel numero dei files')
    #size=size//nvars
    
    #####attenzione linea di debug!!!!
    #size = 4
    #################
    
    log(f'size={size}',0)
    
    eigenvalues=np.empty([n_points2d,n_eigh])
    change_base=np.zeros([n_points2d,n_eigh,size])
    flat2D2flat3D=np.repeat(two2flat[...,None], bmask.shape[-1], -1)[bmask]
    P=np.empty([size,size,n_points2d])
    
    bmasks=mpi.Allgatherv(bmask)
    two2flats=mpi.Allgatherv(two2flat, shapes=[bm[...,0].shape for bm in bmasks])
    three2flats=mpi.Allgatherv(three2flat, shapes=[bm.shape for bm in bmasks])
    
    log("Allgather",0)
    
    io=IO(cells, mask_shape)
    
    #scaling=np.ones([n_points,nvars])
    scaling=np.zeros([nvars,n_points])
    for var in ['P1l','P2l','P3l','P4l']:
        scaling[varlist.index(var),three2flat[:,:,0]]=1
        
    def transform(tracer, small=10**-8):
        tracer[tracer<small]=small
        return np.log(tracer)
    
    def back(tracer):
        return np.exp(tracer)
    
    tracers=np.empty([size,nvars,n_points])
    for i, member in enumerate(members):
        #filelist=[prefix+f'{i:03}.'+var+'.nc' for var in varlist]
        filelist=[member.format(var=var) for var in varlist]
        for j,(f, var) in enumerate(zip(filelist, varlist)):
            tracers[i,j]=readnc(f,'TRN'+var,cell)[bmask]
            
    mpi.Barrier()
    log('tracers done',0)
    
    std=tracers.std(0)
    if True: 
        for i,var in enumerate(varlist):
            io.write(output+"std."+var+".nc", "TRN"+var, std[i], bmasks, three2flats)
        #io.write_all()
            
        #mpi.Barrier()
        #log('std written',0)
    
    for i in range(size):
        for j in range(nvars):
            tracers[i,j]=transform(tracers[i,j])
    
    tracers_mean=tracers.mean(0)
    
    small=np.log(10**-5)
    for i in range(len(members)):
        for j in range(len(varlist)):
            tracers[i,j][tracers[i,j]<small]=small
    tracers-=tracers.mean(0)
    
    std=tracers.std(0)
    if True: 
        for i,var in enumerate(varlist):
            io.write(output+"logstd."+var+".nc", "TRN"+var, std[i], bmasks, three2flats)
        #io.write_all()
            
        #mpi.Barrier()
        #log('logstd written',0)
    
    ########
    # disattiva questo per non scrivere le medie. Utile per fare debug.
    if True: 
        for i,var in enumerate(varlist):
            #writenc(output+"mean."+var+".nc", "TRN"+var, cell, bmask, mask_shape, tracers_mean[three2flat, i])
            io.write(output+"mean."+var+".nc", "TRN"+var, back(tracers_mean[i]), bmasks, three2flats)
        io.write_all()
            
        #mpi.Barrier()
        #log('mean written',0)
    ############
    
    #P=np.empty([size,size,n_points2d])
    scaling_sum=np.empty([n_points2d])
    for i in range(bmask.shape[0]):
        for j in range(bmask.shape[1]):
            if bmask[i,j,0]:
                start=three2flat[i,j,0]
                stop=three2flat[i,j,-1]+1
                m=tracers[...,start:stop].reshape([size,-1])
                P[:,:,two2flat[i,j]]=np.matmul(m*scaling[:,start:stop].flatten(),m.transpose())
                scaling_sum[two2flat[i,j]]=(scaling[:,start:stop]>0).sum()
                
    mpi.Barrier()
    log('P done',0)
                
    circle_template=np.zeros([local_r*2+1]*2)
    for i in range(-local_r,local_r+1):
        for j in range(-local_r,local_r+1):
            d2=(i**2+j**2)/(local_r+0.5)**2
            if d2<1:
                circle_template[local_r+i,local_r+j]=(1-d2)**2
                #circle_template[local_r+i,local_r+j]=1
            
    #eigenvalues=np.empty([n_points2d,n_eigh])
    for i in range(cell.x_start-cell.x_start_wb-edge if cell.x_start-cell.x_start_wb-edge>0 else 0 ,cell.x_next-cell.x_start_wb+edge if cell.x_next-cell.x_start_wb+edge<cell.x_next_wb-cell.x_start_wb else cell.x_next_wb-cell.x_start_wb):
        for j in range(cell.y_start-cell.y_start_wb-edge if cell.y_start-cell.y_start_wb-edge>0 else 0,cell.y_next-cell.y_start_wb+edge if cell.y_next-cell.y_start_wb+edge<cell.y_next_wb-cell.y_start_wb else cell.y_next_wb-cell.y_start_wb):
            if bmask[i,j,0]:
                circle=np.zeros(bmask.shape[:2])
                
                left=i-local_r
                left, left_template=(0,-left) if left <0 else (left,0)
                right=i+local_r+1
                right, right_template=(bmask.shape[0],bmask.shape[0]-right) if right>bmask.shape[0] else (right,0)
                
                bot=j-local_r
                bot, bot_template=(0,-bot) if bot <0 else (bot,0)
                top=j+local_r+1
                top, top_template=(bmask.shape[1],bmask.shape[1]-top) if top>bmask.shape[1] else (top,0)
                
                circle[left:right,bot:top]=circle_template[left_template:local_r*2+1+right_template,bot_template:local_r*2+1+top_template]
                bcircle=circle.astype(bool)
                circle_and_mask=np.logical_and(bcircle,bmask[...,0])
                P_loc=(P[:,:,two2flat[circle_and_mask]]*circle[circle_and_mask]**2).sum(-1)/((scaling_sum[two2flat[circle_and_mask]]*circle[circle_and_mask]**2).sum()*size)
                
                ########################
                #debug
                ########################
                if np.logical_not(np.isfinite(P_loc)).sum()>1:
                    log(f'i={i}, j={j}')
                    log('P_loc:')
                    log(P_loc)
                    log('Psum:')
                    log((P[:,:,two2flat[circle_and_mask]]*circle[circle_and_mask]).sum(-1))
                    log('denom:')
                    log((scaling_sum[two2flat[circle_and_mask]]*circle[circle_and_mask]).sum())
                    log('scaling_sum:')
                    log(scaling_sum[two2flat[circle_and_mask]])
                    log('circle:')
                    log(circle[circle_and_mask])
                    log('P[0,0]:')
                    log(P[0,two2flat[circle_and_mask]])
                    np.save("circle", bcircle.astype(np.int64))
                    np.save("bmask", bmask[...,0].astype(np.int64))
                    np.set_printoptions(threshold=10**5,linewidth=230)
                    log('bcircle:')
                    log(bcircle.astype(np.int64))
                    log('bmask[...,0]:')
                    log(bmask[...,0].astype(np.int64))
                    mpi.comm.Abort()
                ###############################
                
                eigenvalues[two2flat[i,j]], eigenvectors= eigh(P_loc,subset_by_index=[size-n_eigh,size-1],overwrite_a=True)
                eigenvectors=np.sign(np.take_along_axis(eigenvectors, np.expand_dims(np.abs(eigenvectors).argmax(0), axis=0), axis=0))*eigenvectors
                eigenvectors=eigenvectors[:,::-1]/np.sqrt(size)
                change_base[two2flat[circle_and_mask]]+=eigenvectors.transpose()*circle[circle_and_mask][:, None, None]/circle[circle_and_mask].sum()
    eigenvalues=eigenvalues[:,::-1].transpose().copy()
    
    mpi.Barrier()
    log('eigenvalues done',0)
    
    if True:
        for i in range(n_eigh):
            #writenc(output+f"eig{i:03}.nc", f"eig{i:03}", cell, bmask[...,0], mask_shape, eigenvalues[two2flat, -i-1])
            io.write(output+f"eig{i:03}.nc", f"eig{i:03}", eigenvalues[i], [bm[...,0] for bm in bmasks] , two2flats)
        #io.write_all()
        
        #mpi.Barrier()
        #log('eigenvalues written',0)
        
    if True:
        for i in range(n_eigh-1):
            #writenc(output+f"diff{i:03}_{i+1:03}.nc", f"diff{i:03}", cell, bmask[...,0], mask_shape, eigenvalues[two2flat, -i-1]-eigenvalues[two2flat, -i-2])
            #writenc(output+f"ratio{i:03}.nc", f"ratio{i:03}", cell, bmask[...,0], mask_shape, 1-eigenvalues[two2flat, -i-2]/eigenvalues[two2flat, -i-1])
            
            #io.write(output+f"diff{i:03}_{i+1:03}.nc", f"diff{i:03}", eigenvalues[i]-eigenvalues[i+1], [bm[...,0] for bm in bmasks] , two2flats)
            io.write(output+f"ratio{i:03}.nc", f"ratio{i:03}", 1 - eigenvalues[i+1]/eigenvalues[i], [bm[...,0] for bm in bmasks] , two2flats)
        #io.write_all()
            
        #mpi.Barrier()
        #log('diff and ratio written',0)
    
    if True:
        for i in range(n_eigh):
            for j,var in enumerate(varlist):
                tracer=np.zeros(n_points)
                for k in range(size):
                    #tracer[three2flat]+=tracers[k,j,three2flat]*change_base[two2flat,i,k][...,None]
                    tracer+=tracers[k,j]*change_base[flat2D2flat3D,i,k]
                io.write(output+f"Base{i:03}.{var}.nc", "TRN"+var, tracer, bmasks, three2flats)
        #io.write_all()
        
        #mpi.Barrier()
        #log('Base written',0)
    
    if True:
        ortmatrix=np.zeros([EnsSize-1,EnsSize])
        for i in range(1,EnsSize):
            den=np.sqrt(EnsSize/(i*(i+1)))
            ortmatrix[i-1,:i]=den
            ortmatrix[i-1,i]=-i*den
        ortmatrix=ortmatrix.transpose()
        
        for i in range(EnsSize):
            for j,var in enumerate(varlist):
                tracer=tracers_mean[j]
                change2member=np.matmul(ortmatrix[i],change_base[:,:EnsSize-1,:])
                for k in range(size):
                    #tracer[three2flat]+=tracers[k,j,three2flat]*change2member[two2flat,k][...,None]
                    tracer+=tracers[k,j]*change2member[flat2D2flat3D,k]
                io.write(output+f"Member{i:03}.{var}.nc", "TRN"+var, back(tracer), bmasks, three2flats)
        #io.write_all()
          
        #mpi.Barrier()
        #log('Members written',0)
        
    io.write_all()
    
    mpi.Barrier()
    log('End',0)
                
                
if __name__=='__main__':
    #try:
        main()
    #finally:
    #    mpi.comm.Abort()
                
        
