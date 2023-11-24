 
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

import gc

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
    def Scatterv(sendbuf=None, shape=None, root=0, counters=None):
        t0=time.time()
        
        if mpi.rank==root:
            shapes=[sb.shape for sb in sendbuf]
        else:
            shapes=None
                
        if shape is None:
            shape=mpi.scatter(shapes,root=root)
            
            if mpi.rank==root:
                dt_name=mpi.bcast(sendbuf[0].dtype.name, root=root)
            else:
                dt_name=mpi.bcast(None, root=root)
                
            dt_np=np.dtype(dt_name)
        else:
            dt_np=sendbuf[0].dtype
            
        dt=dtlib.from_numpy_dtype(dt_np)
            
        message=np.empty(shape, dtype=dt_np)
            
        if mpi.rank==root:
            counts=np.prod(shapes,-1)
            msg=np.concatenate([sb.flatten() for sb in sendbuf])
            mpi.comm.Scatterv([msg, counts, dt], [message, message.size, dt] , root)
                               
        else:
            mpi.comm.Scatterv(None, [message, message.size, dt] , root)
            
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
        
        ysum=rect.sum(0)
        i=0
        while ysum[i]==0:
            i+=1
        self.x_start_wb+=i
        i=-1
        while ysum[i]==0:
            i-=1
        self.x_next_wb+=i+1
        
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
        
    def send_to(self,cell):
        xmax=max(self.x_start, cell.x_start_wb)
        xmin=min(self.x_next, cell.x_next_wb)
        if xmin-xmax<=0:
            return None
        ymax=max(self.y_start, cell.y_start_wb)
        ymin=min(self.y_next, cell.y_next_wb)
        if ymin-ymax<=0:
            return None
        return {'x_start':xmax, 'x_next':xmin,'y_start':ymax, 'y_next':ymin}
        
    def recv_from(self,cell):
        return cell.send_to(self)
        
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
    
    def plot(self, draw_border=True, draw_cells=True, border_color='r', cells_color='w'):
        fig,ax=plt.subplots()
        ax.imshow(self.hmask.mask.astype(np.int64), origin='lower')
        ylim=ax.get_ylim()
        xlim=ax.get_xlim()
        if draw_border:
            for cell in self.cells:
                ax.plot((cell.x_start_wb,cell.x_next_wb), (cell.y_start_wb,)*2,border_color)
                ax.plot((cell.x_start_wb,cell.x_next_wb), (cell.y_next_wb,)*2, border_color)
                ax.plot((cell.x_start_wb,)*2, (cell.y_start_wb,cell.y_next_wb),border_color)
                ax.plot((cell.x_next_wb,)*2, (cell.y_start_wb,cell.y_next_wb),border_color)
        if draw_cells:
            for cell in self.cells:
                ax.plot((cell.x_start,)*2, (cell.y_start,cell.y_next),cells_color)
                ax.plot((cell.x_next,)*2, (cell.y_start,cell.y_next),cells_color)
                ax.plot((cell.x_start,cell.x_next), (cell.y_start,)*2,cells_color)
                ax.plot((cell.x_start,cell.x_next), (cell.y_next,)*2, cells_color)
            
        ax.set_ylim(ylim)
        ax.set_xlim(xlim)
        plt.show()
            
def log(msg, node=None):
    if mpi.rank==node or node is None:
        print(f'rank = {mpi.rank} - '+str(datetime.now())+': '+str(msg))
        sys.stdout.flush()
    
def readnc(filename, varname, cell, border=False):
    with Dataset(filename=filename, mode='r') as ds:
        m=ds.variables[varname]
        if border:
            y_start=cell.y_start_wb
            y_next=cell.y_next_wb
            x_start=cell.x_start_wb
            x_next=cell.x_next_wb
        else:
            y_start=cell.y_start
            y_next=cell.y_next
            x_start=cell.x_start
            x_next=cell.x_next
        if len(m.shape) == 4:
            mask=m[0,:,y_start:y_next, x_start:x_next]
        elif len(m.shape) == 3:
            mask=m[:,y_start:y_next, x_start:x_next]
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
        else:
            data=data.transpose()
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
        if mpi.rank==self.index:
            lens=[[dim2flat.reshape([-1])[-1]+1] for dim2flat in dim2flats]            
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
                if dim2flat.shape[:2]==(cell.x_next-cell.x_start,cell.y_next-cell.y_start):
                    x_ref=cell.x_start
                    y_ref=cell.y_start
                elif dim2flat.shape[:2]==(cell.x_next_wb-cell.x_start_wb,cell.y_next_wb-cell.y_start_wb):
                    x_ref=cell.x_start_wb
                    y_ref=cell.y_start_wb
                else:
                    raise Exception(f'bad shape. shape dim2flat: {dim2flat.shape[:2]}, cell: {str(cell)}')
                v[0, ...,
                    cell.y_start:cell.y_next,
                    cell.x_start:cell.x_next] = data[...,
                                                    cell.y_start-y_ref:cell.y_next-y_ref,
                                                    cell.x_start-x_ref:cell.x_next-x_ref]
        log(self.filename + ' written')
        self.index=0
    
def main():
    maskfile="/g100_work/OGS_devC/SSPADA/GHOSH/T1/wrkdir/MODEL/meshmask.nc"
    local_r=30
    n_eigh=100
    EnsSize=24
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
        hmask=hMask(mask.mask.sum(0), local_r=0)
        #hmask=hMask(mask.mask[0,:,:], local_r=local_r)
        #hmask=hMask(mask.mask[0,:,:], local_r=0)
        dec=Decomposition(hmask, nodes=nodes)
        hmask.local_r=local_r
        for cell in dec.cells:
            cell._set_properties(hmask)
            cell.trim(hmask)        
        log(dec)
        if False:
            dec.plot(draw_border=False, draw_cells=True, cells_color='r')
            mpi.comm.Abort()
        
        mpi.bcast(mask_shape)
        cells=mpi.bcast(dec.cells[:mpi.size])
        
        mask=np.array(mask.mask[0,:,:]).transpose()-1
        boundaries=[]
        for i,cell in enumerate(cells):
            a=mask[cell.x_start:cell.x_next,cell.y_start:cell.y_next]
            a[a==0]=i
        for cell in cells:
            boundaries.append(mask[cell.x_start_wb:cell.x_next_wb,cell.y_start_wb:cell.y_next_wb].copy())
        
    else:
        mask_shape=mpi.bcast(None)
        #cell=mpi.scatter(None)
        cells=mpi.bcast(None)
        boundaries=[np.array([1])]
        
    cell=cells[mpi.rank]
    boundary2d=mpi.Scatterv(boundaries, (cell.x_next_wb-cell.x_start_wb, cell.y_next_wb-cell.y_start_wb))
    
    log(f'mask_shape = {mask_shape}:',0)
    mask=readnc(maskfile, 'tmask', cell)
    
    bmask=mask.astype(bool)
    three2flat=np.cumsum(bmask)
    n_points=three2flat[-1]
    three2flat=three2flat.reshape(mask.shape)-1
    
    #two2flat=np.cumsum(bmask[...,0])
    #n_points2d=two2flat[-1]
    #two2flat=two2flat.reshape(mask.shape[:2])-1
    
    bmask2d=boundary2d>-1
    two2flat=np.cumsum(bmask2d)
    n_points2d=two2flat[-1]
    two2flat=two2flat.reshape(bmask2d.shape)-1
    
    bound_send=[]
    bound_recv=[]
    for i in range(mpi.size):
        if i==mpi.rank:
            bound_d = None
        else:
            bound_d=cell.send_to(cells[i])
        if bound_d is None:
            bound_send.append(np.zeros([0], dtype=np.int64))
            bound_recv.append(np.zeros([0], dtype=np.int64))
        else:         
            for _ in range(2):
                bm=np.zeros(bmask2d.shape,dtype=bmask2d.dtype)
                bm[bound_d['x_start']-cell.x_start_wb:bound_d['x_next']-cell.x_start_wb,
                    bound_d['y_start']-cell.y_start_wb:bound_d['y_next']-cell.y_start_wb]=bmask2d[bound_d['x_start']-cell.x_start_wb:bound_d['x_next']-cell.x_start_wb,
                                                                                                bound_d['y_start']-cell.y_start_wb:bound_d['y_next']-cell.y_start_wb]
                bound_send.append(two2flat[bm])
                bound_d=cell.recv_from(cells[i])
            bound_recv.append(bound_send.pop())     
    
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
    log(f'size={size}',0)
    
    #ndiags=1
    #nvars+=ndiags
    
    eigenvalues=np.empty([n_points2d,n_eigh])
    #change_base=np.zeros([n_points2d,n_eigh,size])
    all_eigenvectors=np.zeros([n_points2d,n_eigh,size])
    flat2D2flat3D=np.repeat(two2flat[...,None], bmask.shape[-1], -1)[cell.x_start-cell.x_start_wb:cell.x_next-cell.x_start_wb, 
                                                                     cell.y_start-cell.y_start_wb:cell.y_next-cell.y_start_wb][bmask] 
    P=np.empty([n_points2d, size, size])
    scaling=np.zeros([nvars,n_points])
    
    bmasks=mpi.Allgatherv(bmask)
    bmask2ds=mpi.Allgatherv(bmask2d)
    two2flats=mpi.Allgatherv(two2flat, shapes=[bm.shape for bm in bmask2ds])
    three2flats=mpi.Allgatherv(three2flat, shapes=[bm.shape for bm in bmasks])
    
    mpi.Barrier()
    log("Allgather",0)
    
    io=IO(cells, mask_shape)
        
    def transform(tracer, small=10**-8):
        tracer[tracer<small]=small
        return np.log(tracer)
    
    def back(tracer):
        return np.exp(tracer)
    
    tracers=np.empty([size,nvars,n_points])
    for i, member in enumerate(members):
        filelist=[member.format(var=var) for var in varlist]
        #tracers[i,-1]=0
        for j,(f, var) in enumerate(zip(filelist, varlist)):
            tracers[i,j]=readnc(f,'TRN'+var,cell)[bmask]
            #if var in ['P1l','P2l','P3l','P4l']:
                #tracers[i,-1]+=tracers[i,j]
        
        if i%100==(size-1)%100:
            mpi.Barrier()
    log('tracers done',0)
    
    #varlist.append('P_l')
    
    #scaling=np.zeros([nvars,n_points])
    #for var in ['P1l','P2l','P3l','P4l']:
    for var in ['N1p']:
        scaling[varlist.index(var),three2flat[:,:,0]]=1
    #scaling[-1,three2flat[:,:,0]]=1
    
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
                P[two2flat[i+cell.x_start-cell.x_start_wb,j+cell.y_start-cell.y_start_wb],:,:]=np.matmul(m*scaling[:,start:stop].flatten(),m.transpose())
                scaling_sum[two2flat[i+cell.x_start-cell.x_start_wb,j+cell.y_start-cell.y_start_wb]]=(scaling[:,start:stop]>0).sum()
                
    mpi.Barrier()
    log('P done',0)
    
    def boundary_comm(flat2Darray):
        boundary_msg=[]
        for i in range(mpi.size):
            msg=flat2Darray[bound_send[i]]
            #log(f'size={msg.size}, msg={msg}')
            #log([br.shape for br in bound_recv],i)
            boundary_msg=mpi.Gatherv(msg,shapes=[br.shape+flat2Darray.shape[1:] for br in bound_recv],root=i)
            #log(f'fatto {i}',i)
            #mpi.Barrier()
            if i==mpi.rank:
                for j in range(mpi.size):
                    flat2Darray[bound_recv[j]]=boundary_msg[j]
                    
    boundary_comm(P)
    boundary_comm(scaling_sum)
                
    mpi.Barrier()
    log('P communication done',0)
                
    circle_template=np.zeros([local_r*2+1]*2)
    for i in range(-local_r,local_r+1):
        for j in range(-local_r,local_r+1):
            d2=(i**2+j**2)/(local_r+0.5)**2
            if d2<1:
                circle_template[local_r+i,local_r+j]=(1-d2)**2
                #circle_template[local_r+i,local_r+j]=1
            
    for i in range(cell.x_start-cell.x_start_wb,cell.x_next-cell.x_start_wb):
        for j in range(cell.y_start-cell.y_start_wb,cell.y_next-cell.y_start_wb):
            if bmask2d[i,j]:
                circle=np.zeros(bmask2d.shape)
                
                left=i-local_r
                left, left_template=(0,-left) if left <0 else (left,0)
                right=i+local_r+1
                right, right_template=(bmask2d.shape[0],bmask2d.shape[0]-right) if right>bmask2d.shape[0] else (right,0)
                
                bot=j-local_r
                bot, bot_template=(0,-bot) if bot <0 else (bot,0)
                top=j+local_r+1
                top, top_template=(bmask2d.shape[1],bmask2d.shape[1]-top) if top>bmask2d.shape[1] else (top,0)
                
                circle[left:right,bot:top]=circle_template[left_template:local_r*2+1+right_template,bot_template:local_r*2+1+top_template]
                bcircle=circle.astype(bool)
                circle_and_mask=np.logical_and(bcircle,bmask2d)
                P_loc=(P[two2flat[circle_and_mask]]*circle[circle_and_mask][...,None,None]**2).sum(0)/((scaling_sum[two2flat[circle_and_mask]]*circle[circle_and_mask]**2).sum()*size)
                
                eigenvalues[two2flat[i,j]], eigenvectors= eigh(P_loc,subset_by_index=[size-n_eigh,size-1],overwrite_a=True)
                eigenvectors=np.sign(np.take_along_axis(eigenvectors, np.expand_dims(np.abs(eigenvectors).argmax(0), axis=0), axis=0))*eigenvectors
                eigenvectors=eigenvectors[:,::-1]/np.sqrt(size)
                all_eigenvectors[two2flat[i,j]]=eigenvectors.transpose()
                #change_base[two2flat[circle_and_mask]]+=eigenvectors.transpose()*circle[circle_and_mask][:, None, None]/circle[circle_and_mask].sum()
                
    eigenvalues=eigenvalues[:,::-1].transpose().copy()
    
    mpi.Barrier()
    log('eigenvalues done',0)
                
    boundary_comm(all_eigenvectors)
    
    mpi.Barrier()
    log('eigenvectors communication done',0)
    
    del P
    del scaling_sum
    gc.collect()
    change_base=np.zeros([n_points2d,n_eigh,size])
    
    for i in range(cell.x_start-cell.x_start_wb,cell.x_next-cell.x_start_wb):
        for j in range(cell.y_start-cell.y_start_wb,cell.y_next-cell.y_start_wb):
            if bmask2d[i,j]:
                circle=np.zeros(bmask2d.shape)
                
                left=i-local_r
                left, left_template=(0,-left) if left <0 else (left,0)
                right=i+local_r+1
                right, right_template=(bmask2d.shape[0],bmask2d.shape[0]-right) if right>bmask2d.shape[0] else (right,0)
                
                bot=j-local_r
                bot, bot_template=(0,-bot) if bot <0 else (bot,0)
                top=j+local_r+1
                top, top_template=(bmask2d.shape[1],bmask2d.shape[1]-top) if top>bmask2d.shape[1] else (top,0)
                
                circle[left:right,bot:top]=circle_template[left_template:local_r*2+1+right_template,bot_template:local_r*2+1+top_template]
                bcircle=circle.astype(bool)
                circle_and_mask=np.logical_and(bcircle,bmask2d)
                
                change_base[two2flat[i,j]]=(all_eigenvectors[two2flat[circle_and_mask]]*circle[circle_and_mask][:, None, None]).sum(0)/circle[circle_and_mask].sum()
                
    mpi.Barrier()
    log('change_base done',0)
    
    if True:
        for i in range(n_eigh):
            io.write(output+f"eig{i:03}.nc", f"eig{i:03}", eigenvalues[i], bmask2ds, two2flats)
        #io.write_all()
        
        #mpi.Barrier()
        #log('eigenvalues written',0)
        
    if True:
        for i in range(n_eigh-1):            
            #io.write(output+f"diff{i:03}_{i+1:03}.nc", f"diff{i:03}", eigenvalues[i]-eigenvalues[i+1], bmask2ds, two2flats)
            io.write(output+f"ratio{i:03}.nc", f"ratio{i:03}", 1 - eigenvalues[i+1]/eigenvalues[i], bmask2ds, two2flats)
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
                    tracer+=tracers[k,j]*change2member[flat2D2flat3D,k]
                io.write(output+f"RST{i:03}.20190101-00:00:00.{var}.nc", "TRN"+var, back(tracer), bmasks, three2flats)
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
                
        
