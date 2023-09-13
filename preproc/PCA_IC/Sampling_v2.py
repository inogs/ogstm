import numpy as np
from netCDF4 import Dataset

import glob
import argparse

def readnc(filename, varname):
    with Dataset(filename=filename, mode='r') as ds:
        m=ds.variables[varname]
        v=np.array(m[...]).squeeze()
        
    print(filename + ' read')
        
    return v

def writenc(filename, varname, data, mask, fill_value=1.e+20):
    with Dataset(filename=filename, mode='w', compression='zlib') as nc:
        x = nc.createDimension('x', mask.shape[-1])
        y = nc.createDimension('y', mask.shape[-2])
        dims=('time',)
        if len(mask.shape)>2:
            z = nc.createDimension('z', mask.shape[-3])
            dims=('time','z','y','x')
        else:
            dims=('time','y','x')
        time = nc.createDimension('time', 1)
        v=nc.createVariable(varname, np.float64, dims, fill_value=fill_value)
        data[np.logical_not(mask)]=fill_value
        
        v[0, ...] = data[...]
    print(filename + ' written')
    
def ortmatrix3_old(EnsSize):
    ortmatrix=np.zeros([EnsSize-1,EnsSize])
    n=EnsSize//2
    for i in range(1,n):
        den=np.sqrt(EnsSize/(2*i*(i+1)))
        ortmatrix[n+i-1,:i]=den
        ortmatrix[n+i-1,i]=-i*den
    ortmatrix[:n,:n]=np.eye(n)*np.sqrt(EnsSize/2)
    ortmatrix[:n,n:2*n]=-ortmatrix[:n,:n]
    ortmatrix[n:2*n-1,n:n*2]=ortmatrix[n:2*n-1,:n]
    if EnsSize%2==1:
        ortmatrix[-1,:-1]=np.sqrt(1/(EnsSize-1))
        ortmatrix[-1,-1]=-np.sqrt(EnsSize-1)
    return ortmatrix.transpose()

def ortmatrix3(EnsSize):
    ortmatrix=np.zeros([EnsSize-1,EnsSize])
    n=EnsSize//2
    ortmatrix[0,:n]=np.sqrt(EnsSize/(2*n))
    for i in range(1,n):
        den=np.sqrt(EnsSize/(2*i*(i+1)))
        ortmatrix[i,:i]=den
        ortmatrix[i,i]=-i*den
    ortmatrix[:n,n:2*n]=-ortmatrix[:n,:n]
    
    for i in range(1,n):
        den=np.sqrt(EnsSize/(2*i*(i+1)))
        ortmatrix[n+i-1,:i]=den
        ortmatrix[n+i-1,i]=-i*den
    ortmatrix[n:2*n-1,n:n*2]=ortmatrix[n:2*n-1,:n]
    if EnsSize%2==1:
        ortmatrix[-1,:-1]=np.sqrt(1/(EnsSize-1))
        ortmatrix[-1,-1]=-np.sqrt(EnsSize-1)
    return ortmatrix.transpose()
    
def main():
    maskfile="/g100_work/OGS_devC/SSPADA/GHOSH/T1/wrkdir/MODEL/meshmask.nc"
    template='/g100_scratch/userexternal/sspada00/restart_PCA/test/output_20230908_11x20/Base{member}.{var}.nc'
    output='/g100_scratch/userexternal/sspada00/restart_PCA/test/ordine3_v2/RST{member}.20190101-00:00:00.{var}.nc'
    meanfile='/g100_scratch/userexternal/sspada00/restart_PCA/test/output_20230908_11x20/mean.{var}.nc'
    EnsSize=24
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--var")
    args = parser.parse_args()
    selvar=args.var
    
    ortmatrix=ortmatrix3(EnsSize)
    mask=readnc(maskfile,'tmask').astype(bool)
    base=np.empty((EnsSize-1,)+mask.shape)
    
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
    print(f'nvars={nvars}, varlist={" ".join(varlist)}')
    #varlist=varlist[varlist.index('Z6c'):]
    #nvars=len(varlist)
    #print(f'nvars={nvars}, varlist={varlist}')
    
    #for var in varlist:
    var=selvar
    if True:
        mean=readnc(meanfile.format(var=var),'TRN'+var)
        mean[mask==False]=0
        #mean=np.log(mean)
        for i in range(EnsSize-1):
            base[i]=readnc(template.format(member=f'{i:03}', var=var),'TRN'+var)
            base[i][mask==False]=0
            
        ensemble=mean*np.exp(np.matmul(ortmatrix, base.reshape([EnsSize-1,-1]))).reshape((EnsSize,)+mask.shape)
        
        for i in range(EnsSize):
            writenc(output.format(member=f'{i:03}', var=var), 'TRN'+var, ensemble[i], mask)
            
            
    
if __name__=='__main__':
    main()
    
    
    
