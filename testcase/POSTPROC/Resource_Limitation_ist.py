# plotting Tilman diagram

import os,sys, getopt
import glob
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def compute_R(mydata,depth,mu,loss,bio,alphaN,Q):
    aa     = mydata[:,depth][mu]/mydata[:,depth][bio]
    bb     = mydata[:,depth][loss]/mydata[:,depth][bio]
    return  Q/ alphaN * (aa*bb)/(aa-bb)

from maskload import *

#define derived type to store informations
mydtype = [('N1p','f4')    ,('N3n','f4')    ,('N4n','f4'),
           ('P1c','f4')    ,('P2c','f4')    ,('P3c','f4')    ,('P4c','f4'),
           ('ruPPY1c','f4'),('ruPPY2c','f4'),('ruPPY3c','f4'),('ruPPY4c','f4'),
           ('loPPY1c','f4'),('loPPY2c','f4'),('loPPY3c','f4'),('loPPY4c','f4')]          

filename = 'ave.TEST04.nc'

M    = NC.netcdf_file(filename,"r",mmap=False)
aux  = (M.variables[mydtype[0][0]].data[:,:,2,2]).copy()
ntime= aux.shape[0]
jpk  = aux.shape[1] 

mydata=np.zeros((ntime,jpk),dtype=mydtype)


# Center coordinates
ci=jpi/2
cj=jpj/2

for var in mydata.dtype.names:
    mydata[:,:][var]=M.variables[var].data[:,:,2,2].copy()

alpha_p=[0.025, 0.025, 0.25, 0.025];#m3/mgC/d
alpha_n=[0.025, 0.025, 0.25, 0.025];#m3/mgC/d

Qpcmin =[1.80e-4,   1.80e-4,  1.80e-4, 4.29e-4];#mmolP/mgC 
Qncmin =[4.193e-3, 4.193e-3, 4.193e-3, 0.00687];#mmolN/mgC 

depth=getDepthIndex(nav_lev, 0.)
print(depth)
# R --> PO4 
rp = np.zeros(4,dtype='float32')
rn = np.zeros(4,dtype='float32')

rp[0] =  compute_R(mydata,depth,'ruPPY1c','loPPY1c','P1c',alpha_p[0],Qpcmin[0]).mean(0)
rp[1] =  compute_R(mydata,depth,'ruPPY2c','loPPY2c','P2c',alpha_p[1],Qpcmin[1]).mean(0)
rp[2] =  compute_R(mydata,depth,'ruPPY3c','loPPY3c','P3c',alpha_p[2],Qpcmin[2]).mean(0)
rp[3] =  compute_R(mydata,depth,'ruPPY4c','loPPY4c','P4c',alpha_p[3],Qpcmin[3]).mean(0)

rn[0] =  compute_R(mydata,depth,'ruPPY1c','loPPY1c','P1c',alpha_n[0],Qncmin[0]).mean(0)
rn[1] =  compute_R(mydata,depth,'ruPPY2c','loPPY2c','P2c',alpha_n[1],Qncmin[1]).mean(0)
rn[2] =  compute_R(mydata,depth,'ruPPY3c','loPPY3c','P3c',alpha_n[2],Qncmin[2]).mean(0)
rn[3] =  compute_R(mydata,depth,'ruPPY4c','loPPY4c','P4c',alpha_n[3],Qncmin[3]).mean(0)


plt.plot(mydata[:,depth]['N1p'],(mydata[:,depth]['N3n']+mydata[:,depth]['N4n']))

line_cs=['r--','g--','b--','c--'];
for i in range(4):
    plt.plot([rp[i],rp[i]],[rn[i],    10 ], line_cs[i], lw=2)
    plt.plot([rp[i],1],    [rn[i],  rn[i]], line_cs[i], lw=2)
plt.xlim([0,0.1])
plt.ylim([0,2.5])
#plt.show()
fileout="TILave.png"
plt.savefig(fileout, format='png',dpi=900)

