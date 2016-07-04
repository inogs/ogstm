# plotting Tilman diagram

import os,sys, getopt
import glob
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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

depth=getDepthIndex(nav_lev, 50.)
print(depth)
# R --> PO4 
rp = np.zeros(4,dtype='float32')
rn = np.zeros(4,dtype='float32')

mu = [mydata[:,depth]['ruPPY1c'].mean(0)/mydata[:,depth]['P1c'].mean(0), 
      mydata[:,depth]['ruPPY2c'].mean(0)/mydata[:,depth]['P2c'].mean(0),
      mydata[:,depth]['ruPPY3c'].mean(0)/mydata[:,depth]['P3c'].mean(0),
      mydata[:,depth]['ruPPY4c'].mean(0)/mydata[:,depth]['P4c'].mean(0)]; 

dd = [mydata[:,depth]['loPPY1c'].mean(0)/mydata[:,depth]['P1c'].mean(0),    
      mydata[:,depth]['loPPY2c'].mean(0)/mydata[:,depth]['P2c'].mean(0),
      mydata[:,depth]['loPPY3c'].mean(0)/mydata[:,depth]['P3c'].mean(0),
      mydata[:,depth]['loPPY4c'].mean(0)/mydata[:,depth]['P4c'].mean(0)];

for i in range(4):
     rp[i] = Qpcmin[i]/ alpha_p[i] * (mu[i]*dd[i])/(mu[i]-dd[i])
     rn[i] = Qncmin[i]/ alpha_n[i] * (mu[i]*dd[i])/(mu[i]-dd[i])


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

