# plotting Tilman diagram

import os,sys, getopt
import glob
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.io.netcdf as NC

def compute_R(mydata,depth,mu,loss,bio,alphaN,Q):
    aa     = mydata[:,depth][mu].mean(0)/mydata[:,depth][bio].mean(0)
    bb     = mydata[:,depth][loss].mean(0)/mydata[:,depth][bio].mean(0)
    return  Q/ alphaN * (aa*bb)/(aa-bb)

def plot_RL_NP(test):
#   Domain paramters
    jpi=test['jpi'];
    jpj=test['jpj'];
    jpk=test['jpk'];
    time = 1
    maskfile=test['Dir'] + '/meshmask.nc'

    M=NC.netcdf_file(maskfile,"r")

    Lon     =  M.variables['glamt'].data[0,0,:,:].copy()
    Lat     =  M.variables['gphit'].data[0,0,:,:].copy()
    gdept   =  M.variables['gdept'].data[0,:,0,0].copy()
    gdepw   =  M.variables['gdepw'].data[0,:,0,0].copy()

    M.close()

#define derived type to store informations
    mydtype = [('N1p','f4')    ,('N3n','f4')    ,('N4n','f4'),
               ('P1c','f4')    ,('P2c','f4')    ,('P3c','f4')    ,('P4c','f4'),
               ('ruPPY1c','f4'),('ruPPY2c','f4'),('ruPPY3c','f4'),('ruPPY4c','f4'),
               ('loPPY1c','f4'),('loPPY2c','f4'),('loPPY3c','f4'),('loPPY4c','f4')]          

    filename     = 'POSTPROC/' + test['Area'] + '.nc'
    filename_dia = 'POSTPROC/' + test['Area'] + '_dia.nc'

    M       = NC.netcdf_file(filename,"r",mmap=False)
    aux     = (M.variables[mydtype[0][0]].data[:,:]).copy()

    ntime= aux.shape[0]
    jpk  = aux.shape[1] 

    M_dia   = NC.netcdf_file(filename_dia,"r",mmap=False)

    mydata=np.zeros((ntime,jpk),dtype=mydtype)


# Center coordinates

    for var in mydata.dtype.names[0:7]:
       mydata[:,:][var]=M.variables[var].data[:,:].copy()
    for var in mydata.dtype.names[7:]:
       mydata[:,:][var]=M_dia.variables[var].data[:,:].copy()

    alpha_p=[0.025, 0.025, 0.25, 0.025];#m3/mgC/d
    alpha_n=[0.025, 0.025, 0.25, 0.025];#m3/mgC/d

    Qpcmin =[1.80e-4,   1.80e-4,  1.80e-4, 4.29e-4];#mmolP/mgC 
    Qncmin =[4.193e-3, 4.193e-3, 4.193e-3, 0.00687];#mmolN/mgC 

#   depth=getDepthIndex(nav_lev, 0.)
# R --> PO4 
    fig=plt.figure(figsize=(10, 10))
    for d,depth in enumerate(np.arange(0,100,5)):
        plt.subplot(4, 5, d+1)
        rp = np.zeros(4,dtype='float32')
        rn = np.zeros(4,dtype='float32')

        rp[0] =  max(0.,compute_R(mydata,depth,'ruPPY1c','loPPY1c','P1c',alpha_p[0],Qpcmin[0]))
        rp[1] =  max(0.,compute_R(mydata,depth,'ruPPY2c','loPPY2c','P2c',alpha_p[1],Qpcmin[1]))
        rp[2] =  max(0.,compute_R(mydata,depth,'ruPPY3c','loPPY3c','P3c',alpha_p[2],Qpcmin[2]))
        rp[3] =  max(0.,compute_R(mydata,depth,'ruPPY4c','loPPY4c','P4c',alpha_p[3],Qpcmin[3]))

        rn[0] =  max(0.,compute_R(mydata,depth,'ruPPY1c','loPPY1c','P1c',alpha_n[0],Qncmin[0]))
        rn[1] =  max(0.,compute_R(mydata,depth,'ruPPY2c','loPPY2c','P2c',alpha_n[1],Qncmin[1]))
        rn[2] =  max(0.,compute_R(mydata,depth,'ruPPY3c','loPPY3c','P3c',alpha_n[2],Qncmin[2]))
        rn[3] =  max(0.,compute_R(mydata,depth,'ruPPY4c','loPPY4c','P4c',alpha_n[3],Qncmin[3]))

        plt.plot(mydata[:,depth]['N1p'],(mydata[:,depth]['N3n']+mydata[:,depth]['N4n']),'.')
        plt.scatter(mydata[:,depth]['N1p'],(mydata[:,depth]['N3n']+mydata[:,depth]['N4n']),s=1)

        line_cs=['r--','g--','b--','c--'];

        for i in range(4):
            plt.plot([rp[i],rp[i]],[rn[i],  2.5], line_cs[i], lw=2)
            plt.plot([rp[i],0.15 ],[rn[i],rn[i]], line_cs[i], lw=2)
#       plt.plot([rp[i],rp[i]   ],[rn[i], 10*rn[i]], line_cs[i], lw=2)
#       plt.plot([rp[i],10*rp[i]],[rn[i],    rn[i]], line_cs[i], lw=2)
        plt.xlim([0,0.15])
        plt.ylim([0,2.5])
        plt.yticks(fontsize=7)
        plt.xticks(fontsize=7)
        if depth < 15:
           ax = plt.gca()
           ax.set_xticklabels([])
        if depth > 14:
           plt.xticks(rotation='vertical')
        plt.title(str(depth),fontsize=7)
#       plt.title(str(nav_lev[depth]),fontsize=7)
#plt.show()
    fileout="TIL" + test['Area'] + ".png"
    plt.savefig(fileout, format='png',dpi=900)
    plt.close(fig)

