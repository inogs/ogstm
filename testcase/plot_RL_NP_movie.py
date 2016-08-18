# plotting Tilman diagram

import os,sys, getopt
import glob
import datetime
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.io.netcdf as NC
from matplotlib.lines import Line2D
import matplotlib.animation as animation


def compute_R(mydata,depth,mu,loss,bio,alphaN,Q):
    aa     = mydata[:,depth][mu].mean(0)/mydata[:,depth][bio].mean(0)
    bb     = mydata[:,depth][loss].mean(0)/mydata[:,depth][bio].mean(0)
    return  Q/ alphaN * (aa*bb)/(aa-bb)

def plot_RL_NP_movie(test):
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
    rp = np.zeros((4,20),dtype='float32')
    rn = np.zeros((4,20),dtype='float32')

    for d,depth in enumerate(np.arange(0,100,5)):

       rp[0,d] =  max(0.,compute_R(mydata,depth,'ruPPY1c','loPPY1c','P1c',alpha_p[0],Qpcmin[0]))
       rp[1,d] =  max(0.,compute_R(mydata,depth,'ruPPY2c','loPPY2c','P2c',alpha_p[1],Qpcmin[1]))
       rp[2,d] =  max(0.,compute_R(mydata,depth,'ruPPY3c','loPPY3c','P3c',alpha_p[2],Qpcmin[2]))
       rp[3,d] =  max(0.,compute_R(mydata,depth,'ruPPY4c','loPPY4c','P4c',alpha_p[3],Qpcmin[3]))

       rn[0,d] =  max(0.,compute_R(mydata,depth,'ruPPY1c','loPPY1c','P1c',alpha_n[0],Qncmin[0]))
       rn[1,d] =  max(0.,compute_R(mydata,depth,'ruPPY2c','loPPY2c','P2c',alpha_n[1],Qncmin[1]))
       rn[2,d] =  max(0.,compute_R(mydata,depth,'ruPPY3c','loPPY3c','P3c',alpha_n[2],Qncmin[2]))
       rn[3,d] =  max(0.,compute_R(mydata,depth,'ruPPY4c','loPPY4c','P4c',alpha_n[3],Qncmin[3]))

       line_cs=['r--','g--','c--','b--'];

    class SubplotAnimation(animation.TimedAnimation):
        def __init__(self):
            fig = plt.figure(figsize=(10, 10))
            ax=[] 
            N_P_lines=[]

            self.data=mydata

            for d,depth in enumerate(np.arange(0,100,5)):
                ax.append(fig.add_subplot(4, 5, d+1))
                N_P_lines.append(Line2D([], [], color='black'))
                ax[d].add_line(N_P_lines[d])
                ax[d].set_xlim([0,0.15])
                ax[d].set_ylim([0,2.5])

                for label in (ax[d].get_xticklabels() + ax[d].get_yticklabels()):
                     label.set_fontsize(7) 

                for i in range(4):
                    ax[d].plot([rp[i,d],rp[i,d]],[rn[i,d],    2.5], line_cs[i], lw=2)
                    ax[d].plot([rp[i,d],0.15   ],[rn[i,d],rn[i,d]], line_cs[i], lw=2)

                if depth <= 70:
                   ax[d].set_xticklabels([])  
                else:
                   for label in (ax[d].get_xticklabels()):
                       label.set_fontsize(7) 
                       label.set_rotation('vertical') 

                ax[d].set_title(str(depth),fontsize=7)


            self.lines = N_P_lines

            date = datetime.datetime(2003, 1, 1) + datetime.timedelta(0)  
            mm   = date.strftime('%m')
            dd   = date.strftime('%d')
            
            main_title = test['Area'] + ' Red--> Dia, Green-->Fla, cia -->Cia, Blue-->Dino date: ' + 'm: ' + mm + ' - d: ' + dd
            self.main_title=plt.suptitle(main_title)

            self.ax=ax

            animation.TimedAnimation.__init__(self, fig, interval=50, blit=True,repeat=False)       

        def _draw_frame(self, framedata):
            i = framedata

            for d,depth in enumerate(np.arange(0,100,5)):
                 x=self.data[0:i,depth]['N1p']
                 y=self.data[0:i,depth]['N3n']+self.data[0:i,depth]['N4n']
                 self.lines[d].set_data(x, y)
                 self.lines[d].set_label(str(i))

            date = datetime.datetime(2003, 1, 1) + datetime.timedelta(i)  
            mm   = date.strftime('%m')
            dd   = date.strftime('%d')
            
            main_title = test['Area'] + ' Red--> Dia, Green-->Fla, cia -->Cia, Blue-->Dino date: ' + 'm: ' + mm + ' - d: ' + dd

#           self._drawn_artists = self.lines
            self.main_title = plt.suptitle(main_title)

        def new_frame_seq(self):
            return iter(range(0,365))

        def _init_draw(self):
            for l in self.lines:
                l.set_data([], [])

    ani = SubplotAnimation()
    fileout="POSTPROC/MOVIE/TIL" + test['Area'] + ".mp4"
    ani.save(fileout)
#   plt.show()





