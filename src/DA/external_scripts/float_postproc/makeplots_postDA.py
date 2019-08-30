import scipy.io.netcdf as NC
import numpy as np
import pylab as pl
#from profiler import *
import basins.OGS as OGS
from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
from commons.layer import Layer
from commons.mask import Mask


def plot_floatvsmodel(modelvarname,idate1,AllProfiles,AllChl, NewPres, Float,wmolist):
        TheMask=Mask("/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc")
        nav_lev = TheMask.zlevels
        dz= TheMask.dz

# AllProfiles matrix: len(Float),dimens=201 levelli, 6
# 6 = 0 Model, 1 Ref, 2 Misfit,3 Lat,4 Lon,5 ID FLOAT
        total_index=0
        tot_shape=0
        sigmab=0
        sigmao=0
        #shapecorr=np.matrix([[4,10,14,17,20], [1,2,4,8,16]])
        shapecorr=np.matrix([[0,9,13,17,20],[9,13,17,20,22],[1,2,4,8,16]])

        fig = pl.figure(wmolist)
# EFFETTUA CALCOLI DI VALORE MEDIO DEL MODELLO-OBS
        pmeanmodel =(np.mean(AllProfiles[:,:,0],0))   #mean model
        pmeanobs   =(np.mean(AllProfiles[:,:,1],0))   #mean obs
        pmeanmisfit=(np.mean(AllProfiles[:,:,2],0))   #mean misfit

# EFFETTUA CALCOLI DI STD DEL MODELLO-OBS
        pstdmodel =(np.std(AllProfiles[:,:,0],0))     #std model
        pstdobs   =(np.std(AllProfiles[:,:,1],0))     #std obs
        pstdmisfit=(np.std(AllProfiles[:,:,2],0))     #std misfit

# EFFETTUA CALCOLI DI VALORE MEDIO DELLA CORREZIONE e STD
        chl_mean=(np.mean(AllChl[:,:],0))        #mean correzione
        chl_std =(np.std(AllChl[:,:],0))         #std correzione

#INTERPOLAZIONE SUI 26LEV DEL 3DVAR 
        interp_pmeanmodel =np.interp(nav_lev[0:26],NewPres,pmeanmodel)
        interp_pmeanobs   =np.interp(nav_lev[0:26],NewPres,pmeanobs)
        interp_pmeanmisfit=np.interp(nav_lev[0:26],NewPres,pmeanmisfit)

#CALCOLO ASSIMILATO
        assimilato=interp_pmeanmodel+chl_mean
        sumindex=0
        sumden=0
        sumintcorr=0
        Z=np.sum(dz[0:26])
        #sumden=0
        for i in np.arange(0,26):
            sumindex+=(abs(AllChl[0,i]-interp_pmeanmisfit[i]))*(dz[i]/Z)
            sumden+=(abs(interp_pmeanmisfit[i]))*(dz[i]/Z)

            sumintcorr+=(abs(AllChl[0,i])*dz[i])/Z
            sigmab+=abs(assimilato[i]-interp_pmeanmodel[i])*abs(interp_pmeanobs[i]-interp_pmeanmodel[i])*(dz[i]/Z)
            sigmao+=abs(interp_pmeanobs[i]-assimilato[i])*abs(interp_pmeanobs[i]-interp_pmeanmodel[i])*(dz[i]/Z)

        total_index=1-(sumindex/sumden)   #IIN
 
        for j in np.arange(0,5):
            j1=shapecorr[0,j]
            j2=shapecorr[1,j]
            submean=np.mean(chl_mean[j1:j2])
            if(submean>0):
            #if(chl_mean[jindex]>0):
               tot_shape+=shapecorr[2,j]
               #print chl_mean[jindex],shapecorr[0,j],shapecorr[1,j]
               
        #print wmolist,'IIN=',total_index, 'INT=',sumintcorr, 'ISC',tot_shape
        print idate1,wmolist,total_index,sumintcorr,tot_shape,sigmab,sigmao

        coordinate=str(AllProfiles[0,0,3])[:6]+';'+str(AllProfiles[0,0,4])[:6]
#FACCIO I PLOT DI TUTTI I PROFILI E IL VALOR MEDIO PER MODEL E OBS SOLO PER I FLOAT CON PROFILI

        pl.subplot(1, 2, 1)
        pl.plot(interp_pmeanmodel,nav_lev[0:26],'b',label='Model',linewidth=2)        #MODEL
        pl.plot(interp_pmeanobs  ,nav_lev[0:26],'r',label='Float',linewidth=2)        #OBS
        pl.plot(assimilato       ,nav_lev[0:26],'k',label='Assimilated',linewidth=2)  #ASSIMILATO
        pl.xticks(np.arange(-0.20,0.70,0.10), fontsize = 8)
        #pl.xticks(np.arange(-5.00,5.50,1.00), fontsize = 8)
        pl.ylim([0,200])
        pl.gca().invert_yaxis()
        pl.yticks(fontsize = 8)
        pl.axvline(x=0,color='k',linewidth=0.5)
        floatlabel = 'Float '+ str(wmolist) + ' (' +coordinate+')'
        fig.suptitle(floatlabel)
        pl.title('FLOAT and MODEL',fontsize=13)
        pl.legend(loc='lower right',fontsize=10)
        pl.xlabel('Chl $[mg/m^3]$')
        #pl.xlabel('NO3 $[mg/m^3]$')
        pl.ylabel('Depth $[m]$')

        pl.subplot(1, 2, 2)
        pl.plot(chl_mean ,nav_lev[0:26],'m',label='Correction',linewidth=2)   #CORREZIONE
        pl.plot(interp_pmeanmisfit,nav_lev[0:26],'g',label='Misfit',linewidth=2)  #MISFIT
        pl.legend(loc='lower right',fontsize=10)
        pl.xticks(np.arange(-0.20,0.70,0.10), fontsize = 8)
        #pl.xticks(np.arange(-5.00,5.50,1.00), fontsize = 8)
        pl.ylim([0,200])
        pl.gca().invert_yaxis()
        pl.yticks(fontsize = 8)
        pl.axvline(x=0,color='k',linewidth=0.5)
        pl.title('CORRECTION 3DVAR',fontsize=13)

        fig_name = ''.join(['png_post/'+modelvarname+'_Week_'+idate1+'_'+wmolist,'.png'])
        fig.suptitle(floatlabel)

        pl.savefig(fig_name)
#       pl.show(block=False)
        pl.close(fig)
        return  # END OF THE FUNCTION 
