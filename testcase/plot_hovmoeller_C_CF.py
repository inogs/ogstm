# plotting hovmoeller diagram

import os,sys, getopt
import glob
import scipy.io.netcdf as NC
import netCDF4 as NC4
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
from matplotlib.colors import LogNorm
#from mpl_toolkits.basemap import Basemap
from mydtype import *

def plot_hovmoeller_C_CF(test,DATADIR,plot_mode):
#   Domain paramters
        jpi=test['jpi'];
        jpj=test['jpj'];
        jpk=test['jpk'];
        time = 1
        maskfile=test['Dir'].decode() + '/meshmask.nc'

        plotDepth=14
#       plot_mode=1
        jpk = plotDepth 
#DATADIR='POSTPROC-MLD8'

        M=NC.netcdf_file(maskfile,"r")

        Lon     =  M.variables['glamt'].data[0,0,:,:].copy()
        Lat     =  M.variables['gphit'].data[0,0,:,:].copy()
        gdept   =  M.variables['gdept'].data[0,0:plotDepth,0,0].copy()
        gdept_p1=  M.variables['gdept'].data[0,0:plotDepth+1,0,0].copy()
        gdepw   =  M.variables['gdepw'].data[0,0:plotDepth,0,0].copy()

        M.close()


# Center coordinates
        ci=jpi/2
        cj=jpj/2

#extract data
#filename = 'ave.' + vrn + '.nc'
        filename      = DATADIR+ '/' + 'TEST01' + '.nc'
        filename_phys = DATADIR+ '/' + 'TEST01' + '_phys.nc'
#       filename      = DATADIR+ '/' + test['BIO-FLOAT'] + '.nc'
#       filename_phys = DATADIR+ '/' + test['BIO-FLOAT'] + '_phys.nc'

        M=NC4.Dataset(filename,"r")
        nTimes   = M.dimensions['time'].size
        dataPO4  = (M.variables['N1p'][:,0:plotDepth]).copy()
        dataCH1  = (M.variables['P1l'][:,0:plotDepth]).copy()
        dataCH2  = (M.variables['P2l'][:,0:plotDepth]).copy()
        dataCH3  = (M.variables['P3l'][:,0:plotDepth]).copy()
        dataCH4  = (M.variables['P4l'][:,0:plotDepth]).copy()
        dataBI1  = (M.variables['P1c'][:,0:plotDepth]).copy()
        dataBI2  = (M.variables['P2c'][:,0:plotDepth]).copy()
        dataBI3  = (M.variables['P3c'][:,0:plotDepth]).copy()
        dataBI4  = (M.variables['P4c'][:,0:plotDepth]).copy()
        M.close()
        dataCH   = dataCH1 + dataCH2 + dataCH3 + dataCH4
        dataBI   = dataBI1 + dataBI2 + dataBI3 + dataBI4
        dataCHs  = dataCH.sum(1)
        dataCHm  = dataCH.max(1)

        M=NC4.Dataset(filename_phys,"r")
        dataP      = (M.variables['par'][:,0:plotDepth]).copy()
        dataK      = (M.variables['votkeavt'][:,0:plotDepth]).copy()
        dataCHL_F  = (M.variables['CHL_F'][:,0:plotDepth]).copy()
        dataCHL_Fs = dataCHL_F.sum(1)
        dataCHL_Fm = dataCHL_F.max(1)
#       lon_BGC_float = (M.variables['lon_BGC_float'][:,0]).copy()
#       lat_BGC_float = (M.variables['lat_BGC_float'][:,0]).copy()
        M.close()

#plot the histogram + hovmoeller 

        masked_array_PO4  = np.ma.masked_where(dataPO4>10**19,dataPO4,copy=True)
        masked_array_CH   = np.ma.masked_where(dataCH>10**19,dataCH,copy=True)
        masked_array_P    = np.ma.masked_where(dataP>10**19,dataP,copy=True)
        masked_array_K    = np.ma.masked_where(dataK>10**19,dataK,copy=True)
        masked_array_CHL_F= np.ma.masked_where(dataCHL_F>10**19,dataCHL_F,copy=True)

        eu_depth      = []
        DCM           = []
        DCM_BGC_FLOAT = []
        PO4_DCM       = []
        BIO_DCM       = []
        for i in range(nTimes):

            eu_depth.append( np.max( np.argwhere(dataP[i,:]>dataP[i,0]/100.) )  + 1 )        
            DCM.append(np.argmax(dataCH[i,:])  + 1 )
            DCM_BGC_FLOAT.append(np.argmax(dataCHL_F[i,:])  + 1 )
            PO4_DCM.append(dataPO4[i,DCM[i]-1])
            BIO_DCM.append(dataBI[i,DCM[i]-1])

        eu_depth= np.asarray( eu_depth )
        DCM           = np.asarray( DCM )
        DCM_BGC_FLOAT = np.asarray( DCM_BGC_FLOAT )
        PO4_DCM       = np.asarray( PO4_DCM )
        BIO_DCM       = np.asarray( BIO_DCM )


        fig=plt.figure(figsize=(10, 10))
# PAR
        vrn='PAR'
        vrn_unit = '(umoles/m2/s)'
#       vrn_unit = '(W/m2)'
        ax2=plt.subplot(2, 2, 1)

        t     = np.transpose(np.tile(np.arange(0,nTimes),(jpk, 1)))
        z     = np.tile(np.flipud(gdept),(nTimes, 1))
        data2plot= np.transpose( np.flipud(masked_array_P.T) ) # matrix must be tranposed wr2 t and z, 0.217 units ->watts/m2

        plt.pcolormesh(t[:,:],z[:,:],data2plot[:,:], norm=LogNorm(vmin=data2plot[:,:].max()/1000., vmax=data2plot[:,:].max()),cmap = 'PuBu_r', edgecolors = 'None')
#       plt.plot( t[:,0], eu_depth, 'k',linewidth=4 )
#       plt.plot( t[:,0], eu_depth, 'w',linewidth=1 )
       
#       plt.pcolormesh(t,z,data2plot, cmap = 'PuBu', edgecolors = 'None') # linear scale
#       plt.pcolormesh(t,z,data2plot, cmap = 'PuBu', edgecolors = 'None') # linear scale

        plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
        plt.axis([0, nTimes, plotDepth, 0.])
        ax2.axis('tight')
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
        plt.xlabel('week', fontsize=16)
        plt.ylabel('depth [m]', fontsize=16)

        vrn='Vert. Eddy diff.'
        vrn_unit = '(m2/s)'
        ax2=plt.subplot(2, 2, 2)

        t= np.transpose(np.tile(np.arange(0,nTimes),(jpk, 1)))
        z= np.tile(np.flipud(gdept),(nTimes, 1))

        data2plot= np.transpose( np.flipud(masked_array_K.T) ) # matrix must be tranposed wr2 t and z
        plt.pcolormesh(t[:,:],z[:,:],data2plot[:,:], norm=LogNorm(vmin=data2plot.min(), vmax=data2plot.max()),cmap = 'BuGn', edgecolors = 'None')
#       plt.plot( t[:,0], eu_depth, 'k',linewidth=4 )
#       plt.plot( t[:,0], eu_depth, 'w',linewidth=1 )
        plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
        plt.axis([0, nTimes, plotDepth, 0.])
        ax2.axis('tight')
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
        plt.xlabel('week', fontsize=16)
        plt.ylabel('depth [m]', fontsize=16)

# CHL
        vrn='Chl-Tot Model'
        vrn_unit = '(mg chl/m3)'
        ax2=plt.subplot(2, 2, 3)
        trans = mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)

        t= np.transpose(np.tile(np.arange(0,nTimes),(jpk, 1)))
        z= np.tile(np.flipud(gdept),(nTimes, 1))
        norms = np.transpose(np.tile(dataCHs,(jpk, 1)))
        normm = np.transpose(np.tile(dataCHm,(jpk, 1)))

        data2plot= np.transpose( np.flipud(masked_array_CH.T) ) # matrix must be tranposed wr2 t and z
        norm     = np.ones(data2plot.shape)
        vmax     = 0.7
        if plot_mode == 1:
           norm = norms
           vmax = (data2plot/norm).max()
        if plot_mode == 2:
           norm = normm
           vmax = (data2plot/norm).max()

        plt.pcolormesh(t[:,:],z[:,:],data2plot[:,:]/norm, cmap = 'nipy_spectral', edgecolors = 'None',vmin=0.,vmax=vmax)
#       plt.plot( t[:,0], eu_depth, 'k',linewidth=4 )
#       plt.plot( t[:,0], eu_depth, 'w',linewidth=1 )
#       ax2.fill_between(t[:,0], 199., 0.5, where= (dmax < 50.) , facecolor='gray', alpha=0.75, transform=trans)
        dataMODEL = data2plot/norm
        ax2.set_ylim(199,1)
        ax2.set_xlim(0., nTimes)

        plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
#       plt.axis([0, nTimes, plotDepth, 0.])
#       plt.set_clim(vmin=0.,vmax=0.75)
        ax2.axis('tight')
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
        plt.xlabel('month', fontsize=16)
        plt.ylabel('depth [m]', fontsize=16)

# Nutrients
        vrn='Phosphate Model'
        vrn_unit = '(mmol P/m3)'
        ax2=plt.subplot(2, 2, 4)

        t= np.transpose(np.tile(np.arange(0,nTimes),(jpk, 1)))
        z= np.tile(np.flipud(gdept),(nTimes, 1))

        data2plot= np.transpose( np.flipud(masked_array_PO4.T) ) # matrix must be tranposed wr2 t and z
        plt.pcolormesh(t,z,data2plot, cmap = 'BuGn', edgecolors = 'None')
#       plt.plot( t[:,0], eu_depth, 'k',linewidth=4 )
#       plt.plot( t[:,0], eu_depth, 'w',linewidth=1 )

        plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
        plt.axis([0, nTimes, 400., 0.])
        ax2.axis('tight')
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
        plt.xlabel('month', fontsize=16)
        plt.ylabel('depth [m]', fontsize=16)

        plt.tight_layout()

# Saving  image

        theOutputFile = DATADIR + '/' + 'TEST01' + '_C_CF_m' + str(plot_mode) +'.png'
#       theOutputFile = DATADIR + '/' + test['BIO-FLOAT'] + '_C_CF_m' + str(plot_mode) +'.png'
        fig.savefig(theOutputFile)
# Saving data for other plots
        theOutputFile = DATADIR + '/' + 'TEST01' + '.txt'

        ff = open( theOutputFile, "w" )

        for i in range(nTimes):
            ff.write(str(eu_depth[i]))
            ff.write(" ")
            ff.write(str(DCM[i]))
            ff.write(" ")
            ff.write(str(DCM_BGC_FLOAT[i]))
            ff.write(" ")
            ff.write(str(PO4_DCM[i]))
            ff.write(" ")
            ff.write(str(BIO_DCM[i]))
            ff.write("\n")

        ff.close()
