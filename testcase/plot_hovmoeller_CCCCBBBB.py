# plotting hovmoeller diagram

import os,sys, getopt
import glob
import scipy.io.netcdf as NC
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm

def plot_hovmoeller_CCCCBBBB(test):
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


# Center coordinates
	ci=jpi/2
	cj=jpj/2

#extract data
#filename = 'ave.' + vrn + '.nc'
	filename      = 'POSTPROC/' + test['Area'] + '.nc'
	filename_phys = 'POSTPROC/' + test['Area'] + '_phys.nc'

	M=NC.netcdf_file(filename,"r",mmap=False)
	dataCH1  = (M.variables['P1l'].data[:,:]).copy()
	dataCH2  = (M.variables['P2l'].data[:,:]).copy()
	dataCH3  = (M.variables['P3l'].data[:,:]).copy()
	dataCH4  = (M.variables['P4l'].data[:,:]).copy()
	dataB1   = (M.variables['P1c'].data[:,:]).copy()
	dataB2   = (M.variables['P2c'].data[:,:]).copy()
	dataB3   = (M.variables['P3c'].data[:,:]).copy()
	dataB4   = (M.variables['P4c'].data[:,:]).copy()
        M.close()

#plot the histogram + hovmoeller 

	masked_array_CH1 = np.ma.masked_where(dataCH1>10**19,dataCH1,copy=True)
	masked_array_CH2 = np.ma.masked_where(dataCH2>10**19,dataCH2,copy=True)
	masked_array_CH3 = np.ma.masked_where(dataCH3>10**19,dataCH3,copy=True)
	masked_array_CH4 = np.ma.masked_where(dataCH4>10**19,dataCH4,copy=True)

	masked_array_B1  = np.ma.masked_where(dataB1>10**19,dataB1,copy=True)
	masked_array_B2  = np.ma.masked_where(dataB2>10**19,dataB2,copy=True)
	masked_array_B3  = np.ma.masked_where(dataB3>10**19,dataB3,copy=True)
	masked_array_B4  = np.ma.masked_where(dataB4>10**19,dataB4,copy=True)

	fig=plt.figure(figsize=(10, 10))
# CHL
        vrnLIST=['Chl-Diatoms','Chl-Flagellates','Chl-Phyto','Chl-Dinoflagellates']
        vrn_unitLIST = ['(mg chl/m3)','(mg chl/m3)','(mg chl/m3)','(mg chl/m3)']
        for v,vrn in enumerate(vrnLIST):
            ax2=plt.subplot(2, 4, v+1)

            t= np.transpose(np.tile(np.arange(0,365),(jpk, 1)))
            z= np.tile(np.flipud(gdept),(365, 1))

            if v == 0 : data2plot= np.transpose( np.flipud(masked_array_CH1.T) ) # matrix must be tranposed wr2 t and z
            if v == 1 : data2plot= np.transpose( np.flipud(masked_array_CH2.T) ) 
            if v == 2 : data2plot= np.transpose( np.flipud(masked_array_CH3.T) )
            if v == 3 : data2plot= np.transpose( np.flipud(masked_array_CH4.T) )
            plt.pcolormesh(t,z,data2plot, cmap = 'BuGn', edgecolors = 'None')

            plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
            plt.axis([0, 365, 400., 0.])
            plt.title(vrn + '\n'+ vrn_unitLIST[v], fontsize=20, y=1.1)
            plt.xlabel('month', fontsize=16)
            plt.ylabel('depth [m]', fontsize=16)
            labels=['J','F','M','A','M','J','J','A','S','O','N','D']
            Xl    = np.arange(0.,365.,30)+0.5 # Major tick position
            xl    = np.arange(0.,365.,30) # Minor tick position
            plt.xticks(Xl,labels)
            ax2.set_xticks(xl, minor=True)
            plt.tick_params(axis='x',which='major',length=0)
            plt.tick_params(which='minor',length=3)
            plt.xticks(Xl, labels)

# Carbon
        vrnLIST=['Bio-Diatoms','Bio-Flagellates','Bio-Phyto','Bio-Dinoflagellates']
        vrn_unitLIST = ['(mg C/m3)','(mg C/m3)','(mg C/m3)','(mg C/m3)']
        for v,vrn in enumerate(vrnLIST):
            ax2=plt.subplot(2, 4, v+5)

            t= np.transpose(np.tile(np.arange(0,365),(jpk, 1)))
            z= np.tile(np.flipud(gdept),(365, 1))

            if v == 0 : data2plot= np.transpose( np.flipud(masked_array_B1.T) ) # matrix must be tranposed wr2 t and z
            if v == 1 : data2plot= np.transpose( np.flipud(masked_array_B2.T) ) 
            if v == 2 : data2plot= np.transpose( np.flipud(masked_array_B3.T) ) 
            if v == 3 : data2plot= np.transpose( np.flipud(masked_array_B4.T) ) 
            plt.pcolormesh(t,z,data2plot, cmap = 'Blues', edgecolors = 'None')

            plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
            plt.axis([0, 365, 400., 0.])
            plt.title(vrn + '\n'+ vrn_unitLIST[v], fontsize=20, y=1.1)
            plt.xlabel('month', fontsize=16)
            plt.ylabel('depth [m]', fontsize=16)
            labels=['J','F','M','A','M','J','J','A','S','O','N','D']
            Xl    = np.arange(0.,365.,30)+0.5 # Major tick position
            xl    = np.arange(0.,365.,30) # Minor tick position
            plt.xticks(Xl,labels)
            ax2.set_xticks(xl, minor=True)
            plt.tick_params(axis='x',which='major',length=0)
            plt.tick_params(which='minor',length=3)
            plt.xticks(Xl, labels)

	plt.tight_layout()

# Saving  image

	theOutputFile = 'POSTPROC/' + test['Area'] + '_CCCCBBBB.png'
	fig.savefig(theOutputFile)
