# plotting hovmoeller diagram

import os,sys, getopt
import glob
import scipy.io.netcdf as NC
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm

def plot_hovmoeller_PTSW(test):
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
	filename = 'POSTPROC/' + test['Area'] + '_phys.nc'

	M=NC.netcdf_file(filename,"r",mmap=False)
	dataP    = (M.variables['par'].data[:,:]).copy()
	dataT    = (M.variables['votemper'].data[:,:]).copy()
	dataS    = (M.variables['vosaline'].data[:,:]).copy()
	dataK    = (M.variables['votkeavt'].data[:,:]).copy()

#plot the histogram + hovmoeller 

	masked_array_P = np.ma.masked_where(dataP>10**19,dataP,copy=True)
	masked_array_T = np.ma.masked_where(dataT>10**19,dataT,copy=True)
	masked_array_S = np.ma.masked_where(dataS>10**19,dataS,copy=True)
	masked_array_K = np.ma.masked_where(dataK>10**19,dataK,copy=True)

	fig=plt.figure(figsize=(10, 10))
# PAR
        vrn='PAR'
#       vrn_unit = '(umoles/m2/s)'
        vrn_unit = '(W/m2)'
        ax2=plt.subplot(2, 2, 1)

        t= np.transpose(np.tile(np.arange(0,365),(jpk, 1)))
        z= np.tile(np.flipud(gdept),(365, 1))

        data2plot= np.transpose( np.flipud(masked_array_P.T*0.217) ) # matrix must be tranposed wr2 t and z
#       data2plot= np.transpose( np.flipud(masked_array_P.T) ) # matrix must be tranposed wr2 t and z
	plt.pcolormesh(t,z,data2plot, norm=LogNorm(vmin=data2plot.max()/1000., vmax=data2plot.max()),cmap = 'PuBu', edgecolors = 'None')
#       plt.pcolormesh(t,z,data2plot, cmap = 'PuBu', edgecolors = 'None')

        plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
        plt.axis([0, 365, 400., 0.])
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
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


# Temperature
	vrn='Temperature'
        vrn_unit = '(Degrees Celsius)'
	ax2=plt.subplot(2, 2, 2)

	t= np.transpose(np.tile(np.arange(0,365),(jpk, 1)))
	z= np.tile(np.flipud(gdept),(365, 1))

	data2plot= np.transpose( np.flipud(masked_array_T.T) ) # matrix must be tranposed wr2 t and z
	plt.pcolormesh(t,z,data2plot, cmap = 'Reds', edgecolors = 'None')

	plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
	plt.axis([0, 365, 400., 0.])
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
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

# Salinity
	vrn='Salinity'
        vrn_unit = '(PSU)'
	ax2=plt.subplot(2, 2, 3)

	t= np.transpose(np.tile(np.arange(0,365),(jpk, 1)))
	z= np.tile(np.flipud(gdept),(365, 1))

	data2plot= np.transpose( np.flipud(masked_array_S.T) ) # matrix must be tranposed wr2 t and z
	plt.pcolormesh(t,z,data2plot, cmap = 'Purples', edgecolors = 'None')

	plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
	plt.axis([0, 365, 400., 0.])
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
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

# vertical Eddy Diffusivity
	vrn='Vert. Eddy diff.'
        vrn_unit = '(m2/s)'
	ax2=plt.subplot(2, 2, 4)

	t= np.transpose(np.tile(np.arange(0,365),(jpk, 1)))
	z= np.tile(np.flipud(gdept),(365, 1))

	data2plot= np.transpose( np.flipud(masked_array_K.T) ) # matrix must be tranposed wr2 t and z
	plt.pcolormesh(t,z,data2plot, norm=LogNorm(vmin=data2plot.min(), vmax=data2plot.max()),cmap = 'BuGn', edgecolors = 'None')

	plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
	plt.axis([0, 365, 400., 0.])
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
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

	theOutputFile = 'POSTPROC/' + test['Area'] + '.png'
	fig.savefig(theOutputFile)
