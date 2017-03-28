# plotting hovmoeller diagram

import os,sys, getopt
import glob
import scipy.io.netcdf as NC
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from mpl_toolkits.basemap import Basemap

def plot_hovmoeller_C_CF(test):
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
	filename      = 'POSTPROC/' + test['BIO-FLOAT'] + '.nc'
	filename_phys = 'POSTPROC/' + test['BIO-FLOAT'] + '_phys.nc'

	M=NC.netcdf_file(filename,"r",mmap=False)
        nTimes   = M.dimensions['time']
	dataPO4  = (M.variables['N1p'].data[:,:]).copy()
	dataCH1  = (M.variables['P1l'].data[:,:]).copy()
	dataCH2  = (M.variables['P2l'].data[:,:]).copy()
	dataCH3  = (M.variables['P3l'].data[:,:]).copy()
	dataCH4  = (M.variables['P4l'].data[:,:]).copy()
        M.close()
        dataCH   = dataCH1 + dataCH2 + dataCH3# + dataCH4

	M=NC.netcdf_file(filename_phys,"r",mmap=False)
	dataP    = (M.variables['par'].data[:,:]).copy()
	dataK    = (M.variables['votkeavt'].data[:,:]).copy()
	dataCHL_F= (M.variables['CHL_F'].data[:,:]).copy()
        M.close()

#plot the histogram + hovmoeller 

	masked_array_PO4  = np.ma.masked_where(dataPO4>10**19,dataPO4,copy=True)
	masked_array_CH   = np.ma.masked_where(dataCH>10**19,dataCH,copy=True)
	masked_array_P    = np.ma.masked_where(dataP>10**19,dataP,copy=True)
	masked_array_K    = np.ma.masked_where(dataK>10**19,dataK,copy=True)
	masked_array_CHL_F= np.ma.masked_where(dataCHL_F>10**19,dataCHL_F,copy=True)

	fig=plt.figure(figsize=(10, 10))
# PAR
        vrn='PAR'
#       vrn_unit = '(umoles/m2/s)'
        vrn_unit = '(W/m2)'
        ax2=plt.subplot(2, 3, 1)

        t= np.transpose(np.tile(np.arange(0,nTimes),(jpk, 1)))
        z= np.tile(np.flipud(gdept),(nTimes, 1))

        data2plot= np.transpose( np.flipud(masked_array_P.T*0.217) ) # matrix must be tranposed wr2 t and z, 0.217 units ->watts/m2
#       data2plot= np.transpose( np.flipud(masked_array_P.T) ) # matrix must be tranposed wr2 t and z, -> units umoles/m2/s
	plt.pcolormesh(t,z,data2plot, norm=LogNorm(vmin=data2plot.max()/1000., vmax=data2plot.max()),cmap = 'PuBu', edgecolors = 'None')
#       plt.pcolormesh(t,z,data2plot, cmap = 'PuBu', edgecolors = 'None') # linear scale

        plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
        plt.axis([0, nTimes, 400., 0.])
        ax2.axis('tight')
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
        plt.xlabel('week', fontsize=16)
        plt.ylabel('depth [m]', fontsize=16)

##      labels=['J','F','M','A','M','J','J','A','S','O','N','D']
##      Xl    = np.arange(0.,nTimes.,30)+0.5 # Major tick position
##      xl    = np.arange(0.,nTimes.,30) # Minor tick position
##      plt.xticks(Xl,labels)
##      ax2.set_xticks(xl, minor=True)
##      plt.tick_params(axis='x',which='major',length=0)
##      plt.tick_params(which='minor',length=3)
##      plt.xticks(Xl, labels)

# vertical Eddy Diffusivity
	vrn='Vert. Eddy diff.'
        vrn_unit = '(m2/s)'
	ax2=plt.subplot(2, 3, 2)

	t= np.transpose(np.tile(np.arange(0,nTimes),(jpk, 1)))
	z= np.tile(np.flipud(gdept),(nTimes, 1))

	data2plot= np.transpose( np.flipud(masked_array_K.T) ) # matrix must be tranposed wr2 t and z
	plt.pcolormesh(t,z,data2plot, norm=LogNorm(vmin=data2plot.min(), vmax=data2plot.max()),cmap = 'BuGn', edgecolors = 'None')

	plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
	plt.axis([0, nTimes, 400., 0.])
        ax2.axis('tight')
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
	plt.xlabel('week', fontsize=16)
	plt.ylabel('depth [m]', fontsize=16)

##      labels=['J','F','M','A','M','J','J','A','S','O','N','D']
##	Xl    = np.arange(0.,nTimes.,30)+0.5 # Major tick position
##	xl    = np.arange(0.,nTimes.,30) # Minor tick position
##	plt.xticks(Xl,labels)
##	ax2.set_xticks(xl, minor=True)
##	plt.tick_params(axis='x',which='major',length=0)
##	plt.tick_params(which='minor',length=3)
##	plt.xticks(Xl, labels)

#   Legend
	ax2=plt.subplot(2, 3, 3)

	D_START = 'Init:' + test['Start'][0:8]
	D___END = 'End_:' + test['End'][0:8]
        TITLE   =  D_START + '\n' + D___END
        ax2.set_title(TITLE,fontsize=14, y=1.1)
        map = Basemap(projection='merc',lat_0=38,lon_0=14,\
                                     llcrnrlon = -5.3, \
                                     llcrnrlat = 28.0, \
                                     urcrnrlon = 37, \
                                     urcrnrlat = 46.0, \
                                     resolution='l')
# draw coastlines, country boundaries, fill continents.
        map.drawcoastlines(linewidth=0.25)
        x,y=map(test['lon0'],test['lat0'])
        map.plot(x,y,'ro',markersize=12)


# Ch l -Float
        vrn='CHL-argo-float'
        vrn_unit = '(mg chl /m3)'
        ax2=plt.subplot(2, 3, 4)

        t= np.transpose(np.tile(np.arange(0,nTimes),(jpk, 1)))
        z= np.tile(np.flipud(gdept),(nTimes, 1))

        data2plot = np.transpose( np.flipud(masked_array_CHL_F.T) ) # matrix must be tranposed wr2 t and z
        dataFLOAT = data2plot
        plt.pcolormesh(t,z,data2plot, cmap = 'BuGn', edgecolors = 'None',vmin=0.,vmax=0.75)

        plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
        plt.axis([0, nTimes, 400., 0.])
#       plt.set_clim(vmin=0.,vmax=0.75)
        ax2.axis('tight')
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
        plt.xlabel('week', fontsize=16)
        plt.ylabel('depth [m]', fontsize=16)

##      labels=['J','F','M','A','M','J','J','A','S','O','N','D']
##      Xl    = np.arange(0.,nTimes.,30)+0.5 # Major tick position
##      xl    = np.arange(0.,nTimes.,30) # Minor tick position
##      plt.xticks(Xl,labels)
##      ax2.set_xticks(xl, minor=True)
##      plt.tick_params(axis='x',which='major',length=0)
##      plt.tick_params(which='minor',length=3)
##      plt.xticks(Xl, labels)


# CHL
        vrn='Chl-Tot'
        vrn_unit = '(mg chl/m3)'
        ax2=plt.subplot(2, 3, 5)

        t= np.transpose(np.tile(np.arange(0,nTimes),(jpk, 1)))
        z= np.tile(np.flipud(gdept),(nTimes, 1))

        data2plot= np.transpose( np.flipud(masked_array_CH.T) ) # matrix must be tranposed wr2 t and z
        dataMODEL = data2plot
        plt.pcolormesh(t,z,data2plot, cmap = 'BuGn', edgecolors = 'None',vmin=0.,vmax=0.75)

        plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
        plt.axis([0, nTimes, 400., 0.])
#       plt.set_clim(vmin=0.,vmax=0.75)
        ax2.axis('tight')
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
        plt.xlabel('week', fontsize=16)
        plt.ylabel('depth [m]', fontsize=16)

##      labels=['J','F','M','A','M','J','J','A','S','O','N','D']
##      Xl    = np.arange(0.,nTimes.,30)+0.5 # Major tick position
##      xl    = np.arange(0.,nTimes.,30) # Minor tick position
##      plt.xticks(Xl,labels)
##      ax2.set_xticks(xl, minor=True)
##      plt.tick_params(axis='x',which='major',length=0)
##      plt.tick_params(which='minor',length=3)
##      plt.xticks(Xl, labels)
# CHL_FLOAT - CHL_MODEL
        vrn='FLOAT - MODEL'
        vrn_unit = '(mg chl/m3)'
        ax2=plt.subplot(2, 3, 6)

        z= np.tile(np.flipud(gdept),(nTimes, 1))

        data2plot= dataFLOAT - dataMODEL
        plt.pcolormesh(t,z,data2plot, cmap = 'RdBu', edgecolors = 'None',vmin=-0.75,vmax=0.75)

        plt.colorbar(orientation="vertical",fraction=0.07,pad=0.12)
        plt.axis([0, nTimes, 400., 0.])
#       plt.set_clim(vmin=0.,vmax=0.75)
        ax2.axis('tight')
        plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
        plt.xlabel('week', fontsize=16)
        plt.ylabel('depth [m]', fontsize=16)


	plt.tight_layout()

# Saving  image

	theOutputFile = 'POSTPROC/' + test['BIO-FLOAT'] + '_C_CF.png'
	fig.savefig(theOutputFile)
