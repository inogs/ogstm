# plotting hovmoeller diagram

import os,sys, getopt
import glob
import scipy.io.netcdf as NC
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker as ticker
from mydtype import *

def plot_COMPARE_DCM(test):
#   Domain paramters
	jpi=test['jpi'];
	jpj=test['jpj'];
	jpk=test['jpk'];
	time = 1
	maskfile=test['Dir'] + '/meshmask.nc'

	plotDepth=200

	M=NC.netcdf_file(maskfile,"r")

	Lon     =  M.variables['glamt'].data[0,0,:,:].copy()
	Lat     =  M.variables['gphit'].data[0,0,:,:].copy()
	gdept   =  M.variables['gdept'].data[0,0:plotDepth,0,0].copy()
	gdepw   =  M.variables['gdepw'].data[0,0:plotDepth,0,0].copy()

	M.close()



# Center coordinates
	ci=jpi/2
	cj=jpj/2

#extract data
#filename = 'ave.' + vrn + '.nc'

	FLOAT_LIST=['lovbio015c','lovbio017b','lovbio018c','lovbio035b','lovbio053b','lovbio066d','lovbio067c','lovbio083d','lovbio090d']
	TIME_DCM  =[          65,          20,          80,          30,          55,          15,          60,          15,          17]

	DATADIR   =['POSTPROC_REF', 'POSTPROC_Kd_1e-04', 'POSTPROC_THETA_DIA_0_5','POSTPROC_alpha-10p','POSTPROC_alpha-20p','POSTPROC_p_srsX2']
	LEGEND_LS =['REF','Dv_1E04','THETA-DIA 0.5','alpha-10p','alpha-20p','p_srsX2']

	if  test['BIO-FLOAT'] in FLOAT_LIST:
		
		dataPO4   =[]
		dataCH    =[]
		dataP     =[]
		dataK     =[]
		dataCHL_F =[]

		t = TIME_DCM [ FLOAT_LIST.index(test['BIO-FLOAT']) ]
		
		for i,inDIR in enumerate(DATADIR) :

			filename      = inDIR+ '/' + test['BIO-FLOAT'] + '.nc'
			filename_phys = inDIR+ '/' + test['BIO-FLOAT'] + '_phys.nc'

			M=NC.netcdf_file(filename,"r",mmap=False)
			nTimes   = M.dimensions['time']
			dataPO4.append( (M.variables['N1p'].data[t,0:plotDepth]).copy() )

			dataCH1  = (M.variables['P1l'].data[t,0:plotDepth]).copy()
			dataCH2  = (M.variables['P2l'].data[t,0:plotDepth]).copy()
			dataCH3  = (M.variables['P3l'].data[t,0:plotDepth]).copy()
			dataCH4  = (M.variables['P4l'].data[t,0:plotDepth]).copy()
			M.close()
			dataCH.append(dataCH1 + dataCH2 + dataCH3)#+ dataCH4

			M=NC.netcdf_file(filename_phys,"r",mmap=False)
			dataP.append( (M.variables['par'].data[t,0:plotDepth]).copy() * 0.217) # 0.217 units ->watts/m2
			dataK.append( (M.variables['votkeavt'].data[t,0:plotDepth]).copy()  )
			dataCHL_F.append( (M.variables['CHL_F'].data[t,0:plotDepth]).copy() )
			M.close()

#plot the histogram + hovmoeller 

		fig=plt.figure(figsize=(10, 10))
# PAR
		vrn='PAR'
		vrn_unit = '(W/m2)'
		ax2=plt.subplot(2, 3, 1)
		for myp in range(len(DATADIR)):

			#plt.plot(dataP[myp],gdept, label = DATADIR[myp])
			plt.plot(dataP[myp],gdept, 'k', label = LEGEND_LS[myp])

		plt.gca().invert_yaxis()
		#plt.xlabel('[$\mu W \, cm^{-2} nm^{-1}$]')
		plt.xlabel('[$W m^{-2} }$]')
		plt.ylabel('depth [m]')
		plt.xscale('log')
		plt.xlim([0.01, 1000])
		#plt.ylim([0,200])
		plt.title('Ed PAR')
#           plt.axis([0, nTimes, plotDepth, 0.])
#           ax2.axis('tight')
#           plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
#           plt.xlabel('week', fontsize=16)
#           plt.ylabel('depth [m]', fontsize=16)

# vertical Eddy Diffusivity
		vrn='Vert. Eddy diff.'
		vrn_unit = '(m2/s)'
		ax2=plt.subplot(2, 3, 2)
		for myp in range(len(DATADIR)):

			#plt.plot(dataK[myp],gdept, label = DATADIR[myp])
			plt.plot(dataK[myp],gdept, color = 'rmybck'[myp], label = LEGEND_LS[myp])

		plt.gca().invert_yaxis()
		plt.title('Vertical Eddy Diffusivity')
		plt.xlabel('[$m^{2} s^{-1}$]')
		plt.ylabel('depth [m]')
		plt.legend(loc='lower right')
		plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(4))
#           plt.axis([0, nTimes, plotDepth, 0.])
#           ax2.axis('tight')
#           plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
#    plt.xlabel('week', fontsize=16)
#    plt.ylabel('depth [m]', fontsize=16)

#   Legend
		ax2=plt.subplot(2, 3, 3)

		start_string = test['Start'][0:8]    ;  year_s 	 = start_string[0:4]  ;  	month_s	 = start_string[4:6]; 	day_s		= start_string[6:8]
		end_string   = test['End'][0:8]      ;  year_e 	 = end_string[0:4]    ;  	month_e	 = end_string[4:6]; 	day_e		= end_string[6:8]

		D_START = 'Start: ' + year_s + '-' + month_s + '-' + day_s 
		D___END = 'End:  '   + year_e + '-' + month_e + '-' + day_e 
		#D_START = 'Start:' + test['Start'][0:8]
		#D___END = 'End:' + test['End'][0:8]
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

		lat = []
		lon = []
		# this is a test. be careful :-)
		TEST_LIST_trajectories = np.loadtxt('TEST_LIST_trajectories.dat', dtype=test_conf,skiprows=1,ndmin=1)
		for test_traj in TEST_LIST_trajectories:
			if test_traj['BIO-FLOAT'] == test['BIO-FLOAT']:
				lat.append(test_traj['lat0'])
				lon.append(test_traj['lon0'])
				print test_traj['BIO-FLOAT']
		
		#for myp in range(len(FLOAT_LIST)):	

		x,y=map(lon,lat)
		map.plot(x,y,'ro',markersize=2)

# Ch l -Float
		vrn='CHL-argo-float'
		vrn_unit = '(mg chl /m3)'
		ax2=plt.subplot(2, 3, 4)

		for myp in range(len(DATADIR)):

			find_DCM = np.amax(dataCHL_F)
			DCM_index = np.argmax(dataCHL_F)
			DCM_depth = gdept[DCM_index]
			
			#plt.plot(dataCHL_F[myp],gdept, label = DATADIR[myp])
			plt.plot(dataCHL_F[myp],gdept, 'g')
		plt.axhline(y=DCM_depth,color = 'g', linestyle= '--', label = DCM_depth)
		plt.gca().invert_yaxis()
		plt.title('Chl-Float')
		plt.xlabel(r'[$mg \, m^{-3} $]')
		plt.ylabel('depth [m]')
		plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(4))
		plt.legend(loc='lower right')
#           ax2.axis('tight')
#           plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
#           plt.xlabel('week', fontsize=16)
#           plt.ylabel('depth [m]', fontsize=16)

# CHL
		vrn='Chl-Tot'
		vrn_unit = '(mg chl/m3)'
		ax2=plt.subplot(2, 3, 5)

		for myp in range(len(DATADIR)):

			find_DCM = np.amax(dataCH[myp])
			DCM_index = np.argmax(dataCH[myp])
			DCM_depth = gdept[DCM_index]
			#plt.plot(dataCH[myp],gdept, label = DATADIR[myp])
			#plt.plot(dataCH[myp],gdept, color = 'rmybck'[myp], label = LEGEND_LS[myp])
			plt.plot(dataCH[myp],gdept, color = 'rmybck'[myp])
			plt.axhline(y=DCM_depth, color = 'rmybck'[myp], linestyle= '--', linewidth=0.5, label = DCM_depth)
		plt.gca().invert_yaxis()
		plt.xlabel(r'[$mg \, m^{-3} $]')
		plt.title('Chl-Tot')
		plt.ylabel('depth [m]')
		plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(4))
		plt.legend(loc='best')
#           plt.axis([0, nTimes, plotDepth, 0.])
#           ax2.axis('tight')
#           plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
#           plt.xlabel('week', fontsize=16)
#           plt.ylabel('depth [m]', fontsize=16)

# CHL_FLOAT - CHL_MODEL
		vrn='FLOAT - MODEL'
		vrn_unit = '(mg chl/m3)'
		ax2=plt.subplot(2, 3, 6)

		for myp in range(len(DATADIR)):

			#plt.plot(dataCH[myp]-dataCHL_F[myp],gdept, label = DATADIR[myp])
			plt.plot(dataCH[myp]-dataCHL_F[myp],gdept, label = LEGEND_LS[myp])
		plt.gca().invert_yaxis()
		plt.xlabel(r'[$mg \, m^{-3} $]')
		plt.title('Chl_float - Chl_model')
		plt.ylabel('depth [m]')
		plt.gca().xaxis.set_major_locator(ticker.MaxNLocator(4))

#       plt.set_clim(vmin=0.,vmax=0.75)
#           ax2.axis('tight')
#           plt.title(vrn + '\n'+ vrn_unit, fontsize=20, y=1.1)
#           plt.xlabel('week', fontsize=16)
#           plt.ylabel('depth [m]', fontsize=16)


		plt.tight_layout()

# Saving  image

		theOutputFile =  'COMPARE_DCM/' + test['BIO-FLOAT'] + '_COMPARE_DCM.png'
		fig.savefig(theOutputFile)
