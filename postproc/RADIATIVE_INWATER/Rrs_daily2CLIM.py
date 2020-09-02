#!/bin/env ipython

from __future__ import print_function, division
import numpy as np
import os.path

from basins import V2 as OGS
from commons import timerequestors
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
import netCDF4 as NC4

INPUTDIR   =  '/gpfs/scratch/userexternal/eterzic0/RRS_SAT/ORIG/'    #'/gpfs/scratch/userexternal/eterzic0/KD490_24_DAILY/DAILY/CHECKED/'
 
OUTPUTDIR  =  '/gpfs/scratch/userexternal/eterzic0/RRS_OUT/1KM/DAILY/CLIM_MONTH/'

MASKDIR    =  '/gpfs/scratch/userexternal/eterzic0/KD490_OUT/1KM/'
 
TI         = TimeInterval('20120101','20171231', '%Y%m%d')

varlist = ['RRS412', 'RRS443', 'RRS490', 'RRS510', 'RRS555' , 'RRS670']
 
nMonths    = 12
MONTHS     = np.arange(1, nMonths + 1)

nSub       = len(OGS.med.basin_list)
 
if not os.path.exists(MASKDIR + 'masked_SUB.npy') : 

	ncin       = NC4.Dataset(TL.filelist[0], 'r')
	 
	lon_dim    = ncin.dimensions['lon'].size
	lat_dim    = ncin.dimensions['lat'].size
	 
	lon        = ncin.variables['lon'][:]
	lat        = ncin.variables['lat'][:]

	lon1, lat1 = np.meshgrid(lon, lat)

	masked_W   = np.zeros((lat_dim, lon_dim),        dtype=bool)
	masked_E   = np.zeros((lat_dim, lon_dim),        dtype=bool)
	masked_SUB = np.zeros((lat_dim, lon_dim, nSub),  dtype=bool)


	for i in range(lon_dim):
		for j in range(lat_dim):
			if OGS.wes.is_inside(lon1[j,i], lat1[j,i]):
				masked_W[j,i] = True
			if OGS.eas.is_inside(lon1[j,i], lat1[j,i]):
				masked_E[j,i] = True
			for isub, sub in enumerate(OGS.med):
				if sub.is_inside(lon1[j,i], lat1[j,i]):
					masked_SUB[j,i,isub] = True
			
	ncin.close()

	np.save(MASKDIR + 'masked_E.npy',  masked_E)
	np.save(MASKDIR + 'masked_W.npy',  masked_W)
	np.save(MASKDIR + 'masked_SUB.npy',  masked_SUB)

else:
	masked_W   = np.load(MASKDIR + 'masked_W.npy'  )
	masked_E   = np.load(MASKDIR + 'masked_E.npy'  )
	masked_SUB = np.load(MASKDIR + 'masked_SUB.npy')


for ivar, var in enumerate(varlist):

	TL         = TimeList.fromfilenames(TI, INPUTDIR,"*.nc", prefix='', dateformat="%Y%m%d", filtervar="_d-OC_CNR-L3-" + var + "-MedOC4AD4_SAM_1KM-MED-REP-v02")

	Rrs_MEAN_MED   = np.zeros((nMonths,))
	Rrs_STD_MED    = np.zeros((nMonths,))

	Rrs_MEAN_MED_W = np.zeros((nMonths,))
	Rrs_STD_MED_W  = np.zeros((nMonths,))

	Rrs_MEAN_MED_E = np.zeros((nMonths,))
	Rrs_STD_MED_E  = np.zeros((nMonths,))

	for iMonth, month in enumerate(MONTHS):  # Loop over climatological months

		CLIM_MONTH_req    = timerequestors.Clim_month(month)
	   
		ilist, wlist      = TL.select(CLIM_MONTH_req)
	   
		sum_weight        = 0.
	   
		Rrs_MEAN, Rrs_MEAN_W, Rrs_MEAN_E   = 0. , 0., 0.    
		Rrs_STD , Rrs_STD_W , Rrs_STD_E    = 0. , 0., 0.    
	 
		for ii, w in zip(ilist,wlist):

			ncin          = NC4.Dataset(TL.filelist[ii], 'r')

			Rrs_read    = np.nanmean(ncin.variables[var][0, :, :].filled(fill_value=np.nan))
			Rrs_read_W  = np.nanmean(ncin.variables[var][0,:,:][masked_W].filled(fill_value=np.nan))
			Rrs_read_E  = np.nanmean(ncin.variables[var][0,:,:][masked_E].filled(fill_value=np.nan))

			if np.isnan(Rrs_read)   : Rrs_read = 0.
			if np.isnan(Rrs_read_W) : Rrs_read_W = 0.
			if np.isnan(Rrs_read_E) : Rrs_read_E = 0.

			ncin.close()

			sum_weight  += w

			meanval_old  = Rrs_MEAN                                
			Rrs_MEAN     += w/sum_weight*(Rrs_read-meanval_old)
			Rrs_STD      += w*(Rrs_read-meanval_old)*(Rrs_read-Rrs_MEAN)

			meanval_old  = Rrs_MEAN_W                                
			Rrs_MEAN_W   += w/sum_weight*(Rrs_read_W-meanval_old)
			Rrs_STD_W    += w*(Rrs_read_W-meanval_old)*(Rrs_read_W-Rrs_MEAN_W)

			meanval_old  = Rrs_MEAN_E                                
			Rrs_MEAN_E   += w/sum_weight*(Rrs_read_E-meanval_old)
			Rrs_STD_E    += w*(Rrs_read_E-meanval_old)*(Rrs_read_E-Rrs_MEAN_E)

		Rrs_STD   = np.sqrt(Rrs_STD/sum_weight)
		Rrs_STD_W = np.sqrt(Rrs_STD_W/sum_weight)
		Rrs_STD_E = np.sqrt(Rrs_STD_E/sum_weight)

		Rrs_MEAN_MED[iMonth]   = Rrs_MEAN
		Rrs_STD_MED[iMonth]    = Rrs_STD

		Rrs_MEAN_MED_W[iMonth] = Rrs_MEAN_W
		Rrs_STD_MED_W[iMonth]  = Rrs_STD_W

		Rrs_MEAN_MED_E[iMonth] = Rrs_MEAN_E
		Rrs_STD_MED_E[iMonth]  = Rrs_STD_E

	np.save(OUTPUTDIR + var + '_MEAN_MED.npy'  ,  Rrs_MEAN_MED)
	np.save(OUTPUTDIR + var + '_STD_MED.npy'   ,  Rrs_STD_MED)
	np.save(OUTPUTDIR + var + '_MEAN_MED_W.npy',  Rrs_MEAN_MED_W)
	np.save(OUTPUTDIR + var + '_STD_MED_W.npy' ,  Rrs_STD_MED_W)
	np.save(OUTPUTDIR + var + '_MEAN_MED_E.npy',  Rrs_MEAN_MED_E)
	np.save(OUTPUTDIR + var + '_STD_MED_E.npy' ,  Rrs_STD_MED_E)

