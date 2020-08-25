#!/bin/env ipython

from __future__ import print_function, division
import numpy as np

from basins import V2 as OGS
from commons import timerequestors
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
import netCDF4 as NC4

INPUTDIR   =  '/gpfs/scratch/userexternal/eterzic0/KD490_24_DAILY/DAILY/CHECKED/'
 
OUTPUTDIR  =  '/gpfs/scratch/userexternal/eterzic0/KD490_OUT/1KM/DAILY/CLIM_MONTH/'
 
TI         = TimeInterval('20120101','20171231', '%Y%m%d')
 
TL         = TimeList.fromfilenames(TI, INPUTDIR,"*.nc", prefix='', dateformat="%Y%m%d", filtervar="_d-OC_CNR-L3-KD490-MedOC4AD4_SAM_1KM-MED-REP-v02")
 
nMonths    = 12
MONTHS     = np.arange(1, nMonths + 1)
 
ncin       = NC4.Dataset(TL.filelist[0], 'r')
 
lon_dim    = ncin.dimensions['lon'].size
lat_dim    = ncin.dimensions['lat'].size
 
lon        = ncin.variables['lon'][:]
lat        = ncin.variables['lat'][:]

lon1, lat1 = np.meshgrid(lon, lat)

masked_W   = np.zeros((lat_dim, lon_dim), dtype=bool)
masked_E   = np.zeros((lat_dim, lon_dim), dtype=bool)


for i in range(lon_dim):
	for j in range(lat_dim):
		if OGS.wes.is_inside(lon1[j,i], lat1[j,i]):
			masked_W[j,i] = True
		if OGS.eas.is_inside(lon1[j,i], lat1[j,i]):
			masked_E[j,i] = True
		
ncin.close()

Kd_MEAN_MED = np.zeros((nMonths,))
Kd_STD_MED  = np.zeros((nMonths,))

for iMonth, month in enumerate(MONTHS):  # Loop over climatological months

	CLIM_MONTH_req    = timerequestors.Clim_month(month)
   
	ilist, wlist      = TL.select(CLIM_MONTH_req)
   
	sum_weight        = 0.
   
	Kd_MEAN           = 0. #np.zeros((lat_dim, lon_dim))
	Kd_STD            = 0. #np.zeros((lat_dim, lon_dim))
 
	for ii, w in zip(ilist,wlist):

		ncin          = NC4.Dataset(TL.filelist[ii], 'r')

		Kd490_read    = np.mean(ncin.variables['KD490'][0, :, :])
		Kd490_read_W  = np.mean(ncin.variables['KD490'][0, masked_W])
		Kd490_read_E  = np.mean(ncin.variables['KD490'][0, masked_E])

		ncin.close()

		sum_weight  += w

		meanval_old  = Kd490_read                                #.copy()
		Kd_MEAN     += w/sum_weight*(Kd490_read-meanval_old)
		Kd_STD      += w*(Kd490_read-meanval_old)*(Kd490_read-Kd_MEAN)

	Kd_STD = np.sqrt(Kd_STD/sum_weight)

	Kd_MEAN_MED[iMonth] = Kd_MEAN
	Kd_STD_MED[iMonth]  = Kd_STD

np.save(OUTPUTDIR + 'Kd_MEAN_MED.npy', Kd_MEAN_MED)
np.save(OUTPUTDIR + 'Kd_MEAN_STD.npy', Kd_MEAN_STD)








