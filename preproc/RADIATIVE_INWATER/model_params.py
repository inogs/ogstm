from commons.time_interval import TimeInterval
import os
command = 'source ../../compilers/machine_modules/galileo.intel'
os.system(command)

wl =np.array( [250., 325., 350., 375., 400.,   425.,  450.,  475.,  500.,  525.,  550.,  575.,  600.,  625.,  650.,  675., 700.,
               725., 775., 850., 950., 1050., 1150., 1250., 1350., 1450., 1550., 1650., 1750., 1900., 2200., 2900., 3700.])
wl_int = wl.astype(np.int64)

wl_str = [np.str(wl_int[i]) for i in range(len(wl_int))] 
str_Ed = ['Ed_' + np.str(wl_str[i]) for i in range(len(wl_str))]
str_Es = ['Es_' + np.str(wl_str[i]) for i in range(len(wl_str))]

maskfile    = '/galileo/home/userexternal/eterzic0/OASIM_POSTPROC/ARGO_MATCHUP/CODES/meshmask.nc'
ncIN=NC.netcdf_file(maskfile,"r")
nav_lev = ncIN.variables['nav_lev'].data.copy()
ncIN.close()