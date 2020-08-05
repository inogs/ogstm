import numpy as np
from scipy import interpolate


INPUTFILE = 'bcs/aPHY.csv'

READFILE  = np.genfromtxt(INPUTFILE, skip_header=2, skip_footer=2, delimiter='\t', missing_values='', filling_values=np.nan)

wl        = np.array( [250., 325., 350., 375., 400.,   425.,  450.,  475.,  500.,  525.,  550.,  575.,  600.,  625.,  650.,  675., 700.,
                       725., 775., 850., 950., 1050., 1150., 1250., 1350., 1450., 1550., 1650., 1750., 1900., 2200., 2900., 3700.])

wl_spec  = READFILE[:,0] 

f = interpolate.interp1d(wl_spec, READFILE[:,1:], axis=0, bounds_error=False, fill_value=0.) 

SPEC_INT = f(wl)

DINO_10  = SPEC_INT[:,0] 
DINO_100 = SPEC_INT[:,1] 
DINO_300 = SPEC_INT[:,2] 

CRYP_10  = SPEC_INT[:,3] 
CRYP_100 = SPEC_INT[:,4] 
CRYP_300 = SPEC_INT[:,5]

COCC_10  = SPEC_INT[:,6] 
COCC_100 = SPEC_INT[:,7] 
COCC_300 = SPEC_INT[:,8]

DIAT_10  = SPEC_INT[:,9] 
DIAT_100 = SPEC_INT[:,10] 
DIAT_300 = SPEC_INT[:,11]

PROC_10  = SPEC_INT[:,12] 
PROC_25  = SPEC_INT[:,13] 
PROC_100 = SPEC_INT[:,14]

SYNE_10  = SPEC_INT[:,15] 
SYNE_100 = SPEC_INT[:,16] 
SYNE_300 = SPEC_INT[:,17]

GREE_10  = SPEC_INT[:,18] 
GREE_100 = SPEC_INT[:,19] 
GREE_300 = SPEC_INT[:,20]

PICO_50  = SPEC_INT[:,21]

# DI CICCO 6 PFT GROUPS

aPFT1 = DIAT_100
aPFT2 = DINO_100
aPFT3 = CRYP_100
aPFT4 = (GREE_100 + PROC_100)*0.5
aPFT5 = (PROC_100 + SYNE_100)*0.5
aPFT6 = COCC_100

