#!/bin/env python

'''
Plotting the abw25_morel.dat values for absorption.
Comparison with aw values from the cited papers, as well as with Mason et al., 2016
updated absorption coefficients.
Spectral seawater absorption and total scattering coefficients in 
units of /m.  Derived from Smith and Baker 1981 (200-300 nm), and 
(730-800 nm), Morel et al 2007 (325-475), Pope and Fry 1997 (500-720), 
Curcio and Petty 1951 (800nm-2.5um), and Maul 1985 (2.5-4um). 
'''

import numpy as np

wl       = np.array([250.    , 325.     , 350.    , 375.    , 400.    , 425.     , 450.     , 475.        , 
	                 500.    , 525.     , 550.    , 575.    , 600.    , 625.     , 650.     , 675.        , 
	                 700.    , 725.     , 775.    , 850.    , 950.    , 1050.    , 1150.    , 1250.       , 
	                 1350.   , 1450.    , 1550.   , 1650.   , 1750.   , 1900.    , 2200.    , 2900.       , 3700.  ] )
 
aw_OASIM = np.array([0.6112  , 0.0218   , 0.0081  , 0.0057  , 0.0047  , 0.0049   , 0.0085   , 0.0117      , 
	                 0.0215  , 0.0407   , 0.0550  , 0.0849  , 0.1995  , 0.2850   , 0.3512   , 0.4559      , 
	                 0.6433  , 1.4449   , 2.3900  , 3.7382  , 27.4805 , 19.3470  , 67.1800  , 94.9976     ,
                   363.1256  , 1118.6070, 944.8757, 519.5995, 646.7179, 3768.5610, 2628.0830, 437623.0000 , 1338404.0000])
 
# Directly from literature 
aw_LIT   = np.array([ 0.559  , 0.01875  , 0.00910 , 0.00561 , 0.00460 , 0.00475  , 0.00922  , 0.01140     , 
	                 0.02040 , 0.0417   , 0.0565  , 0.0772  , 0.2224  , 0.2834   , 0.340    , 0.448       ,
	                  0.624  , 1.489    , 2.3900  , 3.7382  , 27.4805 , 19.3470  , 67.1800  , 94.9976     ,
                   363.1256  , 1118.6070, 944.8757, 519.5995, 646.7179, 3768.5610, 2628.0830, 437623.0000 , 1338404.0000])


aw_MASON = np.array([0.05871 , 0.00121  , 0.00089 , 0.00131 , 0.00222 , 0.00375  , 0.00808  , 0.01119     ,
	                 0.02073 , 0.040525 , 0.05629 , 0.0772  , 0.2224  , 0.2834   , 0.340    , 0.448       ,
	                  0.624  , 1.489    , 2.3900  , 3.7382  , 27.4805 , 19.3470  , 67.1800  , 94.9976     ,
                   363.1256  , 1118.6070, 944.8757, 519.5995, 646.7179, 3768.5610, 2628.0830, 437623.0000 , 1338404.0000])

# Now create the total scattering bw

bw = np.array([0.0567  , 0.0162, 0.0117, 0.0089, 0.0069, 0.0054, 0.0043, 0.0034, 0.0029, 0.0023, 0.0019, 0.0016,
               0.0014  , 0.0012, 0.0009, 0.0007, 0.0007, 0.0006, 0.0004, 0.0002, 0.0000, 0.0000, 0.0000, 0.0000, 
               0.0000  , 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000])


OUTDIR = 'bcs/'
header = '''Spectral seawater absorption and total scattering coefficients in
units of /m.  Derived from Smith and Baker 1981 (200-300 nm), and
(730-800 nm), Morel et al 2007 (325-475), Pope and Fry 1997 (500-720),
Circio and Petty 1951 (800nm-2.5um), and Maul 1985 (2.5-4um).
Format i5,f15.4,f10.4 \n'''

file_OASIM = open(OUTDIR + 'abw25_OASIM.dat', 'w')
file_LIT   = open(OUTDIR + 'abw25_LIT.dat',   'w')
file_MASON = open(OUTDIR + 'abw25_MASON.dat', 'w')


file_OASIM.write(header)
file_LIT.write(header)
file_MASON.write(header)

for i in range(len(wl)):
	file_OASIM.write('%5d%15.4f%10.4f\n'%(wl[i], aw_OASIM[i], bw[i]))
	file_LIT.write(  '%5d%15.4f%10.4f\n'%(wl[i], aw_LIT[i], bw[i]))
	file_MASON.write('%5d%15.4f%10.4f\n'%(wl[i], aw_MASON[i], bw[i]))


file_OASIM.close()
file_LIT.close()
file_MASON.close()
