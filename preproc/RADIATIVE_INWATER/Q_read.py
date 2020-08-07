import numpy as np
from scipy import interpolate


Qo_input  = 'Q_calc/Qo.csv'
SQn_input = 'Q_calc/SQn.csv'

CHL_vals = [0.03, 0.1, 0.3, 1., 3., 10.]

wavel    = [412.5, 442.5, 490., 510., 560., 620., 660.]

Qo_read  = np.genfromtxt(Qo_input,  skip_header=1, usecols=(1,2,3,4,5,6,7), delimiter='\t')
SQn_read = np.genfromtxt(SQn_input, skip_header=1, usecols=(1,2,3,4,5,6,7), delimiter='\t')  

Qo_morel  = interpolate.interp2d(wavel, CHL_vals, Qo_read)
SQn_morel = interpolate.interp2d(wavel, CHL_vals, SQn_read) 