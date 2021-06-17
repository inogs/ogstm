import scipy.io.netcdf as NC
import numpy as np
import netCDF4

def test(intobs,wmo,lat1,lon1):
    """
    Observation Quality check 
    """
    undef_val=0.0001   #for P_i
    #undef_val=0.005   #for N3n
    kk=()
    kk = np.where(intobs<0.)
    intobs[kk]=undef_val
    #write an output file with point below 0.
    if (len(kk[0]>1)):
       f = open("QCoutput.dat",'a')
       out1 ="%i\t%7.5f\t%7.5f\t" %(int(wmo),lat1,lon1)
       f.writelines(out1+"\n")
       out2 =str(kk[:])
       f.writelines(out2+"\n"+"\n")
       f.close

    return intobs
