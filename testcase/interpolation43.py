import numpy as np
def interpolate(var43, nPoints):
   x = np.arange(43)
   xNew=np.arange(nPoints,dtype=nfloat)/nPoints*43
   return np.interp(xNew,x,var43)
