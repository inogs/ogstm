############################################################
# generates test inputs for gen_3dvar_satfloats_namelists.py
# ##########################################################
from commons import genUserDateList as DL
import numpy as np
dateDAfloat = [d.strftime("%Y%m%d-%H:%M:%S") for d in DL.getTimeList("20130101-00:00:00","20140101-00:00:00", "days=3")]
dateDAsat   = [d.strftime("%Y%m%d-%H:%M:%S") for d in DL.getTimeList("20130101-00:00:00","20140101-00:00:00", "days=7")]



dateDA = list(set(dateDAsat + dateDAfloat)) # merge
dateDA.sort()

np.savetxt('daFloat',dateDAfloat,fmt='%s')
np.savetxt('daSat'  ,dateDAsat  ,fmt='%s')
np.savetxt('daTimes',dateDA     ,fmt='%s')

