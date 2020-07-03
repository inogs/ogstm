import datetime
import numpy as np 

dt = datetime.datetime(2012, 01, 01, 0, 0, 0)#create from 2003
end = datetime.datetime(2017, 12, 31, 23, 59, 59)
step = datetime.timedelta(days=1)

result = []

while dt < end:
    result.append(dt.strftime('%Y%m%d'))
    dt += step

np.savetxt("lista_date.txt",result,fmt='%s')
