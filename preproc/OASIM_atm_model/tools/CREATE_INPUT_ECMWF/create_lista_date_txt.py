import datetime
import numpy as np 

dt = datetime.datetime(2003, 1, 1)
end = datetime.datetime(2006, 12, 31, 23, 59, 59)
step = datetime.timedelta(days=1)

result = []

while dt < end:
    result.append(dt.strftime('%Y%m%d'))
    dt += step

np.savetxt("lista_date.txt",result,fmt='%s')
