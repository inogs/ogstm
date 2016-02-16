#! /usr/bin/python

import numpy as np

test_conf=np.dtype([('jpi'   ,np.int)  ,('jpj',np.int)     ,('jpk',np.int),\
                    ('nprocx',np.int)  ,('nprocy',np.int)  ,\
                    ('lon0'  ,np.float),('lat0',np.float)  ,\
                    ('dx'    ,np.float),('dy'    ,np.float),\
                    ('Start' ,'S100')  ,('End'   ,'S100')  ,\
                    ('Dir','S100'),('Code','S100')])

ext_data =np.dtype([('date','S17'),('kext',np.float)])


def file2stringlist(filename):
    LIST=[]
    filein=file(filename)
    for line in filein:
        LIST.append(line[:-1])
    filein.close()
    return LIST



