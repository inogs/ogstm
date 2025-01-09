#! /usr/bin/python

import numpy as np

test_conf=np.dtype([('jpi'   ,int)  ,('jpj',int)     ,('jpk',int),\
                    ('nprocx',int)  ,('nprocy',int)  ,\
                    ('lon0'  ,float),('lat0',float)  ,\
                    ('dx'    ,float),('dy'    ,float),\
                    ('Start' ,'S100')  ,('End'   ,'S100')  ,\
                    ('Dir','S100'),('Code','S100')])

ext_data =np.dtype([('date','S17'),('kext',float)])


def file2stringlist(filename):
    LIST=[]
    filein=open(filename)
    for line in filein:
        LIST.append(line[:-1])
    filein.close()
    return LIST
