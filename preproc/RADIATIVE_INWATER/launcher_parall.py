#!/bin/env python
from __future__ import print_function

import datetime, subprocess

NPROCS   = 36

dt_start = datetime.datetime(2012, 01, 01, 0, 0, 0)
dt_end   = datetime.datetime(2018, 01, 01, 0, 0, 0)

pids = []
'''
We define the number of days, the step and the residual of
the step in case the division is not exact.
'''

ndays = (dt_end - dt_start).days
step  = ndays//NPROCS
rstep = ndays%NPROCS

'''
Set the dates
'''

for ii in range(0,NPROCS):
	dt_step = datetime.timedelta(days=step + 1  if rstep > 0 else step)

	int1 = dt_start.strftime('%Y%m%d')
	int2 = (dt_start + dt_step).strftime('%Y%m%d')

	dt_start += dt_step
	rstep    -= 1 # We take one from the residual

	pids.append( subprocess.Popen(['bash','-c','ipython create_input.py %s %s' % (int1,int2)]) )
#	os.system( 'ipython create_input.py %s %s &' % (int1,int2) )

'''
Wait for everyone
'''

for pid in pids:
    pid.wait()
 
print('Ciao!')
