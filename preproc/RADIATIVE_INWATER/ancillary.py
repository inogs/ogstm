#!/bin/env python

''' Here you put all the functions, also for the IOPs , T-S corrections, etc. '''

def findVars(Varlist):
    allvars=[' CHL', 'IRR_380', 'IRR_412', 'IRR_490', 'PAR']
    if len(Varlist)==0: return False    
    
    for var in allvars:
    	if not var in Varlist:
    		return False
    return True
