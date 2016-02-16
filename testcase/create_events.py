import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

def create_events(test):

    os.system("mkdir -p " + test['Dir'])

    filename = test['Dir'] + '/Start_End_Times'
    f01 = open(filename,'wb')
    f01.write(test['Start'])
    f01.write("\n")
    f01.write(test['End'])
    f01.close()

    filename = test['Dir'] + '/1.aveTimes'
    f01 = open(filename,'wb')
    f01.close()

    filename = test['Dir'] + '/2.aveTimes'
    f01 = open(filename,'wb')
    f01.close()

    
    filename = test['Dir'] + '/restartTimes'
    f01 = open(filename,'wb')
    f01.close()  

    os.system("./genInputsDatelists.sh " + test['Dir'])



    
    
