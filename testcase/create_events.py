import os,sys

import numpy as np

from mydtype import *

import scipy.io.netcdf as NC

import pickle

def create_events(test):

    os.system("mkdir -p " + test['Dir'].decode())

    filename = test['Dir'].decode() + '/Start_End_Times'
    f01 = open(filename,'w')
    f01.write(test['Start'].decode())
    f01.write("\n")
    f01.write(test['End'].decode())
    f01.close()

    filename = test['Dir'].decode() + '/1.aveTimes'
    f01 = open(filename,'w')
    f01.close()

    filename = test['Dir'].decode() + '/2.aveTimes'
    f01 = open(filename,'w')
    f01.close()


    filename = test['Dir'].decode() + '/restartTimes'
    f01 = open(filename,'w')
    f01.close()

    os.system("./genInputsDatelists.sh " + test['Dir'].decode())
