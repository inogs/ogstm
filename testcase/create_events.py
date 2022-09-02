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
    f01.write('20000201-00:00:00')
    f01.write("\n")
    f01.write('20000301-00:00:00')
    f01.write("\n")
    f01.write('20000401-00:00:00')
    f01.write("\n")
    f01.write('20000501-00:00:00')
    f01.write("\n")
    f01.write('20000601-00:00:00')
    f01.write("\n")
    f01.write('20000701-00:00:00')
    f01.write("\n")
    f01.write('20000801-00:00:00')
    f01.write("\n")
    f01.write('20000901-00:00:00')
    f01.write("\n")
    f01.write('20001001-00:00:00')
    f01.write("\n")
    f01.write('20001101-00:00:00')
    f01.write("\n")
    f01.write('20001201-00:00:00')
    f01.write("\n")
    f01.write('20010101-00:00:00')
    f01.write("\n")
    f01.close()

    filename = test['Dir'].decode() + '/2.aveTimes'
    f01 = open(filename,'w')
    f01.write('20000201-00:00:00')
    f01.write("\n")
    f01.write('20000301-00:00:00')
    f01.write("\n")
    f01.write('20000401-00:00:00')
    f01.write("\n")
    f01.write('20000501-00:00:00')
    f01.write("\n")
    f01.write('20000601-00:00:00')
    f01.write("\n")
    f01.write('20000701-00:00:00')
    f01.write("\n")
    f01.write('20000801-00:00:00')
    f01.write("\n")
    f01.write('20000901-00:00:00')
    f01.write("\n")
    f01.write('20001001-00:00:00')
    f01.write("\n")
    f01.write('20001101-00:00:00')
    f01.write("\n")
    f01.write('20001201-00:00:00')
    f01.write("\n")
    f01.write('20010101-00:00:00')
    f01.write("\n")
    f01.close()

    filename = test['Dir'].decode() + '/restartTimes'
    f01 = open(filename,'w')
    f01.close()

    filename = test['Dir'].decode() + '/daTimes'
    f01 = open(filename,'w')
    f01.close()

    os.system("./genInputsDatelists.sh " + test['Dir'].decode())
