#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
import numpy as np

# 10 metre U wind component, 10 metre V wind component, 2 metre dewpoint temperature, 2 metre temperature, 
# Mean sea level pressure, Surface pressure, Total cloud cover, Total column ozone, Vertical integral of cloud liquid water

server = ECMWFDataServer()

mydate=np.loadtxt("lista_date.txt",dtype=np.int)

for dd in mydate:

    yyyymmdd=str(dd)
    yyyy_mm_dd = yyyymmdd[0:4] + "-" +  yyyymmdd[4:6] + "-" + yyyymmdd[6:8]
    fileout="ERAINTERIM_" + yyyymmdd + ".nc"
    
    server.retrieve({
        "class": "ei",
        "dataset": "interim",
        "date": yyyy_mm_dd,
    #   "date": "2011-01-01/to/2011-01-02",
        "expver": "1",
        "area": "90/-180/-90/179",
        "grid": "1.0/1.0",
    #   "grid": "0.75/0.75",
        "levtype": "sfc",
        "param": "56.162/134.128/151.128/164.128/165.128/166.128/167.128/168.128/206.128",
    #   "param": "55.162/56.162/58.162/137.128/151.128/164.128/167.128/168.128/206.128/207.128",
        "step": "0",
        "stream": "oper",
        "time": "00:00:00/06:00:00/12:00:00/18:00:00",
        "type": "an",
        "format": "netcdf",
        "target": fileout,
    })
