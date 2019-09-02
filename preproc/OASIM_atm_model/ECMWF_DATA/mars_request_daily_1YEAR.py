#!/usr/bin/env python
import os,sys
from ecmwfapi import ECMWFDataServer
import numpy as np

# 10 metre U wind component, 10 metre V wind component, 2 metre dewpoint temperature, 2 metre temperature, 
# Mean sea level pressure, Surface pressure, Total cloud cover, Total column ozone, Vertical integral of cloud liquid water

server        = ECMWFDataServer()

year          = sys.argv[1]

time_interval = year + "-01-01/to/" + year + "-12-31"

outfile       = 'ERAINTERIM_' + year + '.nc'

server.retrieve({
        "class": "ei",
        "dataset": "interim",
        "date": time_interval,
        "expver": "1",
        "area": "90/-180/-90/179",
        "grid": "1.0/1.0",
        "levtype": "sfc",
        "param": "56.162/134.128/151.128/164.128/165.128/166.128/167.128/168.128/206.128",
        "step": "0",
        "stream": "oper",
        "time": "00:00:00/06:00:00/12:00:00/18:00:00",
        "type": "an",
        "format": "netcdf",
        "target": outfile, })
