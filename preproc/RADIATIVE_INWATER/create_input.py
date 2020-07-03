#!/bin/env python

''' From saved .pkl files read the Profilelist and then create the input for the radiative transfer model '''

from __future__ import print_function
import pprint, pickle

pkl_file = open('Profilelist.pkl', 'rb')

Profilelist = pickle.load(pkl_file)
#pprint.pprint(Profilelist)

Floatlist = pickle.load(pkl_file)
#pprint.pprint(Floatlist)

pkl_file.close()
