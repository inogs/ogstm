# plotting hovmoeller diagram

import os,sys, getopt
import glob
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from maskload import *

# Center coordinates
ci=jpi/2
cj=jpj/2

#extract data
vrn      = 'P1l'
filename = 'ave.' + vrn + '.nc'

M=NC.netcdf_file(filename,"r",mmap=False)
data     = (M.variables[vrn].data[108:120,:,1,1]).copy()

#plot the histogram + hovmoeller 

masked_array = np.ma.masked_where(data>10**19,data,copy=True)

fig=plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 1,
                   width_ratios=[1,1],
                   height_ratios=[1,5]
                   )
gs.update(wspace = 0, hspace = 0,bottom=0.38)

ax1=plt.subplot(gs[0]) # Histogram
data2plot= (masked_array*e3t[0,:,0,0]).sum(1) # compute integral
ax1.bar(np.arange(0,12), data2plot, alpha=0.4, color='black')
ax1.set_xlim((0, 12))
ax1.axis('off')
ax1.set_xticklabels(())
ax1.set_yticklabels(())

ax2=plt.subplot(gs[1]) # Hovmoeller

t= np.transpose(np.tile(np.arange(0,13),(jpk, 1)))
z= np.tile(np.flipud(nav_lev),(13, 1))

data2plot= np.transpose( np.flipud(masked_array.T) ) # matrix must be tranposed wr2 t and z
plt.pcolormesh(t,z,data2plot, cmap = 'BuGn', edgecolors = 'None')

plt.colorbar(orientation="horizontal",fraction=0.07,pad=0.12)
plt.axis([0, 12, 250., 0.])
plt.title(vrn, fontsize=20, y=1.1)
plt.xlabel('month', fontsize=16)
plt.ylabel('depth [m]', fontsize=16)
labels=['J','F','M','A','M','J','J','A','S','O','N','D']
Xl    = np.arange(0.,12.)+0.5 # Major tick position
xl    = np.arange(0.,12) # Minor tick position
plt.xticks(Xl,labels)
ax2.set_xticks(xl, minor=True)
plt.tick_params(axis='x',which='major',length=0)
plt.tick_params(which='minor',length=3)
plt.xticks(Xl, labels)

theOutputFile= "prova.png"
fig.savefig(theOutputFile)
