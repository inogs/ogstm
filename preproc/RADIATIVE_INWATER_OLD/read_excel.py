import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


INPUTFILE  = 'INPUT/Phyto_spec_Organelli2017_CSV_merged.csv'  # Original csv file from spectrophotometric measurements
OASIM_abs  = 'INPUT/absorption.csv'                      # A csv file used in the OASIM model - absorption
OASIM_scat = 'INPUT/scattering.csv'                      # A csv file used in the OASIM model - scattering

# Use the first two rows as headers, the first column for indexing
dfM        =  pd.read_csv(INPUTFILE, header=[0,1], index_col=0,  delimiter='\t')  

dfA        =  pd.read_csv(OASIM_abs,               index_col=1,  delimiter='\t') 
dfB        =  pd.read_csv(OASIM_scat,              index_col=1,  delimiter='\t') 

del dfA['Unnamed: 0']
del dfB['Unnamed: 0']


# PFT groups
PFT_spec  = dfM.columns.levels[0] 
PFT_OASIM = dfA.columns

'''
 1. First plot the spectrophotometric curves
'''
nrow = 4
ncol = 2
fig, axes = plt.subplots(nrow, ncol,  gridspec_kw = {'wspace':0.15, 'hspace':0.35})
fig.set_size_inches(11.69,8.27)

count = 0   #plot counter
for r in range(nrow):
	for c in range(ncol):
		dfM[PFT_spec[count]].plot(ax=axes[r,c])

		axes[r,c].set_title(PFT_spec[count])
		axes[r,c].legend(loc='upper right', framealpha=0.5)

		if r == 3:
			axes[r,c].set_xlabel('$\lambda [nm]$')
		if c == 0:
			axes[r,c].set_ylabel('$a^{}_{\phi} [m^{2} \, mg^{-1}(Chl)]$')
		count+=1

fig.savefig('Spec.png', dpi=300)


'''
 2. Plot the OASIM curves - both absorption and scattering ones
'''
# Now you want to bin your data set on wavebands of the OASIM model
# Specifiy your desired dz step size
wl_bin = [337.5 ,362.5, 387.5, 412.5, 437.5, 462.5, 487.5, 512.5, 537.5, 562.5, 587.5, 612.5, 637.5, 662.5, 687.5 ,712.5, 737.5, 762.5 ] 
wl_out = [350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 700.0, 725.0, 750.0, 775.0] 		

# Rebin dataframe
dfM_1 = dfM.groupby(pd.cut(dfM.index, wl_bin, labels=False), as_index=False).mean()

# Refill 'depth' column
dfM_1.index = wl_out[:-1]

# Plot OASIM absorption
nrow = 3
ncol = 2

fig, axes = plt.subplots(nrow, ncol,  gridspec_kw = {'wspace':0.15, 'hspace':0.35})
fig.set_size_inches(11.69,8.27)

colors = cm.rainbow(np.linspace(0, 1, len(PFT_OASIM))) 
#plot counter
count = 0
for r in range(nrow):
	for c in range(ncol):

		if count == 5:
			axes[r,c].set_visible(False)
		else:

			dfA[PFT_OASIM[count]].plot(ax=axes[r,c],linestyle=' ', marker='o', markerfacecolor=colors[count])

			axes[r,c].set_title(PFT_OASIM[count])
			if r == 2:
				axes[r,c].set_xlabel('$\lambda [nm]$')
			if c == 0:
				axes[r,c].set_ylabel('$a^{}_{\phi} [m^{2} \, mg^{-1}(Chl)]$')

		count+=1

fig.savefig('OASIM_abs.png', dpi=300)


# Plot OASIM scattering
nrow = 3
ncol = 2

fig, axes = plt.subplots(nrow, ncol,  gridspec_kw = {'wspace':0.15, 'hspace':0.35})
fig.set_size_inches(11.69,8.27)

colors = cm.rainbow(np.linspace(0, 1, len(PFT_OASIM))) 
#plot counter
count = 0
for r in range(nrow):
	for c in range(ncol):

		if count == 5:
			axes[r,c].set_visible(False)
		else:
			dfB[PFT_OASIM[count]].plot(ax=axes[r,c],linestyle=' ', marker='o', markerfacecolor=colors[count])


			axes[r,c].set_title(PFT_OASIM[count])
			if r == 2:
				axes[r,c].set_xlabel('$\lambda [nm]$')
			if c == 0:
				axes[r,c].set_ylabel('$b^{}_{p} [m^{2} \, mg^{-1}(Chl)]$')
		count+=1

fig.savefig('OASIM_scat.png', dpi=300)


'''
 4. Plot the spectrophotometric curves binned on OASIM wavebands (25 nm)
'''
nrow = 4
ncol = 2

fig, axes = plt.subplots(nrow, ncol, gridspec_kw = {'wspace':0.15, 'hspace':0.45})
fig.set_size_inches(11.69,8.27)
#plot counter
count = 0

for r in range(nrow):
	for c in range(ncol):
		dfM_1[PFT_spec[count]].plot(ax=axes[r,c], style='o', xlim=[325,775])
		axes[r,c].set_title(PFT_spec[count])
		if r == 3:
			axes[r,c].set_xlabel('$\lambda [nm]$')
		if c == 0:
			axes[r,c].set_ylabel('$a^{}_{\phi} [m^{2} \, mg^{-1}(Chl)]$')
		count+=1

fig.savefig('Spec_25nm.png', dpi=300)