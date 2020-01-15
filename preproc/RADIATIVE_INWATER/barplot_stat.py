
import numpy as np
import matplotlib.pyplot as plt



runs = ['RUN01/', 'RUN02a/', 'RUN02b/', 'RUN03/', 'RUN04a/', 'RUN04b/', 'RUN04c/', 'RUN04d/']#, 'RUN04e/']
run_name = ['pure_water', 'CDOM min', 'CDOM max','NAP', 'PFT1', 'PFT2', 'PFT3', 'PFT4'] #, 'PFT5']

nWl   = 3
nRuns = len(runs)
nStat = 6
#wl, run, stat
OUTFILE_aux = np.load('OUTFILE.npy')
OUTFILE = OUTFILE_aux[:,:-1,1:]

width = 0.2

ind = np.arange(nRuns)

colors = ['indigo', 'darkcyan', 'navy']
wl_list = ['Ed 380', 'Ed 412', 'Ed 490']

statlist = ['N', 'bias', 'RMSE', 'r', 'slope', 'int']

#[$\mu W \, cm^{-2} \, nm^{-1}$]
fig, axs = plt.subplots(nStat-1,1)


for iStat in range(nStat-1):
    for iWl in range(nWl):
        axs[iStat].bar(ind + width*iWl, OUTFILE[iWl, :, iStat], width, color = colors[iWl], label=wl_list[iWl]) 

    axs[iStat].set_ylabel(statlist[iStat])
    axs[iStat].set_xticks([])
    axs[iStat].set_xticklabels([])
    
    if iStat == nStat-2:
        axs[iStat].set_xticks(ind + width*1.5)
        axs[iStat].set_xticklabels(tuple(run_name[1:]), rotation=45, fontsize=8)
    
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), fancybox=True, framealpha=0.5, ncol=3)


