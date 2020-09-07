#!/bin/env python

from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np

# Take Ed statistiscs to cluster your simulation improvements

wl_ls = ['380', '412', '490']

NAME_aw    = ['aw_LIT', 'aw_Mason', 'aw_Mason_TS']
RMSE_aw    = np.array([[0.33, 0.42, 0.22],[0.46, 0.50, 0.23],[0.43, 0.46, 0.23]])
BIAS_aw    = np.array([[0.30, 0.37, 0.17],[0.42, 0.45, 0.18],[0.40, 0.41, 0.17]])
CORR_aw    = np.array([[0.65, 0.68, 0.89],[0.40, 0.58, 0.89],[0.45, 0.63, 0.89]])

NAME_aNAP  = ['aNAP_Case1_CHL', 'aNAP_Case1_BBP', 'aNAP_Babin_CHL', 'aNAP_Babin_BBP', 'aNAP_Babin_CORR']
RMSE_aNAP  = np.array([[0.32, 0.35, 0.20],[0.28, 0.31, 0.20],[0.26, 0.30, 0.20],[0.14, 0.20, 0.19],[0.15, 0.19, 0.17]])
BIAS_aNAP  = np.array([[0.27, 0.29, 0.15],[0.24, 0.26, 0.14],[0.20, 0.23, 0.15],[0.11, 0.14, 0.13],[0.12, 0.14, 0.11]])
CORR_aNAP  = np.array([[0.63, 0.74, 0.90],[0.70, 0.77, 0.91],[0.70, 0.77, 0.90],[0.85, 0.85, 0.90],[0.85, 0.85, 0.91]])


NAME_aCDOM = ['aCDOM_Case1_CHL', 'aCDOM_Case1_fDOM', 'aCDOM_Kbio_Morel', 'aCDOM_Kbio_Mason', 'aCDOM_Kbio_Mason_CORR']
RMSE_aCDOM = np.array([[0.11, 0.13, 0.14],[0.10, 0.12, 0.14],[0.07, 0.09, 0.13],[0.06, 0.09, 0.12],[0.06, 0.09, 0.12]])
BIAS_aCDOM = np.array([[0.07, 0.07, 0.08],[0.06, 0.06, 0.08],[0.03, 0.03, 0.06],[0.03, 0.02, 0.06],[0.02, 0.01, 0.05]])
CORR_aCDOM = np.array([[0.88, 0.90, 0.93],[0.91, 0.91, 0.93],[0.93, 0.93, 0.93],[0.93, 0.93, 0.93],[0.93, 0.93, 0.93]])

# Wrap up with adding PFT and BBP contributions to close the discussion

NAME_remaining = [['aCDOM'],['aPFT'],['bp']] 
RMSE_remaining = np.array([[0.06, 0.09, 0.12],[0.05, 0.08, 0.10],[0.05, 0.08, 0.09]])
BIAS_remaining = np.array([[0.02, 0.01, 0.05],[0.01, 0.0,  0.02],[0.01, -0.01, 0.02]])
CORR_remaining = np.array([[0.93, 0.93, 0.93],[0.94, 0.95, 0.95],[0.94, 0.95, 0.95]])

# Relative contributions of each group of IOPs is demonstrated with Kd scatter plot skill.

NAME_relative  = ['REF', 'aNAP', 'aCDOM', 'aPFT', 'bp']
RMSE_relative  = np.array([[0.03, 0.02, 0.01],[0.02, 0.01, 0.01],[0.03, 0.02, 0.01],[0.02, 0.02, 0.02],[0.03, 0.02, 0.01]])
BIAS_relative  = np.array([[0.02, 0.01, 0.],[0., 0., 0.],[-0.03, -0.02, 0.],[0.01, 0.,-0.01],[0.02, 0.01, 0.]])
CORR_relative  = np.array([[0.84, 0.89, 0.90],[0.85, 0.90, 0.90],[0.75, 0.85, 0.89],[0.83, 0.83, 0.80],[0.84, 0.89, 0.91]])
SLOPE_relative = np.array([[1.15, 1.08, 0.88],[1.14, 1.08,0.87],[0.53, 0.70, 0.73],[0.90, 0.56, 0.24],[1.15, 1.06, 0.83]])

# Now plot it
fig1, ax1 = plt.subplots(1,3)#, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
fig1.set_size_inches(18, 6)

fig2, ax2 = plt.subplots(1,3)#, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
fig2.set_size_inches(18, 6)

fig3, ax3 = plt.subplots(1,3)#, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
fig3.set_size_inches(18, 6)

fig4, ax4 = plt.subplots(1,3)#, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
fig4.set_size_inches(18, 6)

fig5, ax5 = plt.subplots(1,4)#, gridspec_kw = {'wspace':0.25, 'hspace':0.5})
fig5.set_size_inches(24, 6)


plot_barplot(ax1,  NAME_aw,    RMSE_aw,    BIAS_aw,    CORR_aw,    wl_ls, 30, None)
plot_barplot(ax2,  NAME_aNAP,  RMSE_aNAP,  BIAS_aNAP,  CORR_aNAP,  wl_ls, 30, None)
plot_barplot(ax3,  NAME_aCDOM, RMSE_aCDOM, BIAS_aCDOM, CORR_aCDOM, wl_ls, 30, None)
plot_barplot(ax4,  NAME_remaining, RMSE_remaining, BIAS_remaining, CORR_remaining, wl_ls, 30, None)
plot_barplot(ax5,  NAME_relative, RMSE_relative, BIAS_relative, CORR_relative, wl_ls, 30, 'slope')


fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()
fig4.tight_layout()
fig5.tight_layout()


fig1.savefig('aw.png', dpi=300)
fig2.savefig('aNAP.png', dpi=300)
fig3.savefig('aCDOM.png', dpi=300)
fig4.savefig('aPHY.png', dpi=300)
fig5.savefig('IOP_REL.png', dpi=300)


