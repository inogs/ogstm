from basins import V2 as OGS
from instruments import optbio_float_2019
from instruments import var_conversions
from instruments.matchup_manager import Matchup_Manager
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
from profiler import *

def plot_profile(x,y, color, variable, title):
    plt.clf()
    fig = plt.figure()
    fig.set_size_inches(8,12)
    plt.scatter(x, y, c=color)
    if variable == 'CHL':
        plt.xlabel(r'[$mg \, m^{-3} $]', fontsize=24)
        
    else:
        plt.xlabel('[$\mu W \, cm^{-2} nm^{-1}$]', fontsize=24)
        plt.xlim([0.01, 100])
        plt.xscale('log')
    plt.ylabel('depth [m]', fontsize=24)
    plt.ylim([0, 200.])
    plt.gca().invert_yaxis()
    plt.gca().tick_params(axis='both', which='major', labelsize=20)
    plt.title(title, fontsize = 32)
    plt.savefig(variable + '.png', dpi=400, transparent=True)
    return

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

TI = TimeInterval('20120101', '20171231',"%Y%m%d")

variable='P_l'
varname=var_conversions.FLOAT_OPT_VARS_2019[variable]

Profilelist=optbio_float_2019.FloatSelector(varname,TI , OGS.med)

p = Profilelist[256]

PresCHL, CHLz,    Qc = p.read('CHL')
Pres380, Ed380,    Qc = p.read('IRR_380')
Pres412, Ed412,    Qc = p.read('IRR_412')
Pres490, Ed490,    Qc = p.read('IRR_490')

plot_profile(CHLz, PresCHL, 'darkgreen', 'CHL', 'Chlorophyll a')
plot_profile(Ed380, Pres380, 'indigo', 'Ed_380', '$Ed _{\lambda=380}$')
plot_profile(Ed412, Pres412, 'royalblue', 'Ed_412','$Ed _{\lambda=412}$')
plot_profile(Ed490, Pres490, 'navy', 'Ed_490','$Ed _{\lambda=490}$')