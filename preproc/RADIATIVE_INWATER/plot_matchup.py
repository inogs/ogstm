import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

def plot_matchup(scale, L, ax, color, index, titlestr, xpos, ypos, legendBool, basin):
    
    if scale == 'lin':
        x = L.Ref
        y = L.Model
        
    elif scale == 'log':
        x = np.log(L.Ref)
        y = np.log(L.Model)
        
    '''Mask values in case of any NaNs'''
    mask = ~np.isnan(x) & ~np.isnan(y)
    
    if legendBool == True:
        ax[index].scatter(x, y, marker='o', s=0.05, c=color, label=basin)
    else:
        ax[index].scatter(x, y, marker='o', s=0.05, c=color)

    ax[index].set_xlabel('BGC-Argo float [$\mu W \, cm^{-2} \, nm^{-1}$]')
    if index == 0:
        ax[index].set_ylabel('BIOPTIMOD [$\mu W \, cm^{-2} \, nm^{-1}$]')

    count = L.number()
    corr_coeff = L.correlation()
    bias = L.bias()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask],y[mask]) 
    sigma = L.RMSE()
    
    a = intercept
    b = slope
    
    x_max = x.max() *  1.1
    x_reg = np.arange(0., x_max)
    
    ax[index].plot(x_reg,a+b*x_reg,color)
    ax[index].plot(x_reg,x_reg,'k--')
    
    if scale == 'log':
        ax[index].set_xscale('log')
        ax[index].set_yscale('log')
    
    textstr='$\mathrm{RMS}=%.2f$\n$\mathrm{Bias}=%.2f$\n$\mathrm{r}=%.2f$\n$\mathrm{Slope}=%.2f$\n$\mathrm{Y-int}=%.2f$\n$\mathrm{N}=%.2i$'%(sigma, bias, corr_coeff,b,a,count)
    
    if legendBool == True:
        ax[index].legend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=2, fancybox=True, shadow=True)
    
    ax[index].text(xpos, ypos, textstr, transform=ax[index].transAxes, fontsize=8, color = color, verticalalignment='top',bbox=dict(facecolor='white', alpha = 0.5, edgecolor=color))
    ax[index].set_title(titlestr)
    return ax

def save_stat(L):
    
    x = L.Ref
    y = L.Model
    
    '''Mask values in case of any NaNs'''
    mask = ~np.isnan(x) & ~np.isnan(y)
    
    count = L.number()
    corr_coeff = L.correlation()
    bias = L.bias()
    slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask],y[mask]) 
    sigma = L.RMSE()
    
    a = intercept
    b = slope

    return count, bias, sigma, r_value, b, a

