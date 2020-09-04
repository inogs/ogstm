import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import string

def plot_matchup_scatter(scale, L, ax, color, index, units, titlestr, xpos, ypos, legendBool, basin):
	
	if scale == 'lin':
		x = L.Ref
		y = L.Model
		
	elif scale == 'log':
		x = np.log(L.Ref)
		y = np.log(L.Model)
	
	if legendBool == True:
		ax[index].scatter(x, y, marker='o', s=0.05, c=color, label=basin)
	else:
		ax[index].scatter(x, y, marker='o', s=0.05, c=color)

	if units == 'Ed':

		ax[index].set_xlabel('BGC-Argo float [$W \, m^{-2} \, nm^{-1}$]', fontsize=12)
		if index == 0:
			ax[index].set_ylabel('BIOPTIMOD [$W \, m^{-2} \, nm^{-1}$]', fontsize=12)

	if units == 'Kd':

		ax[index].set_xlabel('BGC-Argo float [$m^{-1}$]', fontsize=12)
		if index == 0:
			ax[index].set_ylabel('BIOPTIMOD [$m^{-1}$]', fontsize=12)


	count      = L.number()
	corr_coeff = L.correlation()
	bias       = L.bias()
	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	sigma      = L.RMSE()
	
	a          = intercept
	b          = slope
	
	x_max      = max(x.max(), y.max())*1.1          
	x_reg      = np.linspace(0., x_max, 50)

	ax[index].plot(x_reg, a + b*x_reg, color)
	ax[index].plot(x_reg,x_reg,'k--')
	
	if scale == 'log':
		ax[index].set_xscale('log')
		ax[index].set_yscale('log')
	
	textstr='$\mathrm{RMS}=%.2f$\n$\mathrm{Bias}=%.2f$\n$\mathrm{r}=%.2f$\n$\mathrm{Slope}=%.2f$\n$\mathrm{Y-int}=%.2f$\n$\mathrm{N}=%.2i$'%(sigma, bias, corr_coeff,b,a,count)
	
	if legendBool == True:
		ax[index].legend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=2, fancybox=True, shadow=True)
	
	ax[index].text(xpos, ypos, textstr, transform=ax[index].transAxes, fontsize=12, color = color, verticalalignment='top',bbox=dict(facecolor='white', alpha = 0.5, edgecolor=color))
	ax[index].set_title(titlestr, fontsize=24)
	ax[index].tick_params(axis='both', which='major', labelsize=14)
	ax[index].set_aspect('equal', adjustable='box')
	ax[index].set_xlim([0., x_max])
	ax[index].set_ylim([0., x_max])
		
	return ax[index]

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


def plot_pcolor(fig, ax, MODEL_mean, FLOAT_mean, BIAS, RMSE, titlestr, strname, basin_list_abbrev, months_str):

	x_ticks_pos = np.arange(len(basin_list_abbrev)) + 0.5
	y_ticks_pos = np.arange(len(months_str)) + 0.5

	vmax_AB = max(np.nanmax(MODEL_mean), np.nanmax(FLOAT_mean))
	vmax_C  = max(abs(np.nanmin(BIAS/FLOAT_mean)), abs(np.nanmax(BIAS/FLOAT_mean)))
	vmax_D  = np.nanmax(RMSE/FLOAT_mean)

	# Plot the model values
	cmap = plt.get_cmap('BuPu', 8) 
	c0 = ax[0,0].pcolormesh(MODEL_mean.T, cmap=cmap, vmin=0., vmax=vmax_AB)
	cmap.set_bad('lightgrey',1.)
	ax[0,0].set_title(r'$' + titlestr + '_{\lambda=' + strname +'}$ M (MEAN) $[m^{-1}]$ ')
	ax[0,0].set_xticks(x_ticks_pos)
	ax[0,0].set_yticks(y_ticks_pos)
	ax[0,0].set_xticklabels(tuple(basin_list_abbrev), rotation=45, fontsize=8, ha='center')
	ax[0,0].set_yticklabels(tuple(months_str), rotation=0, fontsize=8, va='center')
	ax[0,0].text(0.0, 1.03, string.ascii_lowercase[0]+ ')', transform=ax[0,0].transAxes, 
				size=14, weight='bold')
	cbar = fig.colorbar(c0, ax=ax[0,0])
	cbar.locator   = matplotlib.ticker.LinearLocator(numticks=9)
	cbar.formatter = matplotlib.ticker.FormatStrFormatter("%.2f")
	cbar.update_ticks()

	
	# Plot float values
	cmap = plt.get_cmap('BuPu', 8)   
	c1 = ax[0,1].pcolormesh(FLOAT_mean.T, cmap=cmap, vmin=0., vmax=vmax_AB)
	cmap.set_bad('lightgrey',1.)
	ax[0,1].set_title(r'$' + titlestr + '_{\lambda='+ strname +'}$ O (MEAN) $[m^{-1}]$')
	ax[0,1].set_xticks(x_ticks_pos)
	ax[0,1].set_yticks(y_ticks_pos)
	ax[0,1].set_xticklabels(tuple(basin_list_abbrev), rotation=45, fontsize=8, ha='center')
	ax[0,1].set_yticklabels(tuple(months_str), rotation=0, fontsize=8, va='center')
	ax[0,1].text(0.0, 1.03, string.ascii_lowercase[1] + ')', transform=ax[0,1].transAxes, 
				size=14, weight='bold')
	cbar1 = fig.colorbar(c1, ax=ax[0,1])
	cbar1.locator = matplotlib.ticker.LinearLocator(numticks=9)
	cbar1.formatter = matplotlib.ticker.FormatStrFormatter("%.2f")
	cbar1.update_ticks()

	
	# Plot bias, 
	cmap = plt.get_cmap('bwr',8)
	cmap.set_bad('lightgrey',1.)
	c2 = ax[1,0].pcolormesh(BIAS.T/FLOAT_mean.T , cmap=cmap, vmin=-vmax_C, vmax=vmax_C)
	ax[1,0].set_title(r'$' + titlestr + '_{\lambda='+  strname +'}$ BIAS (normalized) $[-]$')
	ax[1,0].set_xticks(x_ticks_pos)
	ax[1,0].set_yticks(y_ticks_pos)
	ax[1,0].set_xticklabels(tuple(basin_list_abbrev), rotation=45, fontsize=8, ha='center')
	ax[1,0].set_yticklabels(tuple(months_str), rotation=0, fontsize=8, va='center')
	ax[1,0].text(0.0, 1.05, string.ascii_lowercase[2] + ')', transform=ax[1,0].transAxes, 
				size=14, weight='bold')
	cbar2 = fig.colorbar(c2, ax=ax[1,0])
	cbar2.locator = matplotlib.ticker.LinearLocator(numticks=9)
	cbar2.formatter = matplotlib.ticker.FormatStrFormatter("%.2f")
	cbar2.update_ticks()


	# Plot RMSE   
	cmap = plt.get_cmap('BuPu',8)
	cmap.set_bad('lightgrey',1.)
	c3 = ax[1,1].pcolormesh(RMSE.T/FLOAT_mean.T , cmap=cmap, vmin=0., vmax=vmax_D)
	ax[1,1].set_title(r'$' + titlestr + '_{\lambda='+  strname +'}$ RMSD (normalized) $[-]$')
	ax[1,1].set_xticks(x_ticks_pos)
	ax[1,1].set_yticks(y_ticks_pos)
	ax[1,1].set_xticklabels(tuple(basin_list_abbrev), rotation=45, fontsize=8, ha='center')
	ax[1,1].set_yticklabels(tuple(months_str), rotation=0, fontsize=8, va='center')
	ax[1,1].text(0.0, 1.05, string.ascii_lowercase[3]+ ')', transform=ax[1,1].transAxes, 
				size=14, weight='bold')
	cbar3 = fig.colorbar(c3, ax=ax[1,1])
	cbar3.locator = matplotlib.ticker.LinearLocator(numticks=9)
	cbar3.formatter = matplotlib.ticker.FormatStrFormatter("%.2f")
	cbar3.update_ticks()
	
	return cbar, cbar1, cbar2, cbar3

