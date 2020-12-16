import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import string

#from matchup.statistics import *


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

	if units == 'Rrs':

		ax[index].set_xlabel('BGC-Argo float [$sr^{-1}$]', fontsize=12)
		if index == 0:
			ax[index].set_ylabel('BIOPTIMOD [$sr^{-1}$]', fontsize=12)


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

	if units == 'Ed':
	
		textstr='$\mathrm{RMSE}=%.2f$\n$\mathrm{Bias}=%.2f$\n$\mathrm{r}=%.2f$\n$\mathrm{Slope}=%.2f$\n$\mathrm{Y-int}=%.2f$\n$\mathrm{N}=%.2i$'%(sigma, bias, corr_coeff,b,a,count)
	
	if units == 'Kd':
		textstr='$\mathrm{RMSE}=%.3f$\n$\mathrm{Bias}=%.3f$\n$\mathrm{r}=%.2f$\n$\mathrm{Slope}=%.2f$\n$\mathrm{Y-int}=%.2f$\n$\mathrm{N}=%.2i$'%(sigma, bias, corr_coeff,b,a,count)

	if units == 'Rrs':
		textstr='$\mathrm{RMSE}=%.3f$\n$\mathrm{Bias}=%.3f$\n$\mathrm{r}=%.2f$\n$\mathrm{Slope}=%.2f$\n$\mathrm{Y-int}=%.2f$\n$\mathrm{N}=%.2i$'%(sigma, bias, corr_coeff,b,a,count)


	if legendBool == True:
		ax[index].legend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=2, fancybox=True, shadow=True)
	
	ax[index].text(xpos, ypos, textstr, transform=ax[index].transAxes, fontsize=12, color = color, verticalalignment='top',bbox=dict(facecolor='white', alpha = 0.5, edgecolor=color))
	ax[index].set_title(titlestr, fontsize=24)
	ax[index].tick_params(axis='both', which='major', labelsize=14)
	ax[index].set_aspect('equal', adjustable='box')
	ax[index].set_xlim([0., x_max])
	ax[index].set_ylim([0., x_max])
		
	return ax[index]

def save_stat(Model, Data, filename):

	L = matchup(Model, Data)
	
	x = L.Ref
	y = L.Model
	
	'''Mask values in case of any NaNs'''
	mask       = ~np.isnan(x) & ~np.isnan(y)
	
	count      = L.number()
	corr_coeff = L.correlation()
	bias       = L.bias()
	maxdiff    = L.maxdiff()
	slope, intercept, r_value, p_value, std_err = stats.linregress(x[mask],y[mask]) 
	sigma      = L.RMSE()
	
	a = intercept
	b = slope

	if filename is not None:
		file = open(filename, 'w')
		file.write('%f\n%f\n%f\n%f\n%f\n%f\n%f\n' %(count, bias, sigma, r_value, b, a, maxdiff) )
		file.close()

	return count, bias, sigma, r_value, b, a, maxdiff

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
	
	return 

def plot_barplot(ax, NAME, RMSE, BIAS, CORR, wl_ls, rot, slope=None):

	X      = np.arange(len(NAME))
	
	ax[0].bar(X - 0.25, RMSE[:,0], color='indigo',   width=0.25)
	ax[0].bar(X + 0.00, RMSE[:,1], color='darkcyan', width=0.25)
	ax[0].bar(X + 0.25, RMSE[:,2], color='navy',     width=0.25)
	ax[0].set_xticks(X)
	ax[0].set_xticklabels(NAME, rotation=rot, ha='right')
	ax[0].set_ylabel('RMSE [$W \, m^{-2} \, nm^{-1}$]')
	ax[0].set_ylim(bottom=0.)


	ax[1].bar(X - 0.25, BIAS[:,0], color='indigo',  label=wl_ls[0], width=0.25)
	ax[1].bar(X + 0.00, BIAS[:,1], color='darkcyan',label=wl_ls[1], width=0.25)
	ax[1].bar(X + 0.25, BIAS[:,2], color='navy',    label=wl_ls[2], width=0.25)
	ax[1].set_xticks(X)
	ax[1].set_xticklabels(NAME, rotation=rot, ha='right')
	ax[1].set_ylabel('bias [$W \, m^{-2} \, nm^{-1}$]')
	#ax[1].set_ylim(bottom=0.)
	ax[1].set_ylim([0., 0.5])

	ax[2].bar(X - 0.25, CORR[:,0], color='indigo',   label=wl_ls[0], width=0.25)
	ax[2].bar(X + 0.00, CORR[:,1], color='darkcyan', label=wl_ls[1], width=0.25)
	ax[2].bar(X + 0.25, CORR[:,2], color='navy',     label=wl_ls[2], width=0.25)
	ax[2].set_xticks(X)
	ax[2].set_xticklabels(NAME, rotation=rot, ha='right')
	ax[2].set_ylabel('r')
	ax[2].set_ylim([0., 1.])

	ilegend=1

	ax[ilegend].legend(loc='upper right', ncol=3)#, fancybox=True, shadow=True)

	return 

def plot_barplot2(ax, NAME, RMSE, BIAS, wl_ls, rot, slope=None):

	X      = np.arange(len(NAME))
	
	ax[0].bar(X - 0.25, RMSE[:,0], color='indigo',   width=0.25)
	ax[0].bar(X + 0.00, RMSE[:,1], color='darkcyan', width=0.25)
	ax[0].bar(X + 0.25, RMSE[:,2], color='navy',     width=0.25)
	ax[0].set_xticks(X)
	ax[0].set_xticklabels(NAME, rotation=rot, ha='right')
	ax[0].set_ylabel('RMSE [$W \, m^{-2} \, nm^{-1}$]')
	ax[0].set_ylim(bottom=0.)


	ax[1].bar(X - 0.25, BIAS[:,0], color='indigo',  label=wl_ls[0], width=0.25)
	ax[1].bar(X + 0.00, BIAS[:,1], color='darkcyan',label=wl_ls[1], width=0.25)
	ax[1].bar(X + 0.25, BIAS[:,2], color='navy',    label=wl_ls[2], width=0.25)
	ax[1].set_xticks(X)
	ax[1].set_xticklabels(NAME, rotation=rot, ha='right')
	ax[1].set_ylabel('bias [$W \, m^{-2} \, nm^{-1}$]')
	ax[1].set_ylim(bottom=0.)

	ilegend=1

	ax[ilegend].legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.0), framealpha=0.75, fancybox=True, shadow=False)

	return 


def plot_lineplot(ax, NAME, RMSE, BIAS, CORR, wl_ls, rot, slope=None):

	X      = np.arange(len(NAME))
	
	ax[0].plot(X, RMSE[:,0], color='indigo',   marker='o')
	ax[0].plot(X, RMSE[:,1], color='darkcyan', marker='o')
	ax[0].plot(X, RMSE[:,2], color='navy',     marker='o')
	ax[0].set_xticks(X)
	ax[0].set_xticklabels(NAME, rotation=rot, ha='right')
	ax[0].set_ylabel('RMSE [$W \, m^{-2} \, nm^{-1}$]')
	ax[0].set_ylim(bottom=0.)

	ax[1].plot(X, BIAS[:,0], color='indigo',   marker='o')
	ax[1].plot(X, BIAS[:,1], color='darkcyan', marker='o')
	ax[1].plot(X, BIAS[:,2], color='navy',     marker='o')
	ax[1].set_xticks(X)
	ax[1].set_xticklabels(NAME, rotation=rot, ha='right')
	ax[1].set_ylabel('BIAS [$W \, m^{-2} \, nm^{-1}$]')
	ax[1].set_ylim(bottom=0.)


	ilegend = 2

	ax[2].plot(X, CORR[:,0], color='indigo',   label=wl_ls[0], marker='o')
	ax[2].plot(X, CORR[:,1], color='darkcyan', label=wl_ls[1], marker='o')
	ax[2].plot(X, CORR[:,2], color='navy',     label=wl_ls[2], marker='o')
	ax[2].set_xticks(X)
	ax[2].set_xticklabels(NAME, rotation=rot, ha='right')
	ax[2].set_ylabel('Correlation')
	ax[2].set_ylim([0., 1.])


	if slope is not None:
		ilegend = 3
		ax[3].plor(X, CORR[:,0], color='indigo',   label=wl_ls[0], marker='o')
		ax[3].plor(X, CORR[:,1], color='darkcyan', label=wl_ls[1], marker='o')
		ax[3].plor(X, CORR[:,2], color='navy',     label=wl_ls[2], marker='o')
		ax[3].set_xticks(X)
		ax[3].set_xticklabels(NAME, rotation=rot, ha='right')
		ax[3].set_ylabel('Slope')

	ax[ilegend].legend(loc='center left', bbox_to_anchor=(1, 0.5))


	return 

