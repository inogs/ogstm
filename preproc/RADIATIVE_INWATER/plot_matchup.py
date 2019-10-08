import matplotlib.pyplot as plt
from matchup.statistics import matchup
from scipy import stats

#x = Ed380_float
#y = Ed380_model

#M380 = matchup(x,y)

def plot_matchup(L):

    fig,ax=plt.subplots()
    ax.scatter(L.Ref, L.Model, marker='o', s=6**2, c='k')

    ax.set_ylabel('BIOPTIMOD [$\mu W \, cm^{-2} \, nm^{-1}$]')
    ax.set_xlabel('BGC-Argo float [$\mu W \, cm^{-2} \, nm^{-1}$]')

    
    count = L.number()
    corr_coeff = L.correlation()
    bias = L.bias()
    slope, intercept, r_value, p_value, std_err = stats.linregress(L.Ref,L.Model)
    sigma = L.RMSE()
    
    a = intercept
    b = slope
    
    x_max = L.Ref.max() *  1.1
    x_reg = np.arange(0., x_max)
    
    ax.plot(x_reg,a+b*x_reg,'k')
    ax.plot(x_reg,x_reg,'r')
    
    textstr='$\mathrm{RMS}=%.2f$\n$\mathrm{Bias}=%.2f$\n$\mathrm{r}=%.2f$\n$\mathrm{Slope}=%.2f$\n$\mathrm{Y-int}=%.2f$\n$\mathrm{N}=%.2f$'%(sigma, bias, corr_coeff,b,a,count)
    
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=8, verticalalignment='top',bbox=dict(facecolor='white', edgecolor='black'))
    return fig, ax