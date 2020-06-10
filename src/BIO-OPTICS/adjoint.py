import numpy as np
import netCDF4 as NC4
import matplotlib.pyplot as plt

wl_mod = [250, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 
          650, 675, 700, 725, 775, 850, 950, 1050, 1150, 1250, 1350, 1450, 1550, 
          1650, 1750, 1900, 2200, 2900, 3700 ]
wl_idx = 6

ji=260
jj=200

dir_in='/gpfs/scratch/userexternal/plazzari/REA_16_CDOM/TEST_02/wrkdir/POSTPROC/output/AVE_FREQ_2/SURF_IRR/'

ed_file= dir_in + 'ave.20170416-00:00:00.Ed.nc'
ncin=NC4.Dataset(ed_file,"r")
Ed0 = ncin.variables['Ed'][wl_idx,jj,ji]
ncin.close()

es_file= dir_in + 'ave.20170416-00:00:00.Es.nc'
ncin=NC4.Dataset(es_file,"r")
Es0 = ncin.variables['Es'][wl_idx,jj,ji]
ncin.close()

eu_file= dir_in + 'ave.20170416-00:00:00.Eu.nc'
ncin=NC4.Dataset(eu_file,"r")
Eu0 = ncin.variables['Eu'][wl_idx,jj,ji]
ncin.close()


### Calculations at 450 un diatoms

#Ed0=1.0

#Es0=1.0

Rrs=0.004 #0.07

Q = 4.0

Eu_sat=Rrs*(Ed0+Es0)*Q

chl_v=np.arange(0.8,4.0,0.05)

Dl_Dcd = np.zeros(chl_v.shape[0])
cd_v   = np.zeros(chl_v.shape[0])
err_v  = np.zeros(chl_v.shape[0])

for chl_idx, chl_p in enumerate(chl_v):

    chl = 1.0
    
    #450
    aw   =  0.0085
    abio =  0.0379 * chl
    abio_p =  0.0379 * chl_p
    
    bw   = 0.0043
    bbio = 0.2738 * chl 
    bbio_p = 0.2738 * chl_p 
    
    #400
    #aw   =  0.0047
    #abio =  0.034 * chl
    
    #bw   = 0.0069
    #bbio = 0.2992 * chl 
    
    #550
    #aw   =  0.3512
    #abio =  0.0104 * chl
    
    #bw   = 0.0009
    #bbio = 0.2154 * chl 
    
    
    a = aw + abio
    a_p = aw + abio_p
    
    b = bw + bbio
    b_p = bw + bbio_p
    
    vd   = 1.0
    rmud = 1.0/vd
    
    cd = (a_p+b_p)*rmud
    cd_v[chl_idx] = cd
    
    print(cd)
    
    bbw = 0.5
    bbc = 0.002
    
    bb=bbw * bw + bbc * bbio
    
    rs=1.5
    vs =0.83
    
    ru=3.0
    vu =0.4
    
    rd=1.0
    
    gdepw =[0., 3.00062532507582, 6.23463155969512, 9.71806351694977, 13.4680658263969, 17.5029578589165, 21.8423136902711, 26.5070474319218,
     31.5195042790438, 36.9035576466558, 42.6847127863584, 48.8902173007809, 55.5491789959051, 62.6926915377699, 70.3539684071875, 78.5684856720181, 87.374134127298,
     96.8113813818782, 106.923444501721, 117.756473850692, 129.359748802613, 141.785886031008, 155.091061115687, 169.335244239759, 184.582450782764, 200.901007649911] 
    
    gdepw = np.asarray(gdepw)
    
    z=np.arange(0,100.)#gdepw
    z=gdepw
    dz=np.arange(len(z)) 
    
    Nz = len(z)
    Ed=np.zeros(z.shape[0])
    Es=np.zeros(z.shape[0])
    Eu=np.zeros(z.shape[0])
    
    Ed_ave=np.zeros(z.shape[0])
    Es_ave=np.zeros(z.shape[0])
    Eu_ave=np.zeros(z.shape[0])
    
    Ed[:]=Ed0*np.exp(-cd*z)
    depth = z
    for ii in range(Nz-1):
        dz     = depth[ii+1] - depth[ii]
        Ed_ave[ii] = ( Ed[ii+1] - Ed[ii] ) / (-cd*dz)
    
    Cs = (a+rs*bb)/vs
    Bu = ru*bb/vu
    Fd = (b-rd*bb)*rmud
    Cu = (a+ru*bb)/vu
    Bs = (rs*bb)/vs
    Bd = rd*bb*rmud
    
    den1 = (cd-Cs)*(cd+Cu) + Bs*Bu
    print(str(1.0/den1))
    
    x    = 1.0/den1 *( -Fd * (cd+Cu) - Bu*Bd  )
    y    = 1.0/den1 *( -Fd*Bs    + (cd-Cs)*Bd )
    
    D    = 0.5*(Cs+Cu + np.sqrt( (Cs+Cu)*(Cs+Cu) - 4.0*Bs*Bu ) )
    
    rp   = Bs/D
    kp   = -(Cu-D)
    cp   = Es0 -x * Ed0
    Rp   = (-Fd*rp-Bd)/(cd-Cs+D)
    
    rm   = Bu/D
    km   = -Cs+D 
    Rm   = (Fd+Bd*rm)/(cd+Cu-D)
    
    print "kp = " + str(kp)
    print "km = " + str(km)
    
    Eu[:] = cp*rp*np.exp(-kp*z[:]) + y * Ed[:]
    Es[:] = cp*np.exp(-kp*z[:])    + x * Ed[:]
    
    lambda_u0      = 2.0*(Eu[0] - Eu_sat)
    err_v[chl_idx] = (Eu[0] - Eu_sat) * (Eu[0] - Eu_sat)
    
    lambda_u = lambda_u0*np.exp(-km*z[:])
    lambda_s = lambda_u0*(-rp)*np.exp(-km*z[:]) 
    lambda_d = lambda_u0*(+Rp)*np.exp(-km*z[:]) 
    
    Dl_Dcd[chl_idx] = 0.0
    
    for ii in range(Nz-1):
        dz      = depth[ii+1] - depth[ii]
        Dl_Dcd[chl_idx] += - 0.5* (lambda_d[ii] * Ed[ii]/vd + lambda_s[ii] * Es[ii]/vs + lambda_u[ii] * Eu[ii]/vu + 
                        lambda_d[ii+1] * Ed[ii+1]/vd + lambda_s[ii+1] * Es[ii+1]/vs + lambda_u[ii+1] * Eu[ii+1]/vu) *dz
    
    
    print 'Eu[0]  = ' + str(Eu[0])
    print 'Eu_sat = ' + str(Eu_sat)
    print 'Dl_Dcd = ' + str(Dl_Dcd[chl_idx])

fig,axs = plt.subplots(2,1, gridspec_kw = {'wspace':0.05, 'hspace':0.05})

X=np.zeros((chl_v.shape[0],chl_v.shape[0]))
Y=np.zeros((chl_v.shape[0],chl_v.shape[0]))
U=np.zeros((chl_v.shape[0],chl_v.shape[0]))
V=np.zeros((chl_v.shape[0],chl_v.shape[0]))
base_0=np.zeros((chl_v.shape[0]))

for idy, cy in enumerate(cd_v):
    for idx, cx in enumerate(cd_v):
        X[idy,idx] = cx
        Y[idy,idx] = cy
for idx, cy in enumerate(cd_v):
    U[idx,idx] = Dl_Dcd[idx]/5.0

axs[0].plot(cd_v, err_v)
axs[0].fill_between(cd_v, base_0,err_v,facecolor='red')
axs[0].set_ylabel('Error')
axs[0].set_xticks([])
axs[1].quiver(X, Y, U, V )
axs[1].set_xlabel('cd')
axs[1].set_ylabel('')
axs[1].set_yticks([])
#Q = axs[1].quiver(X, Y, U, V, units='width')

fileout='lagrangian.png'
plt.savefig(fileout, format='png',dpi=150)
plt.close(fig)    # close the figure window
fig,axs = plt.subplots(1,3, gridspec_kw = {'wspace':0.09, 'hspace':0.05})

axs[0].plot(lambda_d,-z)
axs[0].set_ylabel('depth [m]')
axs[0].set_xlabel('$\lambda_d$')
axs[1].plot(lambda_s,-z)
#axs[1].set_ylabel('depth [m]')
axs[1].set_xlabel('$\lambda_s$')
axs[1].set_yticks([])
axs[2].plot(lambda_u,-z)
#axs[2].set_ylabel('depth [m]')
axs[2].set_xlabel('$\lambda_u$')
axs[2].set_yticks([])

fileout='lagrange_multipliers.png'
plt.savefig(fileout, format='png',dpi=150)
plt.close(fig)    # close the figure window


