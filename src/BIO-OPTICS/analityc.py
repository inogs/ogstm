import numpy as np


### Calculations at 450 un diatoms

Ed0=1.0

Es0=1.0

chl = 1.0

#450
#aw   =  0.0085
#abio =  0.0379 * chl

#bw   = 0.0043
#bbio = 0.2738 * chl 

#400
aw   =  0.0047
abio =  0.034 * chl

bw   = 0.0069
bbio = 0.2992 * chl 

#550
#aw   =  0.3512
#abio =  0.0104 * chl

#bw   = 0.0009
#bbio = 0.2154 * chl 


a = aw + abio

b = bw + bbio

rmud=1.0

cd=(a+b)*rmud

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

Ed[:]=Ed0*np.exp(-cd*rmud*z)
depth = z
for ii in range(Nz-1):
    dz     = depth[ii+1] - depth[ii]
    Ed_ave[ii] = ( Ed[ii+1] - Ed[ii] ) / (-cd*rmud*dz)

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

print kp

Es[:] = cp*np.exp(-kp*z[:])    + x * Ed[:]
Eu[:] = cp*rp*np.exp(-kp*z[:]) + y * Ed[:]

for ii in range(Nz-1):
    dz     = depth[ii+1] - depth[ii]
    Es_ave[ii] = cp * (np.exp(-kp*z[ii+1])      - np.exp(-kp*z[ii]))/(-kp*dz) + x * Ed_ave[ii]
    Eu_ave[ii] = cp * rp * (np.exp(-kp*z[ii+1]) - np.exp(-kp*z[ii]))/(-kp*dz) + y * Ed_ave[ii]


