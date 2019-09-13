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

z=np.arange(0,100.) + 0.5
Ed=np.zeros(z.shape[0])
Es=np.zeros(z.shape[0])
Eu=np.zeros(z.shape[0])

Ed[:]=Ed0*np.exp(-cd*rmud*z)

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

Es[:]=cp*np.exp(-kp*z[:]) + x * Ed[:]
Eu[:]=cp*rp*np.exp(-kp*z[:]) + y * Ed[:]
