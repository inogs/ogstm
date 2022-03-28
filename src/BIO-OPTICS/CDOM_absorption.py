import numpy as np
import matplotlib.pyplot as plt


c_cdom   = 0.18 

s_cdom   = 0.017  #for the Adriatic s_cdom = 0.0192 ; Dutkiewicz: s_cdom = 0.021 

lambda_0 = 443.0  #from Babin et al., 2003 ; Dutkiewicz: lambda_0 = 450.0


lam = [  250.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 
         500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 
         700.0, 725.0, 775.0, 850.0, 950.0, 1050.0, 1150.0, 1250.0,
        1350.0, 1450.0, 1550.0,  1650.0, 1750.0, 1900.0, 2200.0, 2900.0, 3700.0 , 4000.0] 

acdom_bin = np.zeros(len(lam)-1)

lam_1nm = np.arange(250.0, 4000.0, 1.0)

a_cdom = c_cdom * np.exp(-s_cdom*(lam_1nm -lambda_0))

lam_bin = np.zeros(len(lam)-1)
for nl in range(len(lam)-1):
   acdom_bin[nl] = c_cdom * ( np.exp(-s_cdom*(lam[nl+1] -lambda_0)) -np.exp(-s_cdom*(lam[nl] -lambda_0)) )  / (-s_cdom *( lam[nl+1] - lam[nl] ) )   
   lam_bin[nl] = (lam[nl] + lam[nl+1])/2.
   
#acdom_int = np.interp(lam_bin, lam_1nm, a_cdom)

plt.plot(lam_1nm, a_cdom/12., 'm') 
plt.bar(lam_bin, acdom_bin/12., width=25, align='center', color='g')
plt.xlim(350,700) 
plt.ylim(0,2) 
plt.xlabel(r'$\lambda$ [nm]') 
plt.ylabel('$a^{CDOM} [m^{2}.mgC^{-1}]$') 
#plt.show() 

for i in range(len(lam)-1):
    mystring = 'acdom(' + str (i+1) + ')= ' + str(acdom_bin[i]/12)  + 'D0'
    print(mystring)
