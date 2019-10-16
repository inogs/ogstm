#eterzic 07.10.2019
import numpy as np
# Create functions for chl-specific IOP spectra

lam = np.array([  250.0, 325.0, 350.0, 375.0, 400.0, 425.0, 450.0, 475.0, 
         500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0, 
         700.0, 725.0, 775.0, 850.0, 950.0, 1050.0, 1150.0, 1250.0,
        1350.0, 1450.0, 1550.0,  1650.0, 1750.0, 1900.0, 2200.0, 2900.0, 3700.0 , 4000.0]) 
lam = lam.reshape(lam.shape[0],1)

def PFT_calc(CHL, p1, p2, p3, p4):
    #p1, p2, p3 and p4 are relative contributions (0-1)
    # of various PFT to total absorption
    PFT_1 = p1*CHL     # 
    PFT_2 = p2*CHL
    PFT_3 = p3*CHL
    PFT_4 = p4*CHL
    return PFT_1, PFT_2, PFT_3, PFT_4

def NAP_calc(PresCHL, fnap):
    NAP = fnap*np.ones(PresCHL.shape)
    return NAP

def CDOM_calc(PresCHL,fcdom):
    CDOM = fcdom*np.ones(PresCHL.shape)
    return CDOM

def NAP_abs(CHL, Snap, a443):
    a_NAP = a443 * np.exp(-Snap * (lam - 443.))
    a_NAP = a_NAP.reshape((1, a_NAP.shape[0]))
    CHL = CHL.reshape(CHL.shape[0], 1)
    aNAP = a_NAP * CHL / np.max(CHL)
    return aNAP

def CDOM_abs(CHL, Scdom, a440):
    a_cdom = a440 * np.exp(-Scdom*(lam-440.))
    a_cdom = a_cdom.reshape((1, a_cdom.shape[0]))
    CHL = CHL.reshape(CHL.shape[0], 1)
    aCDOM = a_cdom * CHL / np.max(CHL)
    return aCDOM




