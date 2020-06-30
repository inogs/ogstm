#eterzic 07.10.2019
import numpy as np
# Create functions for chl-specific IOP spectra

# Wavelengths of the Radiative Transfer Model
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

# Di Cicco et al., 2017 - both for PFT and PSC
def PFT_MED(CHL):
    x     = np.where(CHL>0, np.log10(CHL), 0.)
    #x = np.log10(CHL)
    MICRO = 0.0667*x**3 + 0.1939*x**2 + 0.2743*x + 0.2994
    NANO  =              -0.1740*x**2 - 0.0851*x + 0.4725
    PICO  = 1 - MICRO - NANO
    DIATOM = 0.0482*x**3 + 0.1877*x**2 + 0.2946*x + 0.2533
    DINOPH = MICRO - DIATOM
    CRYPT  = 0.0171*x**3 + 0.0667*x**2 + 0.1153*x + 0.0952
    GREEN  = (np.exp(-1.5780*x + 2.1841) + 22.6833 *x) ** (-1.)
    PROK   = 0.0664*x**3 + 0.1410*x**2 - 0.2097*x + 0.0979
    HAPT   = 1. - MICRO - CRYPT - GREEN - PROK
    
    PFT_1 = DIATOM * CHL
    PFT_2 = HAPT   * CHL
    PFT_3 = CRYPT  * CHL
    PFT_4 = GREEN  * CHL
    PFT_5 = PROK   * CHL
    PFT_6 = DINOPH  * CHL

    PSC_1 = MICRO  * CHL
    PSC_2 = NANO   * CHL
    PSC_3 = PICO   * CHL
    return PFT_1, PFT_2, PFT_3, PFT_4, PFT_5, PFT_6, PSC_1, PSC_2, PSC_3
    
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




