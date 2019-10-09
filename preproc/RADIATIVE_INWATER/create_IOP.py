#eterzic 07.10.2019

# Create functions for chl-specific IOP spectra

def PFT_calc(CHL, p1, p2, p3, p4):
    #p1, p2, p3 and p4 are relative contributions (0-1)
    # of various PFT to total absorption
    PFT_1 = p1*CHL     # 
    PFT_2 = p2*CHL
    PFT_3 = p3*CHL
    PFT_4 = p4*CHL
    return PFT_1, PFT_2, PFT_3, PFT_4

def NAP_calc(CHL, fnap):
    NAP = fnap*CHL/CHL
    return NAP

def CDOM_calc(CHL,fcdom):
    CDOM = fcdom*CHL/CHL
    return CDOM




