# Generates namelist.optics on stdout
freq_nanom = [250, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625,
       650, 675, 700, 725, 775, 850, 950, 1050, 1150, 1250, 1350, 1450, 1550, 1650,
       1750, 1900, 2200, 2900, 3700]
#   Ednm(1)="Ed_0250"
#   EdWR(1)=1
#   Ed3D(1)=0

for groupvar in ["Ed", "Es", "Eu"]:
    a="&%s_nam" %(groupvar)
    print(a)
    for ie, freq in enumerate(freq_nanom):
        a="   %snm(%d)=\"%s_%04d\"" %(groupvar,ie+1,groupvar,freq)
        print(a)
        a="   %sWR(%d)=1" %(groupvar,ie+1)
        print(a)
        n=0
        if freq in [375,400,425,475,500]: n=1
        a="   %s3D(%d)=%d\n" %(groupvar,ie+1,n)
        print(a)
    print("/\n\n")
    
