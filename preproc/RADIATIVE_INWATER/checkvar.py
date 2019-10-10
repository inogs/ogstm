

def findVars(Varlist):
    allvars=' PAR CHL IRR_380 IRR_412 IRR_490'
    if len(Varlist)==0: return False    
    
    for var in Varlist:
        if allvars.find(var) == -1 :
            return False
    return True

