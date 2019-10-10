

def findVars(Varlist):
    allvars=[' CHL', 'IRR_380', 'IRR_412', 'IRR_490', 'PAR']
    if len(Varlist)==0: return False    
    
    for var in allvars:
        if Varlist.find(var) == -1 :
            return False
    return True

