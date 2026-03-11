import os
import pyfabm

CODEPATH = '../../'
CODEPATH = CODEPATH.replace("~",os.getenv("HOME"))
fabm_yaml=  CODEPATH + "/fabm/extern/ogs/fabm_multispectral_2xDetritus.yaml "
model = pyfabm.Model(fabm_yaml)

#interior_state:
#    O2_o:
#        ctrmax: 1.000000e+03 
#        ctrhf:  1
#        relax:  1

with open('ogstm.yaml', 'w') as f:
    f.write('interior_state:\n')
    for i,variable in enumerate(model.state_variables):
        f.write(f"    {variable.name.replace('/','_')}:\n")
        f.write(f"        ctrmax: 1.000000e+03\n")
        # OGSTM want always at least an high freq output
        if i == 0:
            f.write(f"        ctrhf: 1\n")
        else:
            f.write(f"        ctrhf: 0\n")
        f.write(f"        relax: 0\n")

    f.write('interior_diagnostic:\n')
    for variable in model.interior_diagnostic_variables:
        if variable.output:
            f.write(f"    {variable.name.replace('/','_')}:\n")
            f.write(f"        diahf: 0\n")
            f.write(f"        diaWR: 0\n")

    f.write('horizontal_diagnostic:\n')
    for variable in model.horizontal_diagnostic_variables:
        if variable.output:
            f.write(f"    {variable.name.replace('/','_')}:\n")
            f.write(f"        diahf_2d: 0\n")
            f.write(f"        diaWR_2d: 0\n")
