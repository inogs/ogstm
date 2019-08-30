from commons.utils import file2stringlist

def dump_template(ORIG, outfile,SAT_OBS,ARGO):
    LINES=[]
    for line in ORIG:
        newline=line
        if (line.find("@@SAT_OBS@@") != -1):  newline=line.replace("@@SAT_OBS@@",SAT_OBS)
        if (line.find("@@ARGO@@")   != -1):   newline=line.replace("@@ARGO@@"   ,  ARGO)
        LINES.append(newline + "\n")
    fid=open(filename,"w")
    fid.writelines(LINES)
    fid.close()

ORIG=file2stringlist('satfloat.template')
dateDAfloat = file2stringlist("daFloat")
dateDAsat   = file2stringlist("daSat")
dateDA      = file2stringlist("daTimes")

OUTDIR="OUT/"
for d in dateDA:
    filename=OUTDIR + "satfloat." + d[:8] + ".nml"
    #print filename
    ARGO    = str(int(d in dateDAfloat)) # "1" if True, "0" if False
    SAT_OBS = str(int(d in dateDAsat))
    dump_template(ORIG, filename, SAT_OBS, ARGO)





