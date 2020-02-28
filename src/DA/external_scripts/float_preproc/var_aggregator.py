import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates ave files with aggregated var for aveScan.py. 
    Avescan.py. Creates also chl sup.''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '/some/path/MODEL/AVE_FREQ_1/')

    parser.add_argument(   '--tmpdir', '-t',
                                type = str,
                                default = None,
                                help = """ /some/path/POSTPROC/output/AVE_FREQ_1/TMP/ .  
                                Path to put files with aggregated variables for aveScan.py. 
                                No generation of aggregate files if this parameter is omitted.
                                """)
    parser.add_argument(   '--archivedir', '-a',
                                type = str,
                                default = None,
                                help = '''/some/path/POSTPROC/output/AVE_FREQ_1/Archive/  . 
                                Path to put native vars as they are, in order to compress them.
                                No generation of archived files if this parameter is omitted.
                                ''')
    
    parser.add_argument(   '--chlsupdir', '-c',
                                type = str,
                                default = None,
                                help = """/some/path/POSTPROC/output/AVE_FREQ_1/CHL_SUP.
                                No generation of chl sup if this parameter is omitted.
                                """)    
    parser.add_argument(   '--avelist',"-l",
                                type = str,
                                default = "ave*N1p.nc",
                                help = 'ave*.N1p.nc, they configure the date list')
    parser.add_argument(   '--descriptor',"-d",
                                type = str,
                                required = True,
                                help = 'VarDescriptor_1.xml, or the complete path')
    parser.add_argument(   '--maskfile',"-m",
                                type = str,
                                help = '''Path of the mask file''')

    return parser.parse_args()

args = argument()

import glob
import os
import GB_lib as G
import read_descriptor
from commons.mask import Mask

try :
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks = comm.size 
except:
    rank   = 0
    nranks = 1

def addsep(string):    
    if string[-1] != os.sep:
        return string + os.sep
    else:
        return  string
    
 

AVEDIR     = addsep(args.inputdir)

RD = read_descriptor.read_descriptor(args.descriptor)
PATH_NAME = AVEDIR + args.avelist
if rank==0 : print "INPUT_DIR =", AVEDIR


if args.archivedir:
    ARCHIVEdir = addsep(args.archivedir)
    if rank==0 : print "ARCHIVEDIR=", ARCHIVEdir
    os.system("mkdir -p " + ARCHIVEdir)

if args.tmpdir:
    TMPOUTdir  = addsep(args.tmpdir)
    if rank==0 : print "TMPOUTDIR= ", TMPOUTdir
    os.system("mkdir -p " + TMPOUTdir)

if args.chlsupdir:
    CHLSUPdir  = addsep(args.chlsupdir)
    if rank==0 : print "CHLSUPDIR =", CHLSUPdir
    os.system("mkdir -p " + CHLSUPdir)

SingleVar_filelist=glob.glob(PATH_NAME)
SingleVar_filelist.sort()
TheMask=Mask(args.maskfile)

for N1pfile in SingleVar_filelist[rank::nranks]:
    dailyAve  = os.path.basename(N1pfile)
    print "writing ", dailyAve

    if args.tmpdir:
        G.WriteAggregateAvefiles(TheMask, N1pfile, AVEDIR, TMPOUTdir, TMPOUTdir, RD)
    
    if args.chlsupdir:    
        F = G.filename_manager(N1pfile)
        chl3dfile = TMPOUTdir + F.prefix + "." + F.datestr + ".P_l.nc"
        chl2dfile = CHLSUPdir + "chl."         + F.datestr + ".nc"
        G.writeChlSup(chl3dfile, chl2dfile, 'P_l')
