import numpy as np
from commons.utils import file2stringlist
import argparse
import os


def argument():
    parser = argparse.ArgumentParser(description = '''
    Merges chl and nitrate misfit files in an unique misfit file, input for 3d_var
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--nit', '-n',
                            type = str,
                            required = True,
                            help = 'path of the input nitrate misfit file')
    parser.add_argument(   '--chl', '-c',
                            type = str,
                            required = True,
                            help = 'path of the input chlorophyll misfit file')
    parser.add_argument(   '--outfile', '-o',
                            type = str,
                            help = 'path of the output misfit file')

    return parser.parse_args()

args = argument()


idate0=args.inputdate

arg_misN3n = args.nit
exists = os.path.isfile(arg_misN3n)
if exists:
    merge1 = arg_misN3n
    N3nmis = file2stringlist(arg_misN3n)
    N3nmis0 = N3nmis[1:]
else:
    merge1 = ''
    N3nmis = [0]
    N3nmis0 = []

arg_misP_l = args.chl
exists = os.path.isfile(arg_misP_l)
if exists:
    merge2 = arg_misP_l
    P_lmis = file2stringlist(arg_misP_l)
    P_lmis0 = P_lmis[1:]
else:
    merge2 = ''
    P_lmis = [0]
    P_lmis0 = []

totobs =  np.int(P_lmis[0]) + np.int(N3nmis[0])
allmis0 = N3nmis0 + P_lmis0

allmis = [np.str(totobs)] + allmis0


print 'Merging 1 ' + merge1 + ' and 2 ' + merge2 + ' into ' + outfile

ff = open(outfile,'wb')
ff.writelines( "%s\n" % item for item in allmis )

ff.close()
