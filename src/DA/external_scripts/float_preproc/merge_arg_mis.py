import numpy as np
from commons.utils import file2stringlist
import argparse
import os


def argument():
    parser = argparse.ArgumentParser(description = '''
    read something
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--inputdate', '-t',
                            type = str,
                            required = True,
                            default = '20150101',
                            help = 'Input date')

    return parser.parse_args()

args = argument()


idate0=args.inputdate

arg_misN3n = idate0 + ".N3n_arg_mis.dat"
exists = os.path.isfile(arg_misN3n)
if exists:
    merge1 = arg_misN3n
    N3nmis = file2stringlist(arg_misN3n)
    N3nmis0 = N3nmis[1:]
else:
    merge1 = ''
    N3nmis = [0]
    N3nmis0 = []

arg_misP_l = idate0 + ".P_l_arg_mis.dat"
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


filename = idate0 + '.arg_mis.dat'

print 'Merging 1 ' + merge1 + ' and 2 ' + merge2 + ' into ' + filename 

ff = open(filename,'wb')
ff.writelines( "%s\n" % item for item in allmis )

ff.close()
