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
    parser.add_argument(   '--oxy', '-x',
                            type = str,
                            required = True,
                            help = 'path of the input oxygen misfit file')
    parser.add_argument(   '--outfile', '-o',
                            type = str,
                            help = 'path of the output misfit file')

    return parser.parse_args()

args = argument()


arg_misN3n = args.nit
arg_misP_l = args.chl
arg_misO2o = args.oxy

exists = os.path.isfile(arg_misN3n)
if exists:
    merge1 = arg_misN3n
    N3nmis = file2stringlist(arg_misN3n)
    N3nmis0 = N3nmis[1:]
else:
    merge1 = ''
    N3nmis = [0]
    N3nmis0 = []


exists = os.path.isfile(arg_misP_l)
if exists:
    merge2 = arg_misP_l
    P_lmis = file2stringlist(arg_misP_l)
    P_lmis0 = P_lmis[1:]
else:
    merge2 = ''
    P_lmis = [0]
    P_lmis0 = []

exists = os.path.isfile(arg_misO2o)
if exists:
    merge3 = arg_misO2o
    O2omis = file2stringlist(arg_misO2o)
    O2omis0 = O2omis[1:]
else:
    merge3 = ''
    O2omis = [0]
    O2omis0 = []

totobs =  int(P_lmis[0]) + int(N3nmis[0]) + int(O2omis[0])
allmis0 = N3nmis0 + P_lmis0 + O2omis0

allmis = [str(totobs)] + allmis0


ff = open(args.outfile,'wb')
ff.writelines((item + "\n").encode() for item in allmis )

ff.close()
