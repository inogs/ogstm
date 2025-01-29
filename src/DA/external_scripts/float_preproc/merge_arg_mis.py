from bitsea.commons.utils import file2stringlist
import argparse
from bitsea.utilities.argparse_types import generic_path


def argument():
    parser = argparse.ArgumentParser(description = '''
    Merges chl and nitrate misfit files in an unique misfit file, input for 3d_var
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--nit', '-n',
                            type = generic_path,
                            required = True,
                            help = 'path of the input nitrate misfit file')
    parser.add_argument(   '--chl', '-c',
                            type = generic_path,
                            required = True,
                            help = 'path of the input chlorophyll misfit file')
    parser.add_argument(   '--oxy', '-x',
                            type = generic_path,
                            required = True,
                            help = 'path of the input oxygen misfit file')
    parser.add_argument(   '--outfile', '-o',
                            type = generic_path,
                            help = 'path of the output misfit file')

    return parser.parse_args()

args = argument()


arg_misN3n = args.nit
arg_misP_l = args.chl
arg_misO2o = args.oxy

if arg_misN3n.is_file():
    merge1 = arg_misN3n
    N3nmis = file2stringlist(arg_misN3n)
    N3nmis0 = N3nmis[1:]
else:
    merge1 = ''
    N3nmis = [0]
    N3nmis0 = []


if arg_misP_l.is_file():
    merge2 = arg_misP_l
    P_lmis = file2stringlist(arg_misP_l)
    P_lmis0 = P_lmis[1:]
else:
    merge2 = ''
    P_lmis = [0]
    P_lmis0 = []

if arg_misO2o.is_file():
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


with open(args.outfile,'wb') as ff:
    ff.writelines((item + "\n").encode() for item in allmis )

