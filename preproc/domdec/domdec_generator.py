import argparse
 
def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates domdec.txt file for ogstm-bfm model, a file containing
    this line for each rank:
    rank i_pos j_pos jpi jpj nimpp njmpp nbondi nbondj west east north south neighbors
    ''', formatter_class=argparse.RawTextHelpFormatter)
 
 
    parser.add_argument(   '--mpiprocs','-n',
                                type = int,
                                required = True,
                                help = 'number of MPI ranks of ogstm simulation')
    parser.add_argument(   '--max_proc_i','-i',
                                type = int,
                                required = True,
                                help = 'Maximum number of subdivision along i')
    parser.add_argument(   '--max_proc_j','-j',
                                type = int,
                                required = True,
                                help = 'Maximum number of subdivision along j')
 
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = ''' Path of maskfile''')
 
 
    return parser.parse_args()
 
args = argument()
nproc = args.mpiprocs
max_proc_i = args.max_proc_i
max_proc_j = args.max_proc_j

from domdec import *
from commons.mask import Mask

TheMask = Mask(args.maskfile, dzvarname="e3t_0")
tmask = TheMask.mask_at_level(0)
jpjglo, jpiglo = tmask.shape


USED_PROCS, COMMUNICATION = candidate_decompositions(tmask, max_proc_i, max_proc_j, nproc)

choosen_procs, nproci, nprocj =  get_best_decomposition(USED_PROCS, COMMUNICATION, nproc, jpiglo, jpjglo)

dump_outfile(tmask, choosen_procs,nproci, nprocj)

fig, ax = plot_decomposition(tmask, nproci, nprocj)
fig.set_dpi(150)
outfile='domdec_' + str(choosen_procs) + ".png"
fig.savefig(outfile)
