import numpy as np
from commons.mask import Mask
import pylab as pl
#TheMask = Mask("/Users/gbolzon/Documents/workspace/ogs_bounday_conditions/masks/meshmask.nc")
TheMask= Mask("/gpfs/scratch/userexternal/plazzari/eas_v6/eas_v6_1/wrkdir/MODEL/meshmask.nc")
tmask = TheMask.mask_at_level(0)
jpjglo, jpiglo = tmask.shape

def riparto(lenglo,nprocs):
    ''' works in a single direction '''
    mean_value, remainder = divmod(lenglo,nprocs)
    #print "rem = ", remainder
    JP = np.ones((nprocs),np.int)*mean_value + 2
    JP[ 0] = JP[ 0] - 1
    JP[-1] = JP[-1] - 1
    for r in range(remainder):
        JP[-r] = JP[-r]+1
    return JP

def get_startpoints(JP):
    '''returns indexes in fortran format'''
    startpoint=1
    startpoints=[]
    for j in JP:
        startpoints.append(startpoint)
        startpoint = startpoint + j -2
    return startpoints




def get_wp_matrix(tmask, nprocj, nproci):
    '''
    returns 
    * M * a 2d array (nprocj, nproci) of integers containing the sum of waterpoints
    * C * 
    '''
    jpjglo, jpiglo = tmask.shape
    JPI = riparto(jpiglo,nproci)
    Start_I = get_startpoints(JPI)
    End_I = Start_I + JPI -1
    
    JPJ = riparto(jpjglo,nprocj)
    Start_J = get_startpoints(JPJ)
    End_J = Start_J + JPJ -1

    M = np.zeros((nprocj, nproci),dtype=np.int32)
    C = np.zeros((nprocj, nproci),dtype=np.int32)
    for i in range(nproci): 
        for j in range(nprocj):
            start_i = Start_I[i] -1
            end_i   = End_I[i] -1 
            start_j = Start_J[j] -1
            end_j   = End_J[j] -1
            #print start_i, end_i, start_j, end_j
            m = tmask[start_j:end_j, start_i:end_i]
            M[j,i] = m.sum()
            C[j,i] = m[0,:].sum() + m[:,0].sum()
    return M,C




def candidate_decompositions(tmask, max_proc_i,max_proc_j,nproc):
    M_processes = np.zeros((max_proc_j,max_proc_i),np.int)
    C_processes = np.zeros((max_proc_j,max_proc_i),np.int)
    for i in range(max_proc_i):
        nproci = i+1
        for j in range(max_proc_j):
            nprocj = j+1
            if (nproci * nprocj < nproc)   : continue
            if (nproci * nprocj > nproc*3 ): continue
            M,C = get_wp_matrix(tmask, nprocj, nproci)
            M_processes[j,i] = (M>0).sum()
            C_processes[j,i] = C.sum()
    return M_processes,C_processes
        
nproc = 128
max_proc_i = 20
max_proc_j = 12

USED_PROCS, COMMUNICATION = candidate_decompositions(tmask, max_proc_i, max_proc_j, nproc)
good = USED_PROCS == nproc
J,I = good.nonzero() # poi vanno incrementati di 1
nCandidates = len(I)
HYP_COMMUNICATION_LINE=np.zeros(nCandidates,dtype=np.int)
EFF_COMMUNICATION_LINE=np.zeros(nCandidates,dtype=np.int)
for k in range(nCandidates):
    nproci = I[k]+1
    nprocj = J[k]+1
    line = (nproci -1 )*jpjglo + (nprocj-1)*jpiglo
    HYP_COMMUNICATION_LINE[k]=line
    EFF_COMMUNICATION_LINE[k] = COMMUNICATION[J[k],I[k]]
    JPI = riparto(jpiglo,nproci)
    JPJ = riparto(jpjglo,nprocj)
    print (JPI.mean()==JPI[0]) , (JPJ.mean()==JPJ[0]) 


print I+1,J+1, EFF_COMMUNICATION_LINE

choosen = EFF_COMMUNICATION_LINE.argmin()
nproci  = I[choosen]+1
nprocj  = J[choosen]+1
M,C = get_wp_matrix(tmask, nprocj, nproci)
J,I = M.nonzero()


def neighbors(M,nproc):
    J,I = M.nonzero()
    WEST =np.zeros((nproc,),dtype=np.int)
    SOUTH=np.zeros((nproc,),dtype=np.int)
    EAST =np.zeros((nproc,),dtype=np.int)
    NORTH=np.zeros((nproc,),dtype=np.int)
    
    for rank in range(nproc):
        j = J[rank]
        i = I[rank]
        if i==0 :
            west = -1
        else:
            if M[j,i-1]>0:
                west = np.argwhere((J == j) & ( I == i-1))[0][0]
            else:
                west = -1
        if i==nproci-1 :
            east = -1
        else:
            if M[j,i+1]>0:
                east = np.argwhere((J == j) & ( I == i+1))[0][0]
            else:
                east = -1
    
    
        if j==0 :
            south = -1
        else:
            if M[j-1,i]>0:
                south = np.argwhere((J == j-1) & ( I == i))[0][0]
            else:
                south = -1
        if j==nprocj-1 :
            north = -1
        else:
            if M[j+1,i]>0:
                north = np.argwhere((J == j+1) & ( I == i))[0][0]
            else:
                north = -1
    
        WEST[rank] = west
        SOUTH[rank] = south
        EAST[rank] = east
        NORTH[rank] = north
        
        return WEST, EAST, NORTH, SOUTH
    
WEST, EAST, NORTH, SOUTH = neighbors(M, nproc)

for rank in range(nproc):
    print WEST[rank],rank, EAST[rank], I[rank], J[rank]

def plot_decomposition(nproci, nprocj):
    M,C = get_wp_matrix(tmask, nprocj, nproci)
    J,I = M.nonzero()

    JPI = riparto(jpiglo,nproci)
    Start_I = get_startpoints(JPI)
    End_I = Start_I + JPI -1
     
    JPJ = riparto(jpjglo,nprocj)
    Start_J = get_startpoints(JPJ)
    End_J = Start_J + JPJ -1
    
    fig,ax = pl.subplots()
    ax.imshow(tmask)
    for i in range(1,nproci):
        x=Start_I[i]
        ax.plot([x,x],[0,jpjglo],'w')
    
    for j in range(1,nprocj):
        y=Start_J[j]
        ax.plot([0,jpiglo],[y,y],'w')
    
    for rank in range(nproc):
        j = J[rank]
        i = I[rank]
        x = Start_I[i] + JPI[0]/2
        y = Start_J[j] + JPJ[0]/2
        ax.text(x,y,str(rank), color='w', ha='center', va='center', fontsize=8)
    
    
    ax.invert_yaxis()
    return fig, ax

fig, ax = plot_decomposition(nproci, nprocj)
fig.set_dpi(150)


        
