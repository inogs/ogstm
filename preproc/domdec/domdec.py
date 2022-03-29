import numpy as np
import pylab as pl

def riparto(lenglo,nprocs):
    ''' Uniform decomposition of a 1d array of size lenglo in nprocs subdomains
        The load balancing of advection and hdf depends on that decomposition,
        that is supposed uniform.

    Arguments :
     * lenglo * integer, longitudinal or latitudinal dimension of the global mesh
     * nprocs * integer, the number of processors along a dimension
    Features :
    -  the ghost cell is taken in account
    - subdomains a bit largers (one cell) are the last ones

    Returns a numpy array of integers, called jpi or jpj in ogstm '''
    mean_value, remainder = divmod(lenglo,nprocs)
    #print("rem = ", remainder)
    JP = np.ones((nprocs),np.int)*mean_value + 2
    JP[ 0] = JP[ 0] - 1
    JP[-1] = JP[-1] - 1
    for r in range(remainder):
        JP[-r] = JP[-r]+1
    return JP

def get_startpoints(JP):
    '''
    Calculates startpoints, called nimpp or njmpp in ogstm
    IO/IOnc.f90:      start    = (/nimpp, njmpp,  1,  1/)

    Argument:
     * JP* the list of jpi (or jpj)

    Returns:
    * indexes * list of integers in fortran format'''
    startpoint=1
    startpoints=[]
    for j in JP:
        startpoints.append(startpoint)
        startpoint = startpoint + j -2
    return startpoints




def get_wp_matrix(tmask, nprocj, nproci):
    '''
    Generates a Waterpoint matrix over the domain decomposition obtained by
    the number of processors in each direction
    Arguments:
    * tmask  * a 2d logical array, the surface tmask
    * nprocj * integer, number of latitudinal subdivisions
    * nproci * integer, number of longitudinal subdivisions


    The waterpoint number is calculated on surface, in order to detect
    land processors

    Returns:
    * M * a 2d array (nprocj, nproci) of integers containing the sum of waterpoints
          useful to detect land processors
    * C * a 2d array (nprocj, nproci) of integers containing, the sum of waterpoints on west and south boundary
          useful to have an idea of MPI communication
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
            #print(start_i, end_i, start_j, end_j)
            m = tmask[start_j:end_j, start_i:end_i]
            M[j,i] = m.sum()
            C[j,i] = m[0,:].sum() + m[:,0].sum()
    return M,C




def candidate_decompositions(tmask, max_proc_i,max_proc_j,nproc):
    '''
    Calculates the number of needed ranks for all the possible decompositions
    we can generate by fixing the maximum number of decompositions in each direction.
    A decomposition is considered candidate if nproc < nproci*nprocj < nproc*3

    In general, there are many decompositions for nproc ranks, so we need
     - to find them
     - then to choice the best.

    Arguments:
     * tmask      * a 2d logical array, the surface tmask
     * max_proc_i * integer, a maximum number of longitudinal subdomains
     * max_proc_j * integer, a maximum number of latitudinal subdomains
     * nproc      * the number of processors effectively used in simulation


    Returns:
    * Needed_procs * a 2d integer array (max_proc_j,max_proc_i)
                    Needed_procs[nprocj,nproci] is the number of no-land processors
                    for a (nprocj,nproci) decomposition
                    Needed_procs == nproc will be the next step candidate decomposition.

    * Comm_table * a 2d integer array (max_proc_j,max_proc_i)
                    Comm_table[nproci,nprocj] is the MPI communication,
                    useful to choice between candidates.

    '''
    Needed_procs = np.zeros((max_proc_j,max_proc_i),np.int)
    Comm_table = np.zeros((max_proc_j,max_proc_i),np.int)
    for i in range(max_proc_i):
        nproci = i+1
        for j in range(max_proc_j):
            nprocj = j+1
            if (nproci * nprocj < nproc)   : continue
            if (nproci * nprocj > nproc*3 ): continue
            M,C = get_wp_matrix(tmask, nprocj, nproci)
            Needed_procs[j,i] = (M>0).sum()
            Comm_table[j,i] = C.sum()
    good = Needed_procs == nproc
    if good.sum()==0:
        print("No valid candidate have been found. Try modify max_proc_i and/or max_proc_j.")
        raise ValueError
    return Needed_procs,Comm_table
        
def neighbors(M,nproc,nproci,nprocj):
    '''
    Generates number of neighbors ranks for each rank,
    corresponding to nowe, noea, nono, noso in ogstm.


    Arguments:
    * M     * a 2d array of integers (nproci, nprocj), as provided by get_wp_matrix
    * nproc * integer, the number of MPI ranks
    * nproci* integer, subdivisions along i
    * nprocj* integer, subdivisions along j

    This method is tested for a M waterpoint matrix associated to nproc, M should be the best choice.

    Returns:
    * WEST, EAST, NORTH, SOUTH, NBONDI, NBONDJ *   1d arrays of integers (nproc)
    '''
    J,I = M.nonzero()
    WEST =np.zeros((nproc,),dtype=np.int)
    SOUTH=np.zeros((nproc,),dtype=np.int)
    EAST =np.zeros((nproc,),dtype=np.int)
    NORTH=np.zeros((nproc,),dtype=np.int)
    NBONDI=np.zeros((nproc,),dtype=np.int)
    NBONDJ=np.zeros((nproc,),dtype=np.int)
    
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
        
        nbondi=2
        if (east>  -1) & (west>  -1) : nbondi= 0
        if (east== -1) & (west>  -1) : nbondi= 1
        if (east>  -1) & (west== -1) : nbondi=-1
        nbondj=2
        if (south>  -1) & (north>  -1) : nbondj= 0
        if (south>  -1) & (north== -1) : nbondj= 1
        if (south== -1) & (north>  -1) : nbondj=-1


        NBONDI[rank] = nbondi
        NBONDJ[rank] = nbondj
        WEST[  rank] = west
        SOUTH[ rank] = south
        EAST[  rank] = east
        NORTH[ rank] = north
    return WEST, EAST, NORTH, SOUTH,NBONDI, NBONDJ
    


def plot_decomposition(tmask, nproci, nprocj):
    '''
    Plots the domain decomposition scheme

    Arguments :
    * tmask  * a 2d logical array, the surface tmask
    * nproci * integer, number of longitudinal subdivisions
    * nprocj * integer, number of latitudinal subdivisions

    Returns:
    fig, ax : matplotlib handles
    '''
    jpjglo, jpiglo = tmask.shape
    M,C = get_wp_matrix(tmask, nprocj, nproci)
    J,I = M.nonzero()
    nproc = len(I)

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



def get_best_decomposition(USED_PROCS, COMMUNICATION, max_nproc, jpiglo, jpjglo):
    '''
    Choose of best decomposition on the basis of:
     - fit with the
     - minor MPI communication
    Arguments:
    * USED_PROCS    * output of  candidate_decomposition()
    * COMMUNICATION * idem
    * max_nproc     * maximum number of MPI ranks we can use in simulation
    * jpiglo        * integer, global domain size
    * jpjglo        * integer, global domain size

    Returns:
     * nproci, nprocj * integers
    '''

    CANDIDATES = np.zeros((0,8),np.int)
    min_nproc = max_nproc -8
    iCandidate=0
    for nproc in range(max_nproc,min_nproc,-1):
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
            xm,ym = np.meshgrid(JPI,JPJ)
            load = (xm*ym).max()
            linearray = np.array([nproc, nproci, nprocj, JPI.max(), JPJ.max(), load, EFF_COMMUNICATION_LINE[k], iCandidate],dtype=np.int,ndmin=2)
            CANDIDATES=np.concatenate((CANDIDATES,linearray),axis=0)
            iCandidate = iCandidate+1

    nCandidates,_ = CANDIDATES.shape
    print("There are ", nCandidates, "candidate decompositions")
    print("nproc, nproci, nprocj, jpi, jpj, jpi*jpj, COMM, progr")
    print(CANDIDATES)
    totalwork = CANDIDATES[:,5]+CANDIDATES[:,6] # domain_size+communication
    choosen = totalwork.argmin()
    print ("choosen=", choosen)
    nprocs  = CANDIDATES[choosen,0]
    nproci  = CANDIDATES[choosen,1]
    nprocj  = CANDIDATES[choosen,2]
    return nprocs, nproci, nprocj


def waterpoints_3d(maskobj, nprocj, nproci):
    '''
    Info about load balance of BFM calls
    '''
    jpjglo, jpiglo = tmask.shape
    JPI = riparto(jpiglo,nproci)
    Start_I = get_startpoints(JPI)
    End_I = Start_I + JPI -1

    JPJ = riparto(jpjglo,nprocj)
    Start_J = get_startpoints(JPJ)
    End_J = Start_J + JPJ -1

    M = np.zeros((nprocj, nproci),dtype=np.int32)

    for i in range(nproci):
        for j in range(nprocj):
            start_i = Start_I[i] -1
            end_i   = End_I[i] -1
            start_j = Start_J[j] -1
            end_j   = End_J[j] -1
            #print(start_i, end_i, start_j, end_j)
            m = maskobj.mask[:,start_j:end_j, start_i:end_i]
            M[j,i] = m.sum()
    return M


def dump_outfile(tmask, choosen_procs,nproci, nprocj, filename="domdec.txt"):
    jpjglo, jpiglo = tmask.shape
    M,C = get_wp_matrix(tmask, nprocj, nproci)
    J,I = M.nonzero()
    JPI = riparto(jpiglo,nproci)
    JPJ = riparto(jpjglo,nprocj)
    Start_I = get_startpoints(JPI) #nimpp
    Start_J = get_startpoints(JPJ) #njmpp
    
    
    
    WEST, EAST, NORTH, SOUTH, NBONDI,NBONDJ = neighbors(M, choosen_procs,nproci, nprocj)
    
    OUT = np.zeros((choosen_procs,13), dtype=np.int32)
    for rank in range(choosen_procs):
        i = I[rank]
        j = J[rank]
        jpi = JPI[i]
        jpj = JPJ[j]
        nimpp = Start_I[i]
        njmpp = Start_J[j]
        OUT[rank, 0] = rank
        OUT[rank, 1] = i
        OUT[rank, 2] = j
        OUT[rank, 3] = jpi
        OUT[rank, 4] = jpj
        OUT[rank, 5] = nimpp
        OUT[rank, 6] = njmpp
        OUT[rank, 7] = NBONDI[rank]
        OUT[rank, 8] = NBONDJ[rank]
        OUT[rank, 9] = WEST[rank]
        OUT[rank,10] = EAST[rank]
        OUT[rank,11] = NORTH[rank]
        OUT[rank,12] = SOUTH[rank]
    np.savetxt(filename, OUT, fmt=13*"%5d")

