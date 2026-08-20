import numpy as np
try:
    import pylab as pl
except ImportError:
    pl = None

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
    # estimate how many points must be assigned to each riparto cell
    mean_value, remainder = divmod(lenglo,nprocs)
    print("rem = ", remainder)
    print("mean_value = ", mean_value)
    # generate the riparto
    # the 2 additional term adds lines of halo points (columns/rows)
    JP = np.ones((nprocs),int)*mean_value + 2
    # subtract lines of halo points
    JP[ 0] = JP[ 0] - 1    # south-most/west-most subdomain
    JP[-1] = JP[-1] - 1    # north-most/east-most subdomain

    # distribute the remaining amount of points
    # it starts from the initial subdomain: this avoids the possibility to have
    # adjacent groups of points with a higher-than-1 difference
    for r in range(remainder):
        JP[-r] = JP[-r]+1
    return JP

def riparto_global(lenglo,nprocs,axis):
    ''' Uniform decomposition of a 1d array of size lenglo in nprocs subdomains
        The load balancing of advection and hdf depends on that decomposition,
        that is supposed uniform.
        The decomposition is different for the 'lat' axis: it removes one stripe of ghost cells

    Arguments :
     * lenglo * integer, longitudinal or latitudinal dimension of the global mesh
     * nprocs * integer, the number of processors along a dimension
     * axis * string, 'lat' or 'lon'
    Features :
    -  the ghost cell is taken in account
    - subdomains a bit largers (one cell) are the last ones

    Returns a numpy array of integers, called jpi or jpj in ogstm '''

    if axis == 'lat':
        print("Riparto in the south-north direction")
        # estimate how many points must be assigned to each riparto cell
        mean_value, remainder = divmod(lenglo,nprocs)
        print("rem = ", remainder)
        print("mean_value = ", mean_value)
        # generate the riparto
        # the 2 additional term adds inner lines of halo points (rows)
        JP = np.ones((nprocs),int)*mean_value + 2
        # subtract lines of halo points
        JP[ 0] = JP[ 0] - 1    # south-most subdomain
        JP[-1] = JP[-1] - 1    # north-most subdomain

        # distribute the remaining amount of points
        # it starts from the initial subdomain: this avoids the possibility to have
        # adjacent subdomains of points with a higher-than-1 difference
        for r in range(remainder):
            JP[-r] = JP[-r]+1
        return JP

    elif axis == 'lon':
        print("Riparto in the west-east direction")
        # estimate how many points must be assigned to each riparto cell
        # subtract 1 to exclude the west-most column
        mean_value, remainder = divmod(lenglo,nprocs)
        print("rem = ", remainder)
        print("mean_value = ", mean_value)
        # generate the riparto
        # the 2 additional term adds inner lines of halo points (columns)
        JP = np.ones((nprocs),int)*mean_value + 2
        # subtract lines of halo points
        JP[ 0] = JP[ 0] - 1    # west-most subdomain
        JP[-1] = JP[-1] - 1    # east-most subdomain

        # distribute the even remaining amount of points
        # assumption: the remainder is even, so that the distribution is symmetric
        # it starts from the subdomain in the middle (it is not incremented)
        # there could be adjacent subdomains with a maximum difference of 2 lines
        _, remDIV2 = divmod(nprocs,2)
        if remDIV2 != 0:
            idx_center_sd = int((nprocs-1)/2)
            for r in range(int(remainder/2)):
                idx_shift = r + 1
                idx_w = idx_center_sd - idx_shift
                idx_e = idx_center_sd + idx_shift
                JP[idx_w] = JP[idx_w]+1
                JP[idx_e] = JP[idx_e]+1
        else:
            idx_e_center_start = int(nprocs/2)
            idx_w_center_start = idx_e_center_start - 1
            for r in range(int(remainder/2)):
                idx_shift = r
                idx_w = idx_w_center_start - idx_shift
                idx_e = idx_e_center_start + idx_shift
                JP[idx_w] = JP[idx_w]+1
                JP[idx_e] = JP[idx_e]+1
        return JP

    else:
        print("Wrong value for the axis argument: it must be  'lat' or 'lon'.")
        raise  ValueError()

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


def calculate_riparto(jpjglo, jpiglo, nprocj, nproci, map_is_global):
    '''
    Computes the riparto scheme: it assigns a number of points to each process, considering both directions.
    It computes the start and end indices for each process along the two directions.
    '''
    if map_is_global:
        JPI = riparto_global(jpiglo,nproci,'lon')
        JPJ = riparto_global(jpjglo,nprocj,'lat')
    else:
        JPI = riparto(jpiglo,nproci)
        JPJ = riparto(jpjglo,nprocj)
    Start_I = get_startpoints(JPI)
    Start_J = get_startpoints(JPJ)
    End_I = Start_I + JPI -1
    End_J = Start_J + JPJ -1

    return JPI, JPJ, Start_I, Start_J, End_I, End_J


def get_wp_matrix(tmask, nprocj, nproci, Start_I, End_I, Start_J, End_J):
    '''
    Generates a Waterpoint matrix over the domain decomposition obtained by
    the number of processors in each direction
    Arguments:
    * tmask  * a 2d logical array, the surface tmask
    * nprocj * integer, number of latitudinal subdivisions
    * nproci * integer, number of longitudinal subdivisions
    * map_is_global * boolean, indicates if the decomposition is performed on a global map


    The waterpoint number is calculated on surface, in order to detect
    land processors

    Returns:
    * M * a 2d array (nprocj, nproci) of integers containing the sum of waterpoints
          useful to detect land processors
    * C * a 2d array (nprocj, nproci) of integers containing, the sum of waterpoints on west and south boundary
          useful to have an idea of MPI communication
    '''

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




def candidate_decompositions(tmask, max_proc_i,max_proc_j,nproc,map_is_global):
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
     * map_is_global * boolean, indicates if the decomposition is performed on a global map


    Returns:
    * Needed_procs * a 2d integer array (max_proc_j,max_proc_i)
                    Needed_procs[nprocj,nproci] is the number of no-land processors
                    for a (nprocj,nproci) decomposition
                    Needed_procs == nproc will be the next step candidate decomposition.

    * Comm_table * a 2d integer array (max_proc_j,max_proc_i)
                    Comm_table[nproci,nprocj] is the MPI communication,
                    useful to choice between candidates.

    '''
    Needed_procs = np.zeros((max_proc_j,max_proc_i),int) # this is a 2D matrix where each element corresponds to a pair of subdomain count values (each value corresponds to an axis)
    Comm_table = np.zeros((max_proc_j,max_proc_i),int)

    jpjglo, jpiglo = tmask.shape
    print("Global domain size: ", jpiglo, "x", jpjglo)
    for i in range(max_proc_i):
        nproci = i+1
        if map_is_global:
            # conditions that must be fulfilled
            # 1. even number of subdomains
            # 2. even number of points in the remainder
            # discard all the odd values for nproci
            _, remDIV2   = divmod(nproci,2)
            if remDIV2 != 0: continue

            # check if jpiglo divided by nproci has an even remainder 
            _, remainder = divmod(jpiglo,nproci)
            _, remDIV2   = divmod(remainder,2)
            if (remDIV2 != 0): continue

        for j in range(max_proc_j):
            nprocj = j+1
            if (nproci * nprocj < nproc)   : continue
            if (nproci * nprocj > nproc*3 ): continue
            JPI, JPJ, Start_I, Start_J, End_I, End_J = calculate_riparto(jpjglo, jpiglo, nprocj, nproci, map_is_global)
            print("Candidate decomposition: ", nproci, "x", nprocj, "=", nproci*nprocj, "processes")
            print("  S-N (j) decomposition")
            print(f"    {JPJ}")
            print(f"    {Start_J}")
            print(f"    {End_J}")

            print("  W-E (i) decomposition")
            print(f"    {JPI}")
            print(f"    {Start_I}")
            print(f"    {End_I}")
            M,C = get_wp_matrix(tmask, nprocj, nproci, Start_I, End_I, Start_J, End_J)
            Needed_procs[j,i] = (M>0).sum() # the non-zero values contain the necessary number of processes: each process is associated to a subdomain in which there is at least one waterpoint
            Comm_table[j,i] = C.sum()
    print("needed procs matrix - only not zero elements")
    print("    S-N (j)   ||     E-W (i)     ||  needed processes")
    for i in range(max_proc_i):
        for j in range(max_proc_j):
            if Needed_procs[j,i] != 0:
                print(f"   {j+1}     ||    {i+1}         ||  {Needed_procs[j,i]}")
    good = Needed_procs == nproc # counts the combinations that use all the processes: if there are combinations that use less than nproc?
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
    WEST =np.zeros((nproc,),dtype=int)
    SOUTH=np.zeros((nproc,),dtype=int)
    EAST =np.zeros((nproc,),dtype=int)
    NORTH=np.zeros((nproc,),dtype=int)
    NBONDI=np.zeros((nproc,),dtype=int)
    NBONDJ=np.zeros((nproc,),dtype=int)
    
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


def neighbors_global(M,nproc,nproci,nprocj):
    '''
    Generates number of neighbors ranks for each rank,
    corresponding to west, east, nord, south in ogstm.
    Case of global map covered by tripolar ORCA grid, with a north-boundary folded onto itself.

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
    WEST =np.zeros((nproc,),dtype=int)
    SOUTH=np.zeros((nproc,),dtype=int)
    EAST =np.zeros((nproc,),dtype=int)
    NORTH=np.zeros((nproc,),dtype=int)
    NBONDI=np.zeros((nproc,),dtype=int)
    NBONDJ=np.zeros((nproc,),dtype=int)
    NORTH_BND=np.zeros((nproc,),dtype=int)
    
    for rank in range(nproc):
        j = J[rank]
        i = I[rank]

        if i==0:
            i_west = -1+nproci
        else:
            i_west = i-1
        if M[j,i_west]>0:
            west = np.argwhere((J == j) & ( I == i_west))[0][0]
        else:
            west = -1

        if i==nproci-1:
            i_east = i+1-nproci
        else:
            i_east = i+1
        if M[j,i_east]>0:
            east = np.argwhere((J == j) & ( I == i_east))[0][0]
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
            j_north = j
            i_north = nproci-1-i   # opposite domain, wrt the domain in the middle (with index (nproci-1)/2 )
            north_bnd = 1
        else:
            j_north = j+1
            i_north = i
            north_bnd = 0
        # activate the following 'if block' to lock communication of the process assigned to the domain in the middle with itself
        if j_north == j and i_north == i:
            north = -1
        elif M[j_north,i_north]>0:
            north = np.argwhere((J == j_north) & ( I == i_north))[0][0]
        else:
            north = -1

        # codes to detect the number of boundaries that communicate with adjacent subdomains
        # west-east direction
        nbondi=2
        if (east>  -1) & (west>  -1) : nbondi= 0
        if (east== -1) & (west>  -1) : nbondi= 1
        if (east>  -1) & (west== -1) : nbondi=-1
        # south-north direction
        nbondj=2
        if (south>  -1) & (north>  -1) : nbondj= 0
        if (south>  -1) & (north== -1) : nbondj= 1
        if (south== -1) & (north>  -1) : nbondj=-1
        # store the codes into the rank-based list
        NBONDI   [rank] = nbondi
        NBONDJ   [rank] = nbondj
        WEST     [rank] = west
        SOUTH    [rank] = south
        EAST     [rank] = east
        NORTH    [rank] = north
        NORTH_BND[rank] = north_bnd
    return WEST, EAST, NORTH, SOUTH,NBONDI, NBONDJ, NORTH_BND
    


def plot_decomposition(tmask, nproci, nprocj, map_is_global):
    '''
    Plots the domain decomposition scheme

    Arguments :
    * tmask  * a 2d logical array, the surface tmask
    * nproci * integer, number of longitudinal subdivisions
    * nprocj * integer, number of latitudinal subdivisions
    * map_is_global * boolean, indicates if the decomposition is performed on a global map

    Returns:
    fig, ax : matplotlib handles (None, None if matplotlib is not available)
    '''
    if pl is None:
        print("matplotlib is not available: skipping plot_decomposition().")
        return None, None

    jpjglo, jpiglo = tmask.shape
    JPI, JPJ, Start_I, Start_J, End_I, End_J = calculate_riparto(jpjglo, jpiglo, nprocj, nproci, map_is_global)

    M,C = get_wp_matrix(tmask, nprocj, nproci, Start_I, End_I, Start_J, End_J)
    J,I = M.nonzero()
    nproc = len(I)
    
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



def get_best_decomposition(USED_PROCS, COMMUNICATION, max_nproc, jpiglo, jpjglo, map_is_global):
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
    * map_is_global * boolean, indicates if the decomposition is performed on a global map

    Returns:
     * nproci, nprocj * integers
    '''

    CANDIDATES = np.zeros((0,8),int)
    min_nproc = max_nproc -8
    iCandidate=0
    for nproc in range(max_nproc,min_nproc,-1):
        good = USED_PROCS == nproc     # select the combinations that use a given number of processes - this number is reduced by 1 at each iteration

        J,I = good.nonzero() # poi vanno incrementati di 1  # stores into different arrays the number of domain subdivisions on the two axes
        nCandidates = len(I)

        HYP_COMMUNICATION_LINE=np.zeros(nCandidates,dtype=int) # hypothetical
        EFF_COMMUNICATION_LINE=np.zeros(nCandidates,dtype=int) # effective

        for k in range(nCandidates):
            nproci = I[k]+1      # reconstruct the number of subdivisions by adding 1
            nprocj = J[k]+1      # reconstruct the number of subdivisions by adding 1
            line = (nproci -1 )*jpjglo + (nprocj-1)*jpiglo
            HYP_COMMUNICATION_LINE[k]=line
            EFF_COMMUNICATION_LINE[k] = COMMUNICATION[J[k],I[k]]
            JPI, JPJ, _, _, _, _ = calculate_riparto(jpjglo, jpiglo, nprocj, nproci, map_is_global)
            xm,ym = np.meshgrid(JPI,JPJ)
            load = (xm*ym).max()
            linearray = np.array([nproc, nproci, nprocj, JPI.max(), JPJ.max(), load, EFF_COMMUNICATION_LINE[k], iCandidate],dtype=int,ndmin=2)
            CANDIDATES=np.concatenate((CANDIDATES,linearray),axis=0)
            iCandidate = iCandidate+1

    nCandidates,_ = CANDIDATES.shape
    print("There are ", nCandidates, "candidate decompositions")
    print("nproc, nproci, nprocj, jpi, jpj, jpi*jpj, COMM, progr")
    print(CANDIDATES)
    totalwork = CANDIDATES[:,5]+CANDIDATES[:,6] # domain_size + effective communication
    choosen = totalwork.argmin()
    print ("choosen=", choosen)
    nprocs  = CANDIDATES[choosen,0]
    nproci  = CANDIDATES[choosen,1]
    nprocj  = CANDIDATES[choosen,2]
    return nprocs, nproci, nprocj


def waterpoints_3d(tmask, maskobj, nprocj, nproci, map_is_global):
    '''
    Info about load balance of BFM calls
    '''
    jpjglo, jpiglo = tmask.shape
    _, _, Start_I, Start_J, End_I, End_J = calculate_riparto(jpjglo, jpiglo, nprocj, nproci, map_is_global)

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


def dump_outfile(tmask, choosen_procs, nproci, nprocj, map_is_global, filename="domdec.txt"):
    '''
    Arguments:
    * tmask         *
    * choosen_procs *
    * nproci        *
    * nprocj        *
    * map_is_global * boolean, indicates if the decomposition is performed on a global map
    * filename      *
    '''
    jpjglo, jpiglo = tmask.shape
    JPI, JPJ, Start_I, Start_J, End_I, End_J = calculate_riparto(jpjglo, jpiglo, nprocj, nproci, map_is_global)

    M,C = get_wp_matrix(tmask, nprocj, nproci, Start_I, End_I, Start_J, End_J)
    J,I = M.nonzero()    
    
    if map_is_global:
        WEST, EAST, NORTH, SOUTH, NBONDI,NBONDJ, NORTH_BND = neighbors_global(M, choosen_procs,nproci, nprocj)

        OUT = np.zeros((choosen_procs,14), dtype=np.int32)
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
            OUT[rank,13] = NORTH_BND[rank]
        np.savetxt(filename, OUT, fmt=14*"%5d")
    else:
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

if __name__ == "__main__":
    
    print("Testing calculate_riparto()")
    jpiglo = 1442
    jpjglo = 1021
    nproci = 5
    nprocj = 5
    map_is_global = True
    JPI, JPJ, Start_I, Start_J, End_I, End_J = calculate_riparto(jpjglo, jpiglo, nprocj, nproci, map_is_global)
    print("latitudinal decomposition")
    print(JPJ)
    print(Start_J)
    print(End_J)

    print("longitudinal decomposition")
    print(JPI)
    print(Start_I)
    print(End_I)

