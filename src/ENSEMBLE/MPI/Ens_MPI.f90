! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_MPI
    use mpi

implicit none
    integer :: EnsDebug
    integer :: EnsComm, EnsRank, EnsSize

contains
    subroutine Ens_MPI_nodes_find_and_split
        use NODES_MODULE, &
            ONLY: nodes, writing_procs, ind_col
        USE NODE_NAME, &
            ONLY: local_array
        USE ogstm_mpi_module, &
            ONLY: glcomm, glsize, glrank
        use myalloc, &
            only: lwp, mycomm, myrank, mysize
        
        INTEGER :: i, j, k
        INTEGER :: IERROR
        CHARACTER (len = MPI_MAX_PROCESSOR_NAME), dimension(glsize) :: total_array
        integer, dimension(glsize) :: nodes_array
        integer, dimension(:), allocatable :: ranks_per_node_array
        integer :: glcomm_ordered, nodecomm, noderank, writingcomm

        CALL MPI_GATHER(local_array, MPI_MAX_PROCESSOR_NAME,MPI_CHAR, total_array,MPI_MAX_PROCESSOR_NAME,MPI_CHAR, 0, glcomm, IERROR)
        
        !call mpi_comm_split_type(glcomm, SPLIT_TYPE, KEY, INFO, NEWCOMM, IERROR) Openmpi permette di usare il tipo OMPI_COMM_TYPE_SOCKET. Si potrebbe usarlo per individuare i socket a cui sono bindati i task.
        !Ma intelmpi pare non supporti la cosa e il tutto e' ancora mal documentato. Per il momento lasciamo perdere.

        IF (glrank == 0) THEN
            
            nodes = 0

            DO i=1, glsize
                DO j=i-1, 1, -1
                    IF (total_array(i) == total_array(j)) THEN
                        EXIT
                    END IF
                END DO
                IF (j==0) THEN
                    nodes = nodes + 1 
                    nodes_array(i)=nodes
                else
                    nodes_array(i)=nodes_array(j)
                END IF
            END DO

            write(*,*) 'Total number of nodes is', nodes
            
            if (EnsDebug>1) then
                do i=1,glsize
                    write(*,*) i,': ',total_array(i), ', ', nodes_array(i)
                end do
            end if
        
        end if
        
        call mpi_scatter(nodes_array, 1, MPI_INT, ind_col, 1, MPI_INT, 0, glcomm, ierror)
        
        call mpi_comm_split(glcomm, 0, ind_col, glcomm_ordered, ierror)
        glcomm = glcomm_ordered
        CALL mpi_comm_rank(glcomm,glrank,ierror)
        
        if (mod(glsize,EnsSize)/=0) then
            write(*,*) "The total number of mpi tasks (glsize=,", glsize, ") is not multiple of the ensemble size (EnsSize=", EnsSize, "). Aborting."
        end if
        
        call MPI_Comm_split(glcomm, glrank/(glsize/EnsSize), glrank, mycomm, ierror)
        CALL mpi_comm_rank(mycomm,myrank,ierror)
        CALL mpi_comm_size(mycomm,mysize,ierror)
        
        call MPI_Comm_split(glcomm, myrank, glrank, EnsComm, ierror)
        CALL mpi_comm_rank(EnsComm,EnsRank,ierror)
        !CALL mpi_comm_size(EnsComm,EnsSize,ierror)
        
        
        call mpi_comm_split(mycomm, ind_col, myrank, nodecomm, ierror)
        CALL mpi_comm_rank(nodecomm,noderank,ierror)
        call mpi_comm_split(mycomm, noderank, myrank, writingcomm, ierror)
        
        CALL mpi_comm_rank(writingcomm,ind_col,ierror)
        ind_col=ind_col+1
        CALL MPI_Bcast(ind_col,1,MPI_INT,0,nodecomm, IERROR)
        
        CALL mpi_comm_size(writingcomm,nodes,ierror)
        CALL MPI_Bcast(nodes,1,MPI_INT,0,nodecomm, IERROR)
        
        ALLOCATE (writing_procs(nodes))
        
        if (noderank==0) CALL mpi_allgather(myrank, 1,MPI_INT, writing_procs,1,MPI_INT, writingcomm, IERROR)
        CALL MPI_Bcast(writing_procs,nodes,MPI_INT,0,nodecomm, IERROR)
        
        call mpi_comm_free(nodecomm, IERROR)
        call mpi_comm_free(writingcomm, IERROR)
        
        if (EnsDebug>1) then
            do i=1,glsize
                CALL mpi_barrier(glcomm,ierror)
                if (glrank==i) then
                    write(*,*) 'gl = ',i, ', my = ',myrank ,', ind_col = ',ind_col
                    write(*,*) 'writing_procs = ', writing_procs
                end if
            end do
        end if
    
    end subroutine
    

end module
