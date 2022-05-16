! developed by Simone Spada (sspada@ogs.it) at OGS

module Ens_MPI
    use mpi
    use Ens_Mem, &
        only: EnsDebug, EnsComm, EnsRank, EnsSize

implicit none

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
            
        integer, parameter :: procs_distribution_type=2 !1 is not supported at this moment
            ! 1 (not supported): faster communications inside ogstm communicator (mycomm)
            ! 2 (supported if ranks per node is multiple of EnsSize): faster communications between ensemble members (EnsComm). Less RAM per per node needed.
        
        INTEGER :: i, j, k
        INTEGER :: IERROR
        CHARACTER (len = MPI_MAX_PROCESSOR_NAME), dimension(glsize) :: total_array
        integer, dimension(glsize) :: nodes_array
        integer, dimension(:), allocatable :: ranks_per_node_array
        integer :: glcomm_ordered, nodecomm, noderank, nodesize, writingcomm

        CALL MPI_GATHER(local_array, MPI_MAX_PROCESSOR_NAME,MPI_CHAR, total_array,MPI_MAX_PROCESSOR_NAME,MPI_CHAR, 0, glcomm, IERROR)
        
        !call mpi_comm_split_type(glcomm, SPLIT_TYPE, KEY, INFO, NEWCOMM, IERROR) Openmpi permette di usare il tipo OMPI_COMM_TYPE_SOCKET. Si potrebbe usarlo per individuare i socket a cui sono bindati i task.
        !Ma intelmpi pare non supporti la cosa e il tutto e' ancora mal documentato. Per il momento lasciamo perdere.

        ! counting nodes:
        
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
            
            if (EnsDebug>0) then
                do i=1,glsize
                    write(*,*) 'proc ', i,': node ', nodes_array(i)
                end do
            end if
        
        end if
        
        call mpi_scatter(nodes_array, 1, mpi_integer, ind_col, 1, mpi_integer, 0, glcomm, ierror)
        
        ! ordering:
        
        call mpi_comm_split(glcomm, 0, ind_col, glcomm_ordered, ierror)
        call mpi_comm_free(glcomm, ierror)
        glcomm = glcomm_ordered
        CALL mpi_comm_rank(glcomm,glrank,ierror)
        
        ! creating ensemble communicators:
        
        if (mod(glsize,EnsSize)/=0) then
            write(*,*) "The total number of mpi tasks (glsize=,", glsize, ") is not multiple of the ensemble size (EnsSize=", EnsSize, "). Aborting."
            call MPI_abort(glcomm, 1, ierror)
        end if
        
        SELECT CASE (procs_distribution_type)
            CASE(1)
                call MPI_Comm_split(glcomm, glrank/(glsize/EnsSize), glrank, mycomm, ierror)
                CALL mpi_comm_rank(mycomm,myrank,ierror)
                CALL mpi_comm_size(mycomm,mysize,ierror)
                
                call MPI_Comm_split(glcomm, myrank, glrank, EnsComm, ierror)
                CALL mpi_comm_rank(EnsComm,EnsRank,ierror)
                !CALL mpi_comm_size(EnsComm,EnsSize,ierror)
                
            CASE(2)
                call MPI_Comm_split(glcomm, glrank/EnsSize, glrank, EnsComm, ierror)
                CALL mpi_comm_rank(EnsComm,EnsRank,ierror)
                !CALL mpi_comm_size(EnsComm,EnsSize,ierror)
                
                call MPI_Comm_split(glcomm, EnsRank, glrank, mycomm, ierror)
                CALL mpi_comm_rank(mycomm,myrank,ierror)
                CALL mpi_comm_size(mycomm,mysize,ierror)   
                
            CASE DEFAULT
                write(*,*) 'Wrong procs_distribution_type value. Aborting.'
                call MPI_abort(glcomm, 1, ierror)
        END SELECT
        
        !checks
        call mpi_barrier(glcomm,ierror)
        call mpi_comm_split(EnsComm, ind_col, EnsRank, nodecomm, ierror)
        call mpi_comm_size(nodecomm, nodesize, ierror)
        if (nodesize/=EnsSize) then
            write(*,*) 'The number of task per node is not a multiple of EnsSize. This is not supported at this moment. Aborting.'
            write(*,*) 'glrank: ', glrank, ', node number: ', ind_col, ', Ens-procs in this node: ', nodesize
            call MPI_abort(glcomm, 1, ierror)
        end if
        call mpi_comm_free(nodecomm, ierror)
        
call mpi_barrier(glcomm,ierror)
call mpi_comm_split(glComm, ind_col, glRank, nodecomm, ierror)
CALL mpi_comm_rank(nodecomm, noderank,ierror)  
call mpi_comm_size(nodecomm, nodesize, ierror)
if (noderank==0 .and. nodesize>48) then
    write(*,*) 'WARNING!!! node ', ind_col, ' has ', nodesize, ' procs.'
    !call mpi_abort(glcomm, 1, ierror)
end if
call mpi_comm_free(nodecomm, ierror)
        
        ! computing ind_col and writing_procs:
        
        call mpi_comm_split(mycomm, ind_col, myrank, nodecomm, ierror)
        CALL mpi_comm_rank(nodecomm,noderank,ierror)            
        
        call mpi_comm_split(mycomm, noderank, myrank, writingcomm, ierror)
        
        CALL mpi_comm_rank(writingcomm,ind_col,ierror)
        ind_col=ind_col+1
        CALL MPI_Bcast(ind_col,1,mpi_integer,0,nodecomm, IERROR)
        
        CALL mpi_comm_size(writingcomm,nodes,ierror)
        CALL MPI_Bcast(nodes,1,mpi_integer,0,nodecomm, IERROR)
        
        ALLOCATE (writing_procs(nodes))
        
        if (noderank==0) CALL mpi_allgather(myrank, 1,mpi_integer, writing_procs,1,mpi_integer, writingcomm, IERROR)
        CALL MPI_Bcast(writing_procs,nodes,mpi_integer,0,nodecomm, IERROR)
        
        call mpi_comm_free(nodecomm, IERROR)
        call mpi_comm_free(writingcomm, IERROR)
        
        if (EnsDebug>1) then
            do i=0,glsize-1
                CALL mpi_barrier(glcomm,ierror)
                if (glrank==i) then
                    write(*,*) 'gl = ',i, ', my = ',myrank ,', ind_col = ',ind_col
                    write(*,*) 'writing_procs = ', writing_procs
                end if
            end do
            !CALL mpi_barrier(glcomm,ierror)
            !call MPI_abort(glcomm, 1, ierror)
        end if
    
    end subroutine
    
    subroutine Ens_MPI_Finalize
        use myalloc, &
            only: mycomm
            
        integer ierror
        
        call mpi_comm_free(mycomm, ierror)
        call mpi_comm_free(EnsComm, ierror)
        
    end subroutine

end module
