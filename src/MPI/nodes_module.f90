MODULE NODES_MODULE
        !---------------------------------------------------------------------
        !
        !                       MODULE MATRIX_VARS
        !                     ******************
        !
        !
        !  gcoidessa@inogs.it developing 
        !  Purpose :
        !  ---------
        !  OBTAIN ARRAY OF NODES AND WRITING PROCESSOR
        !  ---------
        !  Subroutines: 
        !       -NODES_FIND()
        !              
        !


        USE calendar
        USE myalloc
        USE IO_mem
        USE FN_mem
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module
        USE NODE_NAME

        IMPLICIT NONE
        


        PUBLIC

        INTEGER :: nodes
        integer :: ind_col
        
        INTEGER, allocatable, dimension (:) :: writing_procs
        
        CONTAINS

!--------------------------------------------------------------
        SUBROUTINE NODES_MODULE_FIND()
        
        INTEGER :: i, j, p, k
        INTEGER :: IERROR        

        CHARACTER (len = MPI_MAX_PROCESSOR_NAME), dimension(mysize) :: total_array
        integer, dimension(mysize) :: nodes_array
        
!        write (*,*) 'allocation rank',myrank, lengt, local_array

        call mppsync()

        CALL MPI_GATHER( local_array, MPI_MAX_PROCESSOR_NAME,MPI_CHAR, total_array,MPI_MAX_PROCESSOR_NAME,MPI_CHAR, 0, mycomm, IERROR)

       
        IF (myrank == 0) THEN
                nodes = 1
                p=1
                k=2
                
                nodes_array(1)=1
                DO i=2, mysize
                        IF (i==1) THEN
                               ! write(*,*)
                        END IF
                        DO j=1, i
                                IF (total_array(i) == total_array(j)) THEN
                                        EXIT
                                END IF
                        END DO
                        IF (i==j) THEN
                                nodes = nodes + 1 
                                nodes_array(i)=nodes
                        else
                                nodes_array(i)=nodes_array(j)
                        END IF
                END DO


        !print number of nodes

                if (lwp) write(*,*) 'Number of nodes is', nodes

        !rank 0 send to all ranks nodes

                DO i=1, mysize - 1
                        CALL MPI_Send(nodes,1,mpi_integer,i,4,mycomm,IERROR)
                END DO

                if (lwp) write(*,*) 'nodes sent'


        !determing how many processor are inside first node, each node, delta numberby counting inside total array
        !creating array of writing nodes

                ALLOCATE (writing_procs(nodes))        
        
                writing_procs(1) = 0
                DO i=2, mysize
                        IF (nodes==1) THEN
                                !if the number of node is one break, no sense to calculate, avoid problems
                                EXIT
                        
                        ELSE IF (nodes_array(p) < nodes_array(i)) THEN
                                writing_procs(k)= i-1
                                p=i
                                k=k+1
                        ELSE
                                CYCLE
                        END IF
                END DO                

        !rank 0 send to all ranks writing_procs

                DO i=1, mysize - 1
                        CALL MPI_Send(writing_procs,nodes,mpi_integer,i,3,mycomm,IERROR)
                END DO

        END IF

        ! broadcast number of nodes from rank 0 to all ranks
        ! all the ranks receive the array of writing procs
        
        IF (myrank >0) THEN
                
                CALL MPI_Recv(nodes,1,mpi_integer,0,4,mycomm,MPI_STATUS_IGNORE, IERROR)
                
                ALLOCATE (writing_procs(nodes))

                CALL MPI_Recv(writing_procs,nodes,mpi_integer,0,3,mycomm,MPI_STATUS_IGNORE, IERROR)

!                DO k=1, nodes
!                        write (*,*) 'writing procs position is ', k, writing_procs(k)
!                END DO
        END IF
        
        call mpi_scatter(nodes_array, 1, mpi_integer, ind_col, 1, mpi_integer, 0, mycomm, ierror)

!check number of nodes

        !write(*,*) 'Number of nodes is', nodes

        END SUBROUTINE NODES_MODULE_FIND

!---------------------------------------------------------------------------------------------------------

END MODULE
                


















                
