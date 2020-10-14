MODULE NODES_MODULE
        !---------------------------------------------------------------------
        !
        !                       MODULE MATRIX_VARS
        !                     ******************
        !
        !
        !  gcoidess developing 
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
        
        !INTEGER IERROR
        !INTEGER :: k
        !INTEGER :: lengt
        !INTEGER :: i, j, p
        !CHARACTER*(MPI_MAX_PROCESSOR_NAME) buff_recv
        !CHARACTER (len=MPI_MAX_PROCESSOR_NAME),dimension(mpi_glcomm_size) :: total_array ! 
        !CHARACTER(len=:),allocatable :: total_array(:)        


        PUBLIC

        INTEGER :: nodes
        !CHARACTER*(MPI_MAX_PROCESSOR_NAME) local_array 
        
        INTEGER, allocatable, dimension (:) :: writing_procs
        !INTEGER :: lengt
        
        CONTAINS
!-------------------------------------------------------------------
        !SUBROUTINE NODES_MODULE_LEN()
        !INTEGER :: lengt
        !INTEGER :: IERROR
        !CALL MPI_GET_PROCESSOR_NAME(local_array, lengt, IERROR)
        !END SUBROUTINE NODES_MODULE_LEN
!------------------------------------------------------------------------------------------------
        !SUBROUTINE NODES_MODULE_DECLARATION()
        !INTEGER :: IERROR
        !CHARACTER (len = lengt), dimension(mpi_glcomm_size) :: total_array
        !END SUBROUTINE NODES_MODULE_DECLARATION
!--------------------------------------------------------------
        SUBROUTINE NODES_MODULE_FIND()
        
        INTEGER :: i, j, p, k
        INTEGER :: IERROR
        !INTEGER :: lengt
        CHARACTER*(MPI_MAX_PROCESSOR_NAME) buff_recv
        
        !CHARACTER(len=lengt),allocatable :: total_array(:)

        !CALL NODES_MODULE_LEN()
        CHARACTER (len = lengt), dimension(mpi_glcomm_size) :: total_array
        
        write (*,*) 'allocation rank',myrank, lengt, local_array

        call mppsync()

        CALL MPI_GATHER( local_array, lengt,MPI_CHAR, total_array,lengt,MPI_CHAR, 0, MPI_COMM_WORLD, IERROR)

       
        IF (myrank == 0) THEN
                nodes = 1
                p=1
                k=2
        
        !check print
                !DO i=0 , mpi_glcomm_size - 1
                !        write(*,*) 'RANK', i , total_array(i+1)
                !END DO
                                                                           
        !define how many strings (alias how many name nodes) are inside the total array
                
                DO i=2, mpi_glcomm_size
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
                        END IF
                END DO


        !print number of nodes

                write(*,*) 'Number of nodes is', nodes

        !rank 0 send to all ranks nodes

                DO i=1, mpi_glcomm_size - 1
                        CALL MPI_Send(nodes,1,MPI_INT,i,4,MPI_COMM_WORLD,IERROR)
                END DO

                write(*,*) 'nodes sent'


        !determing how many processor are inside first node, each node, delta numberby counting inside total array
        !creating array of writing nodes

                ALLOCATE (writing_procs(nodes))        
        
                writing_procs(1) = 0
                DO i=2, mpi_glcomm_size
                        IF (nodes==1) THEN
                                !if the number of node is one break, no sense to calculate, avoid problems
                                EXIT
                        
                        ELSE IF (total_array(p) /= total_array(i)) THEN
                                writing_procs(k)= i-1
                                p=i
                                k=k+1
                        ELSE
                                CYCLE
                        END IF
                END DO                

        !rank 0 send to all ranks writing_procs

                DO i=1, mpi_glcomm_size - 1
                        CALL MPI_Send(writing_procs,nodes,MPI_INT,i,3,MPI_COMM_WORLD,IERROR)
                END DO

        END IF

        ! broadcast number of nodes from rank 0 to all ranks
        
        !CALL MPI_Bcast(nodes,1,MPI_INT,0,MPI_COMM_WORLD, IERROR)

        ! all the ranks receive the array of writing procs
        
        IF (myrank >0) THEN
                
                CALL MPI_Recv(nodes,1,MPI_INT,0,4,MPI_COMM_WORLD,MPI_STATUS_IGNORE, IERROR)
                
                ALLOCATE (writing_procs(nodes))

                CALL MPI_Recv(writing_procs,nodes,MPI_INT,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE, IERROR)

                DO k=1, nodes
                        write (*,*) 'writing procs position is ', k, writing_procs(k)
                END DO
        END IF

!check number of nodes

        !write(*,*) 'Number of nodes is', nodes

        END SUBROUTINE NODES_MODULE_FIND

!---------------------------------------------------------------------------------------------------------

END MODULE
                


















                
