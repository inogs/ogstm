MODULE TREd_var_MP 
        !---------------------------------------------------------------------
        !
        !                       MODULE MATRIX_VARS
        !                     ******************
        !
        !
        !  gcoidessa@inogs.it developing 
        !  Purpose :
        !  ---------
        !  OBTAIN ARRAY PROCESSOR FOR 3D VAR PARALLELISM ON NODES
        !  ---------
        !  Subroutines: 
        !       - ALLOCATE_3D_PARALLEL_PROCS()
        !       - DEFINE_3D_PARALLEL_PROCS()
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
        
        USE NODES_MODULE

        IMPLICIT NONE
        


        PUBLIC
        
        LOGICAL :: V3D_VAR_PARALLEL
        

        V3D_VAR_PARALLEL = .false.
        CONTAINS

!-------------------------------------------------------------
        SUBROUTINE ALLOCATE_3D_PARALLEL_PROCS()
        
        USE DEFINE_3D_PARALLEL_PROCS
        USE NODES_MODULE

        ALLOCATE (3d_procs_per_node_array(nodes * 3d_procs_per_node))

        END SUBROUTINE ALLOCATE_3D_PARALLEL_PROCS
!--------------------------------------------------------------
        SUBROUTINE DEFINE_3D_PARALLEL_PROCS()
        
        USE NODES_MODULE
        
        INTEGER :: 3d_procs_per_node
        INTEGER, allocatable, dimension (:) :: 3d_procs_per_node_array
        INTEGER :: i, j, counter_3d_procs, IERROR

        3d_procs_per_node = 5

        CALL ALLOCATE_3D_PARALLEL_PROCS()
        
        ! only in rank 0
        ! define which processors will be used for 3d var scheme
        IF(myrank == 0)then
                counter_3d_procs = 1
                do i=1, nodes
                        3d_procs_per_node_array(counter_3d_procs)= writing_procs(i)
                        do j=1, 3d_procs_per_node -1 
                                counter_3d_procs = counter_3d_procs + 1
                                3d_procs_per_node_array(counter_3d_procs)=writing_procs(i) + j
                        end do
                end do
                
                ! proc 0 send the array to all processors
                DO i=1, mpi_glcomm_size - 1
                        CALL MPI_Send(3d_procs_per_node_array,nodes * 3d_procs_per_node,MPI_INT,i,5,MPI_COMM_WORLD,IERROR)
                END DO

                !test
                DO i=1, nodes * 3d_procs_per_node

                        write(*,*) '3d_var_proc is',3d_procs_per_node_array (i)

                END DO

        END IF
        
        !all procs =! 0 receive the array
        IF(myrank > 0)then

                CALL MPI_Recv(3d_procs_per_node_array,nodes * 3d_procs_per_node,MPI_INT,0,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE, IERROR)
                
        END IF

        !all processors change the boolean if they are linked to 3d var scheme
        DO i=1, nodes * 3d_procs_per_node

                IF(MYRANK == 3d_procs_per_node_array(i) then
                        
                        V3D_VAR_PARALLEL = .true.
                END IF

        END DO

                
        END SUBROUTINE DEFINE_3D_PARALLEL_PROCS
!--------------------------------------------------------------
END MODULE TREd_var_MP
                


















                
