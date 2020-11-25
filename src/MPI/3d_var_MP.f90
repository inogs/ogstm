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
        
#ifdef ExecDA
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include "petsc/finclude/petscvec.h"
#else
#include <petsc/finclude/petscvecdef.h>
#endif
#endif

#ifdef ExecDA
        use DA_mem
        use mpi_str, only: Var3DCommunicator
        use petscvec, only: PETSC_COMM_WORLD, PETSC_NULL_CHARACTER
#endif

        USE calendar
        USE myalloc
        USE IO_mem
        USE FN_mem
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module
        
        USE NODES_MODULE

        IMPLICIT NONE
        


        PUBLIC
        
        LOGICAL :: V3D_VAR_PARALLEL        

        CONTAINS

!--------------------------------------------------------------
        SUBROUTINE DEFINE_3D_PARALLEL_PROCS()
        
        USE NODES_MODULE        
        
#ifdef ExecDA
        PetscErrorCode :: stat
#endif
        
        INTEGER, allocatable, dimension (:) :: TREd_procs_per_node_array
        INTEGER :: i, j, counter_3d_procs, IERROR
        INTEGER :: k
        INTEGER :: case_selection
        write(*,*) 'startinROUTINE'
        V3D_VAR_PARALLEL = .false.

!define how many procs per node to use for parallel 3d var
!cases:
! 1 node = max procs per 1 node -> 9
! + nodes but nodes*TREd_procs_per_node < DA_nprocs -> 5 procs per node
! + nodes with DA_n procs -> 5 procs per node + distribution        
#ifdef ExecDA
        !TREd_procs_per_node = 5
        !max_procs_per_one_node = 9
        
        if(nodes == 1) then
                case_selection = 1
        else if (nodes*TREd_procs_per_node < DA_Nprocs) then
                case_selection = 2
        else
                case_selection = 3
        end if

        counter_3d_procs = 1

        SELECT CASE (CASE_SELECTION)

                CASE(1)
                        ALLOCATE(TREd_procs_per_node_array(max_procs_per_one_node))
                        do k=1, max_procs_per_one_node
                                TREd_procs_per_node_array(k) = k - 1
                        end do
                        DO i=1, 9
                                write(*,*) '3d_var_proc is',TREd_procs_per_node_array(i)

                        END DO  
                        do i=1,max_procs_per_one_node
                                IF(MYRANK == TREd_procs_per_node_array(i)) then
                                        V3D_VAR_PARALLEL = .true.
                                END IF
                        end do
                CASE(2)
                        ALLOCATE(TREd_procs_per_node_array(nodes*TREd_procs_per_node)) 
                        do i=1, nodes
                                if(counter_3d_procs > nodes*TREd_procs_per_node) then
                                        exit
                                else
                                        TREd_procs_per_node_array(counter_3d_procs)=writing_procs(i)
                                        counter_3d_procs = counter_3d_procs + 1
                                        do j=1, TREd_procs_per_node -1
                                                if(counter_3d_procs > nodes*TREd_procs_per_node)then
                                                        exit
                                                else
                                                        TREd_procs_per_node_array(counter_3d_procs)=writing_procs(i)+ j
                                                        counter_3d_procs = counter_3d_procs + 1
                                                end if
                                        end do
                                end if
                        end do
                        DO i=1, counter_3d_procs - 1
                                write(*,*) '3d_var_procis',TREd_procs_per_node_array(i)
                        END DO
                        do i=1,nodes*TREd_procs_per_node
                                IF(MYRANK == TREd_procs_per_node_array(i)) then
                                        V3D_VAR_PARALLEL = .true.
                                END IF
                        end do
                CASE(3)
                        ALLOCATE (TREd_procs_per_node_array(DA_Nprocs))
                        DO i=1, NODES
                                if(counter_3d_procs > DA_Nprocs)then
                                        exit
                                else
                                        TREd_procs_per_node_array(counter_3d_procs)=writing_procs(i)
                                        counter_3d_procs = counter_3d_procs + 1
                                        do j=1, TREd_procs_per_node -1
                                                if(counter_3d_procs > DA_Nprocs) then
                                                        exit
                                                else
                                                        TREd_procs_per_node_array(counter_3d_procs)=writing_procs(i)+ j
                                                        counter_3d_procs = counter_3d_procs + 1
                                                end if
                                        end do
                                end if
                        end do
                        DO i=1, counter_3d_procs - 1
                                write(*,*)'3d_var_procis',TREd_procs_per_node_array(i)
                        END DO
                        do i=1,DA_Nprocs
                                IF(MYRANK == TREd_procs_per_node_array(i)) then
                                        V3D_VAR_PARALLEL = .true.
                                END IF
                        end do
                END SELECT

!PG parts from ogstm_mpi

        if(V3D_VAR_PARALLEL) then

                SELECT CASE (CASE_SELECTION)
                        CASE(1)
                                call MPI_Comm_split(MPI_COMM_WORLD, max_procs_per_one_node,myrank,Var3DCommunicator, ierror)
                        CASE(2)
                                call MPI_Comm_split(MPI_COMM_WORLD,nodes*TREd_procs_per_node,myrank,Var3DCommunicator, ierror)
                        CASE(3)
                                call MPI_Comm_split(MPI_COMM_WORLD,DA_nprocs,myrank,Var3DCommunicator,ierror)
                END SELECT

                PETSC_COMM_WORLD = Var3DCommunicator
                call PetscInitialize(PETSC_NULL_CHARACTER,stat)
                CHKERRQ(stat)
        else
                call MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED,myrank,Var3DCommunicator, ierror)
        endif
#endif
                
        END SUBROUTINE DEFINE_3D_PARALLEL_PROCS
!--------------------------------------------------------------
END MODULE TREd_var_MP
