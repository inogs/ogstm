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
        

        !V3D_VAR_PARALLEL = .false.
        CONTAINS

!-------------------------------------------------------------
!        SUBROUTINE ALLOCATE_3D_PARALLEL_PROCS()
!        
!        USE DEFINE_3D_PARALLEL_PROCS
!        USE NODES_MODULE
!
!        ALLOCATE (TREd_procs_per_node_array(DA_Nprocs))
!
!        END SUBROUTINE ALLOCATE_3D_PARALLEL_PROCS
!--------------------------------------------------------------
        SUBROUTINE DEFINE_3D_PARALLEL_PROCS()
        
        USE NODES_MODULE        
        
#ifdef ExecDA
        PetscErrorCode :: stat
#endif
        INTEGER :: TREd_procs_per_node
        INTEGER, allocatable, dimension (:) :: TREd_procs_per_node_array
        INTEGER :: i, j, counter_3d_procs, IERROR
        INTEGER :: k, max_procs_per_one_node
        write(*,*) 'startinROUTINE'
        V3D_VAR_PARALLEL = .false.

        TREd_procs_per_node = 5
        max_procs_per_one_node = 9
        
        if (nodes == 1) then
                ALLOCATE (TREd_procs_per_node_array(max_procs_per_one_node))
        else
                if (nodes*TREd_procs_per_node < DA_Nprocs)then
                        ALLOCATE(TREd_procs_per_node_array(nodes*TREd_procs_per_node))
                else
                        ALLOCATE (TREd_procs_per_node_array(DA_Nprocs))
                end if
        end if
        !CALL ALLOCATE_3D_PARALLEL_PROCS()

        ! define which processors will be used for 3d var scheme

        counter_3d_procs = 1
        write(*,*) 'first_write 3dvarparallel', counter_3d_procs, myrank 
                
        if (nodes == 1)then
                do k=1, max_procs_per_one_node
                        TREd_procs_per_node_array(k) = k - 1
                end do
        else if (nodes*TREd_procs_per_node < DA_Nprocs) then
                if(counter_3d_procs > nodes*TREd_procs_per_node) then
                        exit
                else
                        TREd_procs_per_node_array(counter_3d_procs)=writing_procs(i)
                        counter_3d_procs = counter_3d_procs + 1 
                        do j=1, TREd_procs_per_node -1
                                if(counter_3d_procs > nodes*TREd_procs_per_node) then
                                        exit
                                else
                                        TREd_procs_per_node_array(counter_3d_procs)=writing_procs(i) + j
                                        write(*,*) 'fourth_write3dvarparallel',counter_3d_procs,'',writing_procs(i)+j,TREd_procs_per_node_array(counter_3d_procs)
                                        counter_3d_procs = counter_3d_procs + 1
                                        write(*,*) 'third_write3dvarparallel',counter_3d_procs, myrank
                                end if
                        end do
                end if
        else
                TREd_procs_per_node_array(counter_3d_procs)= writing_procs(i)
                        write(*,*) 'second_write3dvarparallel',counter_3d_procs,'',writing_procs(i),TREd_procs_per_node_array(counter_3d_procs)
                        counter_3d_procs = counter_3d_procs + 1
                        write(*,*) 'third_write 3dvarparallel',counter_3d_procs,myrank
                        do j=1, TREd_procs_per_node -1
                                if(counter_3d_procs > DA_Nprocs) then
                                        exit
                                else
                                        TREd_procs_per_node_array(counter_3d_procs)=writing_procs(i) + j
                                        write(*,*) 'fourth_write3dvarparallel',counter_3d_procs,'',writing_procs(i)+j,TREd_procs_per_node_array(counter_3d_procs)
                                        counter_3d_procs = counter_3d_procs + 1
                                        write(*,*) 'third_write3dvarparallel',counter_3d_procs, myrank
                                end if
                        end do
        end if

!                write(*,*) "before send"
!                ! proc 0 send the array to all processors
!                DO i=1, mpi_glcomm_size - 1
!                        CALL MPI_Send(TREd_procs_per_node_array,DA_Nprocs,MPI_INT,i,5,MPI_COMM_WORLD,IERROR)
!                END DO
!                write(*,*) 'after send'
!                !test
        if(myrank==0)then        
                DO i=1, 9

                        write(*,*) '3d_var_proc is',TREd_procs_per_node_array (i)

                END DO
        end if

        
        !all procs =! 0 receive the array
!        IF(myrank > 0)then
!
!                CALL MPI_Recv(TREd_procs_per_node_array,DA_Nprocs,MPI_INT,0,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE, IERROR)
!                
!        END IF

        !all processors change the boolean if they are linked to 3d var scheme
        if (nodes ==1) then
                do i=1,max_procs_per_one_node
                        IF(MYRANK == TREd_procs_per_node_array(i)) then
                                V3D_VAR_PARALLEL = .true.
                        END IF
                end do
        else if (nodes*TREd_procs_per_node < DA_Nprocs) then
                do i=1,nodes*TREd_procs_per_node
                        IF(MYRANK == TREd_procs_per_node_array(i)) then
                                V3D_VAR_PARALLEL = .true.
                        END IF
                end do
        else
                do i=1,DA_Nprocs
                        IF(MYRANK == TREd_procs_per_node_array(i)) then
                                V3D_VAR_PARALLEL = .true.
                        END IF
                end do
        end if


!PG parts from ogstm_mpi

#ifdef ExecDA
        if(V3D_VAR_PARALLEL) then
                if(nodes ==1) then
                        call MPI_Comm_split(MPI_COMM_WORLD, max_procs_per_one_node,myrank,Var3DCommunicator, ierror)
                else if (nodes*TREd_procs_per_node <DA_Nprocs) then
                        call MPI_Comm_split(MPI_COMM_WORLD,nodes*TREd_procs_per_node,myrank,Var3DCommunicator, ierror)
                else
                        call MPI_Comm_split(MPI_COMM_WORLD,DA_nprocs,myrank,Var3DCommunicator,ierror)
                end if
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
                


















                
