MODULE MPI_GATHER_INFO 
        !---------------------------------------------------------------------
        !
        !                       MODULE MPI_GATHER
        !                     ******************
        !
        !
        !  gcoidessa@inogs.it developing 
        !  Purpose :
        !  ---------
        !     Initialise indices for IO,
        !          subroutines: 
        !               - ALLOCATE_MPI_GATHER_INFO()
        !               - INIT_MPI_GATHER_INFO()
        !               - CLEAN_MEMORY_MPI_GATHER_INFO()
        !               - CLEAN_MEMORY_INIT_MPI_GATHERV_COUNTS_INFO()



        USE calendar
        USE myalloc
        USE IO_mem
        USE FN_mem
        USE TIME_MANAGER
        use mpi
        USE ogstm_mpi_module

        USE dtype_procs_string_module        
        USE nodes_module
        USE MATRIX_VARS

        IMPLICIT NONE
        INTEGER IERROR
        INTEGER ::wr_procs

        PUBLIC
        
        !timing

        DOUBLE PRECISION :: communication_start_time_gather_info, communication_finish_time_gather_info
        DOUBLE PRECISION :: communication_proctime_time_gather_info, communication_max_time_gather_info

        !variable do loops

        INTEGER :: loop_ind,loop_ind_2d
        INTEGER :: cont,cont_2d
        INTEGER :: sendcount,sendcount_2d

        LOGICAL :: WRITING_RANK_WR
        !allocatable

        INTEGER, allocatable, dimension (:) :: jpi_rec_a
        INTEGER, allocatable, dimension (:) :: jpj_rec_a
        INTEGER, allocatable, dimension (:) :: istart_a
        INTEGER, allocatable, dimension (:) :: jstart_a
        INTEGER, allocatable, dimension (:) :: iPe_a
        INTEGER, allocatable, dimension (:) :: jPe_a
        INTEGER, allocatable, dimension (:) :: iPd_a
        INTEGER, allocatable, dimension (:) :: jPd_a
        DOUBLE PRECISION, allocatable, dimension(:) :: bufftrn_TOT
        DOUBLE PRECISION, allocatable, dimension(:) :: buffDIA2d_TOT
        DOUBLE PRECISION, allocatable, dimension(:) :: buffDIA_TOT        
        DOUBLE PRECISION, allocatable, dimension(:) :: buffPHYS2d_TOT
        DOUBLE PRECISION, allocatable, dimension(:) :: buffPHYS_TOT

        INTEGER, allocatable, dimension (:) :: jprcv_count,jprcv_count_2d
        INTEGER, allocatable, dimension (:) :: jpdispl_count,jpdispl_count_2d


        double precision, allocatable ::  tottrnIO(:,:,:)


        CHARACTER(len=2) :: n_ranks_per_node_char
        INTEGER :: n_ranks_per_node
        !WRITING_RANK = .FALSE.

        CONTAINS

!------------------------------------------------------------
        SUBROUTINE ALLOCATE_MPI_GATHER_INFO()



        ALLOCATE (jpi_rec_a(mpi_glcomm_size))
        ALLOCATE (jpj_rec_a(mpi_glcomm_size))
        ALLOCATE (istart_a(mpi_glcomm_size))
        ALLOCATE (jstart_a(mpi_glcomm_size))
        ALLOCATE (iPe_a(mpi_glcomm_size))
        ALLOCATE (jPe_a(mpi_glcomm_size))
        ALLOCATE (iPd_a(mpi_glcomm_size))
        ALLOCATE (jPd_a(mpi_glcomm_size))
        ALLOCATE (bufftrn_TOT(jpi_max* jpj_max* jpk* mpi_glcomm_size))
        ALLOCATE (buffDIA2d_TOT(jpi_max* jpj_max*mpi_glcomm_size))
        ALLOCATE (buffDIA_TOT(jpi_max* jpj_max* jpk*mpi_glcomm_size)) 
        ALLOCATE (buffPHYS2d_TOT(jpi_max* jpj_max*mpi_glcomm_size))
        ALLOCATE (buffPHYS_TOT(jpi_max* jpj_max* jpk*mpi_glcomm_size))

        ALLOCATE (jprcv_count(mpi_glcomm_size))
        ALLOCATE (jpdispl_count(mpi_glcomm_size))
        ALLOCATE (jprcv_count_2d(mpi_glcomm_size))
        ALLOCATE (jpdispl_count_2d(mpi_glcomm_size))


        DO wr_procs=1, nodes

                if (myrank==writing_procs(wr_procs))then
                        allocate(tottrnIO(jpk,jpjglo,jpiglo)) 
                        tottrnIO  = huge(tottrnIO(1,1,1)) 
                        allocate(totglamt(jpjglo,jpiglo))
                        totglamt = huge(totglamt(1,1))  
                        bufftrn_TOT = huge(bufftrn_TOT(1)) 

                        WRITING_RANK_WR = .TRUE.     
#ifdef ExecDA                        
                        allocate(tottrnDA(jpiglo, jpjglo, jpk))
                        tottrnDA = huge(tottrnDA(1,1,1))
#endif
                end if
        end do

        DO wr_procs=2, nodes

                if (myrank==writing_procs(wr_procs))then
                        
                        allocate(tottrn(jpk, jpjglo, jpiglo))
                        tottrn = huge(tottrn(1,1,1))
                        allocate(tottrb(jpk, jpjglo, jpiglo))
                        tottrb = huge(tottrb(1,1,1))
                        !allocate(tottrnIO(jpk,jpjglo,jpiglo)) 
                        !tottrnIO  = huge(tottrnIO(1,1,1)) 
                        allocate(tottrbIO(jpk,jpjglo,jpiglo)) 
                        tottrbIO  = huge(tottrbIO(1,1,1)) 
                        allocate(totsnIO (jpk,jpjglo,jpiglo)) 
                        totsnIO   = huge(totsnIO(1,1,1))  
                        allocate(tottnIO (jpk,jpjglo,jpiglo)) 
                        tottnIO   = huge(tottnIO(1,1,1))  
                        allocate(totvatmIO(jpjglo,jpiglo))    
                        totvatmIO = huge(totvatmIO(1,1))  
                        allocate(totempIO(jpjglo,jpiglo))     
                        totempIO  = huge(totempIO(1,1))   
                        allocate(totqsrIO(jpjglo,jpiglo))     
                        totqsrIO  = huge(totqsrIO(1,1))   
                        allocate(totunIO(jpk,jpjglo,jpiglo))  
                        totunIO   = huge(totunIO(1,1,1))  
                        allocate(totvnIO(jpk,jpjglo,jpiglo))  
                        totvnIO   = huge(totvnIO(1,1,1))
                        allocate(totwnIO(jpk,jpjglo,jpiglo))  
                        totwnIO   = huge(totwnIO(1,1,1))  
                        allocate(totavtIO(jpk,jpjglo,jpiglo))
                        totavtIO  = huge(totavtIO(1,1,1))
                        allocate(tote3tIO(jpk,jpjglo,jpiglo))
                        tote3tIO  = huge(tote3tIO(1,1,1))
                        allocate(tottmaIO(jpk,jpjglo,jpiglo))
                        tottmaIO  = huge(tottmaIO(1,1,1))
                        allocate(tottrnIO2d(jpjglo,jpiglo))
                        tottrnIO2d= huge(tottrnIO2d(1,1))
                        !allocate(totglamt(jpjglo,jpiglo)) 
                        !totglamt = huge(totglamt(1,1))
                        allocate(totgphit(jpjglo,jpiglo))
                        totgphit = huge(totgphit(1,1))
                end if
        end do

        END SUBROUTINE ALLOCATE_MPI_GATHER_INFO
!------------------------------------------------------------

        SUBROUTINE INIT_MPI_GATHER_INFO()

        WRITING_RANK_WR = .FALSE.

        call getenv('RANKS_PER_NODE', n_ranks_per_node_char)
        if (n_ranks_per_node_char .eq. '') then
                write(*,*) 'ERROR: RANKS_PER_NODE environment variable not defined'
                write(*,*) 'EXAMPLE: export RANKS_PER_NODE=48'
                stop
        end if

        read(n_ranks_per_node_char , *) n_ranks_per_node

        CALL ALLOCATE_MPI_GATHER_INFO()
        !gather(send+recv from each rank, stored in array of each indices)

        call mppsync()

        communication_start_time_gather_info= MPI_Wtime()

        DO wr_procs=1, nodes



                CALL MPI_GATHER( jpi, 1, MPI_INTEGER, jpi_rec_a,1,MPI_INTEGER, writing_procs(wr_procs), MPI_COMM_WORLD, IERROR)
                CALL MPI_GATHER( jpj, 1, MPI_INTEGER, jpj_rec_a, 1,MPI_INTEGER, writing_procs(wr_procs), MPI_COMM_WORLD, IERROR)
                CALL MPI_GATHER( nimpp, 1, MPI_INTEGER, istart_a, 1,MPI_INTEGER, writing_procs(wr_procs), MPI_COMM_WORLD, IERROR)
                CALL MPI_GATHER( njmpp, 1, MPI_INTEGER, jstart_a, 1,MPI_INTEGER, writing_procs(wr_procs), MPI_COMM_WORLD, IERROR)
                CALL MPI_GATHER( nlei, 1, MPI_INTEGER, iPe_a, 1, MPI_INTEGER,writing_procs(wr_procs), MPI_COMM_WORLD, IERROR)
                CALL MPI_GATHER( nlej, 1, MPI_INTEGER, jPe_a, 1, MPI_INTEGER,writing_procs(wr_procs), MPI_COMM_WORLD, IERROR)
                CALL MPI_GATHER( nldi, 1, MPI_INTEGER, iPd_a, 1, MPI_INTEGER,writing_procs(wr_procs), MPI_COMM_WORLD, IERROR)
                CALL MPI_GATHER( nldj, 1, MPI_INTEGER, jPd_a, 1, MPI_INTEGER,writing_procs(wr_procs), MPI_COMM_WORLD, IERROR)
        
        END DO        
        communication_finish_time_gather_info= MPI_Wtime()
        communication_proctime_time_gather_info= communication_finish_time_gather_info - communication_start_time_gather_info
        
        !CALL MPI_Reduce( communication_proctime_time_gather_info, communication_max_time_gather_info, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, IERROR)

        !if(myrank == 0) then
        !        write(*,*) 'Communication time gather info is', communication_max_time_gather_info
        !end if
        
  
!        END SUBROUTINE INIT_MPI_GATHER_INFO

!---------------------------------------------------
!        SUBROUTINE INIT_MPI_GATHERV_COUNTS_INFO()
                
        sendcount = jpi * jpj * jpk
        sendcount_2d = jpi * jpj

        if(WRITING_RANK_WR)then
        
                cont = 0
                DO loop_ind = 1, mpi_glcomm_size
                        jprcv_count(loop_ind) = jpi_rec_a(loop_ind) * jpj_rec_a(loop_ind) * jpk
                        jpdispl_count(loop_ind) = cont
                        cont = cont + jprcv_count(loop_ind)        
                        !write(*,*) 'do loop',loop_ind,jprcv_count(loop_ind),jpdispl_count(loop_ind),cont
                end DO        
                !write(*,*) 'do loop finished'

                cont_2d = 0
                DO loop_ind_2d = 1, mpi_glcomm_size
                        jprcv_count_2d(loop_ind_2d) = jpi_rec_a(loop_ind_2d) * jpj_rec_a(loop_ind_2d) 
                        jpdispl_count_2d(loop_ind_2d) = cont_2d
                        cont_2d = cont_2d + jprcv_count_2d(loop_ind_2d)
                        !write(*,*) 'do loop',loop_ind_2d,jprcv_count_2d(loop_ind_2d),jpdispl_count_2d(loop_ind_2d),cont_2d
                end DO
                !write(*,*) 'do loop2d finished'
        end if
        !END DO


END SUBROUTINE INIT_MPI_GATHER_INFO


!--------------------------------------------------
        SUBROUTINE CLEAN_MEMORY_MPI_GATHER_INFO()

!to check

        DEALLOCATE (jpi_rec_a)
        DEALLOCATE (jpj_rec_a)
        DEALLOCATE (istart_a)
        DEALLOCATE (jstart_a)
        DEALLOCATE (iPe_a)
        DEALLOCATE (jPe_a)
        DEALLOCATE (iPd_a)
        DEALLOCATE (jPd_a)

        deallocate(tottrnIO)
        deallocate(totglamt)
        if (writing_rank_wr)then
                deallocate(buffDIA2d_TOT)
                deallocate(bufftrn_TOT)
                deallocate(buffDIA_TOT)
                deallocate(buffPHYS2d_TOT)
                deallocate(buffPHYS_TOT)
        end if
        END SUBROUTINE CLEAN_MEMORY_MPI_GATHER_INFO


!--------------------------------------------------
        SUBROUTINE CLEAN_MEMORY_INIT_MPI_GATHERV_COUNTS_INFO()

        DEALLOCATE (jprcv_count)
        DEALLOCATE (jpdispl_count)

        END SUBROUTINE CLEAN_MEMORY_INIT_MPI_GATHERV_COUNTS_INFO


END MODULE
