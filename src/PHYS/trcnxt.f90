      SUBROUTINE trcnxt
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcnxt
!!!                     *******************
!!!
!!!  PURPOSE :
!!!  ---------
!!!    compute the passive tracers fields at the next time
!!!    step from their temporal trends
!!!
!!    METHOD :
!!
!!
!!      Apply lateral boundary conditions (nperio=0,closed  nperio=1,
!!      east-west cyclic nperio=2, symmetric round the equator)
!!      on tra arrays
!!
!!   default:
!!      arrays swap
!!         (trn) = (tra)   (tra) = (0,0)
!!         (trb) = (trn) 
!!
!!   For Arakawa Sheme : IF key_trc_arakawa defined
!!      A Asselin time filter applied on now tracers (trn) to avoid
!!      the divergence of two consecutive time-steps and tr arrays
!!      to prepare the next time_step:
!!         (trb) = (trn) + gamma [ (trb) + (tra) - 2 (trn) ]
!!         (trn) = (tra)   (tra) = (0,0)
!!
!!      array swap for tracers to start the next time step
!!
!!
!!      macrotasking on tracer slab
!!

       USE myalloc
       USE FN_mem
       USE Time_Manager
       use mpi
       USE ogstm_mpi_module

      IMPLICIT NONE

!! local declarations
!! ==================
      
! omp variables


      INTEGER :: jk,jj,ji,jn


       trcnxtparttime = MPI_WTIME() ! cronometer-start

!! 1. fields at the next time
!! --------------------------
!! Tracer slab
!! ===========


      TRACER_LOOP: DO  jn = 1, jptra


!! 1. Lateral boundary conditions on tra (1,1,1,jn)

#ifdef key_mpp

! ... Mpp : export boundary values to neighboring processors

         CALL mpplnk_my(tra(1,1,1,jn))

#  else

! ... T-point, 3D array, full array tra(1,1,1,jn) is initialised

         CALL lbc( tra(1,1,1,jn), 1, 1, 1, 1, jpk, 1 )
#endif


!!!$omp   parallel default(none) private(mytid,jk,jj,ji)
!!!$omp&      shared(jn,jpk,jpj,jpi,trn,trb,tra,tmask,e3t,e3t_back) 




            DO ji = 1,jpi
          DO jj = 1,jpj
        DO jk = 1,jpk
            
            tra(jk,jj,ji,jn  ) = tra(jk,jj,ji,jn  )*e3t_back(jk,jj,ji)/e3t(jk,jj,ji)
            trb(jk,jj,ji,jn  ) = tra(jk,jj,ji,jn  )
            trn(jk,jj,ji,jn  ) = tra(jk,jj,ji,jn  )*tmask(jk,jj,ji)
            tra(jk,jj,ji,jn  ) = 0.e0
           ! print *,jk,jj,ji,trb(jk,jj,ji,jn  ),trn(jk,jj,ji,jn  ),e3t_back(jk,jj,ji),e3t(jk,jj,ji)
            END DO
          END DO
        END DO


     
!!!$omp end parallel


!! END of tracer slab
!! ==================

       END DO TRACER_LOOP

      

       trcnxtparttime = MPI_WTIME() - trcnxtparttime ! cronometer-stop
       trcnxttottime = trcnxttottime + trcnxtparttime


      END SUBROUTINE trcnxt
