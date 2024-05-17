
      SUBROUTINE snutel()
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE snutel
!!!                     *******************
!!!
!!!  PURPOSE :
!!!  ---------
!!!	compensate along the water column to preserve total mass
!!!
!!!
!!   METHOD :
!!   -------
!!    Must be called between trczdf and trcnxt
!!----------------------------------------------------------------------

       USE myalloc
       ! epascolo USE myalloc_mpp
       USE FN_mem
       use mpi

       IMPLICIT NONE

!!----------------------------------------------------------------------
!! local declarations

!      INTEGER :: elements,nelements(6),idx_element(14,6)
! omp variables

      INTEGER :: jk,jj,ji,jn,jv! ,jnn,gji,gjj
      INTEGER :: queue
      double precision :: zfact,zdt

!!----------------------------------------------------------------------
!! statement functions
!! ===================

       queue=1
       !$acc kernels default(present) async(queue)
       TRA_FN = 0.0
       !$acc end kernels
!      SMALL =  0.00000000001

!****************     INIT PHASE  *********************

      IF (dimen_jvsnu .EQ. 0) THEN
         DO jj = 2,jpjm1
           DO  ji = 2,jpim1
              IF(tmask(1,jj,ji) .NE. 0) THEN
                 dimen_jvsnu = dimen_jvsnu + 1
                 jarr_snu(2,dimen_jvsnu) = ji
                 jarr_snu(1,dimen_jvsnu) = jj
              ENDIF
            END DO
         END DO
      ENDIF

! ***************************************************
      snutelparttime = MPI_WTIME()

      TRACER_LOOP: DO  jn = 1, jptra

!! 3. swap of arrays
!! -----------------
!!!$omp   parallel default(none) private(mytid,jk,jj,ji)
!!!$omp&      shared(jn,jpk,jpj,jpi,tra,tmask,tra_FN,SMALL)


         !$acc parallel loop collapse(3) default(present) async(queue)
         DO ji = 1,jpi
            DO jj = 1,jpj
               DO jk = 1,jpk
                  if (tmask(jk,jj,ji).ne.0.0) then
                     if( tra(jk,jj,ji,jn) .GT. 0.  ) then

                     else
                        tra_FN(jk,jj,ji,jn) =  - tra(jk,jj,ji,jn) + SMALL
                        tra(   jk,jj,ji,jn) =  SMALL
                     end if
                  endif
               END DO
            END DO
         END DO

      
!!!$omp end parallel

      END DO TRACER_LOOP

      !$acc wait(queue)


!! Frequency of correction if plus module of kt

!       ******TO BE DEVELOPED *************************
!      CALL OPA_elements(elements,nelements,idx_element)

!      TOT(:,:) =0
!      TOT_FN(:,:)=0

!!! Compute vertical integral
!      DO jv = 1,dimen_jvsnu, ntids ! cicla sui punti terra 2D
!
! $omp   parallel default(none) private(mytid,jk,jj,ji,jn,jnn)
! $omp&      shared(dimen_jvsnu,jarr_snu,jv,jpk,elements,nelements,TOT,TOT_FN,tra,tra_FN,idx_element,e3t,FN_CORR)
!
!#ifdef __OPENMP1
!        mytid = omp_get_thread_num()  ! take the thread ID
!#endif
!      IF( mytid + jv <= dimen_jvsnu ) THEN
!
!          ji = jarr_snu(1,jv+mytid)
!          jj = jarr_snu(2,jv+mytid)
!          TOT(   jv+mytid,:) = 0.
!          TOT_FN(jv+mytid,:) = 0.
!
!
!          DO jn =1, elements
!          DO jnn=1, nelements(jn)
!
!               DO jk = 1,jpk
!                  TOT(   jv+mytid,jn)=TOT(   jv+mytid,jn) +   tra(jk,jj,ji,idx_element(jnn,jn))*e3t(jk,jj,ji)
!                  TOT_FN(jv+mytid,jn)=TOT_FN(jv+mytid,jn) +tra_FN(jk,jj,ji,idx_element(jnn,jn))*e3t(jk,jj,ji)
!               END DO
!           END DO
!           END DO


!! Compute correction factor

!          DO jn=1, elements
!
!
!                FN_CORR(jv+mytid,jn) = 1-TOT_FN(jv+mytid,jn)/TOT(jv+mytid,jn)
!
!                IF(FN_CORR(jv+mytid,jn) .LE. 0) THEN
!
!    !               gji = idxt2glo(ji, jj, 1,1)
!    !               gjj = idxt2glo(ji, jj, 1,2)
!    !       write(*,250) 'SNUTEL --> NEGATIVE MASS PROBLEM: ji=', gji, ' jj= ', gjj, ' element= ',jn
!    !250   FORMAT(A,I5,A,I5,A,I5)
!                   FN_CORR(jv+mytid,jn) = 1.
!    !              STOP 'NEGATIVE MASS PROBLEM'
!
!                END IF
!
!          END DO


!! Apply correction factor
!
!          DO jn=1, elements
!             DO jnn=1, nelements(jn)
!                DO jk = 1,jpk
!                      tra(jk,jj,ji,idx_element(jnn,jn)) = FN_CORR(jv+mytid,jn) * tra(jk,jj,ji,idx_element(jnn,jn))
!                   END DO
!                END DO
!             END DO

!      ENDIF
! $omp end parallel
!      END DO ! jv

      snutelparttime = MPI_WTIME()   - snutelparttime
      snuteltottime  = snuteltottime + snutelparttime

      END SUBROUTINE SNUTEL
