
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
       USE myalloc_mpp
       USE FN_mem

       IMPLICIT NONE

!!----------------------------------------------------------------------
!! local declarations

!      INTEGER :: elements,nelements(6),idx_element(14,6)
! omp variables
      INTEGER :: mytid, ntids

#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      INTEGER ji,jj,jk,jn,jv! ,jnn,gji,gjj
      REAL(8) zfact,zdt

!!----------------------------------------------------------------------
!! statement functions
!! ===================


#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif

       TRA_FN = 0.0
!      SMALL =  0.00000000001

!****************     INIT PHASE  *********************

      IF (dimen_jvsnu .EQ. 0) THEN
         DO jj = 2,jpjm1
           DO  ji = 2,jpim1
              IF(tmask(ji,jj,1) .NE. 0) THEN
                 dimen_jvsnu = dimen_jvsnu + 1
                 jarr_snu(1,dimen_jvsnu) = ji
                 jarr_snu(2,dimen_jvsnu) = jj
              ENDIF
            END DO
         END DO
      ENDIF

! ***************************************************
      snutelparttime = MPI_WTIME()

      TRACER_LOOP: DO  jn = 1, jptra, ntids

!! 3. swap of arrays
!! -----------------
!!!$omp   parallel default(none) private(mytid,ji,jj,jk)
!!!$omp&      shared(jn,jpk,jpj,jpi,tra,tmask,tra_FN,SMALL)

#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF( mytid + jn <= jptra ) THEN


            DO jk = 1,jpk
               DO jj = 1,jpj
                  DO ji = 1,jpi

                     if (tmask(ji,jj,jk).ne.0.0) then
                     if( tra(ji,jj,jk,jn+mytid) .GT. 0.  ) then

                     else
                        tra_FN(ji,jj,jk,jn+mytid) =  - tra(ji,jj,jk,jn+mytid) + SMALL
                        tra(   ji,jj,jk,jn+mytid) =  SMALL
                     end if
                     endif
                  END DO
               END DO
            END DO

         END IF
!!!$omp end parallel

      END DO TRACER_LOOP


!! Frequency of correction if plus module of kt

!       ******TO BE DEVELOPED *************************
!      CALL OPA_elements(elements,nelements,idx_element)

!      TOT(:,:) =0
!      TOT_FN(:,:)=0

!!! Compute vertical integral
!      DO jv = 1,dimen_jvsnu, ntids ! cicla sui punti terra 2D
!
! $omp   parallel default(none) private(mytid,ji,jj,jk,jn,jnn)
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
!                  TOT(   jv+mytid,jn)=TOT(   jv+mytid,jn) +   tra(ji,jj,jk,idx_element(jnn,jn))*e3t(ji,jj,jk)
!                  TOT_FN(jv+mytid,jn)=TOT_FN(jv+mytid,jn) +tra_FN(ji,jj,jk,idx_element(jnn,jn))*e3t(ji,jj,jk)
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
!                      tra(ji,jj,jk,idx_element(jnn,jn)) = FN_CORR(jv+mytid,jn) * tra(ji,jj,jk,idx_element(jnn,jn))
!                   END DO
!                END DO
!             END DO

!      ENDIF
! $omp end parallel
!      END DO ! jv

      snutelparttime = MPI_WTIME()   - snutelparttime
      snuteltottime  = snuteltottime + snutelparttime

      END SUBROUTINE SNUTEL
