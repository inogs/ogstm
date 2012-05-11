
      SUBROUTINE snutel()
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE snutel
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC	compensate along the water column to preserve total mass
CCC
CCC
CC   METHOD :
CC   -------
CC    Must be called between trczdf and trcnxt
CC
C
CC
CC   OUTPUT :
CC   ------
CC
CC
CC   EXTERNAL :
CC   --------
CC
CC   REFERENCES :
CC   ----------
CC
CC   MODIFICATIONS:
CC   --------------

CC----------------------------------------------------------------------

       USE myalloc
       USE myalloc_mpp
       USE FN_mem

       IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================

      INTEGER :: kt,pack_size
      INTEGER :: elements,nelements(6),idx_element(14,6)
! omp variables
      INTEGER :: mytid, ntids, itid

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      INTEGER ji,jj,jk,jn,jv,jnn,gji,gjj
      REAL(8) zfact,zdt

CC----------------------------------------------------------------------
CC statement functions
CC ===================


#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif

      SMALL =  0.00000000001

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


      TRACER_LOOP: DO  jn = 1, jptra, ntids

C 3. swap of arrays
C -----------------
!$omp   parallel default(none) private(mytid,ji,jj,jk)
!$omp&      shared(jn,jpk,jpj,jpi,trn,trb,tra,tmask,tra_FN,SMALL)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF( mytid + jn <= jptra ) THEN


            DO jk = 1,jpk
               DO jj = 1,jpj
                  DO ji = 1,jpi


                     if( tra(ji,jj,jk,jn+mytid) .GT. 0.  ) then

                     else
                        tra_FN(ji,jj,jk,jn+mytid) =  - tra(ji,jj,jk,jn+mytid) + SMALL
                        tra(ji,jj,jk,jn+mytid)    = SMALL
                     end if

                  END DO
               END DO
            END DO

         END IF
!$omp end parallel

      END DO TRACER_LOOP


C Frequency of correction if plus module of kt


      CALL OPA_elements(elements,nelements,idx_element)

      TOT(:,:) =0
      TOT_FN(:,:)=0

C Compute vertical integral
      DO jn=1, elements
         DO jnn=1, nelements(jn)
            DO jv = 1,dimen_jvsnu
               ji = jarr_snu(1,jv)
               jj = jarr_snu(2,jv)
               DO jk = 1,jpk
                  TOT(jv,jn)=TOT(jv,jn)+tra(ji,jj,jk,idx_element(jnn,jn))*e3t(jk)
                  TOT_FN(jv,jn)=TOT_FN(jv,jn)+tra_FN(ji,jj,jk,idx_element(jnn,jn))*e3t(jk)
               END DO
            END DO
         END DO
      END DO

C Compute correction factor

      DO jn=1, elements
         DO jv = 1,dimen_jvsnu

            FN_CORR(jv,jn)=1-TOT_FN(jv,jn)/TOT(jv,jn)

            IF(FN_CORR(jv,jn) .LE. 0) THEN

               ji = jarr_snu(1,jv)
               jj = jarr_snu(2,jv)

               gji = idxt2glo(ji, jj, 1,1)
               gjj = idxt2glo(ji, jj, 1,2)

               write(*,*) 'SNUTEL --> NEGATIVE MASS PROBLEM: ji=', gji, ' jj= ', gjj, 'element= ',jn

               FN_CORR(jv,jn) = 1.

!              STOP 'NEGATIVE MASS PROBLEM'

            END IF

         END DO
      END DO

C Apply correction factor

      DO jn=1, elements
         DO jnn=1, nelements(jn)
            DO jk = 1,jpk
               DO jv = 1,dimen_jvsnu
                  ji = jarr_snu(1,jv)
                  jj = jarr_snu(2,jv)
                  tra(ji,jj,jk,idx_element(jnn,jn)) = FN_CORR(jv,jn) * tra(ji,jj,jk,idx_element(jnn,jn))
               END DO
            END DO
         END DO
      END DO


      RETURN
      END
