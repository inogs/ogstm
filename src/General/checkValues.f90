      SUBROUTINE CHECKVALUES
      use myalloc
      use myalloc_mpp

      IMPLICIT NONE

      INTEGER jn
      INTEGER ji, jj, jk
      REAL(8) maxV
      CHARACTER*3 varname
      CHARACTER*55  STR

! omp variables
      INTEGER :: mytid, ntids


#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif

      checkVparttime = MPI_WTIME()
      STR='TRACER EXCEPTION: Maximum allowed exceeded in [i,j,k]= '





      DO jn=1, jptra, ntids
!!!$omp parallel default(none) private(mytid, ji,jj,jk,varname,maxV)
!!!$omp&         shared(jn, jpk,jpj,jpi, tra, rank, tmask,ctrcnm,STR,ctrmax,isCheckLOG)

#ifdef __OPENMP1
       mytid = omp_get_thread_num()  ! take the thread ID
#endif

       IF (jn+mytid.le.jptra) then

         maxV     = ctrmax(jn+mytid)
         varname  = ctrcnm(jn+mytid)

         DO jk = 1,jpk
         DO jj = 1,jpj
         DO ji = 1,jpi
            if (tra(ji,jj,jk,jn+mytid).gt.maxV) THEN
                tra(ji,jj,jk,jn+mytid) = maxV*0.2 * tmask(ji,jj,jk)
                if (isCheckLOG) write(*,320) STR, ji,jj,jk, '  rank-> ', rank,' tracer ',varname

            endif
         ENDDO
         ENDDO
         ENDDO
      END IF

!!!$omp end parallel
      ENDDO

320   FORMAT (A,I3, I3, I3, A, I3, A,A)
!321   FORMAT (A,I3, I3, I3, A, I3, A,A)

      checkVparttime = MPI_WTIME() - checkVparttime
      checkVtottime = checkVtottime + checkVparttime

      END SUBROUTINE CHECKVALUES

