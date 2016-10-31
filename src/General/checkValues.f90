      SUBROUTINE CHECKVALUES
      use myalloc
      use mpi
      

      IMPLICIT NONE

      INTEGER :: jn
      INTEGER :: ji, jj, jk
      double precision :: maxV
      CHARACTER*3 :: varname
      CHARACTER*55 ::  STR

! omp variables
   

      checkVparttime = MPI_WTIME()
      STR='TRACER EXCEPTION: Maximum allowed exceeded in [i,j,k]= '





      DO jn=1, jptra
!!!$omp parallel default(none) private(mytid, ji,jj,jk,varname,maxV)
!!!$omp&         shared(jn, jpk,jpj,jpi, tra, myrank, tmask,ctrcnm,STR,ctrmax,isCheckLOG)


    
         maxV     = ctrmax(jn)
         varname  = ctrcnm(jn)

         DO ji = 1,jpi
         DO jj = 1,jpj
         DO jk = 1,jpk
            if (tra(jk,jj,ji,jn).gt.maxV) THEN
                tra(jk,jj,ji,jn) = maxV*0.2 * tmask(jk,jj,ji)
                if (isCheckLOG) write(*,320) STR, jk,jj,ji, '  myrank-> ', myrank,' tracer ',varname

            endif
         ENDDO
         ENDDO
         ENDDO
     

!!!$omp end parallel
      ENDDO

320   FORMAT (A,I3, I3, I3, A, I3, A,A)
!321   FORMAT (A,I3, I3, I3, A, I3, A,A)

      checkVparttime = MPI_WTIME() - checkVparttime
      checkVtottime = checkVtottime + checkVparttime

      END SUBROUTINE CHECKVALUES

