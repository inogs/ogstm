      SUBROUTINE trcdtr
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcdtr
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     computes or READ initial DATA for passive tracer

       USE myalloc
       USE stafun
       IMPLICIT NONE


#if defined key_passivetrc
CC local declarations
CC ==================
      INTEGER ji,jj,jk,jn



      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' *** trcdtr initialisation for '
          WRITE(numout,*) '     passive tracers'
          WRITE(numout,*) ' '
      ENDIF

#if defined key_trc_generic

        DO jk = 1,jpk
          DO jj = 1,jpj
            DO ji = 1,jpi
              trn(ji,jj,jk,1)= 0.
              trn(ji,jj,jk,2)= 0.
              trn(ji,jj,jk,3)= 0.
              trn(ji,jj,jk,4)= 0.
              trn(ji,jj,jk,5)= 0.
            ENDDO
          ENDDO
        ENDDO

#elif defined key_trc_npzd

        DO jk = 1,jpk
          DO jj = 1,jpj
            DO ji = 1,jpi
              trn(ji,jj,jk,1)= 1.
              trn(ji,jj,jk,2)= 0.1
              trn(ji,jj,jk,3)= 1.
              trn(ji,jj,jk,4)= 2.
            ENDDO
          ENDDO
        ENDDO

#elif defined key_trc_npzdb
da inserire       
#elif key_trc_bfm

#endif



       DO jn=1,jptra 
         DO  jk=1,jpk
           DO  jj=1,jpj
             DO  ji=1,jpi
               trb(ji,jj,jk,jn)=trn(ji,jj,jk,jn)!  before field :
             END DO
           END DO
         END DO
       END DO



#else
C     no passive tracers
#endif

      END SUBROUTINE trcdtr

