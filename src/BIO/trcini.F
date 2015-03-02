      SUBROUTINE trcini
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcini
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     initialization for passive tracer
CCC     for restart or not


       USE myalloc
       IMPLICIT NONE


#if defined key_passivetrc

      INTEGER ji,jj,jk,jn


#    if defined key_trc_npzd
#    elif defined key_trc_bfm

      DO jk=1,jpk
        DO jj=1,jpj
          DO ji=1,jpi
            xpar(ji,jj,jk)=0.
          END DO
        END DO
      END DO

      call BFM0D_NO_BOXES(1,1,1,1,1)
      call Initialize()
#    elif defined key_trc_nnpzddom
#    elif defined key_trc_hamocc3
#    else
#    endif
C
C 1. initialization of passives tracers field
C -------------------------------------------
      DO jn=1,jptra
        DO jk=1,jpk
          DO jj=1,jpj
            DO ji=1,jpi
              trn(ji,jj,jk,jn)=0.e0
              tra(ji,jj,jk,jn)=0.e0
            END DO
          END DO
        END DO
      END DO
#    if defined key_trc_diaadd
#    endif
#    if defined key_trc_diabio
#    endif

#    if defined key_trc_diatrd
#    endif      
      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' trncini: generic initialisation done '
          WRITE(numout,*) ' '
      ENDIF
C
C
#else
C
C no passive trnacers
C
#endif

      END SUBROUTINE trcini
