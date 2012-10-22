
       SUBROUTINE h3cvfl(kt)
C**********************************************************************
C
C  computes the effect of E-P (evaporation minus precipitation) which
C  changes surface [DIC] and [DIC-14]
C
C  Input:
C  OUTPUT:
C
C  numbers in brackets refer to equation numbers in simulation design
C  document.
C
C**********************************************************************

       USE myalloc

        IMPLICIT NONE

      INTEGER kt
#if defined key_passivetrc && defined key_trc_hamocc3
      INTEGER nspyr,ji,jj,jk,jn



c ------------------------------------------------------------------
C    For Tracer 1: DIC
c ------------------------------------------------------------------
       DO jn=1,jptra
        DO jj=1,jpj
         DO ji=1,jpi
          tra(ji,jj,1,jn)= tra(ji,jj,1,jn)+tmask(ji,jj,1)
     $          *(emp(ji,jj)*1e-3-wn(ji,jj,1))/e3t(1)
     &          *trn(ji,jj,1,jn)
         END DO
        END DO
       END DO

#    if defined key_trc_diaadd
       DO jj=1,jpj
         DO ji=1,jpi
           trc2d(ji,jj,10)=(emp(ji,jj)*1e-3-wn(ji,jj,1))/e3t(1)
     &          *trn(ji,jj,1,jn)
         END DO
       END DO
#    endif
#endif
        RETURN
        END
