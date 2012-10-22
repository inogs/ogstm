      SUBROUTINE h3cice
#if defined key_passivetrc && defined key_trc_hamocc3
c -------------------------------------------------------------------
c
c sous programme de calcul de la fraction de glace
c   si c est le premier passage :
c lecture des fractions Walsh & Zwally deja traites sur la grille
c
c  si pas de temps kt
c  interpolation des champs mensuels 
c
c -------------------------------------------------------------------
c

       USE myalloc
       IMPLICIT NONE

      INTEGER ji, jj, jk,mo
      INTEGER nvg,ios
      INTEGER ntrice1,ntrice2
      REAL(8)    xcice,xcicem(12)     
      REAL(8)    ciceyr(jpi,jpj)
      
      nvg     = 95
      WRITE(numout,*) 'Fraction de glace'
      WRITE(numout,*)
      WRITE(numout,1111) nvg
 1111 FORMAT ('Glace:nvg :',i7)
      WRITE(numout,'('' ---------------------------------- '')')
c
C Lecture  du fichier glace de mer
C--------------------------------
      DO mo=1,12
#    if defined key_vpp
        CALL read2s(nvg,cicemo(1,1,mo),jpi,jpj)  
#    else
        READ (nvg) cicemo(:,:,mo)
#    endif
      END DO  

c Mask monthly sea ice with tracer mask of layer 1
C-------------------------------------------------
      DO mo=1,12
        DO jj=1,jpj
          DO ji=1,jpi
            cicemo(ji,jj,mo) = cicemo(ji,jj,mo) * tmask(ji,jj,1)
          END do
        END do
      END do

c -------------------------------------------------------
c         Any sea ice coverage below 0.2 is considered open ocean
c         (see correspondence with Bonnie Samuels)
c--------------------------------------------------------
      DO mo=1,12
        DO jj=1,jpj
          DO ji=1,jpi
            IF (cicemo(ji,jj,mo) .LT. 0.2) THEN
                cicemo(ji,jj,mo) = 0.0
            ELSE IF (cicemo(ji,jj,mo) .GT. 1.0) THEN
                cicemo(ji,jj,mo) = 1.0
            ENDIF
          END DO
        END DO
      END DO
c
c--------------------------------------------------------
c         Compute mean annual sea ice
c--------------------------------------------------------
      ciceyr = 0.0
      DO jj=1,jpj
        DO ji=1,jpi
          DO mo=1,12
            ciceyr(ji,jj) = ciceyr(ji,jj) + cicemo(ji,jj,mo)/12.
          END DO
        END DO
      END DO
c
c--------------------------------------------------------
c         Any sea ice coverage below 0.2 is considered open ocean
c         (see correspondence with Bonnie Samuels)
c--------------------------------------------------------
      DO jj=1,jpj
        DO ji=1,jpi
          IF (ciceyr(ji,jj) .LT. 0.2) THEN
              ciceyr(ji,jj) = 0.0
          ELSE IF (ciceyr(ji,jj) .GT. 1.0) THEN
              ciceyr(ji,jj) = 1.0
          ENDIF 
        END DO
      END DO
c
c--------------------------------------------------------
c  Determine annual mean (area weighted) ocean sea ice coverage
c--------------------------------------------------------
      xcice = 0.0
      DO jj=1,jpj
        DO ji=2,jpim1
          xcice = xcice + ciceyr(ji,jj) * e1t(ji,jj)*e2t(ji,jj)
     $        *tmask(ji,jj,1)
        END DO
      END DO
c----------------------------------------------------------
c  Determine monthly mean (area weighted) ocean sea ice coverage
c----------------------------------------------------------
      xcicem = 0.0
      DO mo=1,12
        DO jj=1,jpj
          DO ji=2,jpim1
            xcicem(mo) = xcicem(mo) + cicemo(ji,jj,mo) * 
     $          e1t(ji,jj)*e2t(ji,jj)*tmask(ji,jj,1)
          END DO 
        END DO 
      END DO 
      WRITE (numout,*) 'annual mean sea ice coverage = ', xcice
      WRITE (numout,*) 'mean monthly sea ice coverage = ',xcicem

c
#endif
      RETURN
      END
