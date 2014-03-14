      SUBROUTINE h3cadj
#if defined key_passivetrc 
#if defined key_trc_hamocc3
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE h3cadj
CCC                     ******************
CCC
CC
CC     PURPOSE.
CC     --------
CC          *H3CADJ* CORRECT THE TRACERS CONCENTRATION IN THE WATER 
CC                   BY THE MASS LOST IN THE SEDIMENT
CC
CC     METHOD.
CC     -------
CC
CC     EXTERNALS.
CC     ----------
CC          NONE.
CC
CC     REFERENCE.
CC     ----------
CC
CC     MODIFICATIONS:
CC     --------------
CC      original      : 2001    O. Aumont (& EK)
CC
       USE myalloc
       USE stafun

        IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, jk
      REAL(8) ztal, zsil, zpo4, zvol, alka0, silic1, po4

C
C     Total mass of tracers that are sequestered in the sediment
C
      ztal=0.
      zsil=0.
      zpo4=0.
      zvol=0.
      DO jk= 1,jpk-1
       DO jj=2,jpj-1
         DO ji=2,jpi-1
      ztal = ztal + tmask(ji,jj,jk)*
     .       trn(ji,jj,jk,jptal)*e1t(ji,jj)*e2t(ji,jj)*e3t(jk)
      zsil = zsil + tmask(ji,jj,jk)*
     .       trn(ji,jj,jk,jpsil)*e1t(ji,jj)*e2t(ji,jj)*e3t(jk)
      zpo4 = zpo4 + tmask(ji,jj,jk)*
     .       trn(ji,jj,jk,jppo4)*e1t(ji,jj)*e2t(ji,jj)*e3t(jk)
      zvol = zvol + tmask(ji,jj,jk)*
     .              e1t(ji,jj)*e2t(ji,jj)*e3t(jk)
         ENDDO
        ENDDO
       ENDDO
C
C     Correction to the concentration in the water column
C     to conserve the total mass
C
       alka0  = 2.45e-3
       silic1 = 85.e-6
       po4    = 2.1e-6/po4r
C       ztal = alka0/(ztal/zvol)
C       zsil = silic1/(zsil/zvol)
C       zpo4 = po4/(zpo4/zvol)
C
C 12/1/2004 Paolo - casi 0/0 ================================
      if (ztal.eq.0) then
       ztal=0.
	 else
       ztal = alka0/(ztal/zvol)
      endif
      if (zvol.eq.0) then
       zsil=0. 
	else
       zsil = silic1/(zsil/zvol)
      endif
      if (zvol.eq.0) then
       zpo4=0. 
	else
       zpo4 = po4/(zpo4/zvol)
      endif
C=============================================================

       trn(:,:,:,jptal) = trn(:,:,:,jptal)*ztal
       trb(:,:,:,jptal) = trb(:,:,:,jptal)*ztal
       trn(:,:,:,jpsil) = trn(:,:,:,jpsil)*zsil
       trb(:,:,:,jpsil) = trb(:,:,:,jpsil)*zsil
       trn(:,:,:,jppo4) = trn(:,:,:,jppo4)*zpo4
       trb(:,:,:,jppo4) = trb(:,:,:,jppo4)*zpo4
C
#endif
#endif
      RETURN
      END SUBROUTINE h3cadj
