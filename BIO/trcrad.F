      SUBROUTINE trcrad
#if defined key_passivetrc && key_trc_hamocc3
#if defined key_trahdfiso 
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcrad
CCC                     ******************
CCC
CC
CC     PURPOSE.
CC     --------
CC          *TRCRAD* "CRAPPY" ROUTINE TO CORRECT ARTIFICIAL NEGATIVE
CC                   CONCENTRATIONS DUE TO ISOPYCNAL SCHEME
CC
CC     METHOD.
CC     -------
CC          SET NEGATIVE CONCENTRATIONS TO ZERO
CC          COMPUTE THE CORRESPONDING MASS ADDED TO THE TRACERS 
CC          AND REMOVE IT WHEN POSSIBLE
CC
CC     EXTERNALS.
CC     ----------
CC          NONE.
CC
CC     REFERENCE.
CC     ----------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original      : 2001    O. Aumont (& EK)
CC ---------------------------------------------------------------------
       USE myalloc
       USE stafun
        IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================
C
      INTEGER ji, jj, jk, jn
      REAL(8) trcor,trmas,trmas2
CC----------------------------------------------------------------------
CC statement functions
CC ===================



C
        DO jn = 1, jptra
              trcor=0.
              trmas=0.
         DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              trcor = trcor+
     .  min(0.,trn(ji,jj,jk,jn)) * tmask(ji,jj,jk)
     .  *e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk)
#ifdef key_off_degrad
     .  *facvol(ji,jj,jk)
#endif
             trn(ji,jj,jk,jn)=max(0.,trn(ji,jj,jk,jn))
             trmas=trmas+trn(ji,jj,jk,jn)*tmask(ji,jj,jk)
     .  *e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk)
#ifdef key_off_degrad
     .  *facvol(ji,jj,jk)
#endif
             end do
          end do
          end do
CC 12/1/2004 Paolo - casi 0/0 ================================

          do jk=1,jpkm1
            do jj=1,jpj
               do ji=1,jpi

               trn(ji,jj,jk,jn)=max(0.,trn(ji,jj,jk,jn))

              end do
              end do
              end do



         do jk=1,jpkm1
            do jj=1,jpj
               do ji=1,jpi
 
       if (trmas.ne.0) then
 
            trn(ji,jj,jk,jn)=trn(ji,jj,jk,jn)*
     . (1.+trcor/trmas)*tmask(ji,jj,jk)
     
       endif
     
              end do
              end do
              end do
              end do
CC===========================================================
#endif
#endif
      RETURN
      END
