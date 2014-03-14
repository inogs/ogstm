C $Id: trcini.npzd.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC---------------------------------------------------------------------
CCC
CCC                        trcini.npzd.h
CCC                       ***************
CCC
CCC  purpose :
CCC  ---------
CCC     special initialisation for NPZD model
CCC
CCC  modifications :
CC   -------------
CC      original    : 99-09 (M. Levy) 
CC      additions   : 00-12 (E. Kestenare) add sediment computations
CCC
CCC---------------------------------------------------------------------
CCC  opa8, ipsl (11/96)
CCC---------------------------------------------------------------------
CCC local variables
      REAL(8) ztest
C
C 2. initialization of field for optical model
C --------------------------------------------
C
      DO jj=1,jpj
        DO ji=1,jpi
          xze(ji,jj)=5.
        END DO
      END DO

      DO jk=1,jpk
        DO jj=1,jpj
          DO ji=1,jpi
            xpar(ji,jj,jk)=0.
          END DO
        END DO
      END DO
C
C 3. initialization for passive tracer remineralisation-damping  array
C -------------------------------------------------------------------------
C 
      DO jn=1,jptra
        DO jk=1,jpk
          remdmp(jk,jn)=tminr
        END DO 
      END DO
C
      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' trcini: compute remineralisation-damping  '
          WRITE(numout,*) '         arrays for tracers'
      ENDIF
C
C 4. sediments: no martin's remineralisation profile
C -------------------------------------------------------------------------
C
      dminl = 0.
      dmin3 = 0. 
C
C    CALCUL DU MASK DE COTE
C
        cmask=0.
        do ji=2,jpi-1
          do jj=2,jpj-1
            if (tmask(ji,jj,1).eq.1) then
             ztest=tmask(ji+1,jj,1)*tmask(ji-1,jj,1)*tmask(ji,jj+1,1)
     .             *tmask(ji,jj-1,1)
             if (ztest.eq.0) cmask(ji,jj)=1.
             endif
          end do
        end do

        cmask(1,:)=cmask(jpi-1,:)
        cmask(jpi,:)=cmask(2,:)
C
C     CALCUL DE LA SURFACE COTIERE
C
         do ji=2,jpi-1
          do jj=2,jpj-1
          areacot=areacot+e1t(ji,jj)*e2t(ji,jj)*cmask(ji,jj)
          end do
         end do

