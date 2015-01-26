      SUBROUTINE h3copt
#if defined key_passivetrc && defined key_trc_hamocc3
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE h3copt
CCC                     *****************
CCC
CCC  PURPOSE :
CCC  ---------
CCC         Compute regional deviation from optimal productivity
CCC            depending on the latitude and the month
CCC
CC   METHOD :
CC   -------
CC      
CC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC      common
CC              all the common defined in ogstm
CC
CC
CC   OUTPUT :                   : no
CC   ------
CC
CC   WORKSPACE :
CC   ---------
CC
CC   EXTERNAL :
CC   --------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : E. Maier-Reimer (GBC 1993)
CC      additions : C. Le Quere (1999)
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

       USE myalloc
       IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji, jj
      INTEGER iyy, in
      REAL(8) rum, delta, codel, sidel, phi, argu, t1
C
C Get year
C --------
C
      iyy = ndastp/10000
C
#     if defined key_off_tra
      IF (nbisex.EQ.1 .AND. MOD(iyy,4).EQ.0) THEN
          rum = FLOAT(ijulian)/366.
      ELSE
          rum = FLOAT(ijulian)/365.
      ENDIF
#     else
      IF (nbisex.EQ.1 .AND. MOD(iyy,4).EQ.0) THEN
          rum = FLOAT(njulian)/366.
      ELSE
          rum = FLOAT(njulian)/365.
      ENDIF
#     endif
C
      delta = -COS(rum*rpi*2.)*rpi*23.5/180.
      codel = COS(delta)
      sidel = SIN(delta)

      DO jj = 1,jpj
        DO ji = 1,jpi
          phi = gphit(ji,jj)*rpi/180.
          argu = -SIN(phi)*sidel/(COS(phi)*codel+1.e-10)
          strn(ji,jj) = 0.
          IF (argu.LE.-1.) THEN
              str = rpi*SIN(phi)*sidel
          ENDIF
          IF (argu.GE.1.) str=0.
          IF (ABS(argu).LT.1.) THEN
              t1 = ACOS(argu)
              str = t1*sin(phi)*sidel + sin(t1)*COS(phi)*codel
          ENDIF
             strn(ji,jj) = str/(30.5*24.*3600.)
        END DO
      END DO
C
#endif
      RETURN
      END

