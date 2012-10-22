
      SUBROUTINE lbc ( ptab, ktype, ksgn,kdoloop, kjstart, kjpend, kstep )
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE lbc
CCC                     ***************
CCC
CCC  Purpose :
CCC  --------
CCC      Insure the lateral boundary conditions in non mpp configuration
CCC
CC   Method :
CC   -------
CC
CC   Input :
CC   ------
CC      argument
CC              ptab            : variable array
CC              ktype           : define the nature of the grid-point
CC                  at which ptab is defined
CC                  = 1 ,  T- and W-points
CC                  = 2 ,  U-point
CC                  = 3 ,  V-point
CC                  = 4 ,  F-point
CC        ksgn        : control of the sign change
CC                  = 0 , the sign is modified following
CC                  the type of b.c. used
CC                  = 1 , the sign of the field is un-
CC                  changed at the boundaries
CC              kdoloop         : = 1  macrostaked on jk-slabs for
CC                       2D or 3D arrays
CC                                = 2  macrostaked on jj-slabs for
CC                       a 3D array (jpi,jpj,jpk)
CC                                = 3  macrostaked on jj-slabs for
CC                       a 2D array (jpi,jpj)
CC        kjstart        : starting index of multitasked do-loop
CC        kjpend         : ending index of multitasked do-loop
CC        kstep          : step of multitasked do-loop
CC
CC   Output :
CC   -------
CC      argument
CC              ptab            : variable array with disired lateral
CC                  boundary conditions

      USE myalloc
      IMPLICIT NONE

CC local declarations
CC ==================

      INTEGER ktype, ksgn, kdoloop, kjstart, kjpend, kstep
      REAL(8) ptab(jpi,jpj,jpk)
#if ! defined key_mpp





      REAL(8) zsgn


      INTEGER ijt, iju
      INTEGER ji, jj, jk
C 0. Sign setting
C ---------------

      IF (ksgn.EQ.0) THEN
          zsgn = -1.
        ELSE
          zsgn =  1.
      ENDIF
C
      IF (kdoloop.EQ.1) THEN
C
C
C     Horizontal slab
C     ===============
C
          DO jk = kjstart, kjpend, kstep
C
C
C 1. East-West boundary conditions
C --------------------------------
C
            IF ( nperio.EQ.1 .OR. nperio.EQ.4 ) THEN
C ... cyclic
                DO jj = 1, jpj
                  ptab( 1 ,jj,jk) = ptab(jpim1,jj,jk)
                  ptab(jpi,jj,jk) = ptab(  2  ,jj,jk)
                END DO
              ELSE
C ... closed
                IF ( ktype.EQ.4 ) THEN
                    DO jj = 1, jpj
                      ptab(jpi,jj,jk) = 0.e0
                    END DO
                  ELSE
                    DO jj = 1, jpj
                      ptab( 1 ,jj,jk) = 0.e0
                      ptab(jpi,jj,jk) = 0.e0
                    END DO
                ENDIF
            ENDIF
C
C 2. North-South boundary conditions
C ----------------------------------
C
            IF ( nperio.EQ.2 ) THEN
C ... south symmetric
                IF ( ktype.EQ.1 .OR. ktype.EQ.2 ) THEN
                    DO ji = 1, jpi
                      ptab(ji, 1 ,jk) = ptab(ji,3,jk)
                      ptab(ji,jpj,jk) = 0.e0
                    END DO
                  ELSEIF ( ktype.EQ.3 .OR. ktype.EQ.4 ) THEN
                    DO ji = 1, jpi
                      ptab(ji, 1 ,jk) = zsgn * ptab(ji,2,jk)
                      ptab(ji,jpj,jk) = 0.e0
                    END DO
                ENDIF
              ELSEIF ( nperio.EQ.3 .OR. nperio.EQ.4 ) THEN
C ... north fold
                ptab(1,jpj,jk)=0.e0
                ptab(jpi,jpj,jk)=0.e0
                IF ( ktype.EQ.1 ) THEN
                    DO ji = 2, jpi
                      ijt=jpi-ji+2
                      ptab(ji, 1 ,jk) = 0.e0
                      ptab(ji,jpj,jk) = zsgn * ptab(ijt,jpj-2,jk)
                    END DO
                    DO ji = jpi/2+1, jpi
                      ijt=jpi-ji+2
                      ptab(ji,jpjm1,jk) = zsgn * ptab(ijt,jpjm1,jk)
                    END DO
                  ELSEIF ( ktype.EQ.2 ) THEN
                    DO ji = 1, jpi-1
                      iju=jpi-ji+1
                      ptab(ji, 1 ,jk) = 0.e0
                      ptab(ji,jpj,jk) = zsgn * ptab(iju,jpj-2,jk)
                    END DO
                    DO ji = jpi/2, jpi-1
                      iju=jpi-ji+1
                      ptab(ji,jpjm1,jk) = zsgn * ptab(iju,jpjm1,jk)
                    END DO
                  ELSEIF ( ktype.EQ.3 ) THEN
                    DO ji = 2, jpi
                      ijt=jpi-ji+2
                      ptab(ji, 1 ,jk) = 0.e0
                      ptab(ji,jpj-1,jk) = zsgn * ptab(ijt,jpj-2,jk)
                      ptab(ji,jpj  ,jk) = zsgn * ptab(ijt,jpj-3,jk)
                    END DO
                  ELSEIF ( ktype.EQ.4 ) THEN
                    DO ji = 1, jpi-1
                      iju=jpi-ji+1
C                      ptab(ji, 1 ,jk) = 0.e0
                      ptab(ji,jpj-1,jk) = ptab(iju,jpj-2,jk)
                      ptab(ji,jpj  ,jk) = ptab(iju,jpj-3,jk)
                    END DO
                ENDIF
              ELSE
C ... closed
                IF ( ktype.EQ.4 ) THEN
                    DO ji = 1, jpi
                      ptab(ji,jpj,jk) = 0.e0
                    END DO
                  ELSE
                    DO ji = 1, jpi
                      ptab(ji, 1 ,jk) = 0.e0
                      ptab(ji,jpj,jk) = 0.e0
                    END DO
                ENDIF
            ENDIF
C
C
C     End of slab
C     ===========
C
           END DO
C
C
         ELSEIF ( kdoloop.EQ.2 ) THEN
C
C
C     Vertical slab
C     =============
C
          DO jj = kjstart, kjpend, kstep
C
C 1. East-West boundary conditions
C --------------------------------
C
            IF ( nperio.EQ.1 .OR. nperio.EQ.4 ) THEN
C ... cyclic
                DO jk = 1, jpk
                  ptab( 1 ,jj,jk) = ptab(jpim1,jj,jk)
                  ptab(jpi,jj,jk) = ptab(  2  ,jj,jk)
                END DO
              ELSE
C ... closed
                IF ( ktype.EQ.4 ) THEN
                    DO jk = 1, jpk
                      ptab(jpi,jj,jk) = 0.e0
                    END DO
                  ELSE
                    DO jk = 1, jpk
                      ptab( 1 ,jj,jk) = 0.e0
                      ptab(jpi,jj,jk) = 0.e0
                    END DO
                ENDIF
            ENDIF
C
C
C 2. North-South boundary conditions ( only for row jj = 3 or = 2 )
C ----------------------------------
C
            IF ( nperio.EQ.2 ) THEN
C ... south symmetric
                IF ( ktype.EQ.1 .OR. ktype.EQ.2 .AND. jj.EQ.3 ) THEN
                    DO jk = 1, jpk
                      DO ji = 1, jpi
                        ptab(ji, 1 ,jk) = ptab(ji,3,jk)
                        ptab(ji,jpj,jk) = 0.e0
                      END DO
                    END DO
                  ELSEIF ( ktype.EQ.3 .OR. ktype.EQ.4.AND. jj.EQ.2 ) THEN
                    DO jk = 1, jpk
                      DO ji = 1, jpi
                        ptab(ji, 1 ,jk) = zsgn * ptab(ji,2,jk)
                        ptab(ji,jpj,jk) = 0.e0
                      END DO
                    END DO
                ENDIF
              ELSEIF ( nperio.EQ.3 .OR. nperio.EQ.4 ) THEN
C ... north fold
                IF (jj .EQ. jpj-3) THEN
                    DO jk = 1, jpk
                      ptab(1,jpj,jk)=0.e0
                      ptab(jpi,jpj,jk)=0.e0
                    END DO
                ENDIF
                IF ( ktype.EQ.1 .AND. jj.EQ.jpj-2 ) THEN
                     DO jk = 1, jpk
                       DO ji = 2, jpi
                         ijt=jpi-ji+2
                         ptab(ji, 1 ,jk) = 0.e0
                         ptab(ji,jpj,jk) = zsgn * ptab(ijt,jpj-2,jk)
                       END DO
                     END DO
                   ELSEIF ( ktype.EQ.1 .AND. jj.EQ.jpjm1 ) THEN
                     DO jk = 1, jpk
                       DO ji = jpi/2+1, jpi
                         ijt=jpi-ji+2
                         ptab(ji,jpjm1,jk) = zsgn * ptab(ijt,jpjm1,jk)
                       END DO
                     END DO
                   ELSEIF ( ktype.EQ.2 .AND. jj.EQ.jpj-2 ) THEN
                     DO jk = 1, jpk
                       DO ji = 1, jpi-1
                         iju=jpi-ji+1
                         ptab(ji, 1 ,jk) = 0.e0
                         ptab(ji,jpj,jk) = zsgn * ptab(iju,jpj-2,jk)
                       END DO
                     END DO
                   ELSEIF ( ktype.EQ.2 .AND. jj.EQ.jpjm1 ) THEN
                     DO jk = 1, jpk
                       DO ji = jpi/2, jpi-1
                         iju=jpi-ji+1
                         ptab(ji,jpjm1,jk) = zsgn * ptab(iju,jpjm1,jk)
                       END DO
                    END DO
                  ELSEIF ( ktype.EQ.3 .AND. jj.EQ.jpj-3 ) THEN
                    DO jk = 1, jpk
                      DO ji = 2, jpi
                        ijt=jpi-ji+2
                        ptab(ji, 1 ,jk) = 0.e0
                        ptab(ji,jpj,jk) = zsgn * ptab(ijt,jpj-3,jk)
                      END DO
                    END DO
                  ELSEIF ( ktype.EQ.3 .AND. jj.EQ.jpj-2 ) THEN
                    DO jk = 1, jpk                    
                      DO ji = 2, jpi
                        ijt=jpi-ji+2
                        ptab(ji,jpj-1,jk) = zsgn * ptab(ijt,jpj-2,jk)
                      END DO
                    END DO
                  ELSEIF ( ktype.EQ.4 .AND. jj.EQ.jpj-3 ) THEN
                    DO jk = 1, jpk                       
                      DO ji = 1, jpi-1
                        iju=jpi-ji+1
                        ptab(ji,jpj,jk) =  ptab(iju,jpj-3,jk)
                      END DO
                    END DO
                  ELSEIF ( ktype.EQ.4 .AND. jj.EQ.jpj-2 ) THEN
                    DO jk = 1, jpk                     
                      DO ji = 1, jpi-1
                        iju=jpi-ji+1
                        ptab(ji,jpj-1,jk) =  ptab(iju,jpj-2,jk)
                      END DO
                    END DO
                ENDIF
              ELSEIF( jj.EQ.2 ) THEN
C ... closed
                DO jk = 1, jpk
                  DO ji = 1, jpi
                    ptab(ji, 1 ,jk) = 0.e0
                    ptab(ji,jpj,jk) = 0.e0
                  END DO
                END DO
            ENDIF
C
C
C     End of slab
C     ===========
C
           END DO
C
         ELSEIF ( kdoloop.EQ.3 ) THEN
C
C
C     Vertical slab for a 2D array (jpi,jpj)
C     ============================
C
          DO jj = kjstart, kjpend, kstep
C
C 1. East-West boundary conditions
C --------------------------------
C
            IF ( nperio.EQ.1 .OR. nperio.EQ.4 ) THEN
C ... cyclic
                ptab( 1 ,jj,1) = ptab(jpim1,jj,1)
                ptab(jpi,jj,1) = ptab(  2  ,jj,1)
              ELSE
C ... closed
                IF ( ktype.EQ.4 ) THEN
                    ptab(jpi,jj,1) = 0.e0
                  ELSE
                    ptab( 1 ,jj,1) = 0.e0
                    ptab(jpi,jj,1) = 0.e0
                ENDIF
            ENDIF
C
C
C 2. North-South boundary conditions ( only for row jj = 3 or = 2 )
C ----------------------------------
C
            IF ( nperio.EQ.2 ) THEN
C ... south symmetric
                IF ( ktype.EQ.1 .OR. ktype.EQ.2 .AND. jj.EQ.3 ) THEN
                    DO ji = 1, jpi
                      ptab(ji, 1 ,1) = ptab(ji,3,1)
                      ptab(ji,jpj,1) = 0.e0
                    END DO
                  ELSEIF ( ktype.EQ.3 .OR. ktype.EQ.4.AND. jj.EQ.2 ) THEN
                    DO ji = 1, jpi
                      ptab(ji, 1 ,1) = zsgn * ptab(ji,2,1)
                      ptab(ji,jpj,1) = 0.e0
                    END DO
                ENDIF
              ELSEIF ( nperio.EQ.3 .OR. nperio.EQ.4 ) THEN
C ... north fold
                IF ( jj.EQ. jpj-3) THEN
                    ptab(1,jpj,1)=0.e0
                    ptab(jpi,jpj,1)=0.e0
                ENDIF 
                IF ( ktype.EQ.1 .AND. jj.EQ.jpj-2 ) THEN
                    DO ji = 2, jpi
                      ijt=jpi-ji+2
                      ptab(ji, 1 ,1) = 0.e0
                      ptab(ji,jpj,1) = zsgn * ptab(ijt,jpj-2,1)
                    END DO
                  ELSEIF ( ktype.EQ.1 .AND. jj.EQ.jpjm1 ) THEN
                    DO ji = jpi/2+1, jpi
                      ijt=jpi-ji+2
                      ptab(ji,jpjm1,1) = zsgn * ptab(ijt,jpjm1,1)
                    END DO
                  ELSEIF ( ktype.EQ.2 .AND. jj.EQ.jpj-2 ) THEN
                    DO ji = 1, jpi-1
                      iju=jpi-ji+1
                      ptab(ji, 1 ,1) = 0.e0
                      ptab(ji,jpj,1) = zsgn * ptab(iju,jpj-2,1)
                    END DO
                  ELSEIF ( ktype.EQ.2 .AND. jj.EQ.jpjm1 ) THEN
                    DO ji = jpi/2, jpi-1
                      iju=jpi-ji+1
                      ptab(ji,jpjm1,1) = zsgn * ptab(iju,jpjm1,1)
                    END DO
                  ELSEIF ( ktype.EQ.3 .AND. jj.EQ.jpj-3 ) THEN
                    DO ji = 2, jpi
                      ijt=jpi-ji+2
                      ptab(ji, 1 ,1) = 0.e0
                      ptab(ji,jpj,1) = zsgn * ptab(ijt,jpj-3,1)
                    END DO
                  ELSEIF ( ktype.EQ.3 .AND. jj.EQ.jpj-2 ) THEN
                    DO ji = 2, jpi
                      ijt=jpi-ji+2
                      ptab(ji,jpj-1,1) = zsgn * ptab(ijt,jpj-2,1)
                    END DO
                  ELSEIF ( ktype.EQ.4 .AND. jj.EQ.jpj-3 ) THEN
                    DO ji = 1, jpi-1
                      iju=jpi-ji+1
                      ptab(ji,jpj,1) =  ptab(iju,jpj-3,1)
                    END DO
                  ELSEIF ( ktype.EQ.4 .AND. jj.EQ.jpj-2 ) THEN
                    DO ji = 1, jpi-1
                      iju=jpi-ji+1
                      ptab(ji,jpj-1,1) =  ptab(iju,jpj-2,1)
                    END DO
                ENDIF
              ELSEIF( jj.EQ.2 ) THEN
C ... closed
                DO ji = 1, jpi
                  ptab(ji, 1 ,1) = 0.e0
                  ptab(ji,jpj,1) = 0.e0
                END DO
            ENDIF
C
C
C     End of slab
C     ===========
C
           END DO
C
         ELSE
C
C
C 3. E-R-R-O-R
C ------------
C
           IF(lwp) WRITE(numout,*) ' '
           IF(lwp) WRITE(numout,*) ' lbc: e r r o r in kdoloop argument'
           IF(lwp) WRITE(numout,*) '      we stop'
           STOP
       ENDIF
C
#  else
C      mpp computation the lateral boundary conditions is insured by
C      routines mpplnk and mpplnk2
C
#endif
C
      RETURN
      END
