CCC                       trchdf.laplacian.h
CCC                     *********************
CCC
CC
CC   defined key : no (default option)
CC   ==========
CC
CC
CC   Method :
CC   -------
CC    Second order diffusive operator evaluated using before fields
CC      (forward time scheme). The horizontal diffusive trends of 
CC    passive tracer is given by:
CC
CC       * s-coordinate ('key_s_coord' defined), the vertical scale 
CC      factors e3. are inside the derivatives:
CC          difft = 1/(e1t*e2t*e3t) {  di-1[ ahtt e2u*e3u/e1u di(trb) ]
CC                                   + dj-1[ ahtt e1v*e3v/e2v dj(trb) ] }
CC
CC       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
CC          difft = 1/(e1t*e2t) {  di-1[ ahtt e2u/e1u di(trb) ]
CC                               + dj-1[ ahtt e1v/e2v dj(trb) ] }
CC
CC      Add this trend to the general tracer trend (tra):
CC          (tra) = (tra) + ( difftr )
CC
CC      'key_trc_diatrd' defined: the trend is saved for diagnostics.
CC
CC      macro-tasked on each tracer (slab) (jn-loop)
CC
CC   Input :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC            /comcoo/          : scale factors
CC            /comask/          : masks
CC            /cottrp/          : passive tracer fields 
CC            /comhdf/          : horizontal diffusion
CC            /comtsk/          : multitasking
CC
CC   Output :
CC   ------
CC      common
CC            /cottrp/ tra      : general tracer trend increased by the
CC                                before horizontal diffusion trend
CC            /cotrtd/ trtrd    : now horizontal passive tracer
CC                                diffusion trend
CC                                (IF 'key_trc_diatrd' key is activated)
CC
CC   External :                   no
CC   --------
CC
CC   References :                 no
CC   ----------
CC
CC   Modifications:
CC   --------------
CC       original : 87-06 (P. Andrich - D. L Hostis)
CC       additios : 91-11 (G. Madec)
CC       addition : 95-11 (G. Madec) suppress volumetric scale factors
CC       addition : 96-01 (G. Madec) statement function for e3
CC                                   suppression of common work arrays
CC       addition : 95-02 passive tracers (M. Levy)
CC       addition : 98-04 keeps trends in X and Y 
CC      modification : 00-10 (MA Foujols E Kestenare) USE passive tracer
CC                            coefficient
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================


      USE myalloc
#include <mpif.h>
        IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ktask, kt
#if defined key_passivetrc
      INTEGER ji, jj, jk, jn
      REAL(8) zabe1, zabe2, zbtr, ztra
      REAL(8) ztrx(jpi,jpj), ztry(jpi,jpj), zwtr(jpi,jpj)
CC----------------------------------------------------------------------
CC statement functions
CC ===================

       trclaphdfparttime = MPI_WTIME() ! F79 cronometer-start

C Tracer slab
C =============
C
      DO 1000 jn = ktask,jptra,ncpu
C
C 1. Laplacian for passive tracers
C --------------------------------
C
        DO jk = 1,jpkm1
C
C 1.1 First derivative (gradient)
C
        DO jj=1,jpjm1
          DO ji=1,jpim1
#if defined key_s_coord
            zabe1 = fsahtru(ji,jj,jk) * umask(ji,jj,jk)
     $            * e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
            zabe2 = fsahtrv(ji,jj,jk) * vmask(ji,jj,jk)
     $            * e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)
#  else
            zabe1 = fsahtru(ji,jj,jk) * umask(ji,jj,jk)
     $            * e2u(ji,jj) / e1u(ji,jj)
            zabe2 = fsahtrv(ji,jj,jk) * vmask(ji,jj,jk)
     $            * e1v(ji,jj) / e2v(ji,jj)
#endif
              ztrx(ji,jj) = zabe1
     $                    * (trb(ji+1,jj,jk,jn) - trb(ji,jj,jk,jn))
              ztry(ji,jj) = zabe2
     $                    * (trb(ji,jj+1,jk,jn) - trb(ji,jj,jk,jn))

            END DO
          END DO
C
C 2. Second derivative (divergence)
C --------------------
C
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
#if defined key_s_coord
            zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) )
#  else
            zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) )
#endif
              ztra = ( ztrx(ji,jj) - ztrx(ji-1, jj )
     $             + ztry(ji,jj) - ztry(ji, jj-1) ) * zbtr
              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
C
#if defined key_trc_diatrd
              trtrd(ji,jj,jk,jn,4)=(ztrx(ji,jj)-ztrx(ji-1,jj))*zbtr
              trtrd(ji,jj,jk,jn,5)=(ztry(ji,jj)-ztry(ji,jj-1))*zbtr
#endif
C
            END DO
          END DO
C
        END DO
C
C END of slab
C ===========
C
 1000 CONTINUE

CCC 10 11 2004  F79 cronometer-stop

       trclaphdfparttime = MPI_WTIME() - trclaphdfparttime
       trclaphdftottime = trclaphdftottime + trclaphdfparttime

       write(*,*) "F79T:trclaphdfparttime", trclaphdfparttime
       write(*,*) "F79T:trclaphdftottime", trclaphdftottime

CCC

#else
C
C     no passive tracers
C
#endif
