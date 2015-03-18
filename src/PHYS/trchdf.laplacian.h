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
CC      macro-tasked on each tracer (slab) (jn-loop)
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
      USE myalloc_mpp
        IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================

      INTEGER ji, jj, jk, jn
      REAL(8) zabe1, zabe2, zbtr, ztra
      REAL(8) ztrx(jpi,jpj), ztry(jpi,jpj), zwtr(jpi,jpj)
CC----------------------------------------------------------------------
CC statement functions
CC ===================

       trclaphdfparttime = MPI_WTIME() ! cronometer-start

C Tracer slab
C =============
      DO 1000 jn = 1,jptra,ncpu

C 1. Laplacian for passive tracers
C --------------------------------

        DO jk = 1,jpkm1

C 1.1 First derivative (gradient)

        DO jj=1,jpjm1
          DO ji=1,jpim1

            zabe1 = trcrat * ahtu(jk) * umask(ji,jj,jk)
     $            * e2u(ji,jj) / e1u(ji,jj)
            zabe2 = trcrat * ahtv(jk) * vmask(ji,jj,jk)
     $            * e1v(ji,jj) / e2v(ji,jj)
              ztrx(ji,jj) = zabe1
     $                    * (trb(ji+1,jj,jk,jn) - trb(ji,jj,jk,jn))
              ztry(ji,jj) = zabe2
     $                    * (trb(ji,jj+1,jk,jn) - trb(ji,jj,jk,jn))

            END DO
          END DO
C 2. Second derivative (divergence)
C --------------------
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) )
              ztra = ( ztrx(ji,jj) - ztrx(ji-1, jj )
     $             + ztry(ji,jj) - ztry(ji, jj-1) ) * zbtr
              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra

            END DO
          END DO
        END DO

C END of slab
C ===========

 1000 CONTINUE


       trclaphdfparttime = MPI_WTIME() - trclaphdfparttime
       trclaphdftottime = trclaphdftottime + trclaphdfparttime


