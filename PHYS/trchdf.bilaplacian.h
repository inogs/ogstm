CCC        TRCHDF.BILAPLACIAN
CCC      **********************
CCC
CC   define key : 'key_trc_hdfbilap'   but not key_trc_hdfgeop
CC   ==========
CC
CC
CC   METHOD :
CC   -------
CC      4th order diffusive operator along model level surfaces evalu-
CC    ated using before fields (forward time scheme). The horizontal
CC    diffusive trends of passive tracer is given by:
CC
CC       * s-coordinate ('key_s_coord' defined), the vertical scale 
CC    factors e3. are inside the derivatives:
CC      Laplacian of trb:
CC         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(trb) ]
CC                                  + dj-1[ e1v*e3v/e2v dj(trb) ]  }
CC    Multiply by the eddy diffusivity coef. and insure lateral bc:
CC       zlt   = fsahtt * zlt
CC       call to lbc or mpplnk2
CC      Bilaplacian (laplacian of zlt):
CC         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
CC                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
CC
CC       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
CC      Laplacian of trb:
CC         zlt   = 1/(e1t*e2t) {  di-1[ e2u/e1u di(trb) ]
CC                              + dj-1[ e1v/e2v dj(trb) ] }
CC    Multiply by the eddy diffusivity coef. and insure lateral bc:
CC       zlt   = fsahtt * zlt
CC       call to lbc or mpplnk2
CC      Bilaplacian (laplacian of zlt):
CC         difft = 1/(e1t*e2t) {  di-1[ e2u/e1u di(zlt) ]
CC                              + dj-1[ e1v/e2v dj(zlt) ]  }
CC
CC      Add this trend to the general trend (tra):
CC         (tra) = (tra) + ( difftr )
CC
CC      'key_trc_diatrd' defined: the trend is saved for diagnostics.
CC
CC      macro-tasked on tracer slab (jn-loop)
CC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC            /COMCOO/          : scale factors
CC            /COMASK/          : masks
CC            /COTTRP/          : passive tracer fields
CC            /comhdt/          : eddy diffusivity
CC            /COMTSK/          : multitasking
CC
CC   OUTPUT :
CC   ------
CC      common
CC            /COTTRP/ tra      : general passive tracer trend 
CC                                increased by the
CC                                horizontal diffusion trend
CC            /COTRTD/ trtrd    : horizontal tracer diffusion trend
CC                                ('key_trc_diatrd' defined)
CC
CC   EXTERNAL :      lbc, mpplnk2
CC   --------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original : 87-06 (P. Andrich - D. L Hostis)
CC      addition : 91-11 (G. Madec)
CC      addition : 93-03 (M. Guyon) symetrical conditions
CC      addition : 95-11 (G. Madec) suppress volumetric scale factors
CC      addition : 96-01 (G. Madec) statement function for e3
CC                                  suppression of common work arrays
CC      addition : 96-01 (M. Imbard) mpp exchange
CC      addition : 97-02 (M.A. Foujols) passive tracer modification
CC      addition : 97-07 (G. Madec) optimization, and fsahtt
CC      addition : 98-04 keeps trends in X and Y
CC      addition : 99-02 (M.A. Foujols) lbc in conjonction with ORCA
CC      modification : 00-05 (MA Foujols) add lbc for tracer trends
CC      modification : 00-10 (MA Foujols E Kestenare) USE passive tracer
CC                            coefficient
CC----------------------------------------------------------------------
      USE myalloc
      USE myalloc_mpp
      USE HDF_mem

        IMPLICIT NONE



CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ktask,kt
#if defined key_passivetrc 
      INTEGER ji,jj,jk,jn,jv,mytid,ntids,pack_size!,itid
      INTEGER myji,myjj
      INTEGER locsum,jklef,jjlef,jilef,jkrig,jjrig,jirig
#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif
CC----------------------------------------------------------------------
CC statement functions
CC ===================

       trcbilaphdfparttime = MPI_WTIME()
CCC OpenMP
#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif



CC Define auxiliary matrix

C       dimen_jvhdf1=0
       IF (dimen_jvhdf1 .EQ. 0) THEN
          DO jk = 1,jpk
             DO jj = 1,jpj
                DO ji = 1,jpi
                   jklef = -1
                   jjlef = -1
                   jilef = -1
                   jkrig = +1
                   jjrig = +1
                   jirig = +1
                   if(jk .EQ. 1)   jklef = 0
                   if(jj .EQ. 1)   jjlef = 0
                   if(ji .EQ. 1)   jilef = 0
                   if(jk .EQ. jpk) jkrig = 0
                   if(jj .EQ. jpj) jjrig = 0
                   if(ji .EQ. jpi) jirig = 0
                   locsum = 0
C                  DO myjk=jk+jklef, jk+jkrig
                      DO myjj=jj+jjlef, jj+jjrig
                         DO myji=ji+jilef, ji+jirig
                            locsum = locsum + tmask(myji, myjj, jk)
                         END DO
                      END DO
C                  END DO
                   if(locsum .NE. 0) then
                      dimen_jvhdf1 = dimen_jvhdf1 + 1
                      hdfmask(ji,jj,jk) = 1
                   else
                      hdfmask(ji,jj,jk) = 0
                   endif
                END DO
             END DO
          END DO

C 0. Initialization of metric arrays (for z- or s-coordinates)
C ----------------------------------
          DO jk=1,jpkm1
             DO jj = 1, jpjm1
               DO ji = 1, jpim1
#if defined key_s_coord
C   ... s-coordinates, vertical scale factor are used
                  zbtr(ji,jj,jk) = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk)
     $                        )
                  zeeu(ji,jj,jk) = e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
     $                        * umask(ji,jj,jk)
                  zeev(ji,jj,jk) = e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)
     $                        * vmask(ji,jj,jk)
#else
C   ... z-coordinates, no vertical scale factors
                  zbtr(ji,jj) = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
                  zeeu(ji,jj,jk) = e2u(ji,jj) / e1u(ji,jj) * umask(ji,jj,jk)
                  zeev(ji,jj,jk) = e1v(ji,jj) / e2v(ji,jj) * vmask(ji,jj,jk)
#endif
               END DO
             END DO
          END DO

       ENDIF

       IF (dimen_jvhdf2 .EQ. 0) THEN
          DO jk = 1,jpkm1
             DO jj = 1,jpjm1
                DO  ji = 1,jpim1
                   IF(hdfmask(ji,jj,jk) .NE. 0) THEN
                      dimen_jvhdf2 = dimen_jvhdf2 + 1
                      jarr_hdf(1,dimen_jvhdf2,1) = ji
                      jarr_hdf(2,dimen_jvhdf2,1) = jj
                      jarr_hdf(3,dimen_jvhdf2,1) = jk
                   ENDIF
                END DO
             END DO
          END DO
       ENDIF

C       dimen_jvhdf3=0

       IF (dimen_jvhdf3 .EQ. 0) THEN
          DO jk = 1,jpkm1
             DO jj = 2,jpjm1
                DO  ji = 2,jpim1
                   IF(hdfmask(ji,jj,jk) .NE. 0) THEN
                      dimen_jvhdf3 = dimen_jvhdf3 + 1
                      jarr_hdf(1,dimen_jvhdf3,2) = ji
                      jarr_hdf(2,dimen_jvhdf3,2) = jj
                      jarr_hdf(3,dimen_jvhdf3,2) = jk
                   ENDIF
                END DO
             END DO
          END DO
       ENDIF
C
C tracer slab
C =============
C
      TRACER_LOOP: DO  jn = ktask, jptra, ntids

C 1. Laplacian
C ------------
C
C ... First derivative (gradient)
!$omp  parallel default(none) private(mytid,jv,jk,jj,ji)
!$omp&                        shared(jn,dimen_jvhdf2,jarr_hdf,ztu,zeeu,trb,tmask,ztv,zeev,
!$omp&                               dimen_jvhdf3,zlt,zbtr,trcrat,ahtt)
#ifdef __OPENMP
       mytid = omp_get_thread_num()  ! take the thread ID
#endif
         IF( mytid + jn <= jptra ) THEN

          DO jv=1, dimen_jvhdf2

             ji = jarr_hdf(1,jv,1)
             jj = jarr_hdf(2,jv,1)
             jk = jarr_hdf(3,jv,1)

             ztu(ji,jj,jk,mytid+1) = zeeu(ji,jj,jk) *
     $          ( trb(ji+1,jj,jk,jn+mytid) - trb(ji,jj,jk,jn+mytid) )*
     $          tmask(ji+1,jj,jk) * tmask(ji,jj,jk)
             ztv(ji,jj,jk,mytid+1) = zeev(ji,jj,jk) *
     $          ( trb(ji,jj+1,jk,jn+mytid) - trb(ji,jj,jk,jn+mytid) )*
     $          tmask(ji,jj+1,jk) * tmask(ji,jj,jk)

          END DO
C
C ... Second derivative (divergence)
          DO jv=1, dimen_jvhdf3

             ji = jarr_hdf(1,jv,2)
             jj = jarr_hdf(2,jv,2)
             jk = jarr_hdf(3,jv,2)

             zlt(ji,jj,jk,mytid+1) = (  ztu(ji,jj,jk,mytid+1) - ztu(ji-1,jj,jk,mytid+1)
     $             + ztv(ji,jj,jk,mytid+1) - ztv(ji,jj-1,jk,mytid+1)  ) * zbtr(ji,jj)
C ... Multiply by the eddy diffusivity coefficient
             zlt(ji,jj,jk,mytid+1) = fsahtrt(ji,jj,jk) * zlt(ji,jj,jk,mytid+1)

          END DO

       ENDIF
!$omp  end parallel
C
C
C ... Lateral boundary conditions on the laplacian (zlt,zls)

#ifdef key_mpp
C
C   ... Mpp : export boundary values to neighboring processors
C
       IF( ntids - 1 + jn <= jptra ) THEN
          pack_size = ntids
       ELSE
          pack_size = ntids - (ntids - 1 + jn - jptra)
       END IF

       CALL mpplnk_my(zlt(:,:,:,:), pack_size,1,1)

#else

        DO itid = 1, ntids

           IF( itid - 1 + jn <= jptra ) THEN

              CALL lbc( zlt(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )

           END IF

        END DO
C
#endif
C
C
C 2. Bilaplacian
C --------------
C
C ... third derivative (gradient)
!$omp  parallel default(none) private(mytid,jv,jk,jj,ji)
!$omp&                        shared(jn,dimen_jvhdf2,jarr_hdf,ztu,zeeu,zlt,tmask,ztv,zeev,
!$omp&                               dimen_jvhdf3,zta,zbtr,tra)
#ifdef __OPENMP
       mytid = omp_get_thread_num()  ! take the thread ID
#endif
       IF( mytid + jn <= jptra ) THEN

          DO jv=1, dimen_jvhdf2

             ji = jarr_hdf(1,jv,1)
             jj = jarr_hdf(2,jv,1)
             jk = jarr_hdf(3,jv,1)


             ztu(ji,jj,jk,mytid+1) = zeeu(ji,jj,jk) *
     $        ( zlt(ji+1,jj,jk,mytid+1) - zlt(ji,jj,jk,mytid+1) ) *
     $        tmask(ji+1,jj,jk) * tmask(ji,jj,jk)
             ztv(ji,jj,jk,mytid+1) = zeev(ji,jj,jk) *
     $        ( zlt(ji,jj+1,jk,mytid+1) - zlt(ji,jj,jk,mytid+1) ) *
     $        tmask(ji,jj+1,jk) * tmask(ji,jj,jk)

          END DO
C
C ... fourth derivative (divergence) and add to the general tracer trend

          DO jv=1, dimen_jvhdf3

             ji = jarr_hdf(1,jv,2)
             jj = jarr_hdf(2,jv,2)
             jk = jarr_hdf(3,jv,2)

C   ... horizontal diffusive trends
             zta(mytid+1) = (  ztu(ji,jj,jk,mytid+1) - ztu(ji-1,jj,jk,mytid+1)
     $            + ztv(ji,jj,jk,mytid+1) - ztv(ji,jj-1,jk,mytid+1)  ) * zbtr(ji,jj)
C   ... add it to the general tracer trends
              tra(ji,jj,jk,jn+mytid) = tra(ji,jj,jk,jn+mytid) + zta(mytid+1)
#if defined key_trc_diatrd
C   ... save the horizontal diffusive trends in X and Y
              trtrd(ji,jj,jk,jn,4) = (  ztu(ji,jj) - ztu(ji-1,jj) )
     $                         * zbtr(ji,jj)
              trtrd(ji,jj,jk,jn,5) = (  ztv(ji,jj) - ztv(ji-1,jj) )
     $                         * zbtr(ji,jj)
#endif

         END DO

      ENDIF
!$omp  end parallel

C Lateral boundary conditions on trtrd:
#      if defined key_trc_diatrd
#         ifdef key_mpp
        CALL mpplnk( trtrd(1,1,1,jn,5), 1, 1 )
#         else      
        CALL lbc( trtrd(1,1,1,jn,5), 1, 1, 1, 1, jpk, 1 )
#         endif
#      endif
C
C
C End of slab
C ===========

        END DO TRACER_LOOP

       trcbilaphdfparttime = MPI_WTIME() - trcbilaphdfparttime
       trcbilaphdftottime = trcbilaphdftottime + trcbilaphdfparttime

#  else

C       no passive tracers

#endif

