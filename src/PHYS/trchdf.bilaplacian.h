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
CC    Multiply by the eddy diffusivity coef. and insure lateral bc:
CC      Bilaplacian (laplacian of zlt):
CC         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
CC                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
CC
CC       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
CC      Laplacian of trb:
CC         zlt   = 1/(e1t*e2t) {  di-1[ e2u/e1u di(trb) ]
CC                              + dj-1[ e1v/e2v dj(trb) ] }
CC    Multiply by the eddy diffusivity coef. and insure lateral bc:
CC      Bilaplacian (laplacian of zlt):
CC         difft = 1/(e1t*e2t) {  di-1[ e2u/e1u di(zlt) ]
CC                              + dj-1[ e1v/e2v dj(zlt) ]  }
CC
CC      Add this trend to the general trend (tra):
CC         (tra) = (tra) + ( difftr )
CC
CC
CC      macro-tasked on tracer slab (jn-loop)
CC
CC
CC   OUTPUT :
CC   ------
CC    tra      : general passive tracer trend increased by the
CC                                horizontal diffusion trend


CC----------------------------------------------------------------------
      USE myalloc
      USE myalloc_mpp
      USE HDF_mem
      USE DIA_mem
        IMPLICIT NONE
CC----------------------------------------------------------------------
CC local declarations
CC ==================


      LOGICAL l1,l2,l3
      INTEGER ji,jj,jk,jn,jv,jf,mytid,ntids,pack_size,jp
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
      ntids = mpi_pack_size
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
C   ... z-coordinates, no vertical scale factors
                  zbtr(ji,jj) = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
                  zeeu(ji,jj,jk) = e2u(ji,jj) / e1u(ji,jj) * umask(ji,jj,jk)
                  zeev(ji,jj,jk) = e1v(ji,jj) / e2v(ji,jj) * vmask(ji,jj,jk)
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

          jarr_hdf_flx=0

             DO jf=1,Fsize
                DO jv=1, dimen_jvhdf3

                   l1 = flx_ridxt(jf,2) .EQ. jarr_hdf(1,jv,2)
                   l2 = flx_ridxt(jf,3) .EQ. jarr_hdf(2,jv,2)
                   l3 = flx_ridxt(jf,4) .EQ. jarr_hdf(3,jv,2)

                   IF ( l1 .AND. l2 .AND. l3) THEN
                      jarr_hdf_flx(jv)= jf
                   END IF

                END DO
             END DO
       ENDIF

C tracer slab
C =============

      TRACER_LOOP: DO  jn = 1, jptra, ntids

C 1. Laplacian
C ------------

C ... First derivative (gradient)
!$omp  parallel default(none) private(mytid,jv,jk,jj,ji)
!$omp&                        shared(jn,dimen_jvhdf2,jarr_hdf,ztu,zeeu,trb,tmask,ztv,zeev,
!$omp&                               dimen_jvhdf3,zlt,zbtr,trcrat,ahtt)
#ifdef __OPENMP
       mytid = omp_get_thread_num()  ! take the thread ID
#else
      PACK_LOOP1: DO jp=1,ntids
       mytid=jp-1
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
             zlt(ji,jj,jk,mytid+1) = trcrat * ahtt(jk) * zlt(ji,jj,jk,mytid+1)

          END DO

       ENDIF
!$omp  end parallel
#ifdef __OPENMP
#else
      END DO PACK_LOOP1
      mytid =0
#endif
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

#endif

C 2. Bilaplacian
C --------------

C ... third derivative (gradient)
!$omp  parallel default(none) private(mytid,jv,jk,jj,ji,jf)
!$omp&                        shared(jn,dimen_jvhdf2,jarr_hdf,ztu,zeeu,zlt,tmask,ztv,zeev,
!$omp&                               dimen_jvhdf3,zta,zbtr,tra,jarr_hdf_flx,diaflx,Fsize)
#ifdef __OPENMP
       mytid = omp_get_thread_num()  ! take the thread ID
#else
      PACK_LOOP2: DO jp=1,ntids
       mytid=jp-1
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

C ... fourth derivative (divergence) and add to the general tracer trend

          DO jv=1, dimen_jvhdf3

             ji = jarr_hdf(1,jv,2)
             jj = jarr_hdf(2,jv,2)
             jk = jarr_hdf(3,jv,2)
             jf = jarr_hdf_flx(jv)

C   ... horizontal diffusive trends
             zta(mytid+1) = (  ztu(ji,jj,jk,mytid+1) - ztu(ji-1,jj,jk,mytid+1)
     $            + ztv(ji,jj,jk,mytid+1) - ztv(ji,jj-1,jk,mytid+1)  ) * zbtr(ji,jj)
C   ... add it to the general tracer trends
              tra(ji,jj,jk,jn+mytid) = tra(ji,jj,jk,jn+mytid) + zta(mytid+1)

!     Save diffusive fluxes x,y
              IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                 diaflx(jf,jn+mytid,5) = diaflx(jf,jn+mytid,5) + ztu(ji,jj,jk,mytid+1)
                 diaflx(jf,jn+mytid,6) = diaflx(jf,jn+mytid,6) + ztv(ji,jj,jk,mytid+1)
              END IF

         END DO

      ENDIF
!$omp  end parallel
#ifdef __OPENMP
#else
      END DO PACK_LOOP2
      mytid =0
#endif
      
C End of slab
C ===========

        END DO TRACER_LOOP

       trcbilaphdfparttime = MPI_WTIME() - trcbilaphdfparttime
       trcbilaphdftottime = trcbilaphdftottime + trcbilaphdfparttime


