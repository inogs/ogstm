      SUBROUTINE trchdf
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trchdf
!!!                     ******************
!!!
!!!  PURPOSE :
!!!  ---------
!!!     compute the before horizontal passive tracer diffusive trend
!!!     add it to the general trend of tracer equation.
!!!
!!!     harmonic operator (default option)
!!!    -----------------
!!!       z-coordinates (default option):
!!!             default key : operator acting along model level surfaces
!!!                (which are geopotential surfaces).
!!!        'key_trc_hdfiso' : operator acting along neutral surfaces
!!!                (rotation of the diffusive tensor) plus an eddy
!!!            induced advection if 'key_trc_hdfeiv' defined.
!!!       s-coordinates ('key_s_coord'):
!!!        default key : operator acting along model level surfaces
!!!            (which are not geopotential surfaces).
!!!        'key_trc_hdfgeop' : operator acting along geopotential
!!!            surfaces (rotation of the diffusive tensor).
!!!        'key_trc_hdfiso' : operator acting along neutral surfaces
!!!                     (rotation of the diffusive tensor) plus an eddy
!!!                     induced advection if 'key_trc_hdfeiv' defined.
!!!
!!!    biharmonic operator ('key_trc_hdfbilap')
!!!    -------------------
!!!       z-coordinates (default option):
!!!             default key : operator acting along model level surfaces
!!!                     (which are geopotential surfaces).
!!!        s-coordinates ('key_s_coord'):
!!!             default key : operator acting along model level surfaces
!!!                     (which are not geopotential surfaces).
!!!             'key_trc_hdfgeop' : operator acting along geopotential
!!!                     surfaces (rotation of the diffusive tensor).
!!!
!!!
!!!     notes: 
!!!            IF you want no horizontal diffusion you have to switch
!!!            the boolean lhdf to .FALSE.
!!!
!!!     MODIFICATION : 00-11 (M.A.Foujols E. Kestenare)
!!!              differents keys for passive tracer
!!!

!!!        TRCHDF.BILAPLACIAN
!!!      **********************
!!!
!!   define key : 'key_trc_hdfbilap'   but not key_trc_hdfgeop
!!   ==========
!!
!!
!!   METHOD :
!!   -------
!!      4th order diffusive operator along model level surfaces evalu-
!!    ated using before fields (forward time scheme). The horizontal
!!    diffusive trends of passive tracer is given by:
!!    Multiply by the eddy diffusivity coef. and insure lateral bc:
!!      Bilaplacian (laplacian of zlt):
!!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
!!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
!!
!!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
!!      Laplacian of trb:
!!         zlt   = 1/(e1t*e2t) {  di-1[ e2u/e1u di(trb) ]
!!                              + dj-1[ e1v/e2v dj(trb) ] }
!!    Multiply by the eddy diffusivity coef. and insure lateral bc:
!!      Bilaplacian (laplacian of zlt):
!!         difft = 1/(e1t*e2t) {  di-1[ e2u/e1u di(zlt) ]
!!                              + dj-1[ e1v/e2v dj(zlt) ]  }
!!
!!      Add this trend to the general trend (tra):
!!         (tra) = (tra) + ( difftr )
!!
!!
!!      macro-tasked on tracer slab (jn-loop)
!!
!!
!!   OUTPUT :
!!   ------
!!    tra      : general passive tracer trend increased by the
!!                                horizontal diffusion trend


!!----------------------------------------------------------------------
      USE myalloc
      USE ogstm_mpi_module
      ! epascolo USE myalloc_mpp
      USE HDF_mem
      USE DIA_mem
      use mpi
        IMPLICIT NONE
!!----------------------------------------------------------------------
!! local declarations
!! ==================


      LOGICAL :: l1,l2,l3
      INTEGER :: jk,jj,ji,jn,jv,jf,mytid,ntids,pack_size,jp
      INTEGER :: myji,myjj
      INTEGER :: locsum,jklef,jjlef,jilef,jkrig,jjrig,jirig

!!----------------------------------------------------------------------
!! statement functions
!! ===================

       trcbilaphdfparttime = MPI_WTIME()

!! Define auxiliary matrix


       IF (dimen_jvhdf1 .EQ. 0) THEN
                DO ji = 1,jpi
             DO jj = 1,jpj
          DO jk = 1,jpk
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

                         DO myji=ji+jilef, ji+jirig
                      DO myjj=jj+jjlef, jj+jjrig
                            locsum = locsum + tmask(jk,myjj,myji)
                         END DO
                      END DO

                   if(locsum .NE. 0) then
                      dimen_jvhdf1 = dimen_jvhdf1 + 1
                      hdfmask(jk,jj,ji) = 1
                   else
                      hdfmask(jk,jj,ji) = 0
                   endif
                END DO
             END DO
          END DO

!! 0. Initialization of metric arrays (for z- or s-coordinates)
!! ----------------------------------
               DO ji = 1, jpim1
             DO jj = 1, jpjm1
          DO jk=1,jpkm1
!!   ... z-coordinates, no vertical scale factors
                  zbtr(jk,jj,ji) = 1. / ( e1t(jj,ji)*e2t(jj,ji)*e3t(jk,jj,ji) )
                  zeeu(jk,jj,ji) = e2u(jj,ji)*e3u(jk,jj,ji) / e1u(jj,ji) * umask(jk,jj,ji)
                  zeev(jk,jj,ji) = e1v(jj,ji)*e3v(jk,jj,ji) / e2v(jj,ji) * vmask(jk,jj,ji)
               END DO
             END DO
          END DO

       ENDIF

       IF (dimen_jvhdf2 .EQ. 0) THEN
                DO  ji = 1,jpim1
             DO jj = 1,jpjm1
          DO jk = 1,jpkm1
                   IF(hdfmask(jk,jj,ji) .NE. 0) THEN
                      dimen_jvhdf2 = dimen_jvhdf2 + 1
                      jarr_hdf(1,dimen_jvhdf2,1) = jk
                      jarr_hdf(2,dimen_jvhdf2,1) = jj
                      jarr_hdf(3,dimen_jvhdf2,1) = ji
                   ENDIF
                END DO
             END DO
          END DO
       ENDIF

!!       dimen_jvhdf3=0

       IF (dimen_jvhdf3 .EQ. 0) THEN
                DO  ji = 2,jpim1
             DO jj = 2,jpjm1
          DO jk = 1,jpkm1
                   IF(hdfmask(jk,jj,ji) .NE. 0) THEN
                      dimen_jvhdf3 = dimen_jvhdf3 + 1
                      jarr_hdf(1,dimen_jvhdf3,2) = jk
                      jarr_hdf(2,dimen_jvhdf3,2) = jj
                      jarr_hdf(3,dimen_jvhdf3,2) = ji
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

!! tracer slab
!! =============

      TRACER_LOOP: DO  jn = 1, jptra

!! 1. Laplacian
!! ------------

!! ... First derivative (gradient)
!!!&omp  parallel default(none) private(mytid,jv,jk,jj,ji)
!!!&omp&                        shared(jn,dimen_jvhdf2,jarr_hdf,ztu,zeeu,trb,tmask,ztv,zeev,
!!!&omp&                               dimen_jvhdf3,zlt,zbtr,trcrat,ahtt)

          DO jv=1, dimen_jvhdf2

             ji = jarr_hdf(3,jv,1)
             jj = jarr_hdf(2,jv,1)
             jk = jarr_hdf(1,jv,1)

             ztu(jk,jj,ji, 1) = zeeu(jk,jj,ji) * &
     &          ( trb(jk,jj,ji+1,jn ) - trb(jk,jj,ji,jn ) )* &
     &          tmask(jk,jj,ji+1) * tmask(jk,jj,ji)

             ztv(jk,jj,ji, 1) = zeev(jk,jj,ji) * &
     &          ( trb(jk,jj+1,ji,jn ) - trb(jk,jj,ji,jn ) )* &
     &          tmask(jk,jj+1,ji) * tmask(jk,jj,ji)

          END DO
!!
!! ... Second derivative (divergence)
          DO jv=1, dimen_jvhdf3

             ji = jarr_hdf(3,jv,2)
             jj = jarr_hdf(2,jv,2)
             jk = jarr_hdf(1,jv,2)

             zlt(jk,jj,ji, 1) = (  ztu(jk,jj,ji, 1) - ztu(jk,jj,ji-1, 1) &
     &             + ztv(jk,jj,ji, 1) - ztv(jk,jj-1,ji, 1)  ) * zbtr(jk,jj,ji)
!! ... Multiply by the eddy diffusivity coefficient
             zlt(jk,jj,ji, 1) = trcrat * ahtt(jk) * zlt(jk,jj,ji, 1)

          END DO

       

!!
!!
!! ... Lateral boundary conditions on the laplacian (zlt,zls)

! epascolo mpi comment
! #ifdef key_mpp
! !!
! !!   ... Mpp : export boundary values to neighboring processors
! !!
!        IF( ntids - 1 + jn <= jptra ) THEN
!           pack_size = ntids
!        ELSE
!           pack_size = ntids - (ntids - 1 + jn - jptra)
!        END IF

!        CALL mpplnk_my(zlt(:,:,:,:), pack_size,1,1)

! #else

!         DO itid = 1, ntids

!            IF( itid - 1 + jn <= jptra ) THEN

               CALL lbc( zlt(:,:,:,1), 1, 1, 1, 1, jpk, 1 )

!            END IF
!         END DO

! #endif

!! 2. Bilaplacian
!! --------------

!! ... third derivative (gradient)
!!!&omp  parallel default(none) private(mytid,jv,jk,jj,ji,jf)
!!!&omp&                        shared(jn,dimen_jvhdf2,jarr_hdf,ztu,zeeu,zlt,tmask,ztv,zeev,
!!!&omp&                               dimen_jvhdf3,zta,zbtr,tra,jarr_hdf_flx,diaflx,Fsize)


          DO jv=1, dimen_jvhdf2

             ji = jarr_hdf(3,jv,1)
             jj = jarr_hdf(2,jv,1)
             jk = jarr_hdf(1,jv,1)


             ztu(jk,jj,ji, 1) = zeeu(jk,jj,ji) * &
     &        ( zlt(jk,jj,ji+1, 1) - zlt(jk,jj,ji, 1) ) * &
     &        tmask(jk,jj,ji+1) * tmask(jk,jj,ji)
             ztv(jk,jj,ji, 1) = zeev(jk,jj,ji) * & 
     &        ( zlt(jk,jj+1,ji, 1) - zlt(jk,jj,ji, 1) ) * &
     &        tmask(jk,jj+1,ji) * tmask(jk,jj,ji)

          END DO

!! ... fourth derivative (divergence) and add to the general tracer trend

          DO jv=1, dimen_jvhdf3

             ji = jarr_hdf(3,jv,2)
             jj = jarr_hdf(2,jv,2)
             jk = jarr_hdf(1,jv,2)
             jf = jarr_hdf_flx(jv)

!!   ... horizontal diffusive trends
             zta( 1) = (  ztu(jk,jj,ji, 1) - ztu(jk,jj,ji-1, 1) &
     &            + ztv(jk,jj,ji, 1) - ztv(jk,jj-1,ji, 1)  ) * zbtr(jk,jj,ji)
!!   ... add it to the general tracer trends
              tra(jk,jj,ji,jn ) = tra(jk,jj,ji,jn ) + zta( 1)

!     Save diffusive fluxes x,y
              IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                 diaflx(jf,jn ,5) = diaflx(jf,jn ,5) + ztu(jk,jj,ji, 1)
                 diaflx(jf,jn ,6) = diaflx(jf,jn ,6) + ztv(jk,jj,ji, 1)
              END IF

         END DO

 
      
!! End of slab
!! ===========

        END DO TRACER_LOOP

       trcbilaphdfparttime = MPI_WTIME() - trcbilaphdfparttime
       trcbilaphdftottime = trcbilaphdftottime + trcbilaphdfparttime

      END SUBROUTINE trchdf
