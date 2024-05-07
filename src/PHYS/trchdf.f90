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

      USE DIA_mem
      use mpi

      use simple_timer

      IMPLICIT NONE
!!----------------------------------------------------------------------
!! local declarations
!! ==================


      INTEGER(KIND=1), allocatable,dimension(:,:,:),save :: hdfmask
      double precision, allocatable,dimension(:,:,:),save :: zeeu, zeev, zbtr
      INTEGER,save :: dimen_jvhdf1=0

      LOGICAL :: l1,l2,l3
      INTEGER :: jk,jj,ji,jn,jv,jf,jp
      INTEGER :: myji,myjj
      INTEGER :: locsum,jklef,jjlef,jilef,jkrig,jjrig,jirig
      !INTEGER, allocatable :: jarr_hdf(:,:,:),jarr_hdf_flx(:)
      double precision, allocatable,dimension(:,:,:) :: zlt, ztu, ztv
      integer :: queue
      logical :: use_gpu
!!----------------------------------------------------------------------
!! statement functions
!! ===================

      call tstart("trchdf_1")

#ifdef _OPENACC
      use_gpu=.true.
#else
      use_gpu=.false.
#endif

    !  #include "BFM_var_list.h"
       trcbilaphdfparttime = MPI_WTIME()

!! Define auxiliary matrix


         !   print *,"val ",dimen_jvhdf1
       IF (.not.hdf_initialized ) THEN

        
        !    print *,"ENTRATO val ",dimen_jvhdf1

             allocate(hdfmask(jpk,jpj,jpi   ))    
             hdfmask      = huge(hdfmask(1,1,1))
             allocate(zeeu   (jpk,jpj,jpi      )) 
             zeeu         = huge(zeeu(1,1,1))
             allocate(zeev   (jpk,jpj,jpi      )) 
             zeev         = huge(zeev(1,1,1))
             allocate(zbtr   (jpk,jpj,jpi      )) 
             zbtr         = huge(zbtr(1,1,1))

             !$acc enter data create(hdfmask,zeeu,zeev,zbtr)

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

        queue=1
        !$acc update device(hdfmask)
        hdf_initialized=.true.

        ENDIF

        allocate(zlt    (jpk,jpj,jpi))
        allocate(ztu    (jpk,jpj,jpi))
        allocate(ztv    (jpk,jpj,jpi))
        !$acc enter data create(zlt,ztu,ztv)

!  Metric arrays calculated out of the initialisation phase(for z- or s-coordinates)
! !! ----------------------------------
                !$acc parallel loop gang vector collapse(3) default(present) async(queue)
                DO ji = 1, jpi
              DO jj = 1, jpj
           DO jk=1,jpk
! !!   ... z-coordinates, no vertical scale factors
                   zbtr(jk,jj,ji) = 1. / ( e1t(jj,ji)*e2t(jj,ji)*e3t(jk,jj,ji) )
                   zeeu(jk,jj,ji) = e2u(jj,ji)*e3u(jk,jj,ji) / e1u(jj,ji) * umask(jk,jj,ji)
                   zeev(jk,jj,ji) = e1v(jj,ji)*e3v(jk,jj,ji) / e2v(jj,ji) * vmask(jk,jj,ji)
                END DO
              END DO
           END DO
           !$acc end parallel loop



     
!! tracer slab
!! =============

      call tstop("trchdf_1")
      call tstart("trchdf_tracer")

! $omp  taskloop default(none) private(jv,jk,jj,ji) &
! $omp  private(jn,ztu,ztv,zlt) firstprivate(jpi,jpj,jpk,trcrat) &
! $omp  shared(zeeu,trb,tmask,zeev) &
! $omp  shared(hdfmask,zbtr,ahtt,tra)         
      TRACER_LOOP: DO  jn = 1, jptra


             !$acc kernels default(present) async(queue)
             zlt = 0.
             ztu = 0.
             ztv = 0.
             !$acc end kernels
            
!! 1. Laplacian
!! ------------

!! ... First derivative (gradient)

      !     DO jv=1, dimen_jvhdf2ji

            !  ji = jarr_hdf(3,jv,1a)
            !  jj = jarr_hdf(2,jv,1a)
            !  jk = jarr_hdf(1,jv,1a)
            
            ! $OMP TASK default(shared) private(ji,jj,jk)
                  !$acc parallel loop gang vector collapse(3) default(present) async(queue)
                  DO ji = 1,jpi-1
              DO jj = 1,jpj-1
           DO jk = 1,jpk
            !dir$ vector aligned
             ztu(jk,jj,ji) = zeeu(jk,jj,ji) * (trb(jk,jj,ji+1,jn) - trb(jk,jj,ji,jn)) * tmask(jk,jj,ji+1) * tmask(jk,jj,ji) * &
                 hdfmask(jk,jj,ji)

             !ztv(jk,jj,ji)  = zeev(jk,jj,ji) * ( trb(jk,jj+1,ji,jn ) - trb(jk,jj,ji,jn ) )*tmask(jk,jj+1,ji) * tmask(jk,jj,ji)

             END DO
            END DO
          END DO
          ! $OMP END TASK
          !$acc end parallel loop

                ! $OMP TASK default(shared) private(ji,jj,jk)
                !$acc parallel loop gang vector collapse(3) default(present) async(queue)
                DO ji = 1,jpi-1
              DO jj = 1,jpj-1
           DO jk = 1,jpk
             !dir$ vector aligned
             !ztu(jk,jj,ji)  = zeeu(jk,jj,ji) * ( trb(jk,jj,ji+1,jn ) - trb(jk,jj,ji,jn ) )*tmask(jk,jj,ji+1) * tmask(jk,jj,ji)

             ztv(jk,jj,ji) = zeev(jk,jj,ji) * (trb(jk,jj+1,ji,jn) - trb(jk,jj,ji,jn)) * tmask(jk,jj+1,ji) * tmask(jk,jj,ji) * &
                 hdfmask(jk,jj,ji)

             END DO
            END DO
          END DO
          !$acc end parallel loop
          ! $OMP END TASK

          ! $OMP TASKWAIT

!!
!! ... Second derivative (divergence)
      !     DO jv=1, dimen_jvhdf3* tmask(jk,jj,ji)


      !        ji = jarr_hdf(3,jv,2)
      !        jj = jarr_hdf(2,jv,2)
      !        jk = jarr_hdf(1,jv,2)
                  !$acc parallel loop gang vector collapse(3) default(present) async(queue)
                  DO ji = 2,jpi-1
              DO jj = 2,jpj-1
           DO jk = 1,jpk
             !dir$ vector aligned
             zlt(jk,jj,ji) = trcrat * ahtt(jk) * (ztu(jk,jj,ji) - ztu(jk,jj,ji-1) + ztv(jk,jj,ji) - ztv(jk,jj-1,ji)) * &
                 zbtr(jk,jj,ji) * hdfmask(jk,jj,ji)
 
!! ... Multiply by the eddy diffusivity coefficient
             !zlt(jk,jj,ji)  = trcrat * ahtt(jk) * zlt(jk,jj,ji) 
             END DO
            END DO
          END DO
          !$acc end parallel loop



!!
!!
!! ... Lateral boundary conditions on the laplacian (zlt,zls)


        !$acc wait(queue)
#ifdef key_mpp
!  ... Mpp : export boundary values to neighboring processors
        CALL mpplnk_my(zlt,gpu=use_gpu)
#else
#error
               CALL lbc( zlt(:,:,:), 1, 1, 1, 1, jpk, 1 )
#endif

!! 2. Bilaplacian
!! --------------

!! ... third derivative (gradient)
!!!&omp  parallel default(none) private(mytid,jv,jk,jj,ji,jf)
!!!&omp&                        shared(jn,dimen_jvhdf2,jarr_hdf,ztu,zeeu,zlt,tmask,ztv,zeev,
!!!&omp&                               dimen_jvhdf3,zta,zbtr,tra,jarr_hdf_flx,diaflx,Fsize)

            ! $OMP TASK default(shared) private(ji,jj,jk)
                  !$acc parallel loop gang vector collapse(3) default(present) async(queue)
                  DO ji = 1,jpi-1
              DO jj = 1,jpj-1
           DO jk = 1,jpk
             !dir$ vector aligned
             ztu(jk,jj,ji) = zeeu(jk,jj,ji) * (zlt(jk,jj,ji+1) - zlt(jk,jj,ji)) * tmask(jk,jj,ji+1) * tmask(jk,jj,ji) * &
                 hdfmask(jk,jj,ji)

              END DO
            END DO
          END DO
          !$acc end parallel loop
          ! $OMP END TASK

          ! $OMP TASK default(shared) private(ji,jj,jk)
               !$acc parallel loop gang vector collapse(3) default(present) async(queue)
               DO ji = 1,jpi-1
              DO jj = 1,jpj-1
           DO jk = 1,jpk
              !dir$ vector aligned
      !       ztu(jk,jj,ji)  = zeeu(jk,jj,ji) * ( zlt(jk,jj,ji+1)  - zlt(jk,jj,ji)  ) * tmask(jk,jj,ji+1) * tmask(jk,jj,ji)
             ztv(jk,jj,ji) = zeev(jk,jj,ji) * (zlt(jk,jj+1,ji) - zlt(jk,jj,ji)) * tmask(jk,jj+1,ji) * tmask(jk,jj,ji) * &
                 hdfmask(jk,jj,ji)

              END DO
            END DO
          END DO
          !$acc end parallel loop
          ! $OMP END TASK

!! ... fourth derivative (divergence) and add to the general tracer trend

      !     DO jv=1, dimen_jvhdf3

      !        ji = jarr_hdf(3,jv,2)
      !        jj = jarr_hdf(2,jv,2)
      !        jk = jarr_hdf(1,jv,2)
      !        jf = jarr_hdf_flx(jv)
          
      ! $OMP TASKWAIT
                  !$acc parallel loop gang vector collapse(3) default(present) async(queue)
                  DO ji = 2,jpim1
              DO jj = 2,jpjm1
           DO jk = 1,jpk
!!   ... horizontal diffusive trends
             !dir$ vector aligned
             tra(jk,jj,ji,jn) = tra(jk,jj,ji,jn) + (ztu(jk,jj,ji) - ztu(jk,jj,ji-1) + ztv(jk,jj,ji) - ztv(jk,jj-1,ji)) * &
                 zbtr(jk,jj,ji) * hdfmask(jk,jj,ji)

!!   ... add it to the general tracer trends
              !tra(jk,jj,ji,jn ) = tra(jk,jj,ji,jn ) + zta
               END DO
            END DO
         END DO
         !$acc end parallel loop

         !$acc parallel loop gang vector default(present) async(queue)
         DO jf=1,Fsize
             jk = flx_ridxt(jf,2)
             jj = flx_ridxt(jf,3)
             ji = flx_ridxt(jf,4)

             diaflx(5,jf, jn) = diaflx(5,jf, jn) - ztu(jk,jj,ji)*rdt ! in diffusion we invert face orientation
             diaflx(6,jf, jn) = diaflx(6,jf, jn) - ztv(jk,jj,ji)*rdt

          ENDDO
          !$acc end parallel loop

 
      
!! End of slab
!! ===========

        END DO TRACER_LOOP
        ! $OMP END TASKLOOP
        !$acc wait(queue)

      call tstop("trchdf_tracer")
      call tstart("trchdf_2")

        ! deallocate(hdfmask)
        ! deallocate(zeeu)
        ! deallocate(zeev)
        ! deallocate(zbtr)

        !$acc exit data delete(zlt,ztu,ztv)
        deallocate(zlt)
        deallocate(ztu)
        deallocate(ztv)


       trcbilaphdfparttime = MPI_WTIME() - trcbilaphdfparttime
       trcbilaphdftottime = trcbilaphdftottime + trcbilaphdfparttime

      call tstop("trchdf_2")

      END SUBROUTINE trchdf
