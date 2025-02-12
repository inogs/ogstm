      SUBROUTINE trcsed
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcsed
!!!                     *******************
!!!
!!!  PURPOSE :
!!!  ---------
!!!     compute the now trend due to the vertical sedimentation of
!!!     detritus and add it to the general trend of detritus equations.
!!!
!!!
!!   METHOD :
!!   -------
!!      this ROUTINE compute not exactly the advection but the
!!      transport term, i.e.  dz(wt) and dz(ws)., dz(wtr)
!!      using an upstream scheme
!!
!!    the now vertical advection of tracers is given by:
!!
!!          dz(trn wn) = 1/bt dk+1( e1t e2t vsed (trn) )
!!
!!    add this trend now to the general trend of tracer (ta,sa,tra):
!!
!!                     tra = tra + dz(trn wn)
!!
!!      IF 'key_trc_diabio' key is activated, the now vertical advection
!!      trend of passive tracers is saved for futher diagnostics.
!!
!!    multitasked on vertical slab (jj-loop)
!!
!!
!!
!!   OUTPUT :
!!   ------
!!
!!   WORKSPACE :
!!   ---------
!!    local
!!    ze1e2w, ze3tr, ztra
!!      COMMON
!!
!!   EXTERNAL :                   no
!!   --------
!!
!!   REFERENCES :                 no
!!   ----------

       USE myalloc
       USE BIO_mem
       USE SED_mem
       USE DIA_mem
       IMPLICIT NONE


!!----------------------------------------------------------------------
!! local declarations
!! ==================


#ifdef key_trc_bfm

      LOGICAL :: l1,l2,l3
      LOGICAL,save :: first = .true.
      INTEGER :: ji,jj,jk,jv,jf,js,ntx
      INTEGER :: bottom,queue
      double precision :: ze3tr,d2s
! omp variables
    
!!----------------------------------------------------------------------
!! statement functions
!! ===================


      queue=1

      d2s=1./3600./24.  ! speed from (m/day) to  (m/s)

      IF (first) THEN ! initialization phase
           DO  ji = 2,jpim1
         DO jj = 2,jpjm1
              IF(tmask(1,jj,ji) .NE. 0) THEN
                 dimen_jvsed = dimen_jvsed + 1
                 jarr_sed(1,dimen_jvsed) = jj
                 jarr_sed(2,dimen_jvsed) = ji
              ENDIF
            END DO
         END DO
! Cross matrix between transect and local optimized indexing
         jarr_sed_flx=0

         DO jf=1,Fsize
            DO jv=1, dimen_jvsed
               DO jk=1,jpk
                  l1 = flx_ridxt(jf,4) .EQ. jarr_sed(2,jv)
                  l2 = flx_ridxt(jf,3) .EQ. jarr_sed(1,jv)
                  l3 = flx_ridxt(jf,2) .EQ. jk
                  IF ( l1 .AND. l2 .AND. l3) THEN
                     jarr_sed_flx(jk,jv)= jf
                  END IF
               END DO
            END DO
         END DO
         call myalloc_SED_ztra_zwork()
         !$acc update device(jarr_sed,jarr_sed_flx)
         first=.false.
      ENDIF ! End initialization phase (once at the beginning)


! vertical slab
! =============

      !      if( mytid + jv <=  dimen_jvsed) then
      ! 1. sedimentation of detritus  : upstream scheme
      ! -----------------------------------------------
      ! 1.1 initialisation needed for bottom and surface value

      !$acc parallel loop gang vector collapse(3) default(present) async(queue)
      DO jv=1,dimen_jvsed
         DO js = 1, nsed
            DO  jk = 1,jpk
               ntx=jv
               ji = jarr_sed(2,jv)
               jj = jarr_sed(1,jv)
               zwork(jk,js,ntx) = 0.

            END DO
         END DO
      END DO

      ! 1.2 tracer flux at w-point: we use -vsed (downward flux)
      ! with simplification : no e1*e2

      !                Particulate
      !$acc parallel loop gang vector collapse(3) default(present) async(queue)
      DO jv=1,dimen_jvsed
         DO js =1,4
            DO  jk = 2,jpkm1
               ntx=jv
               ji = jarr_sed(2,jv)
               jj = jarr_sed(1,jv)
               zwork(jk,js,ntx) = -vsed * trn(jk-1,jj,ji, sed_idx(js))
            END DO
         END DO
      END DO

      !                Diatoms
      !$acc parallel loop gang vector collapse(3) default(present) async(queue)
      DO jv=1,dimen_jvsed
         DO js =5,9
            DO  jk = 2,jpkm1
               ntx=jv
               ji = jarr_sed(2,jv)
               jj = jarr_sed(1,jv)
               zwork(jk,js,ntx) = -ogstm_sedipi(jk-1,jj,ji,1) * trn(jk-1,jj,ji, sed_idx(js))
            END DO
         END DO
      END DO

      !                Flagellates
      !$acc parallel loop gang vector collapse(3) default(present) async(queue)
      DO jv=1,dimen_jvsed
         DO js =10,13
            DO  jk = 2,jpkm1
               ntx=jv
               ji = jarr_sed(2,jv)
               jj = jarr_sed(1,jv)
               zwork(jk,js,ntx) = -ogstm_sedipi(jk-1,jj,ji,2) * trn(jk-1,jj,ji, sed_idx(js))
            END DO
         END DO
      END DO

      !                Picophytoplankton
      !$acc parallel loop gang vector collapse(3) default(present) async(queue)
      DO jv=1,dimen_jvsed
         DO js =14,17
            DO  jk = 2,jpkm1
               ntx=jv
               ji = jarr_sed(2,jv)
               jj = jarr_sed(1,jv)
               zwork(jk,js,ntx) = -ogstm_sedipi(jk-1,jj,ji,3) * trn(jk-1,jj,ji, sed_idx(js))
            END DO
         END DO
      END DO

      !                Dinoflagellates
      !$acc parallel loop gang vector collapse(3) default(present) async(queue)
      DO jv=1,dimen_jvsed
         DO js =18,21
            DO  jk = 2,jpkm1
               ntx=jv
               ji = jarr_sed(2,jv)
               jj = jarr_sed(1,jv)
               zwork(jk,js,ntx) = -ogstm_sedipi(jk-1,jj,ji,4) * trn(jk-1,jj,ji, sed_idx(js))
            END DO
         END DO
      END DO

      !                Calcite
      !$acc parallel loop gang vector collapse(3) default(present) async(queue)
      DO jv=1,dimen_jvsed
         DO js =22,22
            DO  jk = 2,jpkm1
               ntx=jv
               ji = jarr_sed(2,jv)
               jj = jarr_sed(1,jv)
               zwork(jk,js,ntx) = - vsedO5c * trn(jk-1,jj,ji, sed_idx(js))
            END DO
         END DO
      END DO

      !$acc parallel loop gang vector collapse(2) default(present) async(queue)
      DO jv=1,dimen_jvsed
         DO js = 1,nsed
            ntx=jv
            ji = jarr_sed(2,jv)
            jj = jarr_sed(1,jv)
            bottom = mbathy(jj,ji) + 1
            zwork(bottom,js,ntx) = bottom_flux * zwork(bottom,js,ntx) ! bottom_flux = 0 -> no flux in the sea floor
         END DO
      END DO

      ! 1.3 tracer flux divergence at t-point added to the general trend

      !$acc parallel loop gang vector collapse(2) default(present) async(queue)
      DO jv=1,dimen_jvsed
         DO  jk = 1,jpkm1
            ntx=jv
            ji = jarr_sed(2,jv)
            jj = jarr_sed(1,jv)

            jf=  jarr_sed_flx(jk,jV)

            ze3tr = 1./e3t(jk,jj,ji)

            DO js =1,nsed
               ztra(js,ntx) = -ze3tr * (zwork(jk,js,ntx) - zwork(jk+1,js,ntx))
               IF ((Fsize .GT. 0) .AND. (jf .GT. 0)) THEN
                  diaflx(4,jf,sed_idx(js)) = diaflx(4, jf,sed_idx(js)) + zwork(jk,js,ntx)*rdt
               ENDIF

!!!  d2s convert speed from (m/day) to  (m/s)
               tra(jk,jj,ji,sed_idx(js)) = tra(jk,jj,ji,sed_idx(js)) + ztra(js,ntx)*d2s
            END DO
         END DO
      END DO
      !$acc end parallel loop
      !$acc wait(queue)

#ifdef key_trc_diabio
      DO jv=1,dimen_jvsed
         ji = jarr_sed(2,jv)
         jj = jarr_sed(1,jv)
         DO  jk = 1,jpkm1
            trbio(jk,jj,ji,8) = ztra
         END DO
      END DO
#endif

#endif ! key_trc_bfm

!!!$omp    end parallel do


!#else

!       no Sedimentation

!#endif

      END SUBROUTINE trcsed
