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
      INTEGER :: ji,jj,jk,jv,jf,js
      INTEGER :: bottom
      double precision :: ze3tr,d2s
! omp variables
    
!!----------------------------------------------------------------------
!! statement functions
!! ===================




      d2s=1./3600./24.  ! speed from (m/day) to  (m/s)

      IF (dimen_jvsed .EQ. 0) THEN ! initialization phase
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

      ENDIF ! End initialization phase (once at the beginning)


! vertical slab
! =============



      MAIN_LOOP: DO jv=1,dimen_jvsed

!      if( mytid + jv <=  dimen_jvsed) then
! 1. sedimentation of detritus  : upstream scheme
! -----------------------------------------------
! 1.1 initialisation needed for bottom and surface value

           ji = jarr_sed(2,jv)
           jj = jarr_sed(1,jv)


                 DO js = 1, nsed
              DO  jk = 1,jpk

                    zwork(jk,js,1) = 0.

                 END DO
              END DO

! 1.2 tracer flux at w-point: we use -vsed (downward flux)
! with simplification : no e1*e2

             

!                Particulate
              DO js =1,4
                 DO  jk = 2,jpkm1
                    zwork(jk,js,1) = -vsed * trn(jk-1,jj,ji, sed_idx(js))
                 END DO
              END DO
!                Diatoms
              DO js =5,9
                 DO  jk = 2,jpkm1
                    zwork(jk,js,1) = -ogstm_sedipi(jk-1,jj,ji,1) * trn(jk-1,jj,ji, sed_idx(js))
                 END DO
              END DO
!                Flagellates
              DO js =10,13
                 DO  jk = 2,jpkm1
                    zwork(jk,js,1) = -ogstm_sedipi(jk-1,jj,ji,2) * trn(jk-1,jj,ji, sed_idx(js))
                 END DO
              END DO
!                Picophytoplankton
              DO js =14,17
                 DO  jk = 2,jpkm1
                    zwork(jk,js,1) = -ogstm_sedipi(jk-1,jj,ji,3) * trn(jk-1,jj,ji, sed_idx(js))
                 END DO
              END DO
!                Dinoflagellates
              DO js =18,21
                 DO  jk = 2,jpkm1
                    zwork(jk,js,1) = -ogstm_sedipi(jk-1,jj,ji,4) * trn(jk-1,jj,ji, sed_idx(js))
                 END DO
              END DO

#ifndef BFMv2
!                Calcite
              DO js =22,22
                 DO  jk = 2,jpkm1
                    zwork(jk,js,1) = - vsedO5c * trn(jk-1,jj,ji, sed_idx(js))
                 END DO
              END DO
#endif
               bottom = mbathy(jj,ji) + 1
               zwork(bottom,:,1) = bottom_flux * zwork(bottom,:,1) ! bottom_flux = 0 -> no flux in the sea floor

! 1.3 tracer flux divergence at t-point added to the general trend

              DO  jk = 1,jpkm1
                  jf=  jarr_sed_flx(jk,jV)

                 ze3tr = 1./e3t(jk,jj,ji)

                 DO js =1,nsed
                    ztra(js,1) = -ze3tr * (zwork(jk,js,1) - zwork(jk+1,js,1))
                    IF ((Fsize .GT. 0) .AND. (jf .GT. 0)) THEN
                         diaflx(4,jf,sed_idx(js)) = diaflx(4, jf,sed_idx(js)) + zwork(jk,js,1)*rdt
                    ENDIF
                 END DO

                 DO js =1,nsed
!!!  d2s convert speed from (m/day) to  (m/s)
                    tra(jk,jj,ji,sed_idx(js)) = tra(jk,jj,ji,sed_idx(js)) + ztra(js,1)*d2s
                 END DO

#ifdef key_trc_diabio
                  trbio(jk,jj,ji,8) = ztra
#endif

              END DO
#endif


      END DO MAIN_LOOP

!!!$omp    end parallel do


!#else

!       no Sedimentation

!#endif

      END SUBROUTINE trcsed
