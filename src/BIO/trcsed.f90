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
      REAL(8) :: ze3tr,d2s
! omp variables
      INTEGER :: mytid, ntids

#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif
!!----------------------------------------------------------------------
!! statement functions
!! ===================




      d2s=1./3600./24.  ! speed from (m/day) to  (m/s)

#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif


      IF (dimen_jvsed .EQ. 0) THEN ! initialization phase
         DO jj = 2,jpjm1
           DO  ji = 2,jpim1
              IF(tmask(1,jj,ji) .NE. 0) THEN
                 dimen_jvsed = dimen_jvsed + 1
                 jarr_sed(1,dimen_jvsed) = ji
                 jarr_sed(2,dimen_jvsed) = jj
              ENDIF
            END DO
         END DO
! Cross matrix between transect and local optimized indexing
         jarr_sed_flx=0

         DO jf=1,Fsize
            DO jv=1, dimen_jvsed
               DO jk=1,jpk
                  l1 = flx_ridxt(jf,2) .EQ. jarr_sed(1,jv)
                  l2 = flx_ridxt(jf,3) .EQ. jarr_sed(2,jv)
                  l3 = flx_ridxt(jf,4) .EQ. jk
                  IF ( l1 .AND. l2 .AND. l3) THEN
                     jarr_sed_flx(jk,jv)= jf
                  END IF
               END DO
            END DO
         END DO

      ENDIF ! End initialization phase (once at the beginning)


! vertical slab
! =============


!!!$omp    parallel do default(none) private(jv,ji,jj,jk,js,mytid,jf,bottom)
!!!$omp&                           shared(dimen_jvsed,jarr_sed,jpk, 
!!!$omp&                                  jpkm1,nsed,zwork,vsed,sediPI,
!!!$omp&                                  trn,sed_idx,ze3tr,e3t,ztra,
!!!$omp&                                  tra,d2s,jarr_sed_flx,Fsize,diaflx,
!!!$omp&                                  mbathy,bottom_flux)

      MAIN_LOOP: DO jv=1,dimen_jvsed

#ifdef __OPENMP1
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
!      if( mytid + jv <=  dimen_jvsed) then
! 1. sedimentation of detritus  : upstream scheme
! -----------------------------------------------
! 1.1 initialisation needed for bottom and surface value

           ji = jarr_sed(1,jv)
           jj = jarr_sed(2,jv)

              DO  jk = 1,jpk

                 DO js = 1, nsed

                    zwork(jk,js,mytid+1) = 0.

                 END DO

              END DO

! 1.2 tracer flux at w-point: we use -vsed (downward flux)
! with simplification : no e1*e2

              DO  jk = 2,jpkm1

!                Particulate
                 DO js =1,4
                    zwork(jk,js,mytid+1) = -vsed * trn(jk-1,jj,ji, sed_idx(js))
                 END DO

!                Diatoms
                 DO js =5,9
                    zwork(jk,js,mytid+1) = -sediPI(jk-1,jj,ji,1) * trn(jk-1,jj,ji, sed_idx(js))
                 END DO

!                Flagellates
                 DO js =10,13
                    zwork(jk,js,mytid+1) = -sediPI(jk-1,jj,ji,2) * trn(jk-1,jj,ji, sed_idx(js))
                 END DO

!                Picophytoplankton
                 DO js =14,17
                    zwork(jk,js,mytid+1) = -sediPI(jk-1,jj,ji,3) * trn(jk-1,jj,ji, sed_idx(js))
                 END DO

!                Dinoflagellates
                 DO js =18,21
                    zwork(jk,js,mytid+1) = -sediPI(jk-1,jj,ji,4) * trn(jk-1,jj,ji, sed_idx(js))
                 END DO

              END DO

               bottom = mbathy(ji,jj) + 1
               zwork(:,bottom,mytid+1) = bottom_flux * zwork(:,bottom,mytid+1) ! bottom_flux = 0 -> no flux in the sea floor

! 1.3 tracer flux divergence at t-point added to the general trend

              DO  jk = 1,jpkm1
                 jf=  jarr_sed_flx(jk,jv)

                 ze3tr = 1./e3t(ji,jj,jk)

                 DO js =1,21
                    ztra(js,mytid+1) = -ze3tr * (zwork(jk,js,mytid+1) - zwork(jk+1,js,mytid+1))
                    IF ((Fsize .GT. 0) .AND. (jf .GT. 0)) THEN
                         diaflx(jf,sed_idx(js),4) = diaflx(jf,sed_idx(js),4) + zwork(jk,js,mytid+1)
                    ENDIF
                 END DO

                 DO js =1,21
!!!  d2s convert speed from (m/day) to  (m/s)
                    tra(ji,jj,jk,sed_idx(js)) = tra(ji,jj,jk,sed_idx(js)) + ztra(js,mytid+1)*d2s
                 END DO

#ifdef key_trc_diabio
                  trbio(ji,jj,jk,8) = ztra
#endif


                END DO
#ENDIF


      END DO MAIN_LOOP

!!!$omp    end parallel do


!#else

!       no Sedimentation

!#endif

      END SUBROUTINE trcsed
