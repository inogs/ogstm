      SUBROUTINE trcsed
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcsed
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     compute the now trend due to the vertical sedimentation of
CCC     detritus and add it to the general trend of detritus equations.
CCC
CCC
CC   METHOD :
CC   -------
CC      this ROUTINE compute not exactly the advection but the
CC      transport term, i.e.  dz(wt) and dz(ws)., dz(wtr)
CC      using an upstream scheme
CC
CC    the now vertical advection of tracers is given by:
CC
CC          dz(trn wn) = 1/bt dk+1( e1t e2t vsed (trn) )
CC
CC    add this trend now to the general trend of tracer (ta,sa,tra):
CC
CC                     tra = tra + dz(trn wn)
CC
CC      IF 'key_trc_diabio' key is activated, the now vertical advection
CC      trend of passive tracers is saved for futher diagnostics.
CC
CC    multitasked on vertical slab (jj-loop)
CC
CC
CC
CC   OUTPUT :
CC   ------
CC
CC   WORKSPACE :
CC   ---------
CC    local
CC    ze1e2w, ze3tr, ztra
CC      COMMON
CC
CC   EXTERNAL :                   no
CC   --------
CC
CC   REFERENCES :                 no
CC   ----------

       USE myalloc
       USE BIO_mem
       USE SED_mem
       USE DIA_mem
       IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================


#   if defined key_trc_bfm

      LOGICAL :: l1,l2,l3
      INTEGER :: ji,jj,jk,jv,jf,js
      INTEGER :: bottom
      REAL(8) :: ze3tr,d2s
! omp variables
      INTEGER :: mytid, ntids

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif
CC----------------------------------------------------------------------
CC statement functions
CC ===================




      d2s=1./3600./24.  ! speed from (m/day) to  (m/s)

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif


      IF (dimen_jvsed .EQ. 0) THEN ! initialization phase
         DO jj = 2,jpjm1
           DO  ji = 2,jpim1
              IF(tmask(ji,jj,1) .NE. 0) THEN
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
                     jarr_sed_flx(jv,jk)= jf
                  END IF
               END DO
            END DO
         END DO

      ENDIF ! End initialization phase (once at the beginning)

C
C
C vertical slab
C =============
C
      MAIN_LOOP: DO jv=1,dimen_jvsed,ntids

!$omp    parallel default(none) private(ji,jj,jk,js,mytid,jf,bottom)
!$omp&                           shared(jv,dimen_jvsed,jarr_sed,
!$omp&                                  jpkm1,nsed,zwork,vsed,sediPI,
!$omp&                                  trn,sed_idx,ze3tr,e3t,ztra,
!$omp&                                  tra,d2s,jarr_sed_flx,Fsize,diaflx,
!$omp&                                  mbathy,bottom_flux)

#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
      if( mytid + jv <=  dimen_jvsed) then
C 1. sedimentation of detritus  : upstream scheme
C -----------------------------------------------
C 1.1 initialisation needed for bottom and surface value
C
           ji = jarr_sed(1,jv+mytid)
           jj = jarr_sed(2,jv+mytid)

              DO  jk = 1,jpk

                 DO js = 1, nsed

                    zwork(js,jk,mytid+1) = 0.

                 END DO

              END DO
C
C 1.2 tracer flux at w-point: we use -vsed (downward flux)
C with simplification : no e1*e2
C
              DO  jk = 2,jpkm1

C                Particulate
                 DO js =1,4
                    zwork(js,jk,mytid+1) = -vsed * trn(ji,jj,jk - 1, sed_idx(js))
                 END DO

C                Diatoms
                 DO js =5,9
                    zwork(js,jk,mytid+1) = -sediPI(ji,jj,jk - 1,1) * trn(ji,jj,jk - 1, sed_idx(js))
                 END DO

C                Flagellates
                 DO js =10,13
                    zwork(js,jk,mytid+1) = -sediPI(ji,jj,jk - 1,2) * trn(ji,jj,jk - 1, sed_idx(js))
                 END DO

C                Picophytoplankton
                 DO js =14,17
                    zwork(js,jk,mytid+1) = -sediPI(ji,jj,jk - 1,3) * trn(ji,jj,jk - 1, sed_idx(js))
                 END DO

C                Dinoflagellates
                 DO js =18,21
                    zwork(js,jk,mytid+1) = -sediPI(ji,jj,jk - 1,4) * trn(ji,jj,jk - 1, sed_idx(js))
                 END DO

              END DO

               bottom = mbathy(ji,jj) + 1
               zwork(:,bottom,mytid+1) = bottom_flux * zwork(:,bottom,mytid+1) ! bottom_flux = 0 -> no flux in the sea floor

C 1.3 tracer flux divergence at t-point added to the general trend

              DO  jk = 1,jpkm1
                 jf=  jarr_sed_flx(jv+mytid,jk)

                 ze3tr = 1./e3t(ji,jj,jk)

                 DO js =1,21
                    ztra(js,mytid+1) = -ze3tr * (zwork(js,jk,mytid+1) - zwork(js,jk + 1,mytid+1))
                    IF ((Fsize .GT. 0) .AND. (jf .GT. 0)) THEN
                         diaflx(jf,sed_idx(js),4) = diaflx(jf,sed_idx(js),4) + zwork(js,jk,mytid+1)
                    ENDIF
                 END DO

                 DO js =1,21
CCC  d2s convert speed from (m/day) to  (m/s)
                    tra(ji,jj,jk,sed_idx(js)) = tra(ji,jj,jk,sed_idx(js)) + ztra(js,mytid+1)*d2s
                 END DO

#            if defined key_trc_diabio
                  trbio(ji,jj,jk,8) = ztra
#            endif


                END DO
        ENDIF

!$omp    end parallel

      END DO MAIN_LOOP



#    else

C       no Sedimentation

#    endif

      END SUBROUTINE trcsed
