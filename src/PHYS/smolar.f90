      SUBROUTINE smolar
       USE myalloc
       USE myalloc_mpp
       USE ADV_mem
       USE DIA_mem
          implicit none

!!!                      trcadv.smolar.h
!!!                     ******************
!!!
!!   defined key : 'key_trc_smolar'
!!   ============
!!
!!  PURPOSE :
!!  ---------
!!     compute the now trend due to the advection of passive tracers
!!     and add it to the general trend of tracer equations:
!!     THEN computes both horizontal and
!!      vertical advection of tracer trn
!!
!!
!!   METHOD :
!!   -------
!!      this ROUTINE compute not exactly the advection but the
!!      transport term, i.e.  div(u*tra).
!!
!!      smolarkevisz scheme
!!      *******************
!!
!!      computes the now horizontal and vertical advection with the
!!                       ----------     --------
!!      complete 3d method.
!!
!!      cf reference
!!
!!      note: - sc is an empirical factor to be used with care
!!            - this advection scheme needs an euler-forward time scheme
!!
!!   remarks :
!!   -------
!!
!!      multitasked on tracer (jn-loop)

!!
!!   --------
!!
!!   REFERENCES :                piotr k. smolarkiewicz, 1983,
!!   ----------                  "a simple positive definit advection
!!                               scheme with small IMPLICIT diffusion"
!!                               monthly weather review, pp 479-486
!!
      LOGICAL :: MPI_CHECK,l1,l2,l3
      INTEGER ji,jj,jk,jt,jn,jf,ju
      INTEGER jp,pack_size
      REAL(8) zbtr,zdt
      REAL(8) junk, junki, junkj, junkk
      INTEGER A,B


! omp variables
      INTEGER :: mytid, ntids, itid

#ifdef __OPENMP1
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

!-------------------------------------------------------------------

      MPI_CHECK = .FALSE.

#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = mpi_pack_size
      mytid = 0
#endif

      if(allpoints .EQ. 0) then  ! INIT phase

!!!&omp parallel default(shared) private(mytid,A)
#ifdef __OPENMP1
         mytid = omp_get_thread_num()  ! take the thread ID
#else
      PACK_LOOP0: DO jp=1,ntids
       mytid=jp-1
#endif
      A = mytid+1
         zti(:,:,:,A)  = 0.
         ztj(:,:,:,A)  = 0.
         zkx(:,:,:,A)  = 0.
         zky(:,:,:,A)  = 0.
         zkz(:,:,:,A)  = 0.
         zbuf(:,:,:,A) = 0.
         zx(:,:,:,A)   = 0.
         zy(:,:,:,A)   = 0.
         zz(:,:,:,A)   = 0.
!!!&omp end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP0
      mytid =0
#endif

         zaa  = 0.
         zbb  = 0.
         zcc  = 0.
         inv_eu = 0.
         inv_ev = 0.
         inv_et = 0.
         big_fact_zaa = 0.
         big_fact_zbb = 0.
         big_fact_zcc = 0.
         zbtr_arr = 0.
         jarr  = 0
         jarr1 = 0
         jarr2 = 0
         jarr3 = 0
         jarrt = 0
         write(*,*) "Storing good points ..."
         DO jk = 1,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  zbtr_arr(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*e3t(ji,jj,jk))
               END DO
            END DO
         END DO

         DO jk = 1,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  inv_eu(ji,jj,jk) = 1./(e1u(ji,jj)*e2u(ji,jj)*e3u(ji,jj,jk) )
                  inv_ev(ji,jj,jk) = 1./(e1v(ji,jj)*e2v(ji,jj)*e3v(ji,jj,jk) )
               END DO
            END DO
         END DO

         DO jk = 2,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  inv_et(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*e3w(ji,jj,jk) )
               END DO
            END DO
         END DO

         allpoints = 0
         DO jk = 1,jpk
            DO jj = 1,jpj
               DO ji = 1,jpi
                  allpoints = allpoints + 1
                  if(tmask(ji,jj,jk) .NE. 0) then
                     tpoints = tpoints + 1
                  endif
               END DO
            END DO
         END DO

         write(*,*) 'trcadv: RANK -> ', RANK, ' all_points -> ', allpoints

         goodpoints = 0
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
                  DO myjk=jk+jklef, jk+jkrig
                     DO myjj=jj+jjlef, jj+jjrig
                        DO myji=ji+jilef, ji+jirig
                           locsum = locsum + tmask(myji, myjj, myjk)
                        END DO
                     END DO
                  END DO
                  if(locsum .NE. 0) then
                     goodpoints = goodpoints + 1
                     advmask(ji,jj,jk) = 1
                  else
                     advmask(ji,jj,jk) = 0
                  endif
               END DO
            END DO
         END DO

         write(*,*) 'trcadv: RANK -> ', RANK, ' good_points -> ', goodpoints

         DO jk = 1,jpk
            DO jj = 2,jpjm1
               DO  ji = 2,jpim1
                  if(advmask(ji,jj,jk) .NE. 0) then
                     zbtr_arr(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*e3t(ji,jj,jk))
                     dimen_jarr = dimen_jarr + 1
                     jarr(1,dimen_jarr) = ji
                     jarr(2,dimen_jarr) = jj
                     jarr(3,dimen_jarr) = jk
                  else
                     zbtr_arr(ji,jj,jk) = 0.
                  endif
               END DO
            END DO
         END DO
         DO jk = 2,jpk
            DO jj = 2,jpjm1
               DO  ji = 2,jpim1
                  if(advmask(ji,jj,jk) .NE. 0) then
                     dimen_jarr1 = dimen_jarr1 + 1
                     jarr1(1,dimen_jarr1) = ji
                     jarr1(2,dimen_jarr1) = jj
                     jarr1(3,dimen_jarr1) = jk
                  endif
               END DO
            END DO
         END DO
         DO jk = 2,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  if(advmask(ji,jj,jk) .NE. 0) then
                     dimen_jarr2 = dimen_jarr2 + 1
                     jarr2(1,dimen_jarr2) = ji
                     jarr2(2,dimen_jarr2) = jj
                     jarr2(3,dimen_jarr2) = jk
                  endif
               END DO
            END DO
         END DO
         DO jk = 1,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  if(advmask(ji,jj,jk) .NE. 0) then
                     dimen_jarr3 = dimen_jarr3 + 1
                     jarr3(1,dimen_jarr3) = ji
                     jarr3(2,dimen_jarr3) = jj
                     jarr3(3,dimen_jarr3) = jk
                  endif
               END DO
            END DO
         END DO
         DO jk = 1,jpk
            DO jj = 1,jpj
               DO ji = 1,jpi
                  if(tmask(ji,jj,jk) .NE. 0) then
                     dimen_jarrt = dimen_jarrt + 1
                     jarrt(1,dimen_jarrt) = ji
                     jarrt(2,dimen_jarrt) = jj
                     jarrt(3,dimen_jarrt) = jk
                  endif
               END DO
            END DO
         END DO

      jarr_adv_flx=0

         DO jf=1,Fsize
            DO ju=1, dimen_jarr3
               l1 = flx_ridxt(jf,2) .EQ. jarr3(1,ju)
               l2 = flx_ridxt(jf,3) .EQ. jarr3(2,ju)
               l3 = flx_ridxt(jf,4) .EQ. jarr3(3,ju)
               IF ( l1 .AND. l2 .AND. l3) THEN
                  jarr_adv_flx(ju)= jf
               END IF
            END DO
         END DO

      endif ! end initialization phase

         jk=1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  zbtr_arr(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*e3t(ji,jj,jk))
               END DO
            END DO

            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  inv_eu(ji,jj,jk) = 1./(e1u(ji,jj)*e2u(ji,jj)*e3u(ji,jj,jk) )
                  inv_ev(ji,jj,jk) = 1./(e1v(ji,jj)*e2v(ji,jj)*e3v(ji,jj,jk) )
               END DO
            END DO

            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  inv_et(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*e3w(ji,jj,jk) )
               END DO
            END DO





      zdt = rdt*ndttrc


!!!&omp  parallel default(none) private(A,mytid,jk,jj,ji) shared(jpk,jpj,jpi,ztj,zx,zy,zz)
#ifdef __OPENMP1
       mytid = omp_get_thread_num()  ! take the thread ID
#else
      PACK_LOOP1: DO jp=1,ntids
       mytid=jp-1
#endif
       A = mytid + 1

       DO jk = 1,jpk
          DO jj = 1,jpj
             DO ji = 1,jpi
                ztj(ji,jj,jk,A) = 0.
                 zx(ji,jj,jk,A) = 0.
                 zy(ji,jj,jk,A) = 0.
                 zz(ji,jj,jk,A) = 0.
             END DO
          END DO
       END DO
!!!&omp  end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP1
      mytid =0
#endif

       DO jk = 1,jpk,ntids
!!!&omp parallel default(none) private(mytid, ji,jj)
!!!&omp&         shared(jk, jpk,jpj,jpi,big_fact_zaa,big_fact_zbb,big_fact_zcc,zaa,zbb,zcc,inv_eu,inv_ev,inv_et,
!!!&omp&                un,vn,wn,e2u,e3u,e3v,e1v,e1t,e2t,e3t,zdt )
#ifdef __OPENMP1
       mytid = omp_get_thread_num()  ! take the thread ID
#else
      PACK_LOOP2: DO jp=1,ntids
       mytid=jp-1
#endif

       if (jk+mytid.le.jpk) then
          DO jj = 1,jpj
             DO ji = 1,jpi
                zaa(ji,jj,jk+mytid) = e2u(ji,jj)*e3u(ji,jj,jk+mytid) * un(ji,jj,jk+mytid)
                zbb(ji,jj,jk+mytid) = e1v(ji,jj)*e3v(ji,jj,jk+mytid) * vn(ji,jj,jk+mytid)
                zcc(ji,jj,jk+mytid) = e1t(ji,jj)*e2t(ji,jj)      * wn(ji,jj,jk+mytid)
                big_fact_zaa(ji,jj,jk+mytid) = ( abs(zaa(ji,jj,jk+mytid)) - zdt*zaa(ji,jj,jk+mytid)**2*inv_eu(ji,jj,jk+mytid) )!/(e1u(ji,jj)*e2u(ji,jj)*e3t(ji,jj,jk+mytid) ) )
                big_fact_zbb(ji,jj,jk+mytid) = ( abs(zbb(ji,jj,jk+mytid)) - zdt*zbb(ji,jj,jk+mytid)**2*inv_ev(ji,jj,jk+mytid) )!/(e1v(ji,jj)*e2v(ji,jj)*e3t(ji,jj,jk+mytid) ) )
                big_fact_zcc(ji,jj,jk+mytid) = ( abs(zcc(ji,jj,jk+mytid)) - zdt*zcc(ji,jj,jk+mytid)**2*inv_et(ji,jj,jk+mytid) )!/(e1t(ji,jj)*e2t(ji,jj)*e3w(ji,jj,jk+mytid) ) )
             END DO
          END DO

         endif
!!!&omp end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP2
      mytid =0
#endif
       END DO


!!     tracer loop parallelized (macrotasking)
!!     =======================================

      trcadvparttime = MPI_WTIME()
      TRACER_LOOP: DO  jn = 1, jptra, ntids


!!        1. tracer flux in the 3 directions
!!        ----------------------------------
!!        1.1 mass flux at u v and t-points and initialization
!!       1.2 calcul of intermediate field with an upstream advection scheme
!!           and mass fluxes calculated above
!!       calcul of tracer flux in the i and j direction

!!!&omp   parallel default(none) private(jk,jj,ji,mytid,A,B,junk)
!!!&omp&      shared(jpk,jpj,jpi,zkx,zky,zbb,zaa,zcc,jarr1,dimen_jarr1,jn,zkz,trn,jpjm1,jpim1)

#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#else
        PACK_LOOP4: DO jp=1,ntids
         mytid=jp-1
#endif

      if( mytid + jn <= jptra ) then
      A = mytid + 1
      B = mytid + jn
           DO jk = 1,jpk
             DO jj=1,jpj; zkx(  1,jj,jk,A)=0.  ; END DO
             DO jj=1,jpj; zkx(jpi,jj,jk,A)=0.  ; END DO

             DO ji=1,jpi; zky(ji,  1,jk,A)=0.  ; END DO
             DO ji=1,jpi; zky(ji,jpj,jk,A)=0.  ; END DO
           END DO

!!
!!       calcul of tracer flux in the k direction
!!
           DO jj = 1,jpj
             DO ji = 1,jpi
               zkz(ji,jj,1,A) = 0.
!               zkz(ji,jj,1,A) = fsz(trn(ji,jj,1,B),trn(ji,jj,1,B),zcc(ji,jj,1))
             END DO
           END DO

           jk = 1
           DO jj = 2,jpjm1
              DO ji = 2,jpim1
                 zkx(ji,jj,jk,A) = fsx(trn(ji,jj,jk,B),trn(ji + 1,jj,jk,B),zaa(ji,jj,jk))
                 zky(ji,jj,jk,A) = fsy(trn(ji,jj,jk,B),trn(ji,jj + 1,jk,B),zbb(ji,jj,jk))
              END DO
           END DO

           DO jk = 2,jpk
             jj=   1 ; DO ji = 1,jpi; zkz(ji,jj,jk,A) = fsz(trn(ji,jj,jk,B),trn(ji,jj,jk - 1,B),zcc(ji,jj,jk)); END DO
             jj= jpj ; DO ji = 1,jpi; zkz(ji,jj,jk,A) = fsz(trn(ji,jj,jk,B),trn(ji,jj,jk - 1,B),zcc(ji,jj,jk)); END DO

             DO jj = 2,jpjm1
                ji=1   ; zkz(ji,jj,jk,A) = fsz(trn(ji,jj,jk,B),trn(ji,jj,jk - 1,B),zcc(ji,jj,jk))
                ji=jpi ; zkz(ji,jj,jk,A) = fsz(trn(ji,jj,jk,B),trn(ji,jj,jk - 1,B),zcc(ji,jj,jk))
             END DO
           END DO

           DO ju=1, dimen_jarr1

              ji = jarr1(1, ju)
              jj = jarr1(2, ju)
              jk = jarr1(3, ju)
              junk = trn(ji,jj,jk,B)
              zkx(ji,jj,jk,A) = fsx(junk,trn(ji + 1,jj,jk,B),zaa(ji,jj,jk))
              zky(ji,jj,jk,A) = fsy(junk,trn(ji,jj + 1,jk,B),zbb(ji,jj,jk))
              zkz(ji,jj,jk,A) = fsz(junk,trn(ji,jj,jk - 1,B),zcc(ji,jj,jk))

           END DO

      end if

!!!&omp end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP4
      mytid=0
#endif

!! ... Lateral boundary conditions on zk[xy]
#ifdef key_mpp

!!   ... Mpp : export boundary values to neighboring processors


      IF( ntids - 1 + jn <= jptra ) THEN
       pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
      END IF
        CALL mpplnk_my(zkx, pack_size,1,1)
        CALL mpplnk_my(zky, pack_size,1,1)




#else

!!   ... T-point, 3D array, full local arrays zk[xy] are initialised


        DO itid = 1, ntids

       IF( itid - 1 + jn <= jptra ) THEN

              CALL lbc( zkx(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
              CALL lbc( zky(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )

       END IF

        END DO

#endif


!! 2. calcul of after field using an upstream advection scheme
!! -----------------------------------------------------------


!!!&omp   parallel default(none) private(mytid,A,B,zbtr,ji,jj,jk,ju,jf)
!!!&omp&      shared(zkx,zky,zkz,zti,jpim1,jpjm1,trn,zdt,jn,jpkm1,zbtr_arr,e1t,e2t,ztj,jarr3,ncor,dimen_jarr3,
!!!&omp&             jarr_adv_flx,Fsize,diaflx)

#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#else
        PACK_LOOP5: DO jp=1,ntids
         mytid=jp-1
#endif

      IF( mytid + jn <= jptra ) THEN
      A = mytid +1
      B = mytid + jn
           DO ju=1, dimen_jarr3

              ji = jarr3(1, ju)
              jj = jarr3(2, ju)
              jk = jarr3(3, ju)
              jf = jarr_adv_flx(ju)

              zbtr = zbtr_arr(ji,jj,jk)

              ztj(ji,jj,jk,A) = -zbtr* &
     &          ( zkx(ji,jj,jk,A) - zkx(ji - 1,jj,jk,A) &
     &          + zky(ji,jj,jk,A) - zky(ji,jj - 1,jk,A) &
     &          + zkz(ji,jj,jk,A) - zkz(ji,jj,jk + 1,A) )


              IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                    diaflx(jf,B,1) = diaflx(jf,B,1) + zkx(ji,jj,jk,A)
                    diaflx(jf,B,2) = diaflx(jf,B,2) + zky(ji,jj,jk,A)
                    diaflx(jf,B,3) = diaflx(jf,B,3) + zkz(ji,jj,jk,A)
              END IF

            END DO
      END IF

!!!&omp end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP5
      mytid=0
#endif

!! 2.1 start of antidiffusive correction loop

        ANTIDIFF_CORR: DO jt = 1,ncor

!! 2.2 calcul of intermediary field zti


!!!&omp     parallel default(none) private(mytid,A,B,ji,jj,jk)
!!!&omp&       shared(jt,jn,ncor,jpkm1,jpjm1,jpim1,zti,ztj,trn,zdt,zbuf)

#ifdef __OPENMP1
          mytid = omp_get_thread_num()  ! take the thread ID
#else
        PACK_LOOP6: DO jp=1,ntids
         mytid=jp-1
#endif

      IF( mytid + jn <= jptra ) THEN
      A = mytid +1
      B = mytid + jn
             if(jt .EQ. 1) then

                if(ncor .EQ. 1) then
                   DO jk = 1,jpkm1
                      DO jj = 2,jpjm1
                         DO ji = 2,jpim1
                            zti(ji,jj,jk,A) = trn(ji,jj,jk,B) + zdt*ztj(ji,jj,jk,A)!zbuf(ji,jj,jk) = 0. + ztj(ji,jj,jk,mytid+1)
                         END DO
                      END DO
                   END DO
                else
                    DO jk = 1,jpkm1
                      DO jj = 2,jpjm1
                         DO ji = 2,jpim1
                            zti(ji,jj,jk,A) = trn(ji,jj,jk,B) + zdt*ztj(ji,jj,jk,A)
                            zbuf(ji,jj,jk,A) = ztj(ji,jj,jk,A)
                         END DO
                      END DO
                   END DO
                endif

             else

                DO jk = 1,jpkm1
                   DO jj = 2,jpjm1
                      DO ji = 2,jpim1
                         zti(ji,jj,jk,A) =  zti(ji,jj,jk,A) + zdt*ztj(ji,jj,jk,A)
                         zbuf(ji,jj,jk,A)  = zbuf(ji,jj,jk,A)       + ztj(ji,jj,jk,A)
                      END DO
                   END DO
                END DO

             endif

      END IF

!!!&omp     end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP6
      mytid=0
#endif



!! ... Lateral boundary conditions on zti

#ifdef key_mpp

!!   ... Mpp : export boundary values to neighboring processors



        IF( ntids - 1 + jn <= jptra ) THEN
           pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
        END IF

        CALL mpplnk_my(zti, pack_size,1,1)



#else

!!   ... T-point, 3D array, full local array zti is initialised

          DO itid = 1, ntids
         IF( itid - 1 + jn <= jptra ) THEN
                CALL lbc( zti(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
         END IF
          END DO

#endif


!! 2.3 calcul of the antidiffusive flux

!!!&omp     parallel default(none) private(mytid,A,junk, junki, junkj, junkk, ji,jj,jk)
!!!&omp&       shared(jn,jpkm1,jpjm1,jpim1,zti,ztj,zy,zx,zz,jarr2,big_fact_zbb,
!!!&omp&              big_fact_zaa,big_fact_zcc,dimen_jarr2,rtrn,rsc)

#ifdef __OPENMP1
          mytid = omp_get_thread_num()  ! take the thread ID
#else
        PACK_LOOP7: DO jp=1,ntids
         mytid=jp-1
#endif

      IF( mytid + jn <= jptra ) THEN
      A = mytid +1
            jk = 1
!          DO jk = 1,jpkm1
            DO jj = 2,jpjm1
              DO ji = 2,jpim1
                  junk  = zti(ji,jj,jk,A)
                  junki = zti(ji + 1,jj,jk,A)
                  junkj = zti(ji,jj + 1,jk,A)
!                 if(advmask(ji,jj,jk) .NE. 0) then
!                zx(ji,jj,jk,A) = ( abs(zaa(ji,jj,jk)) - zdt*zaa(ji,jj,jk)**2/(e1u(ji,jj)*e2u(ji,jj)*e3t(jk) ) )*inv_eu(ji,jj,jk) )
!                zy(ji,jj,jk,A) = ( abs(zbb(ji,jj,jk)) - zdt*zbb(ji,jj,jk)**2/(e1v(ji,jj)*e2v(ji,jj)*e3t(jk) ) )*inv_ev(ji,jj,jk) )
                    zx(ji,jj,jk,A) = big_fact_zaa(ji,jj,jk)*(junki - junk)/(junk + junki + rtrn)* rsc
                    zy(ji,jj,jk,A) = big_fact_zbb(ji,jj,jk)*(junkj - junk)/(junk + junkj + rtrn)* rsc


!                endif

              END DO
            END DO
!          END DO
!!!!!!!!!!!!!!!!! To be checked
!!          DO jj = 2,jpjm1
!!            DO ji = 2,jpim1
!!                junk  = zti(ji,jj,2,A)
!!                junkk = zti(ji,jj,1,A)
!!                zz(ji,jj,1,A) = big_fact_zaa(ji,jj,1)*(junk - junkk)/(junk + junkk + rtrn)* rsc*(-1.)

!!            END DO
!!          END DO
!!!!!!!!!!!!!!!!! 

!                 if(advmask(ji,jj,jk) .NE. 0) then

            DO ju=1, dimen_jarr2

               ji = jarr2(1, ju)
               jj = jarr2(2, ju)
               jk = jarr2(3, ju)
               junk  = zti(ji,jj,jk,A)
               junki = zti(ji + 1,jj,jk,A)
               junkj = zti(ji,jj + 1,jk,A)
               junkk = zti(ji,jj,jk - 1,A)
               zx(ji,jj,jk,A) = big_fact_zaa(ji,jj,jk)*(junki - junk)/(junk + junki + rtrn)* rsc
               zy(ji,jj,jk,A) = big_fact_zbb(ji,jj,jk)*(junkj - junk)/(junk + junkj + rtrn)* rsc
               zz(ji,jj,jk,A) = big_fact_zcc(ji,jj,jk)*(junk - junkk)/(junk + junkk + rtrn)* rsc*( -1.)

           END DO
!                 endif

          END IF

!!!&omp end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP7
      mytid=0
#endif




!! ... Lateral boundary conditions on z[xyz]
#ifdef key_mpp

!!   ... Mpp : export boundary values to neighboring processors


      IF( ntids - 1 + jn <= jptra ) THEN
       pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
      END IF

        CALL mpplnk_my(zx, pack_size,1,1)
        CALL mpplnk_my(zy, pack_size,1,1)
        CALL mpplnk_my(zz, pack_size,1,1)


#else

!!   ... T-point, 3D array, full local array z[xyz] are initialised


          DO itid = 1, ntids
         IF( itid - 1 + jn <= jptra ) THEN
                CALL lbc( zx(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zy(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zz(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
         END IF
      END DO

#endif

!! 2.4 reinitialization

!!!&omp     parallel default(none) private(mytid,A,junk,jk,jj,ji)
!!!&omp&      shared(zkx,zky,zkz,zz,zx,zy,zti,jpjm1,jpim1,dimen_jarr1,jarr1,jpi,jpk,jpj,jn)

#ifdef __OPENMP1
          mytid = omp_get_thread_num()  ! take the thread ID
#else
        PACK_LOOP8: DO jp=1,ntids
         mytid=jp-1
#endif


      IF( mytid + jn <= jptra ) THEN
      A=mytid+1
!!            2.5 calcul of the final field:
!!                advection by antidiffusive mass fluxes and an upstream scheme
             jk = 1
             DO jj = 2,jpjm1
                 DO ji = 2,jpim1
                   junk  = zti(ji,jj,jk,A)
                   zkx(ji,jj,jk,A) = fsx(junk,zti(ji + 1,jj,jk,A),zx(ji,jj,jk,A))!zti(ji,jj,jk,A),zti(ji + 1,jj,jk,A),zaa(ji,jj,jk))
                   zky(ji,jj,jk,A) = fsy(junk,zti(ji,jj + 1,jk,A),zy(ji,jj,jk,A))!zti(ji,jj,jk,A),zti(ji,jj + 1,jk,A),zbb(ji,jj,jk))
!!!!!! To be checked   zkz(ji,jj,jk,A) = fsz(junk,junk,zz(ji,jj,jk,A))
                 END DO
              END DO

             DO jk = 2,jpk
               jj =  1; DO ji = 1,jpi ; zkz(ji,jj,jk,A) = fsz(zti(ji,jj,jk,A),zti(ji,jj,jk - 1,A),zz(ji,jj,jk,A)); END DO!zti(ji,jj,jk,A),zti(ji,jj,jk - 1,A),zcc(ji,jj,jk))
               jj = jpj;DO ji = 1,jpi ; zkz(ji,jj,jk,A) = fsz(zti(ji,jj,jk,A),zti(ji,jj,jk - 1,A),zz(ji,jj,jk,A)); END DO!zti(ji,jj,jk,A),zti(ji,jj,jk - 1,A),zcc(ji,jj,jk))


               DO jj = 2,jpjm1
                 ji =   1 ; zkz(ji,jj,jk,A) = fsz(zti(ji,jj,jk,A),zti(ji,jj,jk - 1,A),zz(ji,jj,jk,A))!zti(ji,jj,jk,A),zti(ji,jj,jk - 1,A),zcc(ji,jj,jk))
                 ji = jpi ; zkz(ji,jj,jk,A) = fsz(zti(ji,jj,jk,A),zti(ji,jj,jk - 1,A),zz(ji,jj,jk,A))!zti(ji,jj,jk,A),zti(ji,jj,jk - 1,A),zcc(ji,jj,jk))
               END DO
             END DO

             DO ju=1, dimen_jarr1

                ji = jarr1(1, ju)
                jj = jarr1(2, ju)
                jk = jarr1(3, ju)
                junk  = zti(ji,jj,jk,A)
                zkx(ji,jj,jk,A) = fsx(junk,zti(ji + 1,jj,jk,A),zx(ji,jj,jk,A))!zti(ji,jj,jk,A),zti(ji + 1,jj,jk,A),zaa(ji,jj,jk))
                zky(ji,jj,jk,A) = fsy(junk,zti(ji,jj + 1,jk,A),zy(ji,jj,jk,A))!zti(ji,jj,jk,A),zti(ji,jj + 1,jk,A),zbb(ji,jj,jk))
                zkz(ji,jj,jk,A) = fsz(junk,zti(ji,jj,jk - 1,A),zz(ji,jj,jk,A))!zti(ji,jj,jk,A),zti(ji,jj,jk - 1,A),zcc(ji,jj,jk))

             END DO

      END IF

!!!&omp end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP8
      mytid=0
#endif


!! ... Lateral boundary conditions on zk[xy]
#ifdef key_mpp
!!   ... Mpp : export boundary values to neighboring processors


        IF( ntids - 1 + jn <= jptra ) THEN
           pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
        END IF

        CALL mpplnk_my(zkx, pack_size,1,1)
        CALL mpplnk_my(zky, pack_size,1,1)
#else
!!   ... T-point, 3D array, full local array zk[xy] are initialised
      DO itid = 1, ntids
        IF( itid - 1 + jn <= jptra ) THEN
               CALL lbc( zkx(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
               CALL lbc( zky(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
        END IF
      END DO
#endif

!!!&omp    parallel default(none) private(mytid,A,B,zbtr,ji,jj,jk,ju,jf)
!!!&omp&      shared(zkx,zky,zkz,zbtr_arr,e1t,e2t,ztj,dimen_jarr3,jarr3,ncor,jn,
!!!&omp&             jarr_adv_flx,Fsize,diaflx)

#ifdef __OPENMP1
         mytid = omp_get_thread_num()  ! take the thread ID
#else
        PACK_LOOP9: DO jp=1,ntids
         mytid=jp-1
#endif
!!        2.6. calcul of after field using an upstream advection scheme


         IF( mytid + jn <= jptra ) THEN
         A = mytid+1
         B = mytid+jn
            if(ncor .EQ. 1) then
               DO ju=1, dimen_jarr3

                  ji = jarr3(1, ju)
                  jj = jarr3(2, ju)
                  jk = jarr3(3, ju)
                  jf = jarr_adv_flx(ju)

!                 if(zbtr_arr(ji,jj,jk) .NE. 0) then zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*e3t(jk))
                  zbtr = zbtr_arr(ji,jj,jk)
                  ztj(ji,jj,jk,A) = -zbtr*( zkx(ji,jj,jk,A) - zkx(ji - 1,jj,jk,A) &
     &              + zky(ji,jj,jk,A) - zky(ji,jj - 1,jk,A)+ zkz(ji,jj,jk,A) - zkz(ji,jj,jk + 1,A) )+ ztj(ji,jj,jk,A)

!     Save advective fluxes x,y,z
              IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                 diaflx(jf,B,1) = diaflx(jf,B,1) + zkx(ji,jj,jk,A)
                 diaflx(jf,B,2) = diaflx(jf,B,2) + zky(ji,jj,jk,A)
                 diaflx(jf,B,3) = diaflx(jf,B,3) + zkz(ji,jj,jk,A)
              END IF

              END DO
           else
              DO ju=1, dimen_jarr3

                  ji = jarr3(1, ju)
                  jj = jarr3(2, ju)
                  jk = jarr3(3, ju)
                  jf = jarr_adv_flx(ju)

!                if(zbtr_arr(ji,jj,jk) .NE. 0) then    zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*e3t(jk))
                  zbtr = zbtr_arr(ji,jj,jk)
                  ztj(ji,jj,jk,A) = -zbtr*( zkx(ji,jj,jk,A) - zkx(ji - 1,jj,jk,A) &
     &              + zky(ji,jj,jk,A) - zky(ji,jj - 1,jk,A)+ zkz(ji,jj,jk,A) - zkz(ji,jj,jk + 1,A) )

!     Save advective fluxes x,y,z
                 IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                    diaflx(jf,B,1) = diaflx(jf,B,1) + zkx(ji,jj,jk,A)
                    diaflx(jf,B,2) = diaflx(jf,B,2) + zky(ji,jj,jk,A)
                    diaflx(jf,B,3) = diaflx(jf,B,3) + zkz(ji,jj,jk,A)
                 END IF
              END DO

           endif

      END IF

!!!&omp end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP9
      mytid=0
#endif
        END DO ANTIDIFF_CORR


!!       3. trend due to horizontal and vertical advection of tracer jn
!!!&omp   parallel default(none) private(mytid,A,B,ji,jj,jk,ju) shared(ncor,dimen_jarrt,jarrt,tra,ztj,jn,zbuf)

#ifdef __OPENMP1
        mytid = omp_get_thread_num()  ! take the thread ID
#else
        PACK_LOOP10: DO jp=1,ntids
         mytid=jp-1
#endif
        A=mytid+1
        B=mytid+jn

        IF( mytid + jn <= jptra ) THEN


           if(ncor .EQ. 1) then
              DO ju=1, dimen_jarrt
                 ji = jarrt(1, ju)
                 jj = jarrt(2, ju)
                 jk = jarrt(3, ju)

                 tra(ji,jj,jk,B) = tra(ji,jj,jk,B)+ ztj(ji,jj,jk,A)!+ (zbuf(ji,jj,jk,A) + ztj(ji,jj,jk,A))

              END DO
           else
              DO ju=1, dimen_jarrt
                 ji = jarrt(1, ju)
                 jj = jarrt(2, ju)
                 jk = jarrt(3, ju)

                 tra(ji,jj,jk,B) = tra(ji,jj,jk,B)+ (zbuf(ji,jj,jk,A) + ztj(ji,jj,jk,A))

              END DO
           endif

      END IF
!!!&omp end parallel
#ifdef __OPENMP1
#else
      END DO PACK_LOOP10
      mytid=0
#endif
       END DO TRACER_LOOP

       trcadvparttime = MPI_WTIME() - trcadvparttime
       trcadvtottime = trcadvtottime + trcadvparttime

      END SUBROUTINE smolar
