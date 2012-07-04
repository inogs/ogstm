      SUBROUTINE smolar
       USE myalloc
       USE myalloc_mpp
       USE ADV_mem
       USE stafun
       USE DIA_mem
          implicit none
CCC
CCC                      trcadv.smolar.h
CCC                     ******************
CCC
CC   defined key : 'key_trc_smolar'
CC   ============
CC
CC  PURPOSE :
CC  ---------
CC     compute the now trend due to the advection of passive tracers
CC     and add it to the general trend of tracer equations:
CC     THEN computes both horizontal and
CC      vertical advection of tracer trn
CC
CC
CC   METHOD :
CC   -------
CC      this ROUTINE compute not exactly the advection but the
CC      transport term, i.e.  div(u*tra).
CC
CC      smolarkevisz scheme
CC      *******************
CC
CC      computes the now horizontal and vertical advection with the
CC                       ----------     --------
CC      complete 3d method.
CC
CC      cf reference
CC
CC      note: - sc is an empirical factor to be used with care
CC            - this advection scheme needs an euler-forward time scheme
CC
CC   remarks :
CC   -------
CC      IF 'key_trc_diatrd' key is activated, the
CC      advection trend is saved for futher diagnostics.
CC
CC      multitasked on tracer (jn-loop)
CC
CC      COMMON
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC            /cottrp/          : present and next fields for passive
CC                              : tracers
CC            /comtsk/          : multitasking
CC            /cot3ad/          : smolar advection
CC
CC   OUTPUT :
CC   ------
CC      COMMON
CC            /cottrp/ tra      : general tracer trend increased by the
CC                                now horizontal tracer advection trend
CC            /cottrd/ trtrd    : now horizontal tracer advection trend
CC                                (IF 'key_trc_diatrd' key is activated)
CC
CC
CC   EXTERNAL :                  upstream, forkjoin, mpplnk, lbc
CC   --------
CC
CC   REFERENCES :                piotr k. smolarkiewicz, 1983,
CC   ----------                  "a simple positive definit advection
CC                               scheme with small IMPLICIT diffusion"
CC                               monthly weather review, pp 479-486
CC
CC   MODIFICATIONS:
CC   --------------
CC      original :  87-06 (pa-dl)
CC      additions : 91-11 (G. Madec)
CC      additions : 94-08 (a. czaja)
CC      modifications : 95-09 (M. Levy) passive tracers
CC      modifications : 98-03 (M.A. Foujols) lateral boundary conditions
CC      modifications : 99-02 (M.A. Foujols) lbc in conjonction with ORCA
CC      modifications : 00-05 (MA Foujols) add lbc for tracer trends
CC      modifications : 00-10 (MA Foujols and E.Kestenare) INCLUDE instead of routine
CC      modifications : 01-05 (E.Kestenare) fix bug in trtrd indexes
CC----------------------------------------------------------------------
      LOGICAL :: MPI_CHECK,l1,l2,l3
      INTEGER ji,jj,jk,jt,jn,jf,ju
      INTEGER pack_size
      REAL(8) zbtr,zdt
      REAL(8) junk, junki, junkj, junkk
      INTEGER A,B
!      REAL(8) zgm,zgz


! omp variables
      INTEGER :: mytid, ntids, itid

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

!-------------------------------------------------------------------

      MPI_CHECK = .FALSE.

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = 1
      mytid = 0
#endif

      if(allpoints .EQ. 0) then  ! INIT phase

!$omp parallel default(shared) private(mytid,A)
#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
      A = mytid+1
         zti(:,:,:,A)  = 0.
         ztj(:,:,:,A)  = 0.
         zkx(:,:,:,A)  = 0.
         zky(:,:,:,A)  = 0.
         zkz(:,:,:,A)  = 0.
         zbuf = 0.
         zx(:,:,:,A)   = 0.
         zy(:,:,:,A)   = 0.
         zz(:,:,:,A)   = 0.
!$omp end parallel

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
                  zbtr_arr(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
               END DO
            END DO
         END DO

         DO jk = 1,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  inv_eu(ji,jj,jk) = 1./(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) )
                  inv_ev(ji,jj,jk) = 1./(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) )
C                  inv_etot(1,ji,jj,jk) = inv_eu(ji,jj,jk)
C                  inv_etot(2,ji,jj,jk) = inv_ev(ji,jj,jk)
               END DO
            END DO
         END DO

         DO jk = 2,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  inv_et(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk) )
C                  inv_etot(3,ji,jj,jk) = inv_et(ji,jj,jk)
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
                     zbtr_arr(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
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



      zdt = rdt*float(ndttrc)


!$omp  parallel default(none) private(A,mytid,jk,jj,ji) shared(jpk,jpj,jpi,ztj,zx,zy,zz)
#ifdef __OPENMP
       mytid = omp_get_thread_num()  ! take the thread ID
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
!$omp  end parallel

       DO jk = 1,jpk
          DO jj = 1,jpj
             DO ji = 1,jpi
                zaa(ji,jj,jk) = e2u(ji,jj)*fse3u(ji,jj,jk) * un(ji,jj,jk)
                zbb(ji,jj,jk) = e1v(ji,jj)*fse3v(ji,jj,jk) * vn(ji,jj,jk)
                zcc(ji,jj,jk) = e1t(ji,jj)*e2t(ji,jj)      * wn(ji,jj,jk)
                big_fact_zaa(ji,jj,jk) = ( abs(zaa(ji,jj,jk)) - zdt*zaa(ji,jj,jk)**2*inv_eu(ji,jj,jk) )!/(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) ) )
                big_fact_zbb(ji,jj,jk) = ( abs(zbb(ji,jj,jk)) - zdt*zbb(ji,jj,jk)**2*inv_ev(ji,jj,jk) )!/(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) ) )
                big_fact_zcc(ji,jj,jk) = ( abs(zcc(ji,jj,jk)) - zdt*zcc(ji,jj,jk)**2*inv_et(ji,jj,jk) )!/(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk) ) )
             END DO
          END DO
       END DO

C

C
C     tracer loop parallelized (macrotasking)
C     =======================================
C
      trcadvparttime = MPI_WTIME()
      TRACER_LOOP: DO  jn = 1, jptra, ntids

C
C        1. tracer flux in the 3 directions
C        ----------------------------------
C
C        1.1 mass flux at u v and t-points and initialization
C

#if defined key_trc_diatrd

!$omp    parallel default(none) private(jk,jj,ji,mytid,B) shared(jpk,jpj,jpi,jn,trtrd)

#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif

      if( mytid + jn <= jptra ) then
      B = mytid+jn
           DO jk = 1,jpk
             DO jj = 1,jpj
               DO ji = 1,jpi
                 trtrd(ji,jj,jk,B,1) = 0.
                 trtrd(ji,jj,jk,B,2) = 0.
                 trtrd(ji,jj,jk,B,3) = 0.
               END DO
             END DO
           END DO

      end if

!$omp    end parallel

#endif


C       1.2 calcul of intermediate field with an upstream advection scheme
C           and mass fluxes calculated above

C       calcul of tracer flux in the i and j direction

!$omp   parallel default(none) private(jk,jj,ji,mytid,A,B,junk)
!$omp&      shared(jpk,jpj,jpi,zkx,zky,zbb,zaa,zcc,jarr1,dimen_jarr1,jn,zkz,trn,jpjm1,jpim1)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
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

C
C       calcul of tracer flux in the k direction
C
           DO jj = 1,jpj
             DO ji = 1,jpi
               zkz(ji,jj,1,A) = 0.
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

!$omp end parallel


C ... Lateral boundary conditions on zk[xy]
#ifdef key_mpp

C   ... Mpp : export boundary values to neighboring processors



C       DO itid = 1, ntids

C  IF( itid - 1 + jn <= jptra ) THEN


C           CALL mpplnk( zkx(:,:,:,itid), 1, 1 )
C           CALL mpplnk( zky(:,:,:,itid), 1, 1 )

C  END IF

C       END DO

      IF( ntids - 1 + jn <= jptra ) THEN
       pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
      END IF

        CALL mpplnk_my(zkx, pack_size,1,1)
        CALL mpplnk_my(zky, pack_size,1,1)




#else

C   ... T-point, 3D array, full local arrays zk[xy] are initialised


        DO itid = 1, ntids

       IF( itid - 1 + jn <= jptra ) THEN

              CALL lbc( zkx(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
              CALL lbc( zky(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )

       END IF

        END DO

#endif


C 2. calcul of after field using an upstream advection scheme
C -----------------------------------------------------------


!$omp   parallel default(none) private(mytid,A,B,zbtr,ji,jj,jk,ju,jf)
!$omp&      shared(zkx,zky,zkz,zti,jpim1,jpjm1,trn,zdt,jn,jpkm1,zbtr_arr,e1t,e2t,ztj,jarr3,ncor,dimen_jarr3,
!$omp&             jarr_adv_flx,Fsize,diaflx)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
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

              ztj(ji,jj,jk,A) = -zbtr*
     $          ( zkx(ji,jj,jk,A) - zkx(ji - 1,jj,jk,A)
     $          + zky(ji,jj,jk,A) - zky(ji,jj - 1,jk,A)
     $          + zkz(ji,jj,jk,A) - zkz(ji,jj,jk + 1,A) )
#if defined key_trc_diatrd
              trtrd(ji,jj,jk,B,1) = trtrd(ji,jj,jk,B,1) -zbtr*( zkx(ji,jj,jk,A) - zkx(ji - 1,jj,jk,A) )
              trtrd(ji,jj,jk,B,2) = trtrd(ji,jj,jk,B,2) -zbtr*( zky(ji,jj,jk,A) - zky(ji,jj - 1,jk,A) )
              trtrd(ji,jj,jk,B,3) = trtrd(ji,jj,jk,B,3) -zbtr*( zkz(ji,jj,jk,A) - zkz(ji,jj,jk + 1,A) )
#endif

              IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                    diaflx(jf,B,1) = diaflx(jf,B,1) + zkx(ji,jj,jk,A)
                    diaflx(jf,B,2) = diaflx(jf,B,2) + zky(ji,jj,jk,A)
                    diaflx(jf,B,3) = diaflx(jf,B,3) + zkz(ji,jj,jk,A)
              END IF

            END DO
      END IF

!$omp end parallel


C 2.1 start of antidiffusive correction loop

        ANTIDIFF_CORR: DO jt = 1,ncor

C 2.2 calcul of intermediary field zti


!$omp     parallel default(none) private(mytid,A,B,ji,jj,jk)
!$omp&       shared(jt,jn,ncor,jpkm1,jpjm1,jpim1,zti,ztj,trn,zdt)

#ifdef __OPENMP
          mytid = omp_get_thread_num()  ! take the thread ID
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
                            zbuf(ji,jj,jk) = ztj(ji,jj,jk,A)
                         END DO
                      END DO
                   END DO
                endif

             else

                DO jk = 1,jpkm1
                   DO jj = 2,jpjm1
                      DO ji = 2,jpim1
                         zti(ji,jj,jk,A) =  zti(ji,jj,jk,A) + zdt*ztj(ji,jj,jk,A)
                         zbuf(ji,jj,jk)  = zbuf(ji,jj,jk)       + ztj(ji,jj,jk,A)
                      END DO
                   END DO
                END DO

             endif

      END IF

!$omp     end parallel




C ... Lateral boundary conditions on zti

#ifdef key_mpp

C   ... Mpp : export boundary values to neighboring processors




C        DO itid = 1, ntids
C    IF( itid - 1 + jn <= jptra ) THEN
C              CALL mpplnk( zti(:,:,:,itid), 1, 1 )
C     END IF
C         END DO

        IF( ntids - 1 + jn <= jptra ) THEN
           pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
        END IF

        CALL mpplnk_my(zti, pack_size,1,1)



#else

C   ... T-point, 3D array, full local array zti is initialised

          DO itid = 1, ntids
         IF( itid - 1 + jn <= jptra ) THEN
                CALL lbc( zti(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
         END IF
          END DO

#endif


C 2.3 calcul of the antidiffusive flux

!$omp     parallel default(none) private(mytid,A,junk, junki, junkj, junkk, ji,jj,jk)
!$omp&       shared(jn,jpkm1,jpjm1,jpim1,zti,ztj,zy,zx,zz,jarr2,big_fact_zbb,
!$omp&              big_fact_zaa,big_fact_zcc,dimen_jarr2,rtrn,rsc)

#ifdef __OPENMP
          mytid = omp_get_thread_num()  ! take the thread ID
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
!                zx(ji,jj,jk,A) = ( abs(zaa(ji,jj,jk)) - zdt*zaa(ji,jj,jk)**2/(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) ) )*inv_eu(ji,jj,jk) )
!                zy(ji,jj,jk,A) = ( abs(zbb(ji,jj,jk)) - zdt*zbb(ji,jj,jk)**2/(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) ) )*inv_ev(ji,jj,jk) )
                    zx(ji,jj,jk,A) = big_fact_zaa(ji,jj,jk)*(junki - junk)/(junk + junki + rtrn)* rsc
                    zy(ji,jj,jk,A) = big_fact_zbb(ji,jj,jk)*(junkj - junk)/(junk + junkj + rtrn)* rsc


!                endif

              END DO
            END DO
!          END DO

!                 if(advmask(ji,jj,jk) .NE. 0) then

            DO ju=1, dimen_jarr2

               ji = jarr2(1, ju)
               jj = jarr2(2, ju)
               jk = jarr2(3, ju)
               junk  = zti(ji,jj,jk,A)
               junki = zti(ji + 1,jj,jk,A)
               junkj = zti(ji,jj + 1,jk,A)
               junkk = zti(ji,jj,jk - 1,A)
!                zx(ji,jj,jk,A) = ( abs(zaa(ji,jj,jk)) - zdt *zaa(ji,jj,jk)**2/(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) ) )*inv_eu(ji,jj,jk) )
!                zy(ji,jj,jk,A) = ( abs(zbb(ji,jj,jk)) - zdt*zbb(ji,jj,jk)**2/(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) ) )*inv_ev(ji,jj,jk) )
!                zz(ji,jj,jk) = ( abs(zcc(ji,jj,jk)) - zdt*zcc(ji,jj,jk)**2/(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk) ) )*inv_et(ji,jj,jk) )
               zx(ji,jj,jk,A) = big_fact_zaa(ji,jj,jk)*(junki - junk)/(junk + junki + rtrn)* rsc
               zy(ji,jj,jk,A) = big_fact_zbb(ji,jj,jk)*(junkj - junk)/(junk + junkj + rtrn)* rsc
               zz(ji,jj,jk,A) = big_fact_zcc(ji,jj,jk)*(junk - junkk)/(junk + junkk + rtrn)* rsc*( -1.)

           END DO
!                 endif

          END IF

!$omp end parallel


C---------------------------------------------------------------------------
C 2.4 cross terms
C
CCC          IF (crosster) THEN
CCC              DO jk = 2,jpkm1
CCC                DO jj = 2,jpjm1
CCC                 DO ji = 2,jpim1
CCC                    zx(ji,jj,jk) = zx(ji,jj,jk)
CCC     $                  - 0.5*zdt*rsc*zaa(ji,jj,jk)*0.25*
CCC     $                  (    (zbb(ji  ,jj - 1,jk  ) + zbb(ji + 1,jj - 1
CCC     $                  ,jk  ) + zbb(ji + 1,jj  ,jk  ) + zbb(ji  ,jj
CCC     $                  ,jk))* (zti(ji  ,jj + 1,jk  ) + zti(ji + 1,jj +
CCC     $                  1,jk  ) - zti(ji + 1,jj - 1,jk  ) - zti(ji  ,jj
CCC     $                  - 1,jk  ))/ (zti(ji  ,jj + 1,jk  ) + zti(ji + 1
CCC     $                  ,jj + 1,jk  ) + zti(ji + 1,jj - 1,jk  ) + zti(ji
CCC     $                  ,jj - 1,jk  ) + rtrn) + (zcc(ji  ,jj  ,jk  ) +
CCC     $                  zcc(ji + 1,jj  ,jk  ) + zcc(ji  ,jj  ,jk + 1) +
CCC     $                  zcc(ji + 1,jj  ,jk + 1))* (zti(ji  ,jj  ,jk - 1)
CCC     $                  + zti(ji + 1,jj  ,jk - 1) - zti(ji  ,jj  ,jk + 1
CCC     $                  )- zti(ji + 1,jj  ,jk + 1))/ (zti(ji  ,jj  ,jk -
CCC     $                  1) + zti(ji + 1,jj  ,jk - 1) + zti(ji  ,jj  ,jk
CCC     $                  +1) + zti(ji + 1,jj  ,jk + 1) + rtrn))/(e1u(ji
CCC     $                  ,jj)*e2u(ji,jj)*fse3u(ji,jj,jk))*vmask(ji  ,jj -
CCC     $                  1,jk  )*vmask(ji + 1,jj - 1,jk  )*vmask(ji + 1
CCC     $                  ,jj,jk)*vmask(ji  ,jj  ,jk  )*tmask(ji  ,jj  ,jk
CCC     $                  )*tmask(ji + 1,jj  ,jk  )*tmask(ji  ,jj  ,jk + 1
CCC     $                  )*tmask(ji + 1,jj  ,jk + 1)
CCC                    zy(ji,jj,jk) = zy(ji,jj,jk)
CCC     $                  - 0.5*zdt*rsc*zbb(ji,jj,jk)*0.25*
CCC     $                  (    (zaa(ji - 1,jj  ,jk  ) + zaa(ji - 1,jj + 1
CCC     $                  ,jk  ) + zaa(ji  ,jj  ,jk  ) + zaa(ji  ,jj + 1
CCC     $                  ,jk))* (zti(ji + 1,jj + 1,jk  ) + zti(ji + 1,jj
CCC     $                  ,jk  ) - zti(ji - 1,jj + 1,jk  ) - zti(ji - 1,jj
CCC     $                  ,jk  ))/ (zti(ji + 1,jj + 1,jk  ) + zti(ji + 1
CCC     $                  ,jj  ,jk  ) + zti(ji - 1,jj + 1,jk  ) + zti(ji
CCC     $                  - 1,jj  ,jk  ) + rtrn) + (zcc(ji  ,jj  ,jk  )
CCC     $                  + zcc(ji  ,jj  ,jk + 1) + zcc(ji  ,jj + 1,jk  )
CCC     $                  + zcc(ji  ,jj + 1,jk + 1))* (zti(ji  ,jj  ,jk -
CCC     $                  1) + zti(ji  ,jj + 1,jk - 1) - zti(ji  ,jj  ,jk
CCC     $                  +1) - zti(ji  ,jj + 1,jk + 1))/ (zti(ji  ,jj
CCC     $                  ,jk- 1) + zti(ji  ,jj + 1,jk - 1) + zti(ji  ,jj
CCC     $                  ,jk+ 1) + zti(ji  ,jj + 1,jk + 1) + rtrn))
CCC     $                  /(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk))
CCC     $                  *umask(ji - 1,jj,jk  )*umask(ji - 1,jj + 1,jk  )
CCC     $                  *umask(ji  ,jj,jk  )*umask(ji  ,jj + 1,jk  )
CCC     $                  *tmask(ji  ,jj,jk)*tmask(ji  ,jj  ,jk + 1)
CCC     $                  *tmask(ji  ,jj + 1,jk)*tmask(ji  ,jj + 1,jk + 1)
CCC                    zz(ji,jj,jk) = zz(ji,jj,jk)
CCC    $                  - 0.5*zdt*rsc*zcc(ji,jj,jk)*0.25*
CCC     $                  (    (zaa(ji - 1,jj  ,jk  ) + zaa(ji  ,jj  ,jk
CCC     $                  ) + zaa(ji  ,jj  ,jk - 1) + zaa(ji - 1,jj  ,jk -
CCC     $                  1))*(zti(ji + 1,jj  ,jk - 1) + zti(ji + 1,jj
CCC     $                  ,jk  ) - zti(ji - 1,jj  ,jk  ) - zti(ji - 1,jj
CCC     $                  ,jk - 1))/(zti(ji + 1,jj  ,jk - 1) + zti(ji + 1
CCC     $                  ,jj,jk  ) + zti(ji - 1,jj  ,jk  ) + zti(ji - 1
CCC     $                  ,jj,jk - 1) + rtrn) + (zbb(ji  ,jj - 1,jk  )
CCC     $                  + zbb(ji  ,jj  ,jk  ) + zbb(ji  ,jj  ,jk - 1)
CCC     $                  + zbb(ji  ,jj - 1,jk - 1))*(zti(ji  ,jj + 1,jk -
CCC     $                  1) + zti(ji  ,jj + 1,jk  ) - zti(ji  ,jj - 1,jk
CCC     $                  ) - zti(ji  ,jj - 1,jk - 1))/(zti(ji  ,jj + 1,jk
CCC     $                  - 1) + zti(ji  ,jj + 1,jk  ) + zti(ji  ,jj - 1
CCC     $                  ,jk  ) + zti(ji  ,jj - 1,jk - 1) + rtrn))
CCC     $                  /(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk))
CCC     $                  *umask(ji - 1,jj,jk  )*umask(ji  ,jj  ,jk  )
CCC     $                  *umask(ji  ,jj,jk- 1)*umask(ji - 1,jj  ,jk - 1)
CCC     $                  *vmask(ji  ,jj- 1,jk)*vmask(ji  ,jj  ,jk  )
CCC     $                  *vmask(ji  ,jj  ,jk-1)*vmask(ji  ,jj - 1,jk - 1)
CCC                  END DO
CCC               END DO
CCC              END DO
CCCC
CCC              DO jj = 2,jpjm1
CCC                DO ji = 2,jpim1
CCC                  zx(ji,jj,1) = zx(ji,jj,1)
CCC     $                - 0.5*zdt*rsc*zaa(ji,jj,1)*0.25*
CCC     $                ( (zbb(ji  ,jj - 1,1  ) + zbb(ji + 1,jj - 1,1  )
CCC     $                + zbb(ji + 1,jj  ,1  ) + zbb(ji  ,jj  ,1  ))
CCC     $                *(zti(ji  ,jj + 1,1  ) + zti(ji + 1,jj + 1,1  )
CCC     $                - zti(ji + 1,jj - 1,1  ) - zti(ji  ,jj - 1,1  ))
CCC     $                /(zti(ji  ,jj + 1,1  ) + zti(ji + 1,jj + 1,1  )
CCC     $                + zti(ji + 1,jj - 1,1  ) + zti(ji  ,jj - 1,1  ) +
CCC     $                rtrn))/(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,1))
CCC     $                *vmask(ji  ,jj - 1,1  )*vmask(ji + 1,jj - 1,1  )
CCC     $                *vmask(ji + 1,jj  ,1  )*vmask(ji  ,jj  ,1  )
CCC                  zy(ji,jj,1) = zy(ji,jj,1)
CCC     $                - 0.5*zdt*rsc*zbb(ji,jj,1)*0.25*
CCC     $                ( (zaa(ji - 1,jj  ,1  ) + zaa(ji - 1,jj + 1,1  )
CCC     $                + zaa(ji  ,jj  ,1  ) + zaa(ji  ,jj + 1,1  ))
CCC     $                *(zti(ji + 1,jj + 1,1  ) + zti(ji + 1,jj  ,1  )
CCC     $                - zti(ji - 1,jj + 1,1  ) - zti(ji - 1,jj  ,1  ))
CCC     $                /(zti(ji + 1,jj + 1,1  ) + zti(ji + 1,jj  ,1  )
CCC     $                + zti(ji - 1,jj + 1,1  ) + zti(ji - 1,jj  ,1  ) +
CCC     $                rtrn))/(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,1))
CCC     $                *umask(ji - 1,jj  ,1  )*umask(ji - 1,jj + 1,1  )
CCC     $                *umask(ji  ,jj  ,1  )*umask(ji  ,jj + 1,1  )
CCC                END DO
CCC              END DO
CCC          ENDIF
C-----------------------------------------------------------------------------



C ... Lateral boundary conditions on z[xyz]
#ifdef key_mpp

C   ... Mpp : export boundary values to neighboring processors




C         DO itid = 1, ntids
C     IF( itid - 1 + jn <= jptra ) THEN
C               CALL mpplnk( zx(:,:,:,itid), 1, 1 )
C               CALL mpplnk( zy(:,:,:,itid), 1, 1 )
C               CALL mpplnk( zz(:,:,:,itid), 1, 1 )
C     END IF
C  END DO

      IF( ntids - 1 + jn <= jptra ) THEN
       pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
      END IF

        CALL mpplnk_my(zx, pack_size,1,1)
        CALL mpplnk_my(zy, pack_size,1,1)
        CALL mpplnk_my(zz, pack_size,1,1)


#else
C
C   ... T-point, 3D array, full local array z[xyz] are initialised
C
C
          DO itid = 1, ntids
         IF( itid - 1 + jn <= jptra ) THEN
                CALL lbc( zx(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zy(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zz(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
         END IF
      END DO
C
#endif
C
C
C 2.4 reinitialization
C

!$omp     parallel default(none) private(mytid,A,junk,jk,jj,ji)
!$omp&      shared(zkx,zky,zkz,zz,zx,zy,zti,jpjm1,jpim1,dimen_jarr1,jarr1,jpi,jpk,jpj,jn)

#ifdef __OPENMP
          mytid = omp_get_thread_num()  ! take the thread ID
#endif


      IF( mytid + jn <= jptra ) THEN
      A=mytid+1
C
C            2.5 calcul of the final field:
C                advection by antidiffusive mass fluxes and an upstream scheme
C
             jk = 1
             DO jj = 2,jpjm1
                 DO ji = 2,jpim1
                   junk  = zti(ji,jj,jk,A)
                   zkx(ji,jj,jk,A) = fsx(junk,zti(ji + 1,jj,jk,A),zx(ji,jj,jk,A))!zti(ji,jj,jk,A),zti(ji + 1,jj,jk,A),zaa(ji,jj,jk))
                   zky(ji,jj,jk,A) = fsy(junk,zti(ji,jj + 1,jk,A),zy(ji,jj,jk,A))!zti(ji,jj,jk,A),zti(ji,jj + 1,jk,A),zbb(ji,jj,jk))
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

!$omp end parallel



C ... Lateral boundary conditions on zk[xy]
#ifdef key_mpp

C   ... Mpp : export boundary values to neighboring processors


C        DO itid = 1, ntids
C    IF( itid - 1 + jn <= jptra ) THEN
C              CALL mpplnk( zkx(:,:,:,itid), 1, 1 )
C              CALL mpplnk( zky(:,:,:,itid), 1, 1 )
C    END IF
C END DO

        IF( ntids - 1 + jn <= jptra ) THEN
           pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
        END IF

        CALL mpplnk_my(zkx, pack_size,1,1)
        CALL mpplnk_my(zky, pack_size,1,1)


!        trcadvparttime = MPI_WTIME() - trcadvparttime
!        trcadvtottime = trcadvtottime + trcadvparttime

#else
C
C   ... T-point, 3D array, full local array zk[xy] are initialised
C
C
      DO itid = 1, ntids
        IF( itid - 1 + jn <= jptra ) THEN
               CALL lbc( zkx(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
               CALL lbc( zky(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
        END IF
      END DO
C
#endif

!$omp    parallel default(none) private(mytid,A,B,zbtr,ji,jj,jk,ju,jf)
!$omp&      shared(zkx,zky,zkz,zbtr_arr,e1t,e2t,ztj,dimen_jarr3,jarr3,ncor,jn,
!$omp&             jarr_adv_flx,Fsize,diaflx)

#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif

C        2.6. calcul of after field using an upstream advection scheme


         IF( mytid + jn <= jptra ) THEN
         A = mytid+1
         B = mytid+jn
            if(ncor .EQ. 1) then
               DO ju=1, dimen_jarr3

                  ji = jarr3(1, ju)
                  jj = jarr3(2, ju)
                  jk = jarr3(3, ju)
                  jf = jarr_adv_flx(ju)

!                 if(zbtr_arr(ji,jj,jk) .NE. 0) then zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
                  zbtr = zbtr_arr(ji,jj,jk)
                  ztj(ji,jj,jk,A) = -zbtr*( zkx(ji,jj,jk,A) - zkx(ji - 1,jj,jk,A)
     $              + zky(ji,jj,jk,A) - zky(ji,jj - 1,jk,A)+ zkz(ji,jj,jk,A) - zkz(ji,jj,jk + 1,A) )+ ztj(ji,jj,jk,A)
#if defined key_trc_diatrd
                  trtrd(ji,jj,jk,B,1) = trtrd(ji,jj,jk,B,1) -zbtr*( zkx(ji,jj,jk,A) - zkx(ji - 1,jj,jk,A) )
                  trtrd(ji,jj,jk,B,2) = trtrd(ji,jj,jk,B,2) -zbtr*( zky(ji,jj,jk,A) - zky(ji,jj - 1,jk,A) )
                  trtrd(ji,jj,jk,B,3) = trtrd(ji,jj,jk,B,3) -zbtr*( zkz(ji,jj,jk,A) - zkz(ji,jj,jk + 1,A) )
#endif
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

!                if(zbtr_arr(ji,jj,jk) .NE. 0) then    zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
                  zbtr = zbtr_arr(ji,jj,jk)
                  ztj(ji,jj,jk,A) = -zbtr*( zkx(ji,jj,jk,A) - zkx(ji - 1,jj,jk,A)
     $              + zky(ji,jj,jk,A) - zky(ji,jj - 1,jk,A)+ zkz(ji,jj,jk,A) - zkz(ji,jj,jk + 1,A) )
#if defined key_trc_diatrd
                  trtrd(ji,jj,jk,B,1) = trtrd(ji,jj,jk,B,1) -zbtr*( zkx(ji,jj,jk,A) - zkx(ji - 1,jj,jk,A) )
                  trtrd(ji,jj,jk,B,2) = trtrd(ji,jj,jk,B,2) -zbtr*( zky(ji,jj,jk,A) - zky(ji,jj - 1,jk,A) )
                  trtrd(ji,jj,jk,B,3) = trtrd(ji,jj,jk,B,3) -zbtr*( zkz(ji,jj,jk,A) - zkz(ji,jj,jk + 1,A) )
#endif
!     Save advective fluxes x,y,z
                 IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                    diaflx(jf,B,1) = diaflx(jf,B,1) + zkx(ji,jj,jk,A)
                    diaflx(jf,B,2) = diaflx(jf,B,2) + zky(ji,jj,jk,A)
                    diaflx(jf,B,3) = diaflx(jf,B,3) + zkz(ji,jj,jk,A)
                 END IF
              END DO



           endif

      END IF

!$omp end parallel

        END DO ANTIDIFF_CORR


C       3. trend due to horizontal and vertical advection of tracer jn
C       --------------------------------------------------------------
C
!$omp   parallel default(none) private(mytid,A,B,ji,jj,jk,ju) shared(ncor,dimen_jarrt,jarrt,tra,ztj,jn)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
        A=mytid+1
        B=mytid+jn

        IF( mytid + jn <= jptra ) THEN


           if(ncor .EQ. 1) then
              DO ju=1, dimen_jarrt
                 ji = jarrt(1, ju)
                 jj = jarrt(2, ju)
                 jk = jarrt(3, ju)

                 tra(ji,jj,jk,B) = tra(ji,jj,jk,B)+ ztj(ji,jj,jk,A)!+ (zbuf(ji,jj,jk) + ztj(ji,jj,jk,A))

              END DO
           else
              DO ju=1, dimen_jarrt
                 ji = jarrt(1, ju)
                 jj = jarrt(2, ju)
                 jk = jarrt(3, ju)

                 tra(ji,jj,jk,B) = tra(ji,jj,jk,B)+ (zbuf(ji,jj,jk) + ztj(ji,jj,jk,A))

              END DO
           endif

      END IF


!$omp end parallel





C 4.0 convert the transport trend into advection trend
C ----------------------------------------------------

#if defined key_trc_diatrd

      DO ju=1, dimen_jarr

          ji = jarr(1, ju)
          jj = jarr(2, ju)
          jk = jarr(3, ju)

    !     if(zbtr_arr(ji,jj,jk) .NE. 0) then!              zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
          zbtr =zbtr_arr(ji,jj,jk)
          zgm =zbtr*trn(ji,jj,jk,B)*(un(ji,jj,jk)*e2u(ji,jj)*fse3u(ji,jj,jk)-un(ji-1,jj,jk)*e2u(ji-1,jj)*fse3u(ji-1,jj,jk))
          zgz =zbtr*trn(ji,jj,jk,B)*(vn(ji,jj,jk)*e1v(ji,jj)*fse3v(ji,jj,jk)-vn(ji,jj-1,jk)*e1v(ji,jj-1)*fse3v(ji,jj-1,jk))
          trtrd(ji,jj,jk,B,1) = trtrd(ji,jj,jk,B,1) + zgm
          trtrd(ji,jj,jk,B,2) = trtrd(ji,jj,jk,B,2) + zgz
          trtrd(ji,jj,jk,B,3) = trtrd(ji,jj,jk,B,3)- trn(ji,jj,jk,B) * hdivn(ji,jj,jk)

      END DO


C Lateral boundary conditions on trtrd:
#   ifdef key_mpp



        CALL mpplnk( trtrd(1,1,1,B,1), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,B,2), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,B,3), 1, 1 )



#   else
        CALL lbc( trtrd(1,1,1,B,1), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,B,2), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,B,3), 1, 1, 1, 1, jpk, 1 )
#   endif

#endif

       END DO TRACER_LOOP

       trcadvparttime = MPI_WTIME() - trcadvparttime
       trcadvtottime = trcadvtottime + trcadvparttime


      END SUBROUTINE smolar
