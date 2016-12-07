      SUBROUTINE trcadv
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcadv
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Compute the now trend due to the advection of tracers
!!!     (tr) and add it to the general trend of passive tracer equations.
!!!
!!!
! CC----------------------------------------------------------------------
! CC parameters and commons
! CC ======================

       USE myalloc
       ! epascolo USE myalloc_mpp
       USE ADV_mem
       USE DIA_mem
       use mpi
       use omp_lib
       USE ogstm_mpi_module
       
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
      INTEGER :: jk,jj,ji,jt,jn,jf,ju
      INTEGER :: jp,pack_size
      double precision :: zbtr,zdt
      double precision :: junk, junki, junkj, junkk
      double precision :: timer
      double precision,dimension(:),allocatable :: array 
      double precision,dimension(:,:),allocatable :: surface 
      INTEGER :: A,B

!-------------------------------------------------------------------

      MPI_CHECK = .FALSE.

      if(allpoints .EQ. 0) then  ! INIT phase

         zti(:,:,:,:)  = 0.
         ztj(:,:,:,:)  = 0.
         zkx(:,:,:,:)  = 0.
         zky(:,:,:,:)  = 0.
         zkz(:,:,:,:)  = 0.
         zbuf(:,:,:,:) = 0.
         zx(:,:,:,:)   = 0.
         zy(:,:,:,:)   = 0.
         zz(:,:,:,:)   = 0.

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
               DO ji = 2,jpim1
            DO jj = 2,jpjm1
                  !dir$ vector aligned
         DO jk = 1,jpkm1
                  zbtr_arr(jk,jj,ji) = 1./(e1t(jj,ji)*e2t(jj,ji)*e3t(jk,jj,ji))
               END DO
            END DO
         END DO

               DO ji = 2,jpim1
            DO jj = 2,jpjm1
            !dir$ vector aligned
         DO jk = 1,jpkm1
                  inv_eu(jk,jj,ji) = 1./(e1u(jj,ji)*e2u(jj,ji)*e3u(jk,jj,ji) )
                  inv_ev(jk,jj,ji) = 1./(e1v(jj,ji)*e2v(jj,ji)*e3v(jk,jj,ji) )
               END DO
            END DO
         END DO

               DO ji = 2,jpim1
            DO jj = 2,jpjm1
            !dir$ vector aligned
         DO jk = 2,jpkm1
                  inv_et(jk,jj,ji) = 1./(e1t(jj,ji)*e2t(jj,ji)*e3w(jk,jj,ji) )
               END DO
            END DO
         END DO

         allpoints = 0

               DO ji = 1,jpi
            DO jj = 1,jpj
         DO jk = 1,jpk
                  allpoints = allpoints + 1
                  if(tmask(jk,jj,ji) .NE. 0) then
                     tpoints = tpoints + 1
                  endif
               END DO
            END DO
         END DO

         write(*,*) 'trcadv: RANK -> ', myrank, ' all_points -> ', allpoints

         goodpoints = 0
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
                              DO myjk=jk+jklef, jk+jkrig
                              locsum = locsum + tmask(myjk,myjj,myji)
                              END DO
                        END DO
                        END DO
                        if(locsum .NE. 0) then
                        goodpoints = goodpoints + 1
                        advmask(jk,jj,ji) = 1
                        else
                        advmask(jk,jj,ji) = 0
                        endif
               END DO
            END DO
         END DO

         write(*,*) 'trcadv: RANK -> ', myrank, ' good_points -> ', goodpoints

         DO  ji = 2,jpim1
            DO jj = 2,jpjm1
                  DO jk = 1,jpk
                        if(advmask(jk,jj,ji) .NE. 0) then
                        zbtr_arr(jk,jj,ji) = 1./(e1t(jj,ji)*e2t(jj,ji)*e3t(jk,jj,ji))
                        dimen_jarr = dimen_jarr + 1
                        jarr(3,dimen_jarr) = ji
                        jarr(2,dimen_jarr) = jj
                        jarr(1,dimen_jarr) = jk
                        else
                        zbtr_arr(jk,jj,ji) = 0.
                        endif
               END DO
            END DO
         END DO

      !    DO  ji = 2,jpim1
      !       DO jj = 2,jpjm1
      !             DO jk = 2,jpk
      !                   if(advmask(jk,jj,ji) .NE. 0) then
      !                   dimen_jarr1 = dimen_jarr1 + 1
      !                   jarr1(3,dimen_jarr1) = ji
      !                   jarr1(2,dimen_jarr1) = jj
      !                   jarr1(1,dimen_jarr1) = jk
      !                   endif
      !          END DO
      !       END DO
      !    END DO

      !    DO ji = 2,jpim1
      !       DO jj = 2,jpjm1
      !             DO jk = 2,jpkm1
      !                   if(advmask(jk,jj,ji) .NE. 0) then
      !                   dimen_jarr2 = dimen_jarr2 + 1
      !                   jarr2(3,dimen_jarr2) = ji
      !                   jarr2(2,dimen_jarr2) = jj
      !                   jarr2(1,dimen_jarr2) = jk
      !                   endif
      !          END DO
      !       END DO
      !    END DO

               DO ji = 2,jpim1
            DO jj = 2,jpjm1
         DO jk = 1,jpkm1
                  if(advmask(jk,jj,ji) .NE. 0) then
                     dimen_jarr3 = dimen_jarr3 + 1
                     jarr3(3,dimen_jarr3) = ji
                     jarr3(2,dimen_jarr3) = jj
                     jarr3(1,dimen_jarr3) = jk
                  endif
               END DO
            END DO
         END DO

               DO ji = 1,jpi
            DO jj = 1,jpj
         DO jk = 1,jpk
                  if(tmask(jk,jj,ji) .NE. 0) then
                     dimen_jarrt = dimen_jarrt + 1
                     jarrt(3,dimen_jarrt) = ji
                     jarrt(2,dimen_jarrt) = jj
                     jarrt(1,dimen_jarrt) = jk
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

      endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end initialization phase

         jk=1
            !ALLOCATE(surface(jpj,jpi))
           ! timer = omp_get_wtime()
               DO ji = 2,jpim1
               !dir$ vector aligned
            DO jj = 2,jpjm1
                  zbtr_arr(1,jj,ji) = 1./(e1t(jj,ji)*e2t(jj,ji)*e3t(1,jj,ji))
               END DO
            END DO
            !zbtr_arr(1,:,:) = surface(:,:)

               DO ji = 2,jpim1
               !dir$ vector aligned
            DO jj = 2,jpjm1
                  inv_eu(1,jj,ji) = 1./(e1u(jj,ji)*e2u(jj,ji)*e3u(1,jj,ji) )
               END DO
            END DO

                  DO ji = 2,jpim1
                 !dir$ vector aligned
            DO jj = 2,jpjm1
                  inv_ev(1,jj,ji) = 1./(e1v(jj,ji)*e2v(jj,ji)*e3v(1,jj,ji) )
               END DO
            END DO

               DO ji = 2,jpim1
               !dir$ vector aligned
            DO jj = 2,jpjm1
                  inv_et(1,jj,ji) = 1./(e1t(jj,ji)*e2t(jj,ji)*e3w(1,jj,ji) )
               END DO
            END DO

            !DEALLOCATE(surface)
            ! timer = omp_get_wtime() -timer
            ! print *,"TTT  = ",timer

      zdt = rdt*ndttrc


!!!&omp  parallel default(none) private(A,mytid,jk,jj,ji) shared(jpk,jpj,jpi,ztj,zx,zy,zz)

       A = 1
       ztj(:,:,:,:) = 0
       zx(:,:,:,:) = 0
       zy(:,:,:,:) = 0
       zz(:,:,:,:) = 0

!!!&omp parallel default(none) private(mytid, jj,ji)
!!!&omp&         shared(jk, jpk,jpj,jpi,big_fact_zaa,big_fact_zbb,big_fact_zcc,zaa,zbb,zcc,inv_eu,inv_ev,inv_et,
!!!&omp&                un,vn,wn,e2u,e3u,e3v,e1v,e1t,e2t,e3t,zdt )

             DO ji = 1,jpi
                  DO jj = 1,jpj
                  !dir$ vector aligned
                        DO jk = 1,jpk
                        zaa(jk,jj,ji) = e2u(jj,ji)*e3u(jk,jj,ji) * un(jk,jj,ji)
                        END DO
                  END DO
            END DO

       DO ji = 1,jpi
          DO jj = 1,jpj
          !dir$ vector aligned
             DO jk = 1,jpk
               big_fact_zaa(jk,jj,ji) = ( abs(zaa(jk,jj,ji)) - zdt*zaa(jk,jj,ji)**2*inv_eu(jk,jj,ji) )!/(e1u(jj,ji)*e2u(jj,ji)*e3t(jk,jj,ji) ) )
             END DO
          END DO
       END DO

       DO ji = 1,jpi
          DO jj = 1,jpj
          !dir$ vector aligned
              DO jk = 1,jpk
                zbb(jk,jj,ji) = e1v(jj,ji)*e3v(jk,jj,ji) * vn(jk,jj,ji)
                END DO
          END DO
       END DO
      
       DO ji = 1,jpi
         DO jj = 1,jpj
         !dir$ vector aligned
            DO jk = 1,jpk
                big_fact_zbb(jk,jj,ji) = ( abs(zbb(jk,jj,ji)) - zdt*zbb(jk,jj,ji)**2*inv_ev(jk,jj,ji) )
             END DO
          END DO
       END DO

       DO ji = 1,jpi
         DO jj = 1,jpj
         !dir$ vector aligned
          DO jk = 1,jpk
                zcc(jk,jj,ji) = e1t(jj,ji)*e2t(jj,ji)* wn(jk,jj,ji)
             END DO
          END DO
       END DO

       DO ji = 1,jpi
         DO jj = 1,jpj
         !dir$ vector aligned
            DO jk = 1,jpk
                big_fact_zcc(jk,jj,ji) = ( abs(zcc(jk,jj,ji)) - zdt*zcc(jk,jj,ji)**2*inv_et(jk,jj,ji) )
             END DO
          END DO
       END DO


!!     tracer loop parallelized (macrotasking)
!!     =======================================

      trcadvparttime = MPI_WTIME()

       
      TRACER_LOOP: DO  jn = 1, jptra


!!        1. tracer flux in the 3 directions
!!        ----------------------------------
!!        1.1 mass flux at u v and t-points and initialization
!!       1.2 calcul of intermediate field with an upstream advection scheme
!!           and mass fluxes calculated above
!!       calcul of tracer flux in the i and j direction

!!!&omp   parallel default(none) private(jk,jj,ji,mytid,A,B,junk)
!!!&omp&      shared(jpk,jpj,jpi,zkx,zky,zbb,zaa,zcc,jarr1,dimen_jarr1,jn,zkz,trn,jpjm1,jpim1)


      A =  1
      B =  jn
      
       zkx(  :,:,1,A)=0.  
       zkx(:,:,jpi,A)=0.    
       zky(:,  1,:,A)=0.  
       zky(:,jpj,:,A)=0.  
       zkz(1,:,:,A)  =0.
! loop unfusion

        DO ji = 2,jpim1
           !dir$ vector aligned
           DO jj = 2,jpjm1
                 zkx(1,jj,ji,A) = fsx(trn(1,jj,ji,B),trn(1,jj,ji + 1,B),zaa(1,jj,ji))
              END DO
           END DO
        

        DO ji = 2,jpim1
          !dir$ vector aligned
           DO jj = 2,jpjm1
                zky(1,jj,ji,A) = fsy(trn(1,jj,ji,B),trn(1,jj+1,ji,B),zbb(1,jj,ji))
           END DO
        END DO

       DO ji = 1,jpi
           !dir$ vector aligned
           DO jk = 2,jpk
            zkz(jk,1,ji,A) = fsz(trn(jk,1,ji,B),trn(jk-1,1,ji,B),zcc(jk,1,ji))
            ENDDO
       ENDDO
      
       DO ji = 1,jpi
            !dir$ vector aligned
            DO jk = 2,jpk
            zkz(jk,jpj,ji,A) = fsz(trn(jk,jpj,ji,B),trn(jk-1,jpj,ji,B),zcc(jk,jpj,ji))
            END DO
       END DO
! loop unfusion
      DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
            zkz(jk,jj,1,A) = fsz(trn(jk,jj,1,B),trn(jk-1,jj,1,B),zcc(jk,jj,1))
            END DO
      END DO

      DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
            zkz(jk,jj,jpi,A) = fsz(trn(jk,jj,jpi,B),trn(jk-1,jj,jpi,B),zcc(jk,jj,jpi))
            END DO
      END DO

      DO  ji = 2,jpim1
        DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
              zkx(jk,jj,ji,A) = fsx(trn(jk,jj,ji,B),trn(jk,jj,ji + 1,B),zaa(jk,jj,ji))*advmask(jk,jj,ji)
            END DO
         END DO
      END DO

            DO  ji = 2,jpim1
        DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
              zky(jk,jj,ji,A) = fsy(trn(jk,jj,ji,B),trn(jk,jj + 1,ji,B),zbb(jk,jj,ji))*advmask(jk,jj,ji)
              END DO
         END DO
      END DO

            DO  ji = 2,jpim1
        DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
              zkz(jk,jj,ji,A) = fsz(trn(jk,jj,ji,B),trn(jk-1,jj,ji,B),zcc(jk,jj,ji))*advmask(jk,jj,ji)
            END DO
         END DO
      END DO

! ! epascolo mpi comment
! !! ... Lateral boundary conditions on zk[xy]
! #ifdef key_mpp

! !!   ... Mpp : export boundary values to neighboring processors


!       IF( ntids - 1 + jn <= jptra ) THEN
!        pack_size = ntids
!         ELSE
!            pack_size = ntids - (ntids - 1 + jn - jptra)
!       END IF
!         CALL mpplnk_my(zkx, pack_size,1,1)
!         CALL mpplnk_my(zky, pack_size,1,1)




! #else

! !!   ... T-point, 3D array, full local arrays zk[xy] are initialised


!         DO itid = 1, ntids

!        IF( itid - 1 + jn <= jptra ) THEN

               CALL lbc( zkx(:,:,:,1), 1, 1, 1, 1, jpk, 1 )
               CALL lbc( zky(:,:,:,1), 1, 1, 1, 1, jpk, 1 )

!        END IF

!         END DO

! #endif


!! 2. calcul of after field using an upstream advection scheme
!! -----------------------------------------------------------


!!!&omp   parallel default(none) private(mytid,A,B,zbtr,jk,jj,ji,ju,jf)
!!!&omp&      shared(zkx,zky,zkz,zti,jpim1,jpjm1,trn,zdt,jn,jpkm1,zbtr_arr,e1t,e2t,ztj,jarr3,ncor,dimen_jarr3,
!!!&omp&             jarr_adv_flx,Fsize,diaflx)



      
      A = 1
      B = jn
           DO ju=1, dimen_jarr3

              ji = jarr3(3, ju)
              jj = jarr3(2, ju)
              jk = jarr3(1, ju)
              jf = jarr_adv_flx(ju)

              zbtr = zbtr_arr(jk,jj,ji)

              ztj(jk,jj,ji,A) = -zbtr* &
     &          ( zkx(jk,jj,ji,A) - zkx(jk,jj,ji-1,A) &
     &          + zky(jk,jj,ji,A) - zky(jk,jj- 1,ji,A) &
     &          + zkz(jk,jj,ji,A) - zkz(jk+1,jj,ji,A) )


              IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                    diaflx(jf,B,3) = diaflx(jf,B,3) + zkx(jk,jj,ji,A)
                    diaflx(jf,B,2) = diaflx(jf,B,2) + zky(jk,jj,ji,A)
                    diaflx(jf,B,1) = diaflx(jf,B,1) + zkz(jk,jj,ji,A)
              END IF

            END DO

!! 2.1 start of antidiffusive correction loop

        ANTIDIFF_CORR: DO jt = 1,ncor

!! 2.2 calcul of intermediary field zti


!!!&omp     parallel default(none) private(mytid,A,B,jk,jj,ji)
!!!&omp&       shared(jt,jn,ncor,jpkm1,jpjm1,jpim1,zti,ztj,trn,zdt,zbuf)

      A = 1
      B = jn
             if(jt .EQ. 1) then

                if(ncor .EQ. 1) then

                         DO ji = 2,jpim1
                      DO jj = 2,jpjm1
                      !dir$ vector aligned
                   DO jk = 1,jpkm1
                            zti(jk,jj,ji,A) = trn(jk,jj,ji,B) + zdt*ztj(jk,jj,ji,A)
                        END DO
                      END DO
                   END DO

                else
                         DO ji = 2,jpim1
                      DO jj = 2,jpjm1
                      !dir$ vector aligned
                    DO jk = 1,jpkm1
                            zti(jk,jj,ji,A) = trn(jk,jj,ji,B) + zdt*ztj(jk,jj,ji,A)
                          !  zbuf(jk,jj,ji,A) = ztj(jk,jj,ji,A)
                         END DO
                      END DO
                   END DO

                        DO ji = 2,jpim1
                      DO jj = 2,jpjm1
                      !dir$ vector aligned
                    DO jk = 1,jpkm1
                         !   zti(jk,jj,ji,A) = trn(jk,jj,ji,B) + zdt*ztj(jk,jj,ji,A)
                            zbuf(jk,jj,ji,A) = ztj(jk,jj,ji,A)
                         END DO
                      END DO
                   END DO
                
                
                endif

             else

                      DO ji = 2,jpim1
                   DO jj = 2,jpjm1
                   !dir$ vector aligned
                DO jk = 1,jpkm1
                         zti(jk,jj,ji,A) =  zti(jk,jj,ji,A) + zdt*ztj(jk,jj,ji,A)
                      END DO
                   END DO
                END DO

                      DO ji = 2,jpim1
                   DO jj = 2,jpjm1
                   !dir$ vector aligned
                DO jk = 1,jpkm1
                         zbuf(jk,jj,ji,A) = zbuf(jk,jj,ji,A) + ztj(jk,jj,ji,A)
                      END DO
                   END DO
                END DO

             endif


!! ... Lateral boundary conditions on zti
! epascolo mpi comment
! #ifdef key_mpp

! !!   ... Mpp : export boundary values to neighboring processors



!         IF( ntids - 1 + jn <= jptra ) THEN
!            pack_size = ntids
!         ELSE
!            pack_size = ntids - (ntids - 1 + jn - jptra)
!         END IF

!         CALL mpplnk_my(zti, pack_size,1,1)



! #else

! !!   ... T-point, 3D array, full local array zti is initialised

!           DO itid = 1, ntids
!          IF( itid - 1 + jn <= jptra ) THEN
                 CALL lbc( zti(:,:,:,1), 1, 1, 1, 1, jpk, 1 )
!          END IF
!           END DO

! #endif


!! 2.3 calcul of the antidiffusive flux

!!!&omp     parallel default(none) private(mytid,A,junk, junki, junkj, junkk, jk,jj,ji)
!!!&omp&       shared(jn,jpkm1,jpjm1,jpim1,zti,ztj,zy,zx,zz,jarr2,big_fact_zbb,
!!!&omp&              big_fact_zaa,big_fact_zcc,dimen_jarr2,rtrn,rsc)


      
      A = 1
      !jk = 1
!          DO jk = 1,jpkm1
      DO ji = 2,jpim1
            DO jj = 2,jpjm1
                  junk  = zti(1,jj,ji,A)
                  junki = zti(1,jj,ji+1,A)
                  junkj = zti(1,jj+1,ji,A)
                  zx(1,jj,ji,A) = big_fact_zaa(1,jj,ji)*(junki - junk)/(junk + junki + rtrn)* rsc
                  zy(1,jj,ji,A) = big_fact_zbb(1,jj,ji)*(junkj - junk)/(junk + junkj + rtrn)* rsc
            END DO
      END DO

      !DO ju=1, dimen_jarr2
         DO ji = 2,jpim1
            DO jj = 2,jpjm1
            !dir$ vector aligned
              DO jk = 2,jpkm1

            !    ji = jarr2(3, ju)
            !    jj = jarr2(2, ju)
            !    jk = jarr2(1, ju)
               !junk  = zti(jk,jj,ji,A)
               !junki = zti(jk,jj,ji+1,A)
               !junkj = zti(jk,jj+ 1,ji,A)
               !junkk = zti(jk-1,jj,ji,A)
               zx(jk,jj,ji,A) = advmask(jk,jj,ji)*(big_fact_zaa(jk,jj,ji)*(zti(jk,jj,ji+1,A) - zti(jk,jj,ji,A))/(zti(jk,jj,ji,A) + zti(jk,jj,ji+1,A) + rtrn)* rsc)
               zy(jk,jj,ji,A) = advmask(jk,jj,ji)*(big_fact_zbb(jk,jj,ji)*(zti(jk,jj+ 1,ji,A) - zti(jk,jj,ji,A))/(zti(jk,jj,ji,A) + zti(jk,jj+ 1,ji,A) + rtrn)* rsc)
               zz(jk,jj,ji,A) = advmask(jk,jj,ji)*(big_fact_zcc(jk,jj,ji)*(zti(jk,jj,ji,A) -  zti(jk-1,jj,ji,A))/(zti(jk,jj,ji,A) +  zti(jk-1,jj,ji,A) + rtrn)* rsc*( -1.))

           END DO
           END DO
           END DO
!                 endif

!! epascolo mpi comment 
!! ... Lateral boundary conditions on z[xyz]
! #ifdef key_mpp

! !!   ... Mpp : export boundary values to neighboring processors


!       IF( ntids - 1 + jn <= jptra ) THEN
!        pack_size = ntids
!         ELSE
!            pack_size = ntids - (ntids - 1 + jn - jptra)
!       END IF

!         CALL mpplnk_my(zx, pack_size,1,1)
!         CALL mpplnk_my(zy, pack_size,1,1)
!         CALL mpplnk_my(zz, pack_size,1,1)


! #else

! !!   ... T-point, 3D array, full local array z[xyz] are initialised


!           DO itid = 1, ntids
!          IF( itid - 1 + jn <= jptra ) THEN
                CALL lbc( zx(:,:,:,1), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zy(:,:,:,1), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zz(:,:,:,1), 1, 1, 1, 1, jpk, 1 )
!          END IF
!       END DO

! #endif

!! 2.4 reinitialization

!!!&omp     parallel default(none) private(mytid,A,junk,jk,jj,ji)
!!!&omp&      shared(zkx,zky,zkz,zz,zx,zy,zti,jpjm1,jpim1,dimen_jarr1,jarr1,jpi,jpk,jpj,jn)

      
      A=1
!!            2.5 calcul of the final field:
!!                advection by antidiffusive mass fluxes and an upstream scheme
           
                 DO ji = 2,jpim1
                 !dir$ vector aligned
             DO jj = 2,jpjm1
                   junk  = zti(1,jj,ji,A)
                   zkx(1,jj,ji,A) = fsx(junk,zti(1,jj,ji+1,A),zx(1,jj,ji,A))
                 END DO
             END DO

                 DO ji = 2,jpim1
                 !dir$ vector aligned
             DO jj = 2,jpjm1
                   junk  = zti(1,jj,ji,A)
                   zky(1,jj,ji,A) = fsy(junk,zti(1,jj+ 1,ji,A),zy(1,jj,ji,A))
                 END DO
              END DO

            DO ji = 1,jpi 
            !dir$ vector aligned
                   DO jk = 2,jpk   
                        zkz(jk,1,ji,A) = fsz(zti(jk,1,ji,A),zti(jk-1,1,ji,A),zz(jk,1,ji,A))
                  ENDDO
            ENDDO
               
            DO ji = 1,jpi
            !dir$ vector aligned
                  DO jk = 2,jpk 
                        zkz(jk,jpj,ji,A) = fsz(zti(jk,jpj,ji,A),zti(jk-1,jpj,ji,A),zz(jk,jpj,ji,A))
                  ENDDO
            ENDDO

             DO jj = 2,jpjm1
             !dir$ vector aligned
                  DO jk = 2,jpk
                        zkz(jk,jj,1,A) = fsz(zti(jk,jj,1,A),zti(jk-1,jj,1,A),zz(jk,jj,1,A))
                  ENDDO
            ENDDO   

             DO jj = 2,jpjm1
             !dir$ vector aligned
                  DO jk = 2,jpk
                        zkz(jk,jj,jpi,A) = fsz(zti(jk,jj,jpi,A),zti(jk-1,jj,jpi,A),zz(jk,jj,jpi,A))
                  END DO
             END DO

               DO  ji = 2,jpim1
            DO jj = 2,jpjm1
            !dir$ vector aligned
         DO jk = 2,jpk    
                zkx(jk,jj,ji,A) = fsx(zti(jk,jj,ji,A),zti(jk,jj,ji + 1,A),zx(jk,jj,ji,A))*advmask(jk,jj,ji)
                !zky(jk,jj,ji,A) = fsy(zti(jk,jj,ji,A),zti(jk,jj+ 1,ji,A),zy(jk,jj,ji,A))*advmask(jk,jj,ji)
                !zkz(jk,jj,ji,A) = fsz(zti(jk,jj,ji,A),zti(jk-1,jj,ji,A),zz(jk,jj,ji,A))*advmask(jk,jj,ji)
         END DO
          END DO
           END DO

        DO  ji = 2,jpim1
            DO jj = 2,jpjm1
            !dir$ vector aligned
         DO jk = 2,jpk    
                !zkx(jk,jj,ji,A) = fsx(zti(jk,jj,ji,A),zti(jk,jj,ji + 1,A),zx(jk,jj,ji,A))*advmask(jk,jj,ji)
                zky(jk,jj,ji,A) = fsy(zti(jk,jj,ji,A),zti(jk,jj+ 1,ji,A),zy(jk,jj,ji,A))*advmask(jk,jj,ji)
                !zkz(jk,jj,ji,A) = fsz(zti(jk,jj,ji,A),zti(jk-1,jj,ji,A),zz(jk,jj,ji,A))*advmask(jk,jj,ji)
         END DO
          END DO
           END DO

        DO  ji = 2,jpim1
            DO jj = 2,jpjm1
            !dir$ vector aligned
         DO jk = 2,jpk    
                !zkx(jk,jj,ji,A) = fsx(zti(jk,jj,ji,A),zti(jk,jj,ji + 1,A),zx(jk,jj,ji,A))*advmask(jk,jj,ji)
                !zky(jk,jj,ji,A) = fsy(zti(jk,jj,ji,A),zti(jk,jj+ 1,ji,A),zy(jk,jj,ji,A))*advmask(jk,jj,ji)
                zkz(jk,jj,ji,A) = fsz(zti(jk,jj,ji,A),zti(jk-1,jj,ji,A),zz(jk,jj,ji,A))*advmask(jk,jj,ji)
         END DO
          END DO
           END DO

!! ... Lateral boundary conditions on zk[xy]
! epascolo mpi comment
! #ifdef key_mpp
! !!   ... Mpp : export boundary values to neighboring processors


!         IF( ntids - 1 + jn <= jptra ) THEN
!            pack_size = ntids
!         ELSE
!            pack_size = ntids - (ntids - 1 + jn - jptra)
!         END IF

!         CALL mpplnk_my(zkx, pack_size,1,1)
!         CALL mpplnk_my(zky, pack_size,1,1)
! #else
! !!   ... T-point, 3D array, full local array zk[xy] are initialised
!       DO itid = 1, ntids
!         IF( itid - 1 + jn <= jptra ) THEN
               CALL lbc( zkx(:,:,:,1), 1, 1, 1, 1, jpk, 1 )
               CALL lbc( zky(:,:,:,1), 1, 1, 1, 1, jpk, 1 )
!         END IF
!       END DO
! #endif

!!!&omp    parallel default(none) private(mytid,A,B,zbtr,jk,jj,ji,ju,jf)
!!!&omp&      shared(zkx,zky,zkz,zbtr_arr,e1t,e2t,ztj,dimen_jarr3,jarr3,ncor,jn,
!!!&omp&             jarr_adv_flx,Fsize,diaflx)


!!        2.6. calcul of after field using an upstream advection scheme

         A = 1
         B = jn
            if(ncor .EQ. 1) then
               DO ju=1, dimen_jarr3

                  ji = jarr3(3, ju)
                  jj = jarr3(2, ju)
                  jk = jarr3(1, ju)
                  jf = jarr_adv_flx(ju)

                  zbtr = zbtr_arr(jk,jj,ji)
                  ztj(jk,jj,ji,A) = -zbtr*( zkx(jk,jj,ji,A) - zkx(jk,jj,ji - 1,A) &
     &              + zky(jk,jj,ji,A) - zky(jk,jj- 1,ji,A)+ zkz(jk,jj,ji,A) - zkz(jk+1,jj,ji,A) )+ ztj(jk,jj,ji,A)

!     Save advective fluxes x,y,z
              IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                 diaflx(jf,B,3) = diaflx(jf,B,3) + zkx(jk,jj,ji,A)
                 diaflx(jf,B,2) = diaflx(jf,B,2) + zky(jk,jj,ji,A)
                 diaflx(jf,B,1) = diaflx(jf,B,1) + zkz(jk,jj,ji,A)
              END IF

              END DO
           else
              DO ju=1, dimen_jarr3

                  ji = jarr3(3, ju)
                  jj = jarr3(2, ju)
                  jk = jarr3(1, ju)
                  jf = jarr_adv_flx(ju)

                  zbtr = zbtr_arr(jk,jj,ji)
                  ztj(jk,jj,ji,A) = -zbtr*( zkx(jk,jj,ji,A) - zkx(jk,jj,ji - 1,A) &
     &              + zky(jk,jj,ji,A) - zky(jk,jj- 1,ji,A)+ zkz(jk,jj,ji,A) - zkz(jk+1,jj,ji,A) )

                 !Save advective fluxes x,y,z
                 IF ( (Fsize .GT. 0) .AND. ( jf .GT. 0 ) ) THEN
                    diaflx(jf,B,3) = diaflx(jf,B,3) + zkx(jk,jj,ji,A)
                    diaflx(jf,B,2) = diaflx(jf,B,2) + zky(jk,jj,ji,A)
                    diaflx(jf,B,1) = diaflx(jf,B,1) + zkz(jk,jj,ji,A)
                 END IF
              END DO

           endif

        ENDDO ANTIDIFF_CORR


!!       3. trend due to horizontal and vertical advection of tracer jn
!!!&omp   parallel default(none) private(mytid,A,B,jk,jj,ji,ju) shared(ncor,dimen_jarrt,jarrt,tra,ztj,jn,zbuf)

        A=1
        B=jn

           if(ncor .EQ. 1) then
              DO ju=1, dimen_jarrt
                 ji = jarrt(3, ju)
                 jj = jarrt(2, ju)
                 jk = jarrt(1, ju)

                 tra(jk,jj,ji,B) = tra(jk,jj,ji,B)+ ztj(jk,jj,ji,A)

              END DO
           else
              DO ju=1, dimen_jarrt
                 ji = jarrt(3, ju)
                 jj = jarrt(2, ju)
                 jk = jarrt(1, ju)

                 tra(jk,jj,ji,B) = tra(jk,jj,ji,B)+ (zbuf(jk,jj,ji,A) + ztj(jk,jj,ji,A))

              END DO
           endif

       END DO TRACER_LOOP

       trcadvparttime = MPI_WTIME() - trcadvparttime
       trcadvtottime = trcadvtottime + trcadvparttime


      END SUBROUTINE
