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
      double precision :: zbtr,zdt
      double precision :: junk, junki, junkj, junkk
      double precision :: timer
      
#ifndef ExecEnsNO
      double precision,dimension(:), allocatable :: array 
      double precision,dimension(:,:), allocatable :: surface 
      double precision, allocatable,dimension(:,:,:) :: zti,ztj
      double precision, allocatable,dimension(:,:,:) :: zx,zy,zz,zbuf
      double precision, allocatable,dimension(:,:,:) :: zkx,zky,zkz
#endif

!-------------------------------------------------------------------

      MPI_CHECK = .FALSE.

      if(.not.adv_initialized ) then  ! INIT phase

        
        

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

         !write(*,*) "Storing good points ..."

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

         !write(*,*) 'trcadv: RANK -> ', myrank, ' all_points -> ', allpoints

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

         !write(*,*) 'trcadv: RANK -> ', myrank, ' good_points -> ', goodpoints



      adv_initialized=.true.
      endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end initialization phase

         jk=1

         zdt = rdt*ndttrc
         !$OMP TASK private(ji,jj) firstprivate(jpim1,jpjm1) shared(zbtr_arr,e1t,e2t,e3t) default(none)

         DO ji = 1,jpi
         DO jj = 1,jpj
            !dir$ vector aligned
         DO jk = 1,jpkm1
                  zbtr_arr(jk,jj,ji) = 1./(e1t(jj,ji)*e2t(jj,ji)*e3t(jk,jj,ji))
         END DO
         END DO
         END DO
         !$OMP END TASK
        
          !$OMP TASK private(ji,jj) firstprivate(jpim1,jpjm1,jpi,jpj,jpk) default(none) &
          !$OMP shared(zdt,zaa,inv_eu,e1u,e2u,e3u,un,big_fact_zaa)

         DO ji = 1,jpi
         DO jj = 1,jpj
            !dir$ vector aligned
         DO jk = 1,jpkm1
                  inv_eu(jk,jj,ji) = 1./(e1u(jj,ji)*e2u(jj,ji)*e3u(jk,jj,ji) )
         END DO
         END DO
         END DO




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
             big_fact_zaa(jk,jj,ji) = ( abs(zaa(jk,jj,ji)) - zdt*zaa(jk,jj,ji)**2*inv_eu(jk,jj,ji) )
             !/(e1u(jj,ji)*e2u(jj,ji)*e3t(jk,jj,ji) ) )
            END DO
            END DO
            END DO
      
          !$OMP END TASK
           
          !$OMP TASK private(ji,jj) firstprivate(jpim1,jpjm1,jpi,jpj,jpk)  default(none) &
          !$OMP shared(inv_ev,e1v,e2v,e3v,vn,zdt,zbb,big_fact_zbb)

         DO ji = 1,jpi
         DO jj = 1,jpj
            !dir$ vector aligned
         DO jk = 1,jpkm1
                  inv_ev(jk,jj,ji) = 1./(e1v(jj,ji)*e2v(jj,ji)*e3v(jk,jj,ji) )
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
            !$OMP END TASK
             
            !$OMP TASK private(ji,jj) firstprivate(jpim1,jpjm1,jpi,jpj,jpk) default(none) &
            !$OMP shared(inv_et,e1t,e2t,e3w,wn,zcc,zdt,big_fact_zcc)   

         DO ji = 1,jpi
         DO jj = 1,jpj
            !dir$ vector aligned
         DO jk = 1,jpkm1
                  inv_et(jk,jj,ji) = 1./(e1t(jj,ji)*e2t(jj,ji)*e3w(jk,jj,ji) )
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

            !$OMP END TASK
      
      !$OMP TASKWAIT

     
!!     tracer loop parallelized (macrotasking)
!!     =======================================

      trcadvparttime = MPI_WTIME()
       
!$omp taskloop default(none) private(jf,junk,junki,junkj,junkk,zbtr) &
!$omp private(zkx,zky,zkz,zti,ztj,zx,zy,zz,zbuf) shared(diaflx,jarrt,tra,zdt) &
!$omp shared(big_fact_zaa,big_fact_zbb,big_fact_zcc,zaa,zbb,zcc,inv_eu,inv_ev,inv_et) &
!$omp shared(jpim1,jpjm1,un,vn,wn,e2u,e3u,e3v,e1v,e1t,e2t,e3t,trn,advmask,jarr3,jarr_adv_flx,zbtr_arr) &
!$omp firstprivate(jpkm1,dimen_jarr3,Fsize,ncor,rtrn,rsc,dimen_jarrt,jpj,jpi,jpk) 

       
      TRACER_LOOP: DO  jn = 1, jptra


!!        1. tracer flux in the 3 directions
!!        ----------------------------------
!!        1.1 mass flux at u v and t-points and initialization
!!       1.2 calcul of intermediate field with an upstream advection scheme
!!           and mass fluxes calculated above
!!       calcul of tracer flux in the i and j direction
       
#ifndef ExecEnsNO
       allocate(zy(jpk,jpj,jpi))  
       allocate(zx(jpk,jpj,jpi))
       allocate(zz(jpk,jpj,jpi))
       allocate(ztj(jpk,jpj,jpi)) 
       allocate(zti(jpk,jpj,jpi))    
       allocate(zkx(jpk,jpj,jpi)) 
       allocate(zky(jpk,jpj,jpi)) 
       allocate(zkz(jpk,jpj,jpi)) 
       allocate(zbuf(jpk,jpj,jpi))
#endif

       zy(:,:,:) = 0
       zz(:,:,:) = 0 
       zx(:,:,:) = 0
       ztj(:,:,:)= 0
       zti(:,:,:)= 0
       zbuf(:,:,:) = 0.
       zkx(:,:,:)=0.  
       zky(:,:,:)=0.  
       zkz(:,:,:)=0.

!        zkx(  :,:,1)=0.  
!        zkx(:,:,jpi)=0.    
!        zky(:,  1,:)=0.  
!        zky(:,jpj,:)=0.  
!        zkz(1,:,:)  =0.
! ! loop unfusion

        DO ji = 2,jpim1
           !dir$ vector aligned
           DO jj = 2,jpjm1
                 zkx(1,jj,ji ) = fsx(trn(1,jj,ji, jn),trn(1,jj,ji + 1, jn),zaa(1,jj,ji))
           END DO
        END DO
        

        DO ji = 2,jpim1
          !dir$ vector aligned
           DO jj = 2,jpjm1
                zky(1,jj,ji ) = fsy(trn(1,jj,ji, jn),trn(1,jj+1,ji, jn),zbb(1,jj,ji))
           END DO
        END DO

       DO ji = 1,jpi
           !dir$ vector aligned
           DO jk = 2,jpk
            zkz(jk,1,ji ) = fsz(trn(jk,1,ji, jn),trn(jk-1,1,ji, jn),zcc(jk,1,ji))
            ENDDO
       ENDDO
      
       DO ji = 1,jpi
            !dir$ vector aligned
            DO jk = 2,jpk
            zkz(jk,jpj,ji ) = fsz(trn(jk,jpj,ji, jn),trn(jk-1,jpj,ji, jn),zcc(jk,jpj,ji))
            END DO
       END DO
! loop unfusion
      DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
            zkz(jk,jj,1 ) = fsz(trn(jk,jj,1, jn),trn(jk-1,jj,1, jn),zcc(jk,jj,1))
            END DO
      END DO

      DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
            zkz(jk,jj,jpi ) = fsz(trn(jk,jj,jpi, jn),trn(jk-1,jj,jpi, jn),zcc(jk,jj,jpi))
            END DO
      END DO

      DO  ji = 2,jpim1
        DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
              zkx(jk,jj,ji ) = fsx(trn(jk,jj,ji, jn),trn(jk,jj,ji + 1, jn),zaa(jk,jj,ji))*advmask(jk,jj,ji)
            END DO
         END DO
      END DO

      DO  ji = 2,jpim1
        DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
              zky(jk,jj,ji ) = fsy(trn(jk,jj,ji, jn),trn(jk,jj + 1,ji, jn),zbb(jk,jj,ji))*advmask(jk,jj,ji)
              END DO
         END DO
      END DO

            DO  ji = 2,jpim1
        DO jj = 2,jpjm1
            !dir$ vector aligned
            DO jk = 2,jpk
              zkz(jk,jj,ji ) = fsz(trn(jk,jj,ji, jn),trn(jk-1,jj,ji, jn),zcc(jk,jj,ji))*advmask(jk,jj,ji)
            END DO
         END DO
      END DO


! ... Lateral boundary conditions on zk[xy]
#ifdef key_mpp

!  ... Mpp : export boundary values to neighboring processors

         CALL mpplnk_my(zkx)
         CALL mpplnk_my(zky)

#else

! !!   ... T-point, 3D array, full local arrays zk[xy] are initialised

               CALL lbc( zkx(:,:,:), 1, 1, 1, 1, jpk, 1 )
               CALL lbc( zky(:,:,:), 1, 1, 1, 1, jpk, 1 )
#endif


!! 2. calcul of after field using an upstream advection scheme
!! -----------------------------------------------------------
          DO ji =2,jpim1
          DO jj =2,jpjm1
          DO jk =1,jpkm1
              ztj(jk,jj,ji ) = -zbtr_arr(jk,jj,ji)* &
     &          ( zkx(jk,jj,ji ) - zkx(jk,jj,ji-1 ) &
     &          + zky(jk,jj,ji ) - zky(jk,jj- 1,ji ) &
     &          + zkz(jk,jj,ji ) - zkz(jk+1,jj,ji ) )
          ENDDO
          ENDDO
          ENDDO

          DO jf=1,Fsize
             jk = flx_ridxt(jf,2)
             jj = flx_ridxt(jf,3)
             ji = flx_ridxt(jf,4)

             diaflx(1,jf, jn) = diaflx(1,jf, jn) + zkx(jk,jj,ji )*rdt
             diaflx(2,jf, jn) = diaflx(2,jf, jn) + zky(jk,jj,ji )*rdt
             diaflx(3,jf, jn) = diaflx(3,jf, jn) + zkz(jk,jj,ji )*rdt
          ENDDO



!! 2.1 start of antidiffusive correction loop

        ANTIDIFF_CORR: DO jt = 1,ncor

!! 2.2 calcul of intermediary field zti



       
        
             if(jt .EQ. 1) then

                if(ncor .EQ. 1) then

                         DO ji = 2,jpim1
                      DO jj = 2,jpjm1
                      !dir$ vector aligned
                   DO jk = 1,jpkm1
                            zti(jk,jj,ji ) = trn(jk,jj,ji, jn) + zdt*ztj(jk,jj,ji )
                        END DO
                      END DO
                   END DO

                else
                         DO ji = 2,jpim1
                      DO jj = 2,jpjm1
                      !dir$ vector aligned
                    DO jk = 1,jpkm1
                            zti(jk,jj,ji ) = trn(jk,jj,ji, jn) + zdt*ztj(jk,jj,ji )
                          !  zbuf(jk,jj,ji ) = ztj(jk,jj,ji )
                         END DO
                      END DO
                   END DO

                        DO ji = 2,jpim1
                      DO jj = 2,jpjm1
                      !dir$ vector aligned
                    DO jk = 1,jpkm1
                         !   zti(jk,jj,ji ) = trn(jk,jj,ji, jn) + zdt*ztj(jk,jj,ji )
                            zbuf(jk,jj,ji ) = ztj(jk,jj,ji )
                         END DO
                      END DO
                   END DO
                
                
                endif

             else

                      DO ji = 2,jpim1
                   DO jj = 2,jpjm1
                   !dir$ vector aligned
                DO jk = 1,jpkm1
                         zti(jk,jj,ji ) =  zti(jk,jj,ji ) + zdt*ztj(jk,jj,ji )
                      END DO
                   END DO
                END DO

                      DO ji = 2,jpim1
                   DO jj = 2,jpjm1
                   !dir$ vector aligned
                DO jk = 1,jpkm1
                         zbuf(jk,jj,ji ) = zbuf(jk,jj,ji ) + ztj(jk,jj,ji )
                      END DO
                   END DO
                END DO

             endif


!! ... Lateral boundary conditions on zti
#ifdef key_mpp
! ... Mpp : export boundary values to neighboring processors
         CALL mpplnk_my(zti)
#else
! ... T-point, 3D array, full local array zti is initialised
                 CALL lbc( zti(:,:,:), 1, 1, 1, 1, jpk, 1 )
#endif


!! 2.3 calcul of the antidiffusive flux
      
       
      !jk = 1
!          DO jk = 1,jpkm1
      DO ji = 2,jpim1
            DO jj = 2,jpjm1
                  junk  = zti(1,jj,ji )
                  junki = zti(1,jj,ji+1 )
                  junkj = zti(1,jj+1,ji )
                  zx(1,jj,ji ) = big_fact_zaa(1,jj,ji)*(junki - junk)/(junk + junki + rtrn)* rsc
                  zy(1,jj,ji ) = big_fact_zbb(1,jj,ji)*(junkj - junk)/(junk + junkj + rtrn)* rsc
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
               !junk  = zti(jk,jj,ji )
               !junki = zti(jk,jj,ji+1 )
               !junkj = zti(jk,jj+ 1,ji )
               !junkk = zti(jk-1,jj,ji )
               zx(jk,jj,ji) = advmask(jk,jj,ji)*(big_fact_zaa(jk,jj,ji)*(zti(jk,jj,ji+1) - zti(jk,jj,ji))/(zti(jk,jj,ji) + &
                   zti(jk,jj,ji+1) + rtrn)*rsc)
               zy(jk,jj,ji) = advmask(jk,jj,ji)*(big_fact_zbb(jk,jj,ji)*(zti(jk,jj+1,ji) - zti(jk,jj,ji))/(zti(jk,jj,ji) + &
                   zti(jk,jj+1,ji) + rtrn)*rsc)
               zz(jk,jj,ji) = advmask(jk,jj,ji)*(big_fact_zcc(jk,jj,ji)*(zti(jk,jj,ji) - zti(jk-1,jj,ji))/(zti(jk,jj,ji) + &
                   zti(jk-1,jj,ji) + rtrn)*rsc*(-1.))

           END DO
           END DO
           END DO
!                 endif


! ... Lateral boundary conditions on z[xyz]
#ifdef key_mpp

! ... Mpp : export boundary values to neighboring processors

         CALL mpplnk_my(zx)
         CALL mpplnk_my(zy)
         CALL mpplnk_my(zz)

#else

!  ... T-point, 3D array, full local array z[xyz] are initialised
                CALL lbc( zx(:,:,:), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zy(:,:,:), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zz(:,:,:), 1, 1, 1, 1, jpk, 1 )
#endif

!! 2.4 reinitialization
!!            2.5 calcul of the final field:
!!                advection by antidiffusive mass fluxes and an upstream scheme
           
                 DO ji = 2,jpim1
                 !dir$ vector aligned
             DO jj = 2,jpjm1
                   zkx(1,jj,ji ) = fsx(zti(1,jj,ji ),zti(1,jj,ji+1 ),zx(1,jj,ji ))
                 END DO
             END DO

                 DO ji = 2,jpim1
                 !dir$ vector aligned
             DO jj = 2,jpjm1 
                   zky(1,jj,ji ) = fsy(zti(1,jj,ji ),zti(1,jj+ 1,ji ),zy(1,jj,ji ))
                 END DO
              END DO

            DO ji = 1,jpi 
            !dir$ vector aligned
                   DO jk = 2,jpk   
                        zkz(jk,1,ji ) = fsz(zti(jk,1,ji ),zti(jk-1,1,ji ),zz(jk,1,ji ))
                  ENDDO
            ENDDO
               
            DO ji = 1,jpi
            !dir$ vector aligned
                  DO jk = 2,jpk 
                        zkz(jk,jpj,ji ) = fsz(zti(jk,jpj,ji ),zti(jk-1,jpj,ji ),zz(jk,jpj,ji ))
                  ENDDO
            ENDDO

             DO jj = 2,jpjm1
             !dir$ vector aligned
                  DO jk = 2,jpk
                        zkz(jk,jj,1 ) = fsz(zti(jk,jj,1 ),zti(jk-1,jj,1 ),zz(jk,jj,1 ))
                  ENDDO
            ENDDO   

             DO jj = 2,jpjm1
             !dir$ vector aligned
                  DO jk = 2,jpk
                        zkz(jk,jj,jpi ) = fsz(zti(jk,jj,jpi ),zti(jk-1,jj,jpi ),zz(jk,jj,jpi ))
                  END DO
             END DO

               DO  ji = 2,jpim1
            DO jj = 2,jpjm1
            !dir$ vector aligned
         DO jk = 2,jpk    
                zkx(jk,jj,ji ) = fsx(zti(jk,jj,ji ),zti(jk,jj,ji + 1 ),zx(jk,jj,ji ))*advmask(jk,jj,ji)
                !zky(jk,jj,ji ) = fsy(zti(jk,jj,ji ),zti(jk,jj+ 1,ji ),zy(jk,jj,ji ))*advmask(jk,jj,ji)
                !zkz(jk,jj,ji ) = fsz(zti(jk,jj,ji ),zti(jk-1,jj,ji ),zz(jk,jj,ji ))*advmask(jk,jj,ji)
         END DO
          END DO
           END DO

        DO  ji = 2,jpim1
            DO jj = 2,jpjm1
            !dir$ vector aligned
         DO jk = 2,jpk    
                !zkx(jk,jj,ji ) = fsx(zti(jk,jj,ji ),zti(jk,jj,ji + 1 ),zx(jk,jj,ji ))*advmask(jk,jj,ji)
                zky(jk,jj,ji ) = fsy(zti(jk,jj,ji ),zti(jk,jj+ 1,ji ),zy(jk,jj,ji ))*advmask(jk,jj,ji)
                !zkz(jk,jj,ji ) = fsz(zti(jk,jj,ji ),zti(jk-1,jj,ji ),zz(jk,jj,ji ))*advmask(jk,jj,ji)
         END DO
          END DO
           END DO

        DO  ji = 2,jpim1
            DO jj = 2,jpjm1
            !dir$ vector aligned
         DO jk = 2,jpk    
                !zkx(jk,jj,ji ) = fsx(zti(jk,jj,ji ),zti(jk,jj,ji + 1 ),zx(jk,jj,ji ))*advmask(jk,jj,ji)
                !zky(jk,jj,ji ) = fsy(zti(jk,jj,ji ),zti(jk,jj+ 1,ji ),zy(jk,jj,ji ))*advmask(jk,jj,ji)
                zkz(jk,jj,ji ) = fsz(zti(jk,jj,ji ),zti(jk-1,jj,ji ),zz(jk,jj,ji ))*advmask(jk,jj,ji)
         END DO
          END DO
           END DO

!... Lateral boundary conditions on zk[xy]
#ifdef key_mpp
!  ... Mpp : export boundary values to neighboring processors

         CALL mpplnk_my(zkx)
         CALL mpplnk_my(zky)
#else
! ... T-point, 3D array, full local array zk[xy] are initialised
               CALL lbc( zkx(:,:,:), 1, 1, 1, 1, jpk, 1 )
               CALL lbc( zky(:,:,:), 1, 1, 1, 1, jpk, 1 )
#endif

!!        2.6. calcul of after field using an upstream advection scheme

          
           
            if(ncor .EQ. 1) then
          DO ji =2,jpim1
          DO jj =2,jpjm1
          DO jk =1,jpkm1
                  ztj(jk,jj,ji ) = -zbtr_arr(jk,jj,ji)*( zkx(jk,jj,ji ) - zkx(jk,jj,ji - 1 ) &
     &              + zky(jk,jj,ji ) - zky(jk,jj- 1,ji )+ zkz(jk,jj,ji ) - zkz(jk+1,jj,ji ) )+ ztj(jk,jj,ji )
          ENDDO
          ENDDO
          ENDDO

         DO jf=1,Fsize
             jk = flx_ridxt(jf,2)
             jj = flx_ridxt(jf,3)
             ji = flx_ridxt(jf,4)

             diaflx(1,jf, jn) = diaflx(1,jf, jn) + zkx(jk,jj,ji )*rdt
             diaflx(2,jf, jn) = diaflx(2,jf, jn) + zky(jk,jj,ji )*rdt
             diaflx(3,jf, jn) = diaflx(3,jf, jn) + zkz(jk,jj,ji )*rdt
          ENDDO



           else
         DO ji =2,jpim1
          DO jj =2,jpjm1
          DO jk =1,jpkm1
                  ztj(jk,jj,ji ) = -zbtr_arr(jk,jj,ji)*( zkx(jk,jj,ji ) - zkx(jk,jj,ji - 1 ) &
     &              + zky(jk,jj,ji ) - zky(jk,jj- 1,ji )+ zkz(jk,jj,ji ) - zkz(jk+1,jj,ji ) )
          ENDDO
          ENDDO
          ENDDO

         DO jf=1,Fsize
             jk = flx_ridxt(jf,2)
             jj = flx_ridxt(jf,3)
             ji = flx_ridxt(jf,4)

             diaflx(1,jf, jn) = diaflx(1,jf, jn) + zkx(jk,jj,ji )*rdt
             diaflx(2,jf, jn) = diaflx(2,jf, jn) + zky(jk,jj,ji )*rdt
             diaflx(3,jf, jn) = diaflx(3,jf, jn) + zkz(jk,jj,ji )*rdt
          ENDDO






           endif

        ENDDO ANTIDIFF_CORR


!!       3. trend due to horizontal and vertical advection of tracer jn
         
          

           if(ncor .EQ. 1) then
           do ji=1,jpi
           do jj=1,jpj
           do jk=1,jpk
               tra(jk,jj,ji, jn) = tra(jk,jj,ji, jn)+ ztj(jk,jj,ji )
           enddo
           enddo
           enddo


           else
           do ji=1,jpi
           do jj=1,jpj
           do jk=1,jpk
               tra(jk,jj,ji, jn) = tra(jk,jj,ji, jn)+ (zbuf(jk,jj,ji ) + ztj(jk,jj,ji ))
           enddo
           enddo
           enddo


           endif

#ifndef ExecEnsNO
        deallocate(zy )  
        deallocate(zx )
        deallocate(zz )
        deallocate(ztj ) 
        deallocate(zti )    
        deallocate(zkx ) 
        deallocate(zky ) 
        deallocate(zkz ) 
        deallocate(zbuf )
#endif


       END DO TRACER_LOOP
      !$OMP end taskloop 

        
       trcadvparttime = MPI_WTIME() - trcadvparttime
       trcadvtottime = trcadvtottime + trcadvparttime
!!!!

      contains

         double precision FUNCTION fsx( pfx1, pfx2, pfu )
         !dir$ attributes vector :: fsx
       IMPLICIT NONE
            double precision, INTENT(IN) :: pfx1, pfx2, pfu
            double precision ::  abspfu
            abspfu = abs(pfu)
            fsx = ( ( pfu + abspfu ) * pfx1+( pfu - abspfu ) * pfx2 ) * 0.5
       END FUNCTION fsx


       double precision FUNCTION fsy( pfy1, pfy2, pfv  )
       !dir$ attributes vector :: fsy
       IMPLICIT NONE
            double precision, INTENT(IN) :: pfy1, pfy2, pfv
            double precision :: abspfv
            abspfv = abs(pfv)
            fsy = ( ( pfv + abspfv ) * pfy1 +( pfv - abspfv ) * pfy2 ) * 0.5
       END FUNCTION fsy


       double precision FUNCTION fsz( pfz1, pfz2, pfw )
       !dir$ attributes vector :: fsz
       IMPLICIT NONE
       double precision, INTENT(IN) :: pfz1, pfz2, pfw
       double precision abspfw
            abspfw = abs(pfw)
            fsz = ( ( pfw + abspfw ) * pfz1+( pfw - abspfw ) * pfz2 ) * 0.5
       END FUNCTION fsz       



      END SUBROUTINE
