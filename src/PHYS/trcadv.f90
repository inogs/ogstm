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
  USE mpi

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
  double precision,dimension(:), allocatable :: array
  double precision,dimension(:,:), allocatable :: surface
  double precision, allocatable,dimension(:,:,:) :: zti,ztj
  double precision, allocatable,dimension(:,:,:) :: zx,zy,zz,zbuf
  double precision, allocatable,dimension(:,:,:) :: zkx,zky,zkz
  real(kind=8) :: t1,t2
  real(kind=8),save :: comm_time=0,mpi_time=0,kernel_time=0,init_time=0,cpu_time=0

  !-------------------------------------------------------------------

  MPI_CHECK = .FALSE.

  if(.not.adv_initialized ) then  ! INIT phase




     t1=MPI_WTime()
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

     write(*,*) "Storing good points ..."

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



     adv_initialized=.true.
     t2=MPI_Wtime()
     cpu_time=cpu_time+t2-t1
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! end initialization phase

  jk=1

  zdt = rdt*ndttrc
  !$OMP TASK private(ji,jj) firstprivate(jpim1,jpjm1) shared(zbtr_arr,e1t,e2t,e3t) default(none)


!!!         !$acc enter data create(...)


  t1=MPI_Wtime()
  !$acc enter data create( zaa(1:jpk,1:jpj,1:jpi), zbb(1:jpk,1:jpj,1:jpi), zcc(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( inv_eu(1:jpk,1:jpj,1:jpi), inv_ev(1:jpk,1:jpj,1:jpi), inv_et(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( big_fact_zaa (1:jpk,1:jpj,1:jpi), big_fact_zbb(1:jpk,1:jpj,1:jpi), big_fact_zcc(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( zbtr_arr(1:jpk,1:jpj,1:jpi) )
  t2=MPI_Wtime()
  init_time=init_time+t2-t1


  t1=MPI_Wtime()
  !$acc enter data create( e1t(1:jpj,1:jpi), e2t(1:jpj,1:jpi), e3t(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( e1u(1:jpj,1:jpi), e2u(1:jpj,1:jpi), e3u(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( e1v(1:jpj,1:jpi), e2v(1:jpj,1:jpi), e3v(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( e3w(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( un(1:jpk,1:jpj,1:jpi), vn(1:jpk,1:jpj,1:jpi), wn(1:jpk,1:jpj,1:jpi) )


!!!!!$acc update device

  !$acc update device( zaa(1:jpk,1:jpj,1:jpi), zbb(1:jpk,1:jpj,1:jpi), zcc(1:jpk,1:jpj,1:jpi) )
  !$acc update device( inv_eu(1:jpk,1:jpj,1:jpi), inv_ev(1:jpk,1:jpj,1:jpi), inv_et(1:jpk,1:jpj,1:jpi) )
  !$acc update device( big_fact_zaa (1:jpk,1:jpj,1:jpi), big_fact_zbb(1:jpk,1:jpj,1:jpi), big_fact_zcc(1:jpk,1:jpj,1:jpi) )
  !$acc update device( zbtr_arr(1:jpk,1:jpj,1:jpi) )

  !$acc update device( e1t(1:jpj,1:jpi), e2t(1:jpj,1:jpi), e3t(1:jpk,1:jpj,1:jpi) )
  !$acc update device( e1u(1:jpj,1:jpi), e2u(1:jpj,1:jpi), e3u(1:jpk,1:jpj,1:jpi) )
  !$acc update device( e1v(1:jpj,1:jpi), e2v(1:jpj,1:jpi), e3v(1:jpk,1:jpj,1:jpi) )
  !$acc update device( e3w(1:jpk,1:jpj,1:jpi) )
  !$acc update device( un(1:jpk,1:jpj,1:jpi), vn(1:jpk,1:jpj,1:jpi), wn(1:jpk,1:jpj,1:jpi) )
  t2=MPI_Wtime()
  comm_time=comm_time+t2-t1



!!$       zaa(jpk,jpj,jpi)
!!$       zbb(jpk,jpj,jpi)
!!$       zcc(jpk,jpj,jpi)
!!$       inv_eu(jpk,jpj,jpi)
!!$       inv_ev(jpk,jpj,jpi)
!!$       inv_et(jpk,jpj,jpi)
!!$       big_fact_zaa (jpk,jpj,jpi)
!!$       big_fact_zbb(jpk,jpj,jpi)
!!$       big_fact_zcc(jpk,jpj,jpi)
!!$       zbtr_arr(jpk,jpj,jpi)
!!$
!!$
!!$
!!$       e1t(jpj,jpi)
!!$       e2t(jpj,jpi)
!!$       e3t(jpk,jpj,jpi)
!!$       e1u(jpj,jpi)
!!$       e2u(jpj,jpi)
!!$       e3u(jpk,jpj,jpi)
!!$       e1v(jpj,jpi)
!!$       e2v(jpj,jpi)
!!$       e3v(jpk,jpj,jpi)
!!$       un(jpk,jpj,jpi)
!!$       vn(jpk,jpj,jpi)
!!$       wn(jpk,jpj,jpi)


!!$         zaa
!!$         zbb
!!$         zcc
!!$         inv_eu
!!$         inv_ev
!!$         inv_et
!!$         big_fact_zaa
!!$         big_fact_zbb
!!$         big_fact_zcc
!!$         zbtr_arr




  t1=MPI_Wtime()
  !$acc kernels default(present)
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpkm1
           zbtr_arr(jk,jj,ji) = 1./(e1t(jj,ji)*e2t(jj,ji)*e3t(jk,jj,ji))
        END DO
     END DO
  END DO
  !$OMP END TASK
  !$acc end kernels

  !$OMP TASK private(ji,jj) firstprivate(jpim1,jpjm1,jpi,jpj,jpk) default(none) &
  !$OMP shared(zdt,zaa,inv_eu,e1u,e2u,e3u,un,big_fact_zaa)

  !$acc kernels default(present)
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpkm1
           inv_eu(jk,jj,ji) = 1./(e1u(jj,ji)*e2u(jj,ji)*e3u(jk,jj,ji) )
        END DO
     END DO
  END DO
  !$acc end kernels


  !$acc kernels default(present)
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpk
           zaa(jk,jj,ji) = e2u(jj,ji)*e3u(jk,jj,ji) * un(jk,jj,ji)
        END DO
     END DO
  END DO
  !$acc end kernels

  !$acc kernels default(present)
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpk
           big_fact_zaa(jk,jj,ji) = ( abs(zaa(jk,jj,ji)) - zdt*zaa(jk,jj,ji)**2*inv_eu(jk,jj,ji) )
           !/(e1u(jj,ji)*e2u(jj,ji)*e3t(jk,jj,ji) ) )
        END DO
     END DO
  END DO
  !$acc end kernels

  !$OMP END TASK

  !$OMP TASK private(ji,jj) firstprivate(jpim1,jpjm1,jpi,jpj,jpk)  default(none) &
  !$OMP shared(inv_ev,e1v,e2v,e3v,vn,zdt,zbb,big_fact_zbb)

  !$acc kernels default(present)
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpkm1
           inv_ev(jk,jj,ji) = 1./(e1v(jj,ji)*e2v(jj,ji)*e3v(jk,jj,ji) )
        END DO
     END DO
  END DO
  !$acc end kernels

  !$acc kernels default(present)
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpk
           zbb(jk,jj,ji) = e1v(jj,ji)*e3v(jk,jj,ji) * vn(jk,jj,ji)
        END DO
     END DO
  END DO
  !$acc end kernels

  !$acc kernels
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpk
           big_fact_zbb(jk,jj,ji) = ( abs(zbb(jk,jj,ji)) - zdt*zbb(jk,jj,ji)**2*inv_ev(jk,jj,ji) )
        END DO
     END DO
  END DO
  !$OMP END TASK
  !$acc end kernels

  !$OMP TASK private(ji,jj) firstprivate(jpim1,jpjm1,jpi,jpj,jpk) default(none) &
  !$OMP shared(inv_et,e1t,e2t,e3w,wn,zcc,zdt,big_fact_zcc)

  !$acc kernels default(present)
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpkm1
           inv_et(jk,jj,ji) = 1./(e1t(jj,ji)*e2t(jj,ji)*e3w(jk,jj,ji) )
        END DO
     END DO
  END DO
  !$acc end kernels

  !$acc kernels
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpk
           zcc(jk,jj,ji) = e1t(jj,ji)*e2t(jj,ji)* wn(jk,jj,ji)
        END DO
     END DO
  END DO
  !$acc end kernels

  !$acc kernels default(present)
  DO ji = 1,jpi
     DO jj = 1,jpj
        !dir$ vector aligned
        DO jk = 1,jpk
           big_fact_zcc(jk,jj,ji) = ( abs(zcc(jk,jj,ji)) - zdt*zcc(jk,jj,ji)**2*inv_et(jk,jj,ji) )
        END DO
     END DO
  END DO
  !$acc end kernels
  !$OMP END TASK

  !$OMP TASKWAIT
  t2=MPI_Wtime()
  kernel_time=kernel_time+t2-t1


!!!!!$acc update host



  t1=MPI_Wtime()
  !$acc update host( zaa(1:jpk,1:jpj,1:jpi), zbb(1:jpk,1:jpj,1:jpi), zcc(1:jpk,1:jpj,1:jpi) )
  !$acc update host( inv_eu(1:jpk,1:jpj,1:jpi), inv_ev(1:jpk,1:jpj,1:jpi), inv_et(1:jpk,1:jpj,1:jpi) )
  !$acc update host( big_fact_zaa (1:jpk,1:jpj,1:jpi), big_fact_zbb(1:jpk,1:jpj,1:jpi), big_fact_zcc(1:jpk,1:jpj,1:jpi) )
  !$acc update host( zbtr_arr(1:jpk,1:jpj,1:jpi) )

!!$       !$acc update host( e1t(1:jpj,1:jpi), e2t(1:jpj,1:jpi), e3t(1:jpk,1:jpj,1:jpi) )
!!$       !$acc update host( e1u(1:jpj,1:jpi), e2u(1:jpj,1:jpi), e3u(1:jpk,1:jpj,1:jpi) )
!!$       !$acc update host( e1v(1:jpj,1:jpi), e2v(1:jpj,1:jpi), e3v(1:jpk,1:jpj,1:jpi) )
!!$       !$acc update host( un(1:jpk,1:jpj,1:jpi), vn(1:jpk,1:jpj,1:jpi), wn(1:jpk,1:jpj,1:jpi) )
  t2=MPI_Wtime()
  comm_time=comm_time+t2-t1


  t1=MPI_Wtime()
!!!               !$acc exit data delete finalize
!!!!$acc exit data delete( zaa, zbb, zcc, inv_eu, inv_ev, inv_et, big_fact_zaa , big_fact_zbb, big_fact_zcc, zbtr_arr ) finalize
!!!!$acc exit data delete( e1t, e2t, e3t, e1u, e2u, e3u, e1v, e2v, e3v, e3w, un, vn, wn ) finalize
  t2=MPI_Wtime()
  init_time=init_time+t2-t1





  !!     tracer loop parallelized (macrotasking)
  !!     =======================================

  trcadvparttime = MPI_WTIME()


  !!OpenMP compatibility broken. Possibility to use ifndef OpenMP + rename the file in trcadv.F90 to keep it
  t1=MPI_Wtime()
  allocate(zy(jpk,jpj,jpi))
  allocate(zx(jpk,jpj,jpi))
  allocate(zz(jpk,jpj,jpi))
  allocate(ztj(jpk,jpj,jpi))
  allocate(zti(jpk,jpj,jpi))
  allocate(zkx(jpk,jpj,jpi))
  allocate(zky(jpk,jpj,jpi))
  allocate(zkz(jpk,jpj,jpi))
  allocate(zbuf(jpk,jpj,jpi))

  zy(:,:,:) = 0
  zz(:,:,:) = 0
  zx(:,:,:) = 0
  ztj(:,:,:)= 0
  zti(:,:,:)= 0
  zbuf(:,:,:) = 0.
  zkx(:,:,:)=0.
  zky(:,:,:)=0.
  zkz(:,:,:)=0.
  t2=MPI_Wtime()
  cpu_time=cpu_time+t2-t1

  t1=MPI_Wtime()
  !!trn could be allocate earlier
  !$acc enter data create(trn(1:jpk,1:jpj,1:jpi,1:jptra))
  !$acc enter data create(tra(1:jpk,1:jpj,1:jpi,1:jptra))
  !$acc enter data create(advmask(1:jpk,1:jpj,1:jpi))
  !$acc enter data create(flx_ridxt(1:Fsize,1:4))
  !$acc enter data create( diaflx(1:7, 1:Fsize, 1:jptra))

  !$acc enter data create( zy(1:jpk,1:jpj,1:jpi), zx(1:jpk,1:jpj,1:jpi), zz(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( ztj(1:jpk,1:jpj,1:jpi), zti(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( zkx(1:jpk,1:jpj,1:jpi), zky(1:jpk,1:jpj,1:jpi), zkz(1:jpk,1:jpj,1:jpi) )
  !$acc enter data create( zbuf(1:jpk,1:jpj,1:jpi) )
  t2=MPI_Wtime()
  init_time=init_time+t2-t1


  t1=MPI_Wtime()
  !$acc update device(tra(1:jpk,1:jpj,1:jpi,1:jptra))
  !$acc update device(trn(1:jpk,1:jpj,1:jpi,1:jptra))
  !$acc update device(advmask(1:jpk,1:jpj,1:jpi))
  !$acc update device(flx_ridxt(1:Fsize,1:4))
  !$acc update device( diaflx(1:7, 1:Fsize, 1:jptra))

!!$       !$acc update device( zy(1:jpk,1:jpj,1:jpi), zx(1:jpk,1:jpj,1:jpi), zz(1:jpk,1:jpj,1:jpi) )
!!$       !$acc update device( ztj(1:jpk,1:jpj,1:jpi), zti(1:jpk,1:jpj,1:jpi) )
!!$       !$acc update device( zkx(1:jpk,1:jpj,1:jpi), zky(1:jpk,1:jpj,1:jpi), zkz(1:jpk,1:jpj,1:jpi) )
!!$       !$acc update device( zbuf(1:jpk,1:jpj,1:jpi) )
  t2=MPI_Wtime()
  comm_time=comm_time+t2-t1

  t1=MPI_Wtime()
  !$acc kernels default(present)
  DO ji = 1, jpi
     DO jj = 1, jpj
        DO jk = 1, jpk
           zy(jk,jj,ji) = 0
           zz(jk,jj,ji) = 0
           zx(jk,jj,ji) = 0
           ztj(jk,jj,ji)= 0
           zti(jk,jj,ji)= 0
           zbuf(jk,jj,ji) = 0.
           zkx(jk,jj,ji)=0.
           zky(jk,jj,ji)=0.
           zkz(jk,jj,ji)=0.
        ENDDO
     ENDDO
  ENDDO
  !$acc end kernels
  t2=MPI_Wtime()
  kernel_time=kernel_time+t2-t1

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

!!$       !!OpenMP compatibility broken. Possibility to use ifdef OpenMP + rename the file in trcadv.F90 to keep it
!!$       allocate(zy(jpk,jpj,jpi))
!!$       allocate(zx(jpk,jpj,jpi))
!!$       allocate(zz(jpk,jpj,jpi))
!!$       allocate(ztj(jpk,jpj,jpi))
!!$       allocate(zti(jpk,jpj,jpi))
!!$       allocate(zkx(jpk,jpj,jpi))
!!$       allocate(zky(jpk,jpj,jpi))
!!$       allocate(zkz(jpk,jpj,jpi))
!!$       allocate(zbuf(jpk,jpj,jpi))
!!$
!!$       zy(:,:,:) = 0
!!$       zz(:,:,:) = 0
!!$       zx(:,:,:) = 0
!!$       ztj(:,:,:)= 0
!!$       zti(:,:,:)= 0
!!$       zbuf(:,:,:) = 0.
!!$       zkx(:,:,:)=0.
!!$       zky(:,:,:)=0.
!!$       zkz(:,:,:)=0.

!!$         !$acc kernels default(present)
!!$         DO ji = 1, jpi
!!$            DO jj = 1, jpj
!!$               DO jk = 1, jpk
!!$                  zy(jk,jj,ji) = 0
!!$                  zz(jk,jj,ji) = 0
!!$                  zx(jk,jj,ji) = 0
!!$                  ztj(jk,jj,ji)= 0
!!$                  zti(jk,jj,ji)= 0
!!$                  zbuf(jk,jj,ji) = 0.
!!$                  zkx(jk,jj,ji)=0.
!!$                  zky(jk,jj,ji)=0.
!!$                  zkz(jk,jj,ji)=0.
!!$               ENDDO
!!$            ENDDO
!!$         ENDDO
!!$         !$acc end kernels


     !        zkx(  :,:,1)=0.
     !        zkx(:,:,jpi)=0.
     !        zky(:,  1,:)=0.
     !        zky(:,jpj,:)=0.
     !        zkz(1,:,:)  =0.
     ! ! loop unfusion

     t1=MPI_Wtime()
     !$acc kernels default(present)
     DO ji = 2,jpim1
        !dir$ vector aligned
        DO jj = 2,jpjm1
           zkx(1,jj,ji ) = fsx(trn(1,jj,ji, jn),trn(1,jj,ji + 1, jn),zaa(1,jj,ji))
        END DO
     END DO
     !$acc end kernels

     !$acc kernels default(present)
     DO ji = 2,jpim1
        !dir$ vector aligned
        DO jj = 2,jpjm1
           zky(1,jj,ji ) = fsy(trn(1,jj,ji, jn),trn(1,jj+1,ji, jn),zbb(1,jj,ji))
        END DO
     END DO
     !$acc end kernels

     !$acc kernels default(present)
     DO ji = 1,jpi
        !dir$ vector aligned
        DO jk = 2,jpk
           zkz(jk,1,ji ) = fsz(trn(jk,1,ji, jn),trn(jk-1,1,ji, jn),zcc(jk,1,ji))
        ENDDO
     ENDDO
     !$acc end kernels

     !$acc kernels default(present)
     DO ji = 1,jpi
        !dir$ vector aligned
        DO jk = 2,jpk
           zkz(jk,jpj,ji ) = fsz(trn(jk,jpj,ji, jn),trn(jk-1,jpj,ji, jn),zcc(jk,jpj,ji))
        END DO
     END DO
     !$acc end kernels
     ! loop unfusion
     !$acc kernels default(present)
     DO jj = 2,jpjm1
        !dir$ vector aligned
        DO jk = 2,jpk
           zkz(jk,jj,1 ) = fsz(trn(jk,jj,1, jn),trn(jk-1,jj,1, jn),zcc(jk,jj,1))
        END DO
     END DO
     !$acc end kernels

     !$acc kernels default(present)
     DO jj = 2,jpjm1
        !dir$ vector aligned
        DO jk = 2,jpk
           zkz(jk,jj,jpi ) = fsz(trn(jk,jj,jpi, jn),trn(jk-1,jj,jpi, jn),zcc(jk,jj,jpi))
        END DO
     END DO
     !$acc end kernels

     !$acc kernels default(present)
     DO  ji = 2,jpim1
        DO jj = 2,jpjm1
           !dir$ vector aligned
           DO jk = 2,jpk
              zkx(jk,jj,ji ) = fsx(trn(jk,jj,ji, jn),trn(jk,jj,ji + 1, jn),zaa(jk,jj,ji))*advmask(jk,jj,ji)
           END DO
        END DO
     END DO
     !$acc end kernels

     !$acc kernels default(present)
     DO  ji = 2,jpim1
        DO jj = 2,jpjm1
           !dir$ vector aligned
           DO jk = 2,jpk
              zky(jk,jj,ji ) = fsy(trn(jk,jj,ji, jn),trn(jk,jj + 1,ji, jn),zbb(jk,jj,ji))*advmask(jk,jj,ji)
           END DO
        END DO
     END DO
     !$acc end kernels

     !$acc parallel loop collapse(3) gang vector default(present)
     DO  ji = 2,jpim1
        DO jj = 2,jpjm1
           !dir$ vector aligned
           DO jk = 2,jpk
              zkz(jk,jj,ji ) = fsz(trn(jk,jj,ji, jn),trn(jk-1,jj,ji, jn),zcc(jk,jj,ji))*advmask(jk,jj,ji)
           END DO
        END DO
     END DO
     !$acc end parallel loop
     t2=MPI_Wtime()
     kernel_time=kernel_time+t2-t1




!!!       !$acc update host( zkx(1:jpk,1:jpj,1:jpi), zky(1:jpk,1:jpj,1:jpi) )




     ! ... Lateral boundary conditions on zk[xy]
     t1=MPI_Wtime()
#ifdef key_mpp

     !  ... Mpp : export boundary values to neighboring processors

#ifndef _OPENACC
     CALL mpplnk_my(zkx)
     CALL mpplnk_my(zky)
#else
     CALL mpplnk_my_openacc(zkx)
     CALL mpplnk_my_openacc(zky)

#endif

#else
#error

     ! !!   ... T-point, 3D array, full local arrays zk[xy] are initialised

     CALL lbc( zkx(:,:,:), 1, 1, 1, 1, jpk, 1 )
     CALL lbc( zky(:,:,:), 1, 1, 1, 1, jpk, 1 )
#endif
     t2=MPI_Wtime()
     mpi_time=mpi_time+t2-t1

!!!       !$acc update device( zkx(1:jpk,1:jpj,1:jpi), zky(1:jpk,1:jpj,1:jpi) )


     !! 2. calcul of after field using an upstream advection scheme
     !! -----------------------------------------------------------

     t1=MPI_Wtime()
     !$acc kernels default(present)
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
     !$acc end kernels

     !$acc kernels default(present)
     DO jf=1,Fsize
        jk = flx_ridxt(jf,2)
        jj = flx_ridxt(jf,3)
        ji = flx_ridxt(jf,4)

        diaflx(1,jf, jn) = diaflx(1,jf, jn) + zkx(jk,jj,ji )*rdt
        diaflx(2,jf, jn) = diaflx(2,jf, jn) + zky(jk,jj,ji )*rdt
        diaflx(3,jf, jn) = diaflx(3,jf, jn) + zkz(jk,jj,ji )*rdt
     ENDDO
     !$acc end kernels
     t2=MPI_Wtime()
     kernel_time=kernel_time+t2-t1

     !! !$acc update host(trn(1:jpk,1:jpj,1:jpi,1:jptra))
     !! !$acc update host(flx_ridxt(1:Fsize,1:4))

!!$       !$acc update host( diaflx(1:7, 1:Fsize, 1:jptra))
!!$       !$acc update host( zy(1:jpk,1:jpj,1:jpi), zx(1:jpk,1:jpj,1:jpi), zz(1:jpk,1:jpj,1:jpi) )
!!$       !$acc update host( ztj(1:jpk,1:jpj,1:jpi), zti(1:jpk,1:jpj,1:jpi) )
!!$       !$acc update host( zkx(1:jpk,1:jpj,1:jpi), zky(1:jpk,1:jpj,1:jpi), zkz(1:jpk,1:jpj,1:jpi) )
!!$       !$acc update host( zbuf(1:jpk,1:jpj,1:jpi) )



     !! 2.1 start of antidiffusive correction loop

     ANTIDIFF_CORR: DO jt = 1,ncor

        !! 2.2 calcul of intermediary field zti





        t1=MPI_Wtime()
        if(jt .EQ. 1) then

           if(ncor .EQ. 1) then
              !$acc kernels default(present)
              DO ji = 2,jpim1
                 DO jj = 2,jpjm1
                    !dir$ vector aligned
                    DO jk = 1,jpkm1
                       zti(jk,jj,ji ) = trn(jk,jj,ji, jn) + zdt*ztj(jk,jj,ji )
                    END DO
                 END DO
              END DO
              !$acc end kernels

           else
              !$acc kernels default(present)
              DO ji = 2,jpim1
                 DO jj = 2,jpjm1
                    !dir$ vector aligned
                    DO jk = 1,jpkm1
                       zti(jk,jj,ji ) = trn(jk,jj,ji, jn) + zdt*ztj(jk,jj,ji )
                       !  zbuf(jk,jj,ji ) = ztj(jk,jj,ji )
                    END DO
                 END DO
              END DO
              !$acc end kernels

              !$acc kernels default(present)
              DO ji = 2,jpim1
                 DO jj = 2,jpjm1
                    !dir$ vector aligned
                    DO jk = 1,jpkm1
                       !   zti(jk,jj,ji ) = trn(jk,jj,ji, jn) + zdt*ztj(jk,jj,ji )
                       zbuf(jk,jj,ji ) = ztj(jk,jj,ji )
                    END DO
                 END DO
              END DO
              !$acc end kernels

           endif

        else
           !$acc kernels default(present)
           DO ji = 2,jpim1
              DO jj = 2,jpjm1
                 !dir$ vector aligned
                 DO jk = 1,jpkm1
                    zti(jk,jj,ji ) =  zti(jk,jj,ji ) + zdt*ztj(jk,jj,ji )
                 END DO
              END DO
           END DO
           !$acc end kernels

           !$acc kernels default(present)
           DO ji = 2,jpim1
              DO jj = 2,jpjm1
                 !dir$ vector aligned
                 DO jk = 1,jpkm1
                    zbuf(jk,jj,ji ) = zbuf(jk,jj,ji ) + ztj(jk,jj,ji )
                 END DO
              END DO
           END DO
           !$acc end kernels
        endif
        t2=MPI_Wtime()
        kernel_time=kernel_time+t2-t1




        !!       !$acc update host( zti(1:jpk,1:jpj,1:jpi) )
        !!       !$acc update host( zbuf(1:jpk,1:jpj,1:jpi) )



        !! ... Lateral boundary conditions on zti
        t1=MPI_Wtime()
#ifdef key_mpp
        ! ... Mpp : export boundary values to neighboring processors
#ifndef _OPENACC
        CALL mpplnk_my(zti)
#else
        CALL mpplnk_my_openacc(zti)
#endif
#else
        ! ... T-point, 3D array, full local array zti is initialised
        CALL lbc( zti(:,:,:), 1, 1, 1, 1, jpk, 1 )
#endif
        t2=MPI_Wtime()
        mpi_time=mpi_time+t2-t1


        !! 2.3 calcul of the antidiffusive flux

        !!       !$acc update device( zti(1:jpk,1:jpj,1:jpi) )



        t1=MPI_Wtime()

        !jk = 1
        !          DO jk = 1,jpkm1
        !$acc kernels default(present)
        DO ji = 2,jpim1
           DO jj = 2,jpjm1
              junk  = zti(1,jj,ji )
              junki = zti(1,jj,ji+1 )
              junkj = zti(1,jj+1,ji )
              zx(1,jj,ji ) = big_fact_zaa(1,jj,ji)*(junki - junk)/(junk + junki + rtrn)* rsc
              zy(1,jj,ji ) = big_fact_zbb(1,jj,ji)*(junkj - junk)/(junk + junkj + rtrn)* rsc
           END DO
        END DO
        !$acc end kernels

        !DO ju=1, dimen_jarr2
        !$acc kernels default(present)
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
        !$acc end kernels
        t2=MPI_Wtime()
        kernel_time=kernel_time+t2-t1
        !                 endif

        !!   !$acc update host( zy(1:jpk,1:jpj,1:jpi), zx(1:jpk,1:jpj,1:jpi), zz(1:jpk,1:jpj,1:jpi) )



        ! ... Lateral boundary conditions on z[xyz]
        t1=MPI_Wtime()
#ifdef key_mpp

        ! ... Mpp : export boundary values to neighboring processors
#ifndef _OPENACC
        CALL mpplnk_my(zx)
        CALL mpplnk_my(zy)
        CALL mpplnk_my(zz)
#else
        CALL mpplnk_my_openacc(zx)
        CALL mpplnk_my_openacc(zy)
        CALL mpplnk_my_openacc(zz)
#endif
#else

        !  ... T-point, 3D array, full local array z[xyz] are initialised
        CALL lbc( zx(:,:,:), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( zy(:,:,:), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( zz(:,:,:), 1, 1, 1, 1, jpk, 1 )
#endif
        t2=MPI_Wtime()
        mpi_time=mpi_time+t2-t1

        !!     !$acc update device(zy(1:jpk,1:jpj,1:jpi), zx(1:jpk,1:jpj,1:jpi), zz(1:jpk,1:jpj,1:jpi) )

        !! 2.4 reinitialization
        !!            2.5 calcul of the final field:
        !!                advection by antidiffusive mass fluxes and an upstream scheme

        t1=MPI_Wtime()
        !$acc kernels default(present)
        DO ji = 2,jpim1
           !dir$ vector aligned
           DO jj = 2,jpjm1
              zkx(1,jj,ji ) = fsx(zti(1,jj,ji ),zti(1,jj,ji+1 ),zx(1,jj,ji ))
           END DO
        END DO
        !$acc end kernels

        !$acc kernels default(present)
        DO ji = 2,jpim1
           !dir$ vector aligned
           DO jj = 2,jpjm1
              zky(1,jj,ji ) = fsy(zti(1,jj,ji ),zti(1,jj+ 1,ji ),zy(1,jj,ji ))
           END DO
        END DO
        !$acc end kernels

        !$acc kernels default(present)
        DO ji = 1,jpi
           !dir$ vector aligned
           DO jk = 2,jpk
              zkz(jk,1,ji ) = fsz(zti(jk,1,ji ),zti(jk-1,1,ji ),zz(jk,1,ji ))
           ENDDO
        ENDDO
        !$acc end kernels

        !$acc kernels default(present)
        DO ji = 1,jpi
           !dir$ vector aligned
           DO jk = 2,jpk
              zkz(jk,jpj,ji ) = fsz(zti(jk,jpj,ji ),zti(jk-1,jpj,ji ),zz(jk,jpj,ji ))
           ENDDO
        ENDDO
        !$acc end kernels

        !$acc kernels default(present)
        DO jj = 2,jpjm1
           !dir$ vector aligned
           DO jk = 2,jpk
              zkz(jk,jj,1 ) = fsz(zti(jk,jj,1 ),zti(jk-1,jj,1 ),zz(jk,jj,1 ))
           ENDDO
        ENDDO
        !$acc end kernels

        !$acc kernels default(present)
        DO jj = 2,jpjm1
           !dir$ vector aligned
           DO jk = 2,jpk
              zkz(jk,jj,jpi ) = fsz(zti(jk,jj,jpi ),zti(jk-1,jj,jpi ),zz(jk,jj,jpi ))
           END DO
        END DO
        !$acc end kernels

        !$acc kernels default(present)
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
        !$acc end kernels

        !$acc kernels default(present)
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
        !$acc end kernels

        !$acc kernels default(present)
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
        !$acc end kernels
        t2=MPI_Wtime()
        kernel_time=kernel_time+t2-t1


        !!           !$acc update host(zkx(1:jpk,1:jpj,1:jpi), zky(1:jpk,1:jpj,1:jpi) )

        !... Lateral boundary conditions on zk[xy]
        t1=MPI_Wtime()
#ifdef key_mpp
        !  ... Mpp : export boundary values to neighboring processors
#ifndef _OPENACC
        CALL mpplnk_my(zkx)
        CALL mpplnk_my(zky)
#else
        CALL mpplnk_my_openacc(zkx)
        CALL mpplnk_my_openacc(zky)
#endif
#else
        ! ... T-point, 3D array, full local array zk[xy] are initialised
        CALL lbc( zkx(:,:,:), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( zky(:,:,:), 1, 1, 1, 1, jpk, 1 )
#endif
        t2=MPI_Wtime()
        mpi_time=mpi_time+t2-t1

        !!          !$acc update device(zkx(1:jpk,1:jpj,1:jpi), zky(1:jpk,1:jpj,1:jpi) )

        !!        2.6. calcul of after field using an upstream advection scheme



        t1=MPI_Wtime()
        if(ncor .EQ. 1) then
           !$acc kernels default(present)
           DO ji =2,jpim1
              DO jj =2,jpjm1
                 DO jk =1,jpkm1
                    ztj(jk,jj,ji ) = -zbtr_arr(jk,jj,ji)*( zkx(jk,jj,ji ) - zkx(jk,jj,ji - 1 ) &
                         &              + zky(jk,jj,ji ) - zky(jk,jj- 1,ji )+ zkz(jk,jj,ji ) - zkz(jk+1,jj,ji ) )+ ztj(jk,jj,ji )
                 ENDDO
              ENDDO
           ENDDO
           !$acc end kernels

           !$acc kernels default(present)
           DO jf=1,Fsize
              jk = flx_ridxt(jf,2)
              jj = flx_ridxt(jf,3)
              ji = flx_ridxt(jf,4)

              diaflx(1,jf, jn) = diaflx(1,jf, jn) + zkx(jk,jj,ji )*rdt
              diaflx(2,jf, jn) = diaflx(2,jf, jn) + zky(jk,jj,ji )*rdt
              diaflx(3,jf, jn) = diaflx(3,jf, jn) + zkz(jk,jj,ji )*rdt
           ENDDO
           !$acc end kernels


        else

           !$acc kernels default(present)
           DO ji =2,jpim1
              DO jj =2,jpjm1
                 DO jk =1,jpkm1
                    ztj(jk,jj,ji ) = -zbtr_arr(jk,jj,ji)*( zkx(jk,jj,ji ) - zkx(jk,jj,ji - 1 ) &
                         &              + zky(jk,jj,ji ) - zky(jk,jj- 1,ji )+ zkz(jk,jj,ji ) - zkz(jk+1,jj,ji ) )
                 ENDDO
              ENDDO
           ENDDO
           !$acc end kernels

           !$acc kernels default(present)
           DO jf=1,Fsize
              jk = flx_ridxt(jf,2)
              jj = flx_ridxt(jf,3)
              ji = flx_ridxt(jf,4)

              diaflx(1,jf, jn) = diaflx(1,jf, jn) + zkx(jk,jj,ji )*rdt
              diaflx(2,jf, jn) = diaflx(2,jf, jn) + zky(jk,jj,ji )*rdt
              diaflx(3,jf, jn) = diaflx(3,jf, jn) + zkz(jk,jj,ji )*rdt
           ENDDO
           !$acc end kernels

!!$       !$acc update host( diaflx(1:7, 1:Fsize, 1:jptra))
!!$       !$acc update host( ztj(1:jpk,1:jpj,1:jpi) )


        endif
        t2=MPI_Wtime()
        kernel_time=kernel_time+t2-t1

     ENDDO ANTIDIFF_CORR


     !!       3. trend due to horizontal and vertical advection of tracer jn



     t1=MPI_Wtime()
     if(ncor .EQ. 1) then

        !$acc kernels default(present)
        do ji=1,jpi
           do jj=1,jpj
              do jk=1,jpk
                 tra(jk,jj,ji, jn) = tra(jk,jj,ji, jn)+ ztj(jk,jj,ji )
              enddo
           enddo
        enddo
        !$acc end kernels

     else

        !$acc kernels default(present)
        do ji=1,jpi
           do jj=1,jpj
              do jk=1,jpk
                 tra(jk,jj,ji, jn) = tra(jk,jj,ji, jn)+ (zbuf(jk,jj,ji ) + ztj(jk,jj,ji ))
              enddo
           enddo
        enddo
        !$acc end kernels

     endif
     t2=MPI_Wtime()
     kernel_time=kernel_time+t2-t1

!!$       !!OpenMP compatibility broken. Possibility to use ifdef OpenMP + rename the file in trcadv.F90 to keep it
!!$        deallocate(zy )
!!$        deallocate(zx )
!!$        deallocate(zz )
!!$        deallocate(ztj )
!!$        deallocate(zti )
!!$        deallocate(zkx )
!!$        deallocate(zky )
!!$        deallocate(zkz )
!!$        deallocate(zbuf )



  END DO TRACER_LOOP
  !$OMP end taskloop

  t1=MPI_Wtime()
  !$acc update host( diaflx(1:7, 1:Fsize, 1:jptra) )
  !$acc update host( tra(1:jpk,1:jpj,1:jpi,1:jptra) )

  !$acc update host( zy(1:jpk,1:jpj,1:jpi), zx(1:jpk,1:jpj,1:jpi), zz(1:jpk,1:jpj,1:jpi) )
  !$acc update host( ztj(1:jpk,1:jpj,1:jpi), zti(1:jpk,1:jpj,1:jpi) )
  !$acc update host( zkx(1:jpk,1:jpj,1:jpi), zky(1:jpk,1:jpj,1:jpi), zkz(1:jpk,1:jpj,1:jpi) )
  !$acc update host( zbuf(1:jpk,1:jpj,1:jpi) )
  t2=MPI_Wtime()
  comm_time=comm_time+t2-t1

  t1=MPI_Wtime()
  !$acc exit data delete( tra) finalize
  !$acc exit data delete( trn, advmask ) finalize
  !$acc exit data delete( flx_ridxt, diaflx ) finalize
  !$acc exit data delete( zy, zx, zz, ztj, zti, zkx, zky, zkz, zbuf ) finalize
  t2=MPI_Wtime()
  init_time=init_time+t2-t1

  !!OpenMP compatibility broken. Possibility to use ifndef OpenMP + rename the file in trcadv.F90 to keep it
  t1=MPI_Wtime()
  deallocate(zy )
  deallocate(zx )
  deallocate(zz )
  deallocate(ztj )
  deallocate(zti )
  deallocate(zkx )
  deallocate(zky )
  deallocate(zkz )
  deallocate(zbuf )
  t2=MPI_Wtime()
  cpu_time=cpu_time+t2-t1


  t1=MPI_Wtime()
  !$acc exit data delete( zaa, zbb, zcc, inv_eu, inv_ev, inv_et, big_fact_zaa , big_fact_zbb, big_fact_zcc, zbtr_arr ) finalize
  !$acc exit data delete( e1t, e2t, e3t, e1u, e2u, e3u, e1v, e2v, e3v, e3w, un, vn, wn ) finalize
  t2=MPI_Wtime()
  init_time=init_time+t2-t1

  trcadvparttime = MPI_WTIME() - trcadvparttime
  trcadvtottime = trcadvtottime + trcadvparttime
!!!!

  print *, ">>> comm kernel init (ms) ", comm_time*1000._8, kernel_time*1000._8, init_time*1000._8
  print *, ">>> mpi cpu total (ms) ", mpi_time*1000._8, cpu_time*1000._8, (comm_time+kernel_time+init_time+mpi_time+cpu_time)*1000._8

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



END SUBROUTINE trcadv
