       MODULE ADV_mem 

       USE modul_param 
       USE myalloc
       USE myalloc_mpp
       USE DIA_mem

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif

       IMPLICIT NONE

       public

      INTEGER :: goodpoints, allpoints, tpoints
      INTEGER :: dimen_jarr, dimen_jarr1, dimen_jarr2, dimen_jarr3, dimen_jarrt 
      INTEGER :: myji, myjj, myjk, locsum
      INTEGER :: jilef, jjlef, jklef, jirig, jjrig, jkrig
      INTEGER, allocatable :: advmask(:,:,:)
      INTEGER, allocatable :: jarr(:,:), jarr1(:,:), jarr2(:,:), jarr3(:,:), jarrt(:,:), jarr_adv_flx(:)
      REAL(8), allocatable :: zti(:,:,:,:),ztj(:,:,:,:)
      REAL(8), allocatable :: zaa(:,:,:),zbb(:,:,:),zcc(:,:,:)
      REAL(8), allocatable :: zx(:,:,:,:),zy(:,:,:,:),zz(:,:,:,:)

      REAL(8), allocatable :: zbuf(:,:,:,:)


      REAL(8), allocatable :: zkx(:,:,:,:),zky(:,:,:,:),zkz(:,:,:,:)

      REAL(8), allocatable :: inv_eu(:,:,:), inv_ev(:,:,:), inv_et(:,:,:) 
      REAL(8), allocatable :: big_fact_zaa(:,:,:), big_fact_zbb(:,:,:), big_fact_zcc(:,:,:) 
      REAL(8), allocatable :: zbtr_arr(:,:,:)

CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_ADV()

#ifdef __OPENMP
      INTEGER :: ntids, omp_get_max_threads
      EXTERNAL :: omp_get_max_threads
#else
      INTEGER :: ntids 
#endif

      INTEGER  :: err
      REAL(8)  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
#else
      ntids = mpi_pack_size
#endif

       allocate(advmask(jpi, jpj, jpk))       ; advmask      = huge(advmask(1,1,1))
       allocate(zti(jpi, jpj, jpk, ntids))    ; zti          = huge(zti(1,1,1,1))
       allocate(ztj(jpi, jpj, jpk, ntids))    ; ztj          = huge(ztj(1,1,1,1))
       allocate(zaa(jpi, jpj, jpk))           ; zaa          = huge(zaa(1,1,1))
       allocate(zbb(jpi, jpj, jpk))           ; zbb          = huge(zbb(1,1,1))
       allocate(zcc(jpi, jpj, jpk))           ; zcc          = huge(zcc(1,1,1))
       allocate(zx(jpi, jpj, jpk, ntids))     ; zx           = huge(zx(1,1,1,1))
       allocate(zy(jpi, jpj, jpk, ntids))     ; zy           = huge(zy(1,1,1,1))
       allocate(zz(jpi, jpj, jpk, ntids))     ; zz           = huge(zz(1,1,1,1))

       allocate(zbuf(jpi, jpj, jpk,ntids))    ; zbuf         = huge(zbuf(1,1,1,1))

       allocate(zkx(jpi, jpj, jpk, ntids))    ; zkx          = huge(zkx(1,1,1,1))
       allocate(zky(jpi, jpj, jpk, ntids))    ; zky          = huge(zky(1,1,1,1))
       allocate(zkz(jpi, jpj, jpk, ntids))    ; zkz          = huge(zkz(1,1,1,1))

       allocate(inv_eu(jpi, jpj, jpk))        ; inv_eu       = huge(inv_eu(1,1,1))
       allocate(inv_ev(jpi, jpj, jpk))        ; inv_ev       = huge(inv_ev(1,1,1))
       allocate(inv_et(jpi, jpj, jpk))        ; inv_et       = huge(inv_et(1,1,1))
       allocate(big_fact_zaa (jpi, jpj, jpk)) ; big_fact_zaa = huge(big_fact_zaa(1,1,1))
       allocate(big_fact_zbb(jpi, jpj, jpk))  ; big_fact_zbb = huge(big_fact_zbb(1,1,1))
       allocate(big_fact_zcc(jpi, jpj, jpk))  ; big_fact_zcc = huge(big_fact_zcc(1,1,1))
       allocate(zbtr_arr(jpi, jpj, jpk))      ; zbtr_arr     = huge(zbtr_arr(1,1,1))
       allocate(jarr(3, jpi*jpj*jpk))         ; jarr         = huge(jarr(1,1))
       allocate(jarr1(3, jpi*jpj*jpk))        ; jarr1        = huge(jarr1(1,1))
       allocate(jarr2(3, jpi*jpj*jpk))        ; jarr2        = huge(jarr2(1,1))
       allocate(jarr3(3, jpi*jpj*jpk))        ; jarr3        = huge(jarr3(1,1))
       allocate(jarrt(3, jpi*jpj*jpk))        ; jarrt        = huge(jarrt(1,1))
       allocate(jarr_adv_flx(jpi*jpj*jpk))    ; jarr_adv_flx = huge(jarr_adv_flx(1))


#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_ADV

      REAL(8) FUNCTION fsx( pfx1, pfx2, pfu )
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: pfx1, pfx2, pfu
      REAL(8) abspfu
      abspfu = abs(pfu)
      fsx = ( ( pfu + abspfu ) * pfx1+( pfu - abspfu ) * pfx2 ) * 0.5
      END FUNCTION fsx

      REAL(8) FUNCTION fsy( pfy1, pfy2, pfv  )
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: pfy1, pfy2, pfv
      REAL(8) abspfv
      abspfv = abs(pfv)
      fsy = ( ( pfv + abspfv ) * pfy1 +( pfv - abspfv ) * pfy2 ) * 0.5
      END FUNCTION fsy


      REAL(8) FUNCTION fsz( pfz1, pfz2, pfw )
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: pfz1, pfz2, pfw
      REAL(8) abspfw
      abspfw = abs(pfw)
      fsz = ( ( pfw + abspfw ) * pfz1+( pfw - abspfw ) * pfz2 ) * 0.5
      END FUNCTION fsz

      END MODULE 
