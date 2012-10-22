       MODULE ADV_mem 

       USE modul_param 
       USE myalloc
       USE DIA_mem

       IMPLICIT NONE

       public

CC----------------------------------------------------------------------
CC Common/comcoh/  : ADVection matrix
CC ---------------------------------------------------------------------
C
C
      INTEGER :: goodpoints, allpoints, tpoints
      INTEGER :: dimen_jarr, dimen_jarr1, dimen_jarr2, dimen_jarr3, dimen_jarrt 
      INTEGER :: myji, myjj, myjk, locsum
      INTEGER :: jilef, jjlef, jklef, jirig, jjrig, jkrig
      INTEGER, allocatable :: advmask(:,:,:)
      INTEGER, allocatable :: jarr(:,:), jarr1(:,:), jarr2(:,:), jarr3(:,:), jarrt(:,:), jarr_adv_flx(:)
      REAL(8), allocatable :: zti(:,:,:,:),ztj(:,:,:,:)
      REAL(8), allocatable :: zaa(:,:,:),zbb(:,:,:),zcc(:,:,:)
      REAL(8), allocatable :: zx(:,:,:,:),zy(:,:,:,:),zz(:,:,:,:)

      REAL(8), allocatable :: zbuf(:,:,:)
!$omp threadprivate( zbuf )

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
      INTEGER :: ntids = 1
#endif

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
#endif

       allocate(advmask(jpi, jpj, jpk))
       allocate(zti(jpi, jpj, jpk, ntids))
       allocate(ztj(jpi, jpj, jpk, ntids))
       allocate(zaa(jpi, jpj, jpk))
       allocate(zbb(jpi, jpj, jpk))
       allocate(zcc(jpi, jpj, jpk))
       allocate(zx(jpi, jpj, jpk, ntids))
       allocate(zy(jpi, jpj, jpk, ntids))
       allocate(zz(jpi, jpj, jpk, ntids))
!$omp parallel default(none) shared(jpi, jpj, jpk)
       allocate(zbuf(jpi, jpj, jpk))
!$omp end parallel
       allocate(zkx(jpi, jpj, jpk, ntids))
       allocate(zky(jpi, jpj, jpk, ntids))
       allocate(zkz(jpi, jpj, jpk, ntids))

       allocate(inv_eu(jpi, jpj, jpk))
       allocate(inv_ev(jpi, jpj, jpk))
       allocate(inv_et(jpi, jpj, jpk))
       allocate(big_fact_zaa (jpi, jpj, jpk))
       allocate(big_fact_zbb(jpi, jpj, jpk))
       allocate(big_fact_zcc(jpi, jpj, jpk))
       allocate(zbtr_arr(jpi, jpj, jpk))
       allocate(jarr(3, jpi*jpj*jpk))
       allocate(jarr1(3, jpi*jpj*jpk))
       allocate(jarr2(3, jpi*jpj*jpk))
       allocate(jarr3(3, jpi*jpj*jpk))
       allocate(jarrt(3, jpi*jpj*jpk))
       allocate(jarr_adv_flx(jpi*jpj*jpk))
       mem_all= rea_len*35 * jpi* jpj* jpk+ int_len *jpi* jpj* jpk

      END subroutine myalloc_ADV

      END MODULE 
