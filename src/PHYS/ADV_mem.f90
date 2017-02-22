       MODULE ADV_mem 

       USE modul_param 
       USE myalloc
       ! epascolo USE myalloc_mpp
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
      INTEGER, allocatable,dimension(:) :: jarr_adv_flx
      INTEGER, allocatable,dimension(:,:) :: jarr, jarr1, jarr2, jarr3, jarrt
      INTEGER(kind=1), allocatable,dimension(:,:,:) :: advmask
      ! double precision, allocatable,dimension(:,:,:),save :: zti,ztj
      ! double precision, allocatable,dimension(:,:,:),save :: zx,zy,zz,zbuf
      ! double precision, allocatable,dimension(:,:,:),save :: zkx,zky,zkz
      double precision, allocatable,dimension(:,:,:) :: zaa,zbb,zcc

      double precision, allocatable,dimension(:,:,:) :: inv_eu, inv_ev, inv_et 
      double precision, allocatable,dimension(:,:,:) :: big_fact_zaa, big_fact_zbb, big_fact_zcc 
      double precision, allocatable,dimension(:,:,:) :: zbtr_arr

!!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_ADV()

#ifdef __OPENMP1
      INTEGER :: ntids, omp_get_max_threads
      EXTERNAL :: omp_get_max_threads
#else
      INTEGER :: ntids 
#endif

      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

#ifdef __OPENMP1
      ntids = omp_get_max_threads() ! take the number of threads
#endif

       allocate(advmask(jpk,jpj,jpi))       
       advmask      = huge(advmask(1,1,1))
      !  allocate(zti(jpk,jpj,jpi))    
      !  zti          = huge(zti(1,1,1))
      !  allocate(ztj(jpk,jpj,jpi))    
      !  ztj          = huge(ztj(1,1,1))
       allocate(zaa(jpk,jpj,jpi))           
       zaa          = huge(zaa(1,1,1))
       allocate(zbb(jpk,jpj,jpi))           
       zbb          = huge(zbb(1,1,1))
       allocate(zcc(jpk,jpj,jpi))           
       zcc          = huge(zcc(1,1,1))
      !  allocate(zx(jpk,jpj,jpi))     
       !zx           = huge(zx(1,1,1))
       !allocate(zy(jpk,jpj,jpi))     
       !zy           = huge(zy(1,1,1))
      !  allocate(zz(jpk,jpj,jpi))     
      !  zz           = huge(zz(1,1,1))

      !  allocate(zbuf(jpk,jpj,jpi))    
      !  zbuf         = huge(zbuf(1,1,1))

      !  allocate(zkx(jpk,jpj,jpi))    
      !  zkx          = huge(zkx(1,1,1))
      !  allocate(zky(jpk,jpj,jpi))    
      !  zky          = huge(zky(1,1,1))
      !  allocate(zkz(jpk,jpj,jpi))    
      !  zkz          = huge(zkz(1,1,1))

       allocate(inv_eu(jpk,jpj,jpi))        
       inv_eu       = huge(inv_eu(1,1,1))
       allocate(inv_ev(jpk,jpj,jpi))        
       inv_ev       = huge(inv_ev(1,1,1))
       allocate(inv_et(jpk,jpj,jpi))        
       inv_et       = huge(inv_et(1,1,1))
       allocate(big_fact_zaa (jpk,jpj,jpi)) 
       big_fact_zaa = huge(big_fact_zaa(1,1,1))
       allocate(big_fact_zbb(jpk,jpj,jpi))  
       big_fact_zbb = huge(big_fact_zbb(1,1,1))
       allocate(big_fact_zcc(jpk,jpj,jpi))  
       big_fact_zcc = huge(big_fact_zcc(1,1,1))
       allocate(zbtr_arr(jpk,jpj,jpi))      
       zbtr_arr     = huge(zbtr_arr(1,1,1))
       allocate(jarr(3, jpk*jpj*jpi))         
       jarr         = huge(jarr(1,1))
       allocate(jarr1(3, jpk*jpj*jpi))        
       jarr1        = huge(jarr1(1,1))
       allocate(jarr2(3, jpk*jpj*jpi))        
       jarr2        = huge(jarr2(1,1))
       allocate(jarr3(3, jpk*jpj*jpi))        
       jarr3        = huge(jarr3(1,1))
       allocate(jarrt(3, jpk*jpj*jpi))        
       jarrt        = huge(jarrt(1,1))
       allocate(jarr_adv_flx(jpk*jpj*jpi))    
       jarr_adv_flx = huge(jarr_adv_flx(1))


#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_ADV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
     



      END MODULE 
