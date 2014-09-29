       MODULE ZDF_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public


      INTEGER :: dimen_jvzdf
      INTEGER, allocatable :: jarr_zdf(:,:),jarr_zdf_flx(:,:)
      REAL(8), allocatable :: zwd(:,:), zws(:,:), zwi(:,:)
      REAL(8), allocatable :: zwx(:,:), zwy(:,:), zwz(:,:), zwt(:,:)


CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_ZDF()

#ifdef __OPENMP
      INTEGER :: ntids, omp_get_max_threads
      EXTERNAL :: omp_get_max_threads
#else
      INTEGER :: ntids = 1
#endif

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
#endif
       dimen_jvzdf=0

       allocate(jarr_zdf(2, jpi*jpj))      ; jarr_zdf     = huge(jarr_zdf(1,1))
       allocate(jarr_zdf_flx(jpi*jpj,jpk)) ; jarr_zdf_flx = huge(jarr_zdf_flx(1,1))
       allocate(zwd(jpk, ntids))           ; zwd          = huge(zwd(1,1))
       allocate(zws(jpk, ntids))           ; zws          = huge(zws(1,1))
       allocate(zwi(jpk, ntids))           ; zwi          = huge(zwi(1,1))
       allocate(zwx(jpk, ntids))           ; zwx          = huge(zwx(1,1))
       allocate(zwy(jpk, ntids))           ; zwy          = huge(zwy(1,1))
       allocate(zwz(jpk, ntids))           ; zwz          = huge(zwz(1,1))
       allocate(zwt(jpk, ntids))           ; zwt          = huge(zwt(1,1))
       mem_all= rea_len*jpk *ntids*7+ int_len *2*jpi* jpj

      END subroutine myalloc_ZDF

      END MODULE 
