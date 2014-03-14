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
       allocate(jarr_zdf(2, jpi*jpj))
       allocate(jarr_zdf_flx(jpi*jpj,jpk))
       allocate(zwd(jpk, ntids))
       allocate(zws(jpk, ntids))
       allocate(zwi(jpk, ntids))
       allocate(zwx(jpk, ntids))
       allocate(zwy(jpk, ntids))
       allocate(zwz(jpk, ntids))
       allocate(zwt(jpk, ntids))
       mem_all= rea_len*jpk *ntids*7+ int_len *2*jpi* jpj

      END subroutine myalloc_ZDF

      END MODULE 
