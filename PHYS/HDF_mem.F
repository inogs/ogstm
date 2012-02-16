       MODULE HDF_mem

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

      INTEGER :: dimen_jvhdf1,dimen_jvhdf2,dimen_jvhdf3
      INTEGER, allocatable :: hdfmask(:,:,:)
      INTEGER, allocatable :: jarr_hdf(:,:,:)
      REAL(8), allocatable :: zta(:)
      REAL(8), allocatable :: zeeu(:,:,:), zeev(:,:,:), zbtr(:,:)
      REAL(8), allocatable :: zlt(:,:,:,:), ztu(:,:,:,:), ztv(:,:,:,:)

CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_HDF()

#ifdef __OPENMP
      INTEGER :: ntids, omp_get_max_threads
      EXTERNAL :: omp_get_max_threads
#else
      INTEGER :: ntids = 1
#endif

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
#endif


       dimen_jvhdf1=0
       dimen_jvhdf2=0
       dimen_jvhdf3=0

       allocate(jarr_hdf(3,jpj*jpj*jpk,2))
       allocate(hdfmask(jpi,jpj,jpk   ))
       allocate(zta    (            ntids))
       allocate(zeeu   (jpi,jpj,jpk      ))
       allocate(zeev   (jpi,jpj,jpk      ))
       allocate(zbtr   (jpi,jpj          ))
       allocate(zlt    (jpi,jpj,jpk,ntids))
       allocate(ztu    (jpi,jpj,jpk,ntids))
       allocate(ztv    (jpi,jpj,jpk,ntids))


       mem_all= int_len*7*jpi*jpj*jpk + rea_len*jpi*jpj + 2*rea_len*jpi*jpj*jpk+
     .           3*rea_len*jpi*jpj*jpk*ntids+rea_len*ntids

      END subroutine myalloc_HDF

      END MODULE 
