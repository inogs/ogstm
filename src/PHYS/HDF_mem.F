       MODULE HDF_mem

       USE modul_param 
       USE myalloc
       USE myalloc_mpp

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif Mem_Monitor

       IMPLICIT NONE

       public

      INTEGER :: dimen_jvhdf1,dimen_jvhdf2,dimen_jvhdf3
      INTEGER, allocatable :: hdfmask(:,:,:)
      INTEGER, allocatable :: jarr_hdf(:,:,:),jarr_hdf_flx(:)
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
      INTEGER :: ntids 
#endif
      INTEGER  :: err
      INTEGER  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif Mem_Monitor

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
#else
      ntids = mpi_pack_size
#endif


       dimen_jvhdf1=0
       dimen_jvhdf2=0
       dimen_jvhdf3=0

       allocate(jarr_hdf(3,jpi*jpj*jpk,2))  ; jarr_hdf     = huge(jarr_hdf(1,1,1))
       allocate(jarr_hdf_flx(jpi*jpj*jpk))  ; jarr_hdf_flx = huge(jarr_hdf_flx(1))
       allocate(hdfmask(jpi,jpj,jpk   ))    ; hdfmask      = huge(hdfmask(1,1,1))
       allocate(zta    (            ntids)) ; zta          = huge(zta(1))
       allocate(zeeu   (jpi,jpj,jpk      )) ; zeeu         = huge(zeeu(1,1,1))
       allocate(zeev   (jpi,jpj,jpk      )) ; zeev         = huge(zeev(1,1,1))
       allocate(zbtr   (jpi,jpj          )) ; zbtr         = huge(zbtr(1,1)) 
       allocate(zlt    (jpi,jpj,jpk,ntids)) ; zlt          = huge(zlt(1,1,1,1)) 
       allocate(ztu    (jpi,jpj,jpk,ntids)) ; ztu          = huge(ztu(1,1,1,1))
       allocate(ztv    (jpi,jpj,jpk,ntids)) ; ztv          = huge(ztv(1,1,1,1))


       ztu = 0.
       ztv = 0.

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif Mem_Monitor

      END subroutine myalloc_HDF

      END MODULE 
