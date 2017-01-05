       MODULE HDF_mem

       USE modul_param 
       USE myalloc
       ! epascolo USE myalloc_mpp

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif

       IMPLICIT NONE

       public

      INTEGER :: dimen_jvhdf1,dimen_jvhdf2,dimen_jvhdf3
      INTEGER, allocatable :: hdfmask(:,:,:)
      INTEGER, allocatable :: jarr_hdf(:,:,:),jarr_hdf_flx(:)
      double precision :: zta
      double precision, allocatable :: zeeu(:,:,:), zeev(:,:,:), zbtr(:,:,:)
      double precision, allocatable,dimension(:,:,:) :: zlt, ztu, ztv

!!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_HDF()

      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif



       dimen_jvhdf1=0
       dimen_jvhdf2=0
       dimen_jvhdf3=0

       allocate(jarr_hdf(3,jpk*jpj*jpi,2))  
       
       jarr_hdf     = huge(jarr_hdf(1,1,1))
       allocate(jarr_hdf_flx(jpk*jpj*jpi))  
       jarr_hdf_flx = huge(jarr_hdf_flx(1))
       allocate(hdfmask(jpk,jpj,jpi   ))    
       hdfmask      = huge(hdfmask(1,1,1))
       
       allocate(zeeu   (jpk,jpj,jpi      )) 
       zeeu         = huge(zeeu(1,1,1))
       allocate(zeev   (jpk,jpj,jpi      )) 
       zeev         = huge(zeev(1,1,1))
       allocate(zbtr   (jpk,jpj,jpi      )) 
       zbtr         = huge(zbtr(1,1,1)) 
       allocate(zlt    (jpk,jpj,jpi)) 
       zlt          = huge(zlt(1,1,1)) 
       allocate(ztu    (jpk,jpj,jpi)) 
       ztu          = huge(ztu(1,1,1))
       allocate(ztv    (jpk,jpj,jpi)) 
       ztv          = huge(ztv(1,1,1))


       ztu = 0.
       ztv = 0.

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_HDF

      END MODULE 
