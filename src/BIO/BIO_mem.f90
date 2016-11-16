       MODULE BIO_mem 

       USE modul_param 
       USE myalloc
       USE TIME_MANAGER

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif


       IMPLICIT NONE

       public

      double precision, allocatable :: bfm_trn(:), bfm_tra(:)
!!!$omp  threadprivate(bfm_trn,bfm_tra)
      double precision, allocatable :: surf_mask(:)
      double precision, allocatable :: sediPI(:,:,:,:)
      double precision, allocatable :: PH(:,:,:) ! GUESS for FOLLOWS algorithm
      double precision, allocatable :: co2(:,:), co2_IO(:,:,:)
      double precision:: ice


!!!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_BIO()

      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

!!!$omp parallel default(none)
       allocate(bfm_trn(jptra))        
       bfm_trn   = huge(bfm_trn(1))
       allocate(bfm_tra(jptra))        
       bfm_tra   = huge(bfm_tra(1))
!!!$omp end parallel
       allocate(surf_mask(jpk))        
       surf_mask = huge(surf_mask(1))
       allocate(co2(jpj,jpi))          
       co2       = huge(co2(1,1))
       allocate(co2_IO(jpj,jpi,2))    
        
       co2_IO    = huge(co2_IO(1,1,1))
       allocate(sediPI(jpk,jpj,jpi,4)) 
       sediPI    = huge(sediPI(1,1,1,1))
       allocate(PH(jpk,jpj,jpi))       
       PH        = huge(PH(1,1,1))
       PH=8.0

       ice=0

#ifdef Mem_Monitor
       mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_BIO

      END MODULE 
