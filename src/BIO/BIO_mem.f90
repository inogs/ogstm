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
      double precision, allocatable :: ogstm_sedipi(:,:,:,:)
      double precision, allocatable :: ogstm_ph(:,:,:) ! GUESS for FOLLOWS algorithm
      double precision, allocatable :: NPPF2(:,:,:)
      double precision, allocatable :: ogstm_co2(:,:), co2_IO(:,:,:)
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
       allocate(ogstm_co2(jpj,jpi))          
       ogstm_co2       = huge(ogstm_co2(1,1))
       allocate(co2_IO(jpj,jpi,2))    
        
       co2_IO    = huge(co2_IO(1,1,1))
       allocate(ogstm_sedipi(jpk,jpj,jpi,4)) 
       ogstm_sedipi    = huge(ogstm_sedipi(1,1,1,1))
       allocate(ogstm_ph(jpk,jpj,jpi))       
       ogstm_ph        = huge(ogstm_ph(1,1,1))
       ogstm_ph=8.0
       allocate(NPPF2(jpk,jpj,jpi))
       NPPF2 = 0 ! nut huge, because it will be assigned only in trcBIO in BFMpoints
                 ! and used in hard_tissue_pump.F also in land points
       ice=0

#ifdef Mem_Monitor
       mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_BIO



      subroutine clean_memory_bio()

            deallocate(bfm_trn)
            deallocate(bfm_tra)
            deallocate(surf_mask)
            deallocate(ogstm_co2)
            deallocate(co2_IO)
            deallocate(ogstm_sedipi)
            deallocate(ogstm_ph)
            deallocate(NPPF2)

      end subroutine clean_memory_bio



      END MODULE 
