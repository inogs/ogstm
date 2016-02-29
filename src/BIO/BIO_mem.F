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

      REAL(8), allocatable :: bfm_trn(:), bfm_tra(:)
!$omp  threadprivate(bfm_trn,bfm_tra)
      REAL(8), allocatable :: surf_mask(:)
      REAL(8), allocatable :: sediPI(:,:,:,:)
      REAL(8), allocatable :: PH(:,:,:) ! GUESS for FOLLOWS algorithm
      REAL(8), allocatable :: NPPF2(:,:,:) 
      REAL(8), allocatable :: co2(:,:), co2_IO(:,:,:)
      REAL(8):: ice!, opa_co2, opa_co2_start, opa_co2_end


CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_BIO()

      INTEGER  :: err
      REAL(8)  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

!$omp parallel default(none)
       allocate(bfm_trn(jptra))         
       bfm_trn   = huge(bfm_trn(1))
       allocate(bfm_tra(jptra))         
       bfm_tra   = huge(bfm_tra(1))
!$omp end parallel
       allocate(surf_mask(jpk))         
       surf_mask = huge(surf_mask(1))  
       allocate(co2(jpi,jpj))           
       co2       = huge(co2(1,1))  
       allocate(co2_IO(jpi,jpj,2))      
       co2_IO    = huge(co2_IO(1,1,1)) 
       allocate(sediPI(jpi,jpj,jpk,4))  
       sediPI    = huge(sediPI(1,1,1,1)) 
       allocate(PH(jpi,jpj,jpk))        
       PH        = huge(PH(1,1,1)) 
       PH=8.0
       allocate(NPPF2(jpi,jpj,jpk))
       NPPF2 = 0.0 ! nut huge, because it will be assigned only in trcBIO in BFMpoints 
                   ! and used in hard_tissue_pump.F also in land points

       ice=0

#ifdef Mem_Monitor
       mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_BIO

      END MODULE 
