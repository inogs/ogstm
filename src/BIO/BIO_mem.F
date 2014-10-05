       MODULE BIO_mem 

       USE modul_param 
       USE myalloc
       USE TIME_MANAGER

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif Mem_Monitor


       IMPLICIT NONE

       public

      REAL(8), allocatable :: bfm_trn(:), bfm_tra(:)
!$omp  threadprivate(bfm_trn,bfm_tra)
      REAL(8), allocatable :: surf_mask(:)
      REAL(8), allocatable :: sediPI(:,:,:,:)
      REAL(8), allocatable :: PH(:,:,:) ! GUESS for FOLLOWS algorithm
      REAL(8), allocatable :: opa_co2(:,:), co2_IO(:,:,:)
      REAL(8):: opa_ice!, opa_co2, opa_co2_start, opa_co2_end


CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_BIO()

      INTEGER  :: err
      INTEGER  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif Mem_Monitor

!$omp parallel default(none)
       allocate(bfm_trn(jptra))        ; bfm_trn   = huge(bfm_trn(1))
       allocate(bfm_tra(jptra))        ; bfm_tra   = huge(bfm_tra(1))
!$omp end parallel
       allocate(surf_mask(jpk))        ; surf_mask = huge(surf_mask(1))
       allocate(opa_co2(jpi,jpj))      ; opa_co2   = huge(opa_co2(1,1))
       allocate(co2_IO(jpi,jpj,2))     ; co2_IO    = huge(co2_IO(1,1,1))
       allocate(sediPI(jpi,jpj,jpk,4)) ; sediPI    = huge(sediPI(1,1,1,1))
       allocate(PH(jpi,jpj,jpk))       ; PH        = huge(PH(1,1,1))
       PH=8.0

       opa_ice=0
       !opa_co2_start=290.0
       !opa_co2_end=370.0

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif Mem_Monitor

      END subroutine myalloc_BIO

      END MODULE 
