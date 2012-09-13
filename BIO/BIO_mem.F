       MODULE BIO_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

      REAL(8), allocatable :: bfm_trn(:), bfm_tra(:)
!$omp  threadprivate(bfm_trn,bfm_tra)
      REAL(8), allocatable :: surf_mask(:)
      REAL(8), allocatable :: sediPI(:,:,:,:)
      REAL(8), allocatable :: PH(:,:,:) ! GUESS for FOLLOWS algorithm
      REAL(8):: opa_ice, opa_co2, opa_co2_start, opa_co2_end

CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_BIO()

!$omp parallel default(none)
       allocate(bfm_trn(jptra))
       allocate(bfm_tra(jptra))
!$omp end parallel
       allocate(surf_mask(jpk))
       allocate(sediPI(jpi,jpj,jpk,4))
       allocate(PH(jpi,jpj,jpk))
       PH=8.0

       opa_ice=0
       opa_co2_start=290.0
       opa_co2_end=370.0
!       sediPI = 0 ! la biologia non assegna valori sui punti di terra
       sediPI = 1000000 ! la biologia non assegna valori sui punti di terra
!       write(*,*) 'sediPI =', sediPI
       mem_all= rea_len*(2*jptra+jpk+ jpi*jpj*jpk*4)

      END subroutine myalloc_BIO

      END MODULE 
