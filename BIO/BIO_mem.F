       MODULE BIO_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

CC----------------------------------------------------------------------
CC Common/comcoh/  : BIO aux matrix matrix
CC ---------------------------------------------------------------------
C
C
      REAL(8), allocatable :: bfm_trn(:), bfm_tra(:)
!$omp  threadprivate(bfm_trn,bfm_tra)
      REAL(8), allocatable :: surf_mask(:)
      REAL(8), allocatable :: sediPI(:,:,:,:)

CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_BIO()

!$omp parallel default(none)
       allocate(bfm_trn(jptra))
       allocate(bfm_tra(jptra))
!$omp end parallel
       allocate(surf_mask(jpk))
       allocate(sediPI(jpi,jpj,jpk,4))
       sediPI = 0 ! altrimenti la biologia non dà valori sui punti di terra
       mem_all= rea_len*(2*jptra+jpk+ jpi*jpj*jpk*4)

      END subroutine myalloc_BIO

      END MODULE 
