       MODULE FN_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
! Common/comcoh/  : False Negative matrix
! ---------------------------------------------------------------------


      REAL(8), allocatable :: tra_FN(:,:,:,:)
      REAL(8), allocatable :: cou_FN(:,:,:)
      REAL(8), allocatable :: cor_FN(:,:,:)
      REAL(8), allocatable :: FN_ranking(:)

!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_FN()

      allocate(tra_FN(jpi,jpj,jpk,jptra))
      allocate(cou_FN(jpi,jpj,jpk))
      allocate(cor_FN(jpi,jpj,jpk))
      allocate(FN_ranking(jptra))
      mem_all = rea_len*(jpi*jpj*jpk*jptra+2*jpi*jpj*jpk+ jptra)

      tra_FN=0.
      cou_FN=0.
      cor_FN=0.
      FN_ranking=0.

      END subroutine myalloc_FN


      END MODULE 
