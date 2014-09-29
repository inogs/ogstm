       MODULE FN_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
! Common/comcoh/  : False Negative matrix
! ---------------------------------------------------------------------
      INTEGER              :: dimen_jvsnu
      INTEGER, allocatable :: jarr_snu(:,:)
      REAL(8)              :: SMALL
      REAL(8), allocatable :: tra_FN(:,:,:,:)
      REAL(8), allocatable :: TOT(:,:),TOT_FN(:,:),FN_CORR(:,:)


      REAL(8), allocatable :: cou_FN(:,:,:)
      REAL(8), allocatable :: FN_ranking(:)
      INTEGER :: elements,nelements(6),idx_element(14,6)


!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_FN()


      IMPLICIT NONE


      dimen_jvsnu = 0
      SMALL=0.00000000001

      allocate(jarr_snu(2, jpi*jpj))      ; jarr_snu = huge(jarr_snu(1,1))
      allocate(tra_FN(jpi,jpj,jpk,jptra)) ; tra_FN   = huge(tra_FN(1,1,1,1))

      CALL OPA_elements(elements,nelements,idx_element)

      allocate(TOT    (jpj*jpi,elements)) ; TOT        = huge(TOT(1,1))        ! elements is the number of element considered (C,P,N,...)
      allocate(TOT_FN (jpj*jpi,elements)) ; TOT_FN     = huge(TOT_FN(1,1))
      allocate(FN_CORR(jpj*jpi,elements)) ; FN_CORR    = huge(FN_CORR(1,1))


      allocate(cou_FN(jpi,jpj,jpk))       ; cou_FN     = huge(cou_FN(1,1,1))
      allocate(FN_ranking(jptra))         ; FN_ranking = huge(FN_ranking(1))

      mem_all = rea_len*(jpi*jpj*jpk*jptra+2*jpi*jpj*jpk+ jptra)

      tra_FN=0.
      cou_FN=0.
!     cor_FN=0.
      FN_ranking=0.

      END subroutine myalloc_FN


      END MODULE FN_mem
