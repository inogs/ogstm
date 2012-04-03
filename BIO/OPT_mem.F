       MODULE OPT_mem 

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
! Common/comcoh/  : ADVection matrix
! ---------------------------------------------------------------------


      INTEGER, allocatable :: itabe(:),imaske(:,:) 
      REAL(8), allocatable :: zpar(:,:),xEPS(:,:)
      REAL(8), allocatable :: zpar0m(:),zpar100(:) 
      REAL(8), allocatable :: kef(:,:)
      REAL(8), allocatable :: kextIO(:,:,:)
!      REAL(4), allocatable :: zkef_f (:,:)
!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_OPT()

       allocate(itabe(jpi))
       allocate(imaske(jpi, jpk))
       allocate(zpar(jpi, jpk))
       allocate(xEPS(jpi, jpk))
       allocate(zpar0m(jpi))
       allocate(zpar100(jpi))
       allocate(kef(jpi,jpj))
       allocate(kextIO(jpi,jpj,2))
!       allocate(zkef_f(jpiglo,jpjglo))

       mem_all = rea_len*(4*jpi*jpk+ 2*jpi+jpi + 2*jpi*jpj)
#if ! defined  key_kef
       kef(:,:) = 0.04
#endif
      END subroutine myalloc_OPT

      END MODULE 
