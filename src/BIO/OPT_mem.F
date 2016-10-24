       MODULE OPT_mem 

       USE modul_param 
       USE myalloc

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif

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
      REAL(4), allocatable :: zkef_f (:,:)
!!!$omp threadprivate (zpar0m,zpar100, zpar, xEPS)
!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_OPT()
      INTEGER  :: err
      REAL(8)  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif
       allocate(itabe(jpi))         
      
       itabe   = huge(itabe(1))
       allocate(imaske(jpk,jpi))   
       imaske  = huge(imaske(1,1))
!!!$omp parallel default (none) shared(jpk,jpi)
       allocate(zpar(jpk,jpi))     
       zpar    = huge(zpar(1,1))
       allocate(xEPS(jpk,jpi))     
       xEPS    = huge(xEPS(1,1))
       allocate(zpar0m(jpi))        
       zpar0m  = huge(zpar0m(1))
       allocate(zpar100(jpi))       
       zpar100 = huge(zpar100(1))
!!!$omp end parallel

       allocate(kef(jpj,jpi))       
       kef     = huge(kef(1,1))
       allocate(kextIO(jpj,jpi,2))  
       kextIO  = huge(kextIO(1,1,1))

#if ! defined  key_kef
       kef(:,:) = 0.04
#endif

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_OPT

      END MODULE 
