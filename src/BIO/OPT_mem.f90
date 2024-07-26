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
      ! double precision, allocatable :: zpar(:,:)
      double precision, allocatable :: xEPS_ogstm(:,:)
      ! double precision, allocatable :: zpar0m(:),zpar100(:)
      double precision, allocatable :: kef(:,:)
      double precision, allocatable :: kextIO(:,:,:)
      real, allocatable :: zkef_f (:,:)
!!!$omp threadprivate (zpar0m,zpar100, zpar, xEPS_ogstm)
!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_OPT()
      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif
       allocate(itabe(jpi))         
      
       itabe   = huge(itabe(1))
       allocate(imaske(jpk,jpi))   
       imaske  = huge(imaske(1,1))
!!!$omp parallel default (none) shared(jpk,jpi)
       ! allocate(zpar(jpk,jpi))
       ! zpar    = huge(zpar(1,1))
       allocate(xEPS_ogstm(jpk,jpi))
       !$acc enter data create(xEPS_ogstm)
       !$acc kernels default(present)
       xEPS_ogstm    = huge(xEPS_ogstm(1,1))
       !$acc end kernels
       ! allocate(zpar0m(jpi))
       ! zpar0m  = huge(zpar0m(1))
       ! allocate(zpar100(jpi))
       ! zpar100 = huge(zpar100(1))
!!!$omp end parallel

       allocate(kef(jpj,jpi))
       !$acc enter data create(kef)
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



      subroutine clean_memory_opt

          deallocate(itabe)
          deallocate(imaske)
          ! deallocate(zpar)
          !$acc exit data delete(xEPS_ogstm)
          deallocate(xEPS_ogstm)
          ! deallocate(zpar0m)
          ! deallocate(zpar100)
          !$acc exit data delete(kef)
          deallocate(kef)
          deallocate(kextIO)

      end subroutine clean_memory_opt



      END MODULE 
