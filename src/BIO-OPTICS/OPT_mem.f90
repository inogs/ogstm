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
      double precision, allocatable :: zpar(:,:),xEPS_ogstm(:,:)
      double precision, allocatable :: zpar0m(:),zpar100(:) 
      double precision, allocatable :: kef(:,:)
      double precision, allocatable :: kextIO(:,:,:)
      real, allocatable :: zkef_f (:,:)
  
 
      integer, parameter            :: nwl=33                     
! Radiative transfer model parameter OASIM Native coordinates
      double precision              :: Ed_0m_COARSE(33,12,18,48), Es_0m_COARSE(33,12,18,48) ! lon, lat, day period, wave length
      double precision              :: OASIM_lon(18,48), OASIM_lat(18,48)  
! Radiative transfer model parameter OGSTM coordinates    
      double precision,allocatable  :: Ed_0m(:,:,:), Es_0m(:,:,:) ! lon, lat, day period, wave length
      
      integer                       :: day_RTcheck
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
       allocate(zpar(jpk,jpi))     
       zpar    = huge(zpar(1,1))
       allocate(xEPS_ogstm(jpk,jpi))     
       xEPS_ogstm    = huge(xEPS_ogstm(1,1))
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

! radiative transfer model
       allocate(Ed_0m(nwl,jpj,jpi))
       Ed_0m  =huge(Ed_0m(1,1,1))
       allocate(Es_0m(nwl,jpj,jpi))
       Es_0m  =huge(Es_0m(1,1,1))


#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_OPT



      subroutine clean_memory_opt

          deallocate(itabe)
          deallocate(imaske)
          deallocate(zpar)
          deallocate(xEPS_ogstm)
          deallocate(zpar0m)
          deallocate(zpar100)
          deallocate(kef)
          deallocate(kextIO)

      end subroutine clean_memory_opt



      END MODULE 
