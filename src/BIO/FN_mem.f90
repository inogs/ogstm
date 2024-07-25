       MODULE FN_mem 

       USE modul_param 
       USE myalloc

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif

       IMPLICIT NONE

       public

!----------------------------------------------------------------------
! Common/comcoh/  : False Negative matrix
! ---------------------------------------------------------------------
      INTEGER              :: dimen_jvsnu
      INTEGER, allocatable :: jarr_snu(:,:)
      double precision              :: SMALL
      double precision, allocatable :: tra_FN(:,:,:,:)
      double precision, allocatable :: TOT(:,:),TOT_FN(:,:),FN_CORR(:,:)
      double precision, allocatable :: TOTcalc(:)

  
      double precision, allocatable :: FN_ranking(:)
      INTEGER :: elements,nelements(6),idx_element(14,6)


!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_FN()


      IMPLICIT NONE

      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

      dimen_jvsnu = 0
      SMALL=0.00000000001

      allocate(jarr_snu(2, jpi*jpj))      
      jarr_snu = huge(jarr_snu(1,1))
      allocate(tra_FN(jpk,jpj,jpi,jptra))
      !$acc enter data create(tra_FN)
      !$acc kernels default(present)
      tra_FN   = huge(tra_FN(1,1,1,1))
      !$acc end kernels

      !CALL OPA_elements(elements,nelements,idx_element)

      !allocate(TOT    (jpj*jpi,elements)) 
      ! TOT        = huge(TOT(1,1))        ! elements is the number of element considered (C,P,N,...)
      !allocate(TOT_FN (jpj*jpi,elements)) 
      ! TOT_FN     = huge(TOT_FN(1,1))
      !allocate(FN_CORR(jpj*jpi,elements)) 
      ! FN_CORR    = huge(FN_CORR(1,1))

      allocate(TOTcalc(jpj*jpi))
      TOTcalc    =huge(TOTcalc(1))
     
      allocate(FN_ranking(jptra))         
      FN_ranking = huge(FN_ranking(1))


      !$acc kernels default(present)
      tra_FN=0.
      !$acc end kernels
     
!     cor_FN=0.
      FN_ranking=0.

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_FN



      subroutine clean_memory_fn

          !$acc exit data delete(tra_FN)
          deallocate(jarr_snu)
          deallocate(tra_FN)
          deallocate(TOTcalc)
          deallocate(FN_ranking)

      end subroutine clean_memory_fn



      END MODULE FN_mem
