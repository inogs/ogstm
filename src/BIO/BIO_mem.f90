       MODULE BIO_mem 

       USE modul_param 
       USE myalloc
       USE TIME_MANAGER

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif


       IMPLICIT NONE

       public

      double precision, allocatable :: ogstm_sedipi(:,:,:,:)
      double precision, allocatable :: ogstm_ph(:,:,:) ! GUESS for FOLLOWS algorithm
      double precision, allocatable :: NPPF2(:,:,:)
      double precision, allocatable :: ogstm_co2(:,:), co2_IO(:,:,:)
      double precision, allocatable :: sediPPY(:,:)
      double precision, allocatable :: local_D3DIAGNOS(:,:)
      double precision, allocatable :: local_D2DIAGNOS(:,:)
      double precision, allocatable :: er(:,:)
      double precision:: ice


!!!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_BIO()

      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

       allocate(ogstm_co2(jpj,jpi))          
       ogstm_co2       = huge(ogstm_co2(1,1))
       allocate(co2_IO(jpj,jpi,2))
       !$acc enter data create(co2_IO)
        
       co2_IO    = huge(co2_IO(1,1,1))
       allocate(ogstm_sedipi(jpk,jpj,jpi,4))
       ogstm_sedipi    = huge(ogstm_sedipi(1,1,1,1))
       allocate(ogstm_ph(jpk,jpj,jpi))       
       ogstm_ph        = huge(ogstm_ph(1,1,1))
       ogstm_ph=8.0
       allocate(NPPF2(jpk,jpj,jpi))
       NPPF2 = 0 ! nut huge, because it will be assigned only in trcBIO in BFMpoints
                 ! and used in hard_tissue_pump.F also in land points
       ice=0

       allocate(sediPPY(jpi * jpj * jpk, 4))
       allocate(local_D3DIAGNOS(jpi * jpj * jpk, jptra_dia))
       allocate(local_D2DIAGNOS(jpi * jpj, jptra_dia_2d))
       allocate(er(jpi * jpj * jpk, 11))
       !$acc enter data create(ogstm_co2,ogstm_sedipi,ogstm_ph,sediPPY,local_D3DIAGNOS,local_D2DIAGNOS,er)

#ifdef Mem_Monitor
       mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_BIO



      subroutine clean_memory_bio()

            deallocate(ogstm_co2)
            !$acc exit data delete(co2_IO)
            deallocate(co2_IO)
            deallocate(ogstm_sedipi)
            deallocate(ogstm_ph)
            deallocate(NPPF2)
            deallocate(sediPPY)
            deallocate(local_D3DIAGNOS)
            deallocate(local_D2DIAGNOS)
            deallocate(er)
            !$acc exit data delete(ogstm_co2,ogstm_sedipi,ogstm_ph,sediPPY,local_D3DIAGNOS,local_D2DIAGNOS,er)

      end subroutine clean_memory_bio



      END MODULE 
