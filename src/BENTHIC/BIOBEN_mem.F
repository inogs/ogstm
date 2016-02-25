       MODULE BIOBEN_mem 

       USE modul_param 
       USE myalloc
       USE TIME_MANAGER

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif


       IMPLICIT NONE

       public

      REAL(8), allocatable :: trn_ben(:,:,:,:), tra_ben(:,:,:,:),trb_ben(:,:,:,:)
      REAL(8), allocatable :: e3t_ben(:,:,:), e3w_ben(:,:,:)


CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_BIOBEN()

      INTEGER  :: err
      REAL(8)  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

       allocate(trb_ben(jpi,jpj,jpk_ben,jptra_b))        ; trb_ben   = huge(trb_ben(1,1,1,1))
       allocate(trn_ben(jpi,jpj,jpk_ben,jptra_b))        ; trn_ben   = huge(trn_ben(1,1,1,1))
       allocate(tra_ben(jpi,jpj,jpk_ben,jptra_b))        ; tra_ben   = huge(tra_ben(1,1,1,1))

       allocate(e3t_ben(jpi,jpj,jpk_ben))                ; e3t_ben   = huge(e3t_ben(1,1,1))
       allocate(e3w_ben(jpi,jpj,jpk_ben))                ; e3w_ben   = huge(e3w_ben(1,1,1))


#ifdef Mem_Monitor
       mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_BIOBEN

      END MODULE 
