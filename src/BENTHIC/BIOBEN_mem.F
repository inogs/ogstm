       MODULE BIOBEN_mem 

       USE modul_param 
       USE myalloc
       USE IO_mem
       USE TIME_MANAGER

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif


       IMPLICIT NONE

       public

      INTEGER              :: jpk_ben
      REAL(8), allocatable :: trn_ben(:,:,:,:), tra_ben(:,:,:,:),trb_ben(:,:,:,:)
      REAL(8), allocatable :: e3t_ben(:,:,:),   e3w_ben(:,:,:) 
! Domain
      REAL(8), allocatable :: gdept_ben(:), gdepw_ben(:)

! Namelist parameters

      CHARACTER(LEN=20) :: cbennm(jptra_b)
      CHARACTER(LEN=12) :: cbenun(jptra_b)
      INTEGER           :: cbe_hf(jptra_b)
      CHARACTER(LEN=20) :: dibnm(jptra_dia_b)
      CHARACTER(LEN=20) :: dibun(jptra_dia_b)
      INTEGER           :: dibhf(jptra_dia_b)
      REAL(8)           :: cbemax(jptra_b)
      CHARACTER(LEN=20) :: dibnm_2d(jptra_dia_b_2d)
      CHARACTER(LEN=20) :: dibun_2d(jptra_dia_b_2d)
      INTEGER           :: dibhf_2d(jptra_dia_b_2d)

! IO
      REAL(8), allocatable :: totbenIO(:,:,:), benIO_HIGH(:,:,:,:),benIO(:,:,:,:),buffben(:)
      REAL(4), allocatable :: d2f3d_ben(:,:,:)


CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_BIOBEN()

      INTEGER  :: err
      REAL(8)  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif
       call getDimension('meshmask.nc','z_ben',jpk_ben)

       allocate(trb_ben(jpi,jpj,jpk_ben,jptra_b))        ; trb_ben    = huge(trb_ben(1,1,1,1))
       allocate(trn_ben(jpi,jpj,jpk_ben,jptra_b))        ; trn_ben    = huge(trn_ben(1,1,1,1))
       allocate(tra_ben(jpi,jpj,jpk_ben,jptra_b))        ; tra_ben    = huge(tra_ben(1,1,1,1))

       allocate(e3t_ben(jpi,jpj,jpk_ben))                ; e3t_ben    = huge(e3t_ben(1,1,1))
       allocate(e3w_ben(jpi,jpj,jpk_ben))                ; e3w_ben    = huge(e3w_ben(1,1,1))

       allocate(gdept_ben(jpk_ben))                      ; gdept_ben  = huge(gdept_ben(1))
       allocate(gdepw_ben(jpk_ben))                      ; gdepw_ben  = huge(gdepw_ben(1))

       allocate(totbenIO(jpi,jpj,jpk_ben))               ; totbenIO   = huge(totbenIO(1,1,1))

       allocate(benIO_HIGH(jpi,jpj,jpk_ben,jptra_b))     ; benIO_HIGH = huge(benIO_HIGH(1,1,1,1))
       allocate(benIO(jpi,jpj,jpk_ben,jptra_b))          ; benIO      = huge(benIO(1,1,1,1))

       allocate(buffben(jpi_max* jpj_max* jpk_ben))      ; buffben    = huge(buffben(1))

       allocate(d2f3d_ben  (jpiglo,jpjglo,jpk_ben))      ; d2f3d_ben  = huge(d2f3d_ben(1,1,1))

#ifdef Mem_Monitor
       mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_BIOBEN

      END MODULE 
