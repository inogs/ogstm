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
      INTEGER, allocatable :: calcbenmask(:,:)

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
      INTEGER           :: jptra_high_ben,jptra_dia_high_ben,jptra_dia2d_high_ben
      INTEGER, allocatable, dimension(:) :: highfreq_table_ben,highfreq_table_dia_ben
      INTEGER, allocatable, dimension(:) :: highfreq_table_dia2d_ben

! IO
      REAL(8), allocatable :: totbenIO(:,:,:), benIO_HIGH(:,:,:,:),benIO(:,:,:,:)
      REAL(8), allocatable :: tottrn_ben(:,:,:),tottrb_ben(:,:,:)
      REAL(8), allocatable :: buffben(:),buffben3D(:,:,:)
      REAL(8), allocatable :: bufftrn_ben(:),bufftrb_ben(:)
      REAL(8), allocatable :: ben_DIA_IO(:,:,:,:),  ben_DIA_IO_HIGH(:,:,:,:)
      REAL(8), allocatable :: ben_DIA(:,:,:,:),  benDIA_HIGH(:,:,:,:)
      REAL(8), allocatable :: ben_DIA_2d(:,:,:)
      REAL(8), allocatable :: ben_DIA_2d_IO(:,:,:)
      REAL(8), allocatable :: ben_DIA_2d_IO_HIGH(:,:,:)
      REAL(8), allocatable :: buffben2D(:),buffben_dia(:)
      REAL(8), allocatable :: totbenIO2d(:,:)
      REAL(4), allocatable :: d2f2d_ben(:,:),d2f3d_ben(:,:,:)


CC----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_BIOBEN()

      INTEGER  :: err
      REAL(8)  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif
       call getDimension('benthic_mask.nc','z_ben',jpk_ben)

       allocate(trb_ben(jpi,jpj,jpk_ben,jptra_b))        ; trb_ben    = huge(trb_ben(1,1,1,1))
       allocate(trn_ben(jpi,jpj,jpk_ben,jptra_b))        ; trn_ben    = huge(trn_ben(1,1,1,1))
       allocate(tra_ben(jpi,jpj,jpk_ben,jptra_b))        ; tra_ben    = huge(tra_ben(1,1,1,1))

       allocate(e3t_ben(jpi,jpj,jpk_ben))                ; e3t_ben    = huge(e3t_ben(1,1,1))
       allocate(e3w_ben(jpi,jpj,jpk_ben))                ; e3w_ben    = huge(e3w_ben(1,1,1))
       allocate(calcbenmask(jpi,jpj))                    ; calcbenmask= huge(calcbenmask(1,1))


       allocate(gdept_ben(jpk_ben))                      ; gdept_ben  = huge(gdept_ben(1))
       allocate(gdepw_ben(jpk_ben))                      ; gdepw_ben  = huge(gdepw_ben(1))

      if (lwp)  allocate(totbenIO(jpiglo,jpjglo,jpk_ben))               ; totbenIO   = huge(totbenIO(1,1,1))

       allocate(benIO(jpi,jpj,jpk_ben,jptra_b))          ; benIO      = huge(benIO(1,1,1,1))
       allocate(benIO_HIGH(jpi,jpj,jpk_ben,jptra_high_ben)); benIO_HIGH = huge(benIO_HIGH(1,1,1,1))

       allocate(bufftrn_ben(jpi_max* jpj_max* jpk_ben))      ; bufftrn_ben    = huge(bufftrn_ben(1))
       allocate(bufftrb_ben(jpi_max* jpj_max* jpk_ben))      ; bufftrb_ben    = huge(bufftrb_ben(1))

       if (lwp) allocate(tottrn_ben(jpiglo, jpjglo, jpk_ben))  ;  if (lwp) tottrn_ben = huge(tottrn_ben(1,1,1)) 
       if (lwp) allocate(tottrb_ben(jpiglo, jpjglo, jpk_ben))  ;  if (lwp) tottrb_ben = huge(tottrb_ben(1,1,1))

       allocate(ben_DIA(jpi,jpj,jpk_ben,jptra_b))        ; ben_DIA  = huge(ben_DIA(1,1,1,1))
       allocate(ben_DIA_IO(jpi,jpj,jpk_ben,jptra_b))     ; ben_DIA_IO = huge(ben_DIA_IO(1,1,1,1))
       allocate(ben_DIA_IO_HIGH(jpi,jpj,jpk_ben,jptra_high_ben)); ben_DIA_IO_HIGH = huge(ben_DIA_IO_HIGH(1,1,1,1))

       allocate(ben_DIA_2d(jpi,jpj,jptra_dia_b_2d))     ; ben_DIA_2d    = huge(ben_DIA_2d(1,1,1))
       allocate(ben_DIA_2d_IO(jpi,jpj,jptra_dia_b_2d))  ; ben_DIA_2d_IO = huge(ben_DIA_2d_IO(1,1,1))
       allocate(ben_DIA_2d_IO_HIGH(jpi,jpj,jptra_dia_b_2d))  ; ben_DIA_2d_IO_HIGH = huge(ben_DIA_2d_IO_HIGH(1,1,1))
       if (lwp) allocate(totbenIO2d(jpiglo,jpjglo))                    ; totbenIO2d  = huge(totbenIO2d(1,1))

       allocate(buffben_dia(jpi_max*jpj_max*jpk_ben)) ; buffben_dia= huge(buffben_dia(1))
       allocate(buffben(jpi_max* jpj_max* jpk_ben))               ; buffben    = huge(buffben(1))
       allocate(buffben2D(jpi_max* jpj_max))                      ; buffben2D  = huge(buffben2D(1))
       allocate(buffben3D(jpi, jpj, jpk_ben))                     ; buffben3D  = huge(buffben3D(1,1,1))

       if (lwp) allocate(d2f3d_ben  (jpiglo,jpjglo,jpk_ben))      ; d2f3d_ben  = huge(d2f3d_ben(1,1,1))
       if (lwp) allocate(d2f2d_ben  (jpiglo,jpjglo))              ; d2f2d_ben  = huge(d2f2d_ben(1,1))

#ifdef Mem_Monitor
       mem_all=get_mem(err) - aux_mem
#endif

      END subroutine myalloc_BIOBEN

      END MODULE 
