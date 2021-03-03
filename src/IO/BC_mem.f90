       MODULE BC_mem

       USE modul_param
       USE myalloc
       ! epascolo USE myalloc_mpp

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

       use bc_set_mod

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

       IMPLICIT NONE

       public

! ----------------------------------------------------------------------
!  Variables for nutrient Nudging
!  ---------------------------------------------------------------------
      INTEGER jn_gib, jn_riv, jn_atm
      INTEGER(4) Gsizeglo, Rsizeglo, Asizeglo
      INTEGER(4) Gsize, Rsize, Asize, lat, lon

      INTEGER(4), allocatable, DIMENSION(:) ::  tra_matrix_gib, tra_matrix_riv, tra_matrix_atm


      INTEGER(4), allocatable, DIMENSION(:)   ::  gib_idxtglo, riv_idxtglo  ! <--domrea
      INTEGER(4), allocatable, DIMENSION(:,:) ::  gib_ridxt, riv_ridxt,atm_idxtglo   ! <--domrea.ReIndexing()
      INTEGER(4), allocatable, DIMENSION(:,:,:) ::  atm_ridxt
!     queste matrici ridotte hanno nel primo indice l'indice globale, negli altri 3 il corrispondente indice i,j,k locale (per cpu)


      double precision, allocatable, DIMENSION(:)       :: gib_aux,riv_aux
      double precision, allocatable, DIMENSION(:)       :: restocorr
      double precision, allocatable, DIMENSION(:,:)     :: gib,riv,atm_aux
      double precision, allocatable, DIMENSION(:,:,:)   :: gib_dtatrc,riv_dtatrc,atm       ! <--bc_gib.load_gib(), bc_tin.
      double precision, allocatable, DIMENSION(:,:,:,:) :: resto,restotr,atm_dtatrc        !array of restoring coeff. for passive tracers

! ----------------------------------------------------------------------
!  BEGIN BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

      type(bc_set), pointer :: boundaries => null()

! ----------------------------------------------------------------------
!  END BC_REFACTORING SECTION
!  ---------------------------------------------------------------------

! ----------------------------------------------------------------------
      CONTAINS
!     *******************************************************************
!     SUBROUTINE alloc_DTATRC
!     *******************************************************************
      SUBROUTINE alloc_DTATRC
      USE TIME_MANAGER
      ! epascolo USE myalloc_mpp

      INTEGER  :: err
      double precision  :: aux_mem

      CHARACTER(LEN=50) atmfile, rivfile, gibfile

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif
!      CHARACTER(LEN=9)  nomedim0
!      CHARACTER(LEN=13) nomedim1
!
!        nomedim0='riv_idxt'//CHAR(0)
!        nomedim1 = 'gib_idxt_N1p'//CHAR(0)

        atmfile='BC/ATM_'//TC_ATM%TimeStrings(1)//'.nc'!//char(0)
      
        ! gibfile='BC/GIB_'//TC_GIB%TimeStrings(1)//'.nc'!//char(0)
      
        ! rivfile='BC/TIN_'//TC_TIN%TimeStrings(1)//'.nc'!//char(0)
      

         call getDimension(atmfile,'lat',lat)
         call getDimension(atmfile,'lon',lon)
         !print *,lat,lon
         ! call getDimension(gibfile,'gib_idxt_N1p',Gsizeglo)
         ! call getDimension(rivfile,'riv_idxt'    ,Rsizeglo)

!           nomedim1(1:12)='gib_idxt_N1p'
       !CALL ioogsnc_idx(gibfile,nomedim1,Gsizeglo)
!           nomedim0(1:8)='riv_idxt'
       
       !CALL ioogsnc_idx(rivfile,nomedim0,Rsizeglo)
!           nomedim0(1:8)='atm_idxt'
       !CALL ioogsnc_idx(atmfile,nomedim0,Asizeglo)


               ! write(*,*) 'Size of vector Gib to allocate -->', Gsizeglo, ' myrank = ',myrank
           if (lwp) then
               ! write(*,*) 'Size of vector Riv to allocate -->', Rsizeglo
               write(*,*) 'Size of vector Atm to allocate -->', lat,lon
           endif

!     We consider nudging of N1p O2o N3n Slca
!     O2o jn = 1 Dissolved Oxygen (seasonal)
!     N1p jn = 2 phophates        (seasonal)
!     N3n jn = 3 nitrates         (climatological)
!     N5s jn = 4 Silicates        (seasonal)

      jn_gib  = 6
      ! jn_riv  = 6
      jn_atm  = 2

       ! resto is kept just to provide compliance with bfmv2, but should be removed with bfmv5
       allocate(resto(jpk,jpj,jpi,jn_gib))
       resto   = huge(resto(1,1,1,1))
       ! allocate(restotr(jpk,jpj,jpi,jptra))
       ! restotr = huge(restotr(1,1,1,1))

       ! IF (Gsizeglo .NE. 0) THEN

           ! allocate(tra_matrix_gib(jn_gib))
       ! tra_matrix_gib = huge(tra_matrix_gib(1))
           ! allocate(restocorr(jn_gib))
       ! restocorr      = huge(restocorr(1))                  !Correction to restoration for O3c and O3h
       !     allocate(gib_aux(       Gsizeglo))
       ! gib_aux        = huge(gib_aux(1))
       !     allocate(gib_idxtglo(   Gsizeglo))
       ! gib_idxtglo    = huge(gib_idxtglo(1))

          ! tra_matrix_gib(1) = ppO2o
       ! restocorr(1)=1. ! dissolved Oxygen
          ! tra_matrix_gib(2) = ppN1p
       ! restocorr(2)=1. ! phosphates
          ! tra_matrix_gib(3) = ppN3n
       ! restocorr(3)=1. ! nitrates
          ! tra_matrix_gib(4) = ppN5s
       ! restocorr(4)=1. ! silicates
          ! tra_matrix_gib(5) = ppO3c
       ! restocorr(5)=2. ! Dic
          ! tra_matrix_gib(6) = ppO3h
       ! restocorr(6)=2. ! Alk
          ! tra_matrix_gib(7) = ppN6r
       ! restocorr(7)=2. ! N6r
       ! ENDIF

       ! IF (Rsizeglo .NE. 0) THEN
           ! allocate(tra_matrix_riv(jn_riv))
       ! tra_matrix_riv = huge(tra_matrix_riv(1))
       !     allocate(riv_aux(       Rsizeglo))
       ! riv_aux        = huge(riv_aux(1))
       !     allocate(riv_idxtglo(   Rsizeglo))
       ! riv_idxtglo    = huge(riv_idxtglo(1))

          ! tra_matrix_riv(1) = ppN1p ! phosphates
          ! tra_matrix_riv(2) = ppN3n ! nitrates
          ! tra_matrix_riv(3) = ppN5s ! silicates
          ! tra_matrix_riv(4) = ppO3c ! Dic
          ! tra_matrix_riv(5) = ppO3h ! Alk
          ! tra_matrix_riv(6) = ppO2o ! Oxygen
       ! ENDIF

       IF ((lat .NE. 0) .AND. (lon .NE. 0)) THEN
           allocate(tra_matrix_atm(jn_atm))    
       tra_matrix_atm = huge(tra_matrix_atm(1))
           allocate(atm_aux(jpj,jpi))  
       atm_aux        = huge(atm_aux(1,1))
           allocate(atm_idxtglo(   jpj,jpi))  
       atm_idxtglo    = huge(atm_idxtglo(1,1))

          tra_matrix_atm(1) = ppN1p ! phosphates
          tra_matrix_atm(2) = ppN3n ! nitrates
       ENDIF


#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END SUBROUTINE alloc_DTATRC

!     ****************************************
!     SUBROUTINE alloc_DTATRC_local_gib
!     called by domrea
!     ****************************************
      SUBROUTINE alloc_DTATRC_local_gib()
#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif
      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif


       ! if(lwp) write(*,*) 'BC_mem -> Gsize : ', Gsize

       ! allocate(gib_ridxt (4, Gsize        ))
       ! gib_ridxt  = huge(gib_ridxt(1,1))
       ! allocate(gib_dtatrc(Gsize, 2, jn_gib))
       ! gib_dtatrc = huge(gib_dtatrc(1,1,1))
       ! allocate(gib       (Gsize,    jn_gib))
       ! gib        = huge(gib(1,1))


#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END SUBROUTINE alloc_DTATRC_local_gib



!     ****************************************
!     SUBROUTINE alloc_DTATRC_local_riv
!     called by domrea
!     ****************************************
      SUBROUTINE alloc_DTATRC_local_riv()

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif
      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif


       ! allocate(riv_ridxt (4, Rsize       ))
       ! riv_ridxt  = huge(riv_ridxt(1,1))
       ! allocate(riv_dtatrc(2,Rsize,jn_riv))
       ! riv_dtatrc = huge(riv_dtatrc(1,1,1))
       ! allocate(riv       (Rsize,   jn_riv))
       ! riv        = huge(riv(1,1))


#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

      END SUBROUTINE alloc_DTATRC_local_riv


!     ****************************************
!     SUBROUTINE alloc_DTATRC_local_atm
!     called by domrea
!     ****************************************
      SUBROUTINE alloc_DTATRC_local_atm()
#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif
      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

      allocate(atm_dtatrc(jpj,jpi, 2, jn_atm)) 
       atm_dtatrc = huge(atm_dtatrc(1,1,1,1))
      allocate(atm       (jpj,jpi,    jn_atm)) 
       atm        = huge(atm(1,1,1))


#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif
      END SUBROUTINE alloc_DTATRC_local_atm



      subroutine clean_memory_bc()

          ! resto is kept just to provide compliance with bfmv2, but should be removed with bfmv5
          deallocate(resto)

          if ((lat /= 0) .and. (lon /= 0)) then
              deallocate(tra_matrix_atm)
              deallocate(atm_aux)
              deallocate(atm_idxtglo)
          endif

          deallocate(atm_dtatrc)
          deallocate(atm)

      end subroutine clean_memory_bc



      END MODULE
