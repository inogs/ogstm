       MODULE BC_mem

       USE modul_param 
       USE myalloc

       IMPLICIT NONE

       public

! ----------------------------------------------------------------------
!  Variables for nutrient Nudging
!  ---------------------------------------------------------------------
      INTEGER jn_gib, jn_riv, jn_atm
      INTEGER(4) Gsizeglo, Rsizeglo, Asizeglo
      INTEGER(4) Gsize, Rsize, Asize


      REAL(4), allocatable,DIMENSION(:,:,:)   :: ztra_nudg00_f
      REAL(4), allocatable,DIMENSION(:,:,:,:) :: ztra_nudg01_f

      INTEGER(4), allocatable, DIMENSION(:)   ::  gib_idxtglo, riv_idxtglo, atm_idxtglo ! <--domrea
!     INTEGER(4), allocatable, DIMENSION(:)   ::  gib_idxt   , riv_idxt   , atm_idxt
      INTEGER(4), allocatable, DIMENSION(:,:) ::  gib_ridxt  , riv_ridxt  , atm_ridxt   ! <--domrea.ReIndexing()
!     queste matrici ridotte hanno nel primo indice l'indice globale, negli altri 3 il corrispondente indice i,j,k locale (per cpu)


      REAL(8), allocatable, DIMENSION(:)     :: gib_aux,riv_aux,atm_aux
      REAL(8), allocatable, DIMENSION(:,:)   :: gib,riv, atm
      REAL(8), allocatable, DIMENSION(:,:,:) :: gib_dtatrc,riv_dtatrc,atm_dtatrc       ! <--bc_gib.load_gib(), bc_tin.

! ----------------------------------------------------------------------
      CONTAINS
!     *******************************************************************
!     SUBROUTINE alloc_DTATRC
!     called by myOPA()
!     *******************************************************************
      SUBROUTINE alloc_DTATRC
      USE TIME_MANAGER

      IMPLICIT NONE
      CHARACTER*50 atmfile, rivfile, gibfile
      CHARACTER*9 nomedim

      INTEGER jn_nudg

        nomedim='12345678'
        nomedim(9:9) = CHAR(0)

        atmfile='BC/ATM_'//TC_ATM%TimeStrings(1)//'.nc'//char(0);
        gibfile='BC/GIB_'//TC_GIB%TimeStrings(1)//'.nc'//char(0);
        rivfile='BC/TIN_'//TC_TIN%TimeStrings(1)//'.nc'//char(0);


           nomedim(1:8)='gib_idxt'; CALL ioogsnc_idx(gibfile,nomedim,Gsizeglo)
           nomedim(1:8)='riv_idxt'; CALL ioogsnc_idx(rivfile,nomedim,Rsizeglo)
           nomedim(1:8)='atm_idxt'; CALL ioogsnc_idx(atmfile,nomedim,Asizeglo)

           if (lwp) then
               write(*,*) 'Size of vector Gib to allocate -->', Gsizeglo
               write(*,*) 'Size of vector Riv to allocate -->', Rsizeglo
               write(*,*) 'Size of vector Atm to allocate -->', Asizeglo
           endif

!     We consider nudging of N1p O2o N3n Slca
!     O2o jn = 1 Dissolved Oxygen (seasonal)
!     N1p jn = 2 phophates        (seasonal)
!     N3n jn = 3 nitrates         (climatological)
!     N5s jn = 4 Silicates        (seasonal)
      
      jn_nudg = 4
      jn_gib  = 4
      jn_riv  = 4
      jn_atm  = 4

       allocate(ztra_nudg00_f(jpiglo,jpjglo,jpk        ))
       allocate(ztra_nudg01_f(jpiglo,jpjglo,jpk,jn_nudg)) 
       
       IF (Gsizeglo .NE. 0) THEN
           allocate(gib_aux(       Gsizeglo))
           allocate(gib_idxtglo(   Gsizeglo))
           mem_all = 4*2*Gsizeglo
       ENDIF
       
       IF (Rsizeglo .NE. 0) THEN
           allocate(riv_aux(       Rsizeglo))
           allocate(riv_idxtglo(   Rsizeglo))
           mem_all = 4*2*Rsizeglo
       ENDIF
       
       IF (Asizeglo .NE. 0) THEN
           allocate(atm_aux(       Asizeglo))
           allocate(atm_idxtglo(   Asizeglo))
           mem_all = 4*2*Asizeglo
       ENDIF

      END SUBROUTINE alloc_DTATRC

!     ****************************************
!     SUBROUTINE alloc_DTATRC_local_gib
!     called by domrea
!     ****************************************
      SUBROUTINE alloc_DTATRC_local_gib()

       if(lwp) write(*,*) 'BC_mem -> Gsize : ', Gsize

       allocate(gib_ridxt (Gsize, 4        ))
       allocate(gib_dtatrc(Gsize, 2, jn_gib))
       allocate(gib       (Gsize,    jn_gib))

      END SUBROUTINE alloc_DTATRC_local_gib



!     ****************************************
!     SUBROUTINE alloc_DTATRC_local_riv
!     called by domrea
!     ****************************************
      SUBROUTINE alloc_DTATRC_local_riv()

       allocate(riv_ridxt (Rsize,4        ))
       allocate(riv_dtatrc(Rsize,2, jn_riv))
       allocate(riv       (Rsize,   jn_riv))

      END SUBROUTINE alloc_DTATRC_local_riv


!     ****************************************
!     SUBROUTINE alloc_DTATRC_local_atm
!     called by domrea
!     ****************************************
      SUBROUTINE alloc_DTATRC_local_atm()

       allocate(atm_ridxt (Asize, 4        ))
       allocate(atm_dtatrc(Asize, 2, jn_atm))
       allocate(atm       (Asize,    jn_atm))

      END SUBROUTINE alloc_DTATRC_local_atm

      END MODULE 
