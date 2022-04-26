       MODULE myalloc

       USE modul_param
       USE timers

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif

       IMPLICIT NONE

       public

!!----------------------------------------------------------------------
!!            domain parameters
!! ---------------------------------------------------------------------
!!      nizoom, njzoom   : (i,j) indexes for the domain left bottom zoom
!!      nperio           : lateral boundary condition type 
!!      nimpp,njmpp      : (i,j) indexes for mpp-subdomain left bottom
!!      nreci,nrecj      : overlap region in i and j
!!      nproc            : number for local processor
!!      narea            : number for local area
!!      nbondi, nbondj   : mark of i- and j-direction local boundaries
!!      nlci, nlcj       : i, j dimensions of the local subdomain
!!      nldi, nlei,      : first and last indoor i- and j-indexes
!!      nldj, nlej   
!!      noea, nowe,      : index of the local neighboring processors in
!!      noso, nono         east, west, south and north directions
!!      nimppt,njmppt(): i-, j-indexes for each processor
!!      nlcit, nlcjt() : dimensions of every subdomain
!!      nldit, nldjt() : first, last indoor index for each i-domain
!!      nleit, nlejt() : first, last indoor index for each j-domain
!!      mindi, mindj() : indexes array of the subdomain

      INTEGER, PARAMETER :: nizoom=1,  njzoom=1
      INTEGER :: mpi_glcomm_size,myrank
      INTEGER nimpp,njmpp
      INTEGER nperio, narea, nlci, nlcj
      INTEGER nbondi, nbondj, nproc, noea, nowe, noso, nono
      INTEGER nreci, nrecj, nldi, nlei, nldj, nlej
      INTEGER, allocatable :: ilcit(:,:), ilcjt(:,:)
      INTEGER, allocatable :: mindi(:), mindj(:)
      INTEGER, allocatable :: nimppt(:), njmppt(:), nlcit(:), nlcjt(:)
      INTEGER, allocatable ::  nldit(:),  nldjt(:), nleit(:), nlejt(:)
      double precision  mem_all

      INTEGER npolj

      double precision, PARAMETER ::  g =9.80665  ! gravity





!!----------------------------------------------------------------------
!!       ocean physical parameters (equation of state, ...)
!! ------------------------------------------
!!        neos         : flag of the type of equation of state used
!!      rau0             : reference volumic mass of the ocean (kg/m3)
!!      ralpha, rbeta    : thermique and haline expension coef. used
!!               for linear equation of state (neos=1 or 2)
      double precision Euphotic_lev
      INTEGER jpk_eu
      LOGICAL forcing_phys_initialized
      LOGICAL read_W_from_file, internal_sponging, ingv_files_direct_reading
      INTEGER ingv_lon_shift
      INTEGER    deflate_ave, deflate_level_ave,deflate_rst, deflate_level_rst
      INTEGER neos
      double precision rau0, ralpha, rbeta
      double precision  rdt     ! dynamics time step
      integer imposed_deltaT(2)
      logical variable_rdt

!!----------------------------------------------------------------------
!!        horizontal curvilinear coordinate and scale factors
!! ---------------------------------------------------------------------
!!      glamt          : longitude of t-point (degre)
!!      glamu          : longitude of u-point (degre)
!!      glamv          : longitude of v-point (degre)
!!      glamf          : longitude of f-point (degre)
!!      gphit          : latitude  of t-point (degre)
!!      gphiu          : latitude  of u-point (degre)
!!      gphiv          : latitude  of v-point (degre)
!!      gphif          : latitude  of f-point (degre)
!!      e1t,e2t        : horizontal scale factors at t-point (m)
!!      e1u,e2u        : horizontal scale factors at u-point (m)
!!      e1v,e2v        : horizontal scale factors at v-point (m)
!!      e1f,e2f        : horizontal scale factors at f-point (m)
!!        ff             : coriolis factor


      double precision, allocatable, dimension(:,:) :: totglamt, glamu, glamv,glamf , glamt
      double precision, allocatable, dimension(:,:) :: totgphit, gphiu, gphiv,gphif , gphit
      double precision, allocatable, dimension(:,:) :: e1t, e1u, e1v, e1f
      double precision, allocatable, dimension(:,:) :: e2t, e2u, e2v, e2f, ff

!!----------------------------------------------------------------------
!!       vertical coordinate and scale factors
!! -------------------------------------------------------

!!                  z-coordinate (default option)
!!                  ------------------------------
!!      gdept, gdepw() : depth of t- and w-points (m)
!!      e3t_0, e3w_0()     : vertical scale factors at t- and w-points (m)
!!
      !dir$ attributes align:64 :: e3t

#ifdef gdept1d
      double precision, allocatable :: gdept(:), gdepw(:)
#else
      double precision, allocatable :: gdept(:,:,:), gdepw(:)
#endif

      double precision, allocatable,dimension(:,:,:), save :: e3t, e3t_back, e3u, e3v, e3w
      double precision, allocatable,dimension(:,:,:), save :: e3t_0, e3u_0, e3v_0, e3w_0
      double precision, allocatable :: spongeT(:,:) , spongeVel(:,:,:)

!!----------------------------------------------------------------------
!!        masks, bathymetry
!! -----------------------------------
!!      mbathy         : number of ocean level (=0, 1, ... , jpk-1)
!!      tmask, umask() : land/ocean mask at t-, u-, v- and f-points
!!      vmask, fmask()

      INTEGER, allocatable :: mbathy(:,:)
      double precision, allocatable, dimension(:,:) :: h_column


      INTEGER(kind = 1), allocatable, dimension(:,:,:) :: tmask,umask, vmask

      ! logical bfmmask is derived from bfmmask.nc
      ! it is supposed to map those points for which bfmv5 has to be resolved
      integer(1), allocatable, dimension(:, :, :) :: bfmmask

      INTEGER NBFMPOINTS, NBFMPOINTS_SUP, NWATERPOINTS
      INTEGER, allocatable, dimension(:,:) :: BFMpoints


!! II. DYNAMICS AND TRACERS
!! ========================
!!----------------------------------------------------------------------
!!       previous fields (before)
!! -----------------------------------------

      


!!----------------------------------------------------------------------
!!      present fields (now)
!! -------------------------------------
!!       un, vn(), wn() : horizontal and vertical velocity (m s-1)
!!        tn,   sn()     : pot. temperature (celsius), salinity (psu)
!!      rdn            : in situ density anomalie rdn=(rho-rau0)/rau0
!!                         (no units)
!!        rhopn          : potential volumic mass (kg m-3)
!!      bn2n           : brunt-vaisala frequency (s-2)
!!
      double precision, allocatable, dimension(:,:,:) :: un, vn, wn
      double precision, allocatable, dimension(:,:,:) :: tn, sn,rdn,rhopn,rho



!! III. OCEAN PHYSICS
!! ==================


!!----------------------------------------------------------------------
!!      lateral diffusivity (tracers)
!! ------------------------------------
!!      aht0             : lateral diffusivity coefficient    (namelist)
!!      ahtu, ahtv()   : lateral diffusivity coef. at u-, v-, w- t-pts
!!      ahtw, ahtt()     (harmonic operator: no rotation, use of u-
!!                          and v-points rotation, use of u-, v- w-pts)
!!                         (biharmonic operator: rotation or not, use of
!!                          t-point only)
!!                         (the arrays used are 3D, 2D, 1D or 0D depen-
!!                          ding on 'key_trahdfcoef.d' )

      double precision aht0
      double precision, allocatable :: ahtu(:), ahtv(:), ahtw(:), ahtt(:)


!!----------------------------------------------------------------------
!!         vertical diffusion
!! ----------------------------------------------------------------------
!!     avt             : vertical diffusivity coeff. at w-point
!!     avtb            : background profile of avm and avt

      double precision, allocatable :: avt(:,:,:)
      double precision, allocatable ::  avtb(:)


!! IV. SURFACE FORCING AND DATA
!! ============================
!!  surface wind stress at givem time_step
!!    taux, tauy()   : wind stress components in (i,j) referential
      double precision, allocatable, dimension(:,:) :: taux, tauy, vatm, freeze

!!----------------------------------------------------------------------
!!     surface fluxes
!! -------------------------------
!!      qt             : total surface heat flux (w m-2)
!!      q              : surface heat flux (w m-2)
!!      emp            : evaporation minus precipitation (mm day-1)
!!      runoff         : annual run off (mm/day)

      double precision, allocatable, dimension(:,:) :: qt, q, emp,runoff
      double precision, allocatable :: qsr(:,:) ! penetrative solar radiation (w m-2)

      INTEGER nsptint ! YPE of spatial interpolation (NAMELIST)


!!      udta,vdta()  : horizontal velocity data array
!!      wdta         : vertical velocity data array
!!      avtdta       : avt data array
!!      flxdta         : additional fluxes

      double precision, allocatable, dimension(:,:,:,:) :: udta,vdta,wdta,avtdta,flxdta
      double precision, allocatable, dimension(:,:,:)   :: flx
      double precision, allocatable, dimension(:,:,:,:) :: tdta,sdta ! : temperature and salinity data array
      double precision, allocatable, dimension(:,:,:,:) :: e3tdta,e3udta,e3vdta,e3wdta ! : temperature and salinity data array

!! V. DIAGNOSTICS
!! ==============


      INTEGER calendarType ! leap years calendar (0/1)

      double precision, ALLOCATABLE, DIMENSION(:,:) :: DAY_LENGTH

      INTEGER            :: numnam = 208 ! unit for namelist
      INTEGER, PARAMETER :: numout = 2   ! unit for output print
      LOGICAL lwp                        ! boolean term for stdout
      INTEGER, PARAMETER :: numnat =80   ! the number of the passive tracer NAMELIST
      INTEGER, PARAMETER :: numphys=3   ! the number for physycs tracers NAMELIST_PHYS

#if defined key_mpp
      INTEGER :: EAST_count_send, WEST_count_send, SOUTH_count_send, NORTH_COUNT_send
      INTEGER :: EAST_count_recv, WEST_count_recv, SOUTH_count_recv, NORTH_COUNT_recv
      double precision, allocatable, dimension(:) :: te_send, tw_send, tn_send, ts_send
      double precision, allocatable, dimension(:) :: te_recv, tw_recv, tn_recv, ts_recv
      INTEGER, allocatable, dimension(:,:) :: EASTpoints_send, WESTpoints_send,NORTHpoints_send, SOUTHpoints_send
      INTEGER, allocatable, dimension(:,:) :: EASTpoints_recv, WESTpoints_recv,NORTHpoints_recv, SOUTHpoints_recv
#endif



!! PASSIVE TRACER MODEL
          CHARACTER(LEN=20) :: ctrcnm(jptra)
          CHARACTER(LEN=12) :: ctrcun(jptra)
          CHARACTER(LEN=20) :: dianm(jptra_dia)
          CHARACTER(LEN=20) :: diaun(jptra_dia)
          INTEGER           :: diahf(jptra_dia)
          INTEGER           :: diaWR(jptra_dia)
          CHARACTER(LEN=20) :: dianm_2d(jptra_dia_2d)
          CHARACTER(LEN=20) :: diaun_2d(jptra_dia_2d)
          INTEGER           :: diahf_2d(jptra_dia_2d)
          INTEGER           :: diaWR_2d(jptra_dia_2d)
          CHARACTER(LEN=17) :: COMMON_DATESTRING
!physical tracers
      INTEGER :: jptra_phys, jptra_phys_2d
      INTEGER :: freq_ave_phys
      CHARACTER(LEN=20), allocatable, dimension(:) :: physnm, physun, physnm_2d, physun_2d
      INTEGER, allocatable, dimension(:) :: physWR,physWR_2d



!!    parameters for the control of passive tracers



      double precision ::  ctrmax(jptra)
      LOGICAL :: isCheckLOG
      LOGICAL :: save_bkp_group2 ! we can avoid to dump bkp of a lot of variables
      INTEGER :: jptra_high, jptra_dia_high, jptra_dia2d_high
      INTEGER :: ctr_hf(jptra)

      INTEGER ave_freq_phys, freq_flux_dump



      INTEGER flagSMS_Dyn                    ! Flag time advance SMS or Dyn
      double precision, allocatable ::  trn(:,:,:,:)
      double precision, allocatable ::  tra(:,:,:,:)
      double precision, allocatable ::  tra_DIA(:,:,:,:)
      double precision, allocatable ::  tra_DIA_2d(:,:,:)
      double precision, allocatable ::  traIO(:,:,:,:)
      double precision, allocatable ::  traIO_HIGH(:,:,:,:)
      double precision, allocatable ::  snIO(:,:,:) 
      double precision, allocatable ::  tnIO(:,:,:) 
      double precision, allocatable ::  vatmIO(:,:) 
      double precision, allocatable ::  empIO(:,:) 
      double precision, allocatable ::  qsrIO(:,:) 
      double precision, allocatable ::  unIO(:,:,:) 
      double precision, allocatable ::  vnIO(:,:,:) 
      double precision, allocatable ::  wnIO(:,:,:) 
      double precision, allocatable ::  avtIO(:,:,:) 
      double precision, allocatable ::  e3tIO(:,:,:) 
      double precision, allocatable ::  tra_DIA_IO(:,:,:,:)
      double precision, allocatable ::  tra_DIA_IO_HIGH(:,:,:,:)
      double precision, allocatable ::  tra_DIA_2d_IO(:,:,:)
      double precision, allocatable ::  tra_DIA_2d_IO_HIGH(:,:,:)
      double precision, allocatable ::  tra_PHYS_IO(:,:,:,:)
      double precision, allocatable ::  tra_PHYS_IO_HIGH(:,:,:,:)
      double precision, allocatable ::  tra_PHYS_2d_IO(:,:,:)
      double precision, allocatable ::  tra_PHYS_2d_IO_HIGH(:,:,:)
      double precision, allocatable :: tottrn(:,:,:)
      double precision, allocatable :: tottrb(:,:,:)

!      double precision, allocatable ::  tottrnIO(:,:,:) ! matrix for i/o writing(trcdit.F)
      double precision, allocatable ::  tottrnIO2d(:,:)
      double precision, allocatable ::  tottrbIO(:,:,:)
      double precision, allocatable ::  totsnIO(:,:,:) 
      double precision, allocatable ::  tottnIO(:,:,:) 
      double precision, allocatable ::  totvatmIO(:,:) 
      double precision, allocatable ::  totempIO(:,:) 
      double precision, allocatable ::  totqsrIO(:,:) 
      double precision, allocatable ::  totunIO(:,:,:) 
      double precision, allocatable ::  totvnIO(:,:,:) 
      double precision, allocatable ::  totwnIO(:,:,:) 
      double precision, allocatable ::  totavtIO(:,:,:) 
      double precision, allocatable ::  tote3tIO(:,:,:) 
      double precision, allocatable ::  tottmaIO(:,:,:) 


      double precision, allocatable ::  trb(:,:,:,:)
      double precision, allocatable ::  buf(:,:,:)
      double precision, allocatable ::  buf2(:,:)
      INTEGER, allocatable, dimension(:) :: highfreq_table,highfreq_table_dia, highfreq_table_dia2d

!!----------------------------------------------------------------------
!!
!! COMMON /cot3ad/ non-centered advection scheme (smolarkiewicz)
!! -------------------------------------------------------------
!!      rsc         : tuning coefficient for anti-diffusion (NAMELIST)
!!      rtrn        : value for truncation (NAMELIST)

      double precision rsc,rtrn


!!----------------------------------------------------------------------
!!
!! COMMON /cit3ad/ non-centered advection scheme (smolarkiewicz)
!! -------------------------------------------------------------
!!      lhdf        : logical if true CALL trchdf (NAMELIST) 
!!      ncor        : number of corrective phases (NAMELIST)
!!      ndttrc      : frequency of step on passive tracers (NAMELIST)

      INTEGER ncor
      double precision ndttrc
      LOGICAL lhdf, ladv, lzdf, lsnu, lsbc
      LOGICAL hdf_initialized, adv_initialized 


!      isopycnal sheme for passive tracers
!! -----------------------------------------------------------
!!      ahtrb0    : background diffusivity coefficient (m2/s)
!!                  for passive tracer
!!      trcrat    : ratio between passive and active tracer coeff
!!                  for diffusion
!!      ahtrc0    : horizontal eddy diffusivity for passive tracers (m2/s)
!!    aeivtr0   : eddy induced velocity coefficient (m2/s)

      double precision ahtrb0,trcrat,ahtrc0,aeivtr0

      INTEGER nwritetrc ! time step frequency for concentration outputs (NAMELIST)


#    if defined key_trc_dmp 
      INTEGER(4), allocatable ::  idxt(:,:,:),idxt2glo(:,:,:,:)
#    endif

      LOGICAL IS_FREE_SURFACE
      LOGICAL lbfm      ! activates bfm model
      LOGICAL latmosph  ! activates atmospheric deposition


!!    Photoperiod formulation
      LOGICAL photop       ! Photoperiod formulation if false daylength is 24 h
      LOGICAL atlantic_bfm ! atlantic buffer biology activation

#    if defined key_trc_bfm
      double precision vsed                        ! sedimentation speed (NAMELIST)
      double precision vsedO5c                     ! sedimentation speed of calcite(NAMELIST)
      double precision bottom_flux                 ! (NAMELIST)

!!     optical parameters
      double precision, allocatable :: xpar(:,:,:) !par (photosynthetic available radiation)

#     endif
!!----------------------------------------------------------------------

      INTEGER ncpu

      integer :: nsubset
      integer,allocatable,dimension(:) :: subset_idx


      CONTAINS

! *******************************************************************
#ifdef key_mpp
       subroutine myalloc_sendrecv()
      INTEGER  :: err
      REAL(8)  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif

       allocate(EASTpoints_send( 2,EAST_count_send )) ; EASTpoints_send  = huge(EASTpoints_send( 1,1))
       allocate(WESTpoints_send( 2,WEST_count_send )) ; WESTpoints_send  = huge(WESTpoints_send( 1,1))
       allocate(NORTHpoints_send(2,NORTH_count_send)) ; NORTHpoints_send = huge(NORTHpoints_send(1,1))
       allocate(SOUTHpoints_send(2,SOUTH_count_send)) ; SOUTHpoints_send = huge(SOUTHpoints_send(1,1))

       allocate(EASTpoints_recv( 2,EAST_count_recv )) ; EASTpoints_recv  = huge(EASTpoints_recv( 1,1))
       allocate(WESTpoints_recv( 2,WEST_count_recv )) ; WESTpoints_recv  = huge(WESTpoints_recv( 1,1))
       allocate(NORTHpoints_recv(2,NORTH_count_recv)) ; NORTHpoints_recv = huge(NORTHpoints_recv(1,1))
       allocate(SOUTHpoints_recv(2,SOUTH_count_recv)) ; SOUTHpoints_recv = huge(SOUTHpoints_recv(1,1))
     
       allocate(te_send(EAST_count_send )) ; te_send = huge(te_send(1))
       allocate(tw_send(WEST_count_send )) ; tw_send = huge(tw_send(1))
       allocate(tn_send(NORTH_count_send)) ; tn_send = huge(tn_send(1))
       allocate(ts_send(SOUTH_count_send)) ; ts_send = huge(ts_send(1))

       allocate(te_recv(EAST_count_recv )) ; te_recv = huge(te_recv(1))
       allocate(tw_recv(WEST_count_recv )) ; tw_recv = huge(tw_recv(1))
       allocate(tn_recv(NORTH_count_recv)) ; tn_recv = huge(tn_recv(1))
       allocate(ts_recv(SOUTH_count_recv)) ; ts_recv = huge(ts_recv(1))

       !$acc enter data create( EASTpoints_send(1:2, 1:EAST_count_send), WESTpoints_send(1:2, 1:WEST_count_send),  NORTHpoints_send(1:2, 1:NORTH_count_send), SOUTHpoints_send(1:2, 1:SOUTH_count_send) )
       !$acc enter data create( EASTpoints_recv(1:2, 1:EAST_count_recv), WESTpoints_recv(1:2, 1:WEST_count_recv),  NORTHpoints_recv(1:2, 1:NORTH_count_recv), SOUTHpoints_recv(1:2, 1:SOUTH_count_recv) )

       !$acc enter data create( te_send(1:EAST_count_send ), tw_send(1:WEST_count_send ), tn_send(1:NORTH_count_send), ts_send(1:SOUTH_count_send) )
       !$acc enter data create( te_recv(1:EAST_count_recv ), tw_recv(1:WEST_count_recv ), tn_recv(1:NORTH_count_recv), ts_recv(1:SOUTH_count_recv) )

       
       
#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

       end subroutine myalloc_sendrecv
#endif
!*******************************************************************

subroutine myalloc_BFM()
      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err) 
#endif

       allocate(BFMpoints(3,NBFMPOINTS)) 
       BFMpoints = huge(BFMpoints(1,1))

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif

end subroutine myalloc_BFM
!*******************************************************************

subroutine alloc_tot()

      INTEGER  :: err
      double precision  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err) 
#endif
     

!      ncpu = 1
!      nimpp=1
!      njmpp=1

      allocate(mindi(jpi))    
      mindi  = huge(mindi(1))
      allocate(mindj(jpj))    
      mindj  = huge(mindj(1))   
      allocate(nimppt(jpnij)) 
      nimppt = huge(nimppt(1))
      allocate(njmppt(jpnij)) 
      njmppt = huge(njmppt(1))
      allocate(nlcit(jpnij))  
      nlcit  = huge(nlcit(1)) 
      allocate(nlcjt(jpnij))  
      nlcjt  = huge(nlcjt(1))
      allocate(nldit(jpnij))  
      nldit  = huge(nldit(1))
      allocate(nldjt(jpnij))  
      nldjt  = huge(nldjt(1))
      allocate(nleit(jpnij))  
      nleit  = huge(nleit(1))
      allocate(nlejt(jpnij))  
      nlejt  = huge(nlejt(1))


      
      allocate(glamt(jpj,jpi))
      glamt    = huge(glamt(1,1))
      allocate(glamu(jpj,jpi))          
      glamu    = huge(glamu(1,1))
      allocate(glamv(jpj,jpi))          
      glamv    = huge(glamv(1,1))
      allocate(glamf(jpj,jpi))          
      glamf    = huge(glamf(1,1))
      allocate(gphit(jpj,jpi))          
      gphit    = huge(gphit(1,1))
      allocate(gphiu(jpj,jpi))          
      gphiu    = huge(gphiu(1,1))
      allocate(gphiv(jpj,jpi))          
      gphiv    = huge(gphiv(1,1))
      allocate(gphif(jpj,jpi))          
      gphif    = huge(gphif(1,1))
      allocate(e1t(jpj,jpi))            
      e1t      = huge(e1t(1,1))
      allocate(e1u(jpj,jpi))            
      e1u      = huge(e1u(1,1))
      allocate(e1v(jpj,jpi))            
      e1v      = huge(e1v(1,1))
      allocate(e1f(jpj,jpi))            
      e1f      = huge(e1f(1,1))
      allocate(e2t(jpj,jpi))            
      e2t      = huge(e2t(1,1))
      allocate(e2u(jpj,jpi))            
      e2u      = huge(e2u(1,1))
      allocate(e2v(jpj,jpi))            
      e2v      = huge(e2v(1,1))
      allocate(e2f(jpj,jpi))            
      e2f      = huge(e2f(1,1))
      allocate(ff (jpj,jpi))            
      ff       = huge(ff(1,1))

#ifdef gdept1d
      allocate(gdept(jpk)) 
        gdept = huge(gdept(1))
#else
     allocate(gdept(jpk,jpj,jpi))
       gdept = huge(gdept(1,1,1))
#endif

      allocate(gdepw(jpk)) 
        gdepw = huge(gdepw(1))
      allocate(e3t_0(jpk,jpj,jpi))
        e3t_0 = huge(e3t_0(1,1,1))
      allocate(e3u_0(jpk,jpj,jpi))
        e3u_0 = huge(e3u_0(1,1,1))
      allocate(e3v_0(jpk,jpj,jpi))
        e3v_0 = huge(e3v_0(1,1,1))
      allocate(e3w_0(jpk,jpj,jpi))
        e3w_0 = huge(e3w_0(1,1,1))


      allocate(e3t(jpk,jpj,jpi)) 
        e3t = huge(e3t(1,1,1))
      allocate(e3t_back(jpk,jpj,jpi)) 
        e3t_back = huge(e3t_back(1,1,1))
      allocate(e3u(jpk,jpj,jpi)) 
        e3u = huge(e3u(1,1,1))
      allocate(e3v(jpk,jpj,jpi)) 
        e3v = huge(e3v(1,1,1))
      allocate(e3w(jpk,jpj,jpi)) 
        e3w = huge(e3w(1,1,1))

      allocate(mbathy(jpj,jpi)) 
       mbathy = huge(mbathy(1,1))

      allocate(tmask(jpk,jpj,jpi)) 
      tmask = huge(tmask(1,1,1))
      allocate(h_column(jpj,jpi))
      h_column = huge(h_column(1,1))
      
      allocate(spongeT(jpj,jpi))
      spongeT = huge(spongeT(1,1))
      allocate(spongeVel(jpk,jpj,jpi))
      spongeVel = huge(spongeVel(1,1,1))

      allocate(umask(jpk,jpj,jpi)) 
      umask = huge(umask(1,1,1))
      allocate(vmask(jpk,jpj,jpi)) 
      vmask = huge(vmask(1,1,1))

      allocate(bfmmask(jpk, jpj, jpi))
      bfmmask = huge(bfmmask(1, 1, 1))

       allocate(un(jpk,jpj,jpi))    
      un     = huge(un(1,1,1))
       allocate(vn(jpk,jpj,jpi))    
      vn     = huge(vn(1,1,1))
       allocate(wn(jpk,jpj,jpi))    
      wn     = huge(wn(1,1,1))
       allocate(tn(jpk,jpj,jpi))    
      tn     = huge(tn(1,1,1))
       allocate(sn(jpk,jpj,jpi))    
      sn     = huge(sn(1,1,1))
       allocate(rdn(jpk,jpj,jpi))   
      rdn    = huge(rdn(1,1,1))
       allocate(rhopn(jpk,jpj,jpi)) 
      rhopn  = huge(rhopn(1,1,1))
       allocate(rho(jpk,jpj,jpi))   
      rho    = huge(rho(1,1,1))
      
      allocate(ahtu(jpk)) 
      ahtu = huge(ahtu(1))
      allocate(ahtv(jpk)) 
      ahtv = huge(ahtv(1))
      allocate(ahtw(jpk)) 
      ahtw = huge(ahtw(1))
      allocate(ahtt(jpk)) 
      ahtt = huge(ahtt(1))

       allocate(avt (jpk,jpj,jpi)) 
       avt  = huge(avt(1,1,1))
       allocate(avtb(jpk))         
       avtb = huge(avtb(1))


       allocate(taux  (jpj,jpi)) 
       taux   = huge(taux(1,1))
       allocate(tauy  (jpj,jpi)) 
       tauy   = huge(tauy(1,1))
       allocate(vatm  (jpj,jpi)) 
       vatm   = huge(vatm(1,1))
       allocate(freeze(jpj,jpi)) 
       freeze = huge(freeze(1,1))



       allocate(qt    (jpj,jpi)) 
       qt     = huge(qt(1,1))
       allocate(q     (jpj,jpi)) 
       q      = huge(q(1,1))
       allocate(emp   (jpj,jpi)) 
       emp    = huge(emp(1,1))
       allocate(runoff(jpj,jpi)) 
       runoff = huge(runoff(1,1))

       allocate(qsr(jpj,jpi)) 
       qsr = huge(qsr(1,1))


       allocate(udta   (jpk,jpj,jpi,2))   
       udta    = huge(udta(1,1,1,1))
       allocate(vdta   (jpk,jpj,jpi,2))   
       vdta    = huge(vdta(1,1,1,1))
       allocate(wdta   (jpk,jpj,jpi,2))   
       wdta    = huge(wdta(1,1,1,1))
       allocate(avtdta (jpk,jpj,jpi,2))   
       avtdta  = huge(avtdta(1,1,1,1))
       allocate(flxdta (jpj,jpi,jpflx,2)) 
       flxdta  = huge(flxdta(1,1,1,1))
       allocate(flx    (jpj,jpi,jpflx))   
       flx     = huge(flx(1,1,1))
       allocate(tdta(jpk,jpj,jpi,2))      
       tdta    = huge(tdta(1,1,1,1))
       allocate(sdta(jpk,jpj,jpi,2))      
       sdta    = huge(sdta(1,1,1,1))
       allocate(e3tdta(jpk,jpj,jpi,2))    
       e3tdta  = huge(e3tdta(1,1,1,1))
       allocate(e3udta(jpk,jpj,jpi,2))    
       e3udta  = huge(e3udta(1,1,1,1))
       allocate(e3vdta(jpk,jpj,jpi,2))    
       e3vdta  = huge(e3vdta(1,1,1,1))
       allocate(e3wdta(jpk,jpj,jpi,2))    
       e3wdta  = huge(e3wdta(1,1,1,1))

!!----------------------------------------------------------------------


       allocate(trn(jpk,jpj,jpi,jptra))                    
       trn    = huge(trn(1,1,1,1))
       allocate(tra(jpk,jpj,jpi,jptra))                    
       tra    = huge(trn(1,1,1,1))
       allocate(tra_DIA(jptra_dia,jpk,jpj,jpi))            
       tra_DIA= huge(tra_DIA(1,1,1,1))
       allocate(tra_DIA_2d(jptra_dia_2d,jpj,jpi))
       tra_DIA_2d= huge(tra_DIA_2d(1,1,1))
       allocate(traIO(jpk,jpj,jpi,jptra))                  
       traIO  = huge(traIO(1,1,1,1)) 
       allocate(snIO(jpk,jpj,jpi))                         
       snIO   = huge(snIO(1,1,1))
       allocate(tnIO(jpk,jpj,jpi))                         
       tnIO   = huge(tnIO(1,1,1))
       allocate(vatmIO(jpj,jpi))                           
       vatmIO = huge(vatmIO(1,1))
       allocate(empIO(jpj,jpi))                            
       empIO  = huge(empIO(1,1))
       allocate(qsrIO(jpj,jpi))                            
       qsrIO  = huge(qsrIO(1,1))
       allocate(unIO  (jpk,jpj,jpi))                       
       unIO   = huge(unIO(1,1,1))
       allocate(vnIO  (jpk,jpj,jpi))                       
       vnIO   = huge(vnIO(1,1,1))
       allocate(wnIO  (jpk,jpj,jpi))                       
       wnIO          = huge(wnIO(1,1,1))
       allocate(avtIO (jpk,jpj,jpi))                       
       avtIO         = huge(avtIO(1,1,1))
       allocate(e3tIO (jpk,jpj,jpi))                       
       e3tIO         = huge(e3tIO(1,1,1))
       allocate(buf   (jpk,jpj,jpi))                       
       buf           = huge(buf(1,1,1))
       allocate(buf2   (jpj,jpi))                          
       buf2          = huge(buf2(1,1))
       allocate(tra_DIA_IO(jptra_dia,jpk,jpj,jpi))         
       tra_DIA_IO    = huge(tra_DIA_IO(1,1,1,1))
       allocate(traIO_HIGH(   jpk,jpj,jpi,jptra_HIGH))     
       traIO_HIGH    = huge(traIO_HIGH(1,1,1,1))
       allocate(tra_DIA_IO_HIGH(jptra_dia_HIGH,jpk,jpj,jpi))
       tra_DIA_IO_HIGH = huge(tra_DIA_IO_HIGH(1,1,1,1))

       allocate(tra_DIA_2d_IO(jptra_dia_2d, jpj,jpi))
       tra_DIA_2d_IO    = huge(tra_DIA_2d_IO(1,1,1))
       allocate(tra_DIA_2d_IO_HIGH(jptra_dia2d_HIGH,jpj,jpi))
       tra_DIA_2d_IO_HIGH = huge(tra_DIA_2d_IO_HIGH(1,1,1))
       
       allocate(tra_PHYS_IO(jptra_phys,jpk,jpj,jpi))
       tra_PHYS_IO    = huge(tra_PHYS_IO(1,1,1,1))
       allocate(tra_PHYS_IO_HIGH(jptra_phys,jpk,jpj,jpi))
       tra_PHYS_IO_HIGH = huge(tra_PHYS_IO_HIGH(1,1,1,1))
       allocate(tra_PHYS_2d_IO(jptra_phys_2d, jpj,jpi))
       tra_PHYS_2d_IO    = huge(tra_PHYS_2d_IO(1,1,1))
       allocate(tra_PHYS_2d_IO_HIGH(jptra_phys_2d,jpj,jpi))
       tra_PHYS_2d_IO_HIGH = huge(tra_PHYS_2d_IO_HIGH(1,1,1))


      if(lwp) then
       allocate(tottrn(jpk, jpjglo, jpiglo))      
       tottrn = huge(tottrn(1,1,1)) 
       allocate(tottrb(jpk, jpjglo, jpiglo))      
       tottrb = huge(tottrb(1,1,1))
!       allocate(tottrnIO(jpk,jpjglo,jpiglo))
!       tottrnIO  = huge(tottrnIO(1,1,1))
       allocate(tottrbIO(jpk,jpjglo,jpiglo)) 
       tottrbIO  = huge(tottrbIO(1,1,1)) 
       allocate(totsnIO (jpk,jpjglo,jpiglo)) 
       totsnIO   = huge(totsnIO(1,1,1))  
       allocate(tottnIO (jpk,jpjglo,jpiglo)) 
       tottnIO   = huge(tottnIO(1,1,1))  
       allocate(totvatmIO(jpjglo,jpiglo))    
       totvatmIO = huge(totvatmIO(1,1))  
       allocate(totempIO(jpjglo,jpiglo))     
       totempIO  = huge(totempIO(1,1))   
       allocate(totqsrIO(jpjglo,jpiglo))     
       totqsrIO  = huge(totqsrIO(1,1))   
       allocate(totunIO(jpk,jpjglo,jpiglo))  
       totunIO   = huge(totunIO(1,1,1))  
       allocate(totvnIO(jpk,jpjglo,jpiglo))  
       totvnIO   = huge(totvnIO(1,1,1))  
       allocate(totwnIO(jpk,jpjglo,jpiglo))  
       totwnIO   = huge(totwnIO(1,1,1))  
       allocate(totavtIO(jpk,jpjglo,jpiglo)) 
       totavtIO  = huge(totavtIO(1,1,1)) 
       allocate(tote3tIO(jpk,jpjglo,jpiglo)) 
       tote3tIO  = huge(tote3tIO(1,1,1)) 
       allocate(tottmaIO(jpk,jpjglo,jpiglo)) 
       tottmaIO  = huge(tottmaIO(1,1,1)) 
       allocate(tottrnIO2d(jpjglo,jpiglo))   
       tottrnIO2d= huge(tottrnIO2d(1,1))
!       allocate(totglamt(jpjglo,jpiglo))
!       totglamt = huge(totglamt(1,1))
       allocate(totgphit(jpjglo,jpiglo)) 
       totgphit = huge(totgphit(1,1))
      endif


       allocate(trb(jpk,jpj,jpi,jptra))              
       trb        = huge(trb(1,1,1,1))


#    if defined key_trc_dmp 

      allocate(idxt(jpk,jpj,jpi))           
       idxt     = huge(idxt(1,1,1))
      allocate(idxt2glo(jpk,jpj,jpi,4))     
       idxt2glo = huge(idxt2glo(1,1,1,1))

#    endif

#ifdef key_trc_bfm
      allocate(xpar(jpk,jpj,jpi))   
       xpar = huge(xpar(1,1,1))

#endif

!!    photoperiod
        allocate(DAY_LENGTH(jpj,jpi))   
       DAY_LENGTH = huge(DAY_LENGTH(1,1))
       forcing_phys_initialized = .false.
#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif
  
        END subroutine alloc_tot



        subroutine clean_memory()

            ! myalloc (memory.f90)

#ifdef key_mpp

            !$acc exit data delete(te_send, tw_send, tn_send, ts_send) finalize
            !$acc exit data delete(te_recv, tw_recv, tn_recv, ts_recv) finalize
          
            !$acc exit data delete(EASTpoints_send, WESTpoints_send, NORTHpoints_send, SOUTHpoints_send) finalize
            !$acc exit data delete(EASTpoints_recv, WESTpoints_recv, NORTHpoints_recv, SOUTHpoints_recv) finalize
          
            deallocate(EASTpoints_send)
            deallocate(WESTpoints_send)
            deallocate(NORTHpoints_send)
            deallocate(SOUTHpoints_send)
            deallocate(EASTpoints_recv)
            deallocate(WESTpoints_recv)
            deallocate(NORTHpoints_recv)
            deallocate(SOUTHpoints_recv)
            deallocate(te_send)
            deallocate(tw_send)
            deallocate(tn_send)
            deallocate(ts_send)
            deallocate(te_recv)
            deallocate(tw_recv)
            deallocate(tn_recv)
            deallocate(ts_recv)
#endif

            deallocate(BFMpoints)
            
            deallocate(mindi)
            deallocate(mindj)
            deallocate(nimppt)
            deallocate(njmppt)
            deallocate(nlcit)
            deallocate(nlcjt)
            deallocate(nldit)
            deallocate(nldjt)
            deallocate(nleit)
            deallocate(nlejt)
            
            deallocate(glamt)
            deallocate(glamu)
            deallocate(glamv)
            deallocate(glamf)
            deallocate(gphit)
            deallocate(gphiu)
            deallocate(gphiv)
            deallocate(gphif)
            deallocate(e1t)
            deallocate(e1u)
            deallocate(e1v)
            deallocate(e1f)
            deallocate(e2t)
            deallocate(e2u)
            deallocate(e2v)
            deallocate(e2f)
            deallocate(ff)
            
            deallocate(gdept)
            deallocate(gdepw)
            deallocate(e3t_0)
            deallocate(e3u_0)
            deallocate(e3v_0)
            deallocate(e3w_0)
            
            deallocate(e3t)
            deallocate(e3t_back)
            deallocate(e3u)
            deallocate(e3v)
            deallocate(e3w)
            
            deallocate(mbathy)
            
            deallocate(tmask)
            deallocate(h_column)
            
            deallocate(spongeT)
            deallocate(spongeVel)
            
            deallocate(umask)
            deallocate(vmask)
            
            deallocate(bfmmask)
            
            deallocate(un)
            deallocate(vn)
            deallocate(wn)
            deallocate(tn)
            deallocate(sn)
            deallocate(rdn)
            deallocate(rhopn)
            deallocate(rho)
            
            deallocate(ahtu)
            deallocate(ahtv)
            deallocate(ahtw)
            deallocate(ahtt)
            
            deallocate(avt)
            deallocate(avtb)
            
            deallocate(taux)
            deallocate(tauy)
            deallocate(vatm)
            deallocate(freeze)
            
            deallocate(qt)
            deallocate(q)
            deallocate(emp)
            deallocate(runoff)
            
            deallocate(qsr)
            
            deallocate(udta)
            deallocate(vdta)
            deallocate(wdta)
            deallocate(avtdta)
            deallocate(flxdta)
            deallocate(flx)
            deallocate(tdta)
            deallocate(sdta)
            deallocate(e3tdta)
            deallocate(e3udta)
            deallocate(e3vdta)
            deallocate(e3wdta)
            
            deallocate(trn)
            deallocate(tra)
            deallocate(tra_DIA)
            deallocate(tra_DIA_2d)
            deallocate(traIO)
            deallocate(snIO)
            deallocate(tnIO)
            deallocate(vatmIO)
            deallocate(empIO)
            deallocate(qsrIO)
            deallocate(unIO)
            deallocate(vnIO)
            deallocate(wnIO)
            deallocate(avtIO)
            deallocate(e3tIO)
            deallocate(buf)
            deallocate(buf2)
            deallocate(tra_DIA_IO)
            deallocate(tra_PHYS_IO)
            deallocate(traIO_HIGH)
            deallocate(tra_DIA_IO_HIGH)
            deallocate(tra_PHYS_IO_HIGH)
            deallocate(tra_DIA_2d_IO)
            deallocate(tra_DIA_2d_IO_HIGH)
            deallocate(tra_PHYS_2d_IO)
            deallocate(tra_PHYS_2d_IO_HIGH)
            

            if(lwp) then
                deallocate(tottrn)
                deallocate(tottrb)
!                deallocate(tottrnIO)
                deallocate(tottrbIO)
                deallocate(totsnIO)
                deallocate(tottnIO)
                deallocate(totvatmIO)
                deallocate(totempIO)
                deallocate(totqsrIO)
                deallocate(totunIO)
                deallocate(totvnIO)
                deallocate(totwnIO)
                deallocate(totavtIO)
                deallocate(tote3tIO)
                deallocate(tottmaIO)
                deallocate(tottrnIO2d)
!                deallocate(totglamt)
                deallocate(totgphit)
            endif
            
            deallocate(trb)

#ifdef key_trc_dmp
            deallocate(idxt)
            deallocate(idxt2glo)
#endif

#ifdef key_trc_bfm
            deallocate(xpar)
#endif

            deallocate(DAY_LENGTH)

            ! trclec

            deallocate(highfreq_table)
            deallocate(highfreq_table_dia)
            deallocate(highfreq_table_dia2d)

        end subroutine clean_memory



        END MODULE myalloc
