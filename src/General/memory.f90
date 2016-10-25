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
      INTEGER  ::threads_pack_size=18
      INTEGER nimpp,njmpp
      INTEGER nperio, narea, nlci, nlcj
      INTEGER nbondi, nbondj, nproc, noea, nowe, noso, nono
      INTEGER nreci, nrecj, nldi, nlei, nldj, nlej
      INTEGER, allocatable :: ilcit(:,:), ilcjt(:,:)
      INTEGER, allocatable :: mindi(:), mindj(:)
      INTEGER, allocatable :: nimppt(:), njmppt(:), nlcit(:), nlcjt(:)
      INTEGER, allocatable ::  nldit(:),  nldjt(:), nleit(:), nlejt(:)
      REAL(8)  mem_all

      INTEGER npolj

      REAL(8), PARAMETER ::  g =9.80665  ! gravity





!!----------------------------------------------------------------------
!!       ocean physical parameters (equation of state, ...)
!! ------------------------------------------
!!        neos         : flag of the type of equation of state used
!!      rau0             : reference volumic mass of the ocean (kg/m3)
!!      ralpha, rbeta    : thermique and haline expension coef. used
!!               for linear equation of state (neos=1 or 2)

      INTEGER neos
      REAL(8) rau0, ralpha, rbeta
      REAL(8)  rdt     ! dynamics time step

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


      REAL(8), allocatable, dimension(:,:) :: totglamt, glamu, glamv,glamf  !, glamt,
      REAL(8), allocatable, dimension(:,:) :: totgphit, gphiu, gphiv,gphif , gphit
      REAL(8), allocatable, dimension(:,:) :: e1t, e1u, e1v, e1f
      REAL(8), allocatable, dimension(:,:) :: e2t, e2u, e2v, e2f, ff

!!----------------------------------------------------------------------
!!       vertical coordinate and scale factors
!! -------------------------------------------------------

!!                  z-coordinate (default option)
!!                  ------------------------------
!!      gdept, gdepw() : depth of t- and w-points (m)
!!      e3t_0, e3w_0()     : vertical scale factors at t- and w-points (m)
!!
      REAL(8), allocatable :: gdept(:), gdepw(:), e3t_0(:), e3w_0(:)
      REAL(8), allocatable :: e3t(:,:,:), e3t_back(:,:,:), e3u(:,:,:), e3v(:,:,:), e3w(:,:,:)

!!----------------------------------------------------------------------
!!        masks, bathymetry
!! -----------------------------------
!!      mbathy         : number of ocean level (=0, 1, ... , jpk-1)
!!      tmask, umask() : land/ocean mask at t-, u-, v- and f-points
!!      vmask, fmask()

      INTEGER, allocatable :: mbathy(:,:)


      REAL(8), allocatable, dimension(:,:,:) :: tmask, fmask,umask, vmask
      INTEGER NBFMPOINTS, NBFMPOINTS_SUP, NWATERPOINTS
      INTEGER, allocatable, dimension(:,:) :: BFMpoints


!! II. DYNAMICS AND TRACERS
!! ========================
!!----------------------------------------------------------------------
!!       previous fields (before)
!! -----------------------------------------

      REAL(8), allocatable, dimension(:,:,:) :: ub, vb ! horizontal velocity (m s-1)


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
      REAL(8), allocatable, dimension(:,:,:) :: un, vn, wn
      REAL(8), allocatable, dimension(:,:,:) :: tn, sn,rdn,rhopn,rho,bn2n
      REAL(8), allocatable, dimension(:,:,:) :: hdivn


!!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!! III. OCEAN PHYSICS
!! ==================

!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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

      REAL(8) aht0
      REAL(8), allocatable :: ahtu(:), ahtv(:), ahtw(:), ahtt(:)


!!----------------------------------------------------------------------
!!         vertical diffusion
!! ----------------------------------------------------------------------
!!     avt             : vertical diffusivity coeff. at w-point
!!     avtb            : background profile of avm and avt

      REAL(8), allocatable :: avt(:,:,:)
      REAL(8), allocatable ::  avtb(:)


!! IV. SURFACE FORCING AND DATA
!! ============================
!!  surface wind stress at givem time_step
!!    taux, tauy()   : wind stress components in (i,j) referential
      REAL(8), allocatable, dimension(:,:) :: taux, tauy, vatm, freeze

!!----------------------------------------------------------------------
!!     surface fluxes
!! -------------------------------
!!      qt             : total surface heat flux (w m-2)
!!      q              : surface heat flux (w m-2)
!!      emp            : evaporation minus precipitation (mm day-1)
!!      runoff         : annual run off (mm/day)

      REAL(8), allocatable, dimension(:,:) :: qt, q, emp,runoff
      REAL(8), allocatable :: qsr(:,:) ! penetrative solar radiation (w m-2)

      INTEGER nsptint ! YPE of spatial interpolation (NAMELIST)


!!      udta,vdta()  : horizontal velocity data array
!!      wdta         : vertical velocity data array
!!      avtdta       : avt data array
!!      flxdta         : additional fluxes

      REAL(8), allocatable, dimension(:,:,:,:) :: udta,vdta,wdta,avtdta,flxdta
      REAL(8), allocatable, dimension(:,:,:)   :: flx
      REAL(8), allocatable, dimension(:,:,:,:) :: tdta,sdta ! : temperature and salinity data array
      REAL(8), allocatable, dimension(:,:,:,:) :: e3tdta,e3udta,e3vdta,e3wdta ! : temperature and salinity data array

!! V. DIAGNOSTICS
!! ==============


      INTEGER calendarType ! leap years calendar (0/1)

      REAL(8), ALLOCATABLE, DIMENSION(:,:) :: DAY_LENGTH

      INTEGER            :: numnam = 208 ! unit for namelist
      INTEGER, PARAMETER :: numout = 2   ! unit for output print
      LOGICAL lwp                        ! boolean term for stdout
      INTEGER, PARAMETER :: numnat =80   ! the number of the passive tracer NAMELIST


#if defined key_mpp 
!!      t3ew           : 3d message passing arrays east-west    
!!      t3we           : 3d message passing arrays west-east
!!      t3ns           : 3d message passing arrays north-south
!!      t3sn           : 3d message passing arrays south-north 
!!      t2ew           : 2d message passing arrays east-west    
!!      t2we           : 2d message passing arrays west-east
!!      t2ns           : 2d message passing arrays north-south
!!      t2sn           : 2d message passing arrays south-north 

      REAL(8), allocatable :: t3ns (:,:,:,:), t3sn (:,:,:,:)
      REAL(8), allocatable :: t3ew (:,:,:,:), t3we (:,:,:,:)
      REAL(8), allocatable :: t3ew_my1 (:,:,:,:,:), t3we_my1 (:,:,:,:,:)
      REAL(8), allocatable :: t3sn_my1 (:,:,:,:,:), t3ns_my1 (:,:,:,:,:)
      REAL(8), allocatable :: t2ns (:,:,:)  , t2sn (:,:,:)
      REAL(8), allocatable :: t2ew (:,:,:)  , t2we (:,:,:)

#  else
!     no mpp
#endif



!! PASSIVE TRACER MODEL
          CHARACTER(LEN=20) :: ctrcnm(jptra)
          CHARACTER(LEN=12) :: ctrcun(jptra)
          CHARACTER(LEN=20) :: dianm(jptra_dia)
          CHARACTER(LEN=20) :: diaun(jptra_dia)
          INTEGER           :: diahf(jptra_dia)
          CHARACTER(LEN=20) :: dianm_2d(jptra_dia_2d)
          CHARACTER(LEN=20) :: diaun_2d(jptra_dia_2d)
          INTEGER           :: diahf_2d(jptra_dia_2d)

!!    parameters for the control of passive tracers



      REAL(8) ::  ctrmax(jptra)
      LOGICAL :: isCheckLOG
      LOGICAL :: save_bkp_group2 ! we can avoid to dump bkp of a lot of variables
      INTEGER :: jptra_high, jptra_dia_high, jptra_dia2d_high
      INTEGER :: ctr_hf(jptra)

      INTEGER freq_ave_phys



      INTEGER flagSMS_Dyn                    ! Flag time advance SMS or Dyn
      REAL(8), allocatable ::  trn(:,:,:,:)
      REAL(8), allocatable ::  tra(:,:,:,:)
      REAL(8), allocatable ::  tra_DIA(:,:,:,:)
      REAL(8), allocatable ::  tra_DIA_2d(:,:,:)
      REAL(8), allocatable ::  traIO(:,:,:,:)
      REAL(8), allocatable ::  traIO_HIGH(:,:,:,:)
      REAL(8), allocatable ::  snIO(:,:,:) 
      REAL(8), allocatable ::  tnIO(:,:,:) 
      REAL(8), allocatable ::  vatmIO(:,:) 
      REAL(8), allocatable ::  empIO(:,:) 
      REAL(8), allocatable ::  qsrIO(:,:) 
      REAL(8), allocatable ::  unIO(:,:,:) 
      REAL(8), allocatable ::  vnIO(:,:,:) 
      REAL(8), allocatable ::  wnIO(:,:,:) 
      REAL(8), allocatable ::  avtIO(:,:,:) 
      REAL(8), allocatable ::  e3tIO(:,:,:) 
      REAL(8), allocatable ::  tra_DIA_IO(:,:,:,:)
      REAL(8), allocatable ::  tra_DIA_IO_HIGH(:,:,:,:)
      REAL(8), allocatable ::  tra_DIA_2d_IO(:,:,:)
      REAL(8), allocatable ::  tra_DIA_2d_IO_HIGH(:,:,:)
      REAL(8), allocatable :: tottrn(:,:,:)
      REAL(8), allocatable :: tottrb(:,:,:)

      REAL(8), allocatable ::  tottrnIO(:,:,:) ! matrix for i/o writing(trcdit.F)
      REAL(8), allocatable ::  tottrnIO2d(:,:)
      REAL(8), allocatable ::  tottrbIO(:,:,:)
      REAL(8), allocatable ::  totsnIO(:,:,:) 
      REAL(8), allocatable ::  tottnIO(:,:,:) 
      REAL(8), allocatable ::  totvatmIO(:,:) 
      REAL(8), allocatable ::  totempIO(:,:) 
      REAL(8), allocatable ::  totqsrIO(:,:) 
      REAL(8), allocatable ::  totunIO(:,:,:) 
      REAL(8), allocatable ::  totvnIO(:,:,:) 
      REAL(8), allocatable ::  totwnIO(:,:,:) 
      REAL(8), allocatable ::  totavtIO(:,:,:) 
      REAL(8), allocatable ::  tote3tIO(:,:,:) 
      REAL(8), allocatable ::  tottmaIO(:,:,:) 


      REAL(8), allocatable ::  trb(:,:,:,:)
      REAL(8), allocatable ::  buf(:,:,:)
      REAL(8), allocatable ::  buf2(:,:)
      INTEGER, allocatable, dimension(:) :: highfreq_table,highfreq_table_dia, highfreq_table_dia2d

!!----------------------------------------------------------------------
!!
!! COMMON /cot3ad/ non-centered advection scheme (smolarkiewicz)
!! -------------------------------------------------------------
!!      rsc         : tuning coefficient for anti-diffusion (NAMELIST)
!!      rtrn        : value for truncation (NAMELIST)

      REAL(8) rsc,rtrn


!!----------------------------------------------------------------------
!!
!! COMMON /cit3ad/ non-centered advection scheme (smolarkiewicz)
!! -------------------------------------------------------------
!!      lhdf        : logical if true CALL trchdf (NAMELIST) 
!!      ncor        : number of corrective phases (NAMELIST)
!!      ndttrc      : frequency of step on passive tracers (NAMELIST)

      INTEGER ncor
      REAL(8) ndttrc
      LOGICAL lhdf


!      isopycnal sheme for passive tracers
!! -----------------------------------------------------------
!!      ahtrb0    : background diffusivity coefficient (m2/s)
!!                  for passive tracer
!!      trcrat    : ratio between passive and active tracer coeff
!!                  for diffusion
!!      ahtrc0    : horizontal eddy diffusivity for passive tracers (m2/s)
!!    aeivtr0   : eddy induced velocity coefficient (m2/s)

      REAL(8) ahtrb0,trcrat,ahtrc0,aeivtr0

      INTEGER nwritetrc ! time step frequency for concentration outputs (NAMELIST)


#    if defined key_trc_dmp 
      INTEGER(4), allocatable ::  idxt(:,:,:),idxt2glo(:,:,:,:)
#    endif


!!    Photoperiod formulation
      LOGICAL photop       ! Photoperiod formulation if false daylength is 24 h
      LOGICAL atlantic_bfm ! atlantic buffer biology activation

#    if defined key_trc_bfm
      REAL(8) vsed                        ! sedimentation speed (NAMELIST)
      REAL(8) bottom_flux                 ! (NAMELIST)

!!     optical parameters
      REAL(8), allocatable :: xpar(:,:,:) !par (photosynthetic available radiation)

#     endif
!!----------------------------------------------------------------------

      INTEGER ncpu



      CONTAINS


!*******************************************************************

subroutine myalloc_BFM()
      INTEGER  :: err
      REAL(8)  :: aux_mem

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
      REAL(8)  :: aux_mem

#ifdef Mem_Monitor
       aux_mem = get_mem(err) 
#endif
     

      ncpu = 1
      nimpp=1
      njmpp=1

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


      allocate(gdept(jpk)) 
        gdept = huge(gdept(1))
      allocate(gdepw(jpk)) 
        gdepw = huge(gdepw(1))
      allocate(e3t_0(jpk)) 
        e3t_0 = huge(e3t_0(1))
      allocate(e3w_0(jpk)) 
        e3w_0 = huge(e3w_0(1))

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
      allocate(fmask(jpk,jpj,jpi)) 
      fmask = huge(fmask(1,1,1))
      allocate(umask(jpk,jpj,jpi)) 
      umask = huge(umask(1,1,1))
      allocate(vmask(jpk,jpj,jpi)) 
      vmask = huge(vmask(1,1,1))

      allocate(ub(jpk,jpj,jpi)) 
      ub = huge(ub(1,1,1))
      allocate(vb(jpk,jpj,jpi)) 
      vb = huge(vb(1,1,1))

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
       allocate(bn2n(jpk,jpj,jpi))  
      bn2n   = huge(bn2n(1,1,1))
       allocate(hdivn(jpk,jpj,jpi)) 
      hdivn  = huge(hdivn(1,1,1))

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

#if defined key_mpp 
       allocate(t3ns(jpi,jprecj,jpk,2))            
       t3ns     = huge(t3ns(1,1,1,1))
       allocate(t3sn(jpi,jprecj,jpk,2))            
       t3sn     = huge(t3sn(1,1,1,1))
       allocate(t3ew(jpj,jpreci,jpk,2))            
       t3ew     = huge(t3ew(1,1,1,1))
       allocate(t3we(jpj,jpreci,jpk,2))            
       t3we     = huge(t3we(1,1,1,1))
       allocate(t3ew_my1(jpj,jpreci,jpk,jptra,2))  
       t3ew_my1 = huge(t3ew_my1(1,1,1,1,1))
       allocate(t3we_my1(jpj,jpreci,jpk,jptra,2))  
       t3we_my1 = huge(t3we_my1(1,1,1,1,1))
       allocate(t3sn_my1(jpi,jpreci,jpk,jptra,2))  
       t3sn_my1 = huge(t3sn_my1(1,1,1,1,1))
       allocate(t3ns_my1(jpi,jpreci,jpk,jptra,2))  
       t3ns_my1 = huge(t3ns_my1(1,1,1,1,1))
       allocate(t2ns(jpi,jprecj,2))                
       t2ns     = huge(t2ns(1,1,1))
       allocate(t2sn(jpi,jprecj,2))                
       t2sn     = huge(t2sn(1,1,1))
       allocate(t2ew(jpj,jpreci,2))                
       t2ew     = huge(t2ew(1,1,1))
       allocate(t2we(jpj,jpreci,2))                
       t2we     = huge(t2we(1,1,1))

#  else
!     no mpp
#endif

       allocate(trn(jpk,jpj,jpi,jptra))                    
       trn    = huge(trn(1,1,1,1))
       allocate(tra(jpk,jpj,jpi,jptra))                    
       tra    = huge(trn(1,1,1,1))
       allocate(tra_DIA(jpk,jpj,jpi,jptra_dia))            
       tra_DIA= huge(tra_DIA(1,1,1,1))
       allocate(tra_DIA_2d(jpj,jpi,jptra_dia_2d))          
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
       allocate(tra_DIA_IO(jpk,jpj,jpi,jptra_dia))         
       tra_DIA_IO    = huge(tra_DIA_IO(1,1,1,1))
       allocate(traIO_HIGH(   jpk,jpj,jpi,jptra_HIGH))     
       traIO_HIGH    = huge(traIO_HIGH(1,1,1,1))
       allocate(tra_DIA_IO_HIGH(jpk,jpj,jpi,jptra_dia_HIGH))
      tra_DIA_IO_HIGH = huge(tra_DIA_IO_HIGH(1,1,1,1))

       allocate(tra_DIA_2d_IO(jpj,jpi,jptra_dia_2d))       
       tra_DIA_2d_IO    = huge(tra_DIA_2d_IO(1,1,1))
       allocate(tra_DIA_2d_IO_HIGH(jpj,jpi,jptra_dia2d_HIGH))
      tra_DIA_2d_IO_HIGH = huge(tra_DIA_2d_IO_HIGH(1,1,1))



      if(lwp) then
       allocate(tottrn(jpiglo, jpjglo, jpk))      
       tottrn = huge(tottrn(1,1,1)) 
       allocate(tottrb(jpiglo, jpjglo, jpk))      
       tottrb = huge(tottrb(1,1,1))
       allocate(tottrnIO(jpk,jpjglo,jpiglo)) 
       tottrnIO  = huge(tottrnIO(1,1,1)) 
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
       allocate(totglamt(jpjglo,jpiglo)) 
       totglamt = huge(totglamt(1,1))
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

#ifdef Mem_Monitor
      mem_all=get_mem(err) - aux_mem
#endif
  
        END subroutine alloc_tot

        END MODULE 
