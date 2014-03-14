       MODULE myalloc

       USE modul_param
       USE timers

       IMPLICIT NONE

       public

CC----------------------------------------------------------------------
CC            domain parameters
CC ---------------------------------------------------------------------
CC      nizoom, njzoom   : (i,j) indexes for the domain left bottom zoom
CC      nperio           : lateral boundary condition type 
CC      nimpp,njmpp      : (i,j) indexes for mpp-subdomain left bottom
CC      nreci,nrecj      : overlap region in i and j
CC      nproc            : number for local processor
CC      narea            : number for local area
CC      nbondi, nbondj   : mark of i- and j-direction local boundaries
CC      nlci, nlcj       : i, j dimensions of the local subdomain
CC      nldi, nlei,      : first and last indoor i- and j-indexes
CC      nldj, nlej   
CC      noea, nowe,      : index of the local neighboring processors in
CC      noso, nono         east, west, south and north directions
CC      nimppt,njmppt(): i-, j-indexes for each processor
CC      nlcit, nlcjt() : dimensions of every subdomain
CC      nldit, nldjt() : first, last indoor index for each i-domain
CC      nleit, nlejt() : first, last indoor index for each j-domain
CC      mindi, mindj() : indexes array of the subdomain

      INTEGER nizoom, njzoom, nimpp, njmpp, nperio, narea, nlci, nlcj
      INTEGER nbondi, nbondj, nproc, noea, nowe, noso, nono
      INTEGER nreci, nrecj, nldi, nlei, nldj, nlej
      INTEGER, allocatable :: ilcit(:,:), ilcjt(:,:)
      INTEGER, allocatable :: mindi(:), mindj(:)
      INTEGER, allocatable :: nimppt(:), njmppt(:), nlcit(:), nlcjt(:)
      INTEGER, allocatable ::  nldit(:),  nldjt(:), nleit(:), nlejt(:)
      REAL(8) mem_all,rea_len,int_len

      INTEGER npolj


CC----------------------------------------------------------------------
CC        physical constants
CC -----------------------------------
CC      nbiss          : number of days per month (leap-year)
CC      nobis          : number of days per month (365 days a year)
CC      rpi              : pi
CC      rday             : day
CC      rsiyea           : sideral year
CC      rsiday           : sideral day  
CC      omega            : earth rotation parameter
CC      ra               : earth radius (meter)
CC      g                : gravity
CC      rtt              : fusion point
CC      rt0              : zero celsius
CC      rauw             : density of pure water kg/m3
CC      rad              : conversion coeff. from degre into radian
CC      rcs              : conversion coeff. from millimeter per day
CC                         into meter per second
CC      raajj            : number of days in one year
CC      raamo            : number of months in one year
CC      rjjhh            : number of hours in one day
CC      rhhmm            : number of minutes in one hour
CC      rmmss            : number of seconds in one minute
CC      raass            : number of seconds in one year
CC      rmoss            : number of seconds in one month
CC      rjjss            : number of seconds in one day

      INTEGER nbiss(12), nobis(12)


      REAL(8) rpi, rday, rsiyea, rsiday, omega, ra, g, rtt, rt0, rauw, rad,
     $    rcs, raajj, raamo, rjjhh, rhhmm, rmmss, raass, rmoss, rjjss 





CC----------------------------------------------------------------------
CC       ocean physical parameters (equation of state, ...)
CC ------------------------------------------
CC        neos         : flag of the type of equation of state used
CC      rau0             : reference volumic mass of the ocean (kg/m3)
CC      ralpha, rbeta    : thermique and haline expension coef. used
CC               for linear equation of state (neos=1 or 2)
CC      rcp              : ocean specific heat
CC      ro0cpr           : = 1. / ( rau0 * rcp )

      INTEGER neos
      REAL(8) rau0, ralpha, rbeta, rcp, ro0cpr


CC----------------------------------------------------------------------
CC       time step
CC ---------------------------------------------------------------------
CC      nacc             : flag to activate the convergence acceleration
CC      atfp, atfp1      : asselin time filter param. (atfp1= 1-2*atfp)
CC      rdt              : dynamics time step
CC      rdttra         : vertical profile of tracer time step
CC      rdtmin           : minimum time step on tracers
CC      rdtmax           : maximum time step on tracers
CC      rdth             : depth scale for vertically dependent rdttra
CC
      INTEGER nacc

      REAL(8) atfp, atfp1, rdt, rdtmin, rdtmax, rdth
      REAL(8), allocatable :: rdttra(:)

CC----------------------------------------------------------------------
CC        horizontal curvilinear coordinate and scale factors
CC ---------------------------------------------------------------------
CC      nummsh           : logical unit of coordinate/meshes output file
CC      ngrid            : flag to read or compute coordinates
CC      nmsh             : flag to write coordinates and meshes file
CC      njpi, njpj       : read dimensions
CC      rra              : read earth radius
CC      glamt          : longitude of t-point (degre)
CC      glamu          : longitude of u-point (degre)
CC      glamv          : longitude of v-point (degre)
CC      glamf          : longitude of f-point (degre)
CC      gphit          : latitude  of t-point (degre)
CC      gphiu          : latitude  of u-point (degre)
CC      gphiv          : latitude  of v-point (degre)
CC      gphif          : latitude  of f-point (degre)
CC      e1t,e2t        : horizontal scale factors at t-point (m)
CC      e1u,e2u        : horizontal scale factors at u-point (m)
CC      e1v,e2v        : horizontal scale factors at v-point (m)
CC      e1f,e2f        : horizontal scale factors at f-point (m)
CC        ff             : coriolis factor (2.*omega*sin(yphi) ) (s-1)

      INTEGER nummsh, ngrid, nmsh, njpi, njpj


      REAL(8) rra
      REAL(8), allocatable, dimension(:,:) :: totglamt, glamu, glamv,glamf  !, glamt,
      REAL(8), allocatable, dimension(:,:) :: totgphit, gphiu, gphiv,gphif , gphit
      REAL(8), allocatable, dimension(:,:) :: e1t, e1u, e1v, e1f
      REAL(8), allocatable, dimension(:,:) :: e2t, e2u, e2v, e2f, ff

#if defined key_off_degrad
      REAL(8), allocatable :: facvol(:,:,:)
#endif

CC----------------------------------------------------------------------
CC       vertical coordinate and scale factors
CC -------------------------------------------------------


      INTEGER numhba    !  logical unit for bathymetry file (meters)


#ifdef key_s_coord
CC                   s-coordinate ('key_s_coord')
CC                   ----------------------------
CC      hbatt, hbatu() : ocean depth at the vertical of  t-, u-, v-
CC      hbatv, hbatf()   and f-point (m)
CC      gsigt, gsigw() : model level depth coefficient at t-, w-levels
CC      gsi3w          : model level depth coefficient at w-level
CC                         defined as the sum of e3w scale factors
CC      esigt, esigw() : vertical scale factor coef. at t-, w-levels

      REAL(8), allocatable, dimension(:,:) :: hbatt, hbatu, hbatv,hbatf
      REAL(8), allocatable, dimension(:)   :: gsigt, gsigw, gsi3w, esigt, esigw

#  else
CC                  z-coordinate (default option)
CC                  ------------------------------
CC      gdept, gdepw() : depth of t- and w-points (m)
CC      e3t, e3w()     : vertical scale factors at t- and w-points (m)
CC
      REAL(8), allocatable :: gdept(:), gdepw(:), e3t(:), e3w(:)

#endif

CC----------------------------------------------------------------------
CC        masks, bathymetry
CC -----------------------------------
CC      numbat          : logical unit for bathymetry file
CC      mbathy         : number of ocean level (=0, 1, ... , jpk-1)
CC      tmask, umask() : land/ocean mask at t-, u-, v- and f-points
CC      vmask, fmask()

      INTEGER numbat
      INTEGER, allocatable :: mbathy(:,:)


      REAL(8), allocatable, dimension(:,:,:) :: tmask, fmask,umask, vmask
      INTEGER NBFMPOINTS, NBFMPOINTS_SUP
      INTEGER, allocatable, dimension(:,:) :: BFMpoints

CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC
CC II. DYNAMICS AND TRACERS
CC ========================
CC
CC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CC
CC----------------------------------------------------------------------
CC       previous fields (before)
CC -----------------------------------------

      REAL(8), allocatable, dimension(:,:,:) :: ub, vb ! horizontal velocity (m s-1)


CC----------------------------------------------------------------------
CC      present fields (now)
CC -------------------------------------
CC       un, vn(), wn() : horizontal and vertical velocity (m s-1)
CC        hdivn          : horizontal divergence
CC        tn,   sn()     : pot. temperature (celsius), salinity (psu)
CC      rdn            : in situ density anomalie rdn=(rho-rau0)/rau0
CC                         (no units)
CC        rhopn          : potential volumic mass (kg m-3)
CC      bn2n           : brunt-vaisala frequency (s-2)
CC
      REAL(8), allocatable, dimension(:,:,:) :: un, vn, wn, hdivn
      REAL(8), allocatable, dimension(:,:,:) :: tn, sn,rdn,rhopn,rho,bn2n


CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC III. OCEAN PHYSICS
CC ==================

CC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

CC----------------------------------------------------------------------
CC      lateral diffusivity (tracers)
CC ------------------------------------
CC      aht0             : lateral diffusivity coefficient    (namelist)
CC      ahtb0            : background diffusivity coefficient (m2/s)
CC      ahtu, ahtv()   : lateral diffusivity coef. at u-, v-, w- t-pts
CC      ahtw, ahtt()     (harmonic operator: no rotation, use of u-
CC                          and v-points; rotation, use of u-, v- w-pts)
CC                         (biharmonic operator: rotation or not, use of
CC                          t-point only)
CC                         (the arrays used are 3D, 2D, 1D or 0D depen-
CC                          ding on 'key_trahdfcoef.d' )
C
      REAL(8) aht0,ahtb0
CCC    Paolo 22/4/2004 i killed all options
#if defined key_trahdfbilap
C      REAL(8) ahtu, ahtv, ahtw
C      REAL(8) ahtt
# else
#  if defined key_trahdfcoef3d
C      REAL(8), allocatable :: ahtu(:,:,:), ahtv(:,:,:), ahtw(:,:,:)
#   elif key_trahdfcoef2d
C      REAL(8), allocatable :: ahtu(:,:), ahtv(:,:), ahtw(:,:)
#   elif key_trahdfcoef1d
C      REAL(8), allocatable :: ahtu(:), ahtv(:), ahtw(:)
#   else
#     if defined key_off_degrad
C      REAL(8), allocatable :: ahtu(:,:,:), ahtv(:,:,:), ahtw(:,:,:)
#     else
C      REAL(8) ahtu, ahtv, ahtw
#     endif
#  endif
C      REAL(8) ahtt
#endif

      REAL(8), allocatable :: ahtu(:), ahtv(:), ahtw(:), ahtt(:)

CC
CC----------------------------------------------------------------------
CC       direction of lateral diffusion (momentum & tracers)
CC -----------------------------------------------
CC      Slopes ('key_trahdfiso' or 'key_s_coord' defined)
CC      uslp, wslpi()  : i-direction slope at u-, w-points
CC      vslp, wslpj()  : j-direction slope at v-, w-points
CC          wk1, wk2()     : work array used in rotated biharmonic ope-
CC         wk3, wk4()       rator on momentum and/or tracers
C
#if defined key_trahdfiso || defined key_trahdfgeop || defined key_debug
      REAL(8), allocatable :: uslp (:,:,:), vslp (:,:,:)
      REAL(8), allocatable :: wslpi(:,:,:), wslpj(:,:,:)
#  if defined key_trahdfbilap
      REAL(8), allocatable :: wk1(:,:,:), wk2(:,:,:)
      REAL(8), allocatable :: wk3(:,:,:), wk4(:,:,:)
#   else
      REAL(8) wk1, wk2, wk3, wk4
#  endif

# else
C     diffusion along model level surfaces: no arrays
#endif

CC----------------------------------------------------------------------
CC        eddy induced velocity parameterization
CC -------------------------------------------------------
CC      aeiv0            : eddy induced velocity coefficient (m2/s)
CC      aeiu, aeiv()   : eddy induced velocity coef. at u-, v- and
CC      aeiw             w-points (m2/s). The arrays used are 3D
CC                         2D, 1D or 0D, depending on 'key_trahdfcoef.d'

      REAL(8) aeiv0
#if defined key_trahdfeiv
#  if defined key_trahdfcoef3d
      REAL(8), allocatable, dimension(:,:,:) :: aeiu, aeiv, aeiw
#   elif key_trahdfcoef2d
      REAL(8), allocatable, dimension(:,:)   :: aeiu, aeiv, aeiw
#   elif key_trahdfcoef1d
      REAL(8), allocatable, dimension(:)     :: aeiu, aeiv, aeiw
#   else
      REAL(8) aeiu, aeiv, aeiw
#  endif
# else
#     if defined key_off_degrad
      REAL(8), allocatable, dimension(:,:,:) :: aeiu, aeiv, aeiw
#     else
      REAL(8) aeiu, aeiv, aeiw
#     endif
#endif

CC-------------------------------------------------------------------
CC      eddy induced velocities (Gent et Mc William)
CC---------------------------------------------------------------
CC       ugm  : eddy induced velocity
CC       vgm
CC       wgm

#if defined key_trahdfeiv
      REAL(8), allocatable, dimension(:,:,:) :: ugm,vgm,wgm
#endif
CC
CC----------------------------------------------------------------------
CC         vertical diffusion
CC -----------------------------------
CC     avt0            : vertical viscosity and diff. coef. (namelist)
CC      ntrbbl         : bottom boundary layer (namelist)
CC      atrbbl         : bottom boundary layer diffusivity (namelist)
CC     avt             : vertical diffusivity coeff. at w-point
CC      avtb           : background profile of avm and avt
C
C
      INTEGER ntrbbl
      REAL(8) avt0,atrbbl
      REAL(8), allocatable :: avt(:,:,:)
      REAL(8), allocatable ::  avtb(:)
      REAL(8), allocatable ::  bblx(:,:),bbly(:,:)

CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC
CC IV. SURFACE FORCING AND DATA
CC ============================
CC
CC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CC
CC----------------------------------------------------------------------
CC  surface wind stress at givem time_step
CC -------------------------------------------------------
CC    taux, tauy()   : wind stress components in (i,j) referential
CC    tauxg, tauyg() : zonal and meridian wind stress component used
CC                       in output (geographical referential)

      REAL(8), allocatable, dimension(:,:) :: taux, tauy, vatm, freeze

CC----------------------------------------------------------------------
CC     surface fluxes
CC -------------------------------
CC      qt             : total surface heat flux (w m-2)
CC      q              : surface heat flux (w m-2)
CC      emp            : evaporation minus precipitation (mm day-1)
CC      runoff         : annual run off (mm/day)

      REAL(8), allocatable, dimension(:,:) :: qt, q, emp,runoff

CC
CC----------------------------------------------------------------------
CC
CC      penetrative solar radiation
CC --------------------------------------------
CC      qsr     : solar radiation (w m-2)

      REAL(8), allocatable :: qsr(:,:)

CC----------------------------------------------------------------------

CC     dynamic data selection
CC ---------------------------------------
CC      numdyn    : logical unit for dynamic DATA read
CC      numtke    : logical unit for tke DATA read
CC      nsptint   : TYPE of spatial interpolation (NAMELIST)

      INTEGER numdyn,numtke, nsptint,numslp,numvgm


CC      udta,vdta()  : horizontal velocity data array
CC      wdta         : vertical velocity data array
CC      avtdta       : avt data array
CC      flxdta         : additional fluxes

      REAL(8), allocatable, dimension(:,:,:,:) :: udta,vdta,wdta,avtdta,flxdta
      REAL(8), allocatable, dimension(:,:,:)   :: flx, bblxdta, bblydta, zmxldta
      REAL(8), allocatable, dimension(:,:)     :: hmld

#if defined key_off_tra
      REAL(8), allocatable, dimension(:,:,:,:) :: tdta,sdta ! : temperature and salinity data array
#endif

#if defined key_trahdfiso
CC----------------------------------------------------------------------
CC
CC      dynamic data selection
CC ---------------------------------------
CC      uslpdta      : zonal slope (u point) data array
CC      vslpdta      : meridian slope (v point) data array
CC      wslpidta     : zonal slope (t point) data array
CC      wslpjdta     : meridian slope (t point) data array
CC
      REAL(8), allocatable, dimension(:,:,:,:) ::  uslpdta,vslpdta,wslpidta,wslpjdta

#endif

#if defined key_trahdfeiv

CC       dynamic data selection
CC ---------------------------------------

      REAL(8), allocatable :: ugmdta(:,:,:,:) ! zonal    Gent et Mc William speed (u point) data array
      REAL(8), allocatable :: vgmdta(:,:,:,:) ! meridian Gent et Mc William speed (v point) data array
      REAL(8), allocatable :: wgmdta(:,:,:,:) ! vertical Gent et Mc William speed (w point) data array
#endif

CC----------------------------------------------------------------------

CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC
CC V. DIAGNOSTICS
CC ==============
CC
CC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CC
CC----------------------------------------------------------------------
CC       parameters for the control
CC -------------------------------------------
CC      numnam           : logical unit for namelist

CC      numrst           : logical unit for restart file
CC      lrstar           : boolean term for restart i/o
CC      calendarType     : leap years calendar (0/1)
CC      neuler           : restart euler forward option (0=Euler)



      INTEGER numnam, numrst
      INTEGER calendarType
      INTEGER neuler
      LOGICAL lrstar

      REAL(8), ALLOCATABLE, DIMENSION(:,:) :: DAY_LENGTH


CC----------------------------------------------------------------------
CC         opa output
CC ---------------------------
CC      numout           : logical unit for output print




      INTEGER numout




CC----------------------------------------------------------------------
CC           multitasking

      INTEGER, allocatable :: mtask(:,:)   ! array of multitasking monitoring

CC----------------------------------------------------------------------
CC         massively parallel processors
CC ----------------------------------------------
CC      lwp              : boolean term for mpp output
CC      nummpp           : logical unit for mpp test print

      LOGICAL lwp

      INTEGER nummpp

#if defined key_mpp || defined key_debug 
CC      t3ew           : 3d message passing arrays east-west    
CC      t3we           : 3d message passing arrays west-east
CC      t3ns           : 3d message passing arrays north-south
CC      t3sn           : 3d message passing arrays south-north 
CC      t2ew           : 2d message passing arrays east-west    
CC      t2we           : 2d message passing arrays west-east
CC      t2ns           : 2d message passing arrays north-south
CC      t2sn           : 2d message passing arrays south-north 
C
      REAL(8), allocatable :: t3ns (:,:,:,:), t3sn (:,:,:,:)
      REAL(8), allocatable :: t3ew (:,:,:,:), t3we (:,:,:,:)
      REAL(8), allocatable :: t3ew_my1 (:,:,:,:,:), t3we_my1 (:,:,:,:,:)
      REAL(8), allocatable :: t3sn_my1 (:,:,:,:,:), t3ns_my1 (:,:,:,:,:)
      REAL(8), allocatable :: t2ns (:,:,:)  , t2sn (:,:,:)
      REAL(8), allocatable :: t2ew (:,:,:)  , t2we (:,:,:)

#  else
C     no mpp
#endif

CC----------------------------------------------------------------------
CC Common/comwio/ : io workspace
CC -----------------------------
CC      tabdta         : data auxilary array
CC      tabglo         : global auxilary array

      REAL(8), allocatable :: tabdta(:,:), tabglo(:,:)

#if defined key_mpp|| defined key_debug
CC      tabio          : i/o workspace array

      REAL(8), allocatable :: tabio(:,:,:,:)

#endif

CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC
CC PASSIVE TRACER MODEL
CC =============
CC
CC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CC
CC----------------------------------------------------------------------
CC passive tracers Common
CC ----------------------
CC
#if defined key_passivetrc
CCC F79 Namelist does not accept allocatable

          CHARACTER(LEN=3) :: ctrcnm(jptra)
          CHARACTER(LEN=8) :: ctrcun(jptra)
          CHARACTER(LEN=6) :: mycln
          CHARACTER(LEN=80) ::  ctrcnl(jptra)
          CHARACTER(LEN=3) :: dia(jptra_dia)

CC
CC
CC Common/citctl/ : parameters for the control of passive tracers
CC --------------------------------------------------------------
CC
CC      numnat    : the number of the passive tracer NAMELIST
CC      lutini    : initialisation from FILE or not (NAMELIST)
CC      nutini    : FORTRAN LOGICAL UNIT for initialisation file
CC
      INTEGER numnat

      LOGICAL ::  lutini(jptra)
      INTEGER ::  nutini(jptra)
      REAL(8) ::  maxvalue(jptra)
      LOGICAL :: isCheckLOG
      INTEGER :: jptra_high, jptra_dia_high
      INTEGER :: highfreq_tra(jptra)
      INTEGER :: highfreq_tra_dia(jptra_dia)
      INTEGER freq_ave_phys

CC
CC----------------------------------------------------------------------
CC
CC COMMON/cottrc/ : passive tracers fields (before,now,after)
CC ---------------------------------------
CC      trai      : initial total tracer
CC      trn()     : traceur concentration for actual time step
CC      tra()     : traceur concentration for next time step
CC      trb()     : traceur concentration for before time step
CC

      INTEGER flagSMS_Dyn ! Flag time advance SMS or Dyn
      REAL(8) trai
      REAL(8), allocatable ::  trn(:,:,:,:)
      REAL(8), allocatable ::  tra(:,:,:,:)
      REAL(8), allocatable ::  tra_pp(:,:,:,:)
      REAL(8), allocatable ::  traIO(:,:,:,:)
      REAL(8), allocatable ::  traIO_HIGH(:,:,:,:)
      REAL(8), allocatable ::  snIO(:,:,:) 
      REAL(8), allocatable ::  tnIO(:,:,:) 
      REAL(8), allocatable ::  vatmIO(:,:) 
      REAL(8), allocatable ::  empIO(:,:) 
      REAL(8), allocatable ::  qsrIO(:,:) 
      REAL(8), allocatable ::  unIO(:,:,:) 
      REAL(8), allocatable ::  bblxIO(:,:) 
      REAL(8), allocatable ::  vnIO(:,:,:) 
      REAL(8), allocatable ::  bblyIO(:,:) 
      REAL(8), allocatable ::  wnIO(:,:,:) 
      REAL(8), allocatable ::  avtIO(:,:,:) 
      REAL(8), allocatable ::  tra_ppIO(:,:,:,:)
      REAL(8), allocatable ::  tra_ppIO_HIGH(:,:,:,:)
      REAL(8), allocatable :: tottrn(:,:,:)
      REAL(8), allocatable :: tottrb(:,:,:)

      REAL(8), allocatable ::  tottrnIO(:,:,:) ! matrix for i/o writing(trcdit.F)
      REAL(8), allocatable ::  tottrbIO(:,:,:)
      REAL(8), allocatable ::  totsnIO(:,:,:) 
      REAL(8), allocatable ::  tottnIO(:,:,:) 
      REAL(8), allocatable ::  totvatmIO(:,:) 
      REAL(8), allocatable ::  totempIO(:,:) 
      REAL(8), allocatable ::  totqsrIO(:,:) 
      REAL(8), allocatable ::  totunIO(:,:,:) 
      REAL(8), allocatable ::  totbblxIO(:,:) 
      REAL(8), allocatable ::  totvnIO(:,:,:) 
      REAL(8), allocatable ::  totbblyIO(:,:) 
      REAL(8), allocatable ::  totwnIO(:,:,:) 
      REAL(8), allocatable ::  totavtIO(:,:,:) 
      REAL(8), allocatable ::  tottmaIO(:,:,:) 


      REAL(8), allocatable ::  trb(:,:,:,:)
      REAL(8), allocatable ::  buf(:,:,:)
      REAL(8), allocatable ::  buf2(:,:)
      INTEGER, allocatable, dimension(:) :: highfreq_table,highfreq_table_dia

CC----------------------------------------------------------------------
CC
CC COMMON /cot3ad/ non-centered advection scheme (smolarkiewicz)
CC -------------------------------------------------------------
CC      rsc         : tuning coefficient for anti-diffusion (NAMELIST)
CC      rtrn        : value for truncation (NAMELIST)

      REAL(8) rsc,rtrn


CC----------------------------------------------------------------------
CC
CC COMMON /cit3ad/ non-centered advection scheme (smolarkiewicz)
CC -------------------------------------------------------------
CC      crosster    : logical if true computes crossterms (NAMELIST) 
CC      lhdf        : logical if true CALL trchdf (NAMELIST) 
CC      ncor        : number of corrective phases (NAMELIST)
CC      ndttrc      : frequency of step on passive tracers (NAMELIST)

      INTEGER ncor,ndttrc
      LOGICAL lhdf,crosster



CC----------------------------------------------------------------------
CC
CC COMMON/citrst/ : passive tracers restart (input and output)
CC -----------------------------------------------------------
CC      nutwrs    : output FILE for passive tracers restart
CC      lrsttr    : boolean term for restart i/o for passive tracers 
CC                  (NAMELIST)
CC      nutrst    : logical unit for restart FILE for passive tracers
CC      nrsttr    : control of the time step ( 0 or 1 ) for pass. tr.
CC                  (NAMELIST)

      LOGICAL lrsttr
      INTEGER nutwrs,nutrst,nrsttr

COMMON/cotiso/ : isopycnal sheme for passive tracers 
CC -----------------------------------------------------------
CC      ahtrb0    : background diffusivity coefficient (m2/s)
CC                  for passive tracer
CC      trcrat    : ratio between passive and active tracer coeff
CC                  for diffusion
CC      ahtrc0    : horizontal eddy diffusivity for passive tracers (m2/s)
CC    aeivtr0   : eddy induced velocity coefficient (m2/s)

      REAL(8) ahtrb0,trcrat,ahtrc0,aeivtr0

CC
CC----------------------------------------------------------------------
CC
CC      information for outputs
CC ------------------------------------------------------------------
CC
CC      nwritetrc: time step frequency for concentration outputs (NAMELIST)
CC
CC    netcdf files and index common
CC

      INTEGER nwritetrc
      REAL(8) djulian

#    if defined key_trc_diaadd
CC----------------------------------------------------------------------
CC
CC          additional 2D/3D outputs
CC ------------------------------------------------------------------
CC
CC      ctrc3d    : 3d output field name (NAMELIST)
CC      ctrc3l    : 3d output field long name (NAMELIST)
CC      ctrc3u    : 3d output field unit (NAMELIST)
CC      ctrc2d    : 2d output field name (NAMELIST)
CC      ctrc2l    : 2d output field long name (NAMELIST)
CC      ctrc2u    : 2d output field unit (NAMELIST)
CC      trc3d     : additional 3d outputs
CC      trc2d     : additional 2d outputs
CC
      CHARACTER(LEN=8), allocatable :: ctrc3d(:),ctrc2d(:)
      CHARACTER(LEN=8), allocatable :: ctrc3u(:),ctrc2u(:)
      CHARACTER(LEN=80), allocatable :: ctrc3l(:),ctrc2l(:)
      REAL(8), allocatable :: trc3d(:,:,:,:), trc2d(:,:,:)

CC
CC    netcdf files and index common
CC
CC      nwriteadd: frequency of additional arrays outputs (NAMELIST)
CC      nitd     : id for additional array output FILE
CC      ndepitd  : id for depth mesh
CC      nhoritd  : id for horizontal mesh
CC
      INTEGER nwriteadd,nitd,ndepitd,nhoritd
#    endif

#    if defined key_trc_diatrd
CC----------------------------------------------------------------------
CC
CC COMMON/cottrd/ : non conservative trends (biological, ...)
CC ------------------------------------------------------------------
CC
CC      luttrd    : large trends diagnostic to WRITE or not (NAMELIST)
CC
      LOGICAL, allocatable ::   luttrd(:)
CC
CC    dynamical trends
CC
CC    trtrd()   : trends of the tracer equations
CC           1 : X advection
CC           2 : Y advection
CC           3 : Z advection
CC           4 : X diffusion
CC           5 : Y diffusion
CC           6 : Z diffusion
CC           7 : X gent velocity
CC           8 : Y gent velocity
CC           9 : Z gent velocity

      REAL(8), allocatable ::  trtrd(:,:,:,:,:)
CC
CC    netcdf files and index common
CC
CC      nwritetrd: frequency for dynamical trends output (NAMELIST)
CC                 one per tracer

CC
      INTEGER nwritetrd,nhorit6(:)
#    endif 
CC
CC----------------------------------------------------------------------
CC
CC COMMON/cotrda/ : passive tracers DATA READ and at given time_step
CC -----------------------------------------------------------------
CC      nclimr    : switch for passive tracers initialization
CC      ncontr    : variable for time interpolation
CC      numtr1    : LOGICAL UNIT for passive tracers DATA
CC      numtr2    : LOGICAL UNIT for passive tracers DATA created IF
CC                  interpolation is needed (ninttr=1)
CC      nlectr    : switch for reading once
CC      ninttr    : switch for interpolation on model grid
CC      nmldmptr  : : = 0/1/2 type of damping in the mixed layer
CC      trcdat()  : passive tracers DATA array for two value
CC                  needed for time interpolation
CC
CC      trdta()   : passive tracers DATA at given time-step
CC      restotr() : array of restoring coeff. for passive tracers
CC
      INTEGER nclimr,ncontr,numtr1,numtr2,nlectr,ninttr,nmldmptr
CC
CC
#    if defined key_trc_dmp 

      REAL(8), allocatable ::  trdta(:,:,:,:),trcdat(:,:,:,:,:)
!      REAL(8), allocatable ::  restoglo(:,:,:), idxtglo(:,:,:)
      INTEGER(4), allocatable ::  idxt(:,:,:),idxt2glo(:,:,:,:)

CC
#    else
CC      no passive tracers DATA at given time step
#    endif
CC
CC
CC
CC    fotoperiod formulation
CC    photop      : Fotoperiod formulation if false daylength is 24 h
      LOGICAL photop
CC    atlantic buffer biology activation
      LOGICAL atlantic_bfm

#    if defined key_trc_bfm
CC
CC      vsed      : sedimentation speed (NAMELIST)
CC      remine()  : remineralisation trend
CC
CC      tmaxr     : maximum coefficient for passive tracer damping (NAMELIST)
CC      tminr     : minimum coefficient for passive tracer damping (NAMELIST)
CC      remdmp()  : damping coefficient of passive tracers (depth dependant)


      REAL(8) vsed,vsed_dia,tmaxr,tminr,bottom_flux
      REAL(8), allocatable :: remine(:,:,:)
      REAL(8), allocatable :: remdmp(:,:)


CC----------------------------------------------------------------------
CC
CC COMMON/cotopt/ : optical parameters
CC -----------------------------------
CC
CC      xze       : euphotic layer depth
CC      xpar      : par (photosynthetic available radiation)
CC      xkr0      : water coefficient absorption in red (NAMELIST)
CC      xkg0      : water coefficient absorption in green (NAMELIST)
CC      xkrp      : pigment coefficient absorption in red (NAMELIST)
CC      xkgp      : pigment coefficient absorption in green (NAMELIST)
CC      xlr       : exposant for pigment absorption in red (NAMELIST)
CC      xlg       : exposant for pigment absorption in green (NAMELIST)
CC      rpig      : chla/chla+phea ratio (NAMELIST)
CC
      REAL(8) xkr0,xkg0,xkrp,xkgp,xlr,xlg,rpig

      REAL(8), allocatable :: xze(:,:)
      REAL(8), allocatable :: xpar(:,:,:)
CC----------------------------------------------------------------------
CC

#    endif
CC
CC
#else
CC
CC no passive tracer COMMON specification
CC
#endif
CC
CC----------------------------------------------------------------------

#ifdef key_trahdfeiv
      REAL(8), allocatable, dimension(:,:,:) :: zaeiu,zaeiv, zaeiw
#endif
#ifdef key_s_coord
      REAL(8), allocatable ::hbattglo(:,:)    , hbatuglo(:,:)
      REAL(8), allocatable ::hbatvglo(:,:)    , hbatfglo(:,:)
#endif
      INTEGER ncpu



      CONTAINS
! *******************************************************************
      subroutine myalloc_DYN()
      mem_all = 0.
#ifdef key_trahdfeiv
       allocate(zaeiu   (jpiglo,jpjglo,jpk))
       allocate(zaeiv   (jpiglo,jpjglo,jpk))
       allocate(zaeiw   (jpiglo,jpjglo,jpk))
        mem_all = 8*(3*jpiglo*jpjglo*jpk)
#endif
      END subroutine myalloc_DYN
! *******************************************************************
      subroutine myalloc_DR()
      mem_all = 0
#ifdef key_s_coord
      allocate(hbattglo(jpiglo,jpjglo))
      allocate(hbatuglo(jpiglo,jpjglo))
      allocate(hbatvglo(jpiglo,jpjglo))
      allocate(hbatfglo(jpiglo,jpjglo))
      mem_all = rea_len*4*jpiglo*jpjglo
#endif
      END subroutine myalloc_DR
! *******************************************************************

       subroutine myalloc_BFM()

       allocate(BFMpoints(3,NBFMPOINTS))
       mem_all = mem_all + 3*NBFMPOINTS*4
       end subroutine myalloc_BFM



       subroutine alloc_tot()

C I. MODEL SET UP
C ===============
C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CC
CC----------------------------------------------------------------------
CC Common/comdom/  : domain parameters
CC -----------------------------------
CC      nizoom, njzoom   : (i,j) indexes for the domain left bottom zoom
CC      nperio           : lateral boundary condition type 
CC      nimpp,njmpp      : (i,j) indexes for mpp-subdomain left bottom
CC      nreci,nrecj      : overlap region in i and j
CC      nproc            : number for local processor
CC      narea            : number for local area
CC      nbondi, nbondj   : mark of i- and j-direction local boundaries
CC      nlci, nlcj       : i, j dimensions of the local subdomain
CC      nldi, nlei,      : first and last indoor i- and j-indexes
CC      nldj, nlej   
CC      noea, nowe,      : index of the local neighboring processors in
CC      noso, nono         east, west, south and north directions
CC      nimppt(),njmppt(): i-, j-indexes for each processor
CC      nlcit(), nlcjt() : dimensions of every subdomain
CC      nldit(), nldjt() : first, last indoor index for each i-domain
CC      nleit(), nlejt() : first, last indoor index for each j-domain
CC      mindi(), mindj() : indexes array of the subdomain
CC
CC      INTEGER mindi(jpi), mindj(jpj)
CC    INTEGER nimppt(jpnij), njmppt(jpnij), nlcit(jpnij), nlcjt(jpnij)
CC    INTEGER  nldit(jpnij),  nldjt(jpnij), nleit(jpnij), nlejt(jpnij)
CCC F79 18 01 2006 to calculate total memory allocated
      ncpu = 1
      mem_all = 0.0
      int_len = 4.0
      rea_len = 8.0
      allocate(mindi(jpi))
      allocate(mindj(jpj))
      allocate(nimppt(jpnij))
      allocate(njmppt(jpnij))
      allocate(nlcit(jpnij))
      allocate(nlcjt(jpnij))
      allocate(nldit(jpnij))
      allocate(nldjt(jpnij))
      allocate(nleit(jpnij))
      allocate(nlejt(jpnij))
      mem_all = mem_all + int_len*(jpi+jpj+8*jpnij) 
C     write(*,*) 'Mem_allocated:', mem_all
CC
CC----------------------------------------------------------------------
CC Common/comcst/ : physical constants
CC -----------------------------------
CC      nbiss()          : number of days per month (leap-year)
CC      nobis()          : number of days per month (365 days a year)
CC      rpi              : pi
CC      rday             : day
CC      rsiyea           : sideral year
CC      rsiday           : sideral day  
CC      omega            : earth rotation parameter
CC      ra               : earth radius (meter)
CC      g                : gravity
CC      rtt              : fusion point
CC      rt0              : zero celsius
CC      rauw             : density of pure water kg/m3
CC      rad              : conversion coeff. from degre into radian
CC      rcs              : conversion coeff. from millimeter per day
CC                         into meter per second
CC      raajj            : number of days in one year
CC      raamo            : number of months in one year
CC      rjjhh            : number of hours in one day
CC      rhhmm            : number of minutes in one hour
CC      rmmss            : number of seconds in one minute
CC      raass            : number of seconds in one year
CC      rmoss            : number of seconds in one month
CC      rjjss            : number of seconds in one day
C
C
CC
CC
CC----------------------------------------------------------------------
CC Common/my_temp/ : cronometer for function calls
CC -----------------------------------

CC----------------------------------------------------------------------
CC Common/comeos/ : ocean physical parameters (equation of state, ...)
CC ------------------------------------------
CC    neos         : flag of the type of equation of state used
CC      rau0             : reference volumic mass of the ocean (kg/m3)
CC      ralpha, rbeta    : thermique and haline expension coef. used
CC               for linear equation of state (neos=1 or 2)
CC      rcp              : ocean specific heat
CC      ro0cpr           : = 1. / ( rau0 * rcp )
C
CC

      allocate(rdttra(jpk))
      mem_all = mem_all + int_len*jpk 
C     write(*,*) 'Mem_allocated:', mem_all
CC
CC
CC----------------------------------------------------------------------
CC     horizontal curvilinear coordinate and scale factors
CC ---------------------------------------------------------------------
CC      nummsh           : logical unit of coordinate/meshes output file
CC      ngrid            : flag to read or compute coordinates
CC      nmsh             : flag to write coordinates and meshes file
CC      njpi, njpj       : read dimensions
CC      rra              : read earth radius
CC      glamt()          : longitude of t-point (degre)
CC      glamu()          : longitude of u-point (degre)
CC      glamv()          : longitude of v-point (degre)
CC      glamf()          : longitude of f-point (degre)
CC      gphit()          : latitude  of t-point (degre)
CC      gphiu()          : latitude  of u-point (degre)
CC      gphiv()          : latitude  of v-point (degre)
CC      gphif()          : latitude  of f-point (degre)
CC      e1t,e2t()        : horizontal scale factors at t-point (m)
CC      e1u,e2u()        : horizontal scale factors at u-point (m)
CC      e1v,e2v()        : horizontal scale factors at v-point (m)
CC      e1f,e2f()        : horizontal scale factors at f-point (m)
CC    ff()             : coriolis factor (2.*omega*sin(yphi) ) (s-1)
C
CC      REAL(8) glamt(jpi,jpj), glamu(jpi,jpj), glamv(jpi,jpj),glamf(jpi,jpj)
CC      REAL(8) gphit(jpi,jpj), gphiu(jpi,jpj), gphiv(jpi,jpj),gphif(jpi,jpj)
CC      REAL(8) e1t(jpi,jpj), e1u(jpi,jpj), e1v(jpi,jpj), e1f(jpi,jpj)
CC      REAL(8) e2t(jpi,jpj), e2u(jpi,jpj), e2v(jpi,jpj), e2f(jpi,jpj)
CC      REAL(8) ff(jpi,jpj)
        allocate(totglamt(jpiglo,jpjglo))
        allocate(totgphit(jpiglo,jpjglo))
       !allocate(glamt(jpi,jpj))
        allocate(glamu(jpi,jpj)) 
        allocate(glamv(jpi,jpj))
        allocate(glamf(jpi,jpj))
        allocate(gphit(jpi,jpj))
        allocate(gphiu(jpi,jpj)) 
        allocate(gphiv(jpi,jpj))
        allocate(gphif(jpi,jpj))
        allocate(e1t(jpi,jpj)) 
        allocate(e1u(jpi,jpj)) 
        allocate(e1v(jpi,jpj)) 
        allocate(e1f(jpi,jpj))
        allocate(e2t(jpi,jpj)) 
        allocate(e2u(jpi,jpj)) 
        allocate(e2v(jpi,jpj)) 
        allocate(e2f(jpi,jpj))
        allocate(ff (jpi,jpj))
        mem_all = mem_all + rea_len*(2*jpiglo*jpjglo+15*jpi*jpj)
C     write(*,*) 'Mem_allocated:', mem_all

#if defined key_off_degrad
        allocate(facvol(jpi,jpj,jpk))
        mem_all = mem_all + rea_len*(jpi*jpj*jpk) 
C     write(*,*) 'Mem_allocated:', mem_all
#endif
C
CC
CC----------------------------------------------------------------------
CC Common/comcoz/  : vertical coordinate and scale factors
CC -------------------------------------------------------
CC      numhba           : logical unit for bathymetry file (meters)
C
C
CC
#ifdef key_s_coord
CC                   s-coordinate ('key_s_coord')
CC                   ----------------------------
CC      hbatt(), hbatu() : ocean depth at the vertical of  t-, u-, v-
CC      hbatv(), hbatf()   and f-point (m)
CC      gsigt(), gsigw() : model level depth coefficient at t-, w-levels
CC      gsi3w()          : model level depth coefficient at w-level
CC                         defined as the sum of e3w scale factors
CC      esigt(), esigw() : vertical scale factor coef. at t-, w-levels
C
       allocate(hbatt(jpi,jpj))
       allocate(hbatu(jpi,jpj)) 
       allocate(hbatv(jpi,jpj))
       allocate(hbatf(jpi,jpj))
       allocate(gsigt(jpk)) 
       allocate(gsigw(jpk))
       allocate(gsi3w(jpk))
       allocate(esigt(jpk))
       allocate(esigw(jpk))
       mem_all = mem_all + rea_len*(4*jpi*jpj+5*jpk) 
      write(*,*) 'Mem_allocated:', mem_all
C
CC
#  else
CC                  z-coordinate (default option)
CC                  ------------------------------
CC      gdept(), gdepw() : depth of t- and w-points (m)
CC      e3t(), e3w()     : vertical scale factors at t- and w-points (m)
CC
      allocate(gdept(jpk)) 
      allocate(gdepw(jpk))
      allocate(  e3t(jpk))
      allocate(  e3w(jpk))
       mem_all = mem_all + rea_len*(4*jpk) 
C     write(*,*) 'Mem_allocated:', mem_all
C
#endif
CC
CC----------------------------------------------------------------------
CC Common/comask/  : masks, bathymetry
CC -----------------------------------
CC      numbat           : logical unit for bathymetry file
CC      mbathy()         : number of ocean level (=0, 1, ... , jpk-1)
CC      tmask(), umask() : land/ocean mask at t-, u-, v- and f-points
CC      vmask(), fmask()
C
      allocate(mbathy(jpi,jpj))
C
C
      allocate(tmask(jpi,jpj,jpk))
      allocate(fmask(jpi,jpj,jpk))
      allocate(umask(jpi,jpj,jpk))
      allocate(vmask(jpi,jpj,jpk))
      mem_all = mem_all + int_len*jpi*jpj + rea_len*4*jpi*jpj*jpk  


C     write(*,*) 'Mem_allocated:', mem_all
C
CC
CC
CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC
CC II. DYNAMICS AND TRACERS
CC ========================
CC
CC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CC
CC----------------------------------------------------------------------
CC Common/combef/  : previous fields (before)
CC -----------------------------------------
CC    ub(), vb()       : horizontal velocity (m s-1)
C
      allocate(ub(jpi,jpj,jpk)) 
      allocate(vb(jpi,jpj,jpk))
      mem_all = mem_all +  rea_len*2*jpi*jpj*jpk  
C     write(*,*) 'Mem_allocated:', mem_all
C
CC
CC----------------------------------------------------------------------
CC Common/comnow/  : present fields (now)
CC -------------------------------------
CC    un(), vn(), wn() : horizontal and vertical velocity (m s-1)
CC    hdivn()          : horizontal divergence
CC    tn(),   sn()     : pot. temperature (celsius), salinity (psu)
CC      rdn()            : in situ density anomalie rdn=(rho-rau0)/rau0
CC                         (no units)
CC    rhopn()          : potential volumic mass (kg m-3)
CC      bn2n()           : brunt-vaisala frequency (s-2)
CC
       allocate(un(jpi,jpj,jpk))
       allocate(vn(jpi,jpj,jpk))
       allocate(wn(jpi,jpj,jpk))
       allocate(hdivn(jpi,jpj,jpk))
       allocate(tn(jpi,jpj,jpk))
       allocate(sn(jpi,jpj,jpk))
       allocate(rdn(jpi,jpj,jpk))
       allocate(rhopn(jpi,jpj,jpk))
       allocate(rho(jpi,jpj,jpk))
       allocate(bn2n(jpi,jpj,jpk))
       mem_all = mem_all +  rea_len*9*jpi*jpj*jpk  
C     write(*,*) 'Mem_allocated:', mem_all
C
CC
CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC III. OCEAN PHYSICS
CC ==================
CC
CC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CC
CC----------------------------------------------------------------------
CC Common/comhdt/ : lateral diffusivity (tracers)
CC ------------------------------------
CC      aht0             : lateral diffusivity coefficient    (namelist)
CC      ahtb0            : background diffusivity coefficient (m2/s)
CC      ahtu(), ahtv()   : lateral diffusivity coef. at u-, v-, w- t-pts
CC      ahtw(), ahtt()     (harmonic operator: no rotation, use of u-
CC                          and v-points; rotation, use of u-, v- w-pts)
CC                         (biharmonic operator: rotation or not, use of
CC                          t-point only)
CC                         (the arrays used are 3D, 2D, 1D or 0D depen-
CC                          ding on 'key_trahdfcoef.d' )
C
CCC Paolo 22/4/2004 i killed all options
#if defined key_trahdfbilap
# else
#  if defined key_trahdfcoef3d
C    allocate(ahtu(jpi,jpj,jpk))
C    allocate(ahtv(jpi,jpj,jpk))
C    allocate(ahtw(jpi,jpj,jpk))
#   elif key_trahdfcoef2d
C    allocate(ahtu(jpi,jpj)) 
C    allocate(ahtv(jpi,jpj))
C    allocate(ahtw(jpi,jpj))
#   elif key_trahdfcoef1d
C    allocate(ahtu(jpk)) 
C    allocate(ahtv(jpk))
C    allocate(ahtw(jpk))
#   else
#     if defined key_off_degrad
C     allocate(ahtu(jpi,jpj,jpk)) 
C     allocate(ahtv(jpi,jpj,jpk))
C     allocate(ahtw(jpi,jpj,jpk))
#     else
#     endif
#  endif
#endif
CCC  Paolo 22/4/2004 new conditions 1D case attention
      allocate(ahtu(jpk))
      allocate(ahtv(jpk))
      allocate(ahtw(jpk))
      allocate(ahtt(jpk))
      mem_all = mem_all +  rea_len*4*jpk  
C     write(*,*) 'Mem_allocated:', mem_all

CC
CC----------------------------------------------------------------------
CC Common/comdhd/ : direction of lateral diffusion (momentum & tracers) 
CC -----------------------------------------------
CC Slopes ('key_trahdfiso' or 'key_s_coord' defined)
CC      uslp(), wslpi()  : i-direction slope at u-, w-points
CC      vslp(), wslpj()  : j-direction slope at v-, w-points
CC    wk1(), wk2()     : work array used in rotated biharmonic ope-
CC    wk3(), wk4()       rator on momentum and/or tracers
C
#if defined key_trahdfiso || defined key_trahdfgeop || defined key_debug
      allocate(uslp(jpi,jpj,jpk)) 
      allocate(vslp(jpi,jpj,jpk))
      allocate(wslpi(jpi,jpj,jpk))
      allocate(wslpj(jpi,jpj,jpk))
       mem_all = mem_all +  rea_len*4*jpi*jpj*jpk  
C     write(*,*) 'Mem_allocated:', mem_all
#  if defined key_trahdfbilap
      allocate(wk1(jpi,jpj,jpk))
      allocate(wk2(jpi,jpj,jpk))
      allocate(wk3(jpi,jpj,jpk))
      allocate(wk4(jpi,jpj,jpk))
       mem_all = mem_all +  rea_len*4*jpi*jpj*jpk  
C     write(*,*) 'Mem_allocated:', mem_all
#   else
#  endif
C
# else
C     diffusion along model level surfaces: no arrays
#endif
CC
CC----------------------------------------------------------------------
CC Common/comeiv/ : eddy induced velocity parameterization
CC -------------------------------------------------------
CC      aeiv0            : eddy induced velocity coefficient (m2/s)
CC      aeiu(), aeiv()   : eddy induced velocity coef. at u-, v- and
CC      aeiw()             w-points (m2/s). The arrays used are 3D
CC                         2D, 1D or 0D, depending on 'key_trahdfcoef.d'
C
#if defined key_trahdfeiv
#  if defined key_trahdfcoef3d
       allocate(aeiu(jpi,jpj,jpk))
       allocate(aeiv(jpi,jpj,jpk))
       allocate(aeiw(jpi,jpj,jpk))
       mem_all = mem_all +  rea_len*3*jpi*jpj*jpk  
C     write(*,*) 'Mem_allocated:', mem_all
#   elif key_trahdfcoef2d
       allocate(aeiu(jpi,jpj)) 
       allocate(aeiv(jpi,jpj))
       allocate(aeiw(jpi,jpj))
       mem_all = mem_all +  rea_len*3*jpi*jpj*jpk  
C     write(*,*) 'Mem_allocated:', mem_all
#   elif key_trahdfcoef1d
       allocate(aeiu(jpk)) 
       allocate(aeiv(jpk))
       allocate(aeiw(jpk))
       mem_all = mem_all +  rea_len*3*jpi*jpj*jpk  
C     write(*,*) 'Mem_allocated:', mem_all
#   else
#  endif
# else
#     if defined key_off_degrad
      allocate(aeiu(jpi,jpj,jpk))
      allocate(aeiv(jpi,jpj,jpk)) 
      allocate(aeiw(jpi,jpj,jpk))
       mem_all = mem_all +  rea_len*3*jpi*jpj*jpk  
C     write(*,*) 'Mem_allocated:', mem_all
#     else
#     endif
#endif
C
CC-------------------------------------------------------------------
CC Common/comgm/ : eddy induced velocities (Gent et Mc William)
CC---------------------------------------------------------------
CC       ugm  : eddy induced velocity
CC       vgm
CC       wgm
C
#if defined key_trahdfeiv
       allocate(ugm(jpi,jpj,jpk))
       allocate(vgm(jpi,jpj,jpk))
       allocate(wgm(jpi,jpj,jpk))
       mem_all = mem_all +  rea_len*3*jpi*jpj*jpk  
C     write(*,*) 'Mem_allocated:', mem_all
#endif
CC
CC----------------------------------------------------------------------
CC Common/comzdf/ : vertical diffusion
CC -----------------------------------
CC    avt0             : vertical viscosity and diff. coef. (namelist)
CC      ntrbbl         : bottom boundary layer (namelist)
CC      atrbbl         : bottom boundary layer diffusivity (namelist)
CC    avt()             : vertical diffusivity coeff. at w-point
CC    avtb()           : background profile of avm and avt
C
C
       allocate(avt (jpi,jpj,jpk))
       allocate(avtb(jpk))
       allocate(bblx(jpi,jpj))
       allocate(bbly(jpi,jpj))
       mem_all = mem_all +  rea_len*(jpi*jpj*jpk+ jpk +2*jpi*jpj) 


CC
CC
CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC
CC IV. SURFACE FORCING AND DATA
CC ============================
CC
CC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CC
CC----------------------------------------------------------------------
CC Common/comtau/ : surface wind stress at givem time_step
CC -------------------------------------------------------
CC    taux(), tauy()   : wind stress components in (i,j) referential
CC    tauxg(), tauyg() : zonal and meridian wind stress component used
CC                       in output (geographical referential)
C
       allocate(taux  (jpi,jpj))
       allocate(tauy  (jpi,jpj))
       allocate(vatm  (jpi,jpj))
       allocate(freeze(jpi,jpj))
       mem_all = mem_all +  rea_len*(4*jpi*jpj) 
C     write(*,*) 'Mem_allocated:', mem_all
C
CC
CC----------------------------------------------------------------------
CC Common/comflx/ : surface fluxes
CC -------------------------------
CC      qt()             : total surface heat flux (w m-2)
CC      q()              : surface heat flux (w m-2)
CC      emp()            : evaporation minus precipitation (mm day-1)
CC      runoff()         : annual run off (mm/day)

       allocate(qt    (jpi,jpj))
       allocate(q     (jpi,jpj))
       allocate(emp   (jpi,jpj))
       allocate(runoff(jpi,jpj))
       mem_all = mem_all +  rea_len*(4*jpi*jpj) 

CC----------------------------------------------------------------------
CC
CC COMMON/comqsr/ : penetrative solar radiation
CC --------------------------------------------
CC      qsr()     : solar radiation (w m-2)
CC
       allocate(qsr(jpi,jpj))
       mem_all = mem_all +  rea_len*(jpi*jpj) 


       allocate(udta   (jpi,jpj,jpk,2))
       allocate(vdta   (jpi,jpj,jpk,2))
       allocate(wdta   (jpi,jpj,jpk,2))
       allocate(avtdta (jpi,jpj,jpk,2))
       allocate(flxdta (jpi,jpj,jpflx,2))
       allocate(flx    (jpi,jpj,jpflx))
       allocate(bblxdta(jpi,jpj,2))
       allocate(bblydta(jpi,jpj,2))
       allocate(zmxldta(jpi,jpj,2))
       allocate(hmld(jpi,jpj))
       mem_all = mem_all +  rea_len*(17*jpi*jpj+ jpi*jpj*jpflx) 


#if defined key_off_tra
CC----------------------------------------------------------------------
CC
CC Common/comdyt/ : dynamic data selection
CC ---------------------------------------
CC      tdta(),sdta()  : temperature and salinity data array
CC
       allocate(tdta(jpi,jpj,jpk,2))
       allocate(sdta(jpi,jpj,jpk,2))
       mem_all = mem_all +  rea_len*(4*jpi*jpj*jpk) 
C     write(*,*) 'Mem_allocated:', mem_all
CC
#endif
CC
#if defined key_trahdfiso
CC----------------------------------------------------------------------
CC
CC Common/comdyi/ : dynamic data selection
CC ---------------------------------------
CC      uslpdta()      : zonal slope (u point) data array
CC      vslpdta()      : meridian slope (v point) data array
CC      wslpidta()     : zonal slope (t point) data array
CC      wslpjdta()     : meridian slope (t point) data array
CC
       allocate(uslpdta(jpi,jpj,jpk,2))
       allocate(vslpdta(jpi,jpj,jpk,2))
       allocate(wslpidta(jpi,jpj,jpk,2))
       allocate(wslpjdta(jpi,jpj,jpk,2))
       mem_all = mem_all +  rea_len*(8*jpi*jpj*jpk) 
C     write(*,*) 'Mem_allocated:', mem_all
CC
#endif
CC
#if defined key_trahdfeiv
CC
CC Common/comdyd/ : dynamic data selection
CC ---------------------------------------
CC      ugmdta()      : zonal Gent et Mc William speed (u point) data array
CC      vgmdta()      : meridian Gent et Mc William speed (v point) data array
CC      wgmdta()      : vertical Gent et Mc William speed (w point) data array
CC
       allocate(ugmdta(jpi,jpj,jpk,2))
       allocate(vgmdta(jpi,jpj,jpk,2))
       allocate(wgmdta(jpi,jpj,jpk,2))
       mem_all = mem_all +  rea_len*(6*jpi*jpj*jpk) 
C     write(*,*) 'Mem_allocated:', mem_all
#endif
CC
CC----------------------------------------------------------------------

CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC
CC V. DIAGNOSTICS
CC ==============


CC----------------------------------------------------------------------
CC           multitasking

        allocate(mtask(3,jpcpu))
       mem_all = mem_all +  int_len*(3*jpcpu) 

CC----------------------------------------------------------------------
CC        massively parallel processors
CC ----------------------------------------------
CC      lwp              : boolean term for mpp output
CC      nummpp           : logical unit for mpp test print
CC
#if defined key_mpp || defined key_debug 
CC      t3ew()           : 3d message passing arrays east-west    
CC      t3we()           : 3d message passing arrays west-east
CC      t3ns()           : 3d message passing arrays north-south
CC      t3sn()           : 3d message passing arrays south-north 
CC      t2ew()           : 2d message passing arrays east-west    
CC      t2we()           : 2d message passing arrays west-east
CC      t2ns()           : 2d message passing arrays north-south
CC      t2sn()           : 2d message passing arrays south-north 
C
       allocate(t3ns(jpi,jprecj,jpk,2))
       allocate(t3sn(jpi,jprecj,jpk,2))
       allocate(t3ew(jpj,jpreci,jpk,2)) 
       allocate(t3we(jpj,jpreci,jpk,2))
       allocate(t3ew_my1(jpj,jpreci,jpk,jptra,2)) 
       allocate(t3we_my1(jpj,jpreci,jpk,jptra,2))
       allocate(t3sn_my1(jpi,jpreci,jpk,jptra,2)) 
       allocate(t3ns_my1(jpi,jpreci,jpk,jptra,2))
       allocate(t2ns(jpi,jprecj,2))  
       allocate(t2sn(jpi,jprecj,2))
       allocate(t2ew(jpj,jpreci,2))
       allocate(t2we(jpj,jpreci,2))
       mem_all = mem_all +  rea_len*(4*jpi*jprecj*jpk+4*jpj*jpreci*jpk 
     .           + 4*jpi*jprecj+4*jpj*jpreci) 
C     write(*,*) 'Mem_allocated:', mem_all
C
#  else
C     no mpp
#endif
CC
CC----------------------------------------------------------------------
CC          io workspace
CC -----------------------------
CC      tabdta()         : data auxilary array
CC      tabglo()         : global auxilary array
C
       allocate(tabdta(jpidta,jpjdta))
       allocate(tabglo(jpiglo,jpjglo))
       mem_all = mem_all +  rea_len*(jpidta*jpjdta+jpiglo*jpjglo) 
C     write(*,*) 'Mem_allocated:', mem_all
CC
#if defined key_mpp|| defined key_debug
CC      tabio()          : i/o workspace array
C
      allocate(tabio(jpi,jpj,jpnij,jpkmod))
       mem_all = mem_all +  rea_len*(jpi*jpj*jpnij*jpkmod) 
C     write(*,*) 'Mem_allocated:', mem_all
CC
#  else
#endif
CC
CC
CC
CC
CC>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
CC
CC PASSIVE TRACER MODEL
CC =============
CC
CC<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
CC

#if defined key_passivetrc
CCC F79 Namelist does not accept allocatable
CCC         allocate(ctrcnm(jptra))
CCC         allocate(ctrcun(jptra))
CCC         allocate(ctrcnl(jptra))
CC
CC
CC Common/citctl/ : parameters for the control of passive tracers
CC --------------------------------------------------------------
CC
CC      numnat    : the number of the passive tracer NAMELIST
CC      lutini    : initialisation from FILE or not (NAMELIST)
CC      nutini    : FORTRAN LOGICAL UNIT for initialisation file
CC
CCC F79 Namelist does not accept allocatable
CCC       allocate(lutini(jptra))
CCC       allocate(nutini(jptra))
CC
CC----------------------------------------------------------------------
CC
CC COMMON/cottrc/ : passive tracers fields (before,now,after)
CC ---------------------------------------
CC      trai      : initial total tracer
CC      trn()     : traceur concentration for actual time step
CC      tra()     : traceur concentration for next time step
CC      trb()     : traceur concentration for before time step
CC
CCC Paolo 30/4/2004 Flag time advance SMS or Dyn

       allocate(trn(jpi,jpj,jpk,jptra))
       allocate(tra(jpi,jpj,jpk,jptra))
       allocate(tra_pp(jpi,jpj,jpk,jptra_dia))
       if (lwp) allocate(tottrn(jpiglo, jpjglo, jpk))
       if (lwp) allocate(tottrb(jpiglo, jpjglo, jpk))
       allocate(traIO(jpi,jpj,jpk,jptra)) 
       allocate(snIO(jpi,jpj,jpk)) 
       allocate(tnIO(jpi,jpj,jpk)) 
       allocate(vatmIO(jpi,jpj)) 
       allocate(empIO(jpi,jpj)) 
       allocate(qsrIO(jpi,jpj)) 
       allocate(unIO  (jpi,jpj,jpk))
       allocate(bblxIO(jpi,jpj)) 
       allocate(vnIO  (jpi,jpj,jpk))
       allocate(bblyIO(jpi,jpj)) 
       allocate(wnIO  (jpi,jpj,jpk))
       allocate(avtIO (jpi,jpj,jpk))
       allocate(buf   (jpi,jpj,jpk))
       allocate(buf2   (jpi,jpj))
       allocate(tra_ppIO(jpi,jpj,jpk,jptra_dia))
       allocate(traIO_HIGH(   jpi,jpj,jpk,jptra_HIGH))
       allocate(tra_ppIO_HIGH(jpi,jpj,jpk,jptra_dia_HIGH))


       mem_all = mem_all +  rea_len*(3*jpi*jpj*jpk*jptra
     .          + (jpi*jpj*jpk)*(jptra_HIGH + jptra_dia_HIGH)
     .          + 13*jpi*jpj*jpk+2*jpiglo* jpjglo* jpk
     .          + 5*jpi*jpj)


       if (lwp) allocate(tottrnIO(jpiglo,jpjglo,jpk))
       if (lwp) allocate(tottrbIO(jpiglo,jpjglo,jpk))
       if (lwp) allocate(totsnIO (jpiglo,jpjglo,jpk))
       if (lwp) allocate(tottnIO (jpiglo,jpjglo,jpk))
       if (lwp) allocate(totvatmIO(jpiglo,jpjglo))
       if (lwp) allocate(totempIO(jpiglo,jpjglo))
       if (lwp) allocate(totqsrIO(jpiglo,jpjglo))
       if (lwp) allocate(totunIO(jpiglo,jpjglo,jpk))
       if (lwp) allocate(totbblxIO(jpiglo,jpjglo))
       if (lwp) allocate(totvnIO(jpiglo,jpjglo,jpk))
       if (lwp) allocate(totbblyIO(jpiglo,jpjglo))
       if (lwp) allocate(totwnIO(jpiglo,jpjglo,jpk))
       if (lwp) allocate(totavtIO(jpiglo,jpjglo,jpk))
       if (lwp) allocate(tottmaIO(jpiglo,jpjglo,jpk))
       mem_all = mem_all +  rea_len*(12*jpiglo*jpjglo*jpk
     .          + 5*jpiglo*jpjglo)

       allocate(trb(jpi,jpj,jpk,jptra))
       mem_all = mem_all +  rea_len*jpi*jpj*jpk*jptra 


CC----------------------------------------------------------------------
CC
CC COMMON /cot3ad/ non-centered advection scheme (smolarkiewicz)
CC -------------------------------------------------------------
CC      rsc         : tuning coefficient for anti-diffusion (NAMELIST)
CC      rtrn        : value for truncation (NAMELIST)
CC
CC

CC----------------------------------------------------------------------
CC
CC COMMON /cit3ad/ non-centered advection scheme (smolarkiewicz)
CC -------------------------------------------------------------
CC      crosster    : logical if true computes crossterms (NAMELIST) 
CC      lhdf        : logical if true CALL trchdf (NAMELIST) 
CC      ncor        : number of corrective phases (NAMELIST)
CC      ndttrc      : frequency of step on passive tracers (NAMELIST)
CC

CC
CC
CC----------------------------------------------------------------------
CC
CC COMMON/citrst/ : passive tracers restart (input and output)
CC -----------------------------------------------------------
CC      nutwrs    : output FILE for passive tracers restart
CC      lrsttr    : boolean term for restart i/o for passive tracers 
CC                  (NAMELIST)
CC      nutrst    : logical unit for restart FILE for passive tracers
CC      nrsttr    : control of the time step ( 0 or 1 ) for pass. tr.
CC                  (NAMELIST)
CC
CC
CC----------------------------------------------------------------------
CC
CC COMMON/cotiso/ : isopycnal sheme for passive tracers 
CC -----------------------------------------------------------
CC      ahtrb0    : background diffusivity coefficient (m2/s)
CC                  for passive tracer
CC      trcrat    : ratio between passive and active tracer coeff
CC                  for diffusion
CC      ahtrc0    : horizontal eddy diffusivity for passive tracers (m2/s)
CC    aeivtr0   : eddy induced velocity coefficient (m2/s)
CC
CC

#    if defined key_trc_diaadd
CC----------------------------------------------------------------------
CC
CC       : additional 2D/3D outputs
CC ------------------------------------------------------------------
CC
CC      ctrc3d    : 3d output field name (NAMELIST)
CC      ctrc3l    : 3d output field long name (NAMELIST)
CC      ctrc3u    : 3d output field unit (NAMELIST)
CC      ctrc2d    : 2d output field name (NAMELIST)
CC      ctrc2l    : 2d output field long name (NAMELIST)
CC      ctrc2u    : 2d output field unit (NAMELIST)
CC      trc3d     : additional 3d outputs
CC      trc2d     : additional 2d outputs
CC
       allocate(ctrc3d(jpdia3d),allocate(ctrc2d(jpdia2d)
       allocate(ctrc3u(jpdia3d),allocate(ctrc2u(jpdia2d)
       allocate(ctrc3l(jpdia3d),allocate(ctrc2l(jpdia2d)
       allocate(trc3d(jpi,jpj,jpk,jpdia3d), allocate(trc2d(jpi,jpj,jpdia2d)

       mem_all = mem_all +  rea_len*(3*jpdia3d+3*jpdia2d +
     .           + jpi*jpj*jpk*jpdia3d+jpi*jpj*jpdia2d )
C     write(*,*) 'Mem_allocated:', mem_all

CC
CC    netcdf files and index common
CC
CC      nwriteadd: frequency of additional arrays outputs (NAMELIST)
CC      nitd     : id for additional array output FILE
CC      ndepitd  : id for depth mesh
CC      nhoritd  : id for horizontal mesh
CC
#    endif

#    if defined key_trc_diatrd
CC----------------------------------------------------------------------
CC
CC COMMON/cottrd/ : non conservative trends (biological, ...)
CC ------------------------------------------------------------------
CC
CC      luttrd    : large trends diagnostic to WRITE or not (NAMELIST)
CC
      allocate(luttrd(jptra))
       mem_all = mem_all +  rea_len*jptra
C     write(*,*) 'Mem_allocated:', mem_all
CC
CC    dynamical trends
CC
CC    trtrd()   : trends of the tracer equations
CC           1 : X advection
CC           2 : Y advection
CC           3 : Z advection
CC           4 : X diffusion
CC           5 : Y diffusion
CC           6 : Z diffusion
CC           7 : X gent velocity
CC           8 : Y gent velocity
CC           9 : Z gent velocity

      allocate(trtrd(jpi,jpj,jpk,jptra,jpdiatrc))
      mem_all = mem_all +  rea_len*jpi*jpj*jpk*jptra*jpdiatrc
C     write(*,*) 'Mem_allocated:', mem_all
CC
CC    netcdf files and index common
CC
CC      nwritetrd: frequency for dynamical trends output (NAMELIST)
CC                 one per tracer

        mem_all = mem_all +  rea_len*3*jptra
C     write(*,*) 'Mem_allocated:', mem_all
#    endif 
CC
CC----------------------------------------------------------------------
CC
CC     passive tracers DATA READ and at given time_step
CC -----------------------------------------------------------------
CC      nclimr    : switch for passive tracers initialization
CC      ncontr    : variable for time interpolation
CC      numtr1    : LOGICAL UNIT for passive tracers DATA
CC      numtr2    : LOGICAL UNIT for passive tracers DATA created IF
CC                  interpolation is needed (ninttr=1)
CC      nlectr    : switch for reading once
CC      ninttr    : switch for interpolation on model grid
CC      nmldmptr  : : = 0/1/2 type of damping in the mixed layer
CC      trcdat()  : passive tracers DATA array for two value
CC                  needed for time interpolation
CC
CC      trdta()   : passive tracers DATA at given time-step
CC      restotr() : array of restoring coeff. for passive tracers
CC
CC
CC
#    if defined key_trc_dmp 

      allocate(trdta(jpi,jpj,jpk,jptra))
      allocate(trcdat(jpi,jpj,jpk,2,jptra))
!      allocate(restoglo(jpiglo,jpjglo,jpk))
      allocate(idxt(jpi,jpj,jpk))
      allocate(idxt2glo(jpi,jpj,jpk,4))
!      allocate(idxtglo(jpiglo,jpjglo,jpk))
      mem_all = mem_all +  rea_len*(4*jpi*jpj*jpk*jptra
     &                  + 2*jpi*jpj*jpk + 0*jpiglo*jpjglo*jpk)
C     write(*,*) 'Mem_allocated:', mem_all
CC
#    else
CC      no passive tracers DATA at given time step
#    endif
CC
# endif

#ifdef key_trc_bfm
CC
CC      vsed      : sedimentation speed (NAMELIST)
CC      remine()  : remineralisation trend
CC
CC      tmaxr     : maximum coefficient for passive tracer damping (NAMELIST)
CC      tminr     : minimum coefficient for passive tracer damping (NAMELIST)
CC      remdmp()  : damping coefficient of passive tracers (depth dependant)
CC

      allocate(remine(jpi,jpj,jpk))
      allocate(remdmp(jpk,jptra))
      mem_all = mem_all +  rea_len*(jpi*jpj*jpk+jpk*jptra)
C     write(*,*) 'Mem_allocated:', mem_all

CC

CC
CC----------------------------------------------------------------------
CC
CC COMMON/cotopt/ : optical parameters
CC -----------------------------------
CC
CC      xze       : euphotic layer depth
CC      xpar      : par (photosynthetic available radiation)
CC      xkr0      : water coefficient absorption in red (NAMELIST)
CC      xkg0      : water coefficient absorption in green (NAMELIST)
CC      xkrp      : pigment coefficient absorption in red (NAMELIST)
CC      xkgp      : pigment coefficient absorption in green (NAMELIST)
CC      xlr       : exposant for pigment absorption in red (NAMELIST)
CC      xlg       : exposant for pigment absorption in green (NAMELIST)
CC      rpig      : chla/chla+phea ratio (NAMELIST)
CC

      allocate(xze(jpi,jpj))
      allocate(xpar(jpi,jpj,jpk))
      mem_all = mem_all +  rea_len*(jpi*jpj*jpk+jpi*jpj)
C     write(*,*) 'Mem_allocated:', mem_all
CC----------------------------------------------------------------------
CC

#endif

CC    photoperiod
        allocate(DAY_LENGTH(jpi,jpj))
  
        END subroutine alloc_tot

        END MODULE 
