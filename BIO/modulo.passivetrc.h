C $Id: modulo.passivetrc.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
#if defined key_passivetrc
CCC F79 Namelist does not accept allocatable
CCC          CHARACTER(LEN=8), allocatable :: ctrcnm(:),ctrcun(:)
CCC          CHARACTER(LEN=80), allocatable ::  ctrcnl(:)
          CHARACTER(LEN=8) :: ctrcnm(jptra),ctrcun(jptra)
          CHARACTER(LEN=6) :: mycln
          CHARACTER(LEN=80) ::  ctrcnl(jptra)

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
CCC F79 Namelist does not accept allocatable
CCC      LOGICAL, allocatable ::  lutini(:)
CCC      INTEGER, allocatable ::  nutini(:)
      LOGICAL ::  lutini(jptra)
      INTEGER ::  nutini(jptra)
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
      INTEGER flagSMS_Dyn 
      REAL(8) trai
      REAL(8), allocatable ::  trn(:,:,:,:)
      REAL(8), allocatable ::  tra(:,:,:,:)
      REAL(8), allocatable ::  traIO(:,:,:,:) 
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
      REAL(8), allocatable :: tottrn(:,:,:)
      REAL(8), allocatable :: tottrb(:,:,:)
CCC 26 9 2004 Paolo matrix for i/o writing(trcdit.F)
      REAL(8), allocatable ::  tottrnIO(:,:,:)
      REAL(8), allocatable ::  tottrnFN(:,:,:)
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
      REAL(8), allocatable ::  rsttrn(:,:,:,:)
      REAL(8), allocatable ::  rsttrb(:,:,:,:)

CC
      REAL(8), allocatable ::  trb(:,:,:,:)
CC
CC----------------------------------------------------------------------
CC
CC COMMON /cot3ad/ non-centered advection scheme (smolarkiewicz)
CC -------------------------------------------------------------
CC      rsc         : tuning coefficient for anti-diffusion (NAMELIST)
CC      rtrn        : value for truncation (NAMELIST)
CC
      REAL(8) rsc,rtrn
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
      INTEGER ncor,ndttrc
      LOGICAL lhdf,crosster
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
      LOGICAL lrsttr
      INTEGER nutwrs,nutrst,nrsttr
CC
COMMON/cotiso/ : isopycnal sheme for passive tracers 
CC -----------------------------------------------------------
CC      ahtrb0    : background diffusivity coefficient (m2/s)
CC                  for passive tracer
CC      trcrat    : ratio between passive and active tracer coeff
CC                  for diffusion
CC      ahtrc0    : horizontal eddy diffusivity for passive tracers (m2/s)
CC	aeivtr0   : eddy induced velocity coefficient (m2/s)
CC
      REAL(8) ahtrb0,trcrat,ahtrc0,aeivtr0
CC
CC
CC----------------------------------------------------------------------
CC
CC COMMON/citcdf/ : information for outputs
CC ------------------------------------------------------------------
CC
CC      nwritetrc: time step frequency for concentration outputs (NAMELIST)
CC      ndext50   : INTEGER arrays for ocean 3D INDEX 
CC      ndext51   : INTEGER arrays for ocean surface INDEX 
CC
CC    netcdf files and index common

      INTEGER nwritetrc

      INTEGER, allocatable :: ndext50(:) 
     $       ,ndext51(:)
     $       ,ndexttot50(:) 
     $       ,ndexttot51(:) 
            
      REAL(8) djulian

#    if defined key_trc_diaadd
CC----------------------------------------------------------------------
CC
CC COMMON/cot23d/ : additional 2D/3D outputs
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
CC	trtrd()   : trends of the tracer equations
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
      INTEGER nwritetrd
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
CC
      INTEGER nclimr,ncontr,numtr1,numtr2,nlectr,ninttr,nmldmptr
CC
CC
#    if defined key_trc_dmp 

      REAL(8), allocatable ::  trdta(:,:,:,:),trcdat(:,:,:,:,:)
CC
#    else
CC      no passive tracers DATA at given time step
#    endif
CC
CC
CC
#    if defined key_trc_bfm
#    include "modulo.passivetrc.bfm.h"
#    endif
CC
CC
#    if defined key_trc_npzdb
#    include "modulo.passivetrc.npzdb.h"
#    endif
CC
CC
#    if defined key_trc_generic
#    include "modulo.passivetrc.generic.h"
#    endif
#else
CC
CC no passive tracer COMMON specification
CC
#endif
