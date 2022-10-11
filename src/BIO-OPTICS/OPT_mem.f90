       MODULE OPT_mem 

       USE modul_param 
       USE myalloc

#ifdef Mem_Monitor
       USE check_mem
       USE iso_c_binding
#endif

       IMPLICIT NONE

       public


      INTEGER, allocatable :: itabe(:),imaske(:,:) 
      double precision, allocatable :: zpar(:,:),xEPS_ogstm(:,:)
      double precision, allocatable :: zpar0m(:),zpar100(:) 
!      double precision, allocatable :: kef(:,:)
!      double precision, allocatable :: kextIO(:,:,:)

! ######     Oasim inputs
      double precision, allocatable :: spIO(:,:,:)
      double precision, allocatable :: mslIO(:,:,:)
      double precision, allocatable :: t2mIO(:,:,:)
      double precision, allocatable :: d2mIO(:,:,:)
      double precision, allocatable :: tccIO(:,:,:)
      double precision, allocatable :: w10IO(:,:,:)
      double precision, allocatable :: tclwIO(:,:,:)
      double precision, allocatable :: tco3IO(:,:,:)
      double precision, allocatable :: cdremIO(:,:,:)
      double precision, allocatable :: cldtcmIO(:,:,:)
      double precision, allocatable :: tauaIO(:,:,:,:)
      double precision, allocatable :: asympIO(:,:,:,:)
      double precision, allocatable :: ssalbIO(:,:,:,:)

      double precision, allocatable :: sp(:,:)
      double precision, allocatable :: msl(:,:)
      double precision, allocatable :: t2m(:,:)
      double precision, allocatable :: d2m(:,:)
      double precision, allocatable :: tcc(:,:)
      double precision, allocatable :: w10(:,:)
      double precision, allocatable :: tclw(:,:)
      double precision, allocatable :: tco3(:,:)
      double precision, allocatable :: cdrem(:,:)
      double precision, allocatable :: cldtcm(:,:)
      double precision, allocatable :: taua(:,:,:)
      double precision, allocatable :: asymp(:,:,:)
      double precision, allocatable :: ssalb(:,:,:)
! ############
      real, allocatable :: zkef_f (:,:)
      CHARACTER(LEN=7), allocatable,dimension(:) :: Ednm, Esnm, Eunm
      INTEGER, allocatable, dimension(:) :: EdWR,EsWR,EuWR
      INTEGER, allocatable, dimension(:) :: Ed3D,Es3D,Eu3D
      integer Ed_dumped_vars, Es_dumped_vars, Eu_dumped_vars





 
      integer, parameter            :: nchl=4 
      integer, parameter            :: nlt=33                     
      integer                       :: lam(33)
      double precision              :: WtoQ(33)
! Radiative transfer model parameter OASIM Native coordinates
      double precision              :: Ed_0m_COARSE(33,12,18,48), Es_0m_COARSE(33,12,18,48) ! lon, lat, day period, wave length
      double precision              :: OASIM_lon(18,48), OASIM_lat(18,48)  
! Radiative transfer model parameter OGSTM coordinates    
      double precision,allocatable  :: Ed_0m(:,:,:), Es_0m(:,:,:) ! wav, lat, lon
      
      INTEGER                       :: day_RTcheck
! in-water model
      INTEGER                       :: it_check
      double precision              :: aw(33),bw(33)
      double precision              :: ac(4,33),ac_ps(4,33),bc(4,33),bbc(4,33)
      double precision              :: acdom(3,33)
      double precision              :: apoc(33), bpoc(33), bbpoc(33) 
! Variables related to computation of solar zenith angle
      double precision              :: rad 
      integer                       :: imon, nutime 
      double precision              :: dpsi, eps
      double precision,allocatable  :: up(:,:,:),no(:,:,:),ea(:,:,:)

! Avergae cosine computation
      double precision, parameter  :: refrac_idx = 1.341D0
! Constant of aberration
      double precision, parameter  :: xk = 0.0056932D0
!  compute irradiance every hours
      logical                      ::  ifst = .TRUE.
      integer, parameter           ::  nstps1 = 5
      double precision             ::  delmin, hrsec, hrsrt, hrend, delh, delx   
      integer                      ::  nstps

      double precision,allocatable  :: Edaux(:), Esaux(:)
      double precision,allocatable  :: cd(:,:),Cs(:,:),Bu(:,:),Cu(:,:),Bs(:,:),Fd(:,:),Bd(:,:) 
      double precision,allocatable  :: au(:,:),as(:,:),bquad(:,:),cquad(:,:),sqarg(:,:)
      double precision,allocatable  :: inhoD(:,:),inhox(:,:),inhoy(:,:)
      double precision,allocatable  :: D(:,:),a_m(:,:),a_p(:,:)
      double precision,allocatable  :: r_m(:,:),r_p(:,:)
      double precision,allocatable  :: e_m(:,:),e_p(:,:)
      double precision,allocatable  :: zeta0(:),eta0(:)
      double precision,allocatable  :: alpha(:,:),beta(:,:),gamm(:,:),delta(:,:)
      double precision,allocatable  :: epsRT(:,:),zeta(:,:),eta(:,:),theta(:,:)
!     double precision,allocatable  :: vD(:,:),vL(:,:),vU(:,:),WW(:,:),WW1(:,:)
!     double precision,allocatable  :: sol(:,:),sol_p(:,:),sol_m(:,:)
      double precision,allocatable  :: err_RT(:)
! Additional variables for approximate model
      double precision, allocatable        :: a1(:,:),a2(:,:),S(:,:)    
      double precision, allocatable        :: SEdz(:,:),a2ma1(:,:),rM(:,:),rN(:,:),c2(:,:) 
      double precision, allocatable        :: Ta2z(:,:), Eutmp(:,:)

! Outputs of radiative transfer model 
      double precision,allocatable  :: Ed(:,:,:,:), Es(:,:,:,:), Eu(:,:,:,:) ! depth, lat, lon, wave length
      double precision,allocatable  :: PAR(:,:,:,:) ! depth, lat, lon, phyto
      double precision,allocatable  :: RMU(:,:)     ! lat, lon
      double precision,allocatable  :: Ed_IO(:,:,:,:),Es_IO(:,:,:,:),Eu_IO(:,:,:,:)
      logical, allocatable, dimension(:,:,:) :: opt_mask
      INTEGER, allocatable, dimension(:) :: Ed_table, Es_table, Eu_table




!----------------------------------------------------------------------
      CONTAINS

      subroutine myalloc_OPT()
      INTEGER  :: err
      double precision  :: aux_mem

! local variables 
      INTEGER           :: nl
      double precision  :: h, c, hc, oavo, hcoavo, rlamm

#ifdef Mem_Monitor
       aux_mem = get_mem(err)
#endif
       allocate(itabe(jpi))         
      
       itabe   = huge(itabe(1))
       allocate(imaske(jpk,jpi))   
       imaske  = huge(imaske(1,1))
!!!$omp parallel default (none) shared(jpk,jpi)
       allocate(zpar(jpk,jpi))     
       zpar    = huge(zpar(1,1))
       allocate(xEPS_ogstm(jpk,jpi))     
       xEPS_ogstm    = huge(xEPS_ogstm(1,1))
       allocate(zpar0m(jpi))        
       zpar0m  = huge(zpar0m(1))
       allocate(zpar100(jpi))       
       zpar100 = huge(zpar100(1))
!!!$omp end parallel

!       allocate(kef(jpj,jpi))
!       kef     = huge(kef(1,1))
!       allocate(kextIO(jpj,jpi,2))
!       kextIO  = huge(kextIO(1,1,1))

       allocate(sp (jpj,jpi))     ;   sp =huge(sp(1,1))
       allocate(msl(jpj,jpi))     ;   msl=huge(sp(1,1))
       allocate(t2m(jpj,jpi))     ;   t2m=huge(sp(1,1))
       allocate(d2m(jpj,jpi))     ;   d2m=huge(sp(1,1))
       allocate(tcc(jpj,jpi))     ;   tcc=huge(sp(1,1))
       allocate(w10(jpj,jpi))     ;   w10=huge(sp(1,1))
       allocate(tclw(jpj,jpi))    ;  tclw=huge(sp(1,1))
       allocate(tco3(jpj,jpi))    ;  tco3=huge(sp(1,1))
       allocate(cdrem(jpj,jpi))   ; cdrem=huge(sp(1,1))
       allocate(cldtcm(jpj,jpi))  ;cldtcm=huge(sp(1,1))
       allocate(taua(33,jpj,jpi)) ;  taua=huge(sp(1,1))
       allocate(asymp(33,jpj,jpi)); asymp=huge(sp(1,1))
       allocate(ssalb(33,jpj,jpi)); ssalb=huge(sp(1,1))

       allocate(spIO (jpj,jpi,2))     ; spIO   =huge(sp(1,1))
       allocate(mslIO(jpj,jpi,2))     ; mslIO  =huge(sp(1,1))
       allocate(t2mIO(jpj,jpi,2))     ; t2mIO  =huge(sp(1,1))
       allocate(d2mIO(jpj,jpi,2))     ; d2mIO  =huge(sp(1,1))
       allocate(tccIO(jpj,jpi,2))     ; tccIO  =huge(sp(1,1))
       allocate(w10IO(jpj,jpi,2))     ; w10IO = huge(sp(1,1))
       allocate(tclwIO(jpj,jpi,2))    ; tclwIO =huge(sp(1,1))
       allocate(tco3IO(jpj,jpi,2))    ; tco3IO =huge(sp(1,1))
       allocate(cdremIO(jpj,jpi,2))   ; cdremIO=huge(sp(1,1))
       allocate(cldtcmIO(jpj,jpi,2))  ;cldtcmIO=huge(sp(1,1))
       allocate(tauaIO(33,jpj,jpi,2)) ; tauaIO =huge(sp(1,1))
       allocate(asympIO(33,jpj,jpi,2)); asympIO=huge(sp(1,1))
       allocate(ssalbIO(33,jpj,jpi,2)); ssalbIO=huge(sp(1,1))




! radiative transfer model

       rad    = 180.0D0/dacos(-1.0D0) ! initialization of radians
       imon   = 1
       nutime = -99999

       allocate(Ed_0m(nlt,jpj,jpi))
       Ed_0m  =huge(Ed_0m(1,1,1))
       allocate(Es_0m(nlt,jpj,jpi))
       Es_0m  =huge(Es_0m(1,1,1))

       allocate(up(jpj,jpi,3))
       up     =huge(up(1,1,1))

       allocate(no(jpj,jpi,3))
       no     =huge(no(1,1,1))

       allocate(ea(jpj,jpi,3))
       ea     =huge(ea(1,1,1))

       call lidata()

       allocate(Edaux(nlt))
       allocate(Esaux(nlt))
       allocate(cd(jpk,nlt),Cs(jpk,nlt),Bu(jpk,nlt),Cu(jpk,nlt),Bs(jpk,nlt),Fd(jpk,nlt),Bd(jpk,nlt))
       allocate(au(jpk,nlt),as(jpk,nlt),bquad(jpk,nlt),cquad(jpk,nlt),sqarg(jpk,nlt))
       allocate(inhoD(jpk,nlt),inhox(jpk,nlt),inhoy(jpk,nlt))
       allocate(D(jpk,nlt),a_m(jpk,nlt),a_p(jpk,nlt))
       allocate(r_m(jpk,nlt),r_p(jpk,nlt))
       allocate(e_m(jpk,nlt),e_p(jpk,nlt))
       allocate(zeta0(nlt),eta0(nlt))
       allocate(alpha(jpk-1,nlt),beta(jpk-1,nlt),gamm(jpk-1,nlt),delta(jpk-1,nlt))
       allocate(epsRT(jpk-1,nlt),zeta(jpk-1,nlt),eta(jpk-1,nlt),theta(jpk-1,nlt))
!      allocate(vD(2*jpk-1,nlt),vL(2*jpk-1,nlt),vU(2*jpk-1,nlt))
!      allocate(WW(2*jpk-1,nlt), WW1(2*jpk-1,nlt))
!      allocate(sol(2*jpk-1,nlt),sol_p(jpk,nlt),sol_m(jpk,nlt))
       allocate(err_RT(nlt))
! Additional variables for approximate model
       allocate(a1(jpk,nlt),a2(jpk,nlt),S(jpk,nlt))
       allocate(SEdz(jpk,nlt),a2ma1(jpk,nlt),rM(jpk,nlt),rN(jpk,nlt),c2(jpk,nlt))
       allocate(Ta2z(jpk,nlt), Eutmp(jpk,nlt))

! Allocate output variables
       allocate(Ed(jpk,jpj,jpi,nlt),Es(jpk,jpj,jpi,nlt),Eu(jpk,jpj,jpi,nlt))
       allocate(PAR(jpk,jpj,jpi,nchl+1)) ! last index total par
       allocate(RMU(jpj,jpi))
       allocate(Ed_IO(jpk,jpj,jpi,nlt),Es_IO(jpk,jpj,jpi,nlt),Eu_IO(jpk,jpj,jpi,nlt))
       allocate(opt_mask(jpk,jpj,jpi)) ; opt_mask=0



      Ed(:,:,:,:)   = 1.e-4
      Es(:,:,:,:)   = 1.e-4
      Eu(:,:,:,:)   = 1.e-8
      PAR(:,:,:,:)  = 0.0d0
      RMU(:,:)  = 0.0d0
      Ed_IO(:,:,:,:)   = 1.e-4
      Es_IO(:,:,:,:)   = 1.e-4
      Eu_IO(:,:,:,:)   = 1.e-8


      h = 6.6256E-34   !Plancks constant J sec
      c = 2.998E8      !speed of light m/sec
      hc = 1.0D0/(h*c)
      oavo = 1.0D0/6.023E23   ! 1/Avogadros number
      hcoavo = hc*oavo
      do nl = 1,nlt
       rlamm = real(lam(nl),8)*1.0E-9      !lambda in m
       WtoQ(nl) = rlamm*hcoavo*1000000.0D0 !Watts to micro mol quanta conversion
      enddo



#ifdef Mem_Monitor
              mem_all=get_mem(err) - aux_mem
#endif
        
      END subroutine myalloc_OPT

      SUBROUTINE opt_lec

      IMPLICIT NONE
      namelist /ED_nam/      Ednm, EdWR, Ed3D
      namelist /ES_nam/      Esnm, EsWR, Es3D
      namelist /EU_nam/      Eunm, EuWR, Eu3D
      !local
      integer jn

       allocate(Ednm(nlt))
       allocate(Esnm(nlt))
       allocate(Eunm(nlt))

       allocate(EdWR(nlt)); EdWR=huge(EdWR(1))
       allocate(EsWR(nlt)); EsWR=huge(EdWR(1))
       allocate(EuWR(nlt)); EuWR=huge(EdWR(1))

       allocate(Ed3D(nlt)); Ed3D=huge(Ed3D(1))
       allocate(Es3D(nlt)); Es3D=huge(Ed3D(1))
       allocate(Eu3D(nlt)); Eu3D=huge(Ed3D(1))


      OPEN(unit=numnat, file='namelist.optics', status= 'OLD')

      REWIND(numnat)
      READ(numnat,ED_nam)

      REWIND(numnat)
      READ(numnat,ES_nam)

      REWIND(numnat)
      READ(numnat,EU_nam)

      CLOSE(numnat)

      Ed_dumped_vars=0
      Es_dumped_vars=0
      Eu_dumped_vars=0

      do jn=1,nlt
        if (EdWR(jn)==1)  Ed_dumped_vars=Ed_dumped_vars+1
        if (EsWR(jn)==1)  Es_dumped_vars=Es_dumped_vars+1
        if (EuWR(jn)==1)  Eu_dumped_vars=Eu_dumped_vars+1
      enddo



      allocate(Ed_table(Ed_dumped_vars))
      allocate(Es_table(Es_dumped_vars))
      allocate(Eu_table(Eu_dumped_vars))

      Ed_dumped_vars=0
      do jn=1,nlt
         if (EdWR(jn).eq.1) then
            Ed_dumped_vars=Ed_dumped_vars+1
            Ed_table(Ed_dumped_vars) = jn
            if (lwp) WRITE(numout,*) Ednm(jn),' will be dumped'
         endif
      enddo


      Es_dumped_vars=0
      do jn=1,nlt
         if (EsWR(jn).eq.1) then
            Es_dumped_vars=Es_dumped_vars+1
            Es_table(Es_dumped_vars) = jn
            if (lwp) WRITE(numout,*) Esnm(jn),' will be dumped'
         endif
      enddo

      Eu_dumped_vars=0
      do jn=1,nlt
         if (EuWR(jn).eq.1) then
            Eu_dumped_vars=Eu_dumped_vars+1
            Eu_table(Eu_dumped_vars) = jn
            if (lwp) WRITE(numout,*) Eunm(jn),' will be dumped'
         endif
      enddo



      END SUBROUTINE opt_lec

      SUBROUTINE init_OPT
      IMPLICIT NONE
      integer jk,jj,ji, tmask_levels

          do ji=1,jpi
          do jj=1,jpj

              if (bfmmask(1,jj,ji).eq.0 ) cycle
              tmask_levels = mbathy(jj,ji)
              opt_mask(1:tmask_levels+1,jj,ji) = 1

          enddo
          enddo
      END SUBROUTINE init_OPT
        
        

      subroutine clean_memory_opt

          deallocate(itabe)
          deallocate(imaske)
          deallocate(zpar)
          deallocate(xEPS_ogstm)
          deallocate(zpar0m)
          deallocate(zpar100)
!          deallocate(kef)
!          deallocate(kextIO)

      end subroutine clean_memory_opt



      END MODULE 
