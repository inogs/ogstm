       MODULE modulo16


       IMPLICIT NONE

       public
        integer jpj16,jpi16,jpk16
        integer mem_all, rea_len, int_len


!C      glamt16          : longitude of t-point (degre)
!C      glamu16          : longitude of u-point (degre)
!C      glamv16          : longitude of v-point (degre)
!C      glamf16          : longitude of f-point (degre)
!C      gphit16          : latitude  of t-point (degre)
!C      gphiu16          : latitude  of u-point (degre)
!C      gphiv16          : latitude  of v-point (degre)
!C      gphif16          : latitude  of f-point (degre)
!C      e1t16,e2t16        : horizontal scale factors at t-point (m)
!C      e1u16,e2u16        : horizontal scale factors at u-point (m)
!C      e1v16,e2v16        : horizontal scale factors at v-point (m)
!C      e1f16,e2f16        : horizontal scale factors at f-point (m)
!C	ff16             : coriolis factor (2.*omega*sin(yphi) ) (s-1)
!C-CC   dso16,dse16,dno16,dne16 :Relative distances
!
      REAL(4), allocatable ::  glamt16(:,:), glamu16(:,:), glamv16(:,:),glamf16(:,:)
      REAL(4), allocatable :: gphit16(:,:), gphiu16(:,:), gphiv16(:,:),gphif16(:,:)
      REAL(4), allocatable :: e1t16(:,:), e1u16(:,:), e1v16(:,:), e1f16(:,:)
      REAL(4), allocatable :: e2t16(:,:), e2u16(:,:), e2v16(:,:), e2f16(:,:)
      REAL(4), allocatable :: ff16(:,:)
      REAL(4), allocatable :: dso16(:,:),dse16(:,:),dno16(:,:),dne16(:,:)
!C                  z-coordinate (default option)
!C                  ------------------------------
!C      gdept16, gdepw16() : depth of t- and w-points (m)
!C      e3t16, e3w16()     : vertical scale factors at t- and w-points (m)
!C
      REAL(4), allocatable :: gdept16(:), gdepw16(:)
      REAL(4), allocatable :: e3t16(:,:,:), e3u16(:,:,:),e3v16(:,:,:),e3w16(:,:,:)
!C
!C----------------------------------------------------------------------
!C Common/comask/  : masks, bathymetry
!C -----------------------------------
!C      numbat           : logical un16it for bathymetry file
!C      mbathy         : number of ocean level (=0, 1, ... , jpk16-1)
!C      tmask16, umask16() : land/ocean mask at t-, u-, v- and f-points
!C      vmask16, fmask16()
!
      INTEGER numbat
      INTEGER, allocatable :: mbathy(:,:)
!
!
      REAL(4), allocatable :: tmask16(:,:,:), fmask16(:,:,:)
      REAL(4), allocatable :: umask16(:,:,:), vmask16(:,:,:)
!
!C
!C
!C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!C
!C II. DYNAMICS AND TRACERS
!C ========================
!C
!C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!C
!C----------------------------------------------------------------------
!C Common/comnow/  : present fields (now)
!C -------------------------------------
!C	un16, vn16(), wn16() : horizontal and vertical velocity (m s-1)
!C	hdivn16          : horizontal divergence
!C	tn16,   sn16()     : pot. temperature (celsius), salinity (psu)
!C      rdn16            : in situ density anomalie rdn16=(rho-rau0)/rau0
!C                         (no un16its)
!C	rhopn16          : potential volumic mass (kg m-3)
!C      bn2n16           : brun16t-vaisala frequency (s-2)
!C
      REAL(4), allocatable :: un16(:,:,:), vn16(:,:,:), wn16(:,:,:)
      REAL(4), allocatable :: hdivn16(:,:,:)
      REAL(4), allocatable :: tn16(:,:,:), sn16(:,:,:)
      REAL(4), allocatable :: rdn16(:,:,:), rhopn16(:,:,:), bn2n16(:,:,:)
!
      REAL(4), allocatable :: ahtu16(:), ahtv16(:), ahtw16(:), ahtt16(:)

!C
!C----------------------------------------------------------------------
!C Common/comzdf/ : vertical diff16usion
!C -----------------------------------
!C	avt160             : vertical viscosity and diff16. coef. (namelist)
!C      ntrbbl		 : bottom boun16dary layer (namelist)
!C      atrbbl		 : bottom boun16dary layer diff16usivity (namelist)
!C	avt16	         : vertical diff16usivity coeff16. at w-point
!C	avt16b           : backgroun16d profile of avm and avt16
!
!
      INTEGER ntrbbl
      REAL(4) avt016,atrbbl
      REAL(4), allocatable :: avt16(:,:,:)
      REAL(4), allocatable ::  avtb16(:)
      REAL(4), allocatable ::  bblx16(:,:),bbly16(:,:)
!
!C
!C IV. SURFACE FORCING AND DATA
!C ============================
!C
!C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!C
!C----------------------------------------------------------------------
!C Common/comtau/ : surface wind stress at givem time_step
!C -------------------------------------------------------
!C	taux16, tauy16()   : wind stress components in (i,j) referential
!C	taux16g, tauy16g() : zonal and meridian wind stress component used
!C	                   in output (geographical referential)
!
      REAL(4), allocatable :: taux16(:,:), tauy16(:,:), vatm16(:,:), frezze16(:,:)
      REAL(4), allocatable :: aux16(:,:), auy16(:,:)
      REAL(4), allocatable :: aux16_01(:,:), auy16_01(:,:)
      REAL(4), allocatable :: aux16_02(:,:), auy16_02(:,:)
      REAL(4), allocatable :: aux16_03(:,:), auy16_03(:,:)
      REAL(4), allocatable :: aux16_04(:,:), auy16_04(:,:)
      REAL(4), allocatable :: aux16_05(:,:), auy16_05(:,:)
      REAL(4), allocatable :: aux16_06(:,:), auy16_06(:,:)
!
!C
!C----------------------------------------------------------------------
!C Common/comflx/ : surface fluxes
!C -------------------------------
!C      qt16             : total surface heat flux (w m-2)
!C      q16              : surface heat flux (w m-2)
!C      emp16            : evaporation minus precipitation (mm day-1)
!C      runoff16         : annual run16 off16 (mm/day)
!
      REAL(4), allocatable :: qt16(:,:), q16(:,:), emp16(:,:)
      REAL(4), allocatable :: runoff16(:,:)
!
!C
!C----------------------------------------------------------------------
!C
!C COMMON/comqsr/ : penetrative solar radiation
!C --------------------------------------------
!C      qsr     : solar radiation (w m-2)
!C
      REAL(4), allocatable :: qsr16(:,:)
      REAL(4), allocatable :: wsp16(:,:)
      REAL(4), allocatable :: ssh16(:,:)
      double precision, allocatable :: ssh16_R8(:,:),ssh16x_R8(:,:),ssh16y_R8(:,:)
      REAL(4), allocatable :: flp16(:,:)
      CONTAINS

       subroutine alloc_tot_16()
      int_len  = 2
      rea_len = 4
      mem_all  = 0
      mem_all = mem_all + int_len*jpk16 
!     write(*,*) 'Mem_allocated:', mem_all
!
        allocate(glamt16(jpj16,jpi16))
        allocate(glamu16(jpj16,jpi16)) 
        allocate(glamv16(jpj16,jpi16))
        allocate(glamf16(jpj16,jpi16))
        allocate(gphit16(jpj16,jpi16)) 
        allocate(gphiu16(jpj16,jpi16)) 
        allocate(gphiv16(jpj16,jpi16))
        allocate(gphif16(jpj16,jpi16))
        allocate(e1t16(jpj16,jpi16)) 
        allocate(e1u16(jpj16,jpi16)) 
        allocate(e1v16(jpj16,jpi16)) 
        allocate(e1f16(jpj16,jpi16))
        allocate(e2t16(jpj16,jpi16)) 
        allocate(e2u16(jpj16,jpi16)) 
        allocate(e2v16(jpj16,jpi16)) 
        allocate(e2f16(jpj16,jpi16))
        allocate(ff16(jpj16,jpi16))
        allocate(dso16(jpj16,jpi16))
        allocate(dse16(jpj16,jpi16))
        allocate(dno16(jpj16,jpi16))
        allocate(dne16(jpj16,jpi16))
        mem_all = mem_all + rea_len*(+21*jpi16*jpj16) 
!     write(*,*) 'Mem_allocated:', mem_all

!
!C
!C----------------------------------------------------------------------
!C Common/comcoz/  : vertical coordinate and scale factors
!C -------------------------------------------------------
!C
!C                  z-coordinate (default option)
!C                  ------------------------------
!C      gdept16(), gdepw16() : depth of t- and w-points (m)
!C      e3t16(), e3w16()     : vertical scale factors at t- and w-points (m)
!C
      allocate(gdept16(jpk16)) 
      allocate(gdepw16(jpk16))
      allocate(e3t16(jpj16,jpi16,jpk16))
      allocate(e3u16(jpj16,jpi16,jpk16))
      allocate(e3v16(jpj16,jpi16,jpk16))
      allocate(e3w16(jpj16,jpi16,jpk16))
       mem_all = mem_all + rea_len*(4*jpk16) 
!     write(*,*) 'Mem_allocated:', mem_all
!
!C
!C----------------------------------------------------------------------
!C Common/comask/  : masks, bathymetry
!C -----------------------------------
!C      numbat           : logical un16it for bathymetry file
!C      mbathy()         : number of ocean level (=0, 1, ... , jpk16-1)
!C      tmask16(), umask16() : land/ocean mask at t-, u-, v- and f-points
!C      vmask16(), fmask16()
!
      allocate(tmask16(jpj16,jpi16,jpk16))
      allocate(fmask16(jpj16,jpi16,jpk16))
      allocate(umask16(jpj16,jpi16,jpk16))
      allocate(vmask16(jpj16,jpi16,jpk16))
      mem_all = mem_all + int_len*jpi16*jpj16 + rea_len*4*jpi16*jpj16*jpk16  
!     write(*,*) 'Mem_allocated:', mem_all
!
!C
!C
!C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!C
!C----------------------------------------------------------------------
!C Common/comnow/  : present fields (now)
!C -------------------------------------
!C	un16(), vn16(), wn16() : horizontal and vertical velocity (m s-1)
!C	hdivn16()          : horizontal divergence
!C	tn16(),   sn16()     : pot. temperature (celsius), salinity (psu)
!C      rdn16()            : in situ density anomalie rdn16=(rho-rau0)/rau0
!C                         (no un16its)
!C	rhopn16()          : potential volumic mass (kg m-3)
!C      bn2n16()           : brun16t-vaisala frequency (s-2)
!C
       allocate(un16(jpj16,jpi16,jpk16))
       allocate(vn16(jpj16,jpi16,jpk16))
       allocate(wn16(jpj16,jpi16,jpk16))
       allocate(hdivn16(jpj16,jpi16,jpk16))
       allocate(tn16(jpj16,jpi16,jpk16))
       allocate(sn16(jpj16,jpi16,jpk16))
       allocate(rdn16(jpj16,jpi16,jpk16))
       allocate(rhopn16(jpj16,jpi16,jpk16))
       allocate(bn2n16(jpj16,jpi16,jpk16))
       mem_all = mem_all +  rea_len*9*jpi16*jpj16*jpk16  
!     write(*,*) 'Mem_allocated:', mem_all
!
!C
!C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!C III. OCEAN PHYSICS
!C ==================
!C
!C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!C
!C----------------------------------------------------------------------
!C Common/comhdt/ : lateral diff16usivity (tracers)
!C ------------------------------------
!C      aht0             : lateral diff16usivity coeff16icient    (namelist)
!C      ahtb0            : backgroun16d diff16usivity coeff16icient (m2/s)
!C      ahtu16(), ahtv16()   : lateral diff16usivity coef. at u-, v-, w- t-pts
!C      ahtw16(), ahtt16()     (harmonic operator: no rotation, use of u-
!C                          and v-points; rotation, use of u-, v- w-pts)
!C                         (biharmonic operator: rotation or not, use of
!C                          t-point only)
!C                         (the arrays used are 3D, 2D, 1D or 0D depen-
!C                          ding on 'key_trahdfcoef.d' )
!
!CC  Paolo 22/4/2004 new conditions 1D case attention
      allocate(ahtu16(jpk16))
      allocate(ahtv16(jpk16))
      allocate(ahtw16(jpk16))
      allocate(ahtt16(jpk16))
      mem_all = mem_all +  rea_len*4*jpk16  
!     write(*,*) 'Mem_allocated:', mem_all

!C
!C----------------------------------------------------------------------
!C Common/comzdf/ : vertical diff16usion
!C -----------------------------------
!C	avt016             : vertical viscosity and diff16. coef. (namelist)
!C      ntrbbl		 : bottom boun16dary layer (namelist)
!C      atrbbl		 : bottom boun16dary layer diff16usivity (namelist)
!C	avt16()	         : vertical diff16usivity coeff16. at w-point
!C	avtb16()           : backgroun16d profile of avm and avt16
!
!
       allocate(avt16(jpj16,jpi16,jpk16))
       allocate(avtb16(jpk16))
       allocate(bblx16(jpj16,jpi16))
       allocate(bbly16(jpj16,jpi16))
       mem_all = mem_all +  rea_len*(jpi16*jpj16*jpk16+ jpk16 +2*jpi16*jpj16) 
!     write(*,*) 'Mem_allocated:', mem_all
!
!C
!C
!C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!C
!C IV. SURFACE FORCING AND DATA
!C ============================
!C
!C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!C
!C----------------------------------------------------------------------
!C Common/comtau/ : surface wind stress at givem time_step
!C -------------------------------------------------------
!C	taux16(), tauy16()   : wind stress components in (i,j) referential
!C	taux16g(), tauy16g() : zonal and meridian wind stress component used
!C	                   in output (geographical referential)
!
       allocate(taux16(jpj16,jpi16))
       allocate(tauy16(jpj16,jpi16))
       allocate(aux16(jpj16,jpi16))
       allocate(auy16(jpj16,jpi16))
       allocate(aux16_01(jpj16,jpi16))
       allocate(auy16_01(jpj16,jpi16))
       allocate(aux16_02(jpj16,jpi16))
       allocate(auy16_02(jpj16,jpi16))
       allocate(aux16_03(jpj16,jpi16))
       allocate(auy16_03(jpj16,jpi16))
       allocate(aux16_04(jpj16,jpi16))
       allocate(auy16_04(jpj16,jpi16))
       allocate(aux16_05(jpj16,jpi16))
       allocate(auy16_05(jpj16,jpi16))
       allocate(aux16_06(jpj16,jpi16))
       allocate(auy16_06(jpj16,jpi16))
       allocate(vatm16(jpj16,jpi16))
       allocate(frezze16(jpj16,jpi16))
       mem_all = mem_all +  rea_len*(5*jpi16*jpj16) 
!     write(*,*) 'Mem_allocated:', mem_all
!
!C
!C----------------------------------------------------------------------
!C Common/comflx/ : surface fluxes
!C -------------------------------
!C      qt16()             : total surface heat flux (w m-2)
!C      q16()              : surface heat flux (w m-2)
!C      emp16()            : evaporation minus precipitation (mm day-1)
!C      runoff16()         : annual run16 off16 (mm/day)
!
       allocate(qt16(jpj16,jpi16)) 
       allocate(q16(jpj16,jpi16)) 
       allocate(emp16(jpj16,jpi16))
       allocate(runoff16(jpj16,jpi16))
       mem_all = mem_all +  rea_len*(4*jpi16*jpj16) 
!     write(*,*) 'Mem_allocated:', mem_all
!
!C
!C----------------------------------------------------------------------
!C
!C COMMON/comqsr/ : penetrative solar radiation
!C --------------------------------------------
!C      qsr16()     : solar radiation (w m-2)
!C
       allocate(qsr16(jpj16,jpi16))
       allocate(wsp16(jpj16,jpi16))
       allocate(ssh16(jpj16,jpi16))
       allocate(ssh16_R8(jpj16,jpi16))
       allocate(ssh16x_R8(jpj16,jpi16))
       allocate(ssh16y_R8(jpj16,jpi16))
       allocate(flp16(jpj16,jpi16))
       mem_all = mem_all +  rea_len*(4*jpi16*jpj16)
!     write(*,*) 'Mem_allocated:', mem_all
!C

  
        END subroutine alloc_tot_16

        END MODULE 
