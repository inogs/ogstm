      subroutine radmod(bottom,zd,Edtop,Estop,rmud,a,bt,bb,Edz,Esz,Euz,Eu0)
      USE OPT_mem
      IMPLICIT NONE
 
!  Model of irradiance in the water column.  Accounts for three 
!  irradiance streams:
 
!  Edz = direct downwelling irradiance
!  Esz = diffuse downwelling irradiance
!  Euz = diffuse upwelling irradiance
 
!  Uses Ackelson's (1994, JGR) mod's to the Aas (1987, AO) 
!  two-stream model.
 
!  Propagation is done in energy units, tests are done in quanta,
!  final is quanta for phytoplankton growth.
 
!  Commented out terms produce a max error of 
!  0.8% in Esz for a > 0.004 and bb > 0.0001 and
!  3.9% in Euz for a > 0.004 and bb > 0.00063

      integer, INTENT(IN)    :: bottom
      real(8), INTENT(IN)    :: zd(jpk)
      real(8), INTENT(IN)    :: rmud
      real(8), INTENT(IN)    :: a(jpk,nlt),bt(jpk,nlt),bb(jpk,nlt)
      real(8), INTENT(IN)    :: Edtop(nlt),Estop(nlt)
      real(8), INTENT(OUT)   :: Edz(jpk,nlt),Esz(jpk,nlt),Euz(jpk,nlt),Eu0(nlt)     

! local variables
      real(8), parameter     :: rd=1.5   !these are taken from Ackleson, et al. 1994 (JGR)
      real(8), parameter     :: ru=3.0
      integer                :: jk,b
      real(8)                :: rmus, rmuu

!  Constants
      rmus = 1.0/0.83            !avg cosine diffuse down
      rmuu = 1.0/0.4             !avg cosine diffuse up
      b = bottom 
!  Downwelling irradiance: Edz, Esz
!  Compute irradiance components at depth
      cd(1:b,:)     = (a(1:b,:)+bt(1:b,:))*rmud
      au(1:b,:)     = a(1:b,:)*rmuu
      Bu(1:b,:)     = ru*bb(1:b,:)*rmuu
      Cu(1:b,:)     = au(1:b,:)+Bu(1:b,:)
      as(1:b,:)     = a(1:b,:)*rmus
      Bs(1:b,:)     = rd*bb(1:b,:)*rmus
      Cs(1:b,:)     = as(1:b,:)+Bs(1:b,:)
      Bd(1:b,:)     = bb(1:b,:)*rmud
      Fd(1:b,:)     = (bt(1:b,:)-bb(1:b,:))*rmud
      bquad(1:b,:)  = Cs(1:b,:) - Cu(1:b,:)
      cquad(1:b,:)  = Bs(1:b,:)*Bu(1:b,:) - Cs(1:b,:)*Cu(1:b,:)
      sqarg(1:b,:)  = bquad(1:b,:)*bquad(1:b,:) - 4.0*cquad(1:b,:)
!     sqarg(1:b,:) = max(0.0000000001,sqarg(1:b,:))
      a1(1:b,:)     = 0.5*(-bquad(1:b,:) + sqrt(sqarg(1:b,:)))
      a2(1:b,:)     = 0.5*(-bquad(1:b,:) - sqrt(sqarg(1:b,:)))
      S(1:b,:)      = -(Bu(1:b,:)*Bd(1:b,:) + Cu(1:b,:)*Fd(1:b,:))

      Edaux = Edtop
      Esaux = Estop
      do jk=1,bottom
         if (jk > 1) then
            Edaux(:) = Edz(jk-1,:)
            Esaux(:) = Esz(jk-1,:)
         endif
         Edz(jk,:)   = Edaux(:)*exp(-cd(jk,:)*zd(jk))
         SEdz(jk,:)  = S(jk,:)*Edz(jk,:)
         a2ma1(jk,:) = a2(jk,:) - a1(jk,:)
         rM(jk,:)    = SEdz(jk,:)/(a1(jk,:)*a2ma1(jk,:))
         rN(jk,:)    = SEdz(jk,:)/(a2(jk,:)*a2ma1(jk,:))
         c2(jk,:)    = Esaux(:) - rM(jk,:) + rN(jk,:) 
         Ta2z(jk,:)  = exp(a2(jk,:)*zd(jk))
         Esz(jk,:)   = c2(jk,:)*Ta2z(jk,:) + rM(jk,:) - rN(jk,:)
         Esz(jk,:)   = max(Esz(jk,:),0.0)
         Eutmp(jk,:) = ((a2(jk,:)+Cs(jk,:))*c2(jk,:))*Ta2z(jk,:)  &
                  + Cs(jk,:)*rM(jk,:) - Cs(jk,:)*rN(jk,:) - Fd(jk,:)*Edz(jk,:)
         Euz(jk,:)   = Eutmp(jk,:)/Bu(jk,:)
         Euz(jk,:)   = max(Euz(jk,:),0.0)
      enddo

      Eu0(:) = -1.0 

      return
      end
