      subroutine radmod(V_POSITION,bottom,zd,Edtop,Estop,rmud,a,bt,bb,Edz,Esz,Euz,Eu0,PARz,SWRz)
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

      integer, INTENT(IN)             :: bottom
      double precision, INTENT(IN)    :: zd(jpk)
      double precision, INTENT(IN)    :: rmud
      double precision, INTENT(IN)    :: a(jpk,nlt),bt(jpk,nlt),bb(jpk,nlt)
      double precision, INTENT(IN)    :: Edtop(nlt),Estop(nlt)
      character(*), INTENT(IN)        ::  V_POSITION
      double precision, INTENT(OUT)   :: Edz(jpk,nlt),Esz(jpk,nlt),Euz(jpk,nlt),Eu0(nlt)     
      double precision, INTENT(OUT)   :: PARz(jpk,nchl+1)
      double precision, INTENT(OUT)   :: SWRz(jpk)

! local variables
      double precision, parameter     :: rd=1.5d0   !these are taken from Ackleson, et al. 1994 (JGR)
      double precision, parameter     :: ru=3.0d0
      double precision                :: Edzm(jpk,nlt),Eszm(jpk,nlt) ! mid cell irradiance
      double precision                :: Edza(jpk,nlt),Esza(jpk,nlt) ! average irradiance
      integer                         :: jk,b,p,k,nl
      double precision                :: rmus, rmuu

!  Constants
      rmus = 1.0d0/0.83d0            !avg cosine diffuse down
      rmuu = 1.0d0/0.4d0             !avg cosine diffuse up
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
      sqarg(1:b,:)  = bquad(1:b,:)*bquad(1:b,:) - 4.0D0*cquad(1:b,:)
!     sqarg(1:b,:) = max(0.0000000001,sqarg(1:b,:))
      a1(1:b,:)     = 0.5d0*(-bquad(1:b,:) + sqrt(sqarg(1:b,:)))
      a2(1:b,:)     = 0.5d0*(-bquad(1:b,:) - sqrt(sqarg(1:b,:)))
      S(1:b,:)      = -(Bu(1:b,:)*Bd(1:b,:) + Cu(1:b,:)*Fd(1:b,:))

      Edaux = Edtop
      Esaux = Estop

      SELECT CASE (TRIM(V_POSITION))

      CASE ("BOTTOM")

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
            Esz(jk,:)   = max(Esz(jk,:),0.0d0)
            Eutmp(jk,:) = ((a2(jk,:)+Cs(jk,:))*c2(jk,:))*Ta2z(jk,:)  &
                     + Cs(jk,:)*rM(jk,:) - Cs(jk,:)*rN(jk,:) - Fd(jk,:)*Edz(jk,:)
            Euz(jk,:)   = Eutmp(jk,:)/Bu(jk,:)
            Euz(jk,:)   = max(Euz(jk,:),0.0d0)

         enddo

      CASE ("MID")

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
            Esz(jk,:)   = max(Esz(jk,:),0.0d0)
            Eutmp(jk,:) = ((a2(jk,:)+Cs(jk,:))*c2(jk,:))*Ta2z(jk,:)  &
                     + Cs(jk,:)*rM(jk,:) - Cs(jk,:)*rN(jk,:) - Fd(jk,:)*Edz(jk,:)
            Euz(jk,:)   = Eutmp(jk,:)/Bu(jk,:)
            Euz(jk,:)   = max(Euz(jk,:),0.0d0)
! Mid cell irradiances
            Edzm(jk,:)  = Edaux(:)*exp(-0.5D0*cd(jk,:)*zd(jk))
            SEdz(jk,:)  = S(jk,:)*Edzm(jk,:)
            a2ma1(jk,:) = a2(jk,:) - a1(jk,:)
            rM(jk,:)    = SEdz(jk,:)/(a1(jk,:)*a2ma1(jk,:))
            rN(jk,:)    = SEdz(jk,:)/(a2(jk,:)*a2ma1(jk,:))
            c2(jk,:)    = Esaux(:) - rM(jk,:) + rN(jk,:)
            Ta2z(jk,:)  = exp(0.5D0*a2(jk,:)*zd(jk))
            Eszm(jk,:)   = c2(jk,:)*Ta2z(jk,:) + rM(jk,:) - rN(jk,:)
            Eszm(jk,:)   = max(Eszm(jk,:),0.0d0)
            Eutmp(jk,:) = ((a2(jk,:)+Cs(jk,:))*c2(jk,:))*Ta2z(jk,:)  &
                     + Cs(jk,:)*rM(jk,:) - Cs(jk,:)*rN(jk,:) - Fd(jk,:)*Edzm(jk,:)
            Euz(jk,:)   = Eutmp(jk,:)/Bu(jk,:)
            Euz(jk,:)   = max(Euz(jk,:),0.0d0)

         enddo

         do jk=1,bottom
            Edz(jk,:) = Edzm(jk,:)
            Esz(jk,:) = Eszm(jk,:)
         enddo

      CASE ("AVERAGE")

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
            Esz(jk,:)   = max(Esz(jk,:),0.0d0)
            Eutmp(jk,:) = ((a2(jk,:)+Cs(jk,:))*c2(jk,:))*Ta2z(jk,:)  &
                     + Cs(jk,:)*rM(jk,:) - Cs(jk,:)*rN(jk,:) - Fd(jk,:)*Edz(jk,:)
            Euz(jk,:)   = Eutmp(jk,:)/Bu(jk,:)
            Euz(jk,:)   = max(Euz(jk,:),0.0d0)
! cell averaged irradiances
            Edza(jk,:)  = Edaux(:)*(1.0D0 - exp(-cd(jk,:)*zd(jk)) ) / (zd(jk) * cd(jk,:)+0.000001d0)
            SEdz(jk,:)  = S(jk,:)*Edza(jk,:)
!           a2ma1(jk,:) = a2(jk,:) - a1(jk,:)
            rM(jk,:)    = SEdz(jk,:)/(a1(jk,:)*a2ma1(jk,:)) ! take average inhomogenous solution
            rN(jk,:)    = SEdz(jk,:)/(a2(jk,:)*a2ma1(jk,:)) ! take average inhomogenous solution
!           c2(jk,:)    = Esaux(:) - rM(jk,:) + rN(jk,:)
            Ta2z(jk,:)  = (exp(a2(jk,:)*zd(jk)) -1.0D0)/(zd(jk) * a2(jk,:)+0.000001d0)

            Esza(jk,:)   = Esaux(:)*Ta2z(jk,:) - S(jk,:)/ (a1(jk,:)*a2ma1(jk,:)) *  &
                           (exp((-cd(jk,:) +a2(jk,:))*zd(jk)) - 1.0D0)/ (zd(jk) * (a2(jk,:) -cd(jk,:))+0.000001d0) &
                           +  S(jk,:)/ (a2(jk,:)*a2ma1(jk,:)) *  &
                           (exp((-cd(jk,:) +a2(jk,:))*zd(jk)) - 1.0D0)/ (zd(jk) * (a2(jk,:) -cd(jk,:))+0.000001d0) & 
                           + rM(jk,:) - rN(jk,:)

            Esza(jk,:)   = max(Esza(jk,:),0.0d0)
!           Eutmp(jk,:)  = (a2(jk,:)+Cs(jk,:))*(Esaux(:)*Ta2z(jk,:) - S(jk,:)/ (a1(jk,:)*a2ma1(jk,:)) *  &
!                          (exp((-cd(jk,:) +a2(jk,:))*zd(jk)) - 1.0D0)/ (zd(jk) * (a2(jk,:) -cd(jk,:))+0.000001d0) &
!                          +  S(jk,:)/ (a2(jk,:)*a2ma1(jk,:)) *  &
!                          (exp((-cd(jk,:) +a2(jk,:))*zd(jk)) - 1.0D0)/ (zd(jk) * (a2(jk,:) -cd(jk,:))+0.000001d0)) & 
!                          + Cs(jk,:)*rM(jk,:) - Cs(jk,:)*rN(jk,:) - Fd(jk,:)*Edza(jk,:)
            Euz(jk,:)   = 0.0D0 ! not computed in this case seems approximation is not good vs analityc solution

         enddo


         do jk=1,bottom
            Edz(jk,:) = Edza(jk,:)
            Esz(jk,:) = Esza(jk,:)
         enddo

      CASE  DEFAULT 

         STOP
      END SELECT


      Eu0(:) = -1.0d0 


! compute PAR for each PFT using scalar irradiance
      PARz(:,:) = 0.0D0

      do p =1,nchl
         do nl =5,19   ! 400 to 700 nm 
            do jk =1,bottom
               PARz(jk,p) = PARz(jk,p) + WtoQ(nl) * ac(p,nl) * ( Edz(jk,nl) * rmud + Esz(jk,nl) * rmus + Euz(jk,nl) * rmuu )   
            enddo
         enddo
      enddo
        
      do nl =5,19  ! 400 to 700 nm  

         do jk =1,bottom

            PARz(jk,nchl+1) = PARz(jk,nchl+1) + WtoQ(nl) * ( Edz(jk,nl) * rmud + Esz(jk,nl) * rmus + Euz(jk,nl) * rmuu )   

         enddo
      enddo

      do nl =1,nlt  ! 0 to 4 um  
         do jk =1,bottom
            SWRz(jk) = SWRz(jk) +  Edz(jk,nl) * rmud + Esz(jk,nl) * rmus + Euz(jk,nl) * rmuu    
         enddo
      enddo

      return
      end
