      subroutine radmod_ex(V_POSITION, bottom,zd,Edtop,Estop,rmud,a,bt,bb,Edz,Esz,Euz,Eu0,PARz)
      USE OPT_mem
      USE Tridiagonal
      IMPLICIT NONE
 
!  Model of irradiance in the water column.  Accounts for three 
!  irradiance streams:
 
!  Edz = direct downwelling irradiance
!  Esz = diffuse downwelling irradiance
!  Euz = diffuse upwelling irradiance
 
!  Propagation is done in energy units, tests are done in quanta,
!  final is quanta for phytoplankton growth.
 
      integer, INTENT(IN) :: bottom
      double precision, INTENT(IN)    :: zd(jpk)
      double precision, INTENT(IN)    :: rmud
      double precision, INTENT(IN)    :: a(jpk,nlt),bt(jpk,nlt),bb(jpk,nlt)
      CHARACTER(*), INTENT(IN)        :: V_POSITION
      double precision, INTENT(IN)    :: Edtop(nlt),Estop(nlt)
      double precision, INTENT(OUT)   :: Edz(jpk,nlt),Esz(jpk,nlt),Euz(jpk,nlt),Eu0(nlt)     
      double precision, INTENT(OUT)   :: PARz(jpk,nchl+1)

! local variables
      double precision, parameter        :: rd=1.5d0   !these are taken from Ackleson, et al. 1994 (JGR)
      double precision, parameter        :: ru=3.0d0
      integer                            :: k,ii,b,p,nl
      double precision                   :: rmus, rmuu
      double precision                   :: Edzm(jpk,nlt),Edza(jpk,nlt)
      double precision                   :: vD(2*bottom-1,nlt),vL(2*bottom-1,nlt),vU(2*bottom-1,nlt)
      double precision                   :: WW(2*bottom-1,nlt),WW1(2*bottom-1,nlt),sol(2*bottom-1,nlt)
      double precision                   :: sol_p(2*bottom-1,nlt), sol_m(2*bottom-1,nlt)
      double precision                   :: aux1(nlt), aux2(nlt)

!  Constants
      rmus = 1.0d0/0.83d0            !avg cosine diffuse down
      rmuu = 1.0d0/0.4d0             !avg cosine diffuse up
      b =bottom 
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

!     
! Solution of the inhomogenous system (inhox and inhoy)
      inhoD(1:b,:)  = 1.0d0/((cd(1:b,:)-Cs(1:b,:))*(cd(1:b,:)+Cu(1:b,:))+Bs(1:b,:)*Bu(1:b,:))
      inhox(1:b,:)  = inhoD(1:b,:)*(- ( cd(1:b,:)+Cu(1:b,:) )*Fd(1:b,:) -Bu(1:b,:)*Bd(1:b,:) )
      inhoy(1:b,:)  = inhoD(1:b,:)*(- Fd(1:b,:)*Bs(1:b,:) + ( cd(1:b,:)-Cs(1:b,:) )*Bd(1:b,:) )
!     write(*,*) 'inhoD',inhoD(1,1)
!     write(*,*) 'inhox',inhox(1,1)
!     write(*,*) 'inhoy',inhoy(1,1)

! eigenvalues of the homogeneous system
      D(1:b,:)      = 0.5d0*( Cs(1:b,:)+Cu(1:b,:) &
                + sqrt((Cs(1:b,:)+Cu(1:b,:))*(Cs(1:b,:)+Cu(1:b,:))-4.0d0*Bs(1:b,:)*Bu(1:b,:) ) )
      a_m(1:b,:)    = D(1:b,:)  - Cs(1:b,:)
      a_p(1:b,:)    = Cs(1:b,:) - ( Bs(1:b,:) *Bu(1:b,:) )/D(1:b,:)
!     write(*,*) 'D',D(1,1)
!     write(*,*) 'a_m',a_m(1,1)
!     write(*,*) 'a_p',a_p(1,1)

   
! solution of the homogeneous system
      r_m(1:b,:)    = Bu(1:b,:)/D(1:b,:) 
      r_p(1:b,:)    = Bs(1:b,:)/D(1:b,:) 
!     write(*,*) 'r_m',r_m(1,1)
!     write(*,*) 'r_p',r_p(1,1)

      do k=1,bottom
         e_m(k,:)    = exp(-a_m(k,:)*zd(k))  
         e_p(k,:)    = exp(-a_p(k,:)*zd(k))  
      enddo

!     write(*,*) 'e_m',e_m(1,1)
!     write(*,*) 'e_p',e_p(1,1)

      Edaux(:)=Edtop(:)
      do k=1,bottom
         if (k > 1) then
            Edaux(:) = Edz(k-1,:)
         endif
         Edz(k,:)   = Edaux(:)*exp(-cd(k,:)*zd(k))
         Edzm(k,:)  = Edaux(:)*exp(-0.5D0*cd(k,:)*zd(k))   ! mid cell
      enddo

! Imposing boundary conditions

      zeta0(:)      = 1.0d0
      eta0(:)       = r_m(1,:)*e_m(1,:) 

      do k=1,bottom-1

         alpha(k,:) = e_p(k,:)*(1.0d0-r_p(k,:)*r_m(k+1,:) )
         beta(k,:)  = r_m(k,:) - r_m(k+1,:)
         gamm(k,:)  = - (1.0d0 -r_p(k+1,:)*r_m(k+1,:) )
         delta(k,:) = ( inhox(k+1,:) - inhox(k,:) - (inhoy(k+1,:) - inhoy(k,:) ) & 
                    * r_m(k+1,:) ) * Edz(k,:) !! Check if Edz(k) is correct

         epsRT(k,:)   = 1.0d0-r_m(k,:)*r_p(k,:) 
         zeta(k,:)  = -(r_p(k+1,:)-r_p(k,:))
         eta(k,:)   = -e_m(k+1,:)*(1.0d0-r_p(k,:)*r_m(k+1,:) )
         theta(k,:) = ( inhoy(k+1,:) - inhoy(k,:) - (inhox(k+1,:) - inhox(k,:) ) & 
                    * r_p(k,:) ) * Edz(k,:) !! Check if Edz(k) is correct
          
      enddo

! contruction of the vectors of the tri-diagonal matrix

! diagonal vector vD

      vD(1,:)= zeta0(:)

      do k =1,bottom-1 
         ii = 2*k  
         vD(ii,:)   = beta(k,:)
         vD(ii+1,:) = zeta(k,:)
      enddo

! lower vector 

      vL(1,:) = 0.0d0
      do k =1,bottom-1 
         ii = 2*k   
         vL(ii,:)   = alpha(k,:)
         vL(ii+1,:) = epsRT(k,:)
      enddo
  
! Upper vector 

      vU(1,:)=eta0(:)

      do k =1,bottom-2 
         ii = 2*k
         vU(ii,:)   = gamm(k,:)
         vU(ii+1,:) = eta(k,:)
      enddo
 
      vU(2*bottom-2,:)  = gamm(bottom-1,:)
      vU(2*bottom-1,:)  = 0.0d0
 
!  constant vector

      WW(1,:)= Estop(:) - inhox(1,:) * Edtop(:)

      do k =1,bottom-1
          ii = 2*k
          WW(ii,:)   = delta(k,:)
          WW(ii+1,:) = theta(k,:)
      enddo

! matrix inversion

! controls

!     do ii=1,nlt

!         if ((Edtop(ii) + Estop(ii)) < 0.0001) then 
!           write(*,*) 'Edtop(ii) + Estop(ii), nl->', ii
!           vL(:,ii) = 0.
!           vD(:,ii) = 0.0001d0
!           vU(:,ii) = 0.
!           write(*,*) '---------------------'
!        endif

!        if ( sum(abs(vL(:,ii))) < 0.0001) then 
!           write(*,*) 'vL, nl->', ii
!           vL(:,ii) = 0.
!           vD(:,ii) = 1.
!           vU(:,ii) = 0.
!           write(*,*) '---------------------'
!        endif

!        if ( sum(abs(vD(2:2*bottom-1,ii))) < 0.0001d0) then 
!           write(*,*) 'vD, nl->', ii
!           vL(:,ii) = 0.
!           vD(:,ii) = 1.0d0
!           vU(:,ii) = 0.
!           write(*,*) '---------------------'
!        endif

!        if (sum(abs(vU(:,ii))) < 0.0001) then 
!           write(*,*) 'vU, nl->', ii
!           vL(:,ii) = 0.
!           vD(:,ii) = 1.
!           vU(:,ii) = 0.
!           write(*,*) '---------------------'
!        endif

!     enddo

      sol(:,:) = 0.0d0

      call SLAE3diag(2*bottom-1, nlt , vL, vD, vU, WW, sol)

!     call matrix3diag_mult(2*bottom-1, nlt, vL, vD, vU, sol, WW1, WW, err_RT)

!     write(*,*) err
! putting the solution on the c_p and c_m vectors

      do k =1,bottom-1
          ii       = 2*k -1
          sol_p(k,:) = sol(ii,:)
          sol_m(k,:) = sol(ii+1,:)
      enddo
          
      sol_p(bottom,:) = sol(2*bottom-1,:)
      sol_m(bottom,:) = 0.0d0

! output solutions on bottom of layers
!     write(*,*)  'sol_p', sol_p(1,:)
!     write(*,*)  'a_p', a_p(1,:)
!     write(*,*)  'sol_m', sol_m(1,:)
!     write(*,*)  'r_m', r_m(1,:)
!     write(*,*)  'inhox', inhox(1,:)
      SELECT CASE (TRIM(V_POSITION))

         CASE ("BOTTOM")

             do k =1,bottom

                aux1(:)  = exp(-a_p(k,:)*zd(k))

                Esz(k,:) = sol_p(k,:)*aux1(:) + sol_m(k,:)*r_m(k,:) + inhox(k,:)*Edz(k,:)
                Euz(k,:) = sol_p(k,:)*r_p(k,:)*aux1(:) + sol_m(k,:) + inhoy(k,:)*Edz(k,:)

             enddo

         CASE ("MID")

             do k =1,bottom

                aux1(:)  = exp(-a_p(k,:)*zd(k)*0.5D0)
                aux2(:)  = exp(-a_m(k,:)*zd(k)*0.5D0)
                Edz(k,:) = Edzm(k,:) 
                Esz(k,:) = sol_p(k,:) * aux1(:)          + sol_m(k,:)*r_m(k,:)*aux2(:) + inhox(k,:)*Edzm(k,:)
                Euz(k,:) = sol_p(k,:)*r_p(k,:) * aux1(:) + sol_m(k,:)*aux2(:)          + inhoy(k,:)*Edzm(k,:)

             enddo

         CASE ("AVERAGE")

             do k =1,bottom

                if (k .EQ. 1) then
                    Edza(k,:) = (Edtop(:) -  Edz(k,:)) /( zd(k) * cd(k,:) + 0.000001d0 ) 
                else
                    Edza(k,:) = (Edz(k-1,:) -  Edz(k,:)) /( zd(k) * cd(k,:) + 0.000001d0 )
                endif

                aux1(:)  = ( 1.0D0 - exp(-a_p(k,:)*zd(k)) )/( zd(k) * a_p(k,:) + 0.000001d0 )
                aux2(:)  = ( 1.0D0 - exp( -a_m(k,:) *zd(k) ) ) /( zd(k) * a_m(k,:) + 0.000001d0 )

                Esz(k,:) = sol_p(k,:)* aux1(:)  + sol_m(k,:)*r_m(k,:) * aux2(:) + inhox(k,:)*Edza(k,:)
                Euz(k,:) = sol_p(k,:)*r_p(k,:) * aux1(:) + sol_m(k,:) * aux2(:) + inhoy(k,:)*Edza(k,:)

             enddo

             Edz(:,:) = Edza(:,:)

         CASE DEFAULT

             STOP

       END SELECT


! surface upward diffuse radiance at depth = 0 
      Eu0(:) = sol_p(1,:)*r_p(1,:) + sol_m(1,:)*exp(-a_m(1,:)*zd(1)) + inhoy(1,:) * Edtop(:)
!     write(*,*) 'Eu', Eu

! compute PAR for each PFT using scalar irradiance
      PARz(:,:) = 0.0D0

      do p =1,nchl
         do nl =5,19   ! 400 to 700 nm 
            do k =1,bottom

               PARz(k,p) = PARz(k,p) + WtoQ(nl) * ac(p,nl) * ( Edz(k,nl) * rmud + Esz(k,nl) * rmus + Euz(k,nl) * rmuu )   

            enddo
         enddo
      enddo
        
      do nl =5,19   ! 400 to 700 nm 
         do k =1,bottom

            PARz(k,nchl+1) = PARz(k,nchl+1) + WtoQ(nl) * ( Edz(k,nl) * rmud + Esz(k,nl) * rmus + Euz(k,nl) * rmuu )   

         enddo
      enddo
      
      return
      end
