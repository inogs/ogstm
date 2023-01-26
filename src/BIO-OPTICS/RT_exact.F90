      subroutine radmod_ex(V_POSITION, bottom,zd,Edtop,Estop,rmud,a,bt,bb,Edz,Esz,Euz,Eu0,PARz,SWRz)
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
      double precision, INTENT(OUT)   :: PARz(jpk,nchl+1),SWRz(jpk)

! local variables
      double precision, parameter        :: rd=1.5d0   !these are taken from Ackleson, et al. 1994 (JGR)
      double precision, parameter        :: ru=3.0d0
      integer                            :: k,ii,b,p,nl
      logical                            :: is_Singular(jpk,nlt)
      CHARACTER(20)                      :: V_POS
      double precision                   :: rmus, rmuu
      double precision                   :: RT_depth(jpk)
      double precision                   :: Edzm(jpk,nlt),Edza(jpk,nlt)
      double precision                   :: vD(2*bottom-1,nlt),vL(2*bottom-1,nlt),vU(2*bottom-1,nlt)
      double precision                   :: WW(2*bottom-1,nlt),WW1(2*bottom-1,nlt),sol(2*bottom-1,nlt)
      double precision                   :: sol_p(2*bottom-1,nlt), sol_m(2*bottom-1,nlt)
      double precision                   :: aux1(nlt), aux2(nlt)

      Edzm=huge(Edzm(1,1))
      Edza=huge(Edza(1,1))
      vD=huge(vD(1,1))
      vL=huge(vL(1,1))
      vU=huge(vU(1,1))
      WW=huge(WW(1,1))
      WW1=huge(WW1(1,1))
      sol=huge(sol(1,1))
      sol_p=huge(sol_p(1,1))
      sol_m=huge(sol_m(1,1))
      aux1=huge(aux1(1))
      aux2=huge(aux2(1))

      Edz=huge(Edz(1,1))
      Esz=huge(Esz(1,1))
      Euz=huge(Euz(1,1))
      Eu0=huge(Eu0(1))

      V_POS = V_POSITION
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
      do nl=1,nlt
          do k =1,bottom
              inhoD(k,nl)  = ((cd(k,nl)-Cs(k,nl))*(cd(k,nl)+Cu(k,nl))+Bs(k,nl)*Bu(k,nl))
              if (abs(inhoD(k,nl)) .LT. 1.0D-4) then
                 is_Singular(k,nl) = .True.
                 V_POS="MID"
              else
                 is_Singular(k,nl) = .False. 
              endif
          enddo
      enddo
      do nl=1,nlt
          do k =1,bottom
              if (is_Singular(k,nl) ) then
                  inhox(k,nl)  = Fd(k,nl)
                  inhoy(k,nl)  = -Bd(k,nl)
              else
                  inhoD_r(k,nl)  = 1.0d0/inhoD(k,nl)
                  inhox(k,nl)  = inhoD_r(k,nl)*(- ( cd(k,nl)+Cu(k,nl) )*Fd(k,nl) -Bu(k,nl)*Bd(k,nl) )
                  inhoy(k,nl)  = inhoD_r(k,nl)*(- Fd(k,nl)*Bs(k,nl) + ( cd(k,nl)-Cs(k,nl) )*Bd(k,nl) )
              endif
          enddo
      enddo
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
         Edz(k,:)     = Edaux(:)*exp(-cd(k,:)*zd(k))
         Edzm(k,:)    = Edaux(:)*exp(-0.5D0*cd(k,:)*zd(k))   ! mid cell
         if (k .EQ.  1) then
            RT_depth(k)  = zd(k)
         else
            RT_depth(k)  = RT_depth(k-1)+zd(k)
         endif
      enddo

! Imposing boundary conditions

      zeta0(:)      = 1.0d0
      eta0(:)       = r_m(1,:)*e_m(1,:) 

      do k=1,bottom-1

         alpha(k,:) = e_p(k,:)*(1.0d0-r_p(k,:)*r_m(k+1,:) )
         beta(k,:)  = r_m(k,:) - r_m(k+1,:)
         gamm(k,:)  = - (1.0d0 -r_p(k+1,:)*r_m(k+1,:) )
!        delta(k,:) = ( inhox(k+1,:) - inhox(k,:) - (inhoy(k+1,:) - inhoy(k,:) ) & 
!                   * r_m(k+1,:) ) * Edz(k,:) !! Check if Edz(k) is correct

         epsRT(k,:)   = 1.0d0-r_m(k,:)*r_p(k,:) 
         zeta(k,:)  = -(r_p(k+1,:)-r_p(k,:))
         eta(k,:)   = -e_m(k+1,:)*(1.0d0-r_p(k,:)*r_m(k+1,:) )
!        theta(k,:) = ( inhoy(k+1,:) - inhoy(k,:) - (inhox(k+1,:) - inhox(k,:) ) & 
!                   * r_p(k,:) ) * Edz(k,:) !! Check if Edz(k) is correct
          
      enddo  

! Theta and Delta depends on the inhomogenous solution and their determination will be different according
! to presence of singularities

      do nl=1,nlt
          do k =1,bottom-1

!  k   Singular
!  -
!  k+1 Singular
              if (is_Singular(k,nl) .AND. is_Singular(k+1,nl)) then

                  delta(k,nl) = ( inhox(k+1,nl) - inhox(k,nl) - (inhoy(k+1,nl) - inhoy(k,nl) ) & 
                        * r_m(k+1,nl) ) * RT_depth(k) * Edz(k,nl) !! 
                  theta(k,nl) = ( inhoy(k+1,nl) - inhoy(k,nl) - (inhox(k+1,nl) - inhox(k,nl) ) & 
                        * r_p(k,nl) ) * RT_depth(k) * Edz(k,nl) !! 

              endif

!  k   Non Singular
!  -
!  k+1 Singular 

              if( (.NOT. is_Singular(k,nl)) .AND. (is_Singular(k+1,nl))) then

                  delta(k,nl) = ( RT_depth(k) * inhox(k+1,nl) - inhox(k,nl) - (RT_depth(k) *inhoy(k+1,nl) - inhoy(k,nl) ) & 
                        * r_m(k+1,nl) )  * Edz(k,nl) !! 
                  theta(k,nl) = ( RT_depth(k) * inhoy(k+1,nl) - inhoy(k,nl) - (RT_depth(k) *inhox(k+1,nl) - inhox(k,nl) ) & 
                        * r_p(k,nl) ) * Edz(k,nl) !! 
              endif

!  k   Singular
!  -
!  k+1 Non Singular 

              if( (is_Singular(k,nl)) .AND. ( .NOT. is_Singular(k+1,nl))) then

                  delta(k,nl) = ( inhox(k+1,nl) - RT_depth(k) *inhox(k,nl) - (inhoy(k+1,nl) - RT_depth(k) *inhoy(k,nl) ) & 
                        * r_m(k+1,nl) )  * Edz(k,nl) !! 
                  theta(k,nl) = ( inhoy(k+1,nl) - RT_depth(k) *inhoy(k,nl) - (inhox(k+1,nl) - RT_depth(k) *inhox(k,nl) ) & 
                        * r_p(k,nl) ) * Edz(k,nl) !! 

              endif

!  k   Non Singular
!  -
!  k+1 Non Singular 


              if (( .NOT. is_Singular(k,nl)) .AND. ( .NOT. is_Singular(k+1,nl))) then

                  delta(k,nl) = ( inhox(k+1,nl) - inhox(k,nl) - (inhoy(k+1,nl) - inhoy(k,nl) ) & 
                        * r_m(k+1,nl) ) * Edz(k,nl)
                  theta(k,nl) = ( inhoy(k+1,nl) - inhoy(k,nl) - (inhox(k+1,nl) - inhox(k,nl) ) & 
                        * r_p(k,nl) ) * Edz(k,nl) 

              endif

          enddo
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

      do nl=1,nlt
          if (is_Singular(1,nl) ) then
              WW(1,nl)= Estop(nl) 
          else
              WW(1,nl)= Estop(nl) - inhox(1,nl) * Edtop(nl)
          endif
      enddo

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

         do nl=1,nlt
             do k =1,bottom


                aux1(nl)  = exp(-a_p(k,nl)*zd(k)*0.5D0)
                aux2(nl)  = exp(-a_m(k,nl)*zd(k)*0.5D0)
                Edz(k,nl) = Edzm(k,nl) 

                if (is_Singular(k,nl) ) then
                    Esz(k,nl) = sol_p(k,nl)           * aux1(nl) + sol_m(k,nl)*r_m(k,nl)*aux2(nl) + inhox(k,nl)*Edzm(k,nl) * (RT_depth(k)-zd(k)*0.5D0)
                    Euz(k,nl) = sol_p(k,nl)*r_p(k,nl) * aux1(nl) + sol_m(k,nl)          *aux2(nl) + inhoy(k,nl)*Edzm(k,nl) * (RT_depth(k)-zd(k)*0.5D0)
                else
                    Esz(k,nl) = sol_p(k,nl)           * aux1(nl) + sol_m(k,nl)*r_m(k,nl)*aux2(nl) + inhox(k,nl)*Edzm(k,nl) 
                    Euz(k,nl) = sol_p(k,nl)*r_p(k,nl) * aux1(nl) + sol_m(k,nl)          *aux2(nl) + inhoy(k,nl)*Edzm(k,nl) 
                endif

             enddo
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
      do nl=1,nlt             
          if (is_Singular(1,nl) ) then
              Eu0(:) = sol_p(1,:)*r_p(1,:) + sol_m(1,:)*exp(-a_m(1,:)*zd(1)) 
          else
              Eu0(:) = sol_p(1,:)*r_p(1,:) + sol_m(1,:)*exp(-a_m(1,:)*zd(1)) + inhoy(1,:) * Edtop(:)
          endif
      enddo
!     write(*,*) 'Eu', Eu
                  
!     subroutine radmod_ex(V_POSITION, bottom,zd,Edtop,Estop,rmud,a,bt,bb,Edz,Esz,Euz,Eu0,PARz,SWRz)
      do nl=1,nlt             
      if (abs(Eu0(nl)) .GT. 1000.) then
          write(1100,*) 'size PARz', size(PARz) 
          write(1101,*) 'V_POSITION', V_POSITION 
          write(1102,*) 'bottom', bottom
          write(1103,*) 'zd', zd
          write(1104,*) 'Edtop', Edtop
          write(1105,*) 'Estop', Estop
          write(1106,*) 'rmud', rmud
          write(1107,*) 'a', a(:,nl)
          write(1108,*) 'bt', bt(:,nl)
          write(1109,*) 'bb', bb(:,nl)
          write(1110,*) 'Edz', Edz(:,nl)
          write(1111,*) 'Esz', Esz(:,nl)
          write(1112,*) 'Euz', Euz(:,nl)
          write(1113,*) 'Eu0', Eu0(:)
          write(1114,*) 'PARz', PARz
          write(1115,*) 'SWRz', SWRz
          
!     call SLAE3diag(2*bottom-1, nlt , vL, vD, vU, WW, sol)
          write(1116,*) 'vL', vL(:,nl)
          write(1117,*) 'vD', vL(:,nl)
          write(1118,*) 'vU', vL(:,nl)
          write(1119,*) 'WW', WW(:,nl)
          write(1120,*) 'sol', sol(:,nl)
          write(1121,*) 'cd', cd(:,nl)
          write(1122,*) 'Cs', Cs(:,nl)
          write(1123,*) 'Cu', Cu(:,nl)
          write(1124,*) 'Bs', Bs(:,nl)
          write(1125,*) 'Bu', Bu(:,nl)
          write(1126,*) 'Fd', Fd(:,nl)
          write(1127,*) 'Bd', Bd(:,nl)
          write(1128,*) 'inhoD', inhoD(:,nl)
          write(1129,*) 'inhox', inhox(:,nl)
          STOP
      endif
      enddo
! compute PAR for each PFT using scalar irradiance
      PARz(:,:) = 0.0D0
      SWRz(:) = 0.0D0

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
      
      do nl =1,nlt ! from 0 to 4 um
         do k =1,bottom

            SWRz(k) = SWRz(k) +( Edz(k,nl) * rmud + Esz(k,nl) * rmus + Euz(k,nl) * rmuu )   

         enddo
      enddo

      return
      end
