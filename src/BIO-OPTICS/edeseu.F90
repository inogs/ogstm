      subroutine edeseu(MODE,bottom,dzRT,Edtop,Estop,CHLz,CDOMz,NAPz,rmud,Edz,Esz,Euz,Eutop)
      USE myalloc
      USE mpi
      USE OPT_mem
      IMPLICIT NONE
!  Model of irradiance in the water column.  Accounts for three 
!  irradiance streams:
 
!  Edz = direct downwelling irradiance
!  Esz = diffuse downwelling irradiance
!  Euz = diffuse upwelling irradiance
 
!  Propagation is done in energy units, tests are done in quanta,
!  final is quanta for phytoplankton growth.
 
      integer :: n, nl, jk 
      integer :: MODE
      integer :: bottom
      double precision :: rmud
      double precision :: dzRT(jpk)
      double precision :: Edz(jpk,nlt),Esz(jpk,nlt)
      double precision :: Euz(jpk,nlt)
      double precision :: Edtop(nlt),Estop(nlt)
      double precision :: Eutop(nlt)
      double precision :: CHLz(jpk,nchl),CDOMz(jpk),NAPz(jpk)
      double precision :: Etop, Ebot
      double precision :: Plte
      double precision :: actot(jpk,nlt),bctot(jpk,nlt),bbctot(jpk,nlt) 
      double precision :: acdom(nlt),anap(nlt)
      double precision :: a(jpk,nlt), bt(jpk,nlt), bb(jpk,nlt) 
      double precision :: bbc(5)
      data bbc /0.002, 0.00071, 0.0032, 0.00071, 0.0029/
      double precision bbw
      data bbw /0.5/       !backscattering to forward scattering ratio
      double precision rmus
 
!  Constants and initialize
      rmus = 1.0/0.83            !avg cosine diffuse down
 
!  Compute irradiance with depth
       Ebot = 0.0
       Eu   = 0.0

       acdom = 0.0 
       anap  = 0.0

       do nl = 1,nlt

!       Edtop(nl) = Ed(nl)
!       Estop(nl) = Es(nl)
!       Ebot = Ebot + (Ed(nl)+Es(nl))

       enddo

!      Ebot = 0.0

       do jk = 1,bottom

!       Etop = Ebot

        do nl = 1,nlt

         Edz(jk,nl) = 0.0
         Esz(jk,nl) = 0.0
         Euz(jk,nl) = 0.0
         actot(jk,nl) = 0.0
         bctot(jk,nl) = 0.0
         bbctot(jk,nl) = 0.0

         do n = 1,nchl
!              Plte = max(P(k,n),0.0)
               actot(jk,nl)  = actot(jk,nl)  + CHLz(jk,n)*ac(n,nl)
               bctot(jk,nl)  = bctot(jk,nl)  + CHLz(jk,n)*bc(n,nl)
               bbctot(jk,nl) = bbctot(jk,nl) + CHLz(jk,n)*bbc(n)*bc(n,nl)
         enddo

         a(jk,nl)  = aw(nl) + CDOMz(jk) * acdom(nl) + NAPz(jk) * anap(nl) + actot(jk,nl)
         bt(jk,nl) = bw(nl) + bctot(jk,nl)
         bb(jk,nl) = bbw*bw(nl) + bbctot(jk,nl)
         bb(jk,nl) = max(bb(jk,nl),0.0002)

        enddo

       enddo
       
!      write(*,*) 'a', a(bottom,1)
!      write(*,*) 'bt', bt(bottom,1)
!      write(*,*) 'bb', bb(bottom,1)
!      write(*,*) 'rmud', rmud


       if (MODE == 1) then ! condition on level < 200 mt
           call radmod(bottom,dzRT(:),Edtop,Estop,rmud,a,bt,bb,Edz,Esz,Euz,Eutop)
       endif 
! condition on level > 200 mt
       if ( MODE ==0) then
           call radmod_ex(bottom,dzRT(:),Edtop,Estop,rmud,a,bt,bb,Edz,Esz,Euz,Eutop)
       endif


      return
      end
