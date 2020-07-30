      subroutine edeseu_IOP(MODE,V_POSITION,bottom,dzRT,Edtop,Estop,CHLz,CDOMz,NAPz,rmud,Edz,Esz,Euz,Eutop,PARz)
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
      integer, INTENT(IN)           :: MODE
      character(*), INTENT(IN)      :: V_POSITION
      integer, INTENT(IN)           :: bottom
      double precision, INTENT(IN)  :: rmud
      double precision, INTENT(IN)  :: dzRT(jpk)
      double precision, INTENT(IN)  :: CHLz(jpk,nchl),CDOMz(jpk, nlt),NAPz(jpk, nlt)
      double precision, INTENT(IN)  :: Edtop(nlt),Estop(nlt)
      double precision, INTENT(OUT) :: Edz(jpk,nlt),Esz(jpk,nlt)
      double precision, INTENT(OUT) :: Euz(jpk,nlt)
      double precision, INTENT(OUT) :: PARz(jpk,nchl+1)
      double precision, intent(OUT) :: Eutop(nlt)
      double precision              :: lambda, saw, sbw
!     Local variables
      double precision :: Etop
      double precision :: Plte
      double precision :: actot(jpk,nlt),bctot(jpk,nlt),bbctot(jpk,nlt) 
      double precision :: anap(nlt)
      double precision :: a(jpk,nlt), bt(jpk,nlt), bb(jpk,nlt) 
      double precision :: bbc(4)
      data bbc /0.002d0, 0.00071d0, 0.001955d0, 0.0029d0/
      double precision bbw
      data bbw /0.5d0/       !backscattering to forward scattering ratio

!  Water data files
      allocate(aw1(bottom, nlt))
      allocate(bw1(bottom, nlt))

      !cfle = cdir//cabw
      open(4,file='bcs/abw25.dat',status='old',form='formatted')
      !do i = 1,5
      ! read(4,'(a80)')title
      ! if(lwp)  write(6,'(a80)')title
      !enddo
      do jk = 1, bottom
        do nl = 1,nlt
          read(4,20)lambda,saw,sbw
          !if (lwp) write(6,20)lambda,saw,sbw
          lam(nl)     = lambda
          aw1(jk,nl)  = saw
          bw1(jk,nl)  = sbw
      enddo
      enddo
      close(4)
20    format(i5,f15.4,f10.4)
 
!  Compute irradiance with depth
       anap(:)  = 0.0d0

       do nl = 1,nlt
          do jk = 1,bottom

          Edz(jk,nl)    = 0.0d0
          Esz(jk,nl)    = 0.0d0
          Euz(jk,nl)    = 0.0d0
          actot(jk,nl)  = 0.0d0
          bctot(jk,nl)  = 0.0d0
          bbctot(jk,nl) = 0.0d0

          do n = 1,nchl
               actot(jk,nl)  = actot(jk,nl)  + CHLz(jk,n)*ac(n,nl)
               bctot(jk,nl)  = bctot(jk,nl)  + CHLz(jk,n)*bc(n,nl)
               bbctot(jk,nl) = bbctot(jk,nl) + CHLz(jk,n)*bbc(n)*bc(n,nl)
          enddo

          a(jk,nl)  = aw1(jk, nl) + CDOMz(jk, nl) + NAPz(jk,nl) + actot(jk,nl)
          bt(jk,nl) = bw1(jk, nl) + bctot(jk,nl)
          bb(jk,nl) = bbw*bw1(jk, nl) + bbctot(jk,nl)
          bb(jk,nl) = max(bb(jk,nl),0.0002d0)

         enddo
       enddo
       
       if (MODE == 1) then ! condition on level < 200 mt
           call radmod(V_POSITION,bottom,dzRT(:),Edtop,Estop,rmud,a,bt,bb,Edz,Esz,Euz,Eutop,PARz)
       endif 
! condition on level > 200 mt
       if ( MODE ==0) then
           call radmod_ex(V_POSITION,bottom,dzRT(:),Edtop,Estop,rmud,a,bt,bb,Edz,Esz,Euz,Eutop,PARz)
       endif


      return
      end
