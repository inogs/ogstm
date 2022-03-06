      subroutine edeseu2(MODE,V_POSITION,bottom,zgrid,Edtop,Estop,CHLz,CDOMz,POCz,rmud,E,PARz)
      USE myalloc
      USE mpi
      USE OPT_mem
      use adj_3stream, only: solve_direct
      IMPLICIT NONE
!  Model of irradiance in the water column.  Accounts for three 
!  irradiance streams:
 
!  Edz = direct downwelling irradiance
!  Esz = diffuse downwelling irradiance
!  Euz = diffuse upwelling irradiance
 
!  Propagation is done in energy units, tests are done in quanta,
!  final is quanta for phytoplankton growth.
 
      integer :: n, nl, jk 
      integer, INTENT(IN)          :: MODE
      character(*), INTENT(IN)     :: V_POSITION
      integer, INTENT(IN)          :: bottom
      double precision, INTENT(IN) :: rmud
      double precision, INTENT(IN) :: zgrid(jpk+1)
      double precision, INTENT(IN) :: CHLz(jpk,nchl),CDOMz(jpk),POCz(jpk)
      double precision, INTENT(IN) :: Edtop(nlt),Estop(nlt)
      double precision, INTENT(OUT) :: PARz(jpk,nchl+1)
      double precision, INTENT(OUT) :: E(3,jpk+1,nlt)
!     double precision, intent(OUT) :: Eutop(nlt)
!     Local variables
      double precision :: Etop
      double precision :: Plte
      double precision :: actot(jpk,nlt),bctot(jpk,nlt),bbctot(jpk,nlt) 
      double precision :: a(jpk,nlt), bt(jpk,nlt), bb(jpk,nlt) 
      double precision :: bbc(4)
      double precision :: rd, rs, ru, vs, vu
      double precision :: vd(jpk,nlt)
      double precision :: E_ave(3,jpk,nlt)
      data bbc /0.002d0, 0.00071d0, 0.001955d0, 0.0029d0/
      double precision bbw
      data bbw /0.5d0/       !backscattering to forward scattering ratio
 
      vd(:,:)=1.0D0/rmud
      rd=1.0D0
      rs=1.5D0
      ru=3.0D0
      vs=0.83D0
      vu=0.4D0

!  Compute irradiance with depth
      E(:,:,:)    = 0.0d0

       do nl = 1,nlt
          do jk = 1,bottom

          actot(jk,nl)  = 0.0d0
          bctot(jk,nl)  = 0.0d0
          bbctot(jk,nl) = 0.0d0

          do n = 1,nchl
               actot(jk,nl)  = actot(jk,nl)  + CHLz(jk,n)*ac(n,nl)
               bctot(jk,nl)  = bctot(jk,nl)  + CHLz(jk,n)*bc(n,nl)
               bbctot(jk,nl) = bbctot(jk,nl) + CHLz(jk,n)*bbc(n)*bc(n,nl)
          enddo

          a(jk,nl)  = aw(nl) + CDOMz(jk) * acdom(nl) + POCz(jk) * apoc(nl) + actot(jk,nl)
          bt(jk,nl) = bw(nl) + bctot(jk,nl) + POCz(jk) * bpoc(nl)
          bb(jk,nl) = bbw*bw(nl) + bbctot(jk,nl) + POCz(jk) * bbpoc(nl)
          bb(jk,nl) = max(bb(jk,nl),0.0002d0)

         enddo
       enddo
       
       call solve_direct(bottom+1, zgrid(1:bottom+1), bottom, zgrid(1:bottom+1), nlt, a, bt, bb, & 
                          rd, rs, ru, vd, vs, vu, Edtop,Estop, E, E_ave)

      return
      end
