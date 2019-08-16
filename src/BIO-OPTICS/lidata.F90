      subroutine lidata()
!     subroutine lidata(lam,aw,bw,ac,bc)

      USE OPT_mem, ONLY: nlt, nchl, lam,aw,bw,ac,bc
      IMPLICIT NONE 
!  Reads in radiative transfer data: specifically 
!  water data (seawater absorption and total scattering coefficients,
!  and chl-specific absorption and total scattering data for 
!  several phytoplankton groups).  PAR (350-700) begins at index 3, 
!  and ends at index 17.
 
!     parameter(nlt=33)
!     parameter(nchl=5)
      character*80 title
      character*80 cfle
      character cacbc*11,cabw*15
      double precision saw,sbw,sac,sbc
      character*4 cdir
!     integer lam(nlt)
!     double precision aw(nlt),bw(nlt)
!     double precision ac(nchl,nlt),bc(nchl,nlt)
      data cdir /'bcs/'/
      data cacbc,cabw /'acbc25b.dat','abw25_morel.dat'/
      integer    :: i, n, nl,lambda
 
!  Water data files
      cfle = cdir//cabw
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,5
       read(4,'(a80)')title
!       write(6,'(a80)')title
      enddo
      do nl = 1,nlt
       read(4,20)lambda,saw,sbw
!       write(6,20)lambda,saw,sbw
       lam(nl) = lambda
       aw(nl) = saw
       bw(nl) = sbw
      enddo
      close(4)
20    format(i5,f15.4,f10.4)
 
!  Phytoplankton group chl-specific absorption and total scattering 
!  data.  Chl-specific absorption data is normalized to 440 nm; convert
!  here to actual ac*(440)

      cfle = cdir//cacbc
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,6
       read(4,'(a80)')title
      enddo
      do n = 1,nchl
       read(4,'(a80)')title
       do nl = 1,19
        read(4,30)lambda,sac,sbc
        ac(n,nl) = sac
        bc(n,nl) = sbc
       enddo
       do nl = 20,nlt
        ac(n,nl) = 0.0
        bc(n,nl) = 0.0
       enddo
      enddo
      close(4)
30    format(i4,2f10.4)
 
      return
      end
