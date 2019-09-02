      subroutine lidata
c
c  Reads in radiative transfer data: specifically atmospheric data
c  (Extraterrestrial irradiance and atmospheric optical data), and
c  water data (seawater absorption and total scattering coefficients,
c  and chl-specific absorption and total scattering data for 
c  several phytoplankton groups).  PAR (350-700) begins at index 3, 
c  and ends at index 17.
c
#include "comlte.h"
      character*80 title
      character*50 cfle
      character catmo*11,cacbc*11,cabw*10
      real sFobar,sthray,soza,sawv,sao,saco2
      real saw,sbw,sac,sbc
      character*4 cdir
      data cdir /'bcs/'/
      data catmo,cacbc,cabw /'atmo25b.dat','acbc25b.dat','abw25b.dat'/
c
c  Atmospheric data file
      cfle = cdir//catmo
      open(4,file=cfle,status='old',form='formatted')
      read(4,'(a80)')title
      read(4,'(a80)')title
      read(4,'(a80)')title
      read(4,'(a80)')title
      np = 0
      do nl = 1,nlt
       read(4,10)ilam,sFobar,sthray,soza,sawv,sao,saco2
       write(6,10)ilam,sFobar,sthray,soza,sawv,sao,saco2
       lam(nl) = ilam
       Fobar(nl) = sFobar
       thray(nl) = sthray
       oza(nl) = soza
       awv(nl) = sawv
       ao(nl) = sao
       aco2(nl) = saco2
      enddo
      close(4)
c
10    format(i5,6f11.4)
c
c
c  Water data files
      cfle = cdir//cabw
      open(4,file=cfle,status='old',form='formatted')
      do i = 1,6
       read(4,'(a50)')title
c       write(6,'(a50)')title
      enddo
      np = 0
      do nl = 1,nlt
       read(4,20)lambda,saw,sbw
c       write(6,20)lambda,saw,sbw
       aw(nl) = saw
       bw(nl) = sbw
      enddo
      close(4)
20    format(i5,f15.4,f10.4)
c
      return
      end
