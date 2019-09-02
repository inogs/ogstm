      subroutine getarea(area)
c
c  Returns the area in m2 of each grid point.
c
#include "define.h"
      parameter (nrg=13+1)
      real area(jm)
      real*8 pi,rad
c
      pi = dacos(-1.0D0)
      rad = 180.0D0/pi
c
c  Area parameters
      xfac = 1.0/dlon
      yfac = 1.0/dlat
      xinc = 1.0/xfac*0.5
      yinc = 1.0/yfac*0.5
c
c  Earth constants
      f = 1.0/298.257
      Re = 6378.137*1000.0   !m
      Rp = -(f*Re-Re)
      arg = 2.0*f - f*f
      ecc = sqrt(arg)
      omf = 1.0-f
c
      ylnth = (111.138*1000.0)/yfac
      do j = 1,jm
       ylt = abs(float(j)/yfac-(abs(alat0)+yinc))   !latitude in degrees
       y90 = 90.0 - ylt     !rearrange as angle from due north
c   Convert to geocentric latitude
       tanc = omf*omf*tan(ylt/rad)
       rlatc = atan(tanc)
       cosyc = cos(rlatc)
       rnum = Re*omf
       arg = 1.0 - (2.0-f)*f*cosyc*cosyc
       rden = sqrt(arg)
       Reloc = rnum/rden    !Earth radius at local point
       a = Reloc*sin(y90/rad)  !radius of circle at latitude
       xlnth = (2.0*a*pi/360.0)/xfac  !length of longitude
       area(j) = xlnth*ylnth
      enddo
c
      return
      end
