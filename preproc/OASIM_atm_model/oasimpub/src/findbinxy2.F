      subroutine findbin(igrd,jgrd,xlon,ylat,ibin,jbin)
c  
c  Given an Earth grid, finds the bin of a particular lon/lat location
c  and returns the i, j location on the grid.
c  This version re-computes dlat, dlon everytime to allow its use
c  several times by the same program.
c
c
c  Constants
      dlonr = float(igrd)/360.0
      dlatr = float(jgrd)/180.0
c
      ibin = int((xlon+180.0)*dlonr + 1.0)
      jbin = min(int((ylat+90.0)*dlatr + 1.0),jgrd)
      if (ibin .gt. igrd)ibin = 1
c
      return
      end
c
c
      subroutine findxy(igrd,jgrd,ibin,jbin,xlon,ylat)
c  
c  Given an Earth grid, finds lon/lat from the bin of a particular 
c  location.
c
c  Constants
      dlonr = float(igrd)/360.0
      dlatr = float(jgrd)/180.0
c
      xlon = (float(ibin)-1.0)/dlonr - (180.0-0.5/dlonr)
      ylat = (float(jbin)-1.0)/dlatr - (90.0-0.5/dlatr)
c
      return
      end
