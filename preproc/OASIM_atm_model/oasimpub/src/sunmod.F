      subroutine sunmod(radeg,iday,iyr,gmt,il,jl,sunz,rs)
c
c  Given year, day of year, time in hours (GMT) and latitude and 
c  longitude, returns an accurate solar zenith and azimuth angle.  
c  Based on IAU 1976 Earth ellipsoid.  Method for computing solar 
c  vector and local vertical from Patt and Gregg, 1994, Int. J. 
c  Remote Sensing.  Only outputs solar zenith angle.  This version
c  utilizes a pre-calculation of the local up, north, and east 
c  vectors, since the locations where the solar zenith angle are 
c  calculated in the model are fixed.
c
c  Subroutines required: sun2000
c                        gha2000
c                        jd
c
#include "define.h"
      real suni(3),sung(3)
      real no
      common /bvec/ up(im,jm,3),no(im,jm,3),ea(im,jm,3)
c
c  Compute sun vector
c   Compute unit sun vector in geocentric inertial coordinates
      sec = gmt*3600.0D0
      call sun2000(radeg,iyr,iday,sec,suni,rs)

c   Get Greenwich mean sidereal angle
      day = float(iday) + sec/86400.0D0
      call gha2000(radeg,iyr,day,gha)
      ghar = gha/radeg

c   Transform Sun vector into geocentric rotating frame
      sung(1) = suni(1)*cos(ghar) + suni(2)*sin(ghar)
      sung(2) = suni(2)*cos(ghar) - suni(1)*sin(ghar)
      sung(3) = suni(3)
c
c  Compute components of spacecraft and sun vector in the
c  vertical (up), North (no), and East (ea) vectors frame
      sunv = 0.0
      sunn = 0.0
      sune = 0.0
      do j = 1,3
       sunv = sunv + sung(j)*up(il,jl,j)
       sunn = sunn + sung(j)*no(il,jl,j)
       sune = sune + sung(j)*ea(il,jl,j)
      enddo
c
c  Compute the solar zenith 
      sunz = radeg*atan2(sqrt(sunn*sunn+sune*sune),sunv)
c
      return
      end
c
c *****************************************************************
	subroutine sun2000(radeg,iyr,iday,sec,sunvec,rs)
c
c  This subroutine computes the Sun vector in geocentric inertial 
c  (equatorial) coodinates.  It uses the model referenced in The 
c  Astronomical Almanac for 1984, Section S (Supplement) and documented
c  in Exact closed-form geolocation algorithm for Earth survey
c  sensors, by F.S. Patt and W.W. Gregg, Int. Journal of Remote
c  Sensing, 1993.  The accuracy of the Sun vector is approximately 0.1 
c  arcminute.
c
c	Arguments:
c
c	Name	Type	I/O	Description
c	--------------------------------------------------------
c	IYR	I*4	 I	Year, four digits (i.e, 1993)
c	IDAY	I*4	 I	Day of year (1-366)
c	SEC	R*8	 I	Seconds of day 
c	SUN(3)	R*8	 O	Unit Sun vector in geocentric inertial 
c				 coordinates of date
c	RS	R*8	 O	Magnitude of the Sun vector (AU)
c
c	Subprograms referenced:
c
c	JD		Computes Julian day from calendar date
c	EPHPARMS	Computes mean solar longitude and anomaly and
c			 mean lunar lontitude and ascending node
c	NUTATE		Compute nutation corrections to lontitude and 
c			 obliquity
c
c	Coded by:  Frederick S. Patt, GSC, November 2, 1992
c	Modified to include Earth constants subroutine by W. Gregg,
c		May 11, 1993.


	real sunvec(3)
	common /nutcm/dpsi,eps,nutime
	data xk/0.0056932/		!Constant of aberration 
	data imon/1/

c   Compute floating point days since Jan 1.5, 2000 
c    Note that the Julian day starts at noon on the specified date
       rjd = float(jd(iyr,imon,iday))
	t = rjd - 2451545.0D0 + (sec-43200.0D0)/86400.0D0

c  Compute solar ephemeris parameters
	call ephparms(t,xls,gs,xlm,omega)

c  Check if need to compute nutation corrections for this day
	nt = int(t)
	if (nt.ne.nutime) then
	  nutime = nt
	  call nutate(radeg,t,xls,gs,xlm,omega,dpsi,eps)
	end if

c  Compute planet mean anomalies
c   Venus Mean Anomaly 	
	g2 = 50.40828D0 + 1.60213022D0*t
	g2 = mod(g2,360.0)

c   Mars Mean Anomaly 		
	g4 = 19.38816D0 + 0.52402078D0*t
	g4 = mod(g4,360.0)

c  Jupiter Mean Anomaly 
	g5 = 20.35116D0 + 0.08309121D0*t
	g5 = mod(g5,360.0)

c  Compute solar distance (AU)
	rs = 1.00014D0 - 0.01671D0*cos(gs/radeg) 
     *       - 0.00014D0*cos(2.0D0*gs/radeg)

c  Compute Geometric Solar Longitude 
	dls = 	(6893.0D0 - 4.6543463D-4*t)*sin(gs/radeg) 
     *		+ 72.0D0*sin(2.0D0*gs/radeg) 
     *		- 7.0D0*cos((gs - g5)/radeg)
     *		+ 6.0D0*sin((xlm - xls)/radeg) 
     *		+ 5.0D0*sin((4.0D0*gs - 8.0D0*g4 + 3.0D0*g5)/radeg) 
     *		- 5.0D0*cos((2.0D0*gs - 2.0D0*g2)/radeg)
     *		- 4.0D0*sin((gs - g2)/radeg) 
     *		+ 4.0D0*cos((4.0D0*gs - 8.0D0*g4 + 3.0D0*g5)/radeg) 
     *		+ 3.0D0*sin((2.0D0*gs - 2.0D0*g2)/radeg)
     *		- 3.0D0*sin(g5/radeg) 
     *		- 3.0D0*sin((2.0D0*gs - 2.0D0*g5)/radeg)  !arcseconds

	xlsg = xls + dls/3600.0D0

c  Compute Apparent Solar Longitude; includes corrections for nutation 
c   in longitude and velocity aberration
	xlsa = xlsg + dpsi - xk/rs

c   Compute unit Sun vector 
	sunvec(1) = cos(xlsa/radeg)
	sunvec(2) = sin(xlsa/radeg)*cos(eps/radeg)
	sunvec(3) = sin(xlsa/radeg)*sin(eps/radeg)
c	type *,' Sunlon = ',xlsg,xlsa,eps

	return
	end
c
c *********************************************************************
	subroutine gha2000(radeg,iyr,day,gha)

c  This subroutine computes the Greenwich hour angle in degrees for the
c  input time.  It uses the model referenced in The Astronomical Almanac
c  for 1984, Section S (Supplement) and documented in Exact 
c  closed-form geolocation algorithm for Earth survey sensors, by 
c  F.S. Patt and W.W. Gregg, Int. Journal of Remote Sensing, 1993.
c  It includes the correction to mean sideral time for nutation
c  as well as precession.

c  Calling Arguments

c  Name		Type 	I/O	Description
c
c  iyr		I*4	 I	Year (four digits)
c  day		R*8	 I	Day (time of day as fraction)
c  gha		R*8	 O	Greenwich hour angle (degrees)


c	Subprograms referenced:
c
c	JD		Computes Julian day from calendar date
c	EPHPARMS	Computes mean solar longitude and anomaly and
c			 mean lunar lontitude and ascending node
c	NUTATE		Compute nutation corrections to lontitude and 
c			 obliquity
c 	
c
c	Program written by:	Frederick S. Patt
c				General Sciences Corporation
c				November 2, 1992
c
c	Modification History:
c
	common /nutcm/dpsi,eps,nutime
	data imon/1/,nutime/-99999/

c  Compute days since J2000
	iday = int(day)
	fday = day - iday
	jday = jd(iyr,imon,iday)
	t = jday - 2451545.5D0 + fday
	
c  Compute Greenwich Mean Sidereal Time	(degrees)
	gmst = 100.4606184D0 + 0.9856473663D0*t + 2.908D-13*t*t

c  Check if need to compute nutation correction for this day
	nt = int(t)
	if (nt.ne.nutime) then
	  nutime = nt
	  call ephparms(t,xls,gs,xlm,omega)
	  call nutate(radeg,t,xls,gs,xlm,omega,dpsi,eps)
	end if

c  Include apparent time correction and time-of-day
	gha = gmst + dpsi*cos(eps/radeg) + fday*360.0D0
	gha = mod(gha,360.0)
	if (gha.lt.0.0D0) gha = gha + 360.0D0

	return
	end

	subroutine ephparms(t,xls,gs,xlm,omega)

c  This subroutine computes ephemeris parameters used by other Mission
c  Operations routines:  the solar mean longitude and mean anomaly, and
c  the lunar mean longitude and mean ascending node.  It uses the model
c  referenced in The Astronomical Almanac for 1984, Section S 
c  (Supplement) and documented and documented in Exact closed-form 
c  geolocation algorithm for Earth survey sensors, by F.S. Patt and 
c  W.W. Gregg, Int. Journal of Remote Sensing, 1993.  These parameters 
c  are used to compute the solar longitude and the nutation in 
c  longitude and obliquity.

c  Calling Arguments

c  Name		Type 	I/O	Description
c
c  t		R*8	 I	Time in days since January 1, 2000 at 
c				 12 hours UT
c  xls		R*8	 O	Mean solar longitude (degrees)
c  gs		R*8	 O	Mean solar anomaly (degrees)
c  xlm		R*8	 O	Mean lunar longitude (degrees)
c  omega	R*8	 O	Ascending node of mean lunar orbit 
c				 (degrees)
c
c
c	Program written by:	Frederick S. Patt
c				General Sciences Corporation
c				November 2, 1992
c
c	Modification History:
c

c  Sun Mean Longitude 		
	xls = 280.46592D0 + 0.9856473516D0*t
	xls = mod(xls,360.0)
 
c  Sun Mean Anomaly		
	gs = 357.52772D0 + 0.9856002831D0*t 
	gs = mod(gs,360.0)

c  Moon Mean Longitude		
	xlm = 218.31643D0 + 13.17639648D0*t 
	xlm = mod(xlm,360.0)

c  Ascending Node of Moons Mean Orbit 	
	omega = 125.04452D0 - 0.0529537648D0*t 
	omega = mod(omega,360.0)
	
	return
	end

	subroutine nutate(radeg,t,xls,gs,xlm,omega,dpsi,eps)

c  This subroutine computes the nutation in longitude and the obliquity
c  of the ecliptic corrected for nutation.  It uses the model referenced
c  in The Astronomical Almanac for 1984, Section S (Supplement) and 
c  documented in Exact closed-form geolocation algorithm for Earth 
c  survey sensors, by F.S. Patt and W.W. Gregg, Int. Journal of 
c  Remote Sensing, 1993.  These parameters are used to compute the 
c  apparent time correction to the Greenwich Hour Angle and for the 
c  calculation of the geocentric Sun vector.  The input ephemeris 
c  parameters are computed using subroutine ephparms.  Terms are 
c  included to 0.1 arcsecond.

c  Calling Arguments

c  Name		Type 	I/O	Description
c
c  t		R*8	 I	Time in days since January 1, 2000 at 
c				 12 hours UT
c  xls		R*8	 I	Mean solar longitude (degrees)
c  gs		R*8	 I	Mean solar anomaly   (degrees)
c  xlm		R*8	 I	Mean lunar longitude (degrees)
c  Omega	R*8	 I	Ascending node of mean lunar orbit 
c  				 (degrees)
c  dPsi		R*8	 O	Nutation in longitude (degrees)
c  Eps		R*8	 O	Obliquity of the Ecliptic (degrees)
c 				 (includes nutation in obliquity)
c
c
c	Program written by:	Frederick S. Patt
c				General Sciences Corporation
c				October 21, 1992
c
c	Modification History:
c
	
c  Nutation in Longitude
	dpsi = - 17.1996D0*sin(omega/radeg) 
     *	 	+ 0.2062D0*sin(2.0D0*omega/radeg)
     *	     	- 1.3187D0*sin(2.0D0*xls/radeg) 
     *		+ 0.1426D0*sin(gs/radeg) 
     *		- 0.2274D0*sin(2.0D0*xlm/radeg) 
	
c  Mean Obliquity of the Ecliptic	
	epsm = 23.439291D0 - 3.560D-7*t 

c  Nutation in Obliquity 
	deps = 9.2025D0*cos(omega/radeg) + 0.5736D0*cos(2.0D0*xls/radeg)

c  True Obliquity of the Ecliptic 
	eps = epsm + deps/3600.0D0

	dpsi = dpsi/3600.0D0

	return
	end
c
c ************************************************************************
      function jd(i,j,k)
c
c
c    This function converts a calendar date to the corresponding Julian
c    day starting at noon on the calendar date.  The algorithm used is
c    from Van Flandern and Pulkkinen, Ap. J. Supplement Series 41, 
c    November 1979, p. 400.
c
c
c	Arguments
c     
c     	Name    Type 	I/O 	Description
c     	----	---- 	--- 	-----------
c     	i	I*4  	 I 	Year - e.g. 1970
c     	j       I*4  	 I  	Month - (1-12)
c     	k       I*4  	 I  	Day  - (1-31)
c     	jd      I*4  	 O  	Julian day
c
c     external references
c     -------------------
c      none
c
c
c     Written by Frederick S. Patt, GSC, November 4, 1992
c
c
      	jd = 367*i - 7*(i+(j+9)/12)/4 + 275*j/9 + k + 1721014

c  This additional calculation is needed only for dates outside of the 
c   period March 1, 1900 to February 28, 2100
c     	jd = jd + 15 - 3*((i+(j-9)/7)/100+1)/4
      	return
      	end
