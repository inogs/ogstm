CCC $Header: /cvsroot/opatm-bfm/opa_model/OPA/trcbbl.offline.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcbbl.offline
CCC                     ************************
CCC
CCC  Purpose :
CCC  --------
CCC     Compute the before tracer (t & s) trend associated with the
CCC	bottom boundary layer and add it to the general trend of tracer
CCC	equations. The bottom boundary layer is supposed to be a purely
CCC     diffusive bottom boundary layer.
CCC
CC
CC
CC   Method :
CC   -------
CC	When the product grad( rho) * grad(h) < 0 (where grad is a
CC	along bottom slope gradient) an additional lateral 2nd order
CC	diffusion along the bottom slope is added to the general
CC	tracer trend, otherwise the additional trend is set to 0.
CC      Second order operator (laplacian type) with variable coefficient
CC      computed as follow for temperature (idem on s): 
CC         difft = 1/(e1t*e2t*e3t) { di-1[ ahbt e2u*e3u/e1u di[ztb] ]
CC                                 + dj-1[ ahbt e1v*e3v/e2v dj[ztb] ] }
CC	where ztb is a 2D array: the bottom ocean te;perature and ahtb
CC	is a time and space varying diffusive coefficient defined by:
CC	   ahbt = zahbp    if grad(rho).grad(h) < 0
CC	        = 0.       otherwise.
CC	Note that grad(.) is the along bottom slope gradient. grad(rho)
CC	is evaluated using the local density (i.e. referenced at the
CC	local depth). Typical value of ahbt is 2000 m2/s (equivalent to
CC	a downslope velocity of 20 cm/s if the condition for slope
CC	convection is satified)
CC
CC      Add this before trend to the general trend (ta,sa) of the 
CC	botton ocean tracer point:
CC              ta = ta + difft
CC
CC      'key_diatrends' defined, the trends are saved for diagnostics.
CC
CC      not macro-tasked
CC
CC   Input :
CC   ------
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      common
CC            /comcoo/          : mesh and scale factors
CC            /comask/          : masks, bathymetry
CC            /combef/          : previous fields (before)
CC            /comaft/          : next fields (after)
CC            /comtsk/          : multitasking
CC
CC   Output :
CC   -------
CC      common
CC            /comaft/ ta, sa   : general tracer trend increased by the
CC                                before horizontal diffusion trend
CC            /comtrd/ ttrd,strd: now horizontal tracer diffusion trend
CC                                (if 'key_diatrends' defined)
CC
CC   References :                 NO
CC   -----------
CC	Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
CC
CC   Modifications :
CC   --------------
CC      original :  96-06 (L. Mortier)
CC      addition :  97-11 (G. Madec)
CC      addition :  01-05 (O. Aumont, E. Kestenare)
CC      additions : 01-06 add key_off_tra for trcbbl routine (ED,EK)
CC----------------------------------------------------------------------

      INTEGER ji, jj, jn
      INTEGER ik, iku, ikv
      REAL(8) zahu(jpi,jpj), zahv(jpi,jpj)

C
      REAL(8) zbtr, zta
      REAL(8) zki(jpi,jpj), zkj(jpi,jpj)
      REAL(8) zkw(jpi,jpj), zkx(jpi,jpj), zky(jpi,jpj), zkz(jpi,jpj)
      REAL(8) ztbb(jpi,jpj)
C
CC----------------------------------------------------------------------
CC statement functions
CC ===================

!                              #include "stafun.h"

       trcbblparttime = MPI_WTIME() !F79 cronometer-start


C 0. 2D fields of bottom temperature and salinity, and bottom slope
C -----------------------------------------------------------------
C mbathy= number of w-level, minimum value=1 (cf dommsk.F)
C
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            iku = max( min( mbathy(ji+1,jj)-1, mbathy(ji,jj)-1 ), 1 )
            ikv = max( min( mbathy(ji,jj+1)-1, mbathy(ji,jj)-1 ), 1 )
            zahu(ji,jj) = atrbbl*e2u(ji,jj)*fse3u(ji,jj,iku)/e1u(ji,jj)
     $                       *umask(ji,jj,1)
            zahv(ji,jj) = atrbbl*e1v(ji,jj)*fse3v(ji,jj,ikv)/e2v(ji,jj)
     $                       *vmask(ji,jj,1)
         END DO
       END DO

      DO jn=1,jptra 
      DO jj = 1, jpj
        DO ji = 1, jpi
C    ... index of the bottom ocean T-level
          ik = max( mbathy(ji,jj)-1, 1 )
C    ... masked before T and S at the ocean bottom 
          ztbb(ji,jj) = trb(ji,jj,ik,jn) * tmask(ji,jj,1)
C    ... depth of the ocean bottom T-level
        END DO
      END DO
C
C
C
C
C
C 2. Additional second order diffusive trends
C -------------------------------------------
C
C ... first derivative (gradient)
        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zkx(ji,jj) = bblx(ji,jj)*zahu(ji,jj)*
     .        ( ztbb(ji+1,jj) - ztbb(ji,jj) )
            zky(ji,jj) = bbly(ji,jj)*zahv(ji,jj)*
     .        ( ztbb(ji,jj+1) - ztbb(ji,jj) )
        END DO
      END DO
C
C ... second derivative (divergence) and add to the general tracer trend
C
        DO jj = 2,jpjm1
          DO ji = 2,jpim1
            ik = max( mbathy(ji,jj)-1, 1 )
            zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,ik) )
            zta = (  zkx(ji,jj) - zkx(ji-1, jj ) 
     $             + zky(ji,jj) - zky( ji ,jj-1)  ) * zbtr
            tra(ji,jj,ik,jn) = tra(ji,jj,ik,jn) + zta
        END DO
      END DO

      END DO
CCC 10 11 2004  F79 cronometer-stop

       trcbblparttime = MPI_WTIME() - trcbblparttime
       trcbbltottime = trcbbltottime + trcbblparttime

CC-CC       write(*,*) "F79T:trcbblparttime", trcbblparttime
CC-CC       write(*,*) "F79T:trcbbltottime", trcbbltottime

CCC

