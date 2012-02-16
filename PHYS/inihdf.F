CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/inihdf.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE inihdf
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE inihdf
CCC                     ******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Initializations of the horizontal coefficients
CCC
CC   Method :
CC   -------
CC
CC	Eddy diffusivity and eddy induced velocity cefficients:
CC         default option   : constant coef. aht0, aeiv0 (namelist)
CC        'key_trahdfcoef1d': depth dependent coef. defined in 
CC                            in 'inihdf.tra.coef3d.h' file
CC        'key_trahdfcoef2d': latitude and longitude dependent coef.
CC                            defined in 'inihdf.tra.coef3d.h' file
CC        'key_trahdfcoef3d': latitude, longitude, depth dependent coef.
CC                            defined in 'inihdf.tra.coef3d.h' file
CC
CC	Direction of lateral diffusion (tracers and/or momentum)
CC	  'key_trahdfiso'   : initialize the slope arrays to zero
CC	  'key_trahdfgeop'  : initialise the slope arrays to the i- and
CC                            j-slopes of s-surfaces
CC	   
CC	Biharmonic operator:
CC	  'key_trahdfbilap' and 'key_trahdfgeop' or 'key_trahdfiso'
CC	                    : initialise working arrays to zero
CC
CC	N.B. User defined include files.  By default, 3d and 2d coef.
CC	are set to a constant value given in the namelist and the 1d
CC	coefficients are initialized to a hyperbolic tangent vertical
CC	profile.
CC
CC   Input :
CC   ------
CC      common
CC            /comask/          : masks, bathymetry
CC            /comcoh/          : scale factors
CC            /comhdt/          : lateral eddy diffusivity
CC
CC   Output :
CC   -------
CC      common
CC            /comhdt/ ahtt,ahtu, : lateral eddy diffusivity coef.
CC                     ahtv,ahtw
CC	      /comdhd/ uslp ,vslp : i- and j- direction of isopycnal or
CC		       wslpi,wslpj  s-surfaces
CC		       wk1, wk2,  : work arrays used in rotated bihar-
CC		       wk3, wk4     monic momentum or tracer operator
CC
CC   External :        prihre, mpplnk or lbc
CC   ---------
CC
CC   Reference :
CC   ----------
CC      Madec, G. and M. Imbard, 1996, A global ocean mesh to overcome
CC      the North Pole singularity, Climate Dynamics, 12, 381-388.
CC
CC   Modifications :
CC   --------------
CC      original :  90-10 (C. Levy - G. Madec)
CC      additions : 96-01 (MA Foujols) adaptation of inimix for offline
CC      additions : 98-04 (L. Bopp MA Foujols) multi dimensional coeff
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

c#include "parameter.h"
c#include "common.h"
       USE myalloc

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji, jj, jk
C
#ifdef  key_trahdfcoef3d
C    local definition for 3d coefficients
#  elif key_trahdfcoef2d
C    local definition for 2d coefficients
#  elif key_trahdfcoef1d
C    local definition for 1d coefficients
      REAL(8) zkah, zahr, za00, za01, zahf, zahs
      REAL(8) zahtf, zahts
#endif
CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"

CCC---------------------------------------------------------------------
CCC  OPA8.1, LODYC (1998)
CCC---------------------------------------------------------------------
C
C    
C
C 1. Lateral eddy diffusivity and eddy induced velocity coefficients
C ==================================================================
C     user: modify the 'inihdf.tra.coef_d.h' include files following
C           your desiderata
#ifdef key_trahdfcoef3d
C
C 3D coefficients : ahtt, aeiv function of longitude, latitude, depth 
C ----------------
#     include "inihdf.tra.coef3d.h"
C
#  elif key_trahdfcoef2d
C
C 2D coefficients : ahtt and aeiv function of longitude and latitude
C ----------------
#     include "inihdf.tra.coef2d.h"
C
#  elif key_trahdfcoef1d
C
C 1D coefficients : ahtt and aeiv function of depth 
C ----------------
#     include "inihdf.tra.coef1d.h"
C
#  else
C
C Constant coefficients :
C ----------------------
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' inihdf: constant eddy diffusivity coef.'
      IF(lwp)WRITE(numout,*) ' ======  '
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '        aht0 = ',aht0
#    ifdef key_trahdfeiv
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) ' inihdf: constant eddy induced velocity',
     $                       ' coef.'
      IF(lwp)WRITE(numout,*) ' ======  '
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '         aeiv0 = ',aeiv0
#      else
C No tracer eddy induced velocity parameterization
#    endif
C
#endif
C
C
C 3. Direction of lateral diffusion (passive tracers)
C =================================
C
#if defined key_trahdfiso 
C
C 3.1 isopycnal diffusion on tracers
C ----------------------------------
C
C ... set the slope of diffusion to zero (even in s-coordinates)
      DO jk = 1, jpk
        DO jj = 1, jpj
          DO ji = 1, jpi
            uslp (ji,jj,jk) = 0.e0
            vslp (ji,jj,jk) = 0.e0
            wslpi(ji,jj,jk) = 0.e0
            wslpj(ji,jj,jk) = 0.e0
            bn2n(ji,jj,jk)  = 0.e0
            rdn(ji,jj,jk)   = 0.e0
          END DO
        END DO
      END DO
C
#  elif defined key_trahdfgeop 
C
C 3.2 geopotential diffusion in s-coordinates on tracers
C ------------------------------------------------------
C The slopes of s-surfaces are computed once (no call to hdfslp in step)
C The slopes for momentum diffusion are i- or j- averaged of those on
C tracers
C
C ... set the slope of diffusion to the slope of s-surfaces
C    ( c a u t i o n : minus sign as fsdep has positive value )
      DO jk = 1, jpk
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            uslp(ji,jj,jk) =-1. / e1u(ji,jj) * umask(ji,jj,jk)
     $                     * ( fsdept(ji+1,jj,jk) - fsdept(ji,jj,jk) )
            vslp(ji,jj,jk) =-1. / e2v(ji,jj) * vmask(ji,jj,jk)
     $                     * ( fsdept(ji,jj+1,jk) - fsdept(ji,jj,jk) )
            wslpi(ji,jj,jk)=-1. / e1t(ji,jj) * tmask(ji,jj,jk)
     $                     * ( fsdepuw(ji+1,jj,jk) - fsdepuw(ji,jj,jk) )
            wslpj(ji,jj,jk)=-1. / e2t(ji,jj) * tmask(ji,jj,jk)
     $                     * ( fsdepvw(ji,jj+1,jk) - fsdepvw(ji,jj,jk) )
          END DO
        END DO
      END DO
C
C ... Lateral boundary conditions on the slopes
#    ifdef key_mpp
C   ... Mpp: export boundary values to neighbouring processors
      CALL mpplnk( uslp , 2, 0 )
      CALL mpplnk( vslp , 3, 0 )
      CALL mpplnk( wslpi, 1, 0 )
      CALL mpplnk( wslpj, 1, 0 )
#      else
C ... mono- or macro-tasking: U-, V-, W-pts, sgn, 3D , no slab
      CALL lbc( uslp , 2, 0, 1, 1, jpk, 1 )
      CALL lbc( vslp , 3, 0, 1, 1, jpk, 1 )
      CALL lbc( wslpi, 1, 0, 1, 1, jpk, 1 )
      CALL lbc( wslpj, 1, 0, 1, 1, jpk, 1 )
#    endif
C
#  else
C
C tracer and momentum diffusion along model coordinates
C
#endif
#if defined key_trahdfbilap
C
C 4. Biharmonic initializations
C -----------------------------
C
#    if defined key_trahdfgeop || defined key_trahdfiso
        DO jk = 1, jpk
          DO jj = 1, jpj
            DO ji = 1, jpi
              wk1(ji,jj,jk) = 0.e0
              wk2(ji,jj,jk) = 0.e0
              wk3(ji,jj,jk) = 0.e0
              wk4(ji,jj,jk) = 0.e0
            END DO
          END DO
        END DO 
C
#      else 
C        No use of wk1, wk2, wk3 and wk4 arrays: no rotation required 
#    endif 
C
#endif 
C
      RETURN
      END
