      SUBROUTINE comp_ssh16 ()
!CC---------------------------------------------------------------------
!CC
!CC                       ROUTINE div
!CC                     ***************
!CC
!CC  Purpose :
!CC  --------
!CC	compute the now horizontal divergence of the velocity field.
!CC
!C   Method :
!C   -------
!C	The now divergence is given by :
!C       * s-coordinate ('key_s_coord' defined)
!C         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
!C       * z-coordinate (default key)
!C         hdivn = 1/(e1t*e2t) [ di(e2u  un) + dj(e1v  vn) ]
!C
!C      Apply lateral boundary conditions on hdivn through a call to
!C      routine mpplnk ('key_mpp' defined) or lbc.
!C
!C      macro-tasked on horizontal slab (jk-loop :  1, jpk-1)
!C
!C   Input :
!C   ------
!C      argument
!C              ktask           : task identificator
!C              kt              : time step
!C      common
!C            /comcoh/   	: scale factors
!C            /comtsk/   	: multitasking
!C            /comnow/		: present fields (now)
!C
!C   Output :
!C   -------
!C      common
!C	      /comnow/ hdivn	: now horizontal divergence
!C
!C   External : mpplnk or lbc
!C   ---------
!C
!C   Modifications :
!C   --------------
!C      original :  87-06 (P. Andrich, D. L Hostis)
!C      additions : 91-11 (G. Madec)
!C                : 93-03 (M. Guyon) symetrical conditions
!C                : 96-01 (G. Madec) s-coordinates
!C                : 97-06 (G. Madec) lateral boundary cond., lbc
!C----------------------------------------------------------------------
!C parameters and commons
!C ======================
              USE modulo16
              USE PolynomialRoots
!C+CC  Implicit typing is never allowed
        IMPLICIT NONE
!C+CC  Implicit typing is never allowed

!C----------------------------------------------------------------------
!C local declarations
!C ==================
      INTEGER ji16, jj16, jk,i
      INTEGER jpi16m1,jpj16m1,jpk16m1
!
      REAL(8) zwu16, zwv16, zww16, inv_fact
      REAL(8) a,b,c,x_0,y_0
      REAL(8) resx, resy, x_1, y_1, du, dv, dz
      REAL(8) resx_1, resy_1
      REAL(8) aa(5)
      REAL(8) u,v,e3
      COMPLEX(8) z(4),zz(4)
!C----------------------------------------------------------------------
!C statement functions
!C ===================
!CC---------------------------------------------------------------------
!CC  OPA8, LODYC (1997)
!CC---------------------------------------------------------------------
!
! Horizontal slab
! ===============
!
      jpi16m1=jpi16-1
      jpj16m1=jpj16-1

      un16(jpi16,:,1)=0; un16(:,jpj16,1)=0;
      vn16(jpi16,:,1)=0; vn16(:,jpj16,1)=0;
      wn16(jpi16,:,1)=0; wn16(:,jpj16,1)=0;

!
        DO jj16 = 1,jpj16
          DO ji16 = 1,jpi16
               ssh16_R8(ji16,jj16)  = REAL(e3t16(ji16,jj16,1),8)
               ssh16x_R8(ji16,jj16) = REAL(e3t16(ji16,jj16,1),8)
               ssh16y_R8(ji16,jj16) = REAL(e3t16(ji16,jj16,1),8)
          END DO  
        END DO  
!
!
! 1. Horizontal fluxes
! --------------------
!
        DO jj16 = jpj16m1,2,-1
          DO ji16 = jpi16m1,2,-1

! Compute the optimal flux needed
!           zwu16    = REAL(e2u16(ji16,jj16) * ssh16x_R8(ji16+1,jj16) * un16(ji16,jj16,1),8)
!           zwv16    = REAL(e1v16(ji16,jj16) * ssh16y_R8(ji16,jj16+1) * vn16(ji16,jj16,1),8)
            zwu16    = REAL(e2u16(ji16,jj16) * 3. * un16(ji16,jj16,1),8)
            zwv16    = REAL(e1v16(ji16,jj16) * 3. * vn16(ji16,jj16,1),8)
            zww16    = REAL(e1t16(ji16,jj16) * e2t16(ji16,jj16)       * wn16(ji16,jj16,2),8)

            a        = -REAL(e2u16(ji16-1,jj16) ,8)
            b        = -REAL(e1v16(ji16,jj16-1) ,8)
            c        = zwu16 + zwv16 - zww16

            x_0      = 3.0 * REAL( un16(ji16-1,jj16,1),8)! average first cell flux
            y_0      = 3.0 * REAL( vn16(ji16,jj16-1,1),8)! average first cell flux

       if ( ABS(x_0*a+y_0*b) .GT. 0) then
            write(*,"(9999(G12.5,:,','))") 'ji16= ',ji16,'jj16= ',jj16
            write(*,"(9999(G12.5,:,','))") 'fx_0= ',x_0*a,'fy_0= ',y_0*b
            write(*,"(9999(G12.5,:,','))") 'fx_1= ',zwu16,'fy_1= ',zwv16
            write(*,"(9999(G12.5,:,','))") 'fz_0= ',zww16,'df/f= ',ABS(c)/ABS(x_0*a+y_0*b)
            write(*,"(9999(G12.5,:,','))") '++++++++++++++++++++++++++++++++++++++++++++++'
       endif
            resx     = (b*( b*x_0-a*y_0)-a*c)/(a*a+b*b)
            resy     = (a*(-b*x_0+a*y_0)-b*c)/(a*a+b*b)

! Compute the cell height and velocity field
! following 4th order polynomial for velocity
! x^4 - x0 x^3 -y0 Fl x - Fl^2 
            aa(1) =   - resx*resx
            aa(2) =   resx*3.
            aa(3) =   0.
            aa(4) =   - REAL( un16(ji16-1,jj16,1),8)
            aa(5) =   1.

            call QuarticRoots(aa, z) ! Find velocity

            zz = 0.
            do i =1,4
               u = DBLE(z(i))
               if( ( AIMAG(z(i)) .EQ. 0)  .AND. (u .NE. 0.) )  then
                if (u * resx .GT. 0. ) then
                  if (u .GT. 0) then
                      u = max(-aa(4)*0.1,min(-aa(4)*20.,u))
                  else
                      u = min(-aa(4)*0.1,max(-aa(4)*20.,u))
                  endif
                  e3=resx/u
                  un16(ji16-1,jj16,1)    = u
                  ssh16x_R8(ji16,jj16)   = e3
                endif
               endif
            enddo

!y component
! Compute the cell height and velocity field
! following 4th order polynomial for velocity
! x^4 - x0 x^3 -y0 Fl x - Fl^2 
            aa(1) =   - resy*resy
            aa(2) =   resy*3.
            aa(3) =   0.
            aa(4) =   - REAL( vn16(ji16,jj16-1,1),8)
            aa(5) =   1.

            call QuarticRoots(aa, z) ! Find velocity

            zz = 0.
            do i =1,4
               v = DBLE(z(i))
               if( ( AIMAG(z(i)) .EQ. 0)  .AND. (v .NE. 0.) ) then
                if (v * resy .GT. 0. ) then
                  if (v .GT. 0) then
                      v = max(-aa(4)*0.1,min(-aa(4)*100.,v ))
                  else
                      v = min(-aa(4)*0.1,max(-aa(4)*100.,v ))
                  endif
                   e3=resy/v
                  vn16(ji16,jj16-1,1)    = v
                  ssh16y_R8(ji16,jj16)   = e3
                endif
               endif
            enddo

            e3u16(ji16,jj16,1)   = min(12.,max ( 0.1,REAL(ssh16x_R8(ji16+1,jj16),4) ))
            e3v16(ji16,jj16,1)   = min(12.,max ( 0.1,REAL(ssh16y_R8(ji16,jj16+1),4) ))
            e3t16(ji16,jj16,1)   = 0.25*(e3u16(ji16,jj16,1)+e3u16(ji16+1,jj16,1) + e3v16(ji16,jj16,1)+e3v16(ji16,jj16+1,1))
            ssh16(ji16,jj16)     = e3u16(ji16,jj16,1)

!           write(*,*) a*a+b*b, e3u16(ji16,jj16,1),e3v16(ji16,jj16,1)


          END DO  
        END DO  
!
!
!
      RETURN
      END
