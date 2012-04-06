!                       ROUTINE trcbbl
!                     ******************
! /*
!  Purpose :
!  --------
!     Compute the before tracer (t & s) trend associated with the
!    bottom boundary layer and add it to the general trend of tracer
!    equations. The bottom boundary layer is supposed to be a purely
!     diffusive bottom boundary layer.
!
!
!
!   Method :
!   -------
!    When the product grad( rho) * grad(h) < 0 (where grad is a
!    along bottom slope gradient) an additional lateral 2nd order
!    diffusion along the bottom slope is added to the general
!    tracer trend, otherwise the additional trend is set to 0.
!      Second order operator (laplacian type) with variable coefficient
!      computed as follow for temperature (idem on s):
!         difft = 1/(e1t*e2t*e3t) { di-1[ ahbt e2u*e3u/e1u di[ztb] ]
!                                 + dj-1[ ahbt e1v*e3v/e2v dj[ztb] ] }
!    where ztb is a 2D array: the bottom ocean te;perature and ahtb
!    is a time and space varying diffusive coefficient defined by:
!       ahbt = zahbp    if grad(rho).grad(h) < 0
!            = 0.       otherwise.
!    Note that grad(.) is the along bottom slope gradient. grad(rho)
!    is evaluated using the local density (i.e. referenced at the
!    local depth). Typical value of ahbt is 2000 m2/s (equivalent to
!    a downslope velocity of 20 cm/s if the condition for slope
!    convection is satified)
!
!      Add this before trend to the general trend (ta,sa) of the
!    botton ocean tracer point:
!              ta = ta + difft
!
!      'key_diatrends' defined, the trends are saved for diagnostics.
!
!      not macro-tasked
!
!   Input :
!   ------
!      argument
!              ktask           : task identificator
!              kt              : time step
!      common
!            /comcoo/          : mesh and scale factors
!            /comask/          : masks, bathymetry
!            /combef/          : previous fields (before)
!            /comaft/          : next fields (after)
!            /comtsk/          : multitasking
!
!   Output :
!   -------
!      common
!            /comaft/ ta, sa   : general tracer trend increased by the
!                                before horizontal diffusion trend
!            /comtrd/ ttrd,strd: now horizontal tracer diffusion trend
!                                (if 'key_diatrends' defined)
!
!   References :                 NO
!   -----------
!    Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
! */

      USE TIME_MANAGER
      INTEGER ktask, kt
      INTEGER ji, jj, jn
      INTEGER ik, iku, ikv
C
      REAL(8) fsalbt, pft, pfh
      REAL(8) zsign, zt, zs, zh, zalbet, zgdrho
      REAL(8) zbtr, zta, zsa
      REAL(8) zki(jpi,jpj), zkj(jpi,jpj)
      REAL(8) zkx(jpi,jpj), zky(jpi,jpj), ztrb(jpi,jpj)
      REAL(8) ztnb(jpi,jpj), zsnb(jpi,jpj), zdep(jpi,jpj)
      REAL(8) zahu(jpi,jpj), zahv(jpi,jpj)
!----------------------------------------------------------------------
! statement functions
! ===================

C ratio alpha/beta
C ================
C
C  fsalbt: ratio of thermal over saline expension coefficients
C       pft :  potential temperature in degrees celcius
C       pfs :  salinity anomaly (s-35) in psu
C       pfh :  depth in meters
C
      fsalbt( pft, pfs, pfh ) =
     $  ( ( ( -0.255019e-07 * pft + 0.298357e-05 ) * pft
     $                            - 0.203814e-03 ) * pft
     $                            + 0.170907e-01 ) * pft
     $                            + 0.665157e-01
     $ +(-0.678662e-05 * pfs - 0.846960e-04 * pft + 0.378110e-02 ) * pfs
     $ +  ( ( - 0.302285e-13 * pfh
     $        - 0.251520e-11 * pfs
     $        + 0.512857e-12 * pft * pft          ) * pfh
     $                             - 0.164759e-06   * pfs
     $     +(   0.791325e-08 * pft - 0.933746e-06 ) * pft
     $                             + 0.380374e-04 ) * pfh

      IF ( kt.EQ.TimeStepStart) THEN
          IF(lwp)WRITE(numout,*) '   bottom boundary layer file'
      ENDIF
C
C 0. 2D fields of bottom temperature and salinity, and bottom slope
C -----------------------------------------------------------------
C mbathy= number of w-level, minimum value=1 (cf dommsk.F)
C
      DO jj = 1, jpj
        DO ji = 1, jpi
C    ... index of the bottom ocean T-level
          ik = max( mbathy(ji,jj)-1, 1 )
C    ... masked now T and S at the ocean bottom 
          ztnb(ji,jj) = tn(ji,jj,ik) * tmask(ji,jj,1)
          zsnb(ji,jj) = sn(ji,jj,ik) * tmask(ji,jj,1)
C    ... masked before T and S at the ocean bottom 
C    ... depth of the ocean bottom T-level
          zdep(ji,jj) = fsdept(ji,jj,ik)
        END DO
      END DO
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
C
C
C 1. Criteria of additional bottom diffusivity: grad(rho).grad(h)<0
C --------------------------------------------
C Sign of the local density gradient along the i- and j-slopes
C multiplied by the slope of the ocean bottom
C
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
C   ... temperature, salinity anomalie and depth
          zt = 0.5 * ( ztnb(ji,jj) + ztnb(ji+1,jj) )
          zs = 0.5 * ( zsnb(ji,jj) + zsnb(ji+1,jj) ) - 35.0
          zh = 0.5 * ( zdep(ji,jj) + zdep(ji+1,jj) )
C   ... masked ratio alpha/beta
          zalbet = fsalbt( zt, zs, zh )*umask(ji,jj,1)
C   ... local density gradient along i-bathymetric slope
          zgdrho = zalbet*( ztnb(ji+1,jj) - ztnb(ji,jj) )
     $               -    ( zsnb(ji+1,jj) - zsnb(ji,jj) )
C   ... sign of local i-gradient of density multiplied by the i-slope
          zsign = sign( 0.5, -zgdrho * ( zdep(ji+1,jj) - zdep(ji,jj) ) )
          zki(ji,jj) = ( 0.5 - zsign ) * zahu(ji,jj)
        END DO
      END DO
C
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
C   ... temperature, salinity anomalie and depth
          zt = 0.5 * ( ztnb(ji,jj+1) + ztnb(ji,jj) )
          zs = 0.5 * ( zsnb(ji,jj+1) + zsnb(ji,jj) ) - 35.0
          zh = 0.5 * ( zdep(ji,jj+1) + zdep(ji,jj) )
C   ... masked ratio alpha/beta
          zalbet = fsalbt( zt, zs, zh )*vmask(ji,jj,1)
C   ... local density gradient along j-bathymetric slope
          zgdrho = zalbet*( ztnb(ji,jj+1) - ztnb(ji,jj) )
     $               -    ( zsnb(ji,jj+1) - zsnb(ji,jj) )
C   ... sign of local j-gradient of density multiplied by the j-slope
          zsign = sign( 0.5, -zgdrho * ( zdep(ji,jj+1) - zdep(ji,jj) ) )
          zkj(ji,jj) = ( 0.5 - zsign ) * zahv(ji,jj)
        END DO
      END DO
C
C
C 2. Additional second order diffusive trends
C -------------------------------------------
C
C ... first derivative (gradient)
      DO jn=1, jptra  

      DO jj = 1, jpj
        DO ji = 1, jpi
C    ... index of the bottom ocean T-level
          ik = max( mbathy(ji,jj)-1, 1 )
C    ... masked before tracer concentration at the ocean bottom
          ztrb(ji,jj) = trb(ji,jj,ik,jn) * tmask(ji,jj,1)
        END DO
      END DO

        DO jj = 1, jpjm1
          DO ji = 1, jpim1
            zkx(ji,jj) = zki(ji,jj)*( ztrb(ji+1,jj) - ztrb(ji,jj) )
            zky(ji,jj) = zkj(ji,jj)*( ztrb(ji,jj+1) - ztrb(ji,jj) )
        END DO
      END DO
C
#if defined key_orca_r2
C   Gibraltar enhancement of BBL
      DO jj = kfindj0(102), kfindj1(102)
        DO ji = kfindi0(139), kfindi1(140) 
          zkx(ji,jj)=4.*zkx(ji,jj)
          zky(ji,jj)=4.*zky(ji,jj)
        END DO
      END DO 
C   Red Sea enhancement of BBL
      DO jj = kfindj0(88), kfindj1(88)
        DO ji = kfindi0(161), kfindi1(162)
          zkx(ji,jj)=10.*zkx(ji,jj)
          zky(ji,jj)=10.*zky(ji,jj)
        END DO
      END DO 
#endif
#if defined key_orca_r4
C   Gibraltar enhancement of BBL
      DO jj = kfindj0(52), kfindj1(52)
        DO ji = kfindi0(70), kfindi1(71)
          zkx(ji,jj)=4.*zkx(ji,jj)
          zky(ji,jj)=4.*zky(ji,jj)
        END DO
      END DO 
#endif
C
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
C
