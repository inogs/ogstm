      SUBROUTINE bn2 ( ktask, kt )
!---------------------------------------------------------------------
!
!                       ROUTINE bn2
!                     ***************
!
!  Purpose :
!  --------
!     Compute the local Brunt-Vaisala frequency
!
!   Method :
!   -------
!	neos = 0  : UNESCO sea water properties
!	   The brunt-vaisala frequency is computed using the polynomial
!	polynomial expression of McDougall (1987):
!		N^2 = g * beta * ( alpha/beta*dk[ tn ] - dk[ sn ] )/e3w
!
!	neos = 1  : linear equation of state (temperature only)
!		N^2 = g * ralpha * dk[ tn ]/e3w
!
!	neos = 2  : linear equation of state (temperature & salinity)
!		N^2 = g * (ralpha * dk[ tn ] - rbeta * dk[ sn ] ) / e3w
!
!	The use of potential density to compute N^2 introduces e r r o r
!      in the sign of N^2 at great depths. We recommand the use of
!      neos = 0, except for academical studies.
!
!      Macro-tasked on horizontal slab (jk-loop)
!
!      N.B. N^2 is set to zero at the first level (JK=1) in inidtr
!      and is never used at this level.
!
!   Output :
!   -------
!	common
!	      /comnow/ bn2n 	: now brunt-vaisala frequency
!
!   References :
!   -----------
!	McDougall, T. J., J. Phys. Oceanogr., 17, 1950-1964, 1987.


        USE myalloc
        USE stafun
        IMPLICIT NONE
!----------------------------------------------------------------------
! local declarations
! ==================
      INTEGER ktask, kt
      INTEGER ji, jj, jk

      REAL(8) zgde3w, zt, zs, zh, zalbet, zbeta
      REAL(8) fsalbt, fsbeta
      REAL(8) pft, pfh
!----------------------------------------------------------------------
! statement functions
! ===================


! Horizontal slab
! ===============
!
      DO 1000 jk = 1, jpkm1
!
!
! 1. first and last levels
! ------------------------
! bn^2=0. at jk=1 and jpk set in inidtr.F : no computation
!
        IF (jk.EQ.1) GOTO 1000
!
!
! 2. Interior points ( 2=< jk =< jpkm1 )
! ------------------
!
        IF ( neos.EQ.0 ) THEN
!
! ... UNESCO seawater equation of state
!
            DO jj = 1, jpj
              DO ji = 1, jpi
                zgde3w = g/fse3w(ji,jj,jk)
!   ... temperature, salinity anomalie and depth at w-level
                zt = 0.5*( tn(ji,jj,jk) + tn(ji,jj,jk-1) )
                zs = 0.5*( sn(ji,jj,jk) + sn(ji,jj,jk-1) ) - 35.0
                zh = fsdepw(ji,jj,jk)
!   ... ratio alpha/beta and beta
                zalbet = fsalbt( zt, zs, zh )
                zbeta  = fsbeta( zt, zs, zh )
!   ... N^2 = g/e3w * beta * ( alpha/beta * dk[tn] - dk[sn])
                bn2n(ji,jj,jk) = zgde3w * zbeta * tmask(ji,jj,jk)
     $                 * ( zalbet * ( tn(ji,jj,jk-1) - tn(ji,jj,jk) )
     $                            - ( sn(ji,jj,jk-1) - sn(ji,jj,jk) ) )
              END DO
            END DO
!
          ELSEIF ( neos.EQ.1 ) THEN
!
! ... First Linear density formulation (function of tempreature only)
!
            DO jj = 1, jpj
              DO ji = 1, jpi
                zgde3w = g/fse3w(ji,jj,jk)
                bn2n(ji,jj,jk) = zgde3w * ralpha * tmask(ji,jj,jk)
     $                           * ( tn(ji,jj,jk-1) - tn(ji,jj,jk) )
              END DO
            END DO
!
          ELSEIF ( neos.EQ.2 ) THEN
!
! ... Second linear density formulation (function of temp. and salinity)
!
            DO jj = 1, jpj
              DO ji = 1, jpi
                zgde3w = g/fse3w(ji,jj,jk)
                bn2n(ji,jj,jk) = zgde3w * tmask(ji,jj,jk)
     $               * (  ralpha * ( tn(ji,jj,jk-1) - tn(ji,jj,jk) )
     $                  - rbeta  * ( sn(ji,jj,jk-1) - sn(ji,jj,jk) )  )
              END DO
            END DO

        ENDIF

! End of slab
! ===========
!
 1000 continue

      END SUBROUTINE bn2

      ! ratio alpha/beta
! ================
!
!  fsalbt: ratio of thermal over saline expension coefficients
!       pft :  potential temperature in degrees celcius
!       pfs :  salinity anomaly (s-35) in psu
!       pfh :  depth in meters

      REAL(8) FUNCTION fsalbt( pft, pfs, pfh )
      IMPLICIT NONE
      REAL(8), INTENT(IN) ::  pft, pfs, pfh
      fsalbt =
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

      END FUNCTION fsalbt


! beta
! ====
!
!  fsbeta: saline expension coefficients
!       pft :  potential temperature in degrees celcius
!       pfs :  salinity anomaly (s-35) in psu
!       pfh :  depth in meters
!
      REAL(8) FUNCTION fsbeta( pft, pfs, pfh )
      IMPLICIT NONE
      REAL(8), INTENT(IN) ::  pft, pfs, pfh
      fsbeta =
     $  ( ( -0.415613e-09 * pft + 0.555579e-07 ) * pft
     *                          - 0.301985e-05 ) * pft
     $                          + 0.785567e-03
     $ +( 0.515032e-08 * pfs + 0.788212e-08 * pft - 0.356603e-06 ) * pfs
     $ +(  (   0.121551e-17 * pfh
     $       - 0.602281e-15 * pfs
     $       - 0.175379e-14 * pft + 0.176621e-12 ) * pfh
     $                            + 0.408195e-10   * pfs
     $    +( - 0.213127e-11 * pft + 0.192867e-09 ) * pft
     $                            - 0.121555e-07 ) * pfh

      END FUNCTION fsbeta
