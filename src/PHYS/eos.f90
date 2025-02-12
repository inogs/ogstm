      SUBROUTINE eos ()
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE eos
!!!                     ***************
!!!
!!!  Purpose :
!!!  --------
!!!    Compute the in situ density (ratio rho/rau0) and the potential
!!!    volumic mass (Kg/m3) from potential temperature and salinity
!!!    fields after the array swap using an equation of state defined
!!!    through the namelist parameter neos.
!!!    !! a u t i o n : here time filtering technique apply to double
!!!    the stability limit associated with internal gravity waves
!!!    (brown & campana 1978).
!!
!!   Method :
!!   -------
!!    default option : use of (tn,sn) after the array swap, eos. is
!!       called at the end of step routine, after the call of dynhpg.
!!
!!    'key_hpgimplicit' defined : use of (ta,sa) after the array swap,
!!       which contain an average over three time levels (before, now
!!       and after) and reset (ta,sa) to zero. eos is called in step
!!       before the call of dynhpg.
!!
!!    neos = 0 : Jackett and McDougall (1994) equation of state.
!!         the in situ density is computed directly as a function of
!!         potential temperature relative to the surface (the ogstm t
!!         variable), salt and pressure (assuming no pressure variation
!!         along geopotential surfaces, i.e. the pressure p in decibars
!!         is approximated by the depth in meters.
!!              rdn(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
!!              rhop(t,s)  = rho(t,s,0)
!!         with pressure                      p        decibars
!!              potential temperature         t        deg celsius
!!              salinity                      s        psu
!!              reference volumic mass        rau0     kg/m**3
!!              in situ volumic mass          rho      kg/m**3
!!              in situ density anomalie      rdn      no units
!!
!!         Check value: rho = 1059.8204 kg/m**3 for p=10000 dbar,
!!          t = 40 deg celcius, s=40 psu
!!
!!      neos = 1 : linear equation of state function of temperature only
!!              rdn(t) = ( rho(t) - rau0 ) / rau0 = 0.028 - ralpha * t
!!              rhop(t,s)  = rho(t,s)
!!
!!      neos = 2 : linear equation of state function of temperature and
!!           salinity
!!              rdn(t,s) = ( rho(t,s) - rau0 ) / rau0 
!!                 = rbeta * s - ralpha * tn - 1.
!!              rhop(t,s)  = rho(t,s)
!!
!!      macro-tasked on horizontal slab (jk-loop)
!!
!!      Note that no boundary condition problem occurs in this routine
!!      as (tn,sn) or (ta,sa) are defined over the whole domain.
!!
!!
!!   Input :
!!   ------
!!      argument
!!      common
!!            /comcoo/          : scale factors
!!            /comnow/        : present fields (now)
!!
!!   Output :
!!   -------
!!    common
!!            /comnow/ rdn()    : now in situ density (no units)
!!               rhop()    : now potential volumic mass (Kg/m3)
!!
!!   References :
!!   -----------
!!      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
!!    Brown, J. A. and K. A. Campana. Mon. Weather Rev., 1978
!!
!!   Modifications :
!!   --------------
!!      original :  89-03 (o. Marti)
!!      additions : 94-08 (G. Madec)
!!                : 96-01 (G. Madec) statement function for e3
!!          : 97-07 (G. Madec) introduction of neos, OPA8.1
!!          : 97-07 (G. Madec) density instead of volumic mass
!!                : 99-02 (G. Madec, N. Grima) 'key_hdfimplicit'
!!----------------------------------------------------------------------

       USE myalloc
       USE mpi
       ! epascolo USE myalloc_mpp
       IMPLICIT NONE

!!----------------------------------------------------------------------
!! local declarations
!! ==================
      INTEGER jk,jj,ji
      double precision :: zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw
      double precision :: zb, zd, zc, zaw, za, zb1, za1, zkw, zk0


      density_partTime = MPI_WTIME()

        IF ( neos.EQ.0 ) THEN
!!!$omp parallel default(none) private(mytid,jj,ji,
!!!$omp&                 zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw,
!!!$omp&                   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
!!!$omp&                   shared(jpk,jpj,jpi,jk,tn,sn,rau0,rhopn,rdn,rho,tmask,gdept)

           ! eos 0.214075s

           !$acc parallel loop gang vector collapse(3) default(present)
           DO ji = 1, jpi
              DO jj = 1, jpj
                 DO jk = 1, jpkm1

!!   ... now potential temperature and salinity
            zt = tn(jk,jj,ji)
            zs = sn(jk,jj,ji)
!!   ... depth

#ifdef gdept1d
            zh = gdept(jk)
#else
            zh = gdept(jk,jj,ji)
#endif

!!   ... square root salinity
            zsr= sqrt( abs( zs ) )
!!   ... compute volumic mass pure water at atm pressure
            zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt-9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
!!   ... seawater volumic mass atm pressure
            zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt -4.0899e-3 ) *zt+0.824493
            zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
            zr4= 4.8314e-4

!!   ... potential volumic mass (reference to the surface)
            zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1

!!   ... save potential volumic mass
            rhopn(jk,jj,ji) = zrhop

!!   ... add the compression terms
            ze = ( -3.508914e-8*zt-1.248266e-8 ) *zt-2.595994e-6
            zbw= (  1.296821e-6*zt-5.782165e-9 ) *zt+1.045941e-4
            zb = zbw + ze * zs

            zd = -2.042967e-2
            zc =   (-7.267926e-5*zt+2.598241e-3 ) *zt+0.1571896
            zaw= ( ( 5.939910e-6*zt+2.512549e-3 ) *zt-0.1028859 ) *zt -4.721788
            za = ( zd*zsr + zc ) *zs + zaw

            zb1=   (-0.1909078*zt+7.390729 ) *zt-55.87545
            za1= ( ( 2.326469e-3*zt+1.553190)*zt-65.00517 ) *zt+1044.077
            zkw= ( ( (-1.361629e-4*zt-1.852732e-2 ) *zt-30.41638 ) *zt +2098.925 ) *zt+190925.6
            zk0= ( zb1*zsr + za1 )*zs + zkw

!!   ... masked in situ density anomaly
            rdn(jk,jj,ji) = ( zrhop &
     &          / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) & 
     &          - rau0 ) / rau0 * tmask(jk,jj,ji)
!!   ... masked in situ density
            rho(jk,jj,ji)=zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) * tmask(jk,jj,ji)

          END DO
        END DO
      ENDDO
      !$acc end parallel loop
!!!$omp end parallel

     ELSEIF( neos.EQ.1 ) THEN

#ifdef _OPENACC
        print *, "eos: code path not tested on GPU (neos=1)"
        stop
#endif

!! 2. First Linear density formulation (function of tempreature only)
!! -----------------------------------

                  DO ji = 1, jpi
            DO jj = 1, jpj
        DO jk = 1, jpkm1
!!   ... now potential temperature and salinity
            zt = tn(jk,jj,ji)
!!   ... density and potential volumic mass
            rdn(jk,jj,ji) = ( 0.028 - ralpha * zt ) * tmask(jk,jj,ji)
            rhopn(jk,jj,ji) = ( rau0 * rdn(jk,jj,ji) + rau0 )* tmask(jk,jj,ji)
          ENDDO
        ENDDO
      ENDDO

  ELSEIF( neos.EQ.2 ) THEN

#ifdef _OPENACC
        print *, "eos: code path not tested on GPU (neos=2)"
        stop
#endif

!! 3. Second linear density formulation (function of temp. and salinity)
!! ------------------------------------

          DO ji = 1, jpi
        DO jj = 1, jpj
      DO jk = 1, jpkm1
!!   ... now potential temperature and salinity
            zt = tn(jk,jj,ji)
            zs = sn(jk,jj,ji)

!!   ... density and potential volumic mass
            rdn(jk,jj,ji) = (   rbeta  * zs - ralpha * zt - 1. )* tmask(jk,jj,ji)
            rhopn(jk,jj,ji) = ( rau0 * rdn(jk,jj,ji) + rau0 )   * tmask(jk,jj,ji)
          END DO
        END DO
      END DO

  ELSE

          IF(lwp) THEN
          WRITE(numout,*) ' E R R O R in neos flag '
          WRITE(numout,*) ' =========    ===='
          WRITE(numout,*) ' '
          WRITE(numout,*) ' we stop'
          WRITE(numout,*) ' '
          ENDIF
          STOP 'eos.f'

      ENDIF

      density_partTime =    MPI_WTIME()  - density_partTime
      density_TotTime  = density_TotTime + density_partTime


      END SUBROUTINE eos
