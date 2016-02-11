      SUBROUTINE eos ()
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE eos
CCC                     ***************
CCC
CCC  Purpose :
CCC  --------
CCC    Compute the in situ density (ratio rho/rau0) and the potential
CCC    volumic mass (Kg/m3) from potential temperature and salinity
CCC    fields after the array swap using an equation of state defined
CCC    through the namelist parameter neos.
CCC    C a u t i o n : here time filtering technique apply to double
CCC    the stability limit associated with internal gravity waves
CCC    (brown & campana 1978).
CC
CC   Method :
CC   -------
CC    default option : use of (tn,sn) after the array swap, eos. is
CC       called at the end of step routine, after the call of dynhpg.
CC
CC    'key_hpgimplicit' defined : use of (ta,sa) after the array swap,
CC       which contain an average over three time levels (before, now
CC       and after) and reset (ta,sa) to zero. eos is called in step
CC       before the call of dynhpg.
CC
CC    neos = 0 : Jackett and McDougall (1994) equation of state.
CC         the in situ density is computed directly as a function of
CC         potential temperature relative to the surface (the ogstm t
CC         variable), salt and pressure (assuming no pressure variation
CC         along geopotential surfaces, i.e. the pressure p in decibars
CC         is approximated by the depth in meters.
CC              rdn(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
CC              rhop(t,s)  = rho(t,s,0)
CC         with pressure                      p        decibars
CC              potential temperature         t        deg celsius
CC              salinity                      s        psu
CC              reference volumic mass        rau0     kg/m**3
CC              in situ volumic mass          rho      kg/m**3
CC              in situ density anomalie      rdn      no units
CC
CC         Check value: rho = 1059.8204 kg/m**3 for p=10000 dbar,
CC          t = 40 deg celcius, s=40 psu
CC
CC      neos = 1 : linear equation of state function of temperature only
CC              rdn(t) = ( rho(t) - rau0 ) / rau0 = 0.028 - ralpha * t
CC              rhop(t,s)  = rho(t,s)
CC
CC      neos = 2 : linear equation of state function of temperature and
CC           salinity
CC              rdn(t,s) = ( rho(t,s) - rau0 ) / rau0 
CC                 = rbeta * s - ralpha * tn - 1.
CC              rhop(t,s)  = rho(t,s)
CC
CC      macro-tasked on horizontal slab (jk-loop)
CC
CC      Note that no boundary condition problem occurs in this routine
CC      as (tn,sn) or (ta,sa) are defined over the whole domain.
CC
CC
CC   Input :
CC   ------
CC      argument
CC      common
CC            /comcoo/          : scale factors
CC            /comnow/        : present fields (now)
CC
CC   Output :
CC   -------
CC    common
CC            /comnow/ rdn()    : now in situ density (no units)
CC               rhop()    : now potential volumic mass (Kg/m3)
CC
CC   References :
CC   -----------
CC      Jackett, D.R., and T.J. McDougall. J. Atmos. Ocean. Tech., 1994
CC    Brown, J. A. and K. A. Campana. Mon. Weather Rev., 1978
CC
CC   Modifications :
CC   --------------
CC      original :  89-03 (o. Marti)
CC      additions : 94-08 (G. Madec)
CC                : 96-01 (G. Madec) statement function for e3
CC          : 97-07 (G. Madec) introduction of neos, OPA8.1
CC          : 97-07 (G. Madec) density instead of volumic mass
CC                : 99-02 (G. Madec, N. Grima) 'key_hdfimplicit'
CC----------------------------------------------------------------------

       USE myalloc
       USE myalloc_mpp
       IMPLICIT NONE

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER ji, jj, jk
      REAL(8) zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw
      REAL(8) zb, zd, zc, zaw, za, zb1, za1, zkw, zk0

      INTEGER :: mytid, ntids! omp variables

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      density_partTime = MPI_WTIME()

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = 1
      mytid = 0
#endif


      DO jk = 1, jpkm1, ntids

        IF ( neos.EQ.0 ) THEN
!$omp parallel default(none) private(mytid,jj,ji,
!$omp&                 zt, zs, zh, zsr, zr1, zr2, zr3, zr4, zrhop, ze, zbw,
!$omp&                   zb, zd, zc, zaw, za, zb1, za1, zkw, zk0)
!$omp&                   shared(jpi,jpj,jpk,jk,tn,sn,rau0,rhopn,rdn,rho,tmask,gdept)
#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif


C 1. Jackett and McDougall (1994) formulation
C -------------------------------------------
        if (jk+mytid.le.jpk) then
        DO jj = 1, jpj
          DO ji = 1, jpi

C   ... now potential temperature and salinity
            zt = tn(ji,jj,jk+mytid)
            zs = sn(ji,jj,jk+mytid)
C   ... depth
            zh = gdept(jk+mytid)
C   ... square root salinity
            zsr= sqrt( abs( zs ) )
C   ... compute volumic mass pure water at atm pressure
            zr1= ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt-9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594
C   ... seawater volumic mass atm pressure
            zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 ) *zt+7.6438e-5 ) *zt -4.0899e-3 ) *zt+0.824493
            zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3
            zr4= 4.8314e-4

C   ... potential volumic mass (reference to the surface)
            zrhop= ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1

C   ... save potential volumic mass
            rhopn(ji,jj,jk+mytid) = zrhop

C   ... add the compression terms
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

C   ... masked in situ density anomaly
            rdn(ji,jj,jk+mytid) = ( zrhop
     $          / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) 
     $          - rau0 ) / rau0 * tmask(ji,jj,jk+mytid)
C   ... masked in situ density
            rho(ji,jj,jk+mytid)=zrhop / (  1.0 - zh / ( zk0 - zh * ( za - zh * zb ) )  ) * tmask(ji,jj,jk+mytid)

          END DO
        END DO
      endif
!$omp end parallel
          ELSEIF( neos.EQ.1 ) THEN


C 2. First Linear density formulation (function of tempreature only)
C -----------------------------------

        DO jj = 1, jpj
          DO ji = 1, jpi

C   ... now potential temperature and salinity
            zt = tn(ji,jj,jk)
C   ... density and potential volumic mass
            rdn(ji,jj,jk) = ( 0.028 - ralpha * zt ) * tmask(ji,jj,jk)
            rhopn(ji,jj,jk) = ( rau0 * rdn(ji,jj,jk) + rau0 )* tmask(ji,jj,jk)
          END DO
        END DO

          ELSEIF( neos.EQ.2 ) THEN


C 3. Second linear density formulation (function of temp. and salinity)
C ------------------------------------

        DO jj = 1, jpj
          DO ji = 1, jpi

C   ... now potential temperature and salinity
            zt = tn(ji,jj,jk)
            zs = sn(ji,jj,jk)

C   ... density and potential volumic mass
            rdn(ji,jj,jk) = (   rbeta  * zs - ralpha * zt - 1. )* tmask(ji,jj,jk)
            rhopn(ji,jj,jk) = ( rau0 * rdn(ji,jj,jk) + rau0 )   * tmask(ji,jj,jk)
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

      END DO

      density_partTime =    MPI_WTIME()  - density_partTime
      density_TotTime  = density_TotTime + density_partTime


      END SUBROUTINE eos
