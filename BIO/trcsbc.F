      SUBROUTINE trcsbc ( ktask, kt )
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcsbc
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Compute the surface boundary condition trend of passive tracer 
CCC     tracers for default CASE : 2 tracers 1 for t & 2 for s
CCC     in the minimum configuration
CCC
CC   Method :
CC   -------
CC    * flux through the air-sea interface:
CC            - temperature : heat flux q (w/m2).
CC        tra(..1) = tra(..1) + q /(rau0 rcp e3t)  for k=1
CC            - salinity    : no salt flux
CC
CC    * concentration/dilution effect:
CC          - temperature : we assume that the temperature of both
CC       precipitations and runoffs is equal to the SST, thus there
CC       is no additional flux since in this case, the concentration
CC       dilution effect is balanced by the net heat flux associated
CC       to the freshwater exchange:
CC        (Tp P - Te E) + STT (P-E) = 0 when Tp=Te=SST
CC          - salinity    : evaporation, precipitation and runoff
CC       water has a zero salinity, thus
CC        sa = sa + emp * sn / e3t   for k=1
CC       where emp, the surface freshwater budget (evaporation minus
CC       precipitation minus runoff) given in kg/m2/s is divided
CC       by 1000 kg/m3 (density of plain water) to obtain m/s.
CC
CC    * vertical advective flux:
CC        The surface vertical advective flux applied in trazad.F
CC       is zero.
CC
CC      macro-tasked on vertical slab (jj-loop)
CC
CC   Input :
CC   ------
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC
CC   Output :
CC   -------
CC      common
CC            /comaft/ tra      : general temperature trend increased
CC                                by surface boundary condition trend

       USE myalloc
       USE stafun
       IMPLICIT NONE
CC local declarations
CC ==================
      INTEGER ktask, kt
#if defined key_passivetrc
      INTEGER ji, jj
      REAL(8) zta, zsa, zqfh, zsrau, zse3t

#    if !defined key_trc_npzd && !defined key_trc_nnpzddom
C default CASE : temperature and salinity as passive tracers


C Vertical slab
C =============

      zsrau=1. / rauw
      DO 1000 jj = ktask+1, jpjm1, ncpu


C 1. Concentration dillution effect on (t,s)
C -------------------------------------------

        DO ji = 2, jpim1
          zse3t=1./fse3t(ji,jj,1)

          zqfh = ro0cpr * qt(ji,jj)

C ... rigid lid
C   ... temperature : heat flux 
          zta = zqfh * zse3t
C   ... salinity :  concent./dilut. effect
          zsa = emp(ji,jj) * zsrau * trn(ji,jj,1,2) * zse3t

C   ... add the trend to the general tracer trend
          tra(ji,jj,1,1) = tra(ji,jj,1,1) + zta
          tra(ji,jj,1,2) = tra(ji,jj,1,2) + zsa

        END DO

C End of slab
C ===========

 1000 CONTINUE
#    endif

C Surface boundary conditions must be implemented here

#else

C     no passive tracer

#endif

      END SUBROUTINE trcsbc
