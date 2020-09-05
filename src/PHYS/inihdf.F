      SUBROUTINE inihdf
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE inihdf
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Initializations of the horizontal coefficients
!!!
!!   Method :
!!   -------
!!
!!    Eddy diffusivity and eddy induced velocity cefficients:
!!         default option   : constant coef. aht0 (namelist)
!!        'key_trahdfcoef1d': depth dependent coef. defined in 
!!                            in 'inihdf.tra.coef3d.h' file
!!
!!    Direction of lateral diffusion (tracers and/or momentum)
!!
!!    Biharmonic operator:
!!
!!   ---------
!!
!!   Reference :
!!   ----------
!!      Madec, G. and M. Imbard, 1996, A global ocean mesh to overcome
!!      the North Pole singularity, Climate Dynamics, 12, 381-388.
       USE myalloc
        IMPLICIT NONE
!! ==================
      INTEGER  jk
!!

      double precision zkah, zahr, za00, za01, zahf, zahs
      double precision zahtf, zahts


!! 1. Lateral eddy diffusivity and eddy induced velocity coefficients
!! ==================================================================
!!     user: modify the 'inihdf.tra.coef_d.h' include files following
!!           your desiderata

!!
!! 1D coefficients : 
!! ----------------
!!!
!!!             inihdf.tra.coef1d.h
!!!           ***********************
!!!
!!!   defined key :  'key_trahdfcoef1d'
!!!   ============
!!!
!!      1D eddy diffusivity coefficients ( depth )
!!      --------------------------------
!!
!!       biharmonic operator    : ahtt = defined at T-level
!!                                ahtu,ahtv,ahtw never used
!!
!!       harmonic operator      : ahtt never used
!!          -1-  iso-model level: ahtu = ahtv defined at T-level
!!                                ahtw never used
!!          -2-  isopycnal or   : ahtu = ahtv defined at T-level
!!               geopotential     ahtw defined at w-level
!!
!!       eddy induced velocity
!!         always harmonic      : aeiu = aeiv defined at T-level
!!				  aeiw defined at w-level
!!
!!!---------------------------------------------------------------------
!!!  OPA8, LODYC (1997)
!!!---------------------------------------------------------------------
!!
      IF(lwp)WRITE(numout,*)

      IF(lwp)WRITE(numout,*) ' inihdf: 1D eddy diffusivity coefficient '
      IF(lwp)WRITE(numout,*) ' ======  --'
      IF(lwp)WRITE(numout,*)
!!    
!! ... initialization of the profile
!!
!!   ... ahts, ahtf: surface and bottom values
       hdf_initialized = .false.
       adv_initialized = .false.
       zahts = aht0
       zahtf = aht0/2.
!!! 
!!   ... zkah, zahr: depth of the inflection pt and width of inflection
          zkah =  -150.
          zahr =   150.
!!!
!!   ... computation coefficients
      za00 = tanh( (-gdept(1,1,1 )-zkah) / zahr )
      za01 = tanh( (-gdept(jpk,1,1)-zkah) / zahr )
      zahf = ( zahts-zahtf ) / ( za00 - za01 )
      zahs = zahts - zahf * za00
!!
#if defined key_trahdfbilap
!!
!! biharmonic operator : 
!! ==================== 
!!
!! ... set ahtt at T-level
      DO jk = 1, jpk
        ahtt(jk) = aht0
      END DO
!!
!! ... control print
      IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) '         aht profile at T-level : '
          WRITE(numout,*)
          WRITE(numout,9200)
          DO jk = 1, jpk
          END DO
      ENDIF
 9200 FORMAT(' level    aht          depth t-level ' )
 9210 FORMAT(i6,2f12.4)
!!
#  else
!!
!! harmonic operator :
!! ================== 
!!
!! ... set ahtu = ahtv at T-level, and ahtw at w-level
      DO jk = 1, jpk
        ahtu(jk) = zahs + zahf * tanh( (-gdept(jk,1,1)-zkah) / zahr )
        ahtv(jk) = ahtu(jk)
        ahtw(jk) = zahs + zahf * tanh( (-gdepw(jk)-zkah) / zahr )
      END DO
!!
!! ... control print
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '         aht profile at T-level : '
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,9200)
      DO jk = 1, jpk
        IF(lwp)WRITE(numout,9210) jk, ahtu(jk), gdept(jk,1,1)
      END DO
 9200 FORMAT(' level    aht          depth t-level ' )
 9210 FORMAT(i6,2f12.4)
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) '         aht profile at W-level : '
      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,9220)
      DO jk = 1, jpk
        IF(lwp)WRITE(numout,9210) jk, ahtw(jk), gdepw(jk)
      END DO
 9220 FORMAT('  jk      aht          depth w-level ' )


#endif




      END SUBROUTINE inihdf
