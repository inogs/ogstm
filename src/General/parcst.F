
      SUBROUTINE parcst
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE parcst
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Print the model parameters and set and print the constants
CC----------------------------------------------------------------------

       USE myalloc
        IMPLICIT NONE
CC ===================


      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' routine parcst'
      IF(lwp)WRITE(numout,*) ' **************'
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' initialization of opa model'
      IF(lwp)WRITE(numout,*) ' '

C 1. Parameters
C -------------

      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' parameter file'
      IF(lwp)WRITE(numout,*) ' **************'
      IF(lwp)WRITE(numout,*) ' '

C 1.1 Dimensions of domain

      IF(lwp)WRITE(numout,*) ' dimension of model'
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' jpni    : ', jpni
      IF(lwp)WRITE(numout,*) ' jpnj    : ', jpnj
      IF(lwp)WRITE(numout,*) ' jpnij   : ', jpnij
      IF(lwp)WRITE(numout,*) ' jpreci  : ', jpreci
      IF(lwp)WRITE(numout,*) ' jprecj  : ', jprecj
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' jpiglo  : ', jpiglo
      IF(lwp)WRITE(numout,*) ' jpjglo  : ', jpjglo
      IF(lwp)WRITE(numout,*) ' jpk     : ', jpk
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' jpi     : ', jpi
      IF(lwp)WRITE(numout,*) ' jpj     : ', jpj
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' jpidta  : ', jpidta
      IF(lwp)WRITE(numout,*) ' jpjdta  : ', jpjdta
      IF(lwp)WRITE(numout,*) ' jpkdta  : ', jpkdta
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' jpim1   : ', jpim1
      IF(lwp)WRITE(numout,*) ' jpjm1   : ', jpjm1
      IF(lwp)WRITE(numout,*) ' jpkm1   : ', jpkm1
      IF(lwp)WRITE(numout,*) ' jpij    : ', jpij
      IF(lwp)WRITE(numout,*) ' '

C 1.2 Option parameter : island, periodic,...

      IF(lwp)WRITE(numout,*) ' island,',
     $    ' lateral boundary conditions, s-coordinate'
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' jperio  : ', jperio
      IF(lwp)WRITE(numout,*) ' jpisl   : ', jpisl 
      IF(lwp)WRITE(numout,*) ' jpnisl  : ', jpnisl
      IF(lwp)WRITE(numout,*) ' '

C 1.3 Multitasking

      IF(lwp)WRITE(numout,*) ' multitasking'
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' jpcpu   : ', jpcpu

C 2. Define constants ( COMMON/comcst/)
C -------------------
C
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' constants'
      IF(lwp)WRITE(numout,*) ' *********'

C 2.1 Mathematical constants

      rpi = 2. * asin( 1. )
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' mathematical constant rpi = ', rpi

C 2.2 Astronomical constants

      rday   = 24. * 60. * 60.
      rsiyea = 365.25 * rday * 2. * rpi / 6.283076
      rsiday = rday / ( 1. + rday / rsiyea )
      omega  = 2. * rpi / rsiday 
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' day          rday   = ', rday,   ' s'
      IF(lwp)WRITE(numout,*) ' sideral year rsiyea = ', rsiyea, ' s'
      IF(lwp)WRITE(numout,*) ' sideral day  rsiday = ', rsiday, ' s'
      IF(lwp)WRITE(numout,*) ' omega        omega  = ', omega,  ' s-1'

C     year=raajj*days day=rjjhh*hours hour=rhhmm*minutes etc ...

      rjjhh = 24.
      rhhmm = 60.
      rmmss = 60.
      raamo = 12.
      rjjss = rjjhh * rhhmm * rmmss
      IF (calendarType .eq. 1) THEN
          raajj = 366.
          raass = raajj * rjjss
          rmoss = raass/raamo
      ELSE IF  (calendarType .eq. 0) THEN
          raajj = 365.
          raass = raajj * rjjss
          rmoss = raass/raamo
      ELSE 
          raajj = calendarType*raamo
          raass = raajj *rjjss
          rmoss = calendarType*rjjss
      ENDIF 
      IF(lwp) THEN
      WRITE(numout,*) ' '
      WRITE(numout,*) ' nb of days per year      raajj = ',raajj,' days'
      WRITE(numout,*) ' nb of months per year    raamo = ', raamo,' months'
      WRITE(numout,*) ' nb of hours per day      rjjhh = ', rjjhh,' hours'
      WRITE(numout,*) ' nb of minutes per hour   rhhmm = ', rhhmm, ' mn'
      WRITE(numout,*) ' nb of seconds per minute rmmss = ', rmmss, ' s'
      WRITE(numout,*) ' nb of seconds per year   raass = ', raass, ' s'
      WRITE(numout,*) ' nb of seconds per month  rmoss = ', rmoss, ' s'
      WRITE(numout,*) ' nb of seconds per day    rjjss = ', rjjss, ' s'
      ENDIF

C 2.3 Geoide

      ra = 6371229.
      g  = 9.80665
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' earth radius ra = ',ra,' m'
      IF(lwp)WRITE(numout,*) ' gravity      g  = ',g,' m s-2'

C 2.4 Thermodynamical constants

      rtt = 273.16
      rt0 = 273.15
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' fusion point rtt = ',rtt
      IF(lwp)WRITE(numout,*) ' zero celsius rt0 = ',rt0

C 2.5 physical constants

      rauw=1000.
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' density of pure water kg/m3 = ',rauw

C 2.6 Conversion terms

      rad = rpi / 180.
      rcs = rday * 1.e3
      IF(lwp)WRITE(numout,*) ' '
      IF(lwp)WRITE(numout,*) ' conversion: degre ==> radian rad = ', rad
      IF(lwp)WRITE(numout,*) ' conversion: mm/day ==> m/s   rcs = ', rcs

      END SUBROUTINE parcst
