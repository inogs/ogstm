
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
        logical mlwp
CC ===================
        mlwp=.false.

      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' routine parcst'
      IF(mlwp)WRITE(numout,*) ' **************'
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' initialization of opa model'
      IF(mlwp)WRITE(numout,*) ' '

C 1. Parameters
C -------------

      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' parameter file'
      IF(mlwp)WRITE(numout,*) ' **************'
      IF(mlwp)WRITE(numout,*) ' '

C 1.1 Dimensions of domain

      IF(mlwp)WRITE(numout,*) ' dimension of model'
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' jpni    : ', jpni
      IF(mlwp)WRITE(numout,*) ' jpnj    : ', jpnj
      IF(mlwp)WRITE(numout,*) ' jpnij   : ', jpnij
      IF(mlwp)WRITE(numout,*) ' jpreci  : ', jpreci
      IF(mlwp)WRITE(numout,*) ' jprecj  : ', jprecj
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' jpiglo  : ', jpiglo
      IF(mlwp)WRITE(numout,*) ' jpjglo  : ', jpjglo
      IF(mlwp)WRITE(numout,*) ' jpk     : ', jpk
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' jpi     : ', jpi
      IF(mlwp)WRITE(numout,*) ' jpj     : ', jpj
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' jpidta  : ', jpidta
      IF(mlwp)WRITE(numout,*) ' jpjdta  : ', jpjdta
      IF(mlwp)WRITE(numout,*) ' jpkdta  : ', jpkdta
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' jpim1   : ', jpim1
      IF(mlwp)WRITE(numout,*) ' jpjm1   : ', jpjm1
      IF(mlwp)WRITE(numout,*) ' jpkm1   : ', jpkm1
      IF(mlwp)WRITE(numout,*) ' jpij    : ', jpij
      IF(mlwp)WRITE(numout,*) ' '

C 1.2 Option parameter : island, periodic,...

      IF(mlwp)WRITE(numout,*) ' island,',
     $    ' lateral boundary conditions, s-coordinate'
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' jperio  : ', jperio
      IF(mlwp)WRITE(numout,*) ' jpisl   : ', jpisl
      IF(mlwp)WRITE(numout,*) ' jpnisl  : ', jpnisl
      IF(mlwp)WRITE(numout,*) ' '

C 1.3 Multitasking

      IF(mlwp)WRITE(numout,*) ' multitasking'
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' jpcpu   : ', jpcpu

C 2. Define constants ( COMMON/comcst/)
C -------------------
C
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' constants'
      IF(mlwp)WRITE(numout,*) ' *********'

C 2.1 Mathematical constants

      rpi = 2. * asin( 1. )
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' mathematical constant rpi = ', rpi

C 2.2 Astronomical constants

      rday   = 24. * 60. * 60.
      rsiyea = 365.25 * rday * 2. * rpi / 6.283076
      rsiday = rday / ( 1. + rday / rsiyea )
      omega  = 2. * rpi / rsiday 
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' day          rday   = ', rday,   ' s'
      IF(mlwp)WRITE(numout,*) ' sideral year rsiyea = ', rsiyea, ' s'
      IF(mlwp)WRITE(numout,*) ' sideral day  rsiday = ', rsiday, ' s'
      IF(mlwp)WRITE(numout,*) ' omega        omega  = ', omega,  ' s-1'

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
      IF(mlwp) THEN
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
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' earth radius ra = ',ra,' m'
      IF(mlwp)WRITE(numout,*) ' gravity      g  = ',g,' m s-2'

C 2.4 Thermodynamical constants

      rtt = 273.16
      rt0 = 273.15
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' fusion point rtt = ',rtt
      IF(mlwp)WRITE(numout,*) ' zero celsius rt0 = ',rt0

C 2.5 physical constants

      rauw=1000.
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' density of pure water kg/m3 = ',rauw

C 2.6 Conversion terms

      rad = rpi / 180.
      rcs = rday * 1.e3
      IF(mlwp)WRITE(numout,*) ' '
      IF(mlwp)WRITE(numout,*) ' conversion: degre ==> radian rad = ', rad
      IF(mlwp)WRITE(numout,*) ' conversion: mm/day ==> m/s   rcs = ', rcs

      END SUBROUTINE parcst
