      SUBROUTINE trcdia(datemean, datefrom, dateend)

!      CALLED by : step

      USE myalloc
      IMPLICIT NONE

#include <mpif.h>

      CHARACTER(LEN=17), INTENT(IN) :: datemean, datefrom, dateend



      ! local
      INTEGER kt, kindic



       trcdiaparttime = MPI_WTIME() !F79 cronometer-start

!      writes ave files for tracer concentration
      CALL   trcdit(datemean, datefrom, dateend)
      CALL  diadump(datemean, datefrom, dateend)
      CALL fluxdump(datemean, datefrom, dateend)
  
       kindic = 0 ! serve a niente...
! ************************************************ FINISCE QUI


!      prendere da trcdit la parte jn=1, fisica, e metterla in
!     una routine per la fisica
!     una routine per altre variabili diagnostiche

       kt = 1 ! solo per compilare le funzioni che seguono

#if defined key_trc_diatrd
      CALL trcdid(kt,kindic) ! outputs for dynamical trends
#endif

#if defined key_trc_diaadd
      CALL trcdii(kt,kindic) ! outputs for additional arrays
#endif

#if defined key_trc_diabio
      CALL trcdib(kt,kindic) ! outputs for biological trends
#endif



       trcdiaparttime =   MPI_WTIME() - trcdiaparttime  ! F79 cronometer-stop
       trcdiatottime  = trcdiatottime + trcdiaparttime

      END SUBROUTINE trcdia
