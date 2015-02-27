      SUBROUTINE trcdia(datemean, datefrom, dateend,FREQ_GROUP)

      USE myalloc
      IMPLICIT NONE

#include <mpif.h>

      CHARACTER(LEN=17), INTENT(IN) :: datemean, datefrom, dateend
      INTEGER, INTENT(IN) :: FREQ_GROUP


       trcdiaparttime = MPI_WTIME() !F79 cronometer-start

!      writes ave files for tracer concentration
      CALL   trcdit(datemean, datefrom, dateend,FREQ_GROUP)
      CALL  diadump(datemean, datefrom, dateend,FREQ_GROUP)
      CALL fluxdump(datemean, datefrom, dateend)
  

       trcdiaparttime =   MPI_WTIME() - trcdiaparttime  ! F79 cronometer-stop
       trcdiatottime  = trcdiatottime + trcdiaparttime

      END SUBROUTINE trcdia
