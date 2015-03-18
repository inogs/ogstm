      SUBROUTINE trcdia(datemean, datefrom, dateend,FREQ_GROUP)

      USE myalloc
      USE myalloc_mpp
      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, datefrom, dateend
      INTEGER, INTENT(IN) :: FREQ_GROUP


       trcdiaparttime = MPI_WTIME() ! cronometer-start

!      writes ave files for tracer concentration
      CALL   trcdit(datemean, datefrom, dateend,FREQ_GROUP)
      CALL  diadump(datemean, datefrom, dateend,FREQ_GROUP)
      CALL fluxdump(datemean, datefrom, dateend)
  

       trcdiaparttime =   MPI_WTIME() - trcdiaparttime  ! cronometer-stop
       trcdiatottime  = trcdiatottime + trcdiaparttime

      END SUBROUTINE trcdia
