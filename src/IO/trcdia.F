      SUBROUTINE trcdia(datemean, datefrom, dateend,FREQ_GROUP)

      USE myalloc
      use mpi
      IMPLICIT NONE


      CHARACTER(LEN=17), INTENT(IN) :: datemean, datefrom, dateend
      INTEGER, INTENT(IN) :: FREQ_GROUP


       trcdiaparttime = MPI_WTIME() ! cronometer-start

!      writes ave files for tracer concentration
      CALL   trcdit(datemean, datefrom, dateend,FREQ_GROUP)
      CALL  diadump(datemean, datefrom, dateend,FREQ_GROUP)
      CALL fluxdump(datemean, datefrom, dateend,FREQ_GROUP)

      if (FREQ_GROUP.eq.1) then
         CALL Ed_dump(datemean, datefrom, dateend)
         CALL Es_dump(datemean, datefrom, dateend)
         CALL Eu_dump(datemean, datefrom, dateend)
      endif


       trcdiaparttime =   MPI_WTIME() - trcdiaparttime  ! cronometer-stop
       trcdiatottime  = trcdiatottime + trcdiaparttime

      END SUBROUTINE trcdia
