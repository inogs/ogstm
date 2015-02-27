       MODULE myalloc_mpp

       IMPLICIT NONE

       public

CC 
#include <mpif.h>
CC 
CC    size : number of process
CC    rank : process number  [ 0 - size-1 ]
CC
      INTEGER  :: size,rank
      INTEGER  :: mpi_pack_size =18


        END MODULE
