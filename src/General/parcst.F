      SUBROUTINE parcst
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE parcst
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Print the model parameters-
CC----------------------------------------------------------------------

       USE myalloc
        IMPLICIT NONE


      if (lwp) then
      WRITE(numout,*) ' '
      WRITE(numout,*) ' routine parcst'
      WRITE(numout,*) ' **************'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' initialization of ogstm model'
      WRITE(numout,*) ' '

C 1. Parameters
C -------------

      WRITE(numout,*) ' '
      WRITE(numout,*) ' parameter file'
      WRITE(numout,*) ' **************'
      WRITE(numout,*) ' '

C 1.1 Dimensions of domain

      WRITE(numout,*) ' dimension of model'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpni    : ', jpni
      WRITE(numout,*) ' jpnj    : ', jpnj
      WRITE(numout,*) ' jpnij   : ', jpnij
      WRITE(numout,*) ' jpreci  : ', jpreci
      WRITE(numout,*) ' jprecj  : ', jprecj
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpiglo  : ', jpiglo
      WRITE(numout,*) ' jpjglo  : ', jpjglo
      WRITE(numout,*) ' jpk     : ', jpk
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpi     : ', jpi
      WRITE(numout,*) ' jpj     : ', jpj
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpidta  : ', jpidta
      WRITE(numout,*) ' jpjdta  : ', jpjdta
      WRITE(numout,*) ' jpkdta  : ', jpkdta
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jpim1   : ', jpim1
      WRITE(numout,*) ' jpjm1   : ', jpjm1
      WRITE(numout,*) ' jpkm1   : ', jpkm1
      WRITE(numout,*) ' jpij    : ', jpij
      WRITE(numout,*) ' '

C 1.2 Option parameter : island, periodic,...

      WRITE(numout,*) ' island,',' lateral boundary conditions, s-coordinate'
      WRITE(numout,*) ' '
      WRITE(numout,*) ' jperio  : ', jperio
      WRITE(numout,*) ' '


      endif

      END SUBROUTINE parcst
