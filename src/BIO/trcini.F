      SUBROUTINE trcini
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcini
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     initialization for passive tracer
CCC     for restart or not
CCC
CC   METHOD :
CC   -------
CC
CC   INPUT :
CC   -----
CC      argument		: no
CC      COMMON
CC
CC   OUTPUT :
CC   ------
CC      argument                : no
CC      COMMON
CC
CC   WORKSPACE :
CC   ---------
CC      ji jj jk
CC
CC
CC   EXTERNAL :
CC   --------
CC
CC   MODIFICATIONS:
CC   --------------
CC                 04/00 (O. Aumont, M.A. Foujols) HAMOCC3 and P3ZD
CC----------------------------------------------------------------------

       USE myalloc
       IMPLICIT NONE


#if defined key_passivetrc

      INTEGER ji,jj,jk,jn!,jl


#    if defined key_trc_npzd
#    include "trcini.npzd.h"
#    elif defined key_trc_bfm
!#    include "trcini.bfm.h"
      DO jj=1,jpj
        DO ji=1,jpi
          xze(ji,jj)=5.
        END DO
      END DO

      DO jk=1,jpk
        DO jj=1,jpj
          DO ji=1,jpi
            xpar(ji,jj,jk)=0.
          END DO
        END DO
      END DO

      call OPA_NO_BOXES(1,1,1,1,1)
      call Initialize()
#    elif defined key_trc_nnpzddom
#    include "trcini.nnpzddom.h"
#    elif defined key_trc_hamocc3
#    include "trcini.hamocc3.h"
#    else
C
C Specific initialization to do
C
#    endif
C
C 1. initialization of passives tracers field
C -------------------------------------------
      DO jn=1,jptra
        DO jk=1,jpk
          DO jj=1,jpj
            DO ji=1,jpi
              trn(ji,jj,jk,jn)=0.e0
              tra(ji,jj,jk,jn)=0.e0
            END DO
          END DO
        END DO
      END DO
#    if defined key_trc_diaadd
C initialization of output 2d and 3d arrays
      DO jn=1,jpdia2d
        DO jj=1,jpj
          DO ji=1,jpi
            trc2d(ji,jj,jn)=0.e0
          END DO
        END DO
      END DO
      DO jn=1,jpdia3d
        DO jk=1,jpk
          DO jj=1,jpj
            DO ji=1,jpi
              trc3d(ji,jj,jk,jn)=0.e0
            END DO
          END DO
        END DO
      END DO
#    endif
#    if defined key_trc_diabio
C initialization of biological trends
      DO jn=1,jpdiabio
        DO jk=1,jpk
          DO jj=1,jpj
            DO ji=1,jpi
              trbio(ji,jj,jk,jn)=0.e0
            END DO
          END DO
        END DO
      END DO
#    endif

#    if defined key_trc_diatrd
C initialization of tracer trends
      DO jl=1,jpdiatrc
        DO jn=1,jptra
          DO jk=1,jpk
            DO jj=1,jpj
              DO ji=1,jpi
                trtrd(ji,jj,jk,jn,jl)=0.e0
              END DO
            END DO
          END DO
        END DO
      END DO
#    endif      
      IF(lwp) THEN
          WRITE(numout,*) ' '
          WRITE(numout,*) ' trncini: generic initialisation done '
          WRITE(numout,*) ' '
      ENDIF
C
C
#else
C
C no passive trnacers
C
#endif

      END SUBROUTINE trcini
