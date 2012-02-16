CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/vgmrea.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
      SUBROUTINE vgmrea(kenr)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE DYNREA
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     READ dynamics fiels from OPA7 standard output
CCC    
CCC 
CC
CC   METHOD :
CC   -------
CC
CC      READ the kenr records of DATA and store in udta(...,2)
CC      vdta(...,2), wdta(...,2), avtdta(...,2), sdta(...,2), tdta(...,2)
CC
CC   INPUT :
CC   -----
CC      argument
CC              kt              : time step
CC		kindic		: indicator of abnormal termination
CC      common
CC
CC   OUTPUT :
CC   ------
CC      file
CC            /COMDYN/          : dynamic fields 
CC
CC
CC   WORKSPACE :
CC   ---------
CC      local
CC              ji,jj,jk, jrec, jfic, junit
CC
CC   EXTERNAL :
CC   --------
CC
CC   REFERENCES :                 no
CC   ----------
CC
CC   MODIFICATIONS:
CC   --------------
CC      original  : 
CC      additions : 
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================

c#include "parameter.h"
c#include "common.h"
       USE myalloc

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
CC Local declarations
CC ==================

      INTEGER kenr

      INTEGER ji, jj, jk, jrec, junit



CC----------------------------------------------------------------------
CC Statement functions
CC ===================

#include "stafun.h"

CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------



#if defined key_trahdfeiv && defined key_off_degrad
        junit=numvgm

C go to the begining of FILE

        REWIND (junit)


C next 12 month for the second year, the good one

        DO jrec=1,kenr-1
          READ (junit)          !
          READ (junit)          !
          READ (junit)          !
        END DO

C READ good fields

cast        READ (junit) ugm        
cast    READ (junit) vgm       
cast    READ (junit) wgm      
        call read3s(junit,ugm,jpi,jpj,jpk)
        call read3s(junit,vgm,jpi,jpj,jpk)
        call read3s(junit,wgm,jpi,jpj,jpk)


        IF (lwp) THEN
            WRITE (numout,*) ' vgm fields '
        END IF 


#endif
      RETURN
      END
