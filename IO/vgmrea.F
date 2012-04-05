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


       USE myalloc
       IMPLICIT NONE
       INTEGER kenr
#if defined key_trahdfeiv && defined key_off_degrad

      INTEGER jrec, junit

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
