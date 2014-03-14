      SUBROUTINE ctlopn ( knum, cdfile, cdstat, cdform, cdacce,
     $                    klengh, klu, kused, kout, ldwp, krequ )
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE ctlopn
CCC                     ******************
CCC
CCC  Purpose :
CCC  --------
CCC     Open file and check if required file is available.
CCC
CC   Method :
CC   -------
CC      Fortan open
CC
CC   Input :
CC   ------
CC      argument:       knum    : logical unit to open
CC                      cdfile  : file name to open
CC                      cdstat  : disposition specifier
CC                      cdform  : formatting specifier
CC                      cdacce  : access specifier
CC                      klengh  : record length
CC                      klu     : number of logical units used 
CC                                (incremented by 1 if file open)
CC                      kused   : array of logical units used (one more
CC                                element if file open)
CC                      kout    : number of logical units for write
CC                      ldwp    : boolean term for print
CC                      krequ   : =1 file required (stop if not exist)
CC                                =0 file not required (create the file
CC                                   if does not exist)
CC
CC   Output :                     no
CC   -------
CC
CC   Modifications :
CC   --------------
CC       original :  95-12 (G. Madec)
CC----------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER knum, klu, krequ, kout, klengh, iost
      INTEGER kused(1,100)
C
      CHARACTER(LEN=21) cdacce, cdform, cdstat
      CHARACTER cdfile*(*)
C
      LOGICAL ldwp
C
C 1. Required file
C ----------------
C
      IF ( krequ.EQ.1 ) THEN
C
          iost=0
          IF ( cdacce(1:6).EQ.'DIRECT' )  THEN
              OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,
     $              STATUS=cdstat, RECL=klengh, ERR=100, IOSTAT=iost )
            ELSE
              OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,
     $              STATUS=cdstat, ERR=100, IOSTAT=iost)
          ENDIF
          IF ( iost.EQ.0 ) THEN
              IF(ldwp) THEN
                  WRITE(kout,*) '     file   : ', cdfile,' open ok'
                  WRITE(kout,*) '     unit   = ', knum
                  WRITE(kout,*) '     status = ', cdstat
                  WRITE(kout,*) '     form   = ', cdform
                  WRITE(kout,*) '     access = ', cdacce
                  WRITE(kout,*)
              ENDIF
              klu = klu + 1
              kused(1,klu) = knum
          ENDIF
 100      CONTINUE
          IF ( iost.NE.0 ) THEN
              IF(ldwp) THEN
                  WRITE(kout,*)
                  WRITE(kout,*) ' ===>>>> : bad opening file: ', cdfile
                  WRITE(kout,*) ' =======   ===  '
                  WRITE(kout,*) '           unit   = ', knum
                  WRITE(kout,*) '           status = ', cdstat
                  WRITE(kout,*) '           form   = ', cdform
                  WRITE(kout,*) '           access = ', cdacce
                  WRITE(kout,*) '           iostat = ', iost
                  WRITE(kout,*) '           we stop. verify the file '
                  WRITE(kout,*)
              ENDIF
              STOP 'ctlopn'
          ENDIF
C
C
C 2. Not required, file create if not exist
C -----------------------------------------
C
        ELSEIF ( krequ.EQ.0 ) THEN
C
          iost = 0
          IF ( cdacce(1:6).EQ.'DIRECT' ) THEN
              OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,
     $              STATUS=cdstat, RECL=klengh, ERR=200, IOSTAT=iost )
            ELSE
              OPEN( UNIT=knum, FILE=cdfile, FORM=cdform, ACCESS=cdacce,
     $              STATUS=cdstat, ERR=200, IOSTAT=iost )
          ENDIF
          IF (iost.eq.0) THEN
              IF(ldwp) THEN
                  WRITE(kout,*) '     file   : ', cdfile,' open ok'
                  WRITE(kout,*) '     unit   = ', knum
                  WRITE(kout,*) '     status = ', cdstat
                  WRITE(kout,*) '     form   = ', cdform
                  WRITE(kout,*) '     access = ', cdacce
                  WRITE(kout,*)
              ENDIF
              klu = klu + 1
              kused(1,klu) = knum
          ENDIF
 200      CONTINUE
          IF ( iost.NE.0 ) THEN
              iost = 0
              IF(ldwp) THEN
                  WRITE(kout,*)
                  WRITE(kout,*) '     ===>>>> : file ', cdfile,
     $                ' does not exist: it is created'
                  WRITE(kout,*) ' =======   ===  '
              ENDIF
              IF ( cdacce(1:6).EQ.'DIRECT' ) THEN
                  OPEN( UNIT=knum, FILE=cdfile, FORM=cdform,
     $                  ACCESS=cdacce, STATUS=cdstat,
     $                  RECL=klengh, ERR=210, IOSTAT=iost )
                ELSE
                  OPEN( UNIT=knum, FILE=cdfile, FORM=cdform,
     $                  ACCESS=cdacce, STATUS=cdstat, ERR=210,
     $                  IOSTAT=iost )
              ENDIF
              IF(ldwp) THEN
                  WRITE(kout,*) '     file   : ', cdfile,' open ok'
                  WRITE(kout,*) '     unit   = ', knum
                  WRITE(kout,*) '     status = ', cdstat
                  WRITE(kout,*) '     form   = ', cdform
                  WRITE(kout,*) '     access = ', cdacce
                  WRITE(kout,*)
              ENDIF
              klu = klu + 1
              kused(1,klu) = knum
 210          CONTINUE
              IF ( iost.NE.0 ) THEN
                  IF(ldwp) THEN
                      WRITE(kout,*) ' logical unit ',knum,' iostat = ',
     $                    iost
                      WRITE(kout,*) ' we stop. verify the file ', cdfile
                      WRITE(kout,*)
                  ENDIF
                  stop001
              ENDIF
          ENDIF
C
        ELSE
C
          IF(ldwp) THEN
              WRITE(kout,*)
              WRITE(kout,*) ' ctlopn : invalid option, krequ = ', krequ
              WRITE(kout,*) ' ******   call for file ', cdfile
              WRITE(kout,*)
          ENDIF
          STOP 'ctlopn'
      ENDIF
C
C
      RETURN
      END
