C $Id: trcctl.nnpzddom.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC Modifications:
CCC --------------
CCC    00-12 (E. Kestenare): 
CCC           assign a parameter to name individual tracers
CCC
c
#if defined key_trc_nnpzddom
      IF(lwp) THEN
          WRITE(numout,*) ' use NNPZDDOM biological model '
          WRITE(numout,*) ' '
      ENDIF
c
c Check number of tracers
c -----------------------
      IF (jptra.ne.6) THEN 
          IF (lwp) THEN 
              WRITE (numout,*) ' ===>>>> : w a r n i n g '
              WRITE (numout,*) ' =======   ============= '
              WRITE (numout,*)
     $            ' STOP, change jptra to 6 in '
     $            ,'parameter.passivetrc.nnpzddom.h '  
          END IF 
          STOP 'TRCCTL'
      END IF 
c
c Check tracer names
c ------------------
      IF (    (ctrcnm(jpdet).ne.'DET').or.(ctrcnm(jpzoo).ne.'ZOO')
     $    .or.(ctrcnm(jpphy).ne.'PHY').or.(ctrcnm(jpno3).ne.'NO3')
     $    .or.(ctrcnm(jpnh4).ne.'NH4').or.(ctrcnm(jpdom).ne.'DOM') 
     $    .or.(ctrcnl(jpdet).ne.'DETRITUS')
     $    .or.(ctrcnl(jpzoo).ne.'ZOOPLANKTON')
     $    .or.(ctrcnl(jpphy).ne.'PHYTOPLANKTON')
     $    .or.(ctrcnl(jpno3).ne.'NITRATE')
     $    .or.(ctrcnl(jpnh4).ne.'AMMONIUM')
     $    .or.(ctrcnl(jpdom).ne.'DISSOLVED ORGANIC MATTER') ) THEN 
          ctrcnm(jpdet)='DET'
          ctrcnl(jpdet)='DETRITUS'
          ctrcnm(jpzoo)='ZOO'
          ctrcnl(jpzoo)='ZOOPLANKTON'
          ctrcnm(jpphy)='PHY'
          ctrcnl(jpphy)='PHYTOPLANKTON'
          ctrcnm(jpno3)='NO3'
          ctrcnl(jpno3)='NITRATE'
          ctrcnm(jpnh4)='NH4'
          ctrcnl(jpnh4)='AMMONIUM'
          ctrcnm(jpdom)='DOM'
          ctrcnl(jpdom)='DISSOLVED ORGANIC MATTER'          
          IF (lwp) THEN
              WRITE (numout,*) ' ===>>>> : w a r n i n g '
              WRITE (numout,*) ' =======   ============= '
              WRITE (numout,*) ' we force tracer names'
              DO jn=1,jptra
                WRITE(numout,*) ' tracer nb: ',jn,' name = ',ctrcnm(jn)
     $              ,ctrcnl(jn) 
              END DO
              WRITE(numout,*) ' '
          ENDIF 
      ENDIF 
c Check tracer units
      DO jn=1,jptra
        IF (ctrcun(jn).ne.'mmole/m3') THEN
            ctrcun(jn)='mmole/m3'
            IF (lwp) THEN
                WRITE (numout,*) ' ===>>>> : w a r n i n g '
                WRITE (numout,*) ' =======   ============= '
                WRITE (numout,*) ' we force tracer unit'
                WRITE(numout,*) ' tracer  ',ctrcnm(jn), 'UNIT= '
     $              ,ctrcun(jn)
                WRITE(numout,*) ' '
            ENDIF 
        ENDIF 
      END DO              
#endif
