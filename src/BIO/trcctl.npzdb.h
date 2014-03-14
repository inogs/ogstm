C $Id: trcctl.npzdb.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC Modifications:
CCC --------------
CCC    00-12 (E. Kestenare): 
CCC           assign a parameter to name individual tracers
CCC
c
#if defined key_trc_generic
      IF(lwp) THEN
          WRITE(numout,*) ' use npzdb tracer model '
          WRITE(numout,*) ' '
      ENDIF
#endif
