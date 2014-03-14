CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/modele.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
      PROGRAM modele
CCC---------------------------------------------------------------------
CCC
CCC                       PROGRAM MODELE
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Send opa model, adjoint model, linear tangent
CCC
C
      write(*,*) "Starting ..."
      CALL myopa()
C
      STOP
      END    
