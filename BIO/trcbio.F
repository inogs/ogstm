
      SUBROUTINE trcbio(ktask,kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcbio
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     compute the now trend due to biogeochemical processes
CCC     and add it to the general trend of passive tracers equations.
CCC
CCC	Three options:
CCC  	  Default option  : no biological trend
CCC       If 'key_trc_npzd' : NPZD model
CCC       IF 'key_trc_nnpzddom' : NNPZDDOM model
CCC
CC   METHOD :
CC   -------
CC      each now biological flux is calculated  in FUNCTION of now
CC      concentrations of tracers.
CC      depending on the tracer, these fluxes are sources or sinks.
CC      the total of the sources and sinks for each tracer
CC      is added to the general trend.
CC
CC		tra = tra + zf...tra - zftra...
CC                             |         |
CC                             |         |
CC                          source      sink
CC
CC
CC      IF 'key_trc_diabio' key is activated, the biogeochemical
CC	trends for passive tracers are saved for futher diagnostics.
CC
CC      multitasked on vertical slab (jj-loop)
CC
CCC   MODIFICATIONS:
CC   --------------
CC       original : 95-02 (M. Levy, NPZD model)
CC                  99-07 (M. Levy, NNPZDDOM model)

      USE STAFUN

#if defined key_passivetrc 
#  if defined key_trc_npzd 
#    include "trcbio.npzd.h"
#  elif defined key_trc_bfm
#    include "trcbio.bfm.h"
#  elif defined key_trc_npzdb
#    include "trcbio.npzdb.h"
#  elif defined key_trc_generic
#    include "trcbio.generic.h"
#  else
C
C    no biological model
C
#  endif
#endif

      END SUBROUTINE trcbio
