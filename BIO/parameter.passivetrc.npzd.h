C $Id: parameter.passivetrc.npzd.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC---------------------------------------------------------------------
CCC
CCC                         PARAMETER passivetrc.npzd
CCC                       *****************************
CCC
CCC  purpose :
CCC  ---------
CCC     INCLUDE PARAMETER FILE for passive tracer NPZD model
CCC
CCC  modifications:  
CCC  --------------
CCC     00-12 (E. Kestenare): 
CCC            assign a parameter to name individual tracers
CCC     01-02 (E. Kestenare): 
CCC            introduce jpno3 instead of jpnut
CC
CC       jptra  : number of tracers
CC
      INTEGER jptra
      PARAMETER (jptra = 4)
CC
CC productive layer depth
CC ----------------------
CC       jpkb   : first vertical layers where biology is active
CC       jpkbm1 : jpkb - 1
CC
      INTEGER jpkb,jpkbm1
      PARAMETER (jpkb = 12,jpkbm1 = 11)
CC
CC number of biological trends
CC ---------------------------
CC
      INTEGER jpdiabio
      PARAMETER (jpdiabio = 9)
CC
CC    NOW ASSIGN A PARAMETER TO NAME INDIVIDUAL TRACERS
CC
CC    JPDET : detritus (mmoleN/m3)
CC    JPZOO : zooplancton concentration (mmoleN/m3)
CC    JPPHY : phytoplancton concentration (mmoleN/m3)
CC    JPNO3 : nutrate concentration (mmoleN/m3)
CC
      INTEGER jpdet,jpzoo,jpphy,jpno3
      PARAMETER (jpdet=1,jpzoo=2,jpphy=3,jpno3=4)

