CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/parameter.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
C
CC---------------------------------------------------------------------
CCC
CCC                         PARAMETER
CCC                       *************
CCC
CCC  GLOBAL VERSION
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Include parameter file
CCC
CC   MODIFICATIONS :
CC   -------------
CC	original : 91 (Imbard, Levy, Madec)
CC
      IMPLICIT NONE
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
CC
CC Several domain sizes are parameterized
CC                                  global      domain   (jpiglo,jpjglo)
CC                                  computation domain   ( jpi  , jpj  )
CC                                  data        domain   (jpidta,jpjdta)
CC
CC * if we dont use massively parallel computer (parameters jpni=jpnj=1)
CC      so jpiglo=jpi and jpjglo=jpj
CC
CC Large domain matrix size
CC ------------------------
CC      jpni    : number of processors following i 
CC      jpnj    : number of processors following j
CC      jpnij   : < or = jpni x jpnj number of processors 
CC                (i.e. local domains)
CC                One can avoid land processors
CC      jpreci  : number of lines for overlap 
CC      jprecj  : number of lines for overlap 
CC
      INTEGER jpni,jpnj,jpnij,jpreci,jprecj
      PARAMETER(jpni=2,jpnj=2) 
      PARAMETER(jpnij=jpni*jpnj)
      PARAMETER(jpreci=1,jprecj=1)
CC
CC	jpiglo  : first  dimension of global domain --> i	 
CC	jpjglo  : second dimension of global domain --> j
CC	jpk	: number of vertical levels
CC
      INTEGER jpiglo,jpjglo,jpk
CCC 25/2/2004 Paolo  Dimension of Mediterranean Grid
#ifdef key_opal
#    ifdef key_antarctic
CCC      PARAMETER(jpiglo=92,jpjglo=15,jpk=31)
#      elif defined key_arctic
CCC      PARAMETER(jpiglo=32,jpjglo=24,jpk=31)
#      else
CCC      PARAMETER(jpiglo=92,jpjglo=76,jpk=31)
#    endif
#  else
#    ifdef key_antarctic
CCC      PARAMETER(jpiglo=182,jpjglo=50,jpk=31)
#      elif defined key_arctic
CCC      PARAMETER(jpiglo=142,jpjglo=53,jpk=31)
#      else
CCC      PARAMETER(jpiglo=182,jpjglo=149,jpk=31)
#    endif
#endif
CCC Paolo 8/4/2004 med parameters
#if defined key_grid831_mfs
CC
CC MFS831
CC
        PARAMETER(jpiglo=363,jpjglo=113,jpk=32)
#endif

CCC
CC
CC      jpkmod  : modulo jpk and jpnij
CC
      INTEGER jpkmod
      PARAMETER (jpkmod=1+(jpk-1)/jpnij)

CC Matrix size
CC -----------
CC	jpi     : first  dimension of grid --> i	 
CC	jpi     : first  dimension of grid --> i	 
CC
      INTEGER jpi,jpj
      PARAMETER(jpi=(jpiglo-2*jpreci + (jpni-1))/jpni + 2*jpreci)
      PARAMETER(jpj=(jpjglo-2*jprecj + (jpnj-1))/jpnj + 2*jprecj)
CC
CC
CC Other dimension parameters
CC --------------------------
CC      jpim1   :  jpi - 1
CC      jpjm1   :  jpj - 1
CC      jpkm1   :  jpk - 1
CC      jpij    :  jpi x jpj
CC
      INTEGER jpim1,jpjm1,jpkm1,jpij
      PARAMETER(jpim1=jpi-1,jpjm1=jpj-1,jpkm1=jpk-1,jpij=jpi*jpj)
CC
CC Original data size
CC ------------------
CC	jpidta	: first horizontal dimension > or = jpi 
CC	jpjdta	: second                     > or = jpj
CC	jpkdta	: number of levels           > or = jpk
CC
      INTEGER jpidta,jpjdta,jpkdta
#ifdef key_opal
CCC      PARAMETER(jpidta=92,jpjdta=76,jpkdta=31)
#  else
CCC      PARAMETER(jpidta=182,jpjdta=149,jpkdta=31)
#endif
CCC
CCC Paolo 8/4/2004 med parameters
#if defined key_grid831_mfs
CC
CC MFS831
CC
           PARAMETER(jpidta=363,jpjdta=113,jpkdta=32)
#endif

CC
CC
CC Domain characteristics
CC ----------------------
CC	jperio	: lateral cond. type for the global domain (4, 3, 2, 1
CC                or 0)  
CC	jpisl	: number of islands
CC	jpnisl	: maximum number of points per island 
CC
      INTEGER jperio,jpisl,jpnisl
CCC Paolo Mediterranean Conditions
#ifdef key_opal
#    ifdef key_antarctic
CCC      PARAMETER(jperio=1)
CCC      PARAMETER(jpisl =1,jpnisl=400)
#      elif defined key_arctic
CCC      PARAMETER(jperio=0)
CCC      PARAMETER(jpisl =0,jpnisl=200)
#      else
CCC      PARAMETER(jperio=1)
CCC      PARAMETER(jpisl =5,jpnisl=200)
#    endif
#else
#  ifdef key_antarctic
CCC      PARAMETER(jperio=1)
CCC      PARAMETER(jpisl =3,jpnisl=400)
#    elif defined key_arctic
CCC      PARAMETER(jperio=3)
CCC      PARAMETER(jpisl =3,jpnisl=400)
#    else
CCC      PARAMETER(jperio=4)
CCC      PARAMETER(jpisl =14,jpnisl=800)
#  endif
#endif
CCC Paolo 8/4/2004 med parameters
#if defined key_grid831_mfs
CC
CC MFS831
CC
         PARAMETER(jperio=0)
         PARAMETER(jpisl =9,jpnisl=800)
#endif



CC
CC
CC Multitasking parameter
CC ----------------------
CC	jpcpu	: maximum number of tasks
CC      jpbyt   : number of bytes for floating  4 or 8 
CC      jpbytda : number of bytes for floating in data files 4 or 8
CC
      INTEGER jpcpu,jpbyt,jpbi3e,jpbytda
#ifdef key_opal
      PARAMETER(jpcpu=8,jpbyt=8,jpbytda=4)
#else
      PARAMETER(jpcpu=8,jpbyt=8,jpbytda=4)
#endif
      PARAMETER(jpbi3e=4)
CC
CC Global grid characteristics
CC ---------------------------
#ifdef key_opal
CC      jpn     : number of latitude lines onto the north hemisphere
CC                                            (irregular grid)
cc      jpeq    : number of latitude lines onto the south hemisphere
cc                                            (  regular grid)
cc      jpl     : number of latitudes
cc      jpx     : maximum of grid points per latitude
CC
      INTEGER jpeq,jpn,jpl,jpx
      PARAMETER (jpeq=32,jpn=jpeq+4 )
      PARAMETER (jpl=jpeq+jpn,jpx=150)
#endif
CC
CC Offline reading parameter
CC -------------------------
CC
#include "parameter.offline.h"
CC
CC Passive tracers parameter
CC -------------------------
CC
#ifdef key_passivetrc
#    include "parameter.passivetrc.h"
#endif
CC
