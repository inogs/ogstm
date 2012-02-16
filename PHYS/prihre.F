CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/prihre.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE prihre(ptab,ki,kj,kideb,kifin,kind,kjdeb,
     $    kjfin,kjnd,pscale,kumout)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE PRIHRE
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Print a real(8) field with the format 10e12.4 or 20f6.2
CCC
CC   METHOD :
CC   -------
CC      the format of print is selected with the dummy argument pscale
CC
CC   INPUT :
CC   -----
CC      argument ptab		: array for print
CC               ki,kj    	: array dimensions
CC               kideb,kifin	: first and last index for i
CC               kind     	: increment for i
CC               kjdeb,kjfin	: first and last index for j
CC               kjnd		: increment for j
CC               pscale		: if pscale=0 print ptab with the
CC                                format e13.5
CC                                else print pscale*ptab with the
CC                                format f6.2
CC               kumout		: unit for print
CC
CC
CC   OUTPUT :                     no
CC   ------
CC
CC   WORKSPACE : 
CC   ---------
CC	ji,jj,jn,isca,il1,il2
CC
CC   MODIFICATIONS:
CC   --------------
CC      original : 86-01 (P. Andrich)
CC      addition : 89-11 (C. Levy)
CC      addition : 92-02 (M. Imbard)
CC      addition : 92-06 (M. Imbard)
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
C
      IMPLICIT NONE
C
      INTEGER kideb,kifin,kind,kjdeb,kjfin,kjnd
      INTEGER ki,kj,kumout
      INTEGER ji,jj,jn
      INTEGER isca,il1,il2
C
      REAL(8) pscale
      REAL(8) ptab(ki,kj)
C
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
C
      isca=10
      IF(pscale.ne.0.) isca=20
C
      il1=kideb
C
      DO jn=1,(kifin-kideb+1)/(kind*isca) + 1
C
        IF(il1.gt.kifin)THEN
            RETURN
C           ------
        ENDIF
C    
        WRITE(kumout,9100)
C
        il2=il1+kind*(isca-1)
        IF(il2.gt.kifin) il2=kifin
        IF(pscale.eq.0.) THEN
            WRITE(kumout,9101) (ji,ji=il1,il2,kind)
            DO jj=kjfin,kjdeb,-kjnd
              WRITE(kumout,9102) jj,(ptab(ji,jj),ji=il1,il2,kind)
            END DO  
        ELSE
            WRITE(kumout,9103) (ji,ji=il1,il2,kind)
            DO jj=kjfin,kjdeb,-kjnd
              WRITE(kumout,9104) jj,(pscale*ptab(ji,jj),ji=il1,il2,kind)
            END DO  
        ENDIF
        il1=il1+kind*isca
C
      END DO  
C
C formats
C
 9100 FORMAT(/)
 9101 FORMAT(10i12,/)
 9102 FORMAT(1x,i3,10(1pe12.4))
 9103 FORMAT(3x,20i6,/)
 9104 FORMAT(1x,i3,1x,20f6.2)
C
C
      RETURN
      END
