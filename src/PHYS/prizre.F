      SUBROUTINE prizre(ptab,ki,kj,kk,kjcut,kideb,kifin,kind,kkdeb,
     $    kkfin,kknd,pscale,kumout)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE PRIZRE
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Print a vertical slab from a tridimentional real(8) field
CCC
CC   METHOD :
CC   -------
CC      the format of print is selected with the dummy argument pscale
CC
CC   INPUT :
CC   -----
CC      argument ptab           : array for print
CC               ki,kj,kk       : array dimensions
CC               kjcut          : index j for the vertical slab
CC               kideb,kifin    : first and last index for i
CC               kind           : increment for i
CC               kkdeb,kkfin    : first and last index for k
CC               kknd           : increment for k
CC               pscale         : IF pscale=0 print ptab with the
CC                                format 10e12.4
CC                                else print pscale*ptab with the
CC                                format 20f6.2
CC               kumout         : unit for print
CC
CC   OUTPUT :                     no
CC   ------
CC
CC   WORKSPACE : ji,jn,jk,isca,il1,il2 
CC   ---------
CC    ji,jn,jk,isca,il1,il2
CC
CC   MODIFICATIONS:
CC   --------------
CC      original : 86-01 (o. Marti)
CC      addition : 92-02 (M. Imbard)
CC      addition : 92-06 (M. Imbard)
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
C
      IMPLICIT NONE
      INTEGER ki,kj,kk,kjcut,kideb,kifin,kind,kkdeb,kkfin,kknd
      INTEGER kumout
      INTEGER ji,jn,jk
      INTEGER isca,il1,il2
C
      REAL(8) pscale
      REAL(8) ptab(ki,kj,kk)

C
      isca=10
      IF(pscale.ne.0.) isca=20
C
      IF (kind.eq.0) kind=1
      IF (kknd.eq.0) kknd=1
C
      il1=kideb
C
      DO jn=1,(kifin-kideb+1)/(kind*isca) + 1
C
        IF(il1.gt.kifin) RETURN
        WRITE(kumout,9100)
        il2=il1+kind*(isca-1)
        IF(il2.gt.kifin) il2=kifin
C
        IF(pscale.eq.0.) THEN
            WRITE(kumout,9101) (ji,ji=il1,il2,kind)
            DO jk=kkdeb,kkfin,kknd
              WRITE (kumout,9102) jk,(ptab(ji,kjcut,jk),ji=il1,il2,kind)
            END DO  
        ELSE
            WRITE (kumout,9103) (ji,ji=il1,il2,kind)
            DO jk=kkdeb,kkfin,kknd
              WRITE(kumout,9104)jk,(pscale*ptab(ji,kjcut,jk),
     $            ji=il1,il2,kind)
            END DO  
        ENDIF
C
        il1=il1+kind*isca
      END DO      
C
 9100 FORMAT(/)
 9101 FORMAT(10i12,/)
 9102 FORMAT(1x,i3,10(1pe12.4))
 9103 FORMAT(3x,20i6,/)
 9104 FORMAT(1x,i3,1x,20f6.1)
C
C
      RETURN
      END
