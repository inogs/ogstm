CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/prihin.F,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $

      SUBROUTINE prihin(ktab,ki,kj,kideb,kifin,kind,kjdeb,
     $    kjfin,kjnd,kscale,kumout)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE PRIHIN
CCC                     ******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Print an integer field
CCC
CC   METHOD :
CC   -------
CC      the format of print is selected with the dummy argument kscale
CC
CC   INPUT :
CC   -----
CC      argument itab            : integer array for print
CC               ki,kj           : array dimensions
CC               kideb,kifin     : first and last index for i
CC               kind            : increment for i
CC               kjdeb,kjfin     : first and last index for j
CC               kjnd            : increment for j
CC               kscale          : IF kscale=0
CC                                 or IF > 5   print itab with format i8
CC                                 IF kscale=1 print itab with format i1
CC                                 IF kscale=2 print itab with format i2
CC                                 IF kscale=3 print itab with format i3
CC                                 IF kscale=4 print itab with format i4
CC                                 IF kscale=5 print itab with format i5
CC               kumout          : unit for print
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
CC      original : 90-04 (0. Marti)
CC      addition : 92-02 (M. Imbard)
CC      addition : 92-06 (M. Imbard)
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================
C
      IMPLICIT NONE
C
      INTEGER kideb,kifin,kind,kjdeb,kjfin,kjnd
      INTEGER ki,kj,kscale,kumout
      INTEGER ktab(ki,kj)
C
      INTEGER ji,jj,jn
      INTEGER isca,il1,il2
C
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
C
      isca=10
      IF(kscale.eq.0) isca=10
      IF(kscale.eq.1) isca=100
      IF(kscale.eq.2) isca=60
      IF(kscale.eq.3) isca=40
      IF(kscale.eq.4) isca=30
      IF(kscale.eq.5) isca=20
C
      IF (kind.eq.0) kind=1
      IF (kjnd.eq.0) kjnd=1
C
      il1=kideb
C
      DO jn=1,(kifin-kideb+1)/(kind*isca) + 1
C
        IF(il1.gt.kifin) RETURN
        WRITE(kumout,'(/)')
        il2=il1+kind*(isca-1)
        IF(il2.gt.kifin) il2=kifin
C
        IF (kscale.eq.1) THEN
            WRITE(kumout,'(4x,i14," to ",1i4," each ",1i4,/)') il1,
     $          il2,kind
            DO jj=kjfin,kjdeb,-kjnd
              WRITE (kumout,'(1x,i3,100i1)') jj,(ktab(ji,jj),
     $            ji=il1,il2,kind)
            END DO  
        ELSE IF (kscale.eq.2) THEN
            WRITE(kumout,'(4x,i14," to ",1i4," each ",1i4,/)')il1,
     $          il2,kind
            DO jj=kjfin,kjdeb,-kjnd
              WRITE (kumout,'(1x,i3,60i2)') jj,(ktab(ji,jj),
     $            ji=il1,il2,kind)
            END DO  
        ELSE IF (kscale.eq.3) THEN
            WRITE(kumout,'(4x,i14," to ",1i4," each ",1i4,/)')il1,
     $          il2,kind
            DO jj=kjfin,kjdeb,-kjnd
              WRITE (kumout,'(1x,i3,40i3)') jj,(ktab(ji,jj),
     $            ji=il1,il2,kind)
            END DO  
        ELSE IF(kscale.eq.4) THEN
            WRITE(kumout,'(4x,30i4,/)') (ji,ji=il1,il2,kind)
            DO jj=kjfin,kjdeb,-kjnd
              WRITE (kumout,'(1x,i3,30i4)') jj,(ktab(ji,jj),
     $            ji=il1,il2,kind)
            END DO  
        ELSE IF(kscale.eq.5) THEN
            WRITE(kumout,'(4x,20i5,/)') (ji,ji=il1,il2,kind)
            DO jj=kjfin,kjdeb,-kjnd
              WRITE (kumout,'(1x,i3,20i5)') jj,(ktab(ji,jj),
     $            ji=il1,il2,kind)
            END DO  
        ELSE
            WRITE(kumout,'(4x,10i8,/)') (ji,ji=il1,il2,kind)
            DO jj=kjfin,kjdeb,-kjnd
              WRITE (kumout,'(1x,i3,10i8)') jj,(ktab(ji,jj),
     $            ji=il1,il2,kind)
            END DO  
        ENDIF
C
        il1=il1+kind*isca
      END DO  
C
C 
      RETURN
      END
