! $Header: /cvsroot/opatm-bfm/opa_model/OPA/allcray.F,v 1.2 2009-09-11 09:20:55 cvsogs01 Exp $

!  Cray subroutines or functions used by OPA model and possibly
!  not found on other platforms.

!  check their existence

!  ssum
!  sdot
!  wheneq
!  orders
!  isamin
!  isamax
!  isrchne
!  isrchfge
!  isrchfne
!  rnlskip
!  ishell

!---------------------------------------------------------
        FUNCTION SSUM(I,X,J)
        IMPLICIT NONE
        INTEGER I,J,N
        REAL(8) X(I),SSUM
        SSUM=0.
        DO 1 N=1,I
        SSUM=SSUM+X(1+(N-1)*J)
  1     CONTINUE
        RETURN
        END
!---------------------------------------------------------
        FUNCTION SDOT(I,X,J,Y,K)
        IMPLICIT NONE
        INTEGER I,J,K,N
        REAL(8) X(I),Y(I),SDOT
        SDOT=0.
        DO 1 N=1,I
        SDOT=SDOT+X(1+(N-1)*J)*Y(1+(N-1)*K)
  1     CONTINUE
        RETURN
        END
!---------------------------------------------------------
        SUBROUTINE WHENEQ(I,X,J,T,IND,NN,dimenarry)
        IMPLICIT NONE
        INTEGER N,I,J,NN,dimenarry
        INTEGER IND(dimenarry)
        REAL(8) X(dimenarry),T

        NN=0
        DO 1 N=1,I
        IF(X(1+(N-1)*J).EQ.T)THEN
        NN=NN+1
        IND(NN)=1+(N-1)*J
        ENDIF
  1     CONTINUE
        RETURN
        END
!---------------------------------------------------------
        SUBROUTINE ORDERS(IMODE,IWORK,IDATA,INDIX,N,IRE,IK,IRA)
        IMPLICIT NONE
        INTEGER I,IMODE,IWORK,N,IRE,IK,IRA,ITEMP
        INTEGER IDATA(IRE,N),INDIX(N)

           DO 10 I=1,N
             INDIX(I)=I
 10         CONTINUE
             DO 20 I=1,N-1
            IF(IDATA(IRE,INDIX(I+1)).LT.IDATA(IRE,INDIX(I))) THEN
              ITEMP=INDIX(I)
              INDIX(I)=INDIX(I+1)
              INDIX(I+1)=ITEMP
            ENDIF
 20       CONTINUE
            DO 30 I=N,2,-1
            IF(IDATA(IRE,INDIX(I)).LT.IDATA(IRE,INDIX(I-1))) THEN
              ITEMP=INDIX(I)
              INDIX(I)=INDIX(I-1)
              INDIX(I-1)=ITEMP
            ENDIF
 30         CONTINUE
        RETURN
        END
!---------------------------------------------------------
        FUNCTION ISAMIN(I,X,J)
        IMPLICIT NONE
        INTEGER I,J,N
        REAL(8) ISAMIN, XMIN, X(I)
        ISAMIN=0
        XMIN=1.E+20
        DO 1 N=1,I
        IF(ABS(X(N)).LT.XMIN)THEN
        XMIN=X(N)
        ISAMIN=N
        ENDIF
  1     CONTINUE
        RETURN
        END
!-----------------------------------------------------------
        FUNCTION ISAMAX(I,X,J)
        IMPLICIT NONE
        INTEGER I,J,N
        REAL(8) ISAMAX, XMAX, X(I)
        ISAMAX=0
        XMAX=-1E+20
        DO 1 N=1,I
        IF(ABS(X(N)).GT.XMAX)THEN
        XMAX=X(N)
        ISAMAX=N
        ENDIF
  1     CONTINUE
        RETURN
        END
!---------------------------------------------------------
        FUNCTION ISRCHNE(K,X,I,B)
        IMPLICIT NONE
        INTEGER N,K,I,ISRCHNE
        REAL(8) B, X(1+(K-1)*I)
        DO 1 N=1,K
        IF(X(1+(N-1)*I).NE.B)THEN
         ISRCHNE=N
        RETURN
        ELSE
          ISRCHNE=N+1
        ENDIF
 1      CONTINUE
        RETURN
        END
!---------------------------------------------------------
        FUNCTION ISRCHFGE(K,X,I,B)
        IMPLICIT NONE
        INTEGER N,K,I,ISRCHFGE
        REAL(8) B, X(1+(K-1)*I)
        DO 1 N=1,K
        IF(X(1+(N-1)*I).GE.B)THEN
         ISRCHFGE=N
        RETURN
        ELSE
          ISRCHFGE=N+1
        ENDIF
 1      CONTINUE
        RETURN
        END
!---------------------------------------------------------
        FUNCTION ISRCHFLE(K,X,I,B)
        IMPLICIT NONE
        INTEGER N,K,I,ISRCHFLE
        REAL(8) B, X(1+(K-1)*I)
        DO 1 N=1,K
        IF(X(1+(N-1)*I).LE.B)THEN
         ISRCHFLE=N
        RETURN
        ELSE
          ISRCHFLE=N+1
        ENDIF
 1      CONTINUE
        RETURN
        END
!---------------------------------------------------------
        SUBROUTINE RNLSKIP(KFLAG)
        INTEGER KFLAG
        RETURN
        END
!---------------------------------------------------------
        SUBROUTINE ISHELL (CL)
        CHARACTER*80 CL
        RETURN
        END
!---------------------------------------------------------

