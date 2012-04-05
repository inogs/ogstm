CCC $Id: trcadv.ppm.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC---------------------------------------------------------------------------
CCC
CCC                         trcadv.ppm.h
CCC                     ******************
CCC
CC   defined key : 'key_trc_ppm'
CC   ===========
CC  
CC   PURPOSE :
CC   ---------
CC      Compute the now trend due to total advection of tracers
CC      and add it to the general trend of tracer equations.
CC      PPM(Piecewise Parabolic Method) scheme is used.
CC 
CC 
CC   METHOD :
CC   -------
CC      this ROUTINE compute not exactly the advection but the
CC      transport term, i.e.  div(u*tra).
CC
CC 
CC   REFERENCES :                
CC   ----------                  
CC
CC   MODIFICATIONS:
CC   --------------
CC       original :  07-00 (A.Estublier)
CC       
CC----------------------------------------------------------------------
C
      INTEGER ji,jj,jk,jn
      REAL(8) ztrax(jpi,jpj,jpk),zatrax(jpi,jpj,jpk)
      REAL(8) ztray(jpi,jpj,jpk),zatray(jpi,jpj,jpk)
      REAL(8) ztraz(jpi,jpj,jpk),zatraz(jpi,jpj,jpk),zakz(jpi,jpj,jpk) 
      REAL(8) zkx(jpi,jpj,jpk),zky(jpi,jpj,jpk),zkz(jpi,jpj,jpk) 
      REAL(8) zkx1(jpi,jpj,jpk)
      REAL(8) ztlx(jpi,jpj,jpk),ztrx(jpi,jpj,jpk)
      REAL(8) ztly(jpi,jpj,jpk),ztry(jpi,jpj,jpk)
      REAL(8) ztlz(jpi,jpj,jpk),ztrz(jpi,jpj,jpk)
      REAL(8) za6x(jpi,jpj,jpk), za6y(jpi,jpj,jpk), za6z(jpi,jpj,jpk)
      REAL(8) zdelta(jpi,jpj,jpk),zeta(jpi,jpj,jpk),zetat(jpi,jpj,jpk)    
      REAL(8) zigma,zeu,zev,zew,zbt,zstep,zign,zign1,zmask(jpi,jpj,jpk)
      REAL(8) zbig,za,zb,zc,zd,zk1,zk2,zeps
C
CC----------------------------------------------------------------------
CC statement functions
CC ===================

!                                   #include "stafun.h"

CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (01/00)
CCC---------------------------------------------------------------------
C

    zbig = 1.e+40
    zk1  = 8.e-4
    zk2  = -300.
    zeps = 1.e-3

C
C tracer loop parallelized (macrotasking)
C =======================================
C
      DO 1000 jn = ktask,jptra,ncpu        

     

C
C Initialization
C --------------
C

    zstep  = rdt*ndttrc

        ztrax  = 0
        zatrax = 0
        zkx    = 0
    ztlx   = 0
    ztrx   = 0
    za6x   = 0

        ztray  = 0
        zatray = 0
        zky    = 0
    ztly   = 0
    ztry   = 0
    za6y   = 0

        ztraz  = 0
        zatraz = 0
        zkz    = 0
    ztlz   = 0
    ztrz   = 0
    za6z   = 0

    zdelta = 0
    zeta   = 0
    zetat  = 0


C
C 1. Slopes computation
C ---------------------
C

C
C In the i and j direction
C
        DO jk=1,jpkm1
          DO jj = 1,jpjm1      
            DO ji = 1,jpim1
              ztrax(ji,jj,jk)  = trn(ji+1,jj,jk,jn) - trn(ji,jj,jk,jn)
              zatrax(ji,jj,jk) = abs(ztrax(ji,jj,jk))
              ztray(ji,jj,jk)  = trn(ji,jj+1,jk,jn) - trn(ji,jj,jk,jn)
              zatray(ji,jj,jk) = abs(ztray(ji,jj,jk))
            ENDDO
          ENDDO
        ENDDO



C
C In the k direction
C
        DO jk=2,jpk
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1
              ztraz(ji,jj,jk)  = trn(ji,jj,jk-1,jn) - trn(ji,jj,jk,jn)
              zatraz(ji,jj,jk) = abs(ztraz(ji,jj,jk))
            ENDDO
          ENDDO
        ENDDO

C
C In the i and j direction
C

        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1
        zign = 0.5*(sign(1.,ztrax(ji,jj,jk)*ztrax(ji-1,jj,jk))+1)
        zkx(ji,jj,jk) = zign*(e1t(ji,jj)*((2*e1t(ji-1,jj)
     $                      +e1t(ji,jj))
     $                      *(trn(ji+1,jj,jk,jn)-trn(ji,jj,jk,jn))/
     $                      (e1t(ji+1,jj)+e1t(ji,jj))+(2*e1t(ji+1,jj)
     $                      +e1t(ji,jj))*(trn(ji,jj,jk,jn)
     $                      -trn(ji-1,jj,jk,jn))/(e1t(ji-1,jj)
     $                      +e1t(ji,jj)))/(e1t(ji-1,jj)+e1t(ji,jj)
     $                      +e1t(ji+1,jj)))


        zign = 0.5*(sign(1.,ztray(ji,jj,jk)*ztray(ji,jj-1,jk))+1)
        zky(ji,jj,jk) = zign*(e2t(ji,jj)*((2*e2t(ji,jj-1)
     $                      +e2t(ji,jj))
     $                      *(trn(ji,jj+1,jk,jn)-trn(ji,jj,jk,jn))/
     $                      (e2t(ji,jj+1)+e2t(ji,jj))+(2*e2t(ji,jj+1)
     $                      +e2t(ji,jj))*(trn(ji,jj,jk,jn)
     $                      -trn(ji,jj-1,jk,jn))/(e2t(ji,jj-1)
     $                      +e2t(ji,jj)))/(e2t(ji,jj-1)+e2t(ji,jj)
     $                      +e2t(ji,jj+1)))

            ENDDO
          ENDDO
        ENDDO        

C
C In the k direction
C
        DO jk=2,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1
        zign = 0.5*(sign(1.,ztraz(ji,jj,jk)*ztraz(ji,jj,jk+1))+1)
            zkz(ji,jj,jk) = zign*(fse3t(ji,jj,jk)*((2*fse3t(ji,jj,jk+1)
     $                      +fse3t(ji,jj,jk))*(-trn(ji,jj,jk-1,jn)
     $                      +trn(ji,jj,jk,jn))/(fse3t(ji,jj,jk-1)
     $                      +fse3t(ji,jj,jk))+(2*fse3t(ji,jj,jk-1)
     $                      +fse3t(ji,jj,jk))*(trn(ji,jj,jk,jn)
     $                      -trn(ji,jj,jk+1,jn))/(fse3t(ji,jj,jk+1)
     $                      +fse3t(ji,jj,jk)))/(fse3t(ji,jj,jk+1)
     $                      +fse3t(ji,jj,jk)+fse3t(ji,jj,jk-1)))
            ENDDO
          ENDDO
        ENDDO        



C
C 2. Slopes limitation
C --------------------
C

C
C In the i and j direction
C
        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1

              zkx(ji,jj,jk) = tmask(ji-1,jj,jk)*tmask(ji+1,jj,jk)*
     $                        sign(1.,zkx(ji,jj,jk)) * 
     $                        min(abs(zkx(ji,jj,jk)),
     $                        2*zatrax(ji-1,jj,jk),
     $                        2*zatrax(ji,jj,jk))

              zky(ji,jj,jk) = tmask(ji,jj-1,jk)*tmask(ji,jj+1,jk)*
     $                        sign(1.,zky(ji,jj,jk)) * 
     $                        min(abs(zky(ji,jj,jk)),
     $                        2*zatray(ji,jj-1,jk),
     $                        2*zatray(ji,jj,jk))
            ENDDO
          ENDDO
        ENDDO        

C
C In the k direction
C

        DO jk=2,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1
              zkz(ji,jj,jk) = tmask(ji,jj,jk+1)*tmask(ji,jj,jk-1)*
     $                        sign(1.,zkz(ji,jj,jk)) * 
     $                        min(abs(zkz(ji,jj,jk)),
     $                        2*zatraz(ji,jj,jk+1),
     $                        2*zatraz(ji,jj,jk))
            ENDDO
          ENDDO
        ENDDO        


C
C 3. T(i+1/2) step computation in the i direction
C
        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 3,jpi-2
           zdelta(ji,jj,jk)=(((trn(ji+1,jj,jk,jn)-trn(ji,jj,jk,jn))
     $                          /(e1t(ji+1,jj)+e1t(ji,jj)))
     $                          -((trn(ji,jj,jk,jn)-trn(ji-1,jj,jk,jn))
     $                          /(e1t(ji,jj)+e1t(ji-1,jj))))
     $                          /(e1t(ji-1,jj)+e1t(ji,jj)
     $                          +e1t(ji+1,jj))

        ENDDO
          ENDDO
        ENDDO


        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 4,jpi-3
           zign = -zdelta(ji-1,jj,jk)*zdelta(ji+1,jj,jk)
           zign = sign(1.,zign)
           zign1 = abs(trn(ji+1,jj,jk,jn)-trn(ji-1,jj,jk,jn))
     $         -zeps*min(abs(trn(ji+1,jj,jk,jn)),
     $         abs(trn(ji-1,jj,jk,jn)))
           zign1 = sign(1.,zign1)
           zetat(ji,jj,jk) = -0.25*(1+zign)*(1+zign1)
     $                           *((zdelta(ji+1,jj,jk)
     $                           -zdelta(ji-1,jj,jk))/(e1u(ji-1,jj)
     $                           +e1u(ji,jj)))*(((e1u(ji-1,jj)
     $                           )**3+(e1u(ji,jj)
     $                           )**3)/(trn(ji+1,jj,jk,jn)
     $                           -trn(ji-1,jj,jk,jn)+rtrn))

           zeta(ji,jj,jk) = min(zk1*(zetat(ji,jj,jk)-zk2),1.)
           zeta(ji,jj,jk) = max(0.,zeta(ji,jj,jk))


        ENDDO
          ENDDO
        ENDDO


        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 3,jpi-3


    ztrx(ji,jj,jk) = trn(ji,jj,jk,jn)+e1t(ji,jj)*(trn(ji+1,jj,jk,jn)
     $                -trn(ji,jj,jk,jn))/(e1t(ji,jj)+e1t(ji+1,jj))
     $                +(2*e1t(ji,jj)*e1t(ji+1,jj)*(trn(ji+1,jj,jk,jn)
     $                -trn(ji,jj,jk,jn))*((e1t(ji-1,jj)+e1t(ji,jj))
     $                /(2*e1t(ji,jj)+e1t(ji+1,jj))-(e1t(ji+1,jj)
     $                +e1t(ji+2,jj))/(2*e1t(ji+1,jj)+e1t(ji,jj)))
     $                /(e1t(ji,jj)+e1t(ji+1,jj))-e1t(ji,jj)
     $                *(e1t(ji,jj)+e1t(ji-1,jj))*zkx(ji+1,jj,jk)
     $                /(2*e1t(ji,jj)+e1t(ji+1,jj))+e1t(ji+1,jj)
     $                *(e1t(ji+1,jj)+e1t(ji+2,jj))*zkx(ji,jj,jk)
     $                /(2*e1t(ji+1,jj)+e1t(ji,jj)))/(e1t(ji-1,jj)
     $                +e1t(ji,jj)+e1t(ji+1,jj)+e1t(ji+2,jj))



        ztlx(ji+1,jj,jk) = ztrx(ji,jj,jk)*(1-zeta(ji+1,jj,jk))
     $                             +(trn(ji,jj,jk,jn)+0.5
     $                             *zkx(ji,jj,jk))*zeta(ji+1,jj,jk)  
        ztrx(ji,jj,jk) = ztrx(ji,jj,jk)*(1-zeta(ji,jj,jk))
     $                           +(trn(ji+1,jj,jk,jn)+0.5
     $                           *zkx(ji+1,jj,jk))*zeta(ji,jj,jk)


        ENDDO
          ENDDO
    ENDDO

C
C Right and Left TR(i+1/2) computation in the i direction
C
        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 4,jpi-3

        za =  ztrx(ji,jj,jk) - trn(ji,jj,jk,jn)
        zb =  ztlx(ji,jj,jk) - trn(ji,jj,jk,jn)
        zc = (ztrx(ji,jj,jk) - ztlx(ji,jj,jk))*
     $                   (trn(ji,jj,jk,jn)-0.5*(ztrx(ji,jj,jk) 
     $                   +ztlx(ji,jj,jk)))
        zd = ((ztrx(ji,jj,jk)-ztlx(ji,jj,jk))**2)/6.


        zign  = sign(1.,za*zb)
        ztlx(ji,jj,jk) = 0.5*(1+zign)
     $                           *trn(ji,jj,jk,jn)+0.5
     $                           *(1-zign)*ztlx(ji,jj,jk)
        ztrx(ji,jj,jk) = 0.5*(1+zign)
     $                           *trn(ji,jj,jk,jn)+0.5
     $                           *(1-zign)*ztrx(ji,jj,jk) 

        zign  = sign(1.,zc-zd)
        ztlx(ji,jj,jk) = 0.5*(1+zign)*(3*trn(ji,jj,jk,jn)
     $                           -2*ztrx(ji,jj,jk))+0.5*(1-zign)
     $                           *ztlx(ji,jj,jk)

        zign  = sign(1.,zc+zd)
        ztrx(ji,jj,jk) = 0.5*(1-zign)*(3*trn(ji,jj,jk,jn)
     $                           -2*ztlx(ji,jj,jk))+0.5*(1+zign)
     $                           *ztrx(ji,jj,jk)  

        za6x(ji,jj,jk) = 6*(trn(ji,jj,jk,jn)-0.5
     $                            *(ztrx(ji,jj,jk)+ztlx(ji,jj,jk)))



        ENDDO
          ENDDO
    ENDDO

C
C 4. T(i+1/2) step computation in the j direction
C


        DO jk=1,jpkm1
          DO jj = 3,jpj-2      
            DO ji = 2,jpim1
           zdelta(ji,jj,jk)=(((trn(ji,jj+1,jk,jn)-trn(ji,jj,jk,jn))
     $                          /(e2t(ji,jj+1)+e2t(ji,jj)))
     $                          -((trn(ji,jj,jk,jn)-trn(ji,jj-1,jk,jn))
     $                          /(e2t(ji,jj)+e2t(ji,jj-1))))
     $                          /(e2t(ji,jj-1)+e2t(ji,jj)
     $                          +e2t(ji,jj+1))

        ENDDO
          ENDDO
        ENDDO


        DO jk=1,jpkm1
          DO jj = 4,jpj-3      
            DO ji = 2,jpim1

           zign = -zdelta(ji,jj-1,jk)*zdelta(ji,jj+1,jk)
           zign = sign(1.,zign)
           zign1 = abs(trn(ji,jj+1,jk,jn)-trn(ji,jj-1,jk,jn))
     $         -zeps*min(abs(trn(ji,jj+1,jk,jn)),
     $         abs(trn(ji,jj-1,jk,jn)))
           zign1 = sign(1.,zign1)

           zeta(ji,jj,jk) = min(zk1*(zetat(ji,jj,jk)-zk2),1.)
           zeta(ji,jj,jk) = max(0.,zeta(ji,jj,jk))

        ENDDO
          ENDDO
        ENDDO


        DO jk=1,jpkm1
          DO jj = 3,jpj-3      
            DO ji = 2,jpim1

    ztry(ji,jj,jk) = trn(ji,jj,jk,jn)+e2t(ji,jj)*(trn(ji,jj+1,jk,jn)
     $                -trn(ji,jj,jk,jn))/(e2t(ji,jj)+e2t(ji,jj+1))
     $                +(2*e2t(ji,jj)*e2t(ji,jj+1)*(trn(ji,jj+1,jk,jn)
     $                -trn(ji,jj,jk,jn))*((e2t(ji,jj-1)+e2t(ji,jj))
     $                /(2*e2t(ji,jj)+e2t(ji,jj+1))-(e2t(ji,jj+1)
     $                +e2t(ji,jj+2))/(2*e2t(ji,jj+1)+e2t(ji,jj)))
     $                /(e2t(ji,jj)+e2t(ji,jj+1))-e2t(ji,jj)
     $                *(e2t(ji,jj)+e2t(ji,jj-1))*zky(ji,jj+1,jk)
     $                /(2*e2t(ji,jj)+e2t(ji,jj+1))+e2t(ji,jj+1)
     $                *(e2t(ji,jj+1)+e2t(ji,jj+2))*zky(ji,jj,jk)
     $                /(2*e2t(ji,jj+1)+e2t(ji,jj)))/(e2t(ji,jj-1)
     $                +e2t(ji,jj)+e2t(ji,jj+1)+e2t(ji,jj+2))




        ztly(ji,jj+1,jk) = ztry(ji,jj,jk)*(1-zeta(ji,jj+1,jk))
     $                             +(trn(ji,jj,jk,jn)+0.5
     $                             *zky(ji,jj,jk))*zeta(ji,jj+1,jk)  

        ztry(ji,jj,jk) = ztry(ji,jj,jk)*(1-zeta(ji,jj,jk))
     $                           +(trn(ji,jj+1,jk,jn)+0.5
     $                           *zky(ji,jj+1,jk))*zeta(ji,jj,jk)

        ENDDO
          ENDDO
    ENDDO

C
C Right and Left TR(j+1/2) computation in the j direction
C
        DO jk=1,jpkm1
          DO jj = 4,jpj-3      
            DO ji = 2,jpim1


        za =  ztry(ji,jj,jk) - trn(ji,jj,jk,jn)
        zb = ztly(ji,jj,jk) - trn(ji,jj,jk,jn)
        zc = (ztry(ji,jj,jk) - ztly(ji,jj,jk))*
     $               (trn(ji,jj,jk,jn)-0.5*(ztry(ji,jj,jk) 
     $                +ztly(ji,jj,jk)))
        zd = ((ztry(ji,jj,jk)-ztly(ji,jj,jk))**2)/6.


            zign  = sign(1.,za*zb)
        ztly(ji,jj,jk) = 0.5*(1+zign)
     $                           *trn(ji,jj,jk,jn)+0.5
     $                           *(1-zign)*ztly(ji,jj,jk)
        ztry(ji,jj,jk) = 0.5*(1+zign)
     $                           *trn(ji,jj,jk,jn)+0.5
     $                           *(1-zign)*ztry(ji,jj,jk) 

        zign  = sign(1.,zc-zd)
        ztly(ji,jj,jk) = 0.5*(1+zign)*(3*trn(ji,jj,jk,jn)
     $                           -2*ztry(ji,jj,jk))+0.5*(1-zign)
     $                           *ztly(ji,jj,jk) 

        zign  = sign(1.,zc+zd)
        ztry(ji,jj,jk) = 0.5*(1-zign)*(3*trn(ji,jj,jk,jn)
     $                           -2*ztly(ji,jj,jk))+0.5*(1+zign)
     $                           *ztry(ji,jj,jk)  

        za6y(ji,jj,jk) = 6*(trn(ji,jj,jk,jn)
     $                           -0.5*(ztry(ji,jj,jk)+ztly(ji,jj,jk)))

        ENDDO
      ENDDO
    ENDDO


C
C 5. T(i+1/2) step computation in the k direction
C

        DO jk=2,jpk-2
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1
          zdelta(ji,jj,jk)=(((-trn(ji,jj,jk-1,jn)+trn(ji,jj,jk,jn))
     $                         /(fse3t(ji,jj,jk-1)+fse3t(ji,jj,jk)))
     $                         -((-trn(ji,jj,jk,jn)+trn(ji,jj,jk+1,jn))
     $                         /(fse3t(ji,jj,jk)+fse3t(ji,jj,jk+1))))
     $                         /(fse3t(ji,jj,jk+1)+fse3t(ji,jj,jk)
     $                         +fse3t(ji,jj,jk-1))

        ENDDO
          ENDDO
        ENDDO


        DO jk=3,jpk-3
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1

           zign = -zdelta(ji,jj,jk-1)*zdelta(ji,jj,jk+1)
           zign = sign(1.,zign)
           zign1 = abs(trn(ji,jj,jk-1,jn)-trn(ji,jj,jk+1,jn))
     $         -zeps*min(abs(trn(ji,jj,jk-1,jn)),
     $         abs(trn(ji,jj,jk+1,jn)))
           zign1 = sign(1.,zign1)
           zetat(ji,jj,jk) = -0.25*(1+zign)*(1+zign1)
     $                           *((-zdelta(ji,jj,jk-1)
     $                           +zdelta(ji,jj,jk+1))/(fse3u(ji,jj,jk+1)
     $                           +fse3u(ji,jj,jk)))*(((fse3u(ji,jj,jk+1)
     $                           )**3+(fse3u(ji,jj,jk)
     $                           )**3)/(-trn(ji,jj,jk-1,jn)
     $                           +trn(ji,jj,jk+1,jn)+rtrn))

           zeta(ji,jj,jk) = min(zk1*(zetat(ji,jj,jk)-zk2),1.)
           zeta(ji,jj,jk) = max(0.,zeta(ji,jj,jk))
        ENDDO
          ENDDO
        ENDDO


        DO jk=3,jpk-3
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1


    ztrz(ji,jj,jk) = trn(ji,jj,jk,jn)+fse3t(ji,jj,jk)*
     $                   (-trn(ji,jj,jk-1,jn)+trn(ji,jj,jk,jn))
     $                   /(fse3t(ji,jj,jk)+fse3t(ji,jj,jk-1))
     $                   +(2*fse3t(ji,jj,jk)*fse3t(ji,jj,jk-1)
     $                   *(-trn(ji,jj,jk-1,jn)+trn(ji,jj,jk,jn))
     $                   *((fse3t(ji,jj,jk+1)+fse3t(ji,jj,jk))
     $                   /(2*fse3t(ji,jj,jk)+fse3t(ji,jj,jk-1))
     $                   -(fse3t(ji,jj,jk-1)+fse3t(ji,jj,jk-2))
     $                   /(2*fse3t(ji,jj,jk-1)+fse3t(ji,jj,jk)))
     $                   /(fse3t(ji,jj,jk)+fse3t(ji,jj,jk-1))
     $                   -fse3t(ji,jj,jk)*(fse3t(ji,jj,jk)
     $                   +fse3t(ji,jj,jk+1))*zkz(ji,jj,jk-1)
     $                   /(2*fse3t(ji,jj,jk)+fse3t(ji,jj,jk-1))
     $                   +fse3t(ji,jj,jk-1)*(fse3t(ji,jj,jk-1)
     $                   +fse3t(ji,jj,jk-2))*zkz(ji,jj,jk)
     $                   /(2*fse3t(ji,jj,jk-1)+fse3t(ji,jj,jk)))
     $                   /(fse3t(ji,jj,jk+1)+fse3t(ji,jj,jk)
     $                   +fse3t(ji,jj,jk-1)+fse3t(ji,jj,jk-2))




        ztlz(ji,jj,jk-1) = ztrz(ji,jj,jk)*(1-zeta(ji,jj,jk-1))
     $                             +(trn(ji,jj,jk,jn)-0.5
     $                             *zkz(ji,jj,jk))*zeta(ji,jj,jk-1)  

        ztrz(ji,jj,jk) = ztrz(ji,jj,jk)*(1-zeta(ji,jj,jk))
     $                           +(trn(ji,jj,jk-1,jn)-0.5
     $                           *zkz(ji,jj,jk-1))*zeta(ji,jj,jk)
        ENDDO
          ENDDO
    ENDDO

C
C Right and Left TR(k+1/2) computation in the k direction
C
        DO jk=3,jpk-4
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1

        za =  ztrz(ji,jj,jk) - trn(ji,jj,jk,jn)
        zb =  ztlz(ji,jj,jk) - trn(ji,jj,jk,jn)
        zc = (ztrz(ji,jj,jk) - ztlz(ji,jj,jk))*
     $                   (trn(ji,jj,jk,jn)-0.5*(ztrz(ji,jj,jk) 
     $                   +ztlz(ji,jj,jk)))
        zd = ((ztrz(ji,jj,jk)-ztlz(ji,jj,jk))**2)/6.


        zign  = sign(1.,za*zb)
        ztlz(ji,jj,jk) = 0.5*(1+zign)
     $                           *trn(ji,jj,jk,jn)+0.5
     $                           *(1-zign)*ztlz(ji,jj,jk)
        ztrz(ji,jj,jk) = 0.5*(1+zign)
     $                           *trn(ji,jj,jk,jn)+0.5
     $                           *(1-zign)*ztrz(ji,jj,jk) 

        zign  = sign(1.,zc-zd)
        ztlz(ji,jj,jk) = 0.5*(1+zign)*(3*trn(ji,jj,jk,jn)
     $                           -2*ztrz(ji,jj,jk))+0.5*(1-zign)
     $                           *ztlz(ji,jj,jk)

        zign  = sign(1.,zc+zd)
        ztrz(ji,jj,jk) = 0.5*(1-zign)*(3*trn(ji,jj,jk,jn)
     $                           -2*ztlz(ji,jj,jk))+0.5*(1+zign)
     $                           *ztrz(ji,jj,jk)  

        za6z(ji,jj,jk) = 6*(trn(ji,jj,jk,jn)-0.5
     $                            *(ztrz(ji,jj,jk)+ztlz(ji,jj,jk)))



        ENDDO
          ENDDO
    ENDDO

C
C 6. Advection terms
C ------------------
C

C
C Fluxes in the i and j direction
C
        
        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 4,jpi-4

          zigma = un(ji,jj,jk)*zstep/e1u(ji,jj)
              zeu   = e2u(ji,jj)* fse3u(ji,jj,jk)* un(ji,jj,jk)
          zign  = sign(1.,un(ji,jj,jk))


          zkx(ji,jj,jk) = tmask(ji,jj,jk)*(tmask(ji-2,jj,jk)
     $                        *tmask(ji+3,jj,jk)*tmask(ji-1,jj,jk)
     $                        *tmask(ji+1,jj,jk)*tmask(ji+2,jj,jk)
     $                        *(0.5*(1+zign)
     $                        *zeu*(ztrx(ji,jj,jk)-0.5
     $                        *zigma*(ztrx(ji,jj,jk)-ztlx(ji,jj,jk)
     $                        -(1-2.*zigma/3.)*za6x(ji,jj,jk)))+0.5
     $                        *(1-zign)*zeu*(ztlx(ji,jj,jk)-0.5*zigma
     $                        *(ztrx(ji+1,jj,jk)-ztlx(ji+1,jj,jk)
     $                        +(1+2.*zigma/3.)*za6x(ji+1,jj,jk))))
     $                        +((1-tmask(ji-2,jj,jk))
     $                        +(1-tmask(ji-1,jj,jk))
     $                        +(1-tmask(ji+1,jj,jk))
     $                        +(1-tmask(ji+2,jj,jk))
     $                        +(1-tmask(ji+3,jj,jk)))
     $                        *zeu*(0.5*(1+zign)*trn(ji,jj,jk,jn)+0.5
     $                        *(1-zign)*trn(ji+1,jj,jk,jn)
     $                        *tmask(ji+1,jj,jk)))

            ENDDO
          ENDDO
        ENDDO       

        DO jk=1,jpkm1
          DO jj = 4,jpj-4      
            DO ji = 2,jpim1

              zigma = vn(ji,jj,jk)*zstep/e2v(ji,jj)
              zev   = e1v(ji,jj)* fse3v(ji,jj,jk)* vn(ji,jj,jk)
              zign  = sign(1.,vn(ji,jj,jk)) 


          zky(ji,jj,jk) = tmask(ji,jj,jk)*(tmask(ji,jj-2,jk)
     $                        *tmask(ji,jj-1,jk)*tmask(ji,jj+1,jk)
     $                        *tmask(ji,jj+2,jk)
     $                        *tmask(ji,jj+3,jk)*(0.5*(1+zign)
     $                        *zev*(ztry(ji,jj,jk)-0.5
     $                        *zigma*(ztry(ji,jj,jk)-ztly(ji,jj,jk)
     $                        -(1-2.*zigma/3.)*za6y(ji,jj,jk)))+0.5
     $                        *(1-zign)*zev*(ztly(ji,jj,jk)-0.5*zigma
     $                        *(ztry(ji,jj+1,jk)-ztly(ji,jj+1,jk)
     $                        +(1+2.*zigma/3.)*za6y(ji,jj+1,jk))))
     $                        +((1-tmask(ji,jj-2,jk))
     $                        +(1-tmask(ji,jj-1,jk))
     $                        +(1-tmask(ji,jj+1,jk))
     $                        +(1-tmask(ji,jj+2,jk))
     $                        +(1-tmask(ji,jj+3,jk)))
     $                        *zeu*(0.5*(1+zign)*trn(ji,jj,jk,jn)+0.5
     $                        *(1-zign)*trn(ji,jj+1,jk,jn)
     $                        *tmask(ji,jj+1,jk)))
            ENDDO
          ENDDO
        ENDDO        


C
C Fluxes in the k direction
C

        DO jk=4,jpk-4
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1

          zigma = wn(ji,jj,jk)*zstep/fse3t(ji,jj,jk)
              zew   = e1t(ji,jj)* e2t(ji,jj)* wn(ji,jj,jk)
          zign  = sign(1.,wn(ji,jj,jk))


          zkz(ji,jj,jk) = tmask(ji,jj,jk)*(tmask(ji,jj,jk-3)
     $                        *tmask(ji,jj,jk-2)
     $                        *tmask(ji,jj,jk-1)
     $                        *tmask(ji,jj,jk+1)
     $                        *tmask(ji,jj,jk+2)
     $                        *tmask(ji,jj,jk+3)*(0.5*(1+zign)
     $                        *zew*(ztrz(ji,jj,jk)-0.5
     $                        *zigma*(-ztrz(ji,jj,jk)+ztlz(ji,jj,jk)
     $                        -(1-2.*zigma/3.)*za6z(ji,jj,jk)))+0.5
     $                        *(1-zign)*zew*(ztlz(ji,jj,jk)-0.5*zigma
     $                        *(-ztrz(ji,jj,jk-1)+ztlz(ji,jj,jk-1)
     $                        +(1+2.*zigma/3.)*za6z(ji,jj,jk-1))))
     $                        +((1-tmask(ji,jj,jk-3))
     $                        +(1-tmask(ji,jj,jk-2))
     $                        +(1-tmask(ji,jj,jk-1))
     $                        +(1-tmask(ji,jj,jk+1))
     $                        +(1-tmask(ji,jj,jk+2))
     $                        +(1-tmask(ji,jj,jk+3)))
     $                        *zeu*(0.5*(1+zign)*trn(ji,jj,jk,jn)+0.5
     $                        *(1-zign)*trn(ji,jj,jk-1,jn)
     $                        *tmask(ji,jj,jk-1)))
            ENDDO
          ENDDO
        ENDDO       


C
C Boundary conditions
C

          DO JK=1,jpkm1
            DO JJ=2,jpjm1

              zkx(1  ,JJ,JK)= 0.

          zigma = un(2,jj,jk)*zstep/e1u(2,jj)
              zeu   = e2u(2,jj)* fse3u(2,jj,jk)* un(2,jj,jk)
          zign  = sign(1.,un(2,jj,jk))
          zkx(2  ,JJ,JK)= tmask(2,jj,jk)*zeu*(0.5*(1+zign)
     $                        *trn(2,jj,jk,jn)+0.5*(1-zign)
     $                        *trn(3,jj,jk,jn)*tmask(3,jj,jk))

          zigma = un(3,jj,jk)*zstep/e1u(3,jj)
              zeu   = e2u(3,jj)* fse3u(3,jj,jk)* un(3,jj,jk)
          zign  = sign(1.,un(3,jj,jk))
          zkx(3  ,JJ,JK)= tmask(3,jj,jk)*zeu*(0.5*(1+zign)
     $                        *trn(3,jj,jk,jn)+0.5*(1-zign)
     $                        *trn(4,jj,jk,jn)*tmask(4,jj,jk))



          zigma = un(jpi-3,jj,jk)*zstep/e1u(jpi-3,jj)
              zeu   = e2u(jpi-3,jj)* fse3u(jpi-3,jj,jk)*un(jpi-3,jj,jk)
          zign  = sign(1.,un(jpi-3,jj,jk))


          zkx(JPI-3  ,JJ,JK)= tmask(jpi-3,jj,jk)*zeu*(0.5*(1+zign)
     $                        *trn(jpi-3,jj,jk,jn)+0.5*(1-zign)
     $                        *trn(jpi-2,jj,jk,jn)*tmask(jpi-2,jj,jk))

          zigma = un(jpi-2,jj,jk)*zstep/e1u(jpi-2,jj)
              zeu   = e2u(jpi-2,jj)* fse3u(jpi-2,jj,jk)* un(jpi-2,jj,jk)
          zign  = sign(1.,un(jpi-2,jj,jk))
          zkx(JPI-2  ,JJ,JK)= tmask(jpi-2,jj,jk)*zeu*(0.5*(1-zign)
     $                        *trn(jpi-1,jj,jk,jn)+0.5*(1+zign)
     $                        *trn(jpi-2,jj,jk,jn)*tmask(jpi-2,jj,jk))

          zigma = un(jpi-1,jj,jk)*zstep/e1u(jpi-1,jj)
              zeu   = e2u(jpi-1,jj)* fse3u(jpi-1,jj,jk)* un(jpi-1,jj,jk)
          zign  = sign(1.,un(jpi-1,jj,jk))
          zkx(JPI-1  ,JJ,JK)= tmask(jpi-1,jj,jk)*zeu*0.5*(1+zign)
     $                        *trn(jpi-1,jj,jk,jn)

              zkx(JPI,JJ,JK)=0.
              zky(1  ,JJ,JK)= 0.
              zky(JPI,JJ,JK)=0.

            ENDDO

            DO JI=1,JPI

              zkx(ji,1,JK)= 0.
              zkx(JI,JPJ,JK)=0.
              zky(ji,1,JK)= 0.

              zigma = vn(ji,2,jk)*zstep/e2v(ji,2)
              zev   = e1v(ji,2)* fse3v(ji,2,jk)* vn(ji,2,jk)
              zign  = sign(1.,vn(ji,2,jk)) 
          zky(ji,2,JK)= tmask(ji,2,jk)*zev*(0.5*(1+zign)
     $                        *trn(ji,2,jk,jn)+0.5*(1+zign)
     $                        *trn(ji,3,jk,jn)*tmask(ji,3,jk))

              zigma = vn(ji,3,jk)*zstep/e2v(ji,3)
              zev   = e1v(ji,3)* fse3v(ji,3,jk)* vn(ji,3,jk)
              zign  = sign(1.,vn(ji,3,jk)) 
          zky(ji,3,JK)= tmask(ji,3,jk)*zev*(0.5*(1+zign)
     $                        *trn(ji,3,jk,jn)+0.5*(1+zign)
     $                        *trn(ji,4,jk,jn)*tmask(ji,4,jk))

              zigma = vn(ji,jpj-3,jk)*zstep/e2v(ji,jpj-3)
              zev   = e1v(ji,jpj-3)* fse3v(ji,jpj-3,jk)* vn(ji,jpj-3,jk)
              zign  = sign(1.,vn(ji,jpj-3,jk)) 
          zky(ji,jpj-3,JK)= tmask(ji,jpj-3,jk)*zev*(0.5*(1+zign)
     $                        *trn(ji,jpj-3,jk,jn)+0.5*(1+zign)
     $                        *trn(ji,jpj-2,jk,jn)*tmask(ji,jpj-2,jk))

              zigma = vn(ji,jpj-2,jk)*zstep/e2v(ji,jpj-2)
              zev   = e1v(ji,jpj-2)* fse3v(ji,jpj-2,jk)* vn(ji,jpj-2,jk)
              zign  = sign(1.,vn(ji,jpj-2,jk)) 
          zky(ji,jpj-2,JK)= tmask(ji,jpj-2,jk)*zev*(0.5*(1+zign)
     $                        *trn(ji,jpj-2,jk,jn)+0.5*(1+zign)
     $                        *trn(ji,jpj-1,jk,jn)*tmask(ji,jpj-1,jk))

              zigma = vn(ji,jpj-1,jk)*zstep/e2v(ji,jpj-1)
              zev   = e1v(ji,jpj-1)* fse3v(ji,jpj-1,jk)* vn(ji,jpj-1,jk)
              zign  = sign(1.,vn(ji,jpj-1,jk)) 
          zky(ji,jpj-1,JK)= tmask(ji,jpj-1,jk)*zev*0.5*(1+zign)
     $                        *trn(ji,jpj-1,jk,jn)


              zky(JI,JPJ,JK)=0.
            ENDDO
          ENDDO

          DO jj=2,jpjm1
            DO ji=2,jpim1
           zkx(JI,JJ,JPK)=0.
           zky(JI,JJ,JPK)=0.
           zkz(ji,jj,1) = 0.

           zigma = wn(ji,jj,2)*zstep/fse3t(ji,jj,2)
               zew   = e1t(ji,jj)* e2t(ji,jj)* wn(ji,jj,2)
               zign  = sign(1.,wn(ji,jj,2)) 
           zkz(ji,jj,2) = tmask(ji,jj,2)*zew*(0.5*(1+zign)
     $                        *trn(ji,jj,2,jn)+0.5*(1+zign)
     $                        *trn(ji,jj,1,jn)*tmask(ji,jj,1))

           zigma = wn(ji,jj,3)*zstep/fse3t(ji,jj,3)
               zew   = e1t(ji,jj)* e2t(ji,jj)* wn(ji,jj,3)
               zign  = sign(1.,wn(ji,jj,3)) 
           zkz(ji,jj,3) = tmask(ji,jj,3)*zew*(0.5*(1+zign)
     $                        *trn(ji,jj,3,jn)+0.5*(1+zign)
     $                        *trn(ji,jj,2,jn)*tmask(ji,jj,2))

           zigma = wn(ji,jj,jpk-3)*zstep/fse3t(ji,jj,jpk-3)
               zew   = e1t(ji,jj)* e2t(ji,jj)* wn(ji,jj,jpk-3)
               zign  = sign(1.,wn(ji,jj,jpk-3)) 
           zkz(ji,jj,jpk-3) = tmask(ji,jj,jpk-3)*zew*(0.5*(1+zign)
     $                        *trn(ji,jj,jpk-3,jn)+0.5*(1+zign)
     $                        *trn(ji,jj,jpk-4,jn)*tmask(ji,jj,jpk-4))

           zigma = wn(ji,jj,jpk-2)*zstep/fse3t(ji,jj,jpk-2)
               zew   = e1t(ji,jj)* e2t(ji,jj)* wn(ji,jj,jpk-2)
               zign  = sign(1.,wn(ji,jj,jpk-2))    
           zkz(ji,jj,jpk-2) = tmask(ji,jj,jpk-2)*zew*(0.5*(1+zign)
     $                        *trn(ji,jj,jpk-2,jn)+0.5*(1+zign)
     $                        *trn(ji,jj,jpk-3,jn)*tmask(ji,jj,jpk-3))

           zigma = wn(ji,jj,2)*zstep/fse3t(ji,jj,jpk-1)
               zew   = e1t(ji,jj)* e2t(ji,jj)* wn(ji,jj,jpk-1)
               zign  = sign(1.,wn(ji,jj,jpk-1))    
           zkz(ji,jj,jpk-1) = tmask(ji,jj,jpk-1)*zew*(0.5*(1+zign)
     $                        *trn(ji,jj,jpk-1,jn)+0.5*(1+zign)
     $                        *trn(ji,jj,jpk-2,jn)*tmask(ji,jj,jpk-2))
 
           zkz(ji,jj,jpk)=0.
        ENDDO
      ENDDO




       IF (JPERIO.EQ.1) THEN

          DO JK=1,JPK
            DO JJ=1,JPJ
           ztrax(2,jj,jk)  = trn(3,jj,jk,jn) - trn(2,jj,jk,jn)
               zatrax(2,jj,jk) = abs(ztrax(2,jj,jk))

           ztrax(3,jj,jk)  = trn(4,jj,jk,jn) - trn(3,jj,jk,jn)
               zatrax(3,jj,jk) = abs(ztrax(3,jj,jk))

           ztrax(jpi-3,jj,jk)  = trn(jpi-2,jj,jk,jn) -
     $                               trn(jpi-3,jj,jk,jn)
               zatrax(jpi-3,jj,jk) = abs(ztrax(jpi-3,jj,jk))

           ztrax(jpi-2,jj,jk)  = trn(jpi-1,jj,jk,jn) -
     $                               trn(jpi-2,jj,jk,jn)
               zatrax(jpi-2,jj,jk) = abs(ztrax(jpi-2,jj,jk))

           ztrax(jpi-1,jj,jk)  = trn(2,jj,jk,jn) -
     $                               trn(jpi-1,jj,jk,jn)
               zatrax(jpi-1,jj,jk) = abs(ztrax(jpi-1,jj,jk))

           zkx(2,jj,jk) = 0.
           zkx(3,jj,jk) = 0.
           zkx(jpi-3,jj,jk) = 0.
           zkx(jpi-2,jj,jk) = 0.
           zkx(jpi-1,jj,jk) = 0.

           ztlx(2,jj,jk)  = 0.
           ztlx(3,jj,jk)  = 0.
           ztlx(jpi-2,jj,jk)  = 0.
           ztlx(jpi-1,jj,jk)  = 0.

           ztrx(2,jj,jk)  = 0.
           ztrx(3,jj,jk)  = 0.
           ztrx(jpi-3,jj,jk)  = 0.
           ztrx(jpi-2,jj,jk)  = 0.
           ztrx(jpi-1,jj,jk)  = 0.

           za6x(2,jj,jk)  = 0.
           za6x(3,jj,jk)  = 0.
           za6x(jpi-3,jj,jk)  = 0.
           za6x(jpi-2,jj,jk)  = 0.
           za6x(jpi-1,jj,jk)  = 0.
      
           zkx1 = 0
        ENDDO
      ENDDO

          DO JK=1,JPK
            DO JJ=1,JPJ

        zign = 0.5*(sign(1.,ztrax(2,jj,jk)*ztrax(jpi-1,jj,jk))+1)
        zkx(2,jj,jk) = zign*(e1t(2,jj)*((2*e1t(jpi-1,jj)
     $                      +e1t(2,jj))
     $                      *(trn(3,jj,jk,jn)-trn(2,jj,jk,jn))/
     $                      (e1t(3,jj)+e1t(2,jj))+(2*e1t(3,jj)
     $                      +e1t(2,jj))*(trn(2,jj,jk,jn)
     $                      -trn(jpi-1,jj,jk,jn))/(e1t(jpi-1,jj)
     $                      +e1t(2,jj)))/(e1t(jpi-1,jj)+e1t(2,jj)
     $                      +e1t(3,jj)))

        zign = 0.5*(sign(1.,ztrax(3,jj,jk)*ztrax(2,jj,jk))+1)
        zkx(3,jj,jk) = zign*(e1t(3,jj)*((2*e1t(2,jj)
     $                      +e1t(3,jj))
     $                      *(trn(4,jj,jk,jn)-trn(3,jj,jk,jn))/
     $                      (e1t(4,jj)+e1t(3,jj))+(2*e1t(4,jj)
     $                      +e1t(3,jj))*(trn(3,jj,jk,jn)
     $                      -trn(2,jj,jk,jn))/(e1t(2,jj)
     $                      +e1t(3,jj)))/(e1t(2,jj)+e1t(3,jj)
     $                      +e1t(4,jj)))

        zign = 0.5*(sign(1.,ztrax(4,jj,jk)*ztrax(3,jj,jk))+1)
        zkx1(4,jj,jk) = zign*(e1t(4,jj)*((2*e1t(3,jj)
     $                      +e1t(4,jj))
     $                      *(trn(5,jj,jk,jn)-trn(4,jj,jk,jn))/
     $                      (e1t(5,jj)+e1t(4,jj))+(2*e1t(5,jj)
     $                      +e1t(4,jj))*(trn(4,jj,jk,jn)
     $                      -trn(3,jj,jk,jn))/(e1t(3,jj)
     $                      +e1t(4,jj)))/(e1t(3,jj)+e1t(4,jj)
     $                      +e1t(5,jj)))


        zign = 0.5*(sign(1.,ztrax(jpi-3,jj,jk)*
     $                          ztrax(jpi-4,jj,jk))+1)
        zkx(jpi-3,jj,jk) = zign*(e1t(jpi-3,jj)*((2*e1t(jpi-4,jj)
     $                      +e1t(jpi-3,jj))
     $                      *(trn(jpi-2,jj,jk,jn)-trn(jpi-3,jj,jk,jn))/
     $                      (e1t(jpi-2,jj)+e1t(jpi-3,jj))
     $                      +(2*e1t(jpi-2,jj)
     $                      +e1t(jpi-3,jj))*(trn(jpi-3,jj,jk,jn)
     $                      -trn(jpi-4,jj,jk,jn))/(e1t(jpi-4,jj)
     $                      +e1t(jpi-3,jj)))/(e1t(jpi-4,jj)
     $                      +e1t(jpi-3,jj)+e1t(jpi-2,jj)))



        zign = 0.5*(sign(1.,ztrax(jpi-2,jj,jk)*
     $                          ztrax(jpi-3,jj,jk))+1)
        zkx(jpi-2,jj,jk) = zign*(e1t(jpi-2,jj)*((2*e1t(jpi-3,jj)
     $                      +e1t(jpi-2,jj))
     $                      *(trn(jpi-1,jj,jk,jn)-trn(jpi-2,jj,jk,jn))/
     $                      (e1t(jpi-1,jj)+e1t(jpi-2,jj))+
     $                      (2*e1t(jpi-1,jj)
     $                      +e1t(jpi-2,jj))*(trn(jpi-2,jj,jk,jn)
     $                      -trn(jpi-3,jj,jk,jn))/(e1t(jpi-3,jj)
     $                      +e1t(jpi-2,jj)))/(e1t(jpi-3,jj)
     $                      +e1t(jpi-2,jj)+e1t(jpi-1,jj)))

        zign = 0.5*(sign(1.,ztrax(jpi-1,jj,jk)*
     $                          ztrax(jpi-2,jj,jk))+1)
        zkx(jpi-1,jj,jk) = zign*(e1t(jpi-1,jj)*((2*e1t(jpi-2,jj)
     $                      +e1t(jpi-1,jj))
     $                      *(trn(2,jj,jk,jn)-trn(jpi-1,jj,jk,jn))/
     $                      (e1t(2,jj)+e1t(jpi-1,jj))+(2*e1t(2,jj)
     $                      +e1t(jpi-1,jj))*(trn(jpi-1,jj,jk,jn)
     $                      -trn(jpi-2,jj,jk,jn))/(e1t(jpi-2,jj)
     $                      +e1t(jpi-1,jj)))/(e1t(jpi-2,jj)
     $                      +e1t(jpi-1,jj)+e1t(2,jj)))

        ENDDO
      ENDDO


          DO JK=1,JPK
            DO JJ=1,JPJ
              zkx(2,jj,jk) = tmask(jpi-1,jj,jk)*tmask(3,jj,jk)*
     $                        sign(1.,zkx(2,jj,jk)) * 
     $                        min(abs(zkx(2,jj,jk)),
     $                        2*zatrax(jpi-1,jj,jk),
     $                        2*zatrax(2,jj,jk))

              zkx(3,jj,jk) = tmask(2,jj,jk)*tmask(4,jj,jk)*
     $                        sign(1.,zkx(3,jj,jk)) * 
     $                        min(abs(zkx(3,jj,jk)),
     $                        2*zatrax(2,jj,jk),
     $                        2*zatrax(3,jj,jk))

              zkx1(4,jj,jk) = tmask(3,jj,jk)*tmask(5,jj,jk)*
     $                        sign(1.,zkx1(4,jj,jk)) * 
     $                        min(abs(zkx1(4,jj,jk)),
     $                        2*zatrax(3,jj,jk),
     $                        2*zatrax(4,jj,jk))

              zkx(jpi-3,jj,jk) = tmask(jpi-4,jj,jk)*tmask(jpi-2,jj,jk)*
     $                        sign(1.,zkx(jpi-3,jj,jk)) * 
     $                        min(abs(zkx(jpi-3,jj,jk)),
     $                        2*zatrax(jpi-4,jj,jk),
     $                        2*zatrax(jpi-3,jj,jk))

              zkx(jpi-2,jj,jk) = tmask(jpi-3,jj,jk)*tmask(jpi-1,jj,jk)*
     $                        sign(1.,zkx(jpi-2,jj,jk)) * 
     $                        min(abs(zkx(jpi-2,jj,jk)),
     $                        2*zatrax(jpi-3,jj,jk),
     $                        2*zatrax(jpi-2,jj,jk))

              zkx(jpi-1,jj,jk) = tmask(jpi-2,jj,jk)*tmask(2,jj,jk)*
     $                        sign(1.,zkx(jpi-1,jj,jk)) * 
     $                        min(abs(zkx(jpi-1,jj,jk)),
     $                        2*zatrax(jpi-2,jj,jk),
     $                        2*zatrax(jpi-1,jj,jk))


        ENDDO
      ENDDO


          DO JK=1,JPK
            DO JJ=1,JPJ
           zdelta(2,jj,jk)=(((trn(3,jj,jk,jn)-trn(2,jj,jk,jn))
     $                          /(e1t(3,jj)+e1t(2,jj)))
     $                          -((trn(2,jj,jk,jn)-trn(jpi-1,jj,jk,jn))
     $                          /(e1t(2,jj)+e1t(jpi-1,jj))))
     $                          /(e1t(jpi-1,jj)+e1t(2,jj)
     $                          +e1t(3,jj))

           zdelta(3,jj,jk)=(((trn(4,jj,jk,jn)-trn(3,jj,jk,jn))
     $                          /(e1t(4,jj)+e1t(3,jj)))
     $                          -((trn(3,jj,jk,jn)-trn(2,jj,jk,jn))
     $                          /(e1t(3,jj)+e1t(2,jj))))
     $                          /(e1t(2,jj)+e1t(3,jj)
     $                          +e1t(4,jj))

           zdelta(4,jj,jk)=(((trn(5,jj,jk,jn)-trn(4,jj,jk,jn))
     $                          /(e1t(5,jj)+e1t(4,jj)))
     $                          -((trn(4,jj,jk,jn)-trn(3,jj,jk,jn))
     $                          /(e1t(4,jj)+e1t(3,jj))))
     $                          /(e1t(3,jj)+e1t(4,jj)
     $                          +e1t(5,jj))

           zdelta(jpi-4,jj,jk)=(((trn(jpi-3,jj,jk,jn)
     $                          -trn(jpi-4,jj,jk,jn))
     $                          /(e1t(jpi-3,jj)+e1t(jpi-4,jj)))
     $                          -((trn(jpi-4,jj,jk,jn)
     $                          -trn(jpi-5,jj,jk,jn))
     $                          /(e1t(jpi-4,jj)+e1t(jpi-5,jj))))
     $                          /(e1t(jpi-5,jj)+e1t(jpi-4,jj)
     $                          +e1t(jpi-3,jj))

           zdelta(jpi-3,jj,jk)=(((trn(jpi-2,jj,jk,jn)
     $                          -trn(jpi-3,jj,jk,jn))
     $                          /(e1t(jpi-2,jj)+e1t(jpi-3,jj)))
     $                          -((trn(jpi-3,jj,jk,jn)
     $                          -trn(jpi-4,jj,jk,jn))
     $                          /(e1t(jpi-3,jj)+e1t(jpi-4,jj))))
     $                          /(e1t(jpi-1,jj)+e1t(2,jj)
     $                          +e1t(jpi-2,jj))

           zdelta(jpi-2,jj,jk)=(((trn(jpi-1,jj,jk,jn)
     $                          -trn(jpi-2,jj,jk,jn))
     $                          /(e1t(jpi-1,jj)+e1t(jpi-2,jj)))
     $                          -((trn(jpi-2,jj,jk,jn)
     $                          -trn(jpi-3,jj,jk,jn))
     $                          /(e1t(jpi-2,jj)+e1t(jpi-3,jj))))
     $                          /(e1t(jpi-3,jj)+e1t(jpi-2,jj)
     $                          +e1t(jpi-1,jj))

           zdelta(jpi-1,jj,jk)=(((trn(2,jj,jk,jn)
     $                          -trn(jpi-1,jj,jk,jn))
     $                          /(e1t(2,jj)+e1t(jpi-1,jj)))
     $                          -((trn(jpi-1,jj,jk,jn)
     $                          -trn(jpi-2,jj,jk,jn))
     $                          /(e1t(jpi-1,jj)+e1t(jpi-2,jj))))
     $                          /(e1t(jpi-2,jj)+e1t(jpi-1,jj)
     $                          +e1t(2,jj))


           ENDDO
      ENDDO

          DO JK=1,JPK
            DO JJ=1,JPJ
           zign = -zdelta(jpi-1,jj,jk)*zdelta(3,jj,jk)
           zign = sign(1.,zign)
           zign1 = abs(trn(3,jj,jk,jn)-trn(jpi-1,jj,jk,jn))
     $         -zeps*min(abs(trn(3,jj,jk,jn)),
     $         abs(trn(jpi-1,jj,jk,jn)))
           zign1 = sign(1.,zign1)
           zetat(2,jj,jk) = -0.25*(1+zign)*(1+zign1)
     $                           *((zdelta(3,jj,jk)
     $                           -zdelta(jpi-1,jj,jk))/(e1u(jpi-1,jj)
     $                           +e1u(2,jj)))*(((e1u(jpi-1,jj)
     $                           )**3+(e1u(2,jj)
     $                           )**3)/(trn(3,jj,jk,jn)
     $                           -trn(jpi-1,jj,jk,jn)+rtrn))

           zeta(2,jj,jk) = min(zk1*(zetat(2,jj,jk)-zk2),1.)
           zeta(2,jj,jk) = max(0.,zeta(2,jj,jk))

           zign = -zdelta(2,jj,jk)*zdelta(4,jj,jk)
           zign = sign(1.,zign)
           zign1 = abs(trn(4,jj,jk,jn)-trn(2,jj,jk,jn))
     $         -zeps*min(abs(trn(4,jj,jk,jn)),
     $         abs(trn(2,jj,jk,jn)))
           zign1 = sign(1.,zign1)
           zetat(3,jj,jk) = -0.25*(1+zign)*(1+zign1)
     $                           *((zdelta(4,jj,jk)
     $                           -zdelta(2,jj,jk))/(e1u(2,jj)
     $                           +e1u(3,jj)))*(((e1u(2,jj)
     $                           )**3+(e1u(3,jj)
     $                           )**3)/(trn(4,jj,jk,jn)
     $                           -trn(2,jj,jk,jn)+rtrn))

           zeta(3,jj,jk) = min(zk1*(zetat(3,jj,jk)-zk2),1.)
           zeta(3,jj,jk) = max(0.,zeta(3,jj,jk))

           zign = -zdelta(jpi-4,jj,jk)*zdelta(jpi-2,jj,jk)
           zign = sign(1.,zign)
           zign1 = abs(trn(jpi-2,jj,jk,jn)-trn(jpi-4,jj,jk,jn))
     $         -zeps*min(abs(trn(jpi-2,jj,jk,jn)),
     $         abs(trn(jpi-4,jj,jk,jn)))
           zign1 = sign(1.,zign1)
           zetat(jpi-3,jj,jk) = -0.25*(1+zign)*(1+zign1)
     $                           *((zdelta(jpi-2,jj,jk)
     $                           -zdelta(jpi-4,jj,jk))/(e1u(jpi-4,jj)
     $                           +e1u(jpi-3,jj)))*(((e1u(jpi-4,jj)
     $                           )**3+(e1u(jpi-3,jj)
     $                           )**3)/(trn(jpi-2,jj,jk,jn)
     $                           -trn(jpi-4,jj,jk,jn)+rtrn))

           zeta(jpi-3,jj,jk) = min(zk1*(zetat(jpi-3,jj,jk)-zk2),1.)
           zeta(jpi-3,jj,jk) = max(0.,zeta(jpi-3,jj,jk))

           zign = -zdelta(jpi-3,jj,jk)*zdelta(jpi-1,jj,jk)
           zign = sign(1.,zign)
           zign1 = abs(trn(jpi-1,jj,jk,jn)-trn(jpi-3,jj,jk,jn))
     $         -zeps*min(abs(trn(jpi-1,jj,jk,jn)),
     $         abs(trn(jpi-3,jj,jk,jn)))
           zign1 = sign(1.,zign1)
           zetat(jpi-2,jj,jk) = -0.25*(1+zign)*(1+zign1)
     $                           *((zdelta(jpi-1,jj,jk)
     $                           -zdelta(jpi-3,jj,jk))/(e1u(jpi-3,jj)
     $                           +e1u(jpi-2,jj)))*(((e1u(jpi-3,jj)
     $                           )**3+(e1u(jpi-2,jj)
     $                           )**3)/(trn(jpi-1,jj,jk,jn)
     $                           -trn(jpi-3,jj,jk,jn)+rtrn))

           zeta(jpi-2,jj,jk) = min(zk1*(zetat(jpi-2,jj,jk)-zk2),1.)
           zeta(jpi-2,jj,jk) = max(0.,zeta(jpi-2,jj,jk))

           zign = -zdelta(jpi-2,jj,jk)*zdelta(2,jj,jk)
           zign = sign(1.,zign)
           zign1 = abs(trn(2,jj,jk,jn)-trn(jpi-2,jj,jk,jn))
     $         -zeps*min(abs(trn(2,jj,jk,jn)),
     $         abs(trn(jpi-2,jj,jk,jn)))
           zign1 = sign(1.,zign1)
           zetat(jpi-1,jj,jk) = -0.25*(1+zign)*(1+zign1)
     $                           *((zdelta(2,jj,jk)
     $                           -zdelta(jpi-2,jj,jk))/(e1u(jpi-2,jj)
     $                           +e1u(jpi-1,jj)))*(((e1u(jpi-2,jj)
     $                           )**3+(e1u(jpi-1,jj)
     $                           )**3)/(trn(2,jj,jk,jn)
     $                           -trn(jpi-2,jj,jk,jn)+rtrn))

           zeta(jpi-1,jj,jk) = min(zk1*(zetat(jpi-1,jj,jk)-zk2),1.)
           zeta(jpi-1,jj,jk) = max(0.,zeta(jpi-1,jj,jk))


       ENDDO
      ENDDO

          DO JK=1,JPK
            DO JJ=1,JPJ

    ztrx(2,jj,jk) = trn(2,jj,jk,jn)+e1t(2,jj)*(trn(3,jj,jk,jn)
     $                -trn(2,jj,jk,jn))/(e1t(2,jj)+e1t(3,jj))
     $                +(2*e1t(2,jj)*e1t(3,jj)*(trn(3,jj,jk,jn)
     $                -trn(2,jj,jk,jn))*((e1t(jpi-1,jj)+e1t(2,jj))
     $                /(2*e1t(2,jj)+e1t(3,jj))-(e1t(3,jj)
     $                +e1t(4,jj))/(2*e1t(3,jj)+e1t(2,jj)))
     $                /(e1t(2,jj)+e1t(3,jj))-e1t(2,jj)
     $                *(e1t(2,jj)+e1t(jpi-1,jj))*zkx(3,jj,jk)
     $                /(2*e1t(2,jj)+e1t(3,jj))+e1t(3,jj)
     $                *(e1t(3,jj)+e1t(4,jj))*zkx(2,jj,jk)
     $                /(2*e1t(3,jj)+e1t(2,jj)))/(e1t(jpi-1,jj)
     $                +e1t(2,jj)+e1t(3,jj)+e1t(4,jj))




        ztlx(3,jj,jk) = ztrx(2,jj,jk)*(1-zeta(3,jj,jk))
     $                             +(trn(2,jj,jk,jn)+0.5
     $                             *zkx(2,jj,jk))*zeta(3,jj,jk) 

        ztrx(2,jj,jk) = ztrx(2,jj,jk)*(1-zeta(2,jj,jk))
     $                           +(trn(3,jj,jk,jn)+0.5
     $                           *zkx(3,jj,jk))*zeta(2,jj,jk) 

    ztrx(3,jj,jk) = trn(3,jj,jk,jn)+e1t(3,jj)*(trn(4,jj,jk,jn)
     $                -trn(3,jj,jk,jn))/(e1t(3,jj)+e1t(4,jj))
     $                +(2*e1t(3,jj)*e1t(4,jj)*(trn(4,jj,jk,jn)
     $                -trn(3,jj,jk,jn))*((e1t(2,jj)+e1t(3,jj))
     $                /(2*e1t(3,jj)+e1t(4,jj))-(e1t(4,jj)
     $                +e1t(5,jj))/(2*e1t(4,jj)+e1t(3,jj)))
     $                /(e1t(3,jj)+e1t(4,jj))-e1t(3,jj)
     $                *(e1t(3,jj)+e1t(2,jj))*zkx1(4,jj,jk)
     $                /(2*e1t(3,jj)+e1t(4,jj))+e1t(4,jj)
     $                *(e1t(4,jj)+e1t(5,jj))*zkx(3,jj,jk)
     $                /(2*e1t(4,jj)+e1t(3,jj)))/(e1t(2,jj)
     $                +e1t(3,jj)+e1t(4,jj)+e1t(5,jj))


        ztrx(3,jj,jk) = ztrx(3,jj,jk)*(1-zeta(3,jj,jk))
     $                           +(trn(4,jj,jk,jn)+0.5
     $                           *zkx1(4,jj,jk))*zeta(3,jj,jk)


    ztrx(jpi-3,jj,jk) = trn(jpi-3,jj,jk,jn)+e1t(jpi-3,jj)
     $                     *(trn(jpi-2,jj,jk,jn)
     $                     -trn(jpi-3,jj,jk,jn))/(e1t(jpi-3,jj)
     $                     +e1t(jpi-2,jj))
     $                     +(2*e1t(jpi-3,jj)*e1t(jpi-2,jj)
     $                     *(trn(jpi-2,jj,jk,jn)
     $                     -trn(jpi-3,jj,jk,jn))*((e1t(jpi-4,jj)
     $                     +e1t(jpi-3,jj))
     $                     /(2*e1t(jpi-3,jj)+e1t(jpi-2,jj))
     $                     -(e1t(jpi-2,jj)
     $                     +e1t(jpi-1,jj))/(2*e1t(jpi-2,jj)
     $                     +e1t(jpi-3,jj)))
     $                     /(e1t(jpi-3,jj)+e1t(jpi-2,jj))
     $                     -e1t(jpi-3,jj)
     $                     *(e1t(jpi-3,jj)+e1t(jpi-4,jj))
     $                     *zkx(jpi-2,jj,jk)
     $                     /(2*e1t(jpi-3,jj)+e1t(jpi-2,jj))
     $                     +e1t(jpi-2,jj)
     $                     *(e1t(jpi-2,jj)+e1t(jpi-1,jj))
     $                     *zkx(jpi-3,jj,jk)
     $                     /(2*e1t(jpi-2,jj)+e1t(jpi-3,jj)))
     $                     /(e1t(jpi-4,jj)
     $                     +e1t(jpi-3,jj)+e1t(jpi-2,jj)
     $                     +e1t(jpi-1,jj))



        ztlx(jpi-2,jj,jk) = ztrx(jpi-3,jj,jk)
     $                             *(1-zeta(jpi-2,jj,jk))
     $                             +(trn(jpi-3,jj,jk,jn)+0.5
     $                             *zkx(jpi-3,jj,jk))*zeta(jpi-2,jj,jk)  

        ztrx(jpi-3,jj,jk) = ztrx(jpi-3,jj,jk)
     $                           *(1-zeta(jpi-3,jj,jk))
     $                           +(trn(jpi-2,jj,jk,jn)+0.5
     $                           *zkx(jpi-2,jj,jk))*zeta(jpi-3,jj,jk)



    ztrx(jpi-2,jj,jk) = trn(jpi-2,jj,jk,jn)+e1t(jpi-2,jj)
     $                *(trn(jpi-1,jj,jk,jn)
     $                -trn(jpi-2,jj,jk,jn))/(e1t(jpi-2,jj)
     $                +e1t(jpi-1,jj))
     $                +(2*e1t(jpi-2,jj)*e1t(jpi-1,jj)
     $                *(trn(jpi-1,jj,jk,jn)
     $                -trn(jpi-2,jj,jk,jn))*((e1t(jpi-3,jj)
     $                +e1t(jpi-2,jj))
     $                /(2*e1t(jpi-2,jj)+e1t(jpi-1,jj))-(e1t(jpi-1,jj)
     $                +e1t(2,jj))/(2*e1t(jpi-1,jj)+e1t(jpi-2,jj)))
     $                /(e1t(jpi-2,jj)+e1t(jpi-1,jj))-e1t(jpi-2,jj)
     $                *(e1t(jpi-2,jj)+e1t(jpi-3,jj))*zkx(jpi-1,jj,jk)
     $                /(2*e1t(jpi-2,jj)+e1t(jpi-1,jj))+e1t(jpi-1,jj)
     $                *(e1t(jpi-1,jj)+e1t(2,jj))*zkx(jpi-2,jj,jk)
     $                /(2*e1t(jpi-1,jj)+e1t(jpi-2,jj)))/(e1t(jpi-3,jj)
     $                +e1t(jpi-2,jj)+e1t(jpi-1,jj)+e1t(2,jj))




        ztlx(jpi-1,jj,jk) = ztrx(jpi-2,jj,jk)
     $                             *(1-zeta(jpi-1,jj,jk))
     $                             +(trn(jpi-2,jj,jk,jn)+0.5
     $                             *zkx(jpi-2,jj,jk))*zeta(jpi-1,jj,jk)  

        ztrx(jpi-2,jj,jk) = ztrx(jpi-2,jj,jk)
     $                           *(1-zeta(jpi-2,jj,jk))
     $                           +(trn(jpi-1,jj,jk,jn)+0.5
     $                           *zkx(jpi-1,jj,jk))*zeta(jpi-2,jj,jk)

    ztrx(jpi-1,jj,jk) = trn(jpi-1,jj,jk,jn)+e1t(jpi-1,jj)
     $                *(trn(2,jj,jk,jn)
     $                -trn(jpi-1,jj,jk,jn))/(e1t(jpi-1,jj)+e1t(2,jj))
     $                +(2*e1t(jpi-1,jj)*e1t(2,jj)*(trn(2,jj,jk,jn)
     $                -trn(jpi-1,jj,jk,jn))*((e1t(jpi-2,jj)
     $                +e1t(jpi-1,jj))
     $                /(2*e1t(jpi-1,jj)+e1t(2,jj))-(e1t(2,jj)
     $                +e1t(3,jj))/(2*e1t(2,jj)+e1t(jpi-1,jj)))
     $                /(e1t(jpi-1,jj)+e1t(2,jj))-e1t(jpi-1,jj)
     $                *(e1t(jpi-1,jj)+e1t(jpi-2,jj))*zkx(2,jj,jk)
     $                /(2*e1t(jpi-1,jj)+e1t(2,jj))+e1t(2,jj)
     $                *(e1t(2,jj)+e1t(3,jj))*zkx(jpi-1,jj,jk)
     $                /(2*e1t(2,jj)+e1t(jpi-1,jj)))/(e1t(jpi-2,jj)
     $                +e1t(jpi-1,jj)+e1t(2,jj)+e1t(3,jj))



        ztlx(2,jj,jk) = ztrx(jpi-1,jj,jk)*(1-zeta(2,jj,jk))
     $                             +(trn(jpi-1,jj,jk,jn)+0.5
     $                             *zkx(jpi-1,jj,jk))*zeta(2,jj,jk)  

        ztrx(jpi-1,jj,jk) = ztrx(jpi-1,jj,jk)
     $                           *(1-zeta(jpi-1,jj,jk))
     $                           +(trn(2,jj,jk,jn)+0.5
     $                           *zkx(2,jj,jk))*zeta(jpi-1,jj,jk)


           ENDDO
          ENDDO


          DO JK=1,JPK
            DO JJ=1,JPJ
        za =  ztrx(2,jj,jk) - trn(2,jj,jk,jn)
        zb =  ztlx(2,jj,jk) - trn(2,jj,jk,jn)
        zc = (ztrx(2,jj,jk) - ztlx(2,jj,jk))*
     $                   (trn(2,jj,jk,jn)-0.5*(ztrx(2,jj,jk) 
     $                   +ztlx(2,jj,jk)))
        zd = ((ztrx(2,jj,jk)-ztlx(2,jj,jk))**2)/6.
        zign  = sign(1.,za*zb)
        ztlx(2,jj,jk) = 0.5*(1+zign)
     $                           *trn(2,jj,jk,jn)+0.5
     $                           *(1-zign)*ztlx(2,jj,jk)
        ztrx(2,jj,jk) = 0.5*(1+zign)
     $                           *trn(2,jj,jk,jn)+0.5
     $                           *(1-zign)*ztrx(2,jj,jk) 
        zign  = sign(1.,zc-zd)
        ztlx(2,jj,jk) = 0.5*(1+zign)*(3*trn(2,jj,jk,jn)
     $                           -2*ztrx(2,jj,jk))+0.5*(1-zign)
     $                           *ztlx(2,jj,jk)
        zign  = sign(1.,zc+zd)
        ztrx(2,jj,jk) = 0.5*(1-zign)*(3*trn(2,jj,jk,jn)
     $                           -2*ztlx(2,jj,jk))+0.5*(1+zign)
     $                           *ztrx(2,jj,jk)  

        za6x(2,jj,jk) = 6*(trn(2,jj,jk,jn)-0.5
     $                            *(ztrx(2,jj,jk)+ztlx(2,jj,jk)))

        za =  ztrx(3,jj,jk) - trn(3,jj,jk,jn)
        zb =  ztlx(3,jj,jk) - trn(3,jj,jk,jn)
        zc = (ztrx(3,jj,jk) - ztlx(3,jj,jk))*
     $                   (trn(3,jj,jk,jn)-0.5*(ztrx(3,jj,jk) 
     $                   +ztlx(3,jj,jk)))
        zd = ((ztrx(3,jj,jk)-ztlx(3,jj,jk))**2)/6.
        zign  = sign(1.,za*zb)
        ztlx(3,jj,jk) = 0.5*(1+zign)
     $                           *trn(3,jj,jk,jn)+0.5
     $                           *(1-zign)*ztlx(3,jj,jk)
        ztrx(3,jj,jk) = 0.5*(1+zign)
     $                           *trn(3,jj,jk,jn)+0.5
     $                           *(1-zign)*ztrx(3,jj,jk) 
        zign  = sign(1.,zc-zd)
        ztlx(3,jj,jk) = 0.5*(1+zign)*(3*trn(3,jj,jk,jn)
     $                           -2*ztrx(3,jj,jk))+0.5*(1-zign)
     $                           *ztlx(3,jj,jk)
        zign  = sign(1.,zc+zd)
        ztrx(3,jj,jk) = 0.5*(1-zign)*(3*trn(3,jj,jk,jn)
     $                           -2*ztlx(3,jj,jk))+0.5*(1+zign)
     $                           *ztrx(3,jj,jk)  
        za6x(3,jj,jk) = 6*(trn(3,jj,jk,jn)-0.5
     $                            *(ztrx(3,jj,jk)+ztlx(3,jj,jk)))


        za =  ztrx(jpi-3,jj,jk) - trn(jpi-3,jj,jk,jn)
        zb =  ztlx(jpi-3,jj,jk) - trn(jpi-3,jj,jk,jn)
        zc = (ztrx(jpi-3,jj,jk) - ztlx(jpi-3,jj,jk))*
     $                   (trn(jpi-3,jj,jk,jn)-0.5*(ztrx(jpi-3,jj,jk) 
     $                   +ztlx(jpi-3,jj,jk)))
        zd = ((ztrx(jpi-3,jj,jk)-ztlx(jpi-3,jj,jk))**2)/6.
        zign  = sign(1.,za*zb)
        ztlx(jpi-3,jj,jk) = 0.5*(1+zign)
     $                           *trn(jpi-3,jj,jk,jn)+0.5
     $                           *(1-zign)*ztlx(jpi-3,jj,jk)
        ztrx(jpi-3,jj,jk) = 0.5*(1+zign)
     $                           *trn(jpi-3,jj,jk,jn)+0.5
     $                           *(1-zign)*ztrx(jpi-3,jj,jk) 
        zign  = sign(1.,zc-zd)
        ztlx(jpi-3,jj,jk) = 0.5*(1+zign)*(3*trn(jpi-3,jj,jk,jn)
     $                           -2*ztrx(jpi-3,jj,jk))+0.5*(1-zign)
     $                           *ztlx(jpi-3,jj,jk)
        zign  = sign(1.,zc+zd)
        ztrx(jpi-3,jj,jk) = 0.5*(1-zign)*(3*trn(jpi-3,jj,jk,jn)
     $                           -2*ztlx(jpi-3,jj,jk))+0.5*(1+zign)
     $                           *ztrx(jpi-3,jj,jk)  
        za6x(jpi-3,jj,jk) = 6*(trn(jpi-3,jj,jk,jn)-0.5
     $                            *(ztrx(jpi-3,jj,jk)
     $                            +ztlx(jpi-3,jj,jk)))

        za =  ztrx(jpi-2,jj,jk) - trn(jpi-2,jj,jk,jn)
        zb =  ztlx(jpi-2,jj,jk) - trn(jpi-2,jj,jk,jn)
        zc = (ztrx(jpi-2,jj,jk) - ztlx(jpi-2,jj,jk))*
     $                   (trn(jpi-2,jj,jk,jn)-0.5*(ztrx(jpi-2,jj,jk) 
     $                   +ztlx(jpi-2,jj,jk)))
        zd = ((ztrx(jpi-2,jj,jk)-ztlx(jpi-2,jj,jk))**2)/6.
        zign  = sign(1.,za*zb)
        ztlx(jpi-2,jj,jk) = 0.5*(1+zign)
     $                           *trn(jpi-2,jj,jk,jn)+0.5
     $                           *(1-zign)*ztlx(jpi-2,jj,jk)
        ztrx(jpi-2,jj,jk) = 0.5*(1+zign)
     $                           *trn(jpi-2,jj,jk,jn)+0.5
     $                           *(1-zign)*ztrx(jpi-2,jj,jk) 
        zign  = sign(1.,zc-zd)
        ztlx(jpi-2,jj,jk) = 0.5*(1+zign)*(3*trn(jpi-2,jj,jk,jn)
     $                           -2*ztrx(jpi-2,jj,jk))+0.5*(1-zign)
     $                           *ztlx(jpi-2,jj,jk)
        zign  = sign(1.,zc+zd)
        ztrx(jpi-2,jj,jk) = 0.5*(1-zign)*(3*trn(jpi-2,jj,jk,jn)
     $                           -2*ztlx(jpi-2,jj,jk))+0.5*(1+zign)
     $                           *ztrx(jpi-2,jj,jk)  
        za6x(jpi-2,jj,jk) = 6*(trn(jpi-2,jj,jk,jn)-0.5
     $                            *(ztrx(jpi-2,jj,jk)
     $                            +ztlx(jpi-2,jj,jk)))

        za =  ztrx(jpi-1,jj,jk) - trn(jpi-1,jj,jk,jn)
        zb =  ztlx(jpi-1,jj,jk) - trn(jpi-1,jj,jk,jn)
        zc = (ztrx(jpi-1,jj,jk) - ztlx(jpi-1,jj,jk))*
     $                   (trn(jpi-1,jj,jk,jn)-0.5*(ztrx(jpi-1,jj,jk) 
     $                   +ztlx(jpi-1,jj,jk)))
        zd = ((ztrx(2,jj,jk)-ztlx(2,jj,jk))**2)/6.
        zign  = sign(1.,za*zb)
        ztlx(jpi-1,jj,jk) = 0.5*(1+zign)
     $                           *trn(jpi-1,jj,jk,jn)+0.5
     $                           *(1-zign)*ztlx(jpi-1,jj,jk)
        ztrx(jpi-1,jj,jk) = 0.5*(1+zign)
     $                           *trn(jpi-1,jj,jk,jn)+0.5
     $                           *(1-zign)*ztrx(jpi-1,jj,jk) 
        zign  = sign(1.,zc-zd)
        ztlx(jpi-1,jj,jk) = 0.5*(1+zign)*(3*trn(jpi-1,jj,jk,jn)
     $                           -2*ztrx(jpi-1,jj,jk))+0.5*(1-zign)
     $                           *ztlx(jpi-1,jj,jk)
        zign  = sign(1.,zc+zd)
        ztrx(jpi-1,jj,jk) = 0.5*(1-zign)*(3*trn(jpi-1,jj,jk,jn)
     $                           -2*ztlx(jpi-1,jj,jk))+0.5*(1+zign)
     $                           *ztrx(jpi-1,jj,jk)  
        za6x(jpi-1,jj,jk) = 6*(trn(jpi-1,jj,jk,jn)-0.5
     $                            *(ztrx(jpi-1,jj,jk)
     $                            +ztlx(jpi-1,jj,jk)))


           ENDDO
          ENDDO

          DO JK=1,JPK
            DO JJ=1,JPJ


          zigma = un(2,jj,jk)*zstep/e1u(2,jj)
              zeu   = e2u(2,jj)* fse3u(2,jj,jk)* un(2,jj,jk)
          zign  = sign(1.,un(2,jj,jk))
          zkx(2,jj,jk) = tmask(2,jj,jk)*(tmask(jpi-2,jj,jk)
     $                        *tmask(jpi-1,jj,jk)
     $                        *tmask(3,jj,jk)*tmask(4,jj,jk)
     $                        *tmask(5,jj,jk)*(0.5*(1+zign)
     $                        *zeu*(ztrx(2,jj,jk)-0.5
     $                        *zigma*(ztrx(2,jj,jk)-ztlx(2,jj,jk)
     $                        -(1-2.*zigma/3.)*za6x(2,jj,jk)))+0.5
     $                        *(1-zign)*zeu*(ztlx(2,jj,jk)-0.5*zigma
     $                        *(ztrx(3,jj,jk)-ztlx(3,jj,jk)
     $                        +(1+2.*zigma/3.)*za6x(3,jj,jk))))
     $                        +((1-tmask(jpi-2,jj,jk))
     $                        +(1-tmask(jpi-1,jj,jk))
     $                        +(1-tmask(3,jj,jk))
     $                        +(1-tmask(4,jj,jk))
     $                        +(1-tmask(5,jj,jk)))
     $                        *zeu*(0.5*(1+zign)*trn(2,jj,jk,jn)+0.5
     $                        *(1-zign)*trn(3,jj,jk,jn)
     $                        *tmask(3,jj,jk)))

          zigma = un(3,jj,jk)*zstep/e1u(3,jj)
              zeu   = e2u(3,jj)* fse3u(3,jj,jk)* un(3,jj,jk)
          zign  = sign(1.,un(3,jj,jk))
          zkx(3,jj,jk) = tmask(3,jj,jk)*(tmask(jpi-1,jj,jk)
     $                        *tmask(2,jj,jk)
     $                        *tmask(4,jj,jk)*tmask(5,jj,jk)
     $                        *tmask(6,jj,jk)*(0.5*(1+zign)
     $                        *zeu*(ztrx(3,jj,jk)-0.5
     $                        *zigma*(ztrx(3,jj,jk)-ztlx(3,jj,jk)
     $                        -(1-2.*zigma/3.)*za6x(3,jj,jk)))+0.5
     $                        *(1-zign)*zeu*(ztlx(3,jj,jk)-0.5*zigma
     $                        *(ztrx(4,jj,jk)-ztlx(4,jj,jk)
     $                        +(1+2.*zigma/3.)*za6x(4,jj,jk))))
     $                        +((1-tmask(jpi-1,jj,jk))
     $                        +(1-tmask(2,jj,jk))
     $                        +(1-tmask(4,jj,jk))
     $                        +(1-tmask(5,jj,jk))
     $                        +(1-tmask(6,jj,jk)))
     $                        *zeu*(0.5*(1+zign)*trn(3,jj,jk,jn)+0.5
     $                        *(1-zign)*trn(4,jj,jk,jn)
     $                        *tmask(4,jj,jk)))

          zigma = un(jpi-3,jj,jk)*zstep/e1u(jpi-3,jj)
              zeu   = e2u(jpi-3,jj)*fse3u(jpi-3,jj,jk)*un(jpi-3,jj,jk)
          zign  = sign(1.,un(jpi-3,jj,jk))
          zkx(jpi-3,jj,jk) = tmask(jpi-3,jj,jk)*(tmask(jpi-5,jj,jk)
     $                        *tmask(jpi-4,jj,jk)*tmask(jpi-2,jj,jk)
     $                        *tmask(jpi-1,jj,jk)
     $                        *tmask(2,jj,jk)*(0.5*(1+zign)
     $                        *zeu*(ztrx(jpi-3,jj,jk)-0.5
     $                        *zigma*(ztrx(jpi-3,jj,jk)
     $                        -ztlx(jpi-3,jj,jk)
     $                        -(1-2.*zigma/3.)*za6x(jpi-3,jj,jk)))+0.5
     $                        *(1-zign)*zeu*(ztlx(jpi-3,jj,jk)-0.5
     $                        *zigma
     $                        *(ztrx(jpi-2,jj,jk)-ztlx(jpi-2,jj,jk)
     $                        +(1+2.*zigma/3.)*za6x(2,jj,jk))))
     $                        +((1-tmask(jpi-5,jj,jk))
     $                        +(1-tmask(jpi-4,jj,jk))
     $                        +(1-tmask(jpi-2,jj,jk))
     $                        +(1-tmask(jpi-1,jj,jk))
     $                        +(1-tmask(2,jj,jk)))
     $                        *zeu*(0.5*(1+zign)*trn(jpi-3,jj,jk,jn)
     $                        +0.5*(1-zign)*trn(jpi-2,jj,jk,jn)
     $                        *tmask(jpi-2,jj,jk)))

          zigma = un(jpi-2,jj,jk)*zstep/e1u(jpi-2,jj)
              zeu   = e2u(jpi-2,jj)*fse3u(jpi-2,jj,jk)*un(jpi-2,jj,jk)
          zign  = sign(1.,un(jpi-2,jj,jk))
          zkx(jpi-2,jj,jk) = tmask(jpi-2,jj,jk)*(tmask(jpi-4,jj,jk)
     $                        *tmask(jpi-3,jj,jk)*tmask(jpi-1,jj,jk)
     $                        *tmask(2,jj,jk)
     $                        *tmask(3,jj,jk)*(0.5*(1+zign)
     $                        *zeu*(ztrx(jpi-2,jj,jk)-0.5
     $                        *zigma*(ztrx(jpi-2,jj,jk)
     $                        -ztlx(jpi-2,jj,jk)
     $                        -(1-2.*zigma/3.)*za6x(jpi-2,jj,jk)))+0.5
     $                        *(1-zign)*zeu*(ztlx(jpi-2,jj,jk)-0.5
     $                        *zigma
     $                        *(ztrx(jpi-1,jj,jk)-ztlx(jpi-1,jj,jk)
     $                        +(1+2.*zigma/3.)*za6x(jpi-1,jj,jk))))
     $                        +((1-tmask(jpi-4,jj,jk))
     $                        +(1-tmask(jpi-3,jj,jk))
     $                        +(1-tmask(jpi-1,jj,jk))
     $                        +(1-tmask(2,jj,jk))
     $                        +(1-tmask(3,jj,jk)))
     $                        *zeu*(0.5*(1+zign)*trn(jpi-2,jj,jk,jn)
     $                        +0.5*(1-zign)*trn(jpi-1,jj,jk,jn)
     $                        *tmask(jpi-1,jj,jk)))

          zigma = un(jpi-1,jj,jk)*zstep/e1u(jpi-1,jj)
              zeu   = e2u(jpi-1,jj)*fse3u(jpi-1,jj,jk)*un(jpi-1,jj,jk)
          zign  = sign(1.,un(jpi-1,jj,jk))
          zkx(jpi-1,jj,jk) = tmask(jpi-1,jj,jk)*(tmask(jpi-3,jj,jk)
     $                        *tmask(jpi-2,jj,jk)*tmask(2,jj,jk)
     $                        *tmask(3,jj,jk)
     $                        *tmask(4,jj,jk)*(0.5*(1+zign)
     $                        *zeu*(ztrx(jpi-1,jj,jk)-0.5
     $                        *zigma*(ztrx(jpi-1,jj,jk)
     $                        -ztlx(jpi-1,jj,jk)
     $                        -(1-2.*zigma/3.)*za6x(jpi-1,jj,jk)))+0.5
     $                        *(1-zign)*zeu*(ztlx(jpi-1,jj,jk)-0.5
     $                        *zigma
     $                        *(ztrx(2,jj,jk)-ztlx(2,jj,jk)
     $                        +(1+2.*zigma/3.)*za6x(2,jj,jk))))
     $                        +((1-tmask(jpi-3,jj,jk))
     $                        +(1-tmask(jpi-2,jj,jk))
     $                        +(1-tmask(2,jj,jk))
     $                        +(1-tmask(3,jj,jk))
     $                        +(1-tmask(4,jj,jk)))
     $                        *zeu*(0.5*(1+zign)*trn(jpi-1,jj,jk,jn)
     $                        +0.5*(1-zign)*trn(2,jj,jk,jn)
     $                        *tmask(2,jj,jk)))


           ENDDO
          ENDDO

          DO JK=1,JPK
            DO JJ=1,JPJ
              zkx( 1 ,JJ,JK)=zkx(jpi-1,JJ,JK)
          zkx(JPI,JJ,JK)=zkx(  2  ,JJ,JK)
              zky( 1 ,JJ,JK)=zky(jpi-1,JJ,JK)
          zky(JPI,JJ,JK)=zky(  2  ,JJ,JK)
              zkz( 1 ,JJ,JK)=zkz(jpi-1,JJ,JK)    
          zkz(JPI,JJ,JK)=zkz(  2  ,JJ,JK)
            ENDDO
          ENDDO
        ENDIF

      DO JJ=1,JPJ
        DO JI=1,JPI
          zkx(JI,JJ,JPK)=0.
      zky(JI,JJ,JPK)=0.
      zkz(JI,JJ,JPK)=0.
        ENDDO
      ENDDO

    DO jj=1,jpj
       DO ji=1,jpi
          zkz(ji,jj,1)=0.
       ENDDO
    ENDDO


C
C 7. Resolution
C -------------
C

        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1
              zbt = e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk)
              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + tmask(ji,jj,jk)
     $                        * (zkx(ji-1,jj,jk)-zkx(ji,jj,jk)
     $                        +  zky(ji,jj-1,jk)-zky(ji,jj,jk)
     $                        +  zkz(ji,jj,jk+1)-zkz(ji,jj,jk))/zbt 

cek #ifdef diatrctrends
#if defined key_trc_diatrd
              trtrd(ji,jj,jk,1,jn) = trtrd(ji,jj,jk,1,jn) -
     $               ( zkx(ji,jj,jk) - zkx(ji-1,jj,jk) ) / zbt
              trtrd(ji,jj,jk,2,jn) = trtrd(ji,jj,jk,2,jn) -
     $               ( zky(ji,jj,jk) - zky(ji,jj-1,jk) ) / zbt
              trtrd(ji,jj,jk,3,jn) = trtrd(ji,jj,jk,3,jn) -
     $               ( zkz(ji,jj,jk) - zkz(ji,jj,jk+1) ) / zbt
#endif


            ENDDO
          ENDDO
        ENDDO  

C
C END of tracer loop
C ==================
C
 1000 CONTINUE
C
