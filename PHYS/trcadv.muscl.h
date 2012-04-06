CCC                      trcadv.muscl.h
CCC                     ******************
CCC
CC   defined key : 'key_trc_muscl'
CC   ============
CC
CC  PURPOSE :
CC  ---------
CC     Compute the now trend due to total advection of tracers
CC     and add it to the general trend of tracer equations.
CC     MUSCL(Monotone Upstream-centered S?chemes for Conservation
CC     Laws) scheme is used.
CC    
CC     
CC
CC
CC   METHOD :
CC   -------
CC      this ROUTINE compute not exactly the advection but the
CC      transport term, i.e.  div(u*tra).
CC
CC
CC
CC
CC
CC   REFERENCES :                
CC   ----------                  
CC
CC   MODIFICATIONS:
CC   --------------
CC       original :  06-00 (A.Estublier)
CC       
CC----------------------------------------------------------------------
C
      INTEGER ji,jj,jk,jn
      REAL(8) ztrax(jpi,jpj,jpk),zatrax(jpi,jpj,jpk)
      REAL(8) ztray(jpi,jpj,jpk),zatray(jpi,jpj,jpk)
      REAL(8) ztraz(jpi,jpj,jpk),zatraz(jpi,jpj,jpk),zakz(jpi,jpj,jpk) 
      REAL(8) zkx(jpi,jpj,jpk),zky(jpi,jpj,jpk),zkz(jpi,jpj,jpk)
      REAL(8) zkx1(jpi,jpj,jpk)      
      REAL(8) zigma,zeu,zev,zew,zbt,zstep,zign,zmask(jpi,jpj,jpk)
      REAL(8) zbig
C
CC----------------------------------------------------------------------
CC statement functions
CC ===================

	zbig = 1.e+40

C tracer loop parallelized (macrotasking)
C =======================================
C
      DO jn = 1,jptra

        
C
C 1. Slopes computation
C ---------------------
C

C
C Initialization
C --------------
C

	zstep  = rdt*ndttrc

        ztrax  = 0
        zatrax = 0
        zkx    = 0
	zkx1 = 0

        ztray  = 0
        zatray = 0
        zky    = 0

        ztraz  = 0
        zatraz = 0
        zakz   = 0
        zkz    = 0


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
	    zkx(ji,jj,jk) = 0.5*(ztrax(ji,jj,jk)
     $                          +ztrax(ji-1,jj,jk))*zign

	    zign = 0.5*(sign(1.,ztray(ji,jj,jk)*ztray(ji,jj-1,jk))+1)
            zky(ji,jj,jk) = 0.5*(ztray(ji,jj,jk)
     $                          +ztray(ji,jj-1,jk))*zign
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
            zakz(ji,jj,jk) = 0.5*(ztraz(ji,jj,jk)
     $                          +ztraz(ji,jj,jk+1))*zign
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
     $                        sign(1.,zky(ji,jj,jk))* 
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
              zakz(ji,jj,jk) = tmask(ji,jj,jk+1)*tmask(ji,jj,jk-1)*
     $                        sign(1.,zakz(ji,jj,jk)) * 
     $                        min(abs(zakz(ji,jj,jk)),
     $                        2*zatraz(ji,jj,jk+1),
     $                        2*zatraz(ji,jj,jk))
            ENDDO
          ENDDO
        ENDDO        




C
C 3. Advection terms
C ------------------
C

C
C Fluxes in the i and j direction
C
        
        DO jk=1,jpkm1
          DO jj = 1,jpjm1      
            DO ji = 1,jpim1
	    
	      zigma = un(ji,jj,jk)*zstep/e1u(ji,jj)
              zeu   = e2u(ji,jj)* fse3u(ji,jj,jk)* un(ji,jj,jk)
	      zign  = sign(1.,un(ji,jj,jk))            

	      zkx(ji,jj,jk) = zeu*(trn(ji+int(0.5*(1-zign)),jj,jk,jn)
     $                            +0.5*zign*(1-zign*zigma)*
     $                            zkx(ji+int(0.5*(1-zign)),jj,jk))

              zigma = vn(ji,jj,jk)*zstep/e2v(ji,jj)
              zev   = e1v(ji,jj)* fse3v(ji,jj,jk)* vn(ji,jj,jk)
              zign  = sign(1.,vn(ji,jj,jk))            

	      zky(ji,jj,jk) = zev*(trn(ji,jj+int(0.5*(1-zign)),jk,jn)
     $                            +0.5*zign*(1-zign*zigma)*
     $                            zky(ji,jj+int(0.5*(1-zign)),jk))
 

            ENDDO
          ENDDO
        ENDDO        



C
C Fluxes in the k direction
C


	DO jj=2,jpjm1
	   DO ji=2,jpim1
	      zakz(ji,jj,1)=0.
	   ENDDO
	ENDDO

        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1

              zigma = wn(ji,jj,jk+1)*zstep/fse3w(ji,jj,jk+1)
              zew   = e1t(ji,jj)*e2t(ji,jj)*wn(ji,jj,jk+1)
	      zign  = sign(1.,wn(ji,jj,jk+1))

	      zkz(ji,jj,jk+1) = zew*(trn(ji,jj,jk+int(0.5*(1+zign)),jn)
     $                              +0.5*zign*(1-zign*zigma)*
     $                              zakz(ji,jj,jk+int(0.5*(1+zign))))

            ENDDO
          ENDDO
        ENDDO 
       

C
C Boundary conditions
C

	DO jj=1,jpj
	   DO ji=1,jpi
	      zkz(ji,jj,1)=0.
	   ENDDO
	ENDDO



       IF (jperio.eq.1) THEN

       DO jk=1,jpk
	DO jj=1,jpj

C
C Points jpi-2, jpi-1 et 2
C

        ztrax(2,jj,jk)  = trn(3,jj,jk,jn) - trn(2,jj,jk,jn)
        zatrax(2,jj,jk) = abs(ztrax(2,jj,jk))
        ztrax(jpi-2,jj,jk)  = trn(jpi-1,jj,jk,jn) - trn(jpi-2,jj,jk,jn)
        zatrax(jpi-2,jj,jk) = abs(ztrax(jpi-2,jj,jk))
        ztrax(jpi-1,jj,jk)  = trn(2,jj,jk,jn) - trn(jpi-1,jj,jk,jn)
        zatrax(jpi-1,jj,jk) = abs(ztrax(jpi-1,jj,jk))
	
	ENDDO
       ENDDO

       DO jk=1,jpk
        DO jj=1,jpj
 	
	zign = 0.5*(sign(1.,ztrax(2,jj,jk)*ztrax(jpi-1,jj,jk))+1)
	zkx1(2,jj,jk) = 0.5*(ztrax(2,jj,jk)
     $                     +ztrax(jpi-1,jj,jk))*zign
        zign = 0.5*(sign(1.,ztrax(3,jj,jk)*ztrax(2,jj,jk))+1)
	zkx1(3,jj,jk) = 0.5*(ztrax(3,jj,jk)
     $                     +ztrax(2,jj,jk))*zign
	zign = 0.5*(sign(1.,ztrax(jpi-2,jj,jk)*ztrax(jpi-3,jj,jk))+1)
	zkx1(jpi-2,jj,jk) = 0.5*(ztrax(jpi-2,jj,jk)
     $                         +ztrax(jpi-3,jj,jk))*zign
	zign = 0.5*(sign(1.,ztrax(jpi-1,jj,jk)*ztrax(jpi-2,jj,jk))+1)
	zkx1(jpi-1,jj,jk) = 0.5*(ztrax(jpi-1,jj,jk)
     $                         +ztrax(jpi-2,jj,jk))*zign

        ENDDO
       ENDDO
	
       DO jk=1,jpk
        DO jj=1,jpj

              zkx1(2,jj,jk) = tmask(3,jj,jk)*tmask(jpi-1,jj,jk)
     $                        *sign(1.,zkx1(2,jj,jk)) * 
     $                        min(abs(zkx1(2,jj,jk)),
     $                        2*zatrax(jpi-1,jj,jk),
     $                        2*zatrax(2,jj,jk))

              zkx1(3,jj,jk) = tmask(4,jj,jk)*tmask(2,jj,jk)
     $                        *sign(1.,zkx1(3,jj,jk)) * 
     $                        min(abs(zkx1(3,jj,jk)),
     $                        2*zatrax(2,jj,jk),
     $                        2*zatrax(3,jj,jk))

              zkx1(jpi-2,jj,jk) = tmask(jpi-1,jj,jk)*tmask(jpi-3,jj,jk)
     $                        *sign(1.,zkx1(jpi-2,jj,jk)) * 
     $                        min(abs(zkx1(jpi-2,jj,jk)),
     $                        2*zatrax(jpi-3,jj,jk),
     $                        2*zatrax(jpi-2,jj,jk))

              zkx1(jpi-1,jj,jk) = tmask(jpi-2,jj,jk)*tmask(2,jj,jk)
     $                        *sign(1.,zkx1(jpi-1,jj,jk)) * 
     $                        min(abs(zkx1(jpi-1,jj,jk)),
     $                        2*zatrax(jpi-2,jj,jk),
     $                        2*zatrax(jpi-1,jj,jk))

	ENDDO
       ENDDO
	
       DO jk=1,jpk
        DO jj=1,jpj
	
	zigma = un(2,jj,jk)*zstep/e1u(2,jj)
        zeu   = e2u(2,jj)* fse3u(2,jj,jk)* un(2,jj,jk)
	zign  = sign(1.,un(2,jj,jk))            

	zkx(2,jj,jk) = zeu*(trn(2+int(0.5*(1-zign)),jj,jk,jn)
     $                          +0.5*zign*(1-zign*zigma)*
     $                          zkx1(2+int(0.5*(1-zign)),jj,jk))


	zigma = un(jpi-2,jj,jk)*zstep/e1u(jpi-2,jj)
        zeu   = e2u(jpi-2,jj)* fse3u(jpi-2,jj,jk)* un(jpi-2,jj,jk)
        zign  = sign(1.,un(jpi-2,jj,jk))            

        zkx(jpi-2,jj,jk) = zeu*(trn(jpi-2+int(0.5*(1-zign)),jj,jk,jn)
     $                          +0.5*zign*(1-zign*zigma)*
     $                          zkx1(jpi-2+int(0.5*(1-zign)),jj,jk))


	zigma = un(jpi-1,jj,jk)*zstep/e1u(jpi-1,jj)
        zeu   = e2u(jpi-1,jj)* fse3u(jpi-1,jj,jk)* un(jpi-1,jj,jk)
        zign  = sign(1.,un(jpi-1,jj,jk))            

	zkx(jpi-1,jj,jk) = zeu*(0.5*(1+zign)*trn(jpi-1,jj,jk,jn)
     $                         +0.5*(1-zign)*trn(2,jj,jk,jn)
     $                         +0.5*zign*(1-zign*zigma)*
     $                         (0.5*(1+zign)*zkx1(jpi-1,jj,jk)
     $                         +0.5*(1-zign)*zkx1(2,jj,jk)))
        ENDDO
       ENDDO
C
C Points 1 and jpi
C 
	DO jk=1,jpk
         DO jj=1,jpj
	       zkx(1,jj,jk)=zkx(jpi-1,jj,jk)
	       zkx(jpi,jj,jk)=zkx(2,jj,jk)
	       zky(1,jj,jk)=zky(jpi-1,jj,jk)
	       zky(jpi,jj,jk)=zky(2,jj,jk)
	       zkz(1,jj,jk)=zkz(jpi-1,jj,jk)
	       zkz(jpi,jj,jk)=zkz(2,jj,jk)
         ENDDO
        ENDDO
        
	ELSE
          DO jk=1,jpk
            DO jj=1,jpj
              zkx(1  ,jj,jk)=0.
              zkx(jpi,jj,jk)=0.
              zky(1  ,jj,jk)=0.
              zky(jpi,jj,jk)=0.
              zkz(1  ,jj,jk)=0.
              zkz(jpi,jj,jk)=0.
            ENDDO
          ENDDO
      ENDIF
C
      DO JJ=1,JPJ
        DO JI=1,JPI
          zkx(JI,JJ,JPK)=0.
	  zky(JI,JJ,JPK)=0.
	  zkz(JI,JJ,JPK)=0.
        ENDDO
      ENDDO

C
C 5. Resolution
C -------------
C

        DO jk=1,jpkm1
          DO jj = 2,jpjm1      
            DO ji = 2,jpim1

              zbt = e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk)
              tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn)*tmask(ji,jj,jk)
     $                        + (zkx(ji-1,jj,jk)-zkx(ji,jj,jk)
     $                        +  zky(ji,jj-1,jk)-zky(ji,jj,jk)
     $                        +  zkz(ji,jj,jk+1)-zkz(ji,jj,jk))/zbt 
#if defined key_trc_diatrd 
              trtrd(ji,jj,jk,jn,1) =
     $               (zkx(ji-1,jj,jk)-zkx(ji,jj,jk))/zbt
              trtrd(ji,jj,jk,jn,2) =
     $               (zky(ji,jj-1,jk)-zky(ji,jj,jk))/zbt
              trtrd(ji,jj,jk,jn,3) =
     $               (zkz(ji,jj,jk+1)-zkz(ji,jj,jk))/zbt
#endif
            ENDDO
          ENDDO
        ENDDO        

#      if defined key_trc_diatrd
#         ifdef key_mpp
        CALL mpplnk( trtrd(1,1,1,jn,1), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,jn,2), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,jn,3), 1, 1 )
#         else      
        CALL lbc( trtrd(1,1,1,jn,1), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,jn,2), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,jn,3), 1, 1, 1, 1, jpk, 1 )
#         endif
#      endif

C
C END of tracer loop
C ==================
C
       END DO

