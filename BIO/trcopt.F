
      SUBROUTINE trcopt(ktask,kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcopt
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC     computes the light propagation in the water column
CCC     and the euphotic layer depth
CCC
CCC
CC   METHOD :
CC   -------
CC
CC      multitasked on vertical slab (jj-loop)
CC      local par is computed in w layers using light propagation
CC      mean par in t layers are computed by integration
CC
CC
CC   INPUT :
CC   -----
CC      argument
CC              ktask           : task identificator
CC              kt              : time step
CC      COMMON
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC                                depths
CC            /comzdf/          : avt vertical eddy diffusivity
CC            /comqsr/          : solar radiation
CC            /comtsk/          : multitasking
CC            /cotopt/          : optical parameters
CC            /cotbio/          : biological parameters
CC
CC   OUTPUT :
CC   ------
CC      COMMON
CC            /cotopt/          : optical parameters
CC
CC   WORKSPACE :
CC   ---------
CC      local     zparr         : red compound of par
CC                zparg         : green compound of par
CC                zpar0m        : irradiance just below the surface
CC                zpar100       : irradiance at euphotic layer depth
CC                zkr           : total absorption coefficient in red
CC                zkg           : total absorption coefficient in green
CC                zpig          : total pigment
CC                imaske        : euphotic layer mask
CC                itabe         : euphotic layer last k index
CC
CC      COMMON
CC
CC   EXTERNAL :                   no
CC   --------

       USE myalloc        
       USE OPT_mem
             IMPLICIT NONE

#include <mpif.h>


CC local declarations
CC ==================
      INTEGER kt,ktask

#if defined key_passivetrc
#    if defined key_trc_nnpzddom || defined key_trc_npzd || key_trc_bfm 

      INTEGER ji,jj,jk,jn,in

      REAL(8) p_eps0
      REAL(8) p_epsChla01,p_epsChla02,p_epsChla03,p_epsChla04

CC+CC      REAL(8) zpar(jpi,jpk), p_par,p_epsR6
CC+CC      REAL(8) xEPS(jpi,jpk)
CC+CC      REAL(8) zpar0m(jpi),zpar100(jpi)

#if defined key_trc_bfm
      REAL(8) redf,rcchl
#endif 


CCC        DO ji = 1,jpi
CCC           WRITE(*,*) 'trcopt.F',ji,jj,kef(ji,jj)
CCC        ENDDO
CCC 10 11 2004  F79 cronometer-start
       trcoptparttime = MPI_WTIME()
C find Phytoplancton index - test CTRCNM

      in=0
      DO jn = 1,jptra
        IF ((ctrcnm(jn) .EQ. 'PHY') .OR.
     $      (ctrcnm(jn) .EQ. 'PHYTO') ) THEN
            in = jn
        END IF
      END DO

#if defined key_trc_bfm
CCC Paolo 15 10 2004 provisory value for the shadow effect
            redf  = 7.
            rcchl = 55.
            in = ppP1n
CCC
#endif
C
C vertical slab
C ===============
C
CCC F79 06 04 2005  3D  implementation
      DO 1000 jj = ktask,jpj,ncpu
CCCC      DO 1000 site_i = ktask,nsite,ncpu
CCCC      ji = site_10(site_i,1)
CCCC      jj = site_10(site_i,2)
C
C
C 1. determination of surface irradiance
C we provide also BFM set up
C   p_par(bfm)= 0.50
C   p_par(OPA)= 0.43
C    zpar0m(ji) = qsr(ji,jj)*p_par
C --------------------------------------
C
C
        DO ji = 1,jpi
          zpar0m(ji) = qsr(ji,jj)*0.50/0.217
CCC WARNING COnversion Einstein to Watt ??? E2W=0.217
          zpar100(ji) = zpar0m(ji)*0.01
          xpar(ji,jj,1) = zpar0m(ji)
          zpar(ji,1) =  zpar0m(1)
C          zparr(ji,1) = 0.5* zpar0m
C          zparg(ji,1) = 0.5* zpar0m(ji)
        END DO
C
C 2. determination of xpar
C ------------------------
C
C determination of local par in w levels
C bfm : p_eps0 = 0.12
C bfm : p_epsESS = 0.04e-3
C bfm : p_epsChla= 10.0e-3  ,  &  ! Chla-contribution to extinction
C bfm : ESS = Suspended sediment in g/m3 to be introduced
C bfm : p_epsR6 = =0.1e-3 (m2/mgC) # Detritus-contribution to extinction
C bfm : xEPS  =   p_eps0+ p_epsESS* ESS+ p_epsR6* R6c
C bfm : xEPS  =   xEPS+ p_epsChla* PhytoPlankton(i,iiI)
C bfm : xpar  =  EIR
CCC F79 for Med area we should put 0.08 !!!!!!
CCC        p_eps0 = 0.12
CCC        p_eps0 = 0.08
CCC        p_epsChla =10.0e-3  
CCC        p_epsR6 = 0.1e-3 
CCC New Values From J.C Blackford et al 
CCC Journal of Marine Ecosystem 52 (2004) 191-215
CCC We consider a Chl/C = 0.05 to convert
CCC carbon dependent shading to Chlorophyll 
CCC dependent shading
        p_eps0 = 0.04
        p_epsChla01 =10.e-3  
        p_epsChla02 =8.e-3  
        p_epsChla03 =8.e-3  
        p_epsChla04 =8.e-3  
CCC        p_epsR6 = 0.1e-3 

        DO jk = 2,jpk
          DO ji = 1,jpi
CCC      xEPS(ji,jk)= (p_eps0 -0.045*(ji+nimpp)/jpiglo)+
CCC F79 7 12 2005 We introduce variable Background
CCC from Crispi et al. journal of Mar. Sys. 33-34 (2002) 497-521
CCC see valuas in OPT_mem.F
            p_eps0 = kef(ji,jj) 
            xEPS(ji,jk)= p_eps0
CCC            xEPS(ji,jk)= p_eps0+p_epsR6* trn(ji,jj,jk - 1,ppR6c)
CCC      write(*,*) 'trcopt:ji,p_eps',ji+nimpp, xEPS(ji,jk)
CCC
CCC            xEPS(ji,jk)= xEPS(ji,jk)+
CCC     &                    p_epsChla01*trn(ji,jj,jk - 1,ppP1i)+
CCC     &                    p_epsChla02*trn(ji,jj,jk - 1,ppP2i)+
CCC     &                    p_epsChla03*trn(ji,jj,jk - 1,ppP3i)+
CCC     &                    p_epsChla04*trn(ji,jj,jk - 1,ppP4i)
CCC F79 this max to avoid denormalized number
            xEPS(ji,jk) = max(xEPS(ji,jk)
     &          ,1.D-15)

             xpar(ji,jj,jk)= xpar(ji,jj,jk-1) *
     &                    exp(-1. * xEPS(ji,jk-1)* e3t(jk-1)) 

CCC F79 this max to avoid denormalized number
            xpar(ji,jj,jk) = max(xpar(ji,jj,jk)
     &          ,1.D-15)

          END DO
        END DO
        DO jk = 1,jpk
          DO ji = 1,jpi
            xpar(ji,jj,jk) = xpar(ji,jj,jk) * exp(- xEPS(ji,jk)* 
     &      0.5D+00* e3t(jk) )
CCC F79 this max to avoid denormalized number
            xpar(ji,jj,jk) = max(xpar(ji,jj,jk)
     &          ,1.D-15)
          END DO
        END DO
C
C 4. determination of euphotic layer depth 
C ----------------------------------------
C
C imaske equal 1 in the euphotic layer, and 0 without
C

CCC Paolo 15 10 2004 we put a Tmask to cancell undefined points
CCC otherwise we obtain an out of bound

        DO jk = 1,jpk
          DO ji = 1,jpi
            imaske(ji,jk) = 0
CCC            IF (Tmask(ji,jj,jk) .EQ. 1) THEN
             IF (xpar(ji,jj,jk) .GE. zpar100(ji)) imaske(ji,jk) = 1
CCC            ENDIF
          END DO
        END DO

        DO ji = 1,jpi
          itabe(ji) = 0
        END DO

        DO jk = 1,jpk
          DO ji = 1,jpi
            itabe(ji) = itabe(ji) + imaske(ji,jk)
          END DO
        END DO
C
        DO ji = 1,jpi
          itabe(ji) = max(1,itabe(ji))
CCC          write(*,*) 'xze',xze(ji,jj),ji,jj
CCC          xze(ji,jj) = fsdepw(ji,jj,itabe(ji) + 1)
        END DO 


C END of slab
C ===========

 1000 CONTINUE

       trcoptparttime = MPI_WTIME() - trcoptparttime !F79 cronometer-stop
       trcopttottime = trcopttottime + trcoptparttime

#    else

C    No optical model

#    endif
#else

C no passive tracers

#endif

      END SUBROUTINE trcopt
