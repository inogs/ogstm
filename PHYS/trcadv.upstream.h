CCC $Id: trcadv.upstream.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
CCC
CCC                      trcadv.smolar.h
CCC                     ******************
CCC
CC   defined key : 'key_trc_smolar'
CC   ============
CC
CC  PURPOSE :
CC  ---------
CC     compute the now trend due to the advection of passive tracers
CC     and add it to the general trend of tracer equations:
CC     THEN computes both horizontal and
CC      vertical advection of tracer trn
CC
CC
CC   METHOD :
CC   -------
CC      this ROUTINE compute not exactly the advection but the
CC      transport term, i.e.  div(u*tra).
CC
CC      smolarkevisz scheme
CC      *******************
CC
CC      computes the now horizontal and vertical advection with the
CC                       ----------     --------
CC      complete 3d method.
CC
CC      cf reference
CC
CC      note: - sc is an empirical factor to be used with care
CC            - this advection scheme needs an euler-forward time scheme
CC
CC   remarks :
CC   -------
CC      IF 'key_trc_diatrd' key is activated, the
CC      advection trend is saved for futher diagnostics.
CC
CC      multitasked on tracer (jn-loop)
CC
CC      COMMON
CC            /comcoo/          : orthogonal curvilinear coordinates
CC                                and scale factors
CC            /cottrp/          : present and next fields for passive
CC                              : tracers
CC            /comtsk/          : multitasking
CC            /cot3ad/          : smolar advection
CC
CC   OUTPUT :
CC   ------
CC      COMMON
CC            /cottrp/ tra      : general tracer trend increased by the
CC                                now horizontal tracer advection trend
CC            /cottrd/ trtrd    : now horizontal tracer advection trend
CC                                (IF 'key_trc_diatrd' key is activated)
CC
CC
CC   EXTERNAL :                  upstream, forkjoin, mpplnk, lbc
CC   --------
CC
CC   REFERENCES :                piotr k. smolarkiewicz, 1983,
CC   ----------                  "a simple positive definit advection
CC                               scheme with small IMPLICIT diffusion"
CC                               monthly weather review, pp 479-486
CC
CC   MODIFICATIONS:
CC   --------------
CC      original :  87-06 (pa-dl)
CC      additions : 91-11 (G. Madec)
CC      additions : 94-08 (a. czaja)
CC      modifications : 95-09 (M. Levy) passive tracers
CC      modifications : 98-03 (M.A. Foujols) lateral boundary conditions
CC      modifications : 99-02 (M.A. Foujols) lbc in conjonction with ORCA
CC      modifications : 00-05 (MA Foujols) add lbc for tracer trends
CC      modifications : 00-10 (MA Foujols and E.Kestenare) INCLUDE instead of routine
CC      modifications : 01-05 (E.Kestenare) fix bug in trtrd indexes
CC----------------------------------------------------------------------

      INTEGER ji,jj,jk,jt,jn
      REAL(8) zbtr,zdt
      REAL(8) zgm,zgz
CC-CC      REAL(8) zti(jpi,jpj,jpk),ztj(jpi,jpj,jpk)
CC-CC      REAL(8) zaa(jpi,jpj,jpk),zbb(jpi,jpj,jpk),zcc(jpi,jpj,jpk)
CC-CC      REAL(8) zx(jpi,jpj,jpk),zy(jpi,jpj,jpk),zz(jpi,jpj,jpk)
CC-CC      REAL(8) zbuf(jpi,jpj,jpk)
CC-CC      REAL(8) zkx(jpi,jpj,jpk),zky(jpi,jpj,jpk),zkz(jpi,jpj,jpk)
CC-CC      SAVE zti, ztj, zaa, zbb, zcc, zx, zy, zz, zbuf, zkx, zky, zkz
CC----------------------------------------------------------------------
CC statement functions
CC ===================

!                                   #include "stafun.h"

CCC 13012004 Opt new AIX52
CC-CC      integer goodpoints, allpoints, tpoints
CC-CC      save goodpoints, allpoints, tpoints 
CC-CC      data goodpoints /0/
CC-CC      data allpoints /0/
CC-CC      data tpoints /0/
CC-CC      integer advmask(jpi,jpj,jpk)
CC-CC      save advmask
CC-CC      real(8) inv_eu(jpi,jpj,jpk), inv_ev(jpi,jpj,jpk), inv_et(jpi,jpj,jpk)
CC-CC      save inv_eu, inv_ev, inv_et
CC-CC      real(8) big_fact_zaa(jpi,jpj,jpk), big_fact_zbb(jpi,jpj,jpk), big_fact_zcc(jpi,jpj,jpk)
CC-CC      save big_fact_zaa, big_fact_zbb, big_fact_zcc
CC-CC      real(8) zbtr_arr(jpi,jpj,jpk)
CC-CC      save zbtr_arr
CC-CC      integer jarr(3, jpi*jpj*jpk), dimen_jarr
CC-CC      save jarr, dimen_jarr
CC-CC      data dimen_jarr /0/
CC-CC      integer ju
CC-CC      integer jarr1(3, jpi*jpj*jpk), dimen_jarr1
CC-CC      save jarr1, dimen_jarr1
CC-CC      data dimen_jarr1 /0/
CC-CC      integer jarr2(3, jpi*jpj*jpk), dimen_jarr2
CC-CC      save jarr2, dimen_jarr2
CC-CC      data dimen_jarr2 /0/
CC-CC      integer jarr3(3, jpi*jpj*jpk), dimen_jarr3
CC-CC      save jarr3, dimen_jarr3
CC-CC      data dimen_jarr3 /0/
CC-CC      integer jarrt(3, jpi*jpj*jpk), dimen_jarrt
CC-CC      save jarrt, dimen_jarrt
CC-CC      data dimen_jarrt /0/
CC-CC      integer myji, myjj, myjk, locsum
CC-CC      integer jilef, jjlef, jklef, jirig, jjrig, jkrig

CCC 13012004 Opt new AIX52


CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (11/96)
CCC---------------------------------------------------------------------
C
CCC 10 11 2004  F79 cronometer-start

CCC       call mppsync()


      if(allpoints .EQ. 0) then
         zti  = 0.
         ztj  = 0.
         zaa  = 0.
         zbb  = 0.
         zcc  = 0.
         zx   = 0.
         zy   = 0. 
         zz   = 0.
         zbuf = 0.
         zkx  = 0.
         zky  = 0.
         zkz  = 0.
         inv_eu = 0.
         inv_ev = 0.
         inv_et = 0.
         big_fact_zaa = 0.
         big_fact_zbb = 0.
         big_fact_zcc = 0.
         zbtr_arr = 0.
         jarr  = 0
         jarr1 = 0
         jarr2 = 0
         jarr3 = 0
         jarrt = 0
         write(*,*) "Storing good points ..." 
         DO jk = 1,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  
                  zbtr_arr(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
               END DO
            END DO
         END DO

         DO jk = 1,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1

                  inv_eu(ji,jj,jk) = 1./(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) )
C                  inv_etot(1,ji,jj,jk) = inv_eu(ji,jj,jk)
                  inv_ev(ji,jj,jk) = 1./(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) )
C                  inv_etot(2,ji,jj,jk) = inv_ev(ji,jj,jk)
               END DO
            END DO
         END DO

         DO jk = 2,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  
                  inv_et(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk) )
C                  inv_etot(3,ji,jj,jk) = inv_et(ji,jj,jk)
               END DO
            END DO
         END DO

         allpoints = 0
         DO jk = 1,jpk
            DO jj = 1,jpj
               DO ji = 1,jpi
                  allpoints = allpoints + 1
                  if(tmask(ji,jj,jk) .NE. 0) then
                     tpoints = tpoints + 1
                  endif
               END DO
            END DO
         END DO

         goodpoints = 0
         DO jk = 1,jpk
            DO jj = 1,jpj
               DO ji = 1,jpi
                  jklef = -1
                  jjlef = -1
                  jilef = -1
                  jkrig = +1
                  jjrig = +1
                  jirig = +1
                  if(jk .EQ. 1)   jklef = 0
                  if(jj .EQ. 1)   jjlef = 0
                  if(ji .EQ. 1)   jilef = 0
                  if(jk .EQ. jpk) jkrig = 0
                  if(jj .EQ. jpj) jjrig = 0
                  if(ji .EQ. jpi) jirig = 0
                  locsum = 0
                  DO myjk=jk+jklef, jk+jkrig
                     DO myjj=jj+jjlef, jj+jjrig
                        DO myji=ji+jilef, ji+jirig
                           locsum = locsum + tmask(myji, myjj, myjk)
                        END DO
                     END DO
                  END DO
                  if(locsum .NE. 0) then
                     goodpoints = goodpoints + 1
                     advmask(ji,jj,jk) = 1
                  else
                     advmask(ji,jj,jk) = 0
                  endif
               END DO
            END DO
         END DO
         DO jk = 1,jpk
            DO jj = 2,jpjm1
               DO  ji = 2,jpim1
                  if(advmask(ji,jj,jk) .NE. 0) then
                     zbtr_arr(ji,jj,jk) = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
                     dimen_jarr = dimen_jarr + 1
                     jarr(1,dimen_jarr) = ji
                     jarr(2,dimen_jarr) = jj
                     jarr(3,dimen_jarr) = jk
                  else
                     zbtr_arr(ji,jj,jk) = 0.
                  endif
               END DO
            END DO
         END DO
         DO jk = 2,jpk
            DO jj = 2,jpjm1
               DO  ji = 2,jpim1
                  if(advmask(ji,jj,jk) .NE. 0) then
                     dimen_jarr1 = dimen_jarr1 + 1
                     jarr1(1,dimen_jarr1) = ji
                     jarr1(2,dimen_jarr1) = jj
                     jarr1(3,dimen_jarr1) = jk
                  endif
               END DO
            END DO
         END DO
         DO jk = 2,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  if(advmask(ji,jj,jk) .NE. 0) then
                     dimen_jarr2 = dimen_jarr2 + 1
                     jarr2(1,dimen_jarr2) = ji
                     jarr2(2,dimen_jarr2) = jj
                     jarr2(3,dimen_jarr2) = jk
                  endif
               END DO
            END DO
         END DO
         DO jk = 1,jpkm1
            DO jj = 2,jpjm1
               DO ji = 2,jpim1
                  if(advmask(ji,jj,jk) .NE. 0) then
                     dimen_jarr3 = dimen_jarr3 + 1
                     jarr3(1,dimen_jarr3) = ji
                     jarr3(2,dimen_jarr3) = jj
                     jarr3(3,dimen_jarr3) = jk
                  endif
               END DO
            END DO
         END DO
         DO jk = 1,jpk
            DO jj = 1,jpj
               DO ji = 1,jpi
                  if(tmask(ji,jj,jk) .NE. 0) then  
                     dimen_jarrt = dimen_jarrt + 1
                     jarrt(1,dimen_jarrt) = ji
                     jarrt(2,dimen_jarrt) = jj
                     jarrt(3,dimen_jarrt) = jk
                  endif
               END DO
            END DO
         END DO
      endif

      trcadvparttime = MPI_WTIME()

      zdt = rdt*float(ndttrc)

CCC
       DO jk = 1,jpk
          DO jj = 1,jpj
             DO ji = 1,jpi
                zaa(ji,jj,jk) = e2u(ji,jj)*fse3u(ji,jj,jk) * un(ji,jj,jk)
                zbb(ji,jj,jk) = e1v(ji,jj)*fse3v(ji,jj,jk) * vn(ji,jj,jk)
                zcc(ji,jj,jk) = e1t(ji,jj)*e2t(ji,jj)      * wn(ji,jj,jk)
                ztj(ji,jj,jk) = 0.
                zx(ji,jj,jk) = 0.
                zy(ji,jj,jk) = 0.
                zz(ji,jj,jk) = 0.
                big_fact_zaa(ji,jj,jk) = ( abs(zaa(ji,jj,jk)) - zdt
     $              *zaa(ji,jj,jk)**2
CCC     $              /(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) ) )
     $              *inv_eu(ji,jj,jk) )
                big_fact_zbb(ji,jj,jk) = ( abs(zbb(ji,jj,jk)) - zdt
     $              *zbb(ji,jj,jk)**2
CCC     $              /(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) ) )
     $              *inv_ev(ji,jj,jk) )
                big_fact_zcc(ji,jj,jk) = ( abs(zcc(ji,jj,jk)) - zdt
     $              *zcc(ji,jj,jk)**2
CCC     $              /(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk) ) )
     $              *inv_et(ji,jj,jk) )
             END DO
          END DO
       END DO

C

C
C tracer loop parallelized (macrotasking)
C =======================================
C
      DO 1000 jn = ktask,jptra,ncpu
C
C 1. tracer flux in the 3 directions
C ----------------------------------
C
C 1.1 mass flux at u v and t-points and initialization
C

        DO jk = 1,jpk
          DO jj = 1,jpj
            DO ji = 1,jpi
CCC              zaa(ji,jj,jk) = e2u(ji,jj)*fse3u(ji,jj,jk) * un(ji,jj,jk)
CCC              zbb(ji,jj,jk) = e1v(ji,jj)*fse3v(ji,jj,jk) * vn(ji,jj,jk)
CCC              zcc(ji,jj,jk) = e1t(ji,jj)*e2t(ji,jj)      * wn(ji,jj,jk)
CCC              zbuf(ji,jj,jk) = 0.
CCC              ztj(ji,jj,jk) = 0.
CCC              zx(ji,jj,jk) = 0.
CCC              zy(ji,jj,jk) = 0.
CCC              zz(ji,jj,jk) = 0.
CCC              zti(ji,jj,jk) = trn(ji,jj,jk,jn)
#if defined key_trc_diatrd
              trtrd(ji,jj,jk,jn,1) = 0.
              trtrd(ji,jj,jk,jn,2) = 0.
              trtrd(ji,jj,jk,jn,3) = 0.
#endif
            END DO
          END DO
       END DO
C
C 1.2 calcul of intermediate field with an upstream advection scheme
C     and mass fluxes calculated above
C
C calcul of tracer flux in the i and j direction
C
       DO jk = 1,jpk
          DO jj=1,jpj
            zkx(1,jj,jk)=0.
          END DO
          DO jj=1,jpj
            zkx(jpi,jj,jk)=0.
          END DO

          DO ji=1,jpi
            zky(ji,1,jk)=0.
          END DO
          DO ji=1,jpi
            zky(ji,jpj,jk)=0.
          END DO
       END DO

C
C calcul of tracer flux in the k direction
C
        DO jj = 1,jpj
          DO ji = 1,jpi
            zkz(ji,jj,1) = 0.
          END DO
        END DO

        jk = 1
        DO jj = 2,jpjm1
           DO ji = 2,jpim1
              zkx(ji,jj,jk) = fsx(
     $            trn(ji,jj,jk,jn),trn(ji + 1,jj,jk,jn),zaa(ji,jj,jk))
              zky(ji,jj,jk) = fsy(
     $            trn(ji,jj,jk,jn),trn(ji,jj + 1,jk,jn),zbb(ji,jj,jk))
           END DO
        END DO

        DO jk = 2,jpk         
          jj=1
          DO ji = 1,jpi
             zkz(ji,jj,jk) = fsz(
     $            trn(ji,jj,jk,jn),trn(ji,jj,jk - 1,jn),zcc(ji,jj,jk))
          END DO
          jj= jpj
          DO ji = 1,jpi
             zkz(ji,jj,jk) = fsz(
     $            trn(ji,jj,jk,jn),trn(ji,jj,jk - 1,jn),zcc(ji,jj,jk))
          END DO
          DO jj = 2,jpjm1
             ji=1
             zkz(ji,jj,jk) = fsz(
     $            trn(ji,jj,jk,jn),trn(ji,jj,jk - 1,jn),zcc(ji,jj,jk))
             ji=jpi
             zkz(ji,jj,jk) = fsz(
     $            trn(ji,jj,jk,jn),trn(ji,jj,jk - 1,jn),zcc(ji,jj,jk))
          END DO
        END DO

CCC        DO jk = 2,jpk
CCC          DO jj = 2,jpjm1
CCC            DO ji = 2,jpim1

        DO ju=1, dimen_jarr1

           ji = jarr1(1, ju)
           jj = jarr1(2, ju)
           jk = jarr1(3, ju)

              zkx(ji,jj,jk) = fsx(
     $            trn(ji,jj,jk,jn),trn(ji + 1,jj,jk,jn),zaa(ji,jj,jk))
              zky(ji,jj,jk) = fsy(
     $            trn(ji,jj,jk,jn),trn(ji,jj + 1,jk,jn),zbb(ji,jj,jk))
              zkz(ji,jj,jk) = fsz(
     $            trn(ji,jj,jk,jn),trn(ji,jj,jk - 1,jn),zcc(ji,jj,jk))    

        END DO

CCC            END DO
CCC          END DO
CCC        END DO
C     
C
C ... Lateral boundary conditions on zk[xy]
#ifdef key_mpp
C
C   ... Mpp : export boundary values to neighboring processors
C

        CALL mpplnk( zkx, 1, 1 )
        CALL mpplnk( zky, 1, 1 )

#else
C
C   ... T-point, 3D array, full local arrays zk[xy] are initialised
C       
C
        CALL lbc( zkx, 1, 1, 1, 1, jpk, 1 )
        CALL lbc( zky, 1, 1, 1, 1, jpk, 1 )
C
#endif
C
C
C 2. calcul of after field using an upstream advection scheme
C -----------------------------------------------------------
C
CCC        DO jk = 1,jpkm1
CCC          DO jj = 2,jpjm1
CCC            DO ji = 2,jpim1

        DO ju=1, dimen_jarr3

           ji = jarr3(1, ju)
           jj = jarr3(2, ju)
           jk = jarr3(3, ju)

CCC             if(zbtr_arr(ji,jj,jk) .NE. 0) then
CCC                 zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
           zbtr = zbtr_arr(ji,jj,jk)
           tra(ji,jj,jk,jn) = -zbtr*
     $          ( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk)
     $          + zky(ji,jj,jk) - zky(ji,jj - 1,jk)
     $          + zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
#if defined key_trc_diatrd
           trtrd(ji,jj,jk,jn,1) = trtrd(ji,jj,jk,jn,1) -
     $          zbtr*( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk) )
           trtrd(ji,jj,jk,jn,2) = trtrd(ji,jj,jk,jn,2) -
     $          zbtr*( zky(ji,jj,jk) - zky(ji,jj - 1,jk) )
           trtrd(ji,jj,jk,jn,3) = trtrd(ji,jj,jk,jn,3) -
     $          zbtr*( zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
#endif

        END DO

C
C
C 4.0 convert the transport trend into advection trend
C ----------------------------------------------------
C
#if defined key_trc_diatrd
CCC        DO jk = 1,jpk
CCC          DO jj = 2,jpjm1
CCC            DO  ji = 2,jpim1

         DO ju=1, dimen_jarr

           ji = jarr(1, ju)
           jj = jarr(2, ju)
           jk = jarr(3, ju)

CCC               if(zbtr_arr(ji,jj,jk) .NE. 0) then
CCC              zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
           zbtr = zbtr_arr(ji,jj,jk)
           zgm = zbtr * trn(ji,jj,jk,jn) * 
     &          ( un(ji  ,jj,jk) * e2u(ji  ,jj) * fse3u(ji  ,jj,jk)
     &          -un(ji-1,jj,jk) * e2u(ji-1,jj) * fse3u(ji-1,jj,jk))
           zgz = zbtr * trn(ji,jj,jk,jn) * 
     &          ( vn(ji,jj  ,jk) * e1v(ji,jj  ) * fse3v(ji,jj  ,jk)
     &          -vn(ji,jj-1,jk) * e1v(ji,jj-1) * fse3v(ji,jj-1,jk))
           trtrd(ji,jj,jk,jn,1) = trtrd(ji,jj,jk,jn,1) + zgm
           trtrd(ji,jj,jk,jn,2) = trtrd(ji,jj,jk,jn,2) + zgz
           trtrd(ji,jj,jk,jn,3) = trtrd(ji,jj,jk,jn,3)
     $          - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
   
        END DO       

CCC               endif
CCC            END DO
CCC          END DO
CCC        END DO

C Lateral boundary conditions on trtrd:
#   ifdef key_mpp
        CALL mpplnk( trtrd(1,1,1,jn,1), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,jn,2), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,jn,3), 1, 1 )
#   else      
        CALL lbc( trtrd(1,1,1,jn,1), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,jn,2), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,jn,3), 1, 1, 1, 1, jpk, 1 )
#   endif
#endif
C
C
C END of tracer loop
C ==================
C
 1000 CONTINUE
C
C
CCC 10 11 2004  F79 cronometer-stop

       trcadvparttime = MPI_WTIME() - trcadvparttime
       trcadvtottime = trcadvtottime + trcadvparttime
       write(*,*) "F79T:trcadvparttime", trcadvparttime
       write(*,*) "F79T:trcadvtottime", trcadvtottime

CCC       trcadvparttime = MPI_WTIME()
CCC       call mppsync()
CCC       trcadvparttime = MPI_WTIME() - trcadvparttime
CCC       write(*,*) "F79T:trcadvBARttime", trcadvparttime

       write(*,*) "F79T:ADVallpoints = ", allpoints, " ADVgoodpoints = ",
     $            goodpoints, " ADVtpoints = ", tpoints

CCC

C$$$      WRITE (0,*) ' tendance smolar ', tra(42,62,1,3)

CCC
