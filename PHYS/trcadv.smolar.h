CCC $Id: trcadv.smolar.h,v 1.2 2009-09-11 09:20:56 cvsogs01 Exp $
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
      LOGICAL :: MPI_CHECK
      INTEGER ji,jj,jk,jt,jn
      INTEGER pack_size
      REAL(8) zbtr,zdt
      REAL(8) zgm,zgz

CC----------------------------------------------------------------------
CC statement functions
CC ===================

#include "stafun.h"


! omp variables
      INTEGER :: mytid, ntids, itid

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (11/96)
CCC---------------------------------------------------------------------

      MPI_CHECK = .FALSE.

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000

#else
      ntids = 1
      mytid = 0
#endif

      if(allpoints .EQ. 0) then

!$omp parallel default(shared) private(mytid)
#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
         zti(:,:,:,mytid+1)  = 0.
         ztj(:,:,:,mytid+1)  = 0.
         zkx(:,:,:,mytid+1)  = 0.
         zky(:,:,:,mytid+1)  = 0.
         zkz(:,:,:,mytid+1)  = 0.
         zbuf = 0.
         zx(:,:,:,mytid+1)   = 0.
         zy(:,:,:,mytid+1)   = 0. 
         zz(:,:,:,mytid+1)   = 0.
!$omp end parallel

         zaa  = 0.
         zbb  = 0.
         zcc  = 0.
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

         write(*,*) 'trcadv: RANK -> ', RANK, ' all_points -> ', allpoints

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

         write(*,*) 'trcadv: RANK -> ', RANK, ' good_points -> ', goodpoints

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


      zdt = rdt*float(ndttrc)

CCC
!$omp  parallel default(none) private(mytid,jk,jj,ji) shared(jpk,jpj,jpi,ztj,zx,zy,zz)
#ifdef __OPENMP
       mytid = omp_get_thread_num()  ! take the thread ID
#endif
       DO jk = 1,jpk
          DO jj = 1,jpj
             DO ji = 1,jpi
                ztj(ji,jj,jk,mytid+1) = 0.
                zx(ji,jj,jk,mytid+1) = 0.
                zy(ji,jj,jk,mytid+1) = 0.
                zz(ji,jj,jk,mytid+1) = 0.
             END DO
          END DO
       END DO
!$omp  end parallel

       DO jk = 1,jpk
          DO jj = 1,jpj
             DO ji = 1,jpi
                zaa(ji,jj,jk) = e2u(ji,jj)*fse3u(ji,jj,jk) * un(ji,jj,jk)
                zbb(ji,jj,jk) = e1v(ji,jj)*fse3v(ji,jj,jk) * vn(ji,jj,jk)
                zcc(ji,jj,jk) = e1t(ji,jj)*e2t(ji,jj)      * wn(ji,jj,jk)
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
C     tracer loop parallelized (macrotasking)
C     =======================================
C
      trcadvparttime = MPI_WTIME()
      TRACER_LOOP: DO  jn = ktask, jptra, ntids

C
C        1. tracer flux in the 3 directions
C        ----------------------------------
C
C        1.1 mass flux at u v and t-points and initialization
C

#if defined key_trc_diatrd

!$omp    parallel default(none) private(jk,jj,ji,mytid) shared(jpk,jpj,jpi,jn,trtrd)

#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
	 if( mytid + jn <= jptra ) then

           DO jk = 1,jpk
             DO jj = 1,jpj
               DO ji = 1,jpi
                 trtrd(ji,jj,jk,jn+mytid,1) = 0.
                 trtrd(ji,jj,jk,jn+mytid,2) = 0.
                 trtrd(ji,jj,jk,jn+mytid,3) = 0.
               END DO
             END DO
           END DO

	 end if

!$omp    end parallel

#endif

C
C       1.2 calcul of intermediate field with an upstream advection scheme
C           and mass fluxes calculated above
C
C       calcul of tracer flux in the i and j direction
C
!$omp   parallel default(none) private(jk,jj,ji,mytid)
!$omp&      shared(jpk,jpj,jpi,zkx,zky,zbb,zaa,zcc,jarr1,dimen_jarr1,jn,zkz,trn,jpjm1,jpim1)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif

	if( mytid + jn <= jptra ) then

           DO jk = 1,jpk
             DO jj=1,jpj
               zkx(1,jj,jk,mytid+1)=0.
             END DO
             DO jj=1,jpj
               zkx(jpi,jj,jk,mytid+1)=0.
             END DO
   
             DO ji=1,jpi
               zky(ji,1,jk,mytid+1)=0.
             END DO
             DO ji=1,jpi
               zky(ji,jpj,jk,mytid+1)=0.
             END DO
           END DO

C
C       calcul of tracer flux in the k direction
C
           DO jj = 1,jpj
             DO ji = 1,jpi
               zkz(ji,jj,1,mytid+1) = 0.
             END DO
           END DO

           jk = 1
           DO jj = 2,jpjm1
              DO ji = 2,jpim1
                 zkx(ji,jj,jk,mytid+1) = fsx(
     $            trn(ji,jj,jk,jn+mytid),trn(ji + 1,jj,jk,jn+mytid),zaa(ji,jj,jk))
                 zky(ji,jj,jk,mytid+1) = fsy(
     $            trn(ji,jj,jk,jn+mytid),trn(ji,jj + 1,jk,jn+mytid),zbb(ji,jj,jk))
              END DO
           END DO

           DO jk = 2,jpk         
             jj=1
             DO ji = 1,jpi
                zkz(ji,jj,jk,mytid+1) = fsz(
     $            trn(ji,jj,jk,jn+mytid),trn(ji,jj,jk - 1,jn+mytid),zcc(ji,jj,jk))
             END DO
             jj= jpj
             DO ji = 1,jpi
                zkz(ji,jj,jk,mytid+1) = fsz(
     $            trn(ji,jj,jk,jn+mytid),trn(ji,jj,jk - 1,jn+mytid),zcc(ji,jj,jk))
             END DO
             DO jj = 2,jpjm1
                ji=1
                zkz(ji,jj,jk,mytid+1) = fsz(
     $            trn(ji,jj,jk,jn+mytid),trn(ji,jj,jk - 1,jn+mytid),zcc(ji,jj,jk))
                ji=jpi
                zkz(ji,jj,jk,mytid+1) = fsz(
     $            trn(ji,jj,jk,jn+mytid),trn(ji,jj,jk - 1,jn+mytid),zcc(ji,jj,jk))
             END DO
           END DO

           DO ju=1, dimen_jarr1

              ji = jarr1(1, ju)
              jj = jarr1(2, ju)
              jk = jarr1(3, ju)

              zkx(ji,jj,jk,mytid+1) = fsx(
     $            trn(ji,jj,jk,jn+mytid),trn(ji + 1,jj,jk,jn+mytid),zaa(ji,jj,jk))
              zky(ji,jj,jk,mytid+1) = fsy(
     $            trn(ji,jj,jk,jn+mytid),trn(ji,jj + 1,jk,jn+mytid),zbb(ji,jj,jk))
              zkz(ji,jj,jk,mytid+1) = fsz(
     $            trn(ji,jj,jk,jn+mytid),trn(ji,jj,jk - 1,jn+mytid),zcc(ji,jj,jk))    

           END DO

	end if

!$omp end parallel

C
C ... Lateral boundary conditions on zk[xy]
#ifdef key_mpp
C
C   ... Mpp : export boundary values to neighboring processors
C
        IF(MPI_CHECK) CALL mt_trace_start()
CCCCCCCCC        trcadvparttime = MPI_WTIME()

C       DO itid = 1, ntids

C  IF( itid - 1 + jn <= jptra ) THEN


C           CALL mpplnk( zkx(:,:,:,itid), 1, 1 )
C           CALL mpplnk( zky(:,:,:,itid), 1, 1 )

C  END IF

C       END DO

	IF( ntids - 1 + jn <= jptra ) THEN
	   pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
	END IF

        CALL mpplnk_my(zkx(:,:,:,:), pack_size,1,1)
        CALL mpplnk_my(zky(:,:,:,:), pack_size,1,1)

        IF(MPI_CHECK) CALL mt_trace_stop()
CCCCCCCCC        trcadvparttime = MPI_WTIME() - trcadvparttime
CCCCCCCCC        trcadvtottime = trcadvtottime + trcadvparttime


#else
C
C   ... T-point, 3D array, full local arrays zk[xy] are initialised
C       
C
        DO itid = 1, ntids

	   IF( itid - 1 + jn <= jptra ) THEN

              CALL lbc( zkx(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
              CALL lbc( zky(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )

	   END IF

        END DO
C
#endif
C
C
C 2. calcul of after field using an upstream advection scheme
C -----------------------------------------------------------
C

!$omp   parallel default(none) private(mytid,zbtr,ji,jj,jk,ju)
!$omp&      shared(zkx,zky,zkz,zti,jpim1,jpjm1,trn,zdt,jn,jpkm1,zbtr_arr,e1t,e2t,ztj,jarr3,ncor,dimen_jarr3)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif

	IF( mytid + jn <= jptra ) THEN

           DO ju=1, dimen_jarr3

              ji = jarr3(1, ju)
              jj = jarr3(2, ju)
              jk = jarr3(3, ju)

              zbtr = zbtr_arr(ji,jj,jk)
              ztj(ji,jj,jk,mytid+1) = -zbtr*
     $          ( zkx(ji,jj,jk,mytid+1) - zkx(ji - 1,jj,jk,mytid+1)
     $          + zky(ji,jj,jk,mytid+1) - zky(ji,jj - 1,jk,mytid+1)
     $          + zkz(ji,jj,jk,mytid+1) - zkz(ji,jj,jk + 1,mytid+1) )
#if defined key_trc_diatrd
              trtrd(ji,jj,jk,jn+mytid,1) = trtrd(ji,jj,jk,jn+mytid,1) -
     $          zbtr*( zkx(ji,jj,jk,mytid+1) - zkx(ji - 1,jj,jk,mytid+1) )
              trtrd(ji,jj,jk,jn+mytid,2) = trtrd(ji,jj,jk,jn+mytid,2) -
     $          zbtr*( zky(ji,jj,jk,mytid+1) - zky(ji,jj - 1,jk,mytid+1) )
              trtrd(ji,jj,jk,jn+mytid,3) = trtrd(ji,jj,jk,jn+mytid,3) -
     $          zbtr*( zkz(ji,jj,jk,mytid+1) - zkz(ji,jj,jk + 1,mytid+1) )
#endif

           END DO

	END IF

!$omp end parallel 

C
C 2.1 start of antidiffusive correction loop
C
        ANTIDIFF_CORR: DO jt = 1,ncor
C
C 2.2 calcul of intermediary field zti
C

!$omp     parallel default(none) private(mytid,ji,jj,jk)
!$omp&       shared(jt,jn,ncor,jpkm1,jpjm1,jpim1,zti,ztj,trn,zdt)

#ifdef __OPENMP
          mytid = omp_get_thread_num()  ! take the thread ID
#endif

	  IF( mytid + jn <= jptra ) THEN

             if(jt .EQ. 1) then

                if(ncor .EQ. 1) then
                   DO jk = 1,jpkm1
                      DO jj = 2,jpjm1
                         DO ji = 2,jpim1
                            zti(ji,jj,jk,mytid+1) = trn(ji,jj,jk,jn+mytid) + zdt*ztj(ji,jj,jk,mytid+1)
CCC                         zbuf(ji,jj,jk) = 0. + ztj(ji,jj,jk,mytid+1)
                         END DO
                      END DO
                   END DO
                else
                    DO jk = 1,jpkm1
                      DO jj = 2,jpjm1
                         DO ji = 2,jpim1
                            zti(ji,jj,jk,mytid+1) = trn(ji,jj,jk,jn+mytid) + zdt*ztj(ji,jj,jk,mytid+1)
                            zbuf(ji,jj,jk) = 0. + ztj(ji,jj,jk,mytid+1)
                         END DO
                      END DO
                   END DO
                endif

             else

                DO jk = 1,jpkm1
                   DO jj = 2,jpjm1
                      DO ji = 2,jpim1
                         zti(ji,jj,jk,mytid+1) = zti(ji,jj,jk,mytid+1) + zdt*ztj(ji,jj,jk,mytid+1)
                         zbuf(ji,jj,jk) = zbuf(ji,jj,jk) + ztj(ji,jj,jk,mytid+1)
                      END DO
                   END DO
                END DO

             endif

	  END IF

!$omp     end parallel


C
C
C ... Lateral boundary conditions on zti
C
#ifdef key_mpp
C
C   ... Mpp : export boundary values to neighboring processors
C
        IF(MPI_CHECK) CALL mt_trace_start()
CCCCCCCCC        trcadvparttime = MPI_WTIME()

C        DO itid = 1, ntids
C    IF( itid - 1 + jn <= jptra ) THEN
C              CALL mpplnk( zti(:,:,:,itid), 1, 1 )
C     END IF
C         END DO

        IF( ntids - 1 + jn <= jptra ) THEN
           pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
        END IF

        CALL mpplnk_my(zti(:,:,:,:), pack_size,1,1)


        IF(MPI_CHECK) CALL mt_trace_stop()
CCCCCCCCC        trcadvparttime = MPI_WTIME() - trcadvparttime
CCCCCCCCC        trcadvtottime = trcadvtottime + trcadvparttime


#else
C
C   ... T-point, 3D array, full local array zti is initialised
C       
          DO itid = 1, ntids
	     IF( itid - 1 + jn <= jptra ) THEN
                CALL lbc( zti(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
	     END IF
          END DO
C
#endif
C
C
C 2.3 calcul of the antidiffusive flux
C
!$omp     parallel default(none) private(mytid,ji,jj,jk)
!$omp&       shared(jn,jpkm1,jpjm1,jpim1,zti,ztj,zy,zx,zz,jarr2,big_fact_zbb,
!$omp&              big_fact_zaa,big_fact_zcc,dimen_jarr2,rtrn,rsc)

#ifdef __OPENMP
          mytid = omp_get_thread_num()  ! take the thread ID
#endif

	  IF( mytid + jn <= jptra ) THEN

            jk = 1
CCC          DO jk = 1,jpkm1
            DO jj = 2,jpjm1
              DO ji = 2,jpim1

CCC                 if(advmask(ji,jj,jk) .NE. 0) then

CCC                zx(ji,jj,jk,mytid+1) = ( abs(zaa(ji,jj,jk)) - zdt
CCC     $              *zaa(ji,jj,jk)**2
CCCCCC     $              /(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) ) )
CCC     $              *inv_eu(ji,jj,jk) )

                    zx(ji,jj,jk,mytid+1) = big_fact_zaa(ji,jj,jk)
     $                   *(zti(ji + 1,jj,jk,mytid+1) - zti( ji ,jj,jk,mytid+1))
     $                   /(zti( ji ,jj,jk,mytid+1) + zti(ji + 1,jj,jk,mytid+1) + rtrn)
     $                   * rsc

CCC                zy(ji,jj,jk,mytid+1) = ( abs(zbb(ji,jj,jk)) - zdt
CCC     $              *zbb(ji,jj,jk)**2
CCCCCC     $              /(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) ) )
CCC     $              *inv_ev(ji,jj,jk) )

                    zy(ji,jj,jk,mytid+1) = big_fact_zbb(ji,jj,jk)
     $                   *(zti(ji,jj + 1,jk,mytid+1) - zti(ji, jj ,jk,mytid+1))
     $                   /(zti(ji, jj ,jk,mytid+1) + zti(ji,jj + 1,jk,mytid+1) + rtrn)
     $                   * rsc


CCC                endif

              END DO
            END DO
CCC          END DO

CCC                 if(advmask(ji,jj,jk) .NE. 0) then

            DO ju=1, dimen_jarr2

               ji = jarr2(1, ju)
               jj = jarr2(2, ju)
               jk = jarr2(3, ju)

CCC                zx(ji,jj,jk,mytid+1) = ( abs(zaa(ji,jj,jk)) - zdt
CCC     $              *zaa(ji,jj,jk)**2
CCCCCC     $              /(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) ) )
CCC     $              *inv_eu(ji,jj,jk) )
                    zx(ji,jj,jk,mytid+1) = big_fact_zaa(ji,jj,jk)
     $                   *(zti(ji + 1,jj,jk,mytid+1) - zti( ji ,jj,jk,mytid+1))
     $                   /(zti( ji ,jj,jk,mytid+1) + zti(ji + 1,jj,jk,mytid+1) + rtrn)
     $                   * rsc
CCC                zy(ji,jj,jk,mytid+1) = ( abs(zbb(ji,jj,jk)) - zdt
CCC     $              *zbb(ji,jj,jk)**2
CCCCCC     $              /(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) ) )
CCC     $              *inv_ev(ji,jj,jk) )
                    zy(ji,jj,jk,mytid+1) = big_fact_zbb(ji,jj,jk)
     $                   *(zti(ji,jj + 1,jk,mytid+1) - zti(ji, jj ,jk,mytid+1))
     $                   /(zti(ji, jj ,jk,mytid+1) + zti(ji,jj + 1,jk,mytid+1) + rtrn)
     $                   * rsc

CCC                zz(ji,jj,jk) = ( abs(zcc(ji,jj,jk)) - zdt
CCC     $              *zcc(ji,jj,jk)**2
CCCCCC     $              /(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk) ) )
CCC     $              *inv_et(ji,jj,jk) )
                    zz(ji,jj,jk,mytid+1) = big_fact_zcc(ji,jj,jk)
     $                   *(zti(ji,jj,jk,mytid+1) - zti(ji,jj,jk - 1,mytid+1))
     $                   /(zti(ji,jj,jk,mytid+1) + zti(ji,jj,jk - 1,mytid+1) + rtrn)
     $                   * rsc*( -1.)

           END DO
CCC                 endif

          END IF

!$omp end parallel 


C---------------------------------------------------------------------------
C 2.4 cross terms
C
CCC          IF (crosster) THEN
CCC              DO jk = 2,jpkm1
CCC                DO jj = 2,jpjm1
CCC                 DO ji = 2,jpim1
CCC                    zx(ji,jj,jk) = zx(ji,jj,jk)
CCC     $                  - 0.5*zdt*rsc*zaa(ji,jj,jk)*0.25*
CCC     $                  (    (zbb(ji  ,jj - 1,jk  ) + zbb(ji + 1,jj - 1
CCC     $                  ,jk  ) + zbb(ji + 1,jj  ,jk  ) + zbb(ji  ,jj
CCC     $                  ,jk))* (zti(ji  ,jj + 1,jk  ) + zti(ji + 1,jj +
CCC     $                  1,jk  ) - zti(ji + 1,jj - 1,jk  ) - zti(ji  ,jj
CCC     $                  - 1,jk  ))/ (zti(ji  ,jj + 1,jk  ) + zti(ji + 1
CCC     $                  ,jj + 1,jk  ) + zti(ji + 1,jj - 1,jk  ) + zti(ji
CCC     $                  ,jj - 1,jk  ) + rtrn) + (zcc(ji  ,jj  ,jk  ) +
CCC     $                  zcc(ji + 1,jj  ,jk  ) + zcc(ji  ,jj  ,jk + 1) +
CCC     $                  zcc(ji + 1,jj  ,jk + 1))* (zti(ji  ,jj  ,jk - 1)
CCC     $                  + zti(ji + 1,jj  ,jk - 1) - zti(ji  ,jj  ,jk + 1
CCC     $                  )- zti(ji + 1,jj  ,jk + 1))/ (zti(ji  ,jj  ,jk -
CCC     $                  1) + zti(ji + 1,jj  ,jk - 1) + zti(ji  ,jj  ,jk
CCC     $                  +1) + zti(ji + 1,jj  ,jk + 1) + rtrn))/(e1u(ji
CCC     $                  ,jj)*e2u(ji,jj)*fse3u(ji,jj,jk))*vmask(ji  ,jj -
CCC     $                  1,jk  )*vmask(ji + 1,jj - 1,jk  )*vmask(ji + 1
CCC     $                  ,jj,jk)*vmask(ji  ,jj  ,jk  )*tmask(ji  ,jj  ,jk
CCC     $                  )*tmask(ji + 1,jj  ,jk  )*tmask(ji  ,jj  ,jk + 1
CCC     $                  )*tmask(ji + 1,jj  ,jk + 1)
CCC                    zy(ji,jj,jk) = zy(ji,jj,jk)
CCC     $                  - 0.5*zdt*rsc*zbb(ji,jj,jk)*0.25*
CCC     $                  (    (zaa(ji - 1,jj  ,jk  ) + zaa(ji - 1,jj + 1
CCC     $                  ,jk  ) + zaa(ji  ,jj  ,jk  ) + zaa(ji  ,jj + 1
CCC     $                  ,jk))* (zti(ji + 1,jj + 1,jk  ) + zti(ji + 1,jj
CCC     $                  ,jk  ) - zti(ji - 1,jj + 1,jk  ) - zti(ji - 1,jj
CCC     $                  ,jk  ))/ (zti(ji + 1,jj + 1,jk  ) + zti(ji + 1
CCC     $                  ,jj  ,jk  ) + zti(ji - 1,jj + 1,jk  ) + zti(ji
CCC     $                  - 1,jj  ,jk  ) + rtrn) + (zcc(ji  ,jj  ,jk  )
CCC     $                  + zcc(ji  ,jj  ,jk + 1) + zcc(ji  ,jj + 1,jk  )
CCC     $                  + zcc(ji  ,jj + 1,jk + 1))* (zti(ji  ,jj  ,jk -
CCC     $                  1) + zti(ji  ,jj + 1,jk - 1) - zti(ji  ,jj  ,jk
CCC     $                  +1) - zti(ji  ,jj + 1,jk + 1))/ (zti(ji  ,jj
CCC     $                  ,jk- 1) + zti(ji  ,jj + 1,jk - 1) + zti(ji  ,jj
CCC     $                  ,jk+ 1) + zti(ji  ,jj + 1,jk + 1) + rtrn))
CCC     $                  /(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk))
CCC     $                  *umask(ji - 1,jj,jk  )*umask(ji - 1,jj + 1,jk  )
CCC     $                  *umask(ji  ,jj,jk  )*umask(ji  ,jj + 1,jk  )
CCC     $                  *tmask(ji  ,jj,jk)*tmask(ji  ,jj  ,jk + 1)
CCC     $                  *tmask(ji  ,jj + 1,jk)*tmask(ji  ,jj + 1,jk + 1)
CCC                    zz(ji,jj,jk) = zz(ji,jj,jk)
CCC    $                  - 0.5*zdt*rsc*zcc(ji,jj,jk)*0.25*
CCC     $                  (    (zaa(ji - 1,jj  ,jk  ) + zaa(ji  ,jj  ,jk
CCC     $                  ) + zaa(ji  ,jj  ,jk - 1) + zaa(ji - 1,jj  ,jk -
CCC     $                  1))*(zti(ji + 1,jj  ,jk - 1) + zti(ji + 1,jj
CCC     $                  ,jk  ) - zti(ji - 1,jj  ,jk  ) - zti(ji - 1,jj
CCC     $                  ,jk - 1))/(zti(ji + 1,jj  ,jk - 1) + zti(ji + 1
CCC     $                  ,jj,jk  ) + zti(ji - 1,jj  ,jk  ) + zti(ji - 1
CCC     $                  ,jj,jk - 1) + rtrn) + (zbb(ji  ,jj - 1,jk  )
CCC     $                  + zbb(ji  ,jj  ,jk  ) + zbb(ji  ,jj  ,jk - 1)
CCC     $                  + zbb(ji  ,jj - 1,jk - 1))*(zti(ji  ,jj + 1,jk -
CCC     $                  1) + zti(ji  ,jj + 1,jk  ) - zti(ji  ,jj - 1,jk
CCC     $                  ) - zti(ji  ,jj - 1,jk - 1))/(zti(ji  ,jj + 1,jk
CCC     $                  - 1) + zti(ji  ,jj + 1,jk  ) + zti(ji  ,jj - 1
CCC     $                  ,jk  ) + zti(ji  ,jj - 1,jk - 1) + rtrn))
CCC     $                  /(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk))
CCC     $                  *umask(ji - 1,jj,jk  )*umask(ji  ,jj  ,jk  )
CCC     $                  *umask(ji  ,jj,jk- 1)*umask(ji - 1,jj  ,jk - 1)
CCC     $                  *vmask(ji  ,jj- 1,jk)*vmask(ji  ,jj  ,jk  )
CCC     $                  *vmask(ji  ,jj  ,jk-1)*vmask(ji  ,jj - 1,jk - 1)
CCC                  END DO
CCC               END DO
CCC              END DO
CCCC
CCC              DO jj = 2,jpjm1
CCC                DO ji = 2,jpim1
CCC                  zx(ji,jj,1) = zx(ji,jj,1)
CCC     $                - 0.5*zdt*rsc*zaa(ji,jj,1)*0.25*
CCC     $                ( (zbb(ji  ,jj - 1,1  ) + zbb(ji + 1,jj - 1,1  )
CCC     $                + zbb(ji + 1,jj  ,1  ) + zbb(ji  ,jj  ,1  ))
CCC     $                *(zti(ji  ,jj + 1,1  ) + zti(ji + 1,jj + 1,1  )
CCC     $                - zti(ji + 1,jj - 1,1  ) - zti(ji  ,jj - 1,1  ))
CCC     $                /(zti(ji  ,jj + 1,1  ) + zti(ji + 1,jj + 1,1  )
CCC     $                + zti(ji + 1,jj - 1,1  ) + zti(ji  ,jj - 1,1  ) +
CCC     $                rtrn))/(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,1))
CCC     $                *vmask(ji  ,jj - 1,1  )*vmask(ji + 1,jj - 1,1  )
CCC     $                *vmask(ji + 1,jj  ,1  )*vmask(ji  ,jj  ,1  )
CCC                  zy(ji,jj,1) = zy(ji,jj,1)
CCC     $                - 0.5*zdt*rsc*zbb(ji,jj,1)*0.25*
CCC     $                ( (zaa(ji - 1,jj  ,1  ) + zaa(ji - 1,jj + 1,1  )
CCC     $                + zaa(ji  ,jj  ,1  ) + zaa(ji  ,jj + 1,1  ))
CCC     $                *(zti(ji + 1,jj + 1,1  ) + zti(ji + 1,jj  ,1  )
CCC     $                - zti(ji - 1,jj + 1,1  ) - zti(ji - 1,jj  ,1  ))
CCC     $                /(zti(ji + 1,jj + 1,1  ) + zti(ji + 1,jj  ,1  )
CCC     $                + zti(ji - 1,jj + 1,1  ) + zti(ji - 1,jj  ,1  ) +
CCC     $                rtrn))/(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,1))
CCC     $                *umask(ji - 1,jj  ,1  )*umask(ji - 1,jj + 1,1  )
CCC     $                *umask(ji  ,jj  ,1  )*umask(ji  ,jj + 1,1  )
CCC                END DO
CCC              END DO
CCC          ENDIF
C-----------------------------------------------------------------------------


C
C ... Lateral boundary conditions on z[xyz]
#ifdef key_mpp
C
C   ... Mpp : export boundary values to neighboring processors
C
        IF(MPI_CHECK) CALL mt_trace_start()
CCCCCCCCC        trcadvparttime = MPI_WTIME()

C         DO itid = 1, ntids
C     IF( itid - 1 + jn <= jptra ) THEN
C               CALL mpplnk( zx(:,:,:,itid), 1, 1 )
C               CALL mpplnk( zy(:,:,:,itid), 1, 1 )
C               CALL mpplnk( zz(:,:,:,itid), 1, 1 )
C     END IF
C  END DO

	IF( ntids - 1 + jn <= jptra ) THEN
	   pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
	END IF

        CALL mpplnk_my(zx(:,:,:,:), pack_size,1,1)
        CALL mpplnk_my(zy(:,:,:,:), pack_size,1,1)
        CALL mpplnk_my(zz(:,:,:,:), pack_size,1,1)

        IF(MPI_CHECK) CALL mt_trace_stop()
CCCCCCCCC        trcadvparttime = MPI_WTIME() - trcadvparttime
CCCCCCCCC        trcadvtottime = trcadvtottime + trcadvparttime

#else
C
C   ... T-point, 3D array, full local array z[xyz] are initialised
C       
C
          DO itid = 1, ntids
	     IF( itid - 1 + jn <= jptra ) THEN
                CALL lbc( zx(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zy(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
                CALL lbc( zz(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
	     END IF
	  END DO
C
#endif
C
C
C 2.4 reinitialization
C

!$omp     parallel default(none) private(mytid,jk,jj,ji)
!$omp&      shared(zkx,zky,zkz,zz,zx,zy,zti,jpjm1,jpim1,dimen_jarr1,jarr1,jpi,jpk,jpj,jn)

#ifdef __OPENMP
          mytid = omp_get_thread_num()  ! take the thread ID
#endif

	  IF( mytid + jn <= jptra ) THEN

C
C            2.5 calcul of the final field:
C                advection by antidiffusive mass fluxes and an upstream scheme
C
             jk = 1
             DO jj = 2,jpjm1
                 DO ji = 2,jpim1
                   zkx(ji,jj,jk,mytid+1) = fsx(
     $              zti(ji,jj,jk,mytid+1),zti(ji + 1,jj,jk,mytid+1),zx(ji,jj,jk,mytid+1))
CCC     $              zti(ji,jj,jk,mytid+1),zti(ji + 1,jj,jk,mytid+1),zaa(ji,jj,jk))
                   zky(ji,jj,jk,mytid+1) = fsy(
     $              zti(ji,jj,jk,mytid+1),zti(ji,jj + 1,jk,mytid+1),zy(ji,jj,jk,mytid+1))
CCC     $              zti(ji,jj,jk,mytid+1),zti(ji,jj + 1,jk,mytid+1),zbb(ji,jj,jk))
                 END DO
              END DO

             DO jk = 2,jpk
               jj = 1
               DO ji = 1,jpi
                  zkz(ji,jj,jk,mytid+1) = fsz(
     $              zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zz(ji,jj,jk,mytid+1))
CCC   $                 zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zcc(ji,jj,jk))
               END DO
               jj = jpj
               DO ji = 1,jpi
                  zkz(ji,jj,jk,mytid+1) = fsz(
     $              zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zz(ji,jj,jk,mytid+1))
CCC   $                 zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zcc(ji,jj,jk))
               END DO

               DO jj = 2,jpjm1
                 ji = 1
                 zkz(ji,jj,jk,mytid+1) = fsz(
     $             zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zz(ji,jj,jk,mytid+1))
CCC   $              zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zcc(ji,jj,jk))
                 ji = jpi
                 zkz(ji,jj,jk,mytid+1) = fsz(
     $             zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zz(ji,jj,jk,mytid+1))
CCC   $              zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zcc(ji,jj,jk))

               END DO
             END DO

             DO ju=1, dimen_jarr1

                ji = jarr1(1, ju)
                jj = jarr1(2, ju)
                jk = jarr1(3, ju)  

                zkx(ji,jj,jk,mytid+1) = fsx(
     $            zti(ji,jj,jk,mytid+1),zti(ji + 1,jj,jk,mytid+1),zx(ji,jj,jk,mytid+1))
CCC     $              zti(ji,jj,jk,mytid+1),zti(ji + 1,jj,jk,mytid+1),zaa(ji,jj,jk))
                zky(ji,jj,jk,mytid+1) = fsy(
     $            zti(ji,jj,jk,mytid+1),zti(ji,jj + 1,jk,mytid+1),zy(ji,jj,jk,mytid+1))
CCC     $              zti(ji,jj,jk,mytid+1),zti(ji,jj + 1,jk,mytid+1),zbb(ji,jj,jk))
                zkz(ji,jj,jk,mytid+1) = fsz(
     $            zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zz(ji,jj,jk,mytid+1))
CCC     $              zti(ji,jj,jk,mytid+1),zti(ji,jj,jk - 1,mytid+1),zcc(ji,jj,jk))
             
             END DO

	  END IF

!$omp end parallel

C
C
C ... Lateral boundary conditions on zk[xy]
#ifdef key_mpp
C
C   ... Mpp : export boundary values to neighboring processors
C
        IF(MPI_CHECK) CALL mt_trace_start()
CCCCCCCCC        trcadvparttime = MPI_WTIME()

C        DO itid = 1, ntids
C    IF( itid - 1 + jn <= jptra ) THEN
C              CALL mpplnk( zkx(:,:,:,itid), 1, 1 )
C              CALL mpplnk( zky(:,:,:,itid), 1, 1 )
C    END IF
C END DO

        IF( ntids - 1 + jn <= jptra ) THEN
           pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
        END IF

        CALL mpplnk_my(zkx(:,:,:,:), pack_size,1,1)
        CALL mpplnk_my(zky(:,:,:,:), pack_size,1,1)


        IF(MPI_CHECK) CALL mt_trace_stop()
CCCCCCCCC        trcadvparttime = MPI_WTIME() - trcadvparttime
CCCCCCCCC        trcadvtottime = trcadvtottime + trcadvparttime

#else
C
C   ... T-point, 3D array, full local array zk[xy] are initialised
C       
C
         DO itid = 1, ntids
	    IF( itid - 1 + jn <= jptra ) THEN
               CALL lbc( zkx(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
               CALL lbc( zky(:,:,:,itid), 1, 1, 1, 1, jpk, 1 )
	    END IF
	 END DO
C
#endif

!$omp    parallel default(none) private(mytid,zbtr,ji,jj,jk,ju)
!$omp&      shared(zkx,zky,zkz,zbtr_arr,e1t,e2t,ztj,dimen_jarr3,jarr3,ncor,jn)

#ifdef __OPENMP
         mytid = omp_get_thread_num()  ! take the thread ID
#endif
C
C
C        2.6. calcul of after field using an upstream advection scheme
C

         IF( mytid + jn <= jptra ) THEN

            if(ncor .EQ. 1) then
               DO ju=1, dimen_jarr3
               
                  ji = jarr3(1, ju)
                  jj = jarr3(2, ju)
                  jk = jarr3(3, ju)

CCC                if(zbtr_arr(ji,jj,jk) .NE. 0) then
CCC                zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
                  zbtr = zbtr_arr(ji,jj,jk)
                  ztj(ji,jj,jk,mytid+1) = -zbtr*
     $              ( zkx(ji,jj,jk,mytid+1) - zkx(ji - 1,jj,jk,mytid+1)
     $              + zky(ji,jj,jk,mytid+1) - zky(ji,jj - 1,jk,mytid+1)
     $              + zkz(ji,jj,jk,mytid+1) - zkz(ji,jj,jk + 1,mytid+1) )
     $              + ztj(ji,jj,jk,mytid+1)
#if defined key_trc_diatrd
                  trtrd(ji,jj,jk,jn+mytid,1) = trtrd(ji,jj,jk,jn+mytid,1) -
     $              zbtr*( zkx(ji,jj,jk,mytid+1) - zkx(ji - 1,jj,jk,mytid+1) )
                  trtrd(ji,jj,jk,jn+mytid,2) = trtrd(ji,jj,jk,jn+mytid,2) -
     $              zbtr*( zky(ji,jj,jk,mytid+1) - zky(ji,jj - 1,jk,mytid+1) )
                  trtrd(ji,jj,jk,jn+mytid,3) = trtrd(ji,jj,jk,jn+mytid,3) -
     $              zbtr*( zkz(ji,jj,jk,mytid+1) - zkz(ji,jj,jk + 1,mytid+1) )
#endif

              END DO
           else
              DO ju=1, dimen_jarr3
               
                  ji = jarr3(1, ju)
                  jj = jarr3(2, ju)
                  jk = jarr3(3, ju)

CCC                if(zbtr_arr(ji,jj,jk) .NE. 0) then
CCC                zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
                  zbtr = zbtr_arr(ji,jj,jk)
                  ztj(ji,jj,jk,mytid+1) = -zbtr*
     $              ( zkx(ji,jj,jk,mytid+1) - zkx(ji - 1,jj,jk,mytid+1)
     $              + zky(ji,jj,jk,mytid+1) - zky(ji,jj - 1,jk,mytid+1)
     $              + zkz(ji,jj,jk,mytid+1) - zkz(ji,jj,jk + 1,mytid+1) )
#if defined key_trc_diatrd
                  trtrd(ji,jj,jk,jn+mytid,1) = trtrd(ji,jj,jk,jn+mytid,1) -
     $              zbtr*( zkx(ji,jj,jk,mytid+1) - zkx(ji - 1,jj,jk,mytid+1) )
                  trtrd(ji,jj,jk,jn+mytid,2) = trtrd(ji,jj,jk,jn+mytid,2) -
     $              zbtr*( zky(ji,jj,jk,mytid+1) - zky(ji,jj - 1,jk,mytid+1) )
                  trtrd(ji,jj,jk,jn+mytid,3) = trtrd(ji,jj,jk,jn+mytid,3) -
     $              zbtr*( zkz(ji,jj,jk,mytid+1) - zkz(ji,jj,jk + 1,mytid+1) )
#endif

              END DO
           endif

	END IF

!$omp end parallel

        END DO ANTIDIFF_CORR

C
C       3. trend due to horizontal and vertical advection of tracer jn
C       --------------------------------------------------------------
C
!$omp   parallel default(none) private(mytid,ji,jj,jk,ju) shared(ncor,dimen_jarrt,jarrt,tra,ztj,jn)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif

        IF( mytid + jn <= jptra ) THEN

           if(ncor .EQ. 1) then
              DO ju=1, dimen_jarrt
                 ji = jarrt(1, ju)
                 jj = jarrt(2, ju)
                 jk = jarrt(3, ju)

                 tra(ji,jj,jk,jn+mytid) = tra(ji,jj,jk,jn+mytid)
     $             + ztj(ji,jj,jk,mytid+1)
CCC     $             + (zbuf(ji,jj,jk) + ztj(ji,jj,jk,mytid+1))
               
              END DO
           else
              DO ju=1, dimen_jarrt
                 ji = jarrt(1, ju)
                 jj = jarrt(2, ju)
                 jk = jarrt(3, ju)

                 tra(ji,jj,jk,jn+mytid) = tra(ji,jj,jk,jn+mytid)
     $             + (zbuf(ji,jj,jk) + ztj(ji,jj,jk,mytid+1))
   
              END DO
           endif

	END IF

!$omp end parallel
 

C
C
C 4.0 convert the transport trend into advection trend
C ----------------------------------------------------
C
#if defined key_trc_diatrd

        DO ju=1, dimen_jarr

           ji = jarr(1, ju)
           jj = jarr(2, ju)
           jk = jarr(3, ju)

CCC               if(zbtr_arr(ji,jj,jk) .NE. 0) then
CCC              zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
           zbtr = zbtr_arr(ji,jj,jk)
           zgm = zbtr * trn(ji,jj,jk,jn+mytid) * 
     &          ( un(ji  ,jj,jk) * e2u(ji  ,jj) * fse3u(ji  ,jj,jk)
     &          -un(ji-1,jj,jk) * e2u(ji-1,jj) * fse3u(ji-1,jj,jk))
           zgz = zbtr * trn(ji,jj,jk,jn+mytid) * 
     &          ( vn(ji,jj  ,jk) * e1v(ji,jj  ) * fse3v(ji,jj  ,jk)
     &          -vn(ji,jj-1,jk) * e1v(ji,jj-1) * fse3v(ji,jj-1,jk))
           trtrd(ji,jj,jk,jn+mytid,1) = trtrd(ji,jj,jk,jn+mytid,1) + zgm
           trtrd(ji,jj,jk,jn+mytid,2) = trtrd(ji,jj,jk,jn+mytid,2) + zgz
           trtrd(ji,jj,jk,jn+mytid,3) = trtrd(ji,jj,jk,jn+mytid,3)
     $          - trn(ji,jj,jk,jn+mytid) * hdivn(ji,jj,jk)
   
        END DO       


C Lateral boundary conditions on trtrd:
#   ifdef key_mpp
        IF(MPI_CHECK) CALL mt_trace_start()
CCCCCCCCC        trcadvparttime = MPI_WTIME()

        CALL mpplnk( trtrd(1,1,1,jn+mytid,1), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,jn+mytid,2), 1, 1 )
        CALL mpplnk( trtrd(1,1,1,jn+mytid,3), 1, 1 )

        IF(MPI_CHECK) CALL mt_trace_stop()
CCCCCCCCC        trcadvparttime = MPI_WTIME() - trcadvparttime
CCCCCCCCC        trcadvtottime = trcadvtottime + trcadvparttime

#   else      
        CALL lbc( trtrd(1,1,1,jn+mytid,1), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,jn+mytid,2), 1, 1, 1, 1, jpk, 1 )
        CALL lbc( trtrd(1,1,1,jn+mytid,3), 1, 1, 1, 1, jpk, 1 )
#   endif

#endif
C
       END DO TRACER_LOOP
C
       trcadvparttime = MPI_WTIME() - trcadvparttime
       trcadvtottime = trcadvtottime + trcadvparttime
