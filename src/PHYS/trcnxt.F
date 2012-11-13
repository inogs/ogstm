      SUBROUTINE trcnxt(ktask,kt)
CCC---------------------------------------------------------------------
CCC
CCC                       ROUTINE trcnxt
CCC                     *******************
CCC
CCC  PURPOSE :
CCC  ---------
CCC    compute the passive tracers fields at the next time
CCC    step from their temporal trends
CCC
CC    METHOD :
CC
CC
CC      Apply lateral boundary conditions (nperio=0,closed ; nperio=1,
CC      east-west cyclic ; nperio=2, symmetric round the equator)
CC      on tra arrays
CC
CC   default:
CC      arrays swap
CC         (trn) = (tra) ; (tra) = (0,0)
CC         (trb) = (trn) 
CC
CC   For Arakawa Sheme : IF key_trc_arakawa defined
CC      A Asselin time filter applied on now tracers (trn) to avoid
CC      the divergence of two consecutive time-steps and tr arrays
CC      to prepare the next time_step:
CC         (trb) = (trn) + gamma [ (trb) + (tra) - 2 (trn) ]
CC         (trn) = (tra) ; (tra) = (0,0)
CC
CC      array swap for tracers to start the next time step
CC
CC
CC      macrotasking on tracer slab
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
CC            /comask/          : masks, bathymetry
CC
CC   OUTPUT :
CC   ------
CC      argument                : no
CC
CC   WORKSPACE :
CC   ---------
CC   ktask kt ji jj jk jn zfact zdt
CC

       USE myalloc
       USE myalloc_mpp
       USE FN_mem
       USE Time_Manager

      IMPLICIT NONE

CC local declarations
CC ==================
      INTEGER ktask,kt,pack_size
! omp variables
      INTEGER :: mytid, ntids


#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

      INTEGER ji,jj,jk,jn,jw


#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif


       trcnxtparttime = MPI_WTIME() ! cronometer-start

C 1. fields at the next time
C --------------------------
C Tracer slab
C ===========


      TRACER_LOOP: DO  jn = ktask, jptra, ntids


C 1. Lateral boundary conditions on tra (1,1,1,jn)

#ifdef key_mpp

C   ... Mpp : export boundary values to neighboring processors

        IF( ntids - 1 + jn <= jptra ) THEN
           pack_size = ntids
        ELSE
           pack_size = ntids - (ntids - 1 + jn - jptra)
        END IF

        CALL mpplnk_my(tra(1,1,1,jn), pack_size,1,1)

#  else

C   ... T-point, 3D array, full array tra(1,1,1,jn) is initialised

        CALL lbc( tra(1,1,1,jn), 1, 1, 1, 1, jpk, 1 )

#endif
#if defined key_trc_arakawa || defined key_trc_tvd 


C 2.a Time filter for Arakawa or TVD Sheme 
C ----------------------------------------

!$omp   parallel default(none) private(mytid,ji,jj,jk)
!$omp&      shared(jn,neuler,kt,TimeStepStart,jpk,jpj,jpi,trb,trn,tra,tmask)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
      IF( mytid + jn <= jptra ) THEN

        IF ((neuler.eq.0).and.(kt.eq.TimeStepStart)) THEN
            DO jk=1,jpk
              DO jj=1,jpj
                DO ji=1,jpi
                  trb(ji,jj,jk,jn+mytid)=trn(ji,jj,jk,jn+mytid)*tmask(ji,jj,jk)
                END DO
              END DO
            END DO 
        ELSE
            DO jk=1,jpk
             DO jj=1,jpj
                DO ji=1,jpi
                  trb(ji,jj,jk,jn+mytid)=fsass(trb(ji,jj,jk,jn+mytid),
     $                trn(ji,jj,jk,jn+mytid),tra(ji,jj,jk,jn+mytid))
                END DO  
              END DO  
            END DO 
        ENDIF

      END IF
!$omp end parallel

#else

C 2.b TRB for Smolarkiewicz scheme
C -------------------------------
c!$omp   parallel default(none) private(mytid,ji,jj,jk)
c!$omp&      shared(jn,jpk,jpj,jpi,trb,tra)

!#ifdef __OPENMP
!        mytid = omp_get_thread_num()  ! take the thread ID
!#endif
!      IF( mytid + jn <= jptra ) THEN
!
!        DO jk=1,jpk
!          DO jj=1,jpj
!            DO ji=1,jpi
!              trb(ji,jj,jk,jn+mytid)=tra(ji,jj,jk,jn+mytid)
!            END DO
!          END DO
!        END DO
!
!      END IF
c!$omp end parallel
#endif
C
C
C 3. swap of arrays
C -----------------
!$omp   parallel default(none) private(mytid,ji,jj,jk,jw)
!$omp&      shared(jn,jpk,jpj,jpi,trn,trb,tra,tmask) !NWATERPOINTS,WATERPOINTS,JWATERPOINTS,KWATERPOINTS)

#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
      IF( mytid + jn <= jptra ) THEN

C        DO jk = 1,jpk
C          DO jj = 1,jpj
C            DO ji = 1,jpi
C              trn(ji,jj,jk,jn+mytid) = tra(ji,jj,jk,jn+mytid)*tmask(ji,jj,jk)
!CCC F79 We put .le. against .lt. since bfm doesn't like 0 concentration
!              if( tra(ji,jj,jk,jn+mytid).le.0.  ) then
!CCC             WRITE(numout,*) 'conc. exception ( ji , jj , jk , jn)'
!                tra_FN(ji,jj,jk,jn+mytid) = tra_FN(ji,jj,jk,jn+mytid) - tra(ji,jj,jk,jn+mytid) + 0.00000000001
!CCC F79 13 01 2005 WARNING MASS INJECTION
!              trn(ji,jj,jk,jn+mytid) = 0.00000000001
!              trb(ji,jj,jk,jn+mytid) = 0.00000000001
!              end if
C              tra(ji,jj,jk,jn+mytid) = 0.e0 ! tra becomes back trend and not source
C            END DO
C          END DO
C        END DO


        DO jk = 1,jpk
          DO jj = 1,jpj
            DO ji = 1,jpi
!        DO jw = 1, NWATERPOINTS
!            ji = WATERPOINTS(1,jw)
!            jj = WATERPOINTS(2,jw)
!            jk = WATERPOINTS(3,jw)
            trb(ji,jj,jk,jn+mytid) = tra(ji,jj,jk,jn+mytid)
            trn(ji,jj,jk,jn+mytid) = tra(ji,jj,jk,jn+mytid)*tmask(ji,jj,jk)
            tra(ji,jj,jk,jn+mytid) = 0.e0
!         ENDDO
            END DO
          END DO
        END DO



      END IF
!$omp end parallel



C END of tracer slab
C ==================

       END DO TRACER_LOOP


       trcnxtparttime = MPI_WTIME() - trcnxtparttime ! cronometer-stop
       trcnxttottime = trcnxttottime + trcnxtparttime


      END SUBROUTINE trcnxt
