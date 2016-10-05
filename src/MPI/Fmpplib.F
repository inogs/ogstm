
CC mpp routines
CC
CC mynode
CC mpplnk_my
CC mpprecv
CC mppsend
CC mppwait
CC mppsync
CC mppstop
CC
      FUNCTION mynode()
CCC---------------------------------------------------------------------
CCC
CCC                       routine mynode
CCC                     ******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Massively parallel processors
CCC     Find processor unit
CCC
CC   Input :
CC   -----
CC      argument                :
CC
CC   Modifications:
CC   --------------
CC       original  : 93-09 (M. Imbard)
CC       additions : 96-05 (j. Escobar)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
C-----------------------------------------------------------------------


       USE myalloc
       USE myalloc_mpp

C-----------------------------------------------------------------------

        IMPLICIT NONE
C
#ifdef key_mpp_mpi
C
C MPI VERSION
C
      INTEGER mynode,ierr
C         -------------
C         Enroll in MPI
C         -------------
C

      CALL mpi_comm_rank(mpi_comm_world,rank,ierr)
      CALL mpi_comm_size(mpi_comm_world,mpi_size_comm,ierr)

      mynode=rank
      RETURN
#  else
      INTEGER mynode
      mynode=0
      RETURN
#endif
      END

      SUBROUTINE mpplnk_my(ptab,packsize,ktype,ksgn)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpplnk_my
CCC                     ******************
CCC
CCC  Purpose :
CCC  ---------
CCC      Message passing manadgement
CCC
CC   Method :
CC   -------
CC       Use mppsend and mpprecv function for passing mask between
CC       processors following neighboring subdomains.
CC
CC   Input :
CC   -----
CC      argument
CC              ptab            : variable array
CC              ktype           : define the nature of the grid-point
CC                  at which ptab is defined for 0
CC                                initialization
CC                  = 1 ,  T- and W-points
CC                  = 2 ,  U-point
CC                  = 3 ,  V-point
CC                  = 4 ,  F-point
CC                                = 11,  T-point only fold treatment
CC                                = 14,  F-point only fold treatment
CC        ksgn        : control of the sign change
CC                  = 0 , the sign is modified following
CC                  the type of b.c. used
CC                  = 1 , the sign of the field is un-
CC                  changed at the boundaries
CC      common
CC            /COMDOM/ : domain parameters
CC                    nlci   : first dimension of the local subdomain
CC                    nlcj   : second dimension of the local subdomain
CC                    nbondi : mark for "east-west local boundary"
CC                    nbondj : mark for "north-south local boundary"
CC                    noea   : number for local neighboring processors
CC                    nowe   : number for local neighboring processors
CC                    noso   : number for local neighboring processors
CC                    nono   : number for local neighboring processors
CC            /COMMPP/ : massively parallel processors
CC                    t3ew() : message passing arrays east-west
CC                    t3we() : message passing arrays west-east
CC                    t3ns() : message passing arrays north-south
CC                    t3sn() : message passing arrays south-north
CC
CC   Output :
CC   ------
CC      common
CC            /COMMPP/ : massively parallel processors
CC                    t3ew() : message passing arrays east-west
CC                    t3we() : message passing arrays west-east
CC                    t3ns() : message passing arrays north-south
CC                    t3sn() : message passing arrays south-north
CC   Workspace :
CC   ---------
CC             ji,jj,jk,jl,imigr,iihom,ijhom
CC
CC   External :
CC   --------
CC             mppsend,mpprecv
CC       or    shmem_put barrier shmem_udcflush
CC
CC
CC   References :                 no
CC   ----------
CC
CC   Modifications:
CC   --------------
CC       original  : 94-11 (M. Guyon)
CC       additions : 95-04 (j. Escobar, M. Imbard)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
CC----------------------------------------------------------------------

      USE myalloc
      USE myalloc_mpp

        IMPLICIT NONE


CC----------------------------------------------------------------------
C
      INTEGER ktype, ksgn
      INTEGER packsize
      REAL(8) ptab(jpi,jpj,jpk,packsize)

      REAL(8) t3p1_my1(jpi,1,jpk,packsize,2)
      REAL(8) t3p2_my1(jpi,1,jpk,packsize,2)

#ifdef key_mpp_mpi

      INTEGER ji,jj,jk,jl
      INTEGER imigr,iihom,ijhom,iloc,ijt,iju
      REAL(8) zsgn
      INTEGER reqs1, reqs2, reqr1, reqr2
      INTEGER jn
! omp variables
      INTEGER :: mytid, ntids!, itid

#ifdef __OPENMP11
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


C
C 0. Initialization
C -----------------
C
C Sign setting
C ...
      IF (ksgn.EQ.0) THEN
          zsgn = -1.
      ELSE
          zsgn =  1.
      ENDIF
C OPENMP settings
#ifdef __OPENMP11
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif


C      trcadvparttime = MPI_WTIME()



!!!!!$omp   parallel default(none) private(jn,jk,jj,ji,mytid,iihom,ijhom)
!!!!!$omp&      shared(packsize,nbondi,nperio,jpk,jpj,jpi,ptab,jpim1,ktype,nlci,jpreci,nlcj,jprecj)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP1: DO jn=1,packsize
#endif
C 1. standard boundary treatment
C ------------------------------
C
C East-West boundary conditions
C
      IF(nbondi.EQ.2.AND.(nperio.EQ.1.or.nperio.EQ.4)) THEN
C ... cyclic
          DO jk = 1,jpk
            DO jj = 1, jpj
              ptab( 1 ,jj,jk,jn+mytid) = ptab(jpim1,jj,jk,jn+mytid)
              ptab(jpi,jj,jk,jn+mytid) = ptab(  2  ,jj,jk,jn+mytid)
            END DO
          END DO
      ELSE
C ... closed
          IF( ktype .NE. 11 .and. ktype .NE. 14 ) Then
              iihom = nlci-jpreci
              DO ji = iihom+1,jpi
                DO jk = 1,jpk
                  DO jj = 1,jpj
                    ptab(ji,jj,jk,jn+mytid) = 0.e0
                  END DO
                END DO
              END DO
              IF ( ktype.NE.4  ) THEN
                  DO ji = 1,jpreci
                    DO jk = 1,jpk
                      DO jj = 1,jpj
                        ptab(ji,jj,jk,jn+mytid) = 0.e0
                      END DO
                    END DO
                  END DO
              ENDIF
          ENDIF
      ENDIF
C
C North-South boundary conditions
C
      IF( ktype .NE. 11 .and. ktype .NE. 14 ) THEN
          ijhom = nlcj-jprecj
          DO jj = ijhom+1,jpj
            DO jk = 1,jpk
              DO ji = 1,jpi
                ptab(ji,jj,jk,jn+mytid) = 0.e0
              END DO
            END DO
          END DO
          IF ( ktype.NE.4 ) THEN
              DO jj = 1,jprecj
                DO jk = 1,jpk
                  DO ji = 1, jpi
                    ptab(ji,jj,jk,jn+mytid) = 0.e0
                  END DO
                END DO
              END DO
          ENDIF
      ENDIF


#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP1
#endif
C
C
C 2. East and west directions exchange
C ------------------------------------
C
C 2.1 Read Dirichlet lateral conditions
C
!!!!!$omp   parallel default(none) private(jn,jk,jj,jl,mytid,iihom)
!!!!!$omp&      shared(packsize,nbondi,nlci,nreci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP2: DO jn=1,packsize
#endif

      IF(nbondi.ne.2) THEN
          iihom=nlci-nreci
          DO jl=1,jpreci
            DO jk=1,jpk
              DO jj=1,jpj
                t3ew_my1(jj,jl,jk,jn+mytid,1)=ptab(jpreci+jl,jj,jk,jn+mytid)
                t3we_my1(jj,jl,jk,jn+mytid,1)=ptab(iihom +jl,jj,jk,jn+mytid)
              END DO
            END DO
          END DO
      ENDIF
#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP2
#endif

C
C 2.2 Migrations
C
C
      imigr=jpreci*jpj*jpk*packsize
C
      IF(nbondi.eq.-1) THEN
          CALL mppsend(2,t3we_my1(1,1,1,1,1),imigr,noea,0,reqs1)
          CALL mpprecv(1,t3ew_my1(1,1,1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondi.eq.0) THEN
          CALL mppsend(1,t3ew_my1(1,1,1,1,1),imigr,nowe,0,reqs1)
          CALL mppsend(2,t3we_my1(1,1,1,1,1),imigr,noea,0,reqs2)
          CALL mpprecv(1,t3ew_my1(1,1,1,1,2),imigr,reqr1)
          CALL mpprecv(2,t3we_my1(1,1,1,1,2),imigr,reqr2)
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondi.eq.1) THEN
          CALL mppsend(1,t3ew_my1(1,1,1,1,1),imigr,nowe,0,reqs1)
          CALL mpprecv(2,t3we_my1(1,1,1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ENDIF
C

C
C 2.3 Write Dirichlet lateral conditions
C
C      trcadvparttime = MPI_WTIME()
!!!!!$omp   parallel default(none) private(jn,jk,jj,jl,mytid,iihom)
!!!!!$omp&      shared(packsize,nbondi,nlci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP3: DO jn=1,packsize
#endif
      iihom=nlci-jpreci
      IF(nbondi.eq.0.or.nbondi.eq.1) THEN
C
          DO jl=1,jpreci
            DO jk=1,jpk
              DO jj=1,jpj
                ptab(jl,jj,jk,jn+mytid)=t3we_my1(jj,jl,jk,jn+mytid,2)
              END DO
            END DO
          END DO
      ENDIF
C
      IF(nbondi.eq.-1.or.nbondi.eq.0) THEN
          DO jl=1,jpreci
            DO jk=1,jpk
              DO jj=1,jpj
                ptab(iihom+jl,jj,jk,jn+mytid)=t3ew_my1(jj,jl,jk,jn+mytid,2)
              END DO
            END DO
          END DO
      ENDIF
#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP3
#endif

C        trcadvparttime = MPI_WTIME() - trcadvparttime
C        trcadvtottime = trcadvtottime + trcadvparttime

C
C
C 3. North and south directions
C -----------------------------
C
C 3.1 Read Dirichlet lateral conditions
C
C        trcadvparttime = MPI_WTIME()
!!!!!$omp   parallel default(none) private(jn,jk,ji,jl,mytid,ijhom)
!!!!!$omp&      shared(packsize,nbondj,nlcj,nrecj,jprecj,jpk,jpi,t3sn_my1,t3ns_my1,ptab)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP4: DO jn=1,packsize
#endif

      IF(nbondj.ne.2) THEN
          ijhom=nlcj-nrecj
C
          DO jl=1,jprecj
            DO jk=1,jpk
              DO ji=1,jpi
                t3sn_my1(ji,jl,jk,jn+mytid,1)=ptab(ji,ijhom +jl,jk,jn+mytid)
                t3ns_my1(ji,jl,jk,jn+mytid,1)=ptab(ji,jprecj+jl,jk,jn+mytid)
              END DO
            END DO
          END DO
      ENDIF

#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP4
#endif

C
C 3.2 Migrations
C
C
      imigr=jprecj*jpi*jpk*packsize

      IF(nbondj.eq.-1) THEN
          CALL mppsend(4,t3sn_my1(1,1,1,1,1),imigr,nono,0,reqs1)
          CALL mpprecv(3,t3ns_my1(1,1,1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondj.eq.0) THEN
          CALL mppsend(3,t3ns_my1(1,1,1,1,1),imigr,noso,0,reqs1)
          CALL mppsend(4,t3sn_my1(1,1,1,1,1),imigr,nono,0,reqs2)
          CALL mpprecv(3,t3ns_my1(1,1,1,1,2),imigr,reqr1)
          CALL mpprecv(4,t3sn_my1(1,1,1,1,2),imigr,reqr2)
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondj.eq.1) THEN
          CALL mppsend(3,t3ns_my1(1,1,1,1,1),imigr,noso,0,reqs1)
          CALL mpprecv(4,t3sn_my1(1,1,1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ENDIF
C
C
C 3.3 Write Dirichlet lateral conditions
C
C        trcadvparttime = MPI_WTIME()

!!!!!$omp   parallel default(none) private(jn,jk,ji,jl,mytid,ijhom)
!!!!!$omp&      shared(packsize,nbondj,nlcj,jprecj,jpk,jpi,t3sn_my1,t3ns_my1,ptab)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP5: DO jn=1,packsize
#endif
      ijhom=nlcj-jprecj
      IF(nbondj.eq.0.or.nbondj.eq.1) THEN
          DO jl=1,jprecj
            DO jk=1,jpk
              DO ji=1,jpi
                ptab(ji,jl,jk,jn+mytid)=t3sn_my1(ji,jl,jk,jn+mytid,2)
              END DO
            END DO
          END DO
      ENDIF
C
      IF(nbondj.eq.0.or.nbondj.eq.-1) THEN
          DO jl=1,jprecj
            DO jk=1,jpk
              DO ji=1,jpi
                ptab(ji,ijhom+jl,jk,jn+mytid)=t3ns_my1(ji,jl,jk,jn+mytid,2)
              END DO
            END DO
          END DO
      ENDIF

#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP5
#endif

C
C
C 4. north fold treatment
C -----------------------
C
C 4.1 treatment without exchange (jpni odd)
C
C      trcadvparttime = MPI_WTIME()

!!!!!$omp   parallel default(none) private(jn,jk,ji,ijt,iju,mytid,iloc)
!!!!!$omp&      shared(packsize,npolj,jpiglo,nimpp,ktype,jpk,nlci,zsgn,ptab,nlcj)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP6: DO jn=1,packsize
#endif
      IF (npolj.eq.4) THEN
          iloc=jpiglo-2*(nimpp-1)
          IF ( ktype.EQ.1 .OR. ktype.EQ.11 ) THEN
              DO jk = 1, jpk
                DO ji = 2, nlci
                  ijt=iloc-ji+2
                  ptab(ji,nlcj,jk,jn+mytid) = zsgn * ptab(ijt,nlcj-2,jk,jn+mytid)
                END DO
                DO ji = nlci/2+1, nlci
                  ijt=iloc-ji+2
                  ptab(ji,nlcj-1,jk,jn+mytid) = zsgn * ptab(ijt,nlcj-1,jk,jn+mytid)
                END DO
              END DO
          ELSEIF ( ktype.EQ.2 ) THEN
              DO jk = 1, jpk
                DO ji = 1, nlci-1
                  iju=iloc-ji+1
                  ptab(ji,nlcj,jk,jn+mytid) = zsgn * ptab(iju,nlcj-2,jk,jn+mytid)
                END DO
                DO ji = nlci/2, nlci-1
                  iju=iloc-ji+1
                  ptab(ji,nlcj-1,jk,jn+mytid) = zsgn * ptab(iju,nlcj-1,jk,jn+mytid)
                END DO
              END DO
          ELSEIF ( ktype.EQ.3 ) THEN
              DO jk = 1, jpk
                DO ji = 2, nlci
                  ijt=iloc-ji+2
                  ptab(ji,nlcj-1,jk,jn+mytid) = zsgn * ptab(ijt,nlcj-2,jk,jn+mytid)
                  ptab(ji,nlcj  ,jk,jn+mytid) = zsgn * ptab(ijt,nlcj-3,jk,jn+mytid)
                END DO
              END DO
          ELSEIF ( ktype.EQ.4 .OR. ktype.EQ.14 ) THEN
              DO jk = 1, jpk
                DO ji = 1, nlci-1
                  iju=iloc-ji+1
                  ptab(ji,nlcj-1,jk,jn+mytid) = ptab(iju,nlcj-2,jk,jn+mytid)
                  ptab(ji,nlcj  ,jk,jn+mytid) = ptab(iju,nlcj-3,jk,jn+mytid)
                END DO
              END DO
          ENDIF
      ENDIF

#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP6
#endif

C
C 4.1 treatment with exchange (jpni greater than 1)
C
C ... sign ans sort are taken into account in the sender processor
C

!!!!!$omp   parallel default(none) private(jn,jk,ji,ijt,iju,mytid,iloc)
!!!!!$omp&      shared(packsize,npolj,jpiglo,nimpp,nimppt,nono,ktype,jpk,jpi,t3p1_my1,t3p2_my1,zsgn,ptab,nlcj)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP7: DO jn=1,packsize
#endif
      IF (npolj.eq.3) THEN
          iloc=jpiglo-(nimpp-1+nimppt(nono+1)-1)
          IF ( ktype.EQ.1 .OR. ktype.EQ.11 ) THEN
              DO jk=1,jpk
                DO ji=2,jpi
                  ijt=iloc-ji+2
                  if(ijt .ge. 1) then
                     t3p1_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(ijt,nlcj-2,jk,jn+mytid)
                     t3p2_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(ijt,nlcj-1,jk,jn+mytid)
                  else
                     t3p1_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(ijt+jpi,nlcj-2-1,jk,jn+mytid)
                     t3p2_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(ijt+jpi,nlcj-1-1,jk,jn+mytid)
                  endif
                END DO
              END DO
          ELSEIF ( ktype.EQ.2 ) THEN
              DO jk=1,jpk
                DO ji = 1, jpi-1
                  iju=iloc-ji+1
                  if(iju .ge. 1) then
                     t3p1_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(iju,nlcj-2,jk,jn+mytid)
                     t3p2_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(iju,nlcj-1,jk,jn+mytid)
                  else
                     t3p1_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(iju+jpi,nlcj-2-1,jk,jn+mytid)
                     t3p2_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(iju+jpi,nlcj-1-1,jk,jn+mytid)
                  endif
                END DO
              END DO
          ELSEIF ( ktype.EQ.3 ) THEN
              DO jk=1,jpk
                DO ji = 2, jpi
                  ijt=iloc-ji+2
                  if(ijt .ge. 1) then
                     t3p1_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(ijt,nlcj-2,jk,jn+mytid)
                     t3p2_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(ijt,nlcj-3,jk,jn+mytid)
                  else
                     t3p1_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(ijt+jpi,nlcj-2-1,jk,jn+mytid)
                     t3p2_my1(ji,1,jk,jn+mytid,1) = zsgn * ptab(ijt+jpi,nlcj-3-1,jk,jn+mytid)
                  endif
                END DO
              END DO
          ELSEIF ( ktype.EQ.4 .OR. ktype.EQ.14 ) THEN
              DO jk=1,jpk
                DO ji = 1, jpi-1
                  iju=iloc-ji+1
                  if(iju .ge. 1) then
                     t3p1_my1(ji,1,jk,jn+mytid,1) = ptab(iju,nlcj-2,jk,jn+mytid)
                     t3p2_my1(ji,1,jk,jn+mytid,1) = ptab(iju,nlcj-3,jk,jn+mytid)
                  else
                     t3p1_my1(ji,1,jk,jn+mytid,1) = ptab(iju+jpi,nlcj-2-1,jk,jn+mytid)
                     t3p2_my1(ji,1,jk,jn+mytid,1) = ptab(iju+jpi,nlcj-3-1,jk,jn+mytid)
                  endif
                END DO
              END DO
          ENDIF
        ENDIF

#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP7
#endif


C
C 4.2 Migrations
          IF(npolj.eq.3) THEN
C
C
          imigr=jprecj*jpi*jpk*packsize


C
          CALL mppsend(3,t3p1_my1(1,1,1,1,1),imigr,nono,0,reqs1)
          CALL mpprecv(3,t3p1_my1(1,1,1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
          CALL mppsend(4,t3p2_my1(1,1,1,1,1),imigr,nono,0,reqs2)
          CALL mpprecv(4,t3p2_my1(1,1,1,1,2),imigr,reqr2)
          CALL mppwait(reqs2)
          CALL mppwait(reqr2)
C
          ENDIF
C
C 4.3 Write north fold conditions
C

!!!!!$omp   parallel default(none) private(jn,jk,ji,mytid)
!!!!!$omp&      shared(packsize,npolj,ktype,jpk,nlci,ptab,t3p1_my1,nimpp,nlcj,t3p2_my1,jpjglo)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP8: DO jn=1,packsize
#endif
         IF(npolj.eq.3) THEN
          IF ( ktype .EQ. 1  .or.ktype .eq. 11 ) THEN
              DO jk = 1, jpk
                DO ji = 2, nlci
                  ptab(ji,nlcj,jk,jn+mytid) = t3p1_my1(ji,1,jk,jn+mytid,2)
                END DO
              END DO
              IF(nimpp+nlcj/2 .gt. jpjglo/2) THEN
                  DO jk = 1, jpk
                    DO ji = 2, nlci
                      ptab(ji,nlcj-1,jk,jn+mytid) = t3p2_my1(ji,1,jk,jn+mytid,2)
                    END DO
                  END DO
              ENDIF
          ELSEIF ( ktype.EQ.2 ) THEn
              DO jk = 1, jpk
                DO ji = 1, nlci-1
                  ptab(ji,nlcj,jk,jn+mytid) = t3p1_my1(ji,1,jk,jn+mytid,2)
                END DO
              END DO
              IF(nimpp+nlcj/2 .gt. jpjglo/2) THEN
                  DO jk = 1, jpk
                    DO ji = 1, nlci-1
                      ptab(ji,nlcj-1,jk,jn+mytid) = t3p2_my1(ji,1,jk,jn+mytid,2)
                    END DO
                  END DO
              ENDIF
          ELSEIF ( ktype .EQ.3 ) THEN
              DO jk = 1, jpk
                DO ji = 2, nlci
                  ptab(ji,nlcj-1,jk,jn+mytid) = t3p1_my1(ji,1,jk,jn+mytid,2)
                  ptab(ji,nlcj  ,jk,jn+mytid) = t3p2_my1(ji,1,jk,jn+mytid,2)
                END DO
              END DO
          ELSEIF ( ktype .EQ.4 .or. ktype .EQ.14 ) THEN
              DO jk = 1, jpk
                DO ji = 1, nlci-1
                  ptab(ji,nlcj-1,jk,jn+mytid) = t3p1_my1(ji,1,jk,jn+mytid,2)
                  ptab(ji,nlcj  ,jk,jn+mytid) = t3p2_my1(ji,1,jk,jn+mytid,2)
                END DO
              END DO
          ENDIF
      ENDIF
#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP8
#endif

C
C
C 5. East and west directions exchange
C ------------------------------------
C

!!!!!$omp   parallel default(none) private(jn,jl,jk,jj,mytid,iihom)
!!!!!$omp&      shared(packsize,npolj,nbondi,nlci,nreci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP9: DO jn=1,packsize
#endif
      IF (npolj.eq.3.or.npolj.eq.4) THEN
C
C 5.1 Read Dirichlet lateral conditions
C
          IF(nbondi.ne.2) THEN
              iihom=nlci-nreci
              DO jl=1,jpreci
                DO jk=1,jpk
C Check the following
                  DO jj=1,jpj
                    t3ew_my1(jj,jl,jk,jn+mytid,1)=ptab(jpreci+jl,jj,jk,jn+mytid)
                    t3we_my1(jj,jl,jk,jn+mytid,1)=ptab(iihom +jl,jj,jk,jn+mytid)
                  END DO
                END DO
              END DO
          ENDIF
         ENDIF
#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP9
#endif

C
C 5.2 Migrations
      IF (npolj.eq.3.or.npolj.eq.4) THEN
C
C
          imigr=jpreci*jpj*jpk*packsize

C
          IF(nbondi.eq.-1) THEN
              CALL mppsend(2,t3we_my1(1,1,1,1,1),imigr,noea,0,reqs1)
              CALL mpprecv(1,t3ew_my1(1,1,1,1,2),imigr,reqr1)
              CALL mppwait(reqs1)
              CALL mppwait(reqr1)
          ELSE IF(nbondi.eq.0) THEN
              CALL mppsend(1,t3ew_my1(1,1,1,1,1),imigr,nowe,0,reqs1)
              CALL mppsend(2,t3we_my1(1,1,1,1,1),imigr,noea,0,reqs2)
              CALL mpprecv(1,t3ew_my1(1,1,1,1,2),imigr,reqr1)
              CALL mpprecv(2,t3we_my1(1,1,1,1,2),imigr,reqr2)
              CALL mppwait(reqs1)
              CALL mppwait(reqs2)
              CALL mppwait(reqr1)
              CALL mppwait(reqr2)
          ELSE IF(nbondi.eq.1) THEN
              CALL mppsend(1,t3ew_my1(1,1,1,1,1),imigr,nowe,0,reqs1)
              CALL mpprecv(2,t3we_my1(1,1,1,1,2),imigr,reqr1)
              CALL mppwait(reqs1)
              CALL mppwait(reqr1)
          ENDIF
C

         ENDIF
C
C 5.3 Write Dirichlet lateral conditions
C

!!!!!$omp   parallel default(none) private(jn,jl,jk,jj,mytid,iihom)
!!!!!$omp&      shared(packsize,npolj,nbondi,nlci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab)
#ifdef __OPENMP11
        mytid = omp_get_thread_num()  ! take the thread ID
        jn=1
        IF(mytid +1 <= packsize) THEN
#else
      PACK_LOOP10: DO jn=1,packsize
#endif
          IF (npolj.eq.3.or.npolj.eq.4) THEN
          iihom=nlci-jpreci
          IF(nbondi.eq.0.or.nbondi.eq.1) THEN
C
              DO jl=1,jpreci
                DO jk=1,jpk
                  DO jj=1,jpj
                    ptab(jl,jj,jk,jn+mytid)=t3we_my1(jj,jl,jk,jn+mytid,2)
                  END DO
                END DO
              END DO
          ENDIF
C
          IF(nbondi.eq.-1.or.nbondi.eq.0) THEN
              DO jl=1,jpreci
                DO jk=1,jpk
                  DO jj=1,jpj
                    ptab(iihom+jl,jj,jk,jn+mytid)=t3ew_my1(jj,jl,jk,jn+mytid,2)
                  END DO
                END DO
              END DO
          ENDIF
      ENDIF
#ifdef __OPENMP11
      END IF
!!!!!$omp    end parallel
#else
      END DO PACK_LOOP10
#endif

#  else
C
C      No mpp computation
C
#endif
C
C
      RETURN
      END


      SUBROUTINE mpplnk2(ptab,ktype,ksgn)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpplnk2
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC      Message passing manadgement for 2d array
CCC
CC   Method :
CC   -------
CC       Use mppsend and mpprecv function for passing mask between
CC       processors following neighboring subdomains.
CC
CC   Input :
CC   -----
CC      argument
CC              ptab            : variable array
CC              ktype           : define the nature of the grid-point
CC                                at which ptab is defined for 0
CC                                initialization
CC                                initialization
CC                                = 1 ,  T- and W-points
CC                                = 2 ,  U-point
CC                                = 3 ,  V-point
CC                                = 4 ,  F-point
CC                                = 11,  T-point only fold treatment
CC                                = 14,  F-point only fold treatment
CC              ksgn            : control of the sign change
CC                                = 0 , the sign is modified following
CC                                the type of b.c. used
CC                                = 1 , the sign of the field is un-
CC                                changed at the boundaries
CC      common
CC            /COMDOM/ : domain parameters
CC                    nlci   : first dimension of the local subdomain
CC                    nlcj   : second dimension of the local subdomain
CC                    nbondi : mark for "east-west local boundary"
CC                    nbondj : mark for "north-south local boundary"
CC                    noea   : number for local neighboring processors
CC                    nowe   : number for local neighboring processors
CC                    noso   : number for local neighboring processors
CC                    nono   : number for local neighboring processors
CC            /COMMPP/ : massively parallel processors
CC                    t2ew() : message passing arrays east-west
CC                    t2we() : message passing arrays west-east
CC                    t2ns() : message passing arrays north-south
CC                    t2sn() : message passing arrays south-north
CC
CC   Output :
CC   ------
CC      common
CC            /COMMPP/ : massively parallel processors
CC                    t2ew() : message passing arrays east-west
CC                    t2we() : message passing arrays west-east
CC                    t2ns() : message passing arrays north-south
CC                    t2sn() : message passing arrays south-north
CC   Workspace :
CC   ---------
CC      local
CC             ji,jj,jl,imigr,iihom,ijhom
CC
CC   External :
CC   --------
CC             mppsend,mpprecv
CC       or    shmem_put barrier shmem_udcflush
CC
CC   References :                 no
CC   ----------
CC
CC   Modifications:
CC   --------------
CC       original  : 94-11 (M. Guyon)
CC       additions : 95-04 (j. Escobar, M. Imbard)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
CC----------------------------------------------------------------------

      USE myalloc
      USE myalloc_mpp

        IMPLICIT NONE

CC----------------------------------------------------------------------
C
      INTEGER ktype, ksgn
      REAL(8) ptab(jpi,jpj)
      REAL(8) t2p1(jpi,1,2)
      REAL(8) t2p2(jpi,1,2)
#ifdef key_mpp_mpi

      INTEGER ji,jj,jl
      INTEGER imigr,iihom,ijhom,iloc,ijt,iju
      REAL(8) zsgn
      INTEGER reqs1, reqs2, reqr1, reqr2
C
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
C 0. Initialization
C -----------------
C
C Sign setting
C ...
      IF (ksgn.EQ.0) THEN
          zsgn = -1.
      ELSE
          zsgn =  1.
      ENDIF
C
C 1. standard boundary treatment
C ------------------------------
C
C East-West boundary conditions
C
      IF(nbondi.EQ.2.AND.(nperio.EQ.1.or.nperio.EQ.4)) THEN
C ... cyclic

          DO jj = 1, jpj
            ptab( 1 ,jj) = ptab(jpim1,jj)
            ptab(jpi,jj) = ptab(  2  ,jj)
          END DO
      ELSE
C ... closed
          IF( ktype .NE. 11 .and. ktype .NE. 14 ) Then
              iihom = nlci-jpreci
              DO ji = iihom+1,jpi
                DO jj = 1,jpj
                  ptab(ji,jj) = 0.e0
                END DO
              END DO
              IF ( ktype.NE.4  ) THEN
                  DO ji = 1,jpreci
                    DO jj = 1,jpj
                      ptab(ji,jj) = 0.e0
                    END DO
                  END DO
              ENDIF
          ENDIF
      ENDIF
C
C North-South boundary conditions
C
      IF( ktype .NE. 11 .and. ktype .NE. 14 ) THEN
          ijhom = nlcj-jprecj
          DO jj = ijhom+1,jpj
            DO ji = 1,jpi
              ptab(ji,jj) = 0.e0
            END DO
          END DO
          IF ( ktype.NE.4 ) THEN
              DO jj = 1,jprecj
                DO ji = 1, jpi
                  ptab(ji,jj) = 0.e0
                END DO
              END DO
          ENDIF
      ENDIF
C
C
C 2. East and west directions
C ---------------------------
C
C 2.1 Read Dirichlet lateral conditions
C
      IF(nbondi.ne.2) THEN
          iihom=nlci-nreci
C
          DO jl=1,jpreci
            DO jj=1,jpj
              t2ew(jj,jl,1)=ptab(jpreci+jl,jj)
              t2we(jj,jl,1)=ptab(iihom +jl,jj)
            END DO
          END DO
      ENDIF
C
C 2.2 Migrations
C
C
      imigr=jpreci*jpj
C
      IF(nbondi.eq.-1) THEN
          CALL mppsend(2,t2we(1,1,1),imigr,noea,0,reqs1)
          CALL mpprecv(1,t2ew(1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondi.eq.0) THEN
          CALL mppsend(1,t2ew(1,1,1),imigr,nowe,0,reqs1)
          CALL mppsend(2,t2we(1,1,1),imigr,noea,0,reqs2)
          CALL mpprecv(1,t2ew(1,1,2),imigr,reqr1)
          CALL mpprecv(2,t2we(1,1,2),imigr,reqr2)
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondi.eq.1) THEN
          CALL mppsend(1,t2ew(1,1,1),imigr,nowe,0,reqs1)
          CALL mpprecv(2,t2we(1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ENDIF
C
C
C 2.3 Write Dirichlet lateral conditions
C
      iihom=nlci-jpreci
      IF(nbondi.eq.0.or.nbondi.eq.1) THEN
C
          DO jl=1,jpreci
            DO jj=1,jpj
              ptab(jl,jj)=t2we(jj,jl,2)
            END DO
          END DO
      ENDIF
C
      IF(nbondi.eq.-1.or.nbondi.eq.0) THEN
          DO jl=1,jpreci
            DO jj=1,jpj
              ptab(iihom+jl,jj)=t2ew(jj,jl,2)
            END DO
          END DO
      ENDIF
C
C
C 3. North and south directions
C -----------------------------
C
C 3.1 Read Dirichlet lateral conditions
C
      IF(nbondj.ne.2) THEN
          ijhom=nlcj-nrecj
C
          DO jl=1,jprecj
            DO ji=1,jpi
              t2sn(ji,jl,1)=ptab(ji,ijhom +jl)
              t2ns(ji,jl,1)=ptab(ji,jprecj+jl)
            END DO
          END DO
      ENDIF
C
C 3.2 Migrations
C
C
C MPI VERSION
C
      imigr=jprecj*jpi
C
      IF(nbondj.eq.-1) THEN
          CALL mppsend(4,t2sn(1,1,1),imigr,nono,0,reqs1)
          CALL mpprecv(3,t2ns(1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondj.eq.0) THEN
          CALL mppsend(3,t2ns(1,1,1),imigr,noso,0,reqs1)
          CALL mppsend(4,t2sn(1,1,1),imigr,nono,0,reqs2)
          CALL mpprecv(3,t2ns(1,1,2),imigr,reqr1)
          CALL mpprecv(4,t2sn(1,1,2),imigr,reqr2)
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondj.eq.1) THEN
          CALL mppsend(3,t2ns(1,1,1),imigr,noso,0,reqs1)
          CALL mpprecv(4,t2sn(1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ENDIF
C
C
C 3.3 Write Dirichlet lateral conditions
C
      ijhom=nlcj-jprecj
      IF(nbondj.eq.0.or.nbondj.eq.1) THEN
          DO jl=1,jprecj
            DO ji=1,jpi
              ptab(ji,jl)=t2sn(ji,jl,2)
            END DO
          END DO
      ENDIF
C
      IF(nbondj.eq.0.or.nbondj.eq.-1) THEN
          DO jl=1,jprecj
            DO ji=1,jpi
              ptab(ji,ijhom+jl)=t2ns(ji,jl,2)
            END DO
          END DO
      ENDIF
C
C 4. north fold treatment
C -----------------------
C
C 4.1 treatment without exchange (jpni odd)
C
      IF (npolj.eq.4) THEN
          iloc=jpiglo-2*(nimpp-1)
          IF ( ktype.EQ.1 .OR. ktype.EQ.11 ) THEN
              DO ji = 2, nlci
                ijt=iloc-ji+2
                ptab(ji,nlcj) = zsgn * ptab(ijt,nlcj-2)
              END DO
              DO ji = nlci/2+1, nlci
                ijt=iloc-ji+2
                ptab(ji,nlcj-1) = zsgn * ptab(ijt,nlcj-1)
              END DO
          ELSEIF ( ktype.EQ.2 ) THEN
              DO ji = 1, nlci-1
                iju=iloc-ji+1
                ptab(ji,nlcj) = zsgn * ptab(iju,nlcj-2)
              END DO
              DO ji = nlci/2, nlci-1
                iju=iloc-ji+1
                ptab(ji,nlcj-1) = zsgn * ptab(iju,nlcj-1)
              END DO
          ELSEIF ( ktype.EQ.3 ) THEN
              DO ji = 2, nlci
                ijt=iloc-ji+2
                ptab(ji,nlcj-1) = zsgn * ptab(ijt,nlcj-2)
                ptab(ji,nlcj  ) = zsgn * ptab(ijt,nlcj-3)
              END DO
          ELSEIF ( ktype.EQ.4 .OR. ktype.EQ.14 ) THEN
              DO ji = 1, nlci-1
                iju=iloc-ji+1
                ptab(ji,nlcj-1) = ptab(iju,nlcj-2)
                ptab(ji,nlcj  ) = ptab(iju,nlcj-3)
              END DO
          ENDIF
      ENDIF
C
C 4.1 treatment with exchange (jpni greater than 1)
C
C ... sign and sort are taken into account in the sender processor
C
      IF (npolj.eq.3) THEN
          iloc=jpiglo-(nimpp-1+nimppt(nono+1)-1)
          IF ( ktype.EQ.1 .OR. ktype.EQ.11 ) THEN
              DO ji=2,jpi
                ijt=iloc-ji+2
                if(ijt .ge. 1) then
                   t2p1(ji,1,1) = zsgn * ptab(ijt,nlcj-2)
                   t2p2(ji,1,1) = zsgn * ptab(ijt,nlcj-1)
                else
                   t2p1(ji,1,1) = zsgn * ptab(ijt+jpi,nlcj-2-1)
                   t2p2(ji,1,1) = zsgn * ptab(ijt+jpi,nlcj-1-1)
                endif
              END DO
          ELSEIF ( ktype.EQ.2 ) THEN
              DO ji = 1, jpi-1
                iju=iloc-ji+1
                if(iju .ge. 1) then
                   t2p1(ji,1,1) = zsgn * ptab(iju,nlcj-2)
                   t2p2(ji,1,1) = zsgn * ptab(iju,nlcj-1)
                else
                   t2p1(ji,1,1) = zsgn * ptab(iju+jpi,nlcj-2-1)
                   t2p2(ji,1,1) = zsgn * ptab(iju+jpi,nlcj-1-1)
                endif
              END DO
          ELSEIF ( ktype.EQ.3 ) THEN
              DO ji = 2, jpi
                ijt=iloc-ji+2
                if(ijt .ge. 1) then
                   t2p1(ji,1,1) = zsgn * ptab(ijt,nlcj-2)
                   t2p2(ji,1,1) = zsgn * ptab(ijt,nlcj-3)
                else
                   t2p1(ji,1,1) = zsgn * ptab(ijt+jpi,nlcj-2-1)
                   t2p2(ji,1,1) = zsgn * ptab(ijt+jpi,nlcj-3-1)
                endif
              END DO
          ELSEIF ( ktype.EQ.4 .OR. ktype.EQ.14 ) THEN
              DO ji = 1, jpi-1
                iju=iloc-ji+1
                if(iju .ge. 1) then
                   t2p1(ji,1,1) = ptab(iju,nlcj-2)
                   t2p2(ji,1,1) = ptab(iju,nlcj-3)
                else
                   t2p1(ji,1,1) = ptab(iju+jpi,nlcj-2-1)
                   t2p2(ji,1,1) = ptab(iju+jpi,nlcj-3-1)
                endif
              END DO
          ENDIF
C
C 4.2 Migrations
C
C
C MPI VERSION
C
          imigr=jprecj*jpi
C
          CALL mppsend(3,t2p1(1,1,1),imigr,nono,0,reqs1)
          CALL mpprecv(3,t2p1(1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
          CALL mppsend(4,t2p2(1,1,1),imigr,nono,0,reqs2)
          CALL mpprecv(4,t2p2(1,1,2),imigr,reqr2)
          CALL mppwait(reqs2)
          CALL mppwait(reqr2)
C
C
C 4.3 Write north fold conditions
C
          IF ( ktype .EQ. 1 .or .ktype .eq. 11 ) THEN
              DO ji = 2, nlci
                ptab(ji,nlcj) = t2p1(ji,1,2)
              END DO
              IF(nimpp+nlcj/2 .gt. jpjglo/2) THEN
                  DO ji = 2, nlci
                    ptab(ji,nlcj-1) = t2p2(ji,1,2)
                  END DO
              ENDIF
          ELSEIF ( ktype.EQ.2 ) THEN
              DO ji = 1, nlci-1
                ptab(ji,nlcj) = t2p1(ji,1,2)
              END DO
              IF(nimpp+nlcj/2 .gt. jpjglo/2) THEN
                  DO ji = 1, nlci-1
                    ptab(ji,nlcj-1) = t2p2(ji,1,2)
                  END DO
              ENDIF
          ELSEIF ( ktype .EQ.3 ) THEN
              DO ji = 2, nlci
                ptab(ji,nlcj-1) = t2p1(ji,1,2)
                ptab(ji,nlcj  ) = t2p2(ji,1,2)
              END DO
          ELSEIF ( ktype .EQ.4 .or. ktype .EQ.14 ) THEN
              DO ji = 1, nlci-1
                ptab(ji,nlcj-1) = t2p1(ji,1,2)
                ptab(ji,nlcj  ) = t2p2(ji,1,2)
              END DO
          ENDIF
      ENDIF
C
C
C 5. East and west directions
C ---------------------------
C
      IF (npolj.eq.3.or.npolj.eq.4) THEN
C
C 5.1 Read Dirichlet lateral conditions
C
          IF(nbondi.ne.2) THEN
              iihom=nlci-nreci
C
              DO jl=1,jpreci
                DO jj=1,jpj
                  t2ew(jj,jl,1)=ptab(jpreci+jl,jj)
                  t2we(jj,jl,1)=ptab(iihom +jl,jj)
                END DO
              END DO
          ENDIF
C
C 5.2 Migrations
C
C
C MPI VERSION
C
          imigr=jpreci*jpj
C
          IF(nbondi.eq.-1) THEN
              CALL mppsend(2,t2we(1,1,1),imigr,noea,0,reqs1)
              CALL mpprecv(1,t2ew(1,1,2),imigr,reqr1)
              CALL mppwait(reqs1)
              CALL mppwait(reqr1)
          ELSE IF(nbondi.eq.0) THEN
              CALL mppsend(1,t2ew(1,1,1),imigr,nowe,0,reqs1)
              CALL mppsend(2,t2we(1,1,1),imigr,noea,0,reqs2)
              CALL mpprecv(1,t2ew(1,1,2),imigr,reqr1)
              CALL mpprecv(2,t2we(1,1,2),imigr,reqr2)
              CALL mppwait(reqs1)
              CALL mppwait(reqs2)
              CALL mppwait(reqr1)
              CALL mppwait(reqr2)
          ELSE IF(nbondi.eq.1) THEN
              CALL mppsend(1,t2ew(1,1,1),imigr,nowe,0,reqs1)
              CALL mpprecv(2,t2we(1,1,2),imigr,reqr1)
              CALL mppwait(reqs1)
              CALL mppwait(reqr1)
          ENDIF
C
C
C 5.3 Write Dirichlet lateral conditions
C
          iihom=nlci-jpreci
          IF(nbondi.eq.0.or.nbondi.eq.1) THEN
C
              DO jl=1,jpreci
                DO jj=1,jpj
                  ptab(jl,jj)=t2we(jj,jl,2)
                END DO
              END DO
          ENDIF
C
          IF(nbondi.eq.-1.or.nbondi.eq.0) THEN
              DO jl=1,jpreci
                DO jj=1,jpj
                  ptab(iihom+jl,jj)=t2ew(jj,jl,2)
                END DO
              END DO
          ENDIF
      ENDIF
#  else
C
C      No mpp computation
C
#endif
C
C
      RETURN
      END



      SUBROUTINE mppsend(ktyp,pmess,kbytes,kdest,kid,ireqsend)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppsend
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Send messag passing array
CC
CC   Input :
CC   -----
CC      argument                :
CC                   ktyp   -> Tag of the message
CC                   pmess  -> array of real(8) to send
CC                   kbytes -> size of pmess in real(8)
CC                   kdest  -> receive process number
CC                   kid    _> ? (note used)
CC
CC   Modifications:
CC   --------------
CC       original  : 93-09 (M. Imbard)
CC       additions : 96-05 (j. Escobar)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
C-----------------------------------------------------------------------


c
      USE myalloc
      USE myalloc_mpp


        IMPLICIT NONE


C-----------------------------------------------------------------------
C
      REAL(8) pmess(*)
      INTEGER kbytes,kdest,ktyp,kid, ireqsend
C
#ifdef key_mpp_mpi




C
      INTEGER iflag
!      INTEGER itid_dest,info
C
C      write(*,*)  "#### kbytes = ", kbytes
      CALL mpi_isend(pmess,kbytes,mpi_real8,kdest,ktyp,
     $    mpi_comm_world,ireqsend,iflag)
C

#endif
      RETURN
      END

      SUBROUTINE mpprecv(ktyp,pmess,kbytes,ireqrecv)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpprecv
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Receive messag passing array
CC
CC   Input :
CC   -----
CC      argument
CC                   ktyp    -> Tag of the recevied message
CC                   pmess   -> array of real(8)
CC                   kbytes  -> suze of the array pmess


CC
CC   Modifications:
CC   --------------
CC       original  : 93-09 (M. Imbard)
CC       additions : 96-05 (j. Escobar)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
C-----------------------------------------------------------------------


c
      USE myalloc
      USE myalloc_mpp


        IMPLICIT NONE


C-----------------------------------------------------------------------
C
      REAL(8) pmess(*)
      INTEGER   kbytes,ktyp, ireqrecv
C
#ifdef key_mpp_mpi

C
C
C MPI VERSION
C
!      INTEGER istatus(mpi_status_size)
      INTEGER iflag
C
      CALL mpi_irecv(pmess,kbytes,mpi_real8,mpi_any_source,ktyp,mpi_comm_world,ireqrecv,iflag)
c

#endif
      RETURN
      END

      SUBROUTINE mppwait(req)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppwait
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Wait message passing isend/irecv
CC
CC   Input :
CC   -----
CC      argument
C-----------------------------------------------------------------------

      USE myalloc
      USE myalloc_mpp
        IMPLICIT NONE


C-----------------------------------------------------------------------


      integer req
      INTEGER istatus(mpi_status_size), ierr


#ifdef key_mpp_mpi
      call MPI_WAIT(req, istatus, ierr)
#endif
      RETURN
      END

      SUBROUTINE mppsync()
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppsync
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Massively parallel processors, synchroneous
CC
CC   Modifications:
CC   --------------
CC       original  : 93-09 (M. Imbard)
CC       additions : 96-05 (j. Escobar)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
C-----------------------------------------------------------------------


c
      USE myalloc
      USE myalloc_mpp


        IMPLICIT NONE
C-----------------------------------------------------------------------

#ifdef key_mpp_mpi

      INTEGER ierror

      CALL mpi_barrier(mpi_comm_world,ierror)


#endif
      RETURN
      END

      SUBROUTINE mppstop
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppstop
CCC                     *******************
CCC
CCC  purpose :
CCC  --------
CCC     Stop massilively parallel processors method
CC

      USE myalloc
      USE myalloc_mpp
      IMPLICIT NONE
CC local declarations
CC ==================
      INTEGER info
#ifdef key_mpp_mpi

      CALL mppsync



#endif

      RETURN
      END
