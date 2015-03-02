CC
CC ORCA
CC
CC mpp routines
CC
CC mynode
CC mpparent    --> to be removed
CC mppspawn    --> to be removed
CC mppshmem    --> to be removed
CC mpplnk
CC mpplnk2
CC mpplnks
CC mpprecv
CC mppsend
CC mppwait
CC mppscatter  --> to be removed
CC mppgather   --> to be removed
CC mppiisl     --> to be removed
CC mppimin     --> to be removed
CC mppisum     --> to be removed
CC mpprisl     --> to be removed
CC mpprmax     --> to be removed
CC mpprmin     --> to be removed
CC mpprsum     --> to be removed
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
      CALL mpi_init(ierr)
      CALL mpi_comm_rank(mpi_comm_world,rank,ierr)
      CALL mpi_comm_size(mpi_comm_world,size,ierr)

CC-CC      if(size .NE. jpnij) then
CC-CC         write(*,*) "mynode: Error !!!! size = ", size, " .NE. jpnij = ", jpnij
CC-CC         call mpi_abort(mpi_comm_world, -1, ierr)
CC-CC         stop
CC-CC      endif

      mynode=rank
      RETURN
#  else
      INTEGER mynode
      mynode=0
      RETURN
#endif
      END

      SUBROUTINE mpparent(kparent_tid)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpparent
ccc                     ********************
CCC
CCC  Purpose :
CCC  ---------
CCC    If key_mpp_pvm then call pvmfparent fonction 
CCC    else use an pvmfparent routine for T3E
CCC   (default key or key_mpp_shmem)
CCC    or  only RETURN -1 (key_mpp_mpi)
CCC
CC   Input :
CC   -----
CC      argument                :
CC                   kparent_tid
CC
CC       original  : 96-05 (j. Escobar)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versionsC
C-----------------------------------------------------------------------


c
      USE myalloc
      USE myalloc_mpp

        IMPLICIT NONE



C-----------------------------------------------------------------------
C
      INTEGER kparent_tid
C

#ifdef key_mpp_mpi
C
C MPI VERSION
C
      kparent_tid=-1
C
#endif
      RETURN
      END

      SUBROUTINE mppspawn(cdexec,kmod,cdwhere,kproc,ktids,kinfo)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppspawn
CCC                     ********************
CCC
CCC  Purpose :
CCC  ---------
CCC    If key_mpp_pvm then call pvmfspawn fonction 
CCC    else use an pvmfspawn routine for T3E
CCC   (default key or key_mpp_shmem)
CCC    or  only RETURN -1 (key_mpp_mpi)
CCC
CC   INPUT :
CC   -----
CC      argument                :
CC                   cdexec
CC                   kmod
CC                   cdwhere
CC                   kproc
CC                   ktids
CC                   kinfo
CC
C
C-----------------------------------------------------------------------


c
      USE myalloc
      USE myalloc_mpp


        IMPLICIT NONE


C-----------------------------------------------------------------------
C
      CHARACTER  cdexec(*),cdwhere(*)
      INTEGER kmod,kproc,ktids(*),kinfo


C
#ifdef key_mpp_mpi
C
C MPI VERSION
C
      kinfo=-1
C
#endif
C
C
      RETURN
      END

      SUBROUTINE mppshmem
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppshmem
CCC
CCC  PURPOSE :
CCC  ---------
CCC     Massively parallel processors
CCC     SHMEM ROUTINE
CCC
CC
C-----------------------------------------------------------------------



      USE myalloc
      USE myalloc_mpp

        IMPLICIT NONE


      RETURN
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

      !REAL(8) t3p1(jpi,1,jpk,2)
      !REAL(8) t3p2(jpi,1,jpk,2)
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

#ifdef __OPENMP
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
#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
C      write(*,*) 'omp :: ntids --> ', ntids
#else
      ntids = 1
      mytid = 0
#endif
C
C WRITE(*,*) 'trcnxt:::mpplnk_my:::jptra=',jptra,', jpj=',jpj,', jpi=',jpi, ', jpk=',jpk
C WRITE(*,*) 'trcnxt:::mpplnk_my:::RANK=',rank,', nbondi= ',nbondi,', nbondj= ',nbondj
    
    
    
    
C      trcadvparttime = MPI_WTIME()



!$omp   parallel default(none) private(jn,jk,jj,ji,mytid,iihom,ijhom)
!$omp&      shared(packsize,nbondi,nperio,jpk,jpj,jpi,ptab,jpim1,ktype,nlci,jpreci,nlcj,jprecj)
#ifdef __OPENMP
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
      

#ifdef __OPENMP
      END IF
!$omp    end parallel
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
!$omp   parallel default(none) private(jn,jk,jj,jl,mytid,iihom)
!$omp&      shared(packsize,nbondi,nlci,nreci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab)
#ifdef __OPENMP
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
#ifdef __OPENMP
      END IF
!$omp    end parallel
#else
      END DO PACK_LOOP2
#endif

C        trcadvparttime = MPI_WTIME() - trcadvparttime
C        trcadvtottime = trcadvtottime + trcadvparttime

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
C          if (rank==1) then
C             WRITE(*,*) ':::2 mppsends proceeded'
C         else
C         endif
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
!$omp   parallel default(none) private(jn,jk,jj,jl,mytid,iihom)
!$omp&      shared(packsize,nbondi,nlci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab)
#ifdef __OPENMP
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
#ifdef __OPENMP
      END IF
!$omp    end parallel
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
!$omp   parallel default(none) private(jn,jk,ji,jl,mytid,ijhom)
!$omp&      shared(packsize,nbondj,nlcj,nrecj,jprecj,jpk,jpi,t3sn_my1,t3ns_my1,ptab)
#ifdef __OPENMP
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

#ifdef __OPENMP
      END IF
!$omp    end parallel
#else
      END DO PACK_LOOP4
#endif
C       trcadvparttime = MPI_WTIME() - trcadvparttime
C       trcadvtottime = trcadvtottime + trcadvparttime

C
C 3.2 Migrations
C
C
      imigr=jprecj*jpi*jpk*packsize
C
C      if (rank==1) then
C        WRITE(*,*) ':::mpplnk:::3.1:::nbondj=',nbondj
C     else
C     endif

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

!$omp   parallel default(none) private(jn,jk,ji,jl,mytid,ijhom)
!$omp&      shared(packsize,nbondj,nlcj,jprecj,jpk,jpi,t3sn_my1,t3ns_my1,ptab)
#ifdef __OPENMP
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

#ifdef __OPENMP
      END IF
!$omp    end parallel
#else
      END DO PACK_LOOP5
#endif

C        trcadvparttime = MPI_WTIME() - trcadvparttime
C        trcadvtottime = trcadvtottime + trcadvparttime

C
C
C 4. north fold treatment
C -----------------------
C
C 4.1 treatment without exchange (jpni odd)
C   
C      trcadvparttime = MPI_WTIME()

!$omp   parallel default(none) private(jn,jk,ji,ijt,iju,mytid,iloc)
!$omp&      shared(packsize,npolj,jpiglo,nimpp,ktype,jpk,nlci,zsgn,ptab,nlcj)
#ifdef __OPENMP
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

#ifdef __OPENMP
      END IF
!$omp    end parallel
#else
      END DO PACK_LOOP6
#endif

C        trcadvparttime = MPI_WTIME() - trcadvparttime
C        trcadvtottime = trcadvtottime + trcadvparttime

C
C 4.1 treatment with exchange (jpni greater than 1)
C 
C ... sign ans sort are taken into account in the sender processor
C
C        trcadvparttime = MPI_WTIME()

!$omp   parallel default(none) private(jn,jk,ji,ijt,iju,mytid,iloc)
!$omp&      shared(packsize,npolj,jpiglo,nimpp,nimppt,nono,ktype,jpk,jpi,t3p1_my1,t3p2_my1,zsgn,ptab,nlcj)
#ifdef __OPENMP
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
CCC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                  t3p1(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-2,jk)
CCC                  t3p2(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-1,jk)
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
CCC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                  t3p1(ji,1,jk,1) = zsgn * ptab(iju,nlcj-2,jk)
CCC                  t3p2(ji,1,jk,1) = zsgn * ptab(iju,nlcj-1,jk)
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
CCC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                  t3p1(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-2,jk)
CCC                  t3p2(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-3,jk)
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
CCC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                  t3p1(ji,1,jk,1) = ptab(iju,nlcj-2,jk)
CCC                  t3p2(ji,1,jk,1) = ptab(iju,nlcj-3,jk)
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

#ifdef __OPENMP
      END IF
!$omp    end parallel
#else
      END DO PACK_LOOP7
#endif

C        trcadvparttime = MPI_WTIME() - trcadvparttime
C        trcadvtottime = trcadvtottime + trcadvparttime

C
C 4.2 Migrations
          IF(npolj.eq.3) THEN
C
C
          imigr=jprecj*jpi*jpk*packsize
C     if (rank==1) then
C       WRITE(*,*) ':::mpplnk:::4.2:::imigr=',imigr
C     else
C     endif

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

C        trcadvparttime = MPI_WTIME()
C

!$omp   parallel default(none) private(jn,jk,ji,mytid)
!$omp&      shared(packsize,npolj,ktype,jpk,nlci,ptab,t3p1_my1,nimpp,nlcj,t3p2_my1,jpjglo)
#ifdef __OPENMP
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
#ifdef __OPENMP
      END IF
!$omp    end parallel
#else
      END DO PACK_LOOP8
#endif

C        trcadvparttime = MPI_WTIME() - trcadvparttime
C        trcadvtottime = trcadvtottime + trcadvparttime

C
C
C 5. East and west directions exchange
C ------------------------------------
C
C      trcadvparttime = MPI_WTIME()

!$omp   parallel default(none) private(jn,jl,jk,jj,mytid,iihom)
!$omp&      shared(packsize,npolj,nbondi,nlci,nreci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab)
#ifdef __OPENMP
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
#ifdef __OPENMP
      END IF
!$omp    end parallel
#else
      END DO PACK_LOOP9
#endif
C        trcadvparttime = MPI_WTIME() - trcadvparttime
C        trcadvtottime = trcadvtottime + trcadvparttime

C
C 5.2 Migrations
      IF (npolj.eq.3.or.npolj.eq.4) THEN
C
C
          imigr=jpreci*jpj*jpk*packsize
C     if (rank==1) then
C       WRITE(*,*) ':::mpplnk:::5.2:::imigr=',imigr
C     else
C     endif

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
C        trcadvparttime = MPI_WTIME()

!$omp   parallel default(none) private(jn,jl,jk,jj,mytid,iihom)
!$omp&      shared(packsize,npolj,nbondi,nlci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab)
#ifdef __OPENMP
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
#ifdef __OPENMP
      END IF
!$omp    end parallel
#else
      END DO PACK_LOOP10
#endif

C        trcadvparttime = MPI_WTIME() - trcadvparttime
C        trcadvtottime = trcadvtottime + trcadvparttime

#  else
C
C      No mpp computation
C
#endif
C
C
      RETURN
      END

      SUBROUTINE mpplnk(ptab,ktype,ksgn)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpplnk
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
      REAL(8) ptab(jpi,jpj,jpk)

      REAL(8) t3p1(jpi,1,jpk,2)
      REAL(8) t3p2(jpi,1,jpk,2)

#ifdef key_mpp_mpi

      INTEGER ji,jj,jk,jl
      INTEGER imigr,iihom,ijhom,iloc,ijt,iju
      REAL(8) zsgn
      INTEGER reqs1, reqs2, reqr1, reqr2
!      INTEGER jn

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
C
C 1. standard boundary treatment
C ------------------------------
C
C East-West boundary conditions
C
      IF(nbondi.EQ.2.AND.(nperio.EQ.1.or.nperio.EQ.4)) THEN
C ... cyclic
          DO jk = 1,jpk
            DO jj = 1, jpj
              ptab( 1 ,jj,jk) = ptab(jpim1,jj,jk)
              ptab(jpi,jj,jk) = ptab(  2  ,jj,jk)
            END DO
          END DO
      ELSE
C ... closed
          IF( ktype .NE. 11 .and. ktype .NE. 14 ) Then
              iihom = nlci-jpreci
              DO ji = iihom+1,jpi
                DO jk = 1,jpk
                  DO jj = 1,jpj
                    ptab(ji,jj,jk) = 0.e0
                  END DO
                END DO
              END DO
              IF ( ktype.NE.4  ) THEN
                  DO ji = 1,jpreci
                    DO jk = 1,jpk
                      DO jj = 1,jpj
                        ptab(ji,jj,jk) = 0.e0
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
                ptab(ji,jj,jk) = 0.e0
              END DO
            END DO
          END DO
          IF ( ktype.NE.4 ) THEN
              DO jj = 1,jprecj
                DO jk = 1,jpk
                  DO ji = 1, jpi
                    ptab(ji,jj,jk) = 0.e0
                  END DO
                END DO
              END DO
          ENDIF
      ENDIF
      
      

C 
C 
C 2. East and west directions exchange
C ------------------------------------
C
C 2.1 Read Dirichlet lateral conditions
C
      IF(nbondi.ne.2) THEN
          iihom=nlci-nreci
          DO jl=1,jpreci
            DO jk=1,jpk
              DO jj=1,jpj
                t3ew(jj,jl,jk,1)=ptab(jpreci+jl,jj,jk)
                t3we(jj,jl,jk,1)=ptab(iihom +jl,jj,jk)
              END DO
            END DO
          END DO
      ENDIF
C
C 2.2 Migrations
C

C
      imigr=jpreci*jpj*jpk
C
C     if (rank==1) then
C       WRITE(*,*) ':::mpplnk:::2.2:::nbondi=',nbondi
C     else
C     endif
      IF(nbondi.eq.-1) THEN
          CALL mppsend(2,t3we(1,1,1,1),imigr,noea,0,reqs1)
          CALL mpprecv(1,t3ew(1,1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondi.eq.0) THEN
          CALL mppsend(1,t3ew(1,1,1,1),imigr,nowe,0,reqs1)
          CALL mppsend(2,t3we(1,1,1,1),imigr,noea,0,reqs2)
          CALL mpprecv(1,t3ew(1,1,1,2),imigr,reqr1)
          CALL mpprecv(2,t3we(1,1,1,2),imigr,reqr2)
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondi.eq.1) THEN
          CALL mppsend(1,t3ew(1,1,1,1),imigr,nowe,0,reqs1)
          CALL mpprecv(2,t3we(1,1,1,2),imigr,reqr1)
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
            DO jk=1,jpk
              DO jj=1,jpj
                ptab(jl,jj,jk)=t3we(jj,jl,jk,2)
              END DO
            END DO
          END DO
      ENDIF
C
      IF(nbondi.eq.-1.or.nbondi.eq.0) THEN
          DO jl=1,jpreci
            DO jk=1,jpk
              DO jj=1,jpj
                ptab(iihom+jl,jj,jk)=t3ew(jj,jl,jk,2)
              END DO
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
            DO jk=1,jpk
              DO ji=1,jpi
                t3sn(ji,jl,jk,1)=ptab(ji,ijhom +jl,jk)
                t3ns(ji,jl,jk,1)=ptab(ji,jprecj+jl,jk)
              END DO
            END DO
          END DO
      ENDIF
C
C 3.2 Migrations
C
C
      imigr=jprecj*jpi*jpk
C
C     if (rank==1) then
C       WRITE(*,*) ':::mpplnk:::3.1:::nbondj=',nbondj
C     else
C     endif

      IF(nbondj.eq.-1) THEN
          CALL mppsend(4,t3sn(1,1,1,1),imigr,nono,0,reqs1)
          CALL mpprecv(3,t3ns(1,1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondj.eq.0) THEN
          CALL mppsend(3,t3ns(1,1,1,1),imigr,noso,0,reqs1)
          CALL mppsend(4,t3sn(1,1,1,1),imigr,nono,0,reqs2)
          CALL mpprecv(3,t3ns(1,1,1,2),imigr,reqr1)
          CALL mpprecv(4,t3sn(1,1,1,2),imigr,reqr2)
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondj.eq.1) THEN
          CALL mppsend(3,t3ns(1,1,1,1),imigr,noso,0,reqs1)
          CALL mpprecv(4,t3sn(1,1,1,2),imigr,reqr1)
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
            DO jk=1,jpk
              DO ji=1,jpi
                ptab(ji,jl,jk)=t3sn(ji,jl,jk,2)
              END DO
            END DO
          END DO
      ENDIF
C
      IF(nbondj.eq.0.or.nbondj.eq.-1) THEN
          DO jl=1,jprecj
            DO jk=1,jpk
              DO ji=1,jpi
                ptab(ji,ijhom+jl,jk)=t3ns(ji,jl,jk,2)
              END DO
            END DO
          END DO
      ENDIF
C
C
C 4. north fold treatment
C -----------------------
C
C 4.1 treatment without exchange (jpni odd)
C   
      IF (npolj.eq.4) THEN 
          iloc=jpiglo-2*(nimpp-1)
          IF ( ktype.EQ.1 .OR. ktype.EQ.11 ) THEN
              DO jk = 1, jpk
                DO ji = 2, nlci
                  ijt=iloc-ji+2
                  ptab(ji,nlcj,jk) = zsgn * ptab(ijt,nlcj-2,jk)
                END DO
                DO ji = nlci/2+1, nlci
                  ijt=iloc-ji+2
                  ptab(ji,nlcj-1,jk) = zsgn * ptab(ijt,nlcj-1,jk)
                END DO
              END DO 
          ELSEIF ( ktype.EQ.2 ) THEN
              DO jk = 1, jpk
                DO ji = 1, nlci-1
                  iju=iloc-ji+1
                  ptab(ji,nlcj,jk) = zsgn * ptab(iju,nlcj-2,jk)
                END DO 
                DO ji = nlci/2, nlci-1
                  iju=iloc-ji+1
                  ptab(ji,nlcj-1,jk) = zsgn * ptab(iju,nlcj-1,jk)
                END DO
              END DO
          ELSEIF ( ktype.EQ.3 ) THEN
              DO jk = 1, jpk
                DO ji = 2, nlci
                  ijt=iloc-ji+2
                  ptab(ji,nlcj-1,jk) = zsgn * ptab(ijt,nlcj-2,jk)
                  ptab(ji,nlcj  ,jk) = zsgn * ptab(ijt,nlcj-3,jk)
                END DO 
              END DO
          ELSEIF ( ktype.EQ.4 .OR. ktype.EQ.14 ) THEN
              DO jk = 1, jpk
                DO ji = 1, nlci-1
                  iju=iloc-ji+1
                  ptab(ji,nlcj-1,jk) = ptab(iju,nlcj-2,jk)
                  ptab(ji,nlcj  ,jk) = ptab(iju,nlcj-3,jk)
                END DO
              END DO 
          ENDIF
      ENDIF
C
C 4.1 treatment with exchange (jpni greater than 1)
C 
C ... sign ans sort are taken into account in the sender processor
C
      IF (npolj.eq.3) THEN
          iloc=jpiglo-(nimpp-1+nimppt(nono+1)-1)
          IF ( ktype.EQ.1 .OR. ktype.EQ.11 ) THEN
              DO jk=1,jpk
                DO ji=2,jpi
                  ijt=iloc-ji+2
CCC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                  t3p1(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-2,jk)
CCC                  t3p2(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-1,jk)
                  if(ijt .ge. 1) then
                     t3p1(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-2,jk)
                     t3p2(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-1,jk)
                  else
                     t3p1(ji,1,jk,1) = zsgn * ptab(ijt+jpi,nlcj-2-1,jk)
                     t3p2(ji,1,jk,1) = zsgn * ptab(ijt+jpi,nlcj-1-1,jk)
                  endif
                END DO
              END DO
          ELSEIF ( ktype.EQ.2 ) THEN
              DO jk=1,jpk
                DO ji = 1, jpi-1
                  iju=iloc-ji+1
CCC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                  t3p1(ji,1,jk,1) = zsgn * ptab(iju,nlcj-2,jk)
CCC                  t3p2(ji,1,jk,1) = zsgn * ptab(iju,nlcj-1,jk)
                  if(iju .ge. 1) then
                     t3p1(ji,1,jk,1) = zsgn * ptab(iju,nlcj-2,jk)
                     t3p2(ji,1,jk,1) = zsgn * ptab(iju,nlcj-1,jk)
                  else
                     t3p1(ji,1,jk,1) = zsgn * ptab(iju+jpi,nlcj-2-1,jk)
                     t3p2(ji,1,jk,1) = zsgn * ptab(iju+jpi,nlcj-1-1,jk)
                  endif
                END DO
              END DO
          ELSEIF ( ktype.EQ.3 ) THEN
              DO jk=1,jpk
                DO ji = 2, jpi
                  ijt=iloc-ji+2
CCC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                  t3p1(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-2,jk)
CCC                  t3p2(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-3,jk)
                  if(ijt .ge. 1) then
                     t3p1(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-2,jk)
                     t3p2(ji,1,jk,1) = zsgn * ptab(ijt,nlcj-3,jk)
                  else
                     t3p1(ji,1,jk,1) = zsgn * ptab(ijt+jpi,nlcj-2-1,jk)
                     t3p2(ji,1,jk,1) = zsgn * ptab(ijt+jpi,nlcj-3-1,jk)
                  endif
                END DO
              END DO 
          ELSEIF ( ktype.EQ.4 .OR. ktype.EQ.14 ) THEN
              DO jk=1,jpk
                DO ji = 1, jpi-1
                  iju=iloc-ji+1
CCC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                  t3p1(ji,1,jk,1) = ptab(iju,nlcj-2,jk)
CCC                  t3p2(ji,1,jk,1) = ptab(iju,nlcj-3,jk)
                  if(iju .ge. 1) then
                     t3p1(ji,1,jk,1) = ptab(iju,nlcj-2,jk)
                     t3p2(ji,1,jk,1) = ptab(iju,nlcj-3,jk)
                  else
                     t3p1(ji,1,jk,1) = ptab(iju+jpi,nlcj-2-1,jk)
                     t3p2(ji,1,jk,1) = ptab(iju+jpi,nlcj-3-1,jk)
                  endif
                END DO
              END DO 
          ENDIF
C
C 4.2 Migrations
C

C
          imigr=jprecj*jpi*jpk
C     if (rank==1) then
C       WRITE(*,*) ':::mpplnk:::4.2:::imigr=',imigr
C     else
C     endif

C
          CALL mppsend(3,t3p1(1,1,1,1),imigr,nono,0,reqs1)
          CALL mpprecv(3,t3p1(1,1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
          CALL mppsend(4,t3p2(1,1,1,1),imigr,nono,0,reqs2)
          CALL mpprecv(4,t3p2(1,1,1,2),imigr,reqr2)
          CALL mppwait(reqs2)
          CALL mppwait(reqr2)
C

C
C 4.3 Write north fold conditions
C
          IF ( ktype .EQ. 1  .or.ktype .eq. 11 ) THEN
              DO jk = 1, jpk
                DO ji = 2, nlci
                  ptab(ji,nlcj,jk) = t3p1(ji,1,jk,2)
                END DO
              END DO 
              IF(nimpp+nlcj/2 .gt. jpjglo/2) THEN
                  DO jk = 1, jpk
                    DO ji = 2, nlci
                      ptab(ji,nlcj-1,jk) = t3p2(ji,1,jk,2)
                    END DO
                  END DO
              ENDIF  
          ELSEIF ( ktype.EQ.2 ) THEn
              DO jk = 1, jpk
                DO ji = 1, nlci-1
                  ptab(ji,nlcj,jk) = t3p1(ji,1,jk,2)
                END DO
              END DO
              IF(nimpp+nlcj/2 .gt. jpjglo/2) THEN
                  DO jk = 1, jpk
                    DO ji = 1, nlci-1
                      ptab(ji,nlcj-1,jk) = t3p2(ji,1,jk,2)
                    END DO
                  END DO
              ENDIF
          ELSEIF ( ktype .EQ.3 ) THEN
              DO jk = 1, jpk
                DO ji = 2, nlci
                  ptab(ji,nlcj-1,jk) = t3p1(ji,1,jk,2)
                  ptab(ji,nlcj  ,jk) = t3p2(ji,1,jk,2)
                END DO 
              END DO
          ELSEIF ( ktype .EQ.4 .or. ktype .EQ.14 ) THEN
              DO jk = 1, jpk
                DO ji = 1, nlci-1
                  ptab(ji,nlcj-1,jk) = t3p1(ji,1,jk,2)
                  ptab(ji,nlcj  ,jk) = t3p2(ji,1,jk,2)
                END DO
              END DO
          ENDIF
      ENDIF 
C
C
C 5. East and west directions exchange
C ------------------------------------
C
      IF (npolj.eq.3.or.npolj.eq.4) THEN
C
C 5.1 Read Dirichlet lateral conditions
C
          IF(nbondi.ne.2) THEN
              iihom=nlci-nreci
              DO jl=1,jpreci
                DO jk=1,jpk
                  DO jj=1,jpj
                    t3ew(jj,jl,jk,1)=ptab(jpreci+jl,jj,jk)
                    t3we(jj,jl,jk,1)=ptab(iihom +jl,jj,jk)
                  END DO
                END DO
              END DO
          ENDIF
C
C 5.2 Migrations
C
C
          imigr=jpreci*jpj*jpk
C     if (rank==1) then
C       WRITE(*,*) ':::mpplnk:::5.2:::imigr=',imigr
C     else
C     endif

C
          IF(nbondi.eq.-1) THEN
              CALL mppsend(2,t3we(1,1,1,1),imigr,noea,0,reqs1)
              CALL mpprecv(1,t3ew(1,1,1,2),imigr,reqr1)
              CALL mppwait(reqs1)
              CALL mppwait(reqr1)             
          ELSE IF(nbondi.eq.0) THEN
              CALL mppsend(1,t3ew(1,1,1,1),imigr,nowe,0,reqs1)
              CALL mppsend(2,t3we(1,1,1,1),imigr,noea,0,reqs2)
              CALL mpprecv(1,t3ew(1,1,1,2),imigr,reqr1)
              CALL mpprecv(2,t3we(1,1,1,2),imigr,reqr2)
              CALL mppwait(reqs1)
              CALL mppwait(reqs2)
              CALL mppwait(reqr1)
              CALL mppwait(reqr2)
          ELSE IF(nbondi.eq.1) THEN
              CALL mppsend(1,t3ew(1,1,1,1),imigr,nowe,0,reqs1)
              CALL mpprecv(2,t3we(1,1,1,2),imigr,reqr1)
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
                DO jk=1,jpk
                  DO jj=1,jpj
                    ptab(jl,jj,jk)=t3we(jj,jl,jk,2)
                  END DO
                END DO
              END DO
          ENDIF
C
          IF(nbondi.eq.-1.or.nbondi.eq.0) THEN
              DO jl=1,jpreci
                DO jk=1,jpk
                  DO jj=1,jpj
                    ptab(iihom+jl,jj,jk)=t3ew(jj,jl,jk,2)
                  END DO
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
CC                  at which ptab is defined for 0
CC                                initialization
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
CCC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                t2p1(ji,1,1) = zsgn * ptab(ijt,nlcj-2)
CCC                t2p2(ji,1,1) = zsgn * ptab(ijt,nlcj-1)
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
CCC                write(*,*) "---> ",ji, iju, nlcj
CC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                t2p1(ji,1,1) = zsgn * ptab(iju,nlcj-2)
CCC                t2p2(ji,1,1) = zsgn * ptab(iju,nlcj-1)
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
CC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                t2p1(ji,1,1) = zsgn * ptab(ijt,nlcj-2)
CCC                t2p2(ji,1,1) = zsgn * ptab(ijt,nlcj-3)
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
CC 24/11/2003 Cristiano - To eliminate an out of bound - Is it right ????
CCC Original
CCC                t2p1(ji,1,1) = ptab(iju,nlcj-2)
CCC                t2p2(ji,1,1) = ptab(iju,nlcj-3)
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

      SUBROUTINE mpplnks(ptab)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpplnks
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC      Message passing manadgement for add 2d array local boundary
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
      REAL(8) ptab(jpi,jpj)
      REAL(8) t2p1(jpi,1,2)
!      REAL(8) t2p2(jpi,1,2)

#ifdef key_mpp_mpi


      INTEGER ji,jj,jl
      INTEGER imigr,iihom,ijhom
      INTEGER reqs1, reqs2, reqr1, reqr2

C 1. north fold treatment
C -----------------------
C
C 1.1 treatment without exchange (jpni odd)
C   
      IF (npolj.eq.4) THEN 
          DO ji = 1, nlci
            ptab(ji,nlcj-2)=ptab(ji,nlcj-2)+t2p1(ji,1,1)
          END DO 
      ENDIF
C
C 1.2 treatment with exchange (jpni greater than 1)
C 
      IF (npolj.eq.3) THEN
C
C
          imigr=jprecj*jpi
C
          CALL mppsend(3,t2p1(1,1,1),imigr,nono,0,reqs1)
          CALL mpprecv(3,t2p1(1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
C
C
C Write north fold conditions
C
          DO ji = 1, nlci
            ptab(ji,nlcj-2)=ptab(ji,nlcj-2)+t2p1(ji,1,2)
          END DO 
      ENDIF
C
C
C 2. East and west directions
C ---------------------------
C
C 2.1 Read Dirichlet lateral conditions
C
      iihom=nlci-jpreci
C
      IF(nbondi.ne.2) THEN
          DO jl=1,jpreci
            DO jj=1,jpj
              t2ew(jj,jl,1)=ptab(  jl    ,jj)
              t2we(jj,jl,1)=ptab(iihom+jl,jj)
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
      IF(nbondi.eq.0.or.nbondi.eq.1) THEN
          DO jl=1,jpreci
            DO jj=1,jpj
              ptab(jpreci+jl,jj)=ptab(jpreci+jl,jj)+t2we(jj,jl,2)
            END DO
          END DO
      ENDIF
C
      IF(nbondi.eq.-1.or.nbondi.eq.0) THEN
          iihom=nlci-nreci
          DO jl=1,jpreci
            DO jj=1,jpj
              ptab(iihom +jl,jj)=ptab(iihom +jl,jj)+t2ew(jj,jl,2)
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
      ijhom=nlcj-jprecj
C
      IF(nbondj.ne.2) THEN
          DO jl=1,jprecj
            DO ji=1,jpi
              t2sn(ji,jl,1)=ptab(ji,ijhom+jl)
              t2ns(ji,jl,1)=ptab(ji,   jl   )
            END DO
          END DO
      ENDIF
C
C 3.2 Migrations
C
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
      IF(nbondj.eq.0.or.nbondj.eq.1) THEN
          DO jl=1,jprecj
            DO ji=1,jpi
              ptab(ji,jprecj+jl)=ptab(ji,jprecj+jl)+t2sn(ji,jl,2)
            END DO
          END DO
      ENDIF
C
      IF(nbondj.eq.0.or.nbondj.eq.-1) THEN
          ijhom=nlcj-nrecj
          DO jl=1,jprecj
            DO ji=1,jpi
              ptab(ji,ijhom +jl)=ptab(ji,ijhom +jl)+t2ns(ji,jl,2)
            END DO
          END DO
      ENDIF
C
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

      SUBROUTINE mppgather(ptab,kk,kp,pio)
CCC--------------------------------------------------------------Enigmail-------
CCC
CCC                       routine mppgather
CCC                     *********************
CCC
CCC  Purpose :
CCC  ---------
CCC     Transfert between a local subdomain array and a work array 
CCC     which is distributed following the vertical level.
CCC
CC   Method :
CC   -------
CC
CC   Input :
CC   -----
CC      argument
CC           ptab   : subdomain array input
CC           kk     : vertical level
CC           kp     : record length
CC
CC   Output :
CC   ------
CC      argument 
CC           pio    : output array
CC
CC   Workspace :
CC   ---------
CC      imess, ic
CC
CC   EXTERNAL :
CC   --------
CC      mppsend mppsync mpprecv
CC   or shmem_put barrier
CC   or mpi_gather
CC
CC   REFERENCES :                 no
CC   ----------
CC
CC   MODIFICATIONS:
CC   --------------
CC       original  : 93-09 (M. Imbard)
CC       additions : 96-05 (j. Escobar)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
CC----------------------------------------------------------------------


c
      USE myalloc
      USE myalloc_mpp


        IMPLICIT NONE


CCC--------------------------------------------------------------------
CCC  OPA8, IDRIS-LODYC (15/11/96)
CCC--------------------------------------------------------------------
C
      INTEGER kk,kp
      REAL(8) ptab(jpi,jpj),pio(jpi,jpj,jpnij)
C
#ifdef key_mpp_mpi




C

C
      INTEGER itaille,ierror
C
      itaille=jpi*jpj
      CALL mpi_gather(ptab,itaille,mpi_real8,pio,itaille
     $    ,mpi_real8,kp,mpi_comm_world,ierror) 

#endif
C
      RETURN
      END

      SUBROUTINE mppscatter(pio,kk,kp,ptab)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppscatter
CCC                     **********************
CCC
CCC  Purpose :
CCC  ---------
CCC     Transfert between awork array which is distributed following 
CCC     the vertical level and the local subdomain array.
CCC
CC   Method :
CC   -------
CC
CC   Input :
CC   -----
CC      argument
CC           pio    -> output array
CC           kk     -> process number
CC           kp     -> Tag (not used with MPI
CC
CC   Output :
CC   ------
CC      argument 
CC           ptab   : subdomain array input
CC
CC   Workspace :
CC   ---------
CC      imess, ic
CC
CC   External :
CC   --------
CC      pvmfscatter mppsync
CC   or shmem_get barrier
CC
CC   References :                 no
CC   ----------
CC
CC   Modifications:
CC   --------------
CC       original  : 93-09 (M. Imbard)
CC       additions : 96-05 (j. Escobar)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
CC----------------------------------------------------------------------


c
      USE myalloc
      USE myalloc_mpp


        IMPLICIT NONE


CCC--------------------------------------------------------------------
CCC  OPA8, IDRIS-LODYC (15/11/96)
CCC--------------------------------------------------------------------
C
      INTEGER kk,kp
      REAL(8) ptab(jpi,jpj),pio(jpi,jpj,jpnij)
C
#ifdef key_mpp_mpi




C

C
      INTEGER itaille, ierror
C
      itaille=jpi*jpj

      CALL mpi_scatter(pio,itaille,mpi_real8,ptab,itaille,
     $     mpi_real8,kp,mpi_comm_world,ierror)
C

#endif
C
C
      RETURN
      END

      SUBROUTINE mppiisl(ktab,kdim,kwork)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppiisl
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Massively parallel processors
CCC     Find the  non zero value
CCC
CC   Input :
CC   -----
CC      argument                :
CC                   ktab
CC                   kdim
CC                   kwork
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
      INTEGER kdim
      INTEGER ktab(kdim),kwork(kdim)
C
#ifdef key_mpp_mpi




C

C
      EXTERNAL lc_isl
      LOGICAL lcommute
      INTEGER mpi_isl,ierror
C
      lcommute=.true.
      CALL mpi_op_create(lc_isl,lcommute,mpi_isl
     $    ,ierror)
CCC 20/06/2003 Probably it should be in this manner ....
CCC      CALL mpi_allreduce(ptab,pwork,kdim,mpi_integer
      CALL mpi_allreduce(ktab,kwork,kdim,mpi_integer
     $    ,mpi_isl,mpi_comm_world,ierror)
      ktab=kwork
C

#endif
      END

!      SUBROUTINE PvmIsl2(kdtatyp,kx,ky,kdim,knfo)
!      INTEGER kx(kdim),ky(kdim)
!      INTEGER kdim,knfo,kdtatyp,ji
!      DO ji=1,kdim
!        IF (ky(ji) .NE. 0) kx(ji)=ky(ji)
!      END DO
!      RETURN
!      END

      SUBROUTINE mppimin(ktab,kdim,kwork)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppimin
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Massively parallel processors
CCC     Find minimum value in an integer layout array
CCC
CC   Input :
CC   -----
CC                   ktab   ->   input array
CC                   kdim   ->   dimension of the array ktab and kwork
CC                   kwork->   work space array used by mpi 
CC
CC   Modifications:
CC   --------------
CC       original  : 93-09 (M. Imbard)
CC       additions : 96-05 (j. Escobar)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
CC----------------------------------------------------------------------
CC parameters and commons
CC ======================


c
      USE myalloc
      USE myalloc_mpp


        IMPLICIT NONE


CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER kdim
      INTEGER ktab(kdim),kwork(kdim)
C
#ifdef key_mpp_mpi




C
C
      INTEGER ierror
C
      CALL mpi_allreduce(ktab,kwork,kdim,mpi_integer
     $    ,mpi_min,mpi_comm_world,ierror)

      ktab=kwork
C

#endif
C
      RETURN
      END

      SUBROUTINE mppisum(ktab,kdim,kwork)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mppisum
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Massively parallel processors
CCC     Global integer sum
CCC
CC   Input :
CC   -----
CC      argument                
CC                   ktab   ->   input array
CC                   kdim   ->   dimension of the array ktab and kwork
CC                   kwork  ->   work space array used by mpi 
CC   Output :
CC   -------
CC       argument 
CC                   ktab 
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
      INTEGER kdim
      INTEGER ktab(kdim),kwork(kdim)
C
#ifdef key_mpp_mpi





      INTEGER ierror
C
      CALL mpi_allreduce(ktab,kwork,kdim,mpi_integer
     $    ,mpi_sum,mpi_comm_world,ierror)

      ktab=kwork
C

#endif
      RETURN
      END

      SUBROUTINE mpprisl(ptab,kdim,pwork)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpprisl
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Massively parallel processors
CCC     Find the  non zero island barotropic stream function value
CCC
CC   Input :
CC   -----
CC      argument                :
CC                   ptab
CC                   kdim
CC                   pwork
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
      INTEGER kdim
      REAL(8) ptab(kdim),pwork(kdim)
C
#ifdef key_mpp_mpi




C

C
      EXTERNAL lc_isl
      LOGICAL lcommute
      INTEGER mpi_isl,ierror
C
      lcommute=.true.
      CALL mpi_op_create(lc_isl,lcommute,mpi_isl
     $    ,ierror)
      CALL mpi_allreduce(ptab,pwork,kdim,mpi_real8
     $    ,mpi_isl,mpi_comm_world,ierror)
      ptab=pwork
C

#endif
      RETURN
      END
      SUBROUTINE PvmIsl(kdtatyp,px,py,kdim,knfo)
      REAL(8) px(kdim),py(kdim)
      integer kdim,knfo,kdtatyp,ji
      DO ji=1,kdim
        IF (py(ji) .NE. 0.) px(ji)=py(ji)
      END DO
      RETURN
      END
      FUNCTION LC_ISL(py,px,kdim,kdtatyp)
      REAL(8) px(kdim),py(kdim)
      Integer kdim,kdtatyp,ji
      Integer LC_ISL
      DO ji=1,kdim
        IF (py(ji) .NE. 0.) px(ji)=py(ji)
      END DO
      LC_ISL=0
      RETURN
      END

      SUBROUTINE mpprmax(ptab,kdim,pwork)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpprmax
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Massively parallel processors
CCC     Maximum
CC
CC   Input :
CC   -----
CC      argument                :
CC                   ptab
CC                   kdim
CC                   pwork
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
      INTEGER kdim
      REAL(8) ptab(kdim),pwork(kdim)
C
#ifdef key_mpp_mpi




C

C
      INTEGER ierror
C
      CALL mpi_allreduce(ptab,pwork,kdim,mpi_real8
     $    ,mpi_max,mpi_comm_world,ierror)
      ptab=pwork
C

#endif
      RETURN
      END

      SUBROUTINE mpprmin(ptab,kdim,pwork)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpprmin
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Massively parallel processors
CCC     Minimum
CC
CC   Input :
CC   -----
CC      argument                :
CC                   ptab
CC                   kdim
CC                   pwork
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
      INTEGER kdim
      REAL(8) ptab(kdim),pwork(kdim)
C
#ifdef key_mpp_mpi




C

C
      INTEGER ierror
C
      CALL mpi_allreduce(ptab,pwork,kdim,mpi_real8
     $    ,mpi_min,mpi_comm_world,ierror)
      ptab=pwork
C

#endif
      RETURN
      END

      SUBROUTINE mpprsum(ptab,kdim,pwork)
CCC---------------------------------------------------------------------
CCC
CCC                       routine mpprsum
CCC                     *******************
CCC
CCC  Purpose :
CCC  ---------
CCC     Massively parallel processors
CCC     Global sum
CC
CC   Input :
CC   -----
CC      argument                :
CC                   ptab   ->   input array
CC                   kdim   ->   dimension of the array ktab and kwork
CC                   pwork  ->   work space array used by mpi 
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
      INTEGER kdim
      REAL(8) ptab(kdim),pwork(kdim)
C
#ifdef key_mpp_mpi




C
C
      INTEGER ierror
C
      CALL mpi_allreduce(ptab,pwork,kdim,mpi_real8
     $    ,mpi_sum,mpi_comm_world,ierror)
      ptab=pwork
C

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
C
#ifdef key_mpp_mpi




C

C
      INTEGER ierror
C
      CALL mpi_barrier(mpi_comm_world,ierror)
C

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

C
C 2. Mpp synchroneus
C ------------------
C

#ifdef key_mpp_mpi

      CALL mpi_finalize(info)

#endif
C
      RETURN
      END
