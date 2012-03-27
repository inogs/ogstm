CC $Header: /cvsroot/opatm-bfm/opa_model/OPA/Fmpplib.F,v 1.2 2009-09-11 09:20:55 cvsogs01 Exp $
CC
CC ORCA
CC
CC mpp routines
CC
CC mynode
CC mpparent
CC mppspawn
CC mppshmem
CC mpplnk
CC mpplnk2
CC mpplnks
CC mpprecv
CC mppsend
CC mppwait
CC mppscatter
CC mppgather
CC mppiisl
CC mppimin
CC mppisum
CC mpprisl
CC mpprmax
CC mpprmin
CC mpprsum
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

c#include "parameter.h"
c#include "common.h"
       USE myalloc
       USE myalloc_mpp

C-----------------------------------------------------------------------
C
CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

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
C
C PVM or SHMEM VERSION
C
      INTEGER mynode
      INTEGER imypid, imyhost
      INTEGER ji,info
      INTEGER iparent_tid
C
      IF(npvm_ipas .NE. nprocmax) THEN
C
C         ---   first passage in mynode
C
C         -------------
C         enroll in pvm
C         -------------
C
          CALL pvmfmytid( npvm_mytid )
          IF( mynode_print .NE. 0) THEN
              WRITE(nummpp,*) 'mynode, npvm_ipas=',npvm_ipas,
     $            ' nprocmax=',nprocmax
              WRITE(nummpp,*) 'mynode, npvm_mytid=',npvm_mytid,
     $            ' after pvmfmytid'
          END IF
C
C         ---------------------------------------------------------------
C         find out IF i am parent or child spawned processes have parents
C         ---------------------------------------------------------------
C
          CALL mpparent(iparent_tid)
          IF( mynode_print .NE. 0) THEN
              WRITE(nummpp,*) 'mynode, npvm_mytid=',npvm_mytid,
     $            ' after mpparent, npvm_tids(0)=',
     $            npvm_tids(0),' iparent_tid=',
     $            iparent_tid
          END IF
          IF(iparent_tid .LT. 0 )  THEN
              WRITE(nummpp,*) 'mynode, npvm_mytid=',npvm_mytid,
     $            ' after mpparent, npvm_tids(0)=',
     $            npvm_tids(0),' iparent_tid=',
     $            iparent_tid
              npvm_tids(0) = npvm_mytid
              npvm_me = 0
              IF(ndim_mpp .GT. nprocmax) THEN
                  WRITE(nummpp,*) 'npvm_mytid=',npvm_mytid,' too great'
                  STOP  ' mynode '
              ELSE
                  npvm_nproc = ndim_mpp
              END IF
C
C             -------------------------
C             start up copies of myself
C             -------------------------
C
              IF(npvm_nproc .GT. 1) THEN
                  CALL mppspawn(executable,pvmdefault,'*'
     $                ,npvm_nproc-1,npvm_tids(1),info)
                  IF(info .NE. npvm_nproc-1) THEN
                      WRITE(nummpp,*) 'mynode.f, problem in spawn '
     $                    ,' info=',info
     $                    ,' executable=',executable
                      STOP
                  END IF
                  IF(mynode_print .NE. 0) THEN
                      WRITE(nummpp,*) 'mynode, npvm_mytid=',npvm_mytid,
     $                    ' maitre=',executable
     $                    ,' info=',info
     $                    ,' npvm_nproc=',npvm_nproc
                      WRITE(nummpp,*) 'mynode, npvm_mytid=',npvm_mytid,
     $                    ' npvm_tids '
     $                    ,(npvm_tids(ji),ji=0,npvm_nproc-1)
                  ENDIF
C
C                 ---------------------------
C                 multicast tids array to children
C                 ---------------------------
C
                  CALL pvmfinitsend( pvmdefault, info )
                  CALL pvmfpack(jpvmint,npvm_nproc,1,1,info)
                  CALL pvmfpack(jpvmint,npvm_tids,npvm_nproc,1,info)
                  CALL pvmfmcast(npvm_nproc-1,npvm_tids(1),10,info)
              END IF
          ELSE
C
C             ---------------------------------
C             receive the tids array and set me
C             ---------------------------------
C
              IF(mynode_print .NE. 0) THEN
                  WRITE(nummpp,*) 'mynode, npvm_mytid=',npvm_mytid,
     $                ' pvmfrecv'
              ENDIF
              CALL pvmfrecv( iparent_tid, 10, info )
              IF(mynode_print .NE. 0) THEN
                  WRITE(nummpp,*) 'mynode, npvm_mytid=',npvm_mytid,
     $                " fin pvmfrecv"
              ENDIF
              CALL pvmfunpack(jpvmint,npvm_nproc,1,1,info)
              CALL pvmfunpack(jpvmint,npvm_tids,npvm_nproc,1,info)
              IF(mynode_print.ne.0) THEN
                  WRITE(nummpp,*) 'mynode, npvm_mytid=',npvm_mytid,
     $                ' esclave='
     $                ,executable,' info=',info
     $                ,' npvm_nproc=',npvm_nproc
                  WRITE(nummpp,*) 'mynode, npvm_mytid=',npvm_mytid,
     $                'npvm_tids'
     $                ,(npvm_tids(ji),ji=0,npvm_nproc-1)
              ENDIF
              DO ji=0, npvm_nproc-1
                IF(npvm_mytid .EQ. npvm_tids(ji) ) npvm_me = ji
              END DO   
          ENDIF
C
C         ------------------------------------------------------------
C         all nproc tasks are equal now
C         and can address each other by tids(0) thru tids(nproc-1)
C         for each process me => process number [0-(nproc-1)]
C         ------------------------------------------------------------
C
          CALL pvmfjoingroup ("bidon", info)
          CALL pvmfbarrier ("bidon", npvm_nproc, info)
          DO ji=0, npvm_nproc-1
            IF(ji .EQ. npvm_me) THEN
                CALL pvmfjoingroup (opaall, npvm_inum)
                IF(npvm_inum .NE. npvm_me) WRITE(nummpp,*) 'mynode',
     $              ' not arrived in the good order for opaall'
            END IF
            CALL pvmfbarrier("bidon",npvm_nproc,info)
          END DO
          CALL pvmfbarrier(opaall,npvm_nproc,info)
C
      ELSE
C
C         ---   other passage in mynode
C
      ENDIF

      npvm_ipas = nprocmax
      mynode    = npvm_me
      imypid    = npvm_mytid
      imyhost   = npvm_tids(0)
      IF(mynode_print .NE. 0) THEN
          WRITE(nummpp,*)'mynode, npvm_mytid=',npvm_mytid
     $        ,' npvm_me=',npvm_me
     $        ,' npvm_nproc=',npvm_nproc
     $        ,' npvm_ipas=',npvm_ipas
      END IF
      RETURN
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp
CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed


C-----------------------------------------------------------------------
C
      INTEGER kparent_tid
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

#ifdef key_mpp_pvm
C
C PVM VERSION
C
      CALL pvmfparent(kparent_tid)
C
#elif defined key_mpp_mpi
C
C MPI VERSION
C
      kparent_tid=-1
C
#else
C
C SHMEM VERSION or PVM onto T3E
C
      INTEGER it3d_my_pe
      INTEGER leadz
      INTEGER ji,info
C
      CALL pvmfmytid(nt3d_mytid)
      CALL pvmfgetpe(nt3d_mytid,it3d_my_pe)
      IF(mpparent_print .NE. 0) THEN
          WRITE(nummpp,*) 'mpparent, nt3d_mytid=',nt3d_mytid
     $        ,' it3d_my_pe=',it3d_my_pe
      END IF
      IF(it3d_my_pe .EQ. 0) THEN
C
C-----------------------------------------------------------------c
C     process = 0 => receive other tids                           c
C-----------------------------------------------------------------c
C
          kparent_tid = -1
          IF(mpparent_print .NE. 0) THEN
              WRITE(nummpp,*) 'mpparent, nt3d_mytid=',nt3d_mytid
     $            ,' kparent_tid=',kparent_tid
          END IF
C
C          --- END receive dimension ---
C
          IF(ndim_mpp .GT. nprocmax) THEN
              WRITE(nummpp,*) 'mytid=',nt3d_mytid,' too great'
              STOP  ' mpparent '
          ELSE
              nt3d_nproc =  ndim_mpp
          END IF
          IF(mpparent_print .NE. 0) THEN
              WRITE(nummpp,*) 'mpparent, nt3d_mytid=',nt3d_mytid
     $            ,' nt3d_nproc=',nt3d_nproc
          END IF
C
C-------- receive tids from others process --------
C
          DO ji=1,nt3d_nproc-1
            CALL pvmfrecv( ji , 100, info )
            CALL pvmfunpack(jpvmint,nt3d_tids(ji),1,1,info)
            IF(mpparent_print .NE. 0) THEN
                WRITE(nummpp,*) 'mpparent, nt3d_mytid=',nt3d_mytid
     $              ,' receive=',nt3d_tids(ji),' from=',ji
            END IF
          END DO
          nt3d_tids(0)=nt3d_mytid
          IF(mpparent_print .NE. 0) THEN
              WRITE(nummpp,*) 'mpparent, nt3d_mytid=',nt3d_mytid
     $            ,' nt3d_tids(ji) =',(nt3d_tids(ji),
     $            ji=0,nt3d_nproc-1)
          END IF
          IF(mpparent_print .NE. 0) THEN
              WRITE(nummpp,*) 'mpparent, nt3d_mytid=',nt3d_mytid
     $            ,' kparent_tid=',kparent_tid
          END IF

      ELSE
C
C-----------------------------------------------------------------c
C     process <> 0 => send  other tids                            c
C-----------------------------------------------------------------c
C
          kparent_tid = 0
          CALL pvmfinitsend( pvmdataraw, info )
          CALL pvmfpack(jpvmint,nt3d_mytid,1,1,info)
          CALL pvmfsend(kparent_tid,100,info)
      END IF
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
      CHARACTER  cdexec(*),cdwhere(*)
      INTEGER kmod,kproc,ktids(*),kinfo
CC
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_pvm
C
C PVM VERSION
C
      CALL pvmfspawn(cdexec,kmod,cdwhere,kproc,ktids,kinfo)
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
      kinfo=-1
C
#  else
C
C SHMEM VERSION or PVM onto T3E
C
      INTEGER ji
C
      DO ji=1,kproc
        ktids(ji)=nt3d_tids(ji)
      END DO
      kinfo=kproc
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
      nrs1sync_shmem = SHMEM_SYNC_VALUE
      nrs2sync_shmem = SHMEM_SYNC_VALUE
      nis1sync_shmem = SHMEM_SYNC_VALUE
      nis2sync_shmem = SHMEM_SYNC_VALUE
      nil1sync_shmem = SHMEM_SYNC_VALUE
      nil2sync_shmem = SHMEM_SYNC_VALUE
      ni11sync_shmem = SHMEM_SYNC_VALUE
      ni12sync_shmem = SHMEM_SYNC_VALUE
      ni21sync_shmem = SHMEM_SYNC_VALUE
      ni22sync_shmem = SHMEM_SYNC_VALUE
C
      CALL barrier()
#endif
#endif
C
C
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
CC				  at which ptab is defined for 0 
CC                                initialization
CC				  = 1 ,  T- and W-points
CC				  = 2 ,  U-point
CC				  = 3 ,  V-point
CC				  = 4 ,  F-point
CC                                = 11,  T-point only fold treatment
CC                                = 14,  F-point only fold treatment
CC		ksgn		: control of the sign change
CC				  = 0 , the sign is modified following
CC				  the type of b.c. used
CC				  = 1 , the sign of the field is un-
CC				  changed at the boundaries
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

c#include  "parameter.h"
c#include  "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
C
      INTEGER ktype, ksgn
      INTEGER packsize
      REAL(8) ptab(jpi,jpj,jpk,packsize)
CCC 20/06/2003 Cristiano - Otherwise it does not work
      REAL(8) t3p1(jpi,1,jpk,2)
      REAL(8) t3p2(jpi,1,jpk,2)
      REAL(8) t3p1_my1(jpi,1,jpk,packsize,2)
      REAL(8) t3p2_my1(jpi,1,jpk,packsize,2)

#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

      INTEGER ji,jj,jk,jl
      INTEGER imigr,iihom,ijhom,iloc,ijt,iju
      REAL(8) zsgn
      INTEGER reqs1, reqs2, reqr1, reqr2
      INTEGER jn
! omp variables
      INTEGER :: mytid, ntids, itid

#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

C
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
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



!$omp   parallel default(none) private(jk,jj,ji,mytid)
!$omp&      shared(packsize,nbondi,nperio,jpk,jpj,jpi,ptab,jpim1,jn,ktype,iihom,nlci,jpreci,ijhom,nlcj,jprecj)
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
!$omp   parallel default(none) private(jk,jj,jl,mytid)
!$omp&      shared(packsize,nbondi,iihom,nlci,nreci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab,jn)
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      imigr=jpreci*jpj*jpk*packsize
C
      IF(nbondi.eq.-1) THEN
      CALL shmem_put(t3we_my1(1,1,1,1,2),t3we_my1(1,1,1,1,1),imigr,noea)
      ELSE IF(nbondi.eq.0) THEN
          CALL shmem_put(t3ew_my1(1,1,1,1,2),t3ew_my1(1,1,1,1,1),imigr,nowe)
          CALL shmem_put(t3we_my1(1,1,1,1,2),t3we_my1(1,1,1,1,1),imigr,noea)
      ELSE IF(nbondi.eq.1) THEN
          CALL shmem_put(t3ew_my1(1,1,1,1,2),t3ew_my(1,1,1,1,1),imigr,nowe)
      ENDIF
      CALL  barrier()
      CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
      imigr=jpreci*jpj*jpk*jpbyt*packsize
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
#endif
C
C 2.3 Write Dirichlet lateral conditions
C
C      trcadvparttime = MPI_WTIME()
!$omp   parallel default(none) private(jk,jj,jl,mytid)
!$omp&      shared(packsize,nbondi,iihom,nlci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab,jn)
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
!$omp   parallel default(none) private(jk,ji,jl,mytid)
!$omp&      shared(packsize,nbondj,ijhom,nlcj,nrecj,jprecj,jpk,jpi,t3sn_my1,t3ns_my1,ptab,jn)
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      imigr=jprecj*jpi*jpk*packsize
C
      IF(nbondj.eq.-1) THEN
          CALL shmem_put(t3sn_my1(1,1,1,1,2),t3sn_my1(1,1,1,1,1),imigr,nono)
      ELSE IF(nbondj.eq.0) THEN
          CALL shmem_put(t3ns_my1(1,1,1,1,2),t3ns_my1(1,1,1,1,1),imigr,noso)
          CALL shmem_put(t3sn_my1(1,1,1,1,2),t3sn_my1(1,1,1,1,1),imigr,nono)
      ELSE IF(nbondj.eq.1) THEN
          CALL shmem_put(t3ns_my1(1,1,1,1,2),t3ns_my1(1,1,1,1,1),imigr,noso)
      ENDIF
      CALL  barrier()
      CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
      imigr=jprecj*jpi*jpk*jpbyt*packsize
C
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
#endif
C
C 3.3 Write Dirichlet lateral conditions
C
C        trcadvparttime = MPI_WTIME()

!$omp   parallel default(none) private(jk,ji,jl,mytid)
!$omp&      shared(packsize,nbondj,ijhom,nlcj,jprecj,jpk,jpi,t3sn_my1,t3ns_my1,ptab,jn)
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

!$omp   parallel default(none) private(jk,ji,ijt,iju,mytid)
!$omp&      shared(packsize,npolj,iloc,jpiglo,nimpp,ktype,jpk,nlci,zsgn,ptab,nlcj,jn)
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

!$omp   parallel default(none) private(jk,ji,ijt,iju,mytid)
!$omp&      shared(packsize,npolj,iloc,jpiglo,nimpp,nimppt,nono,ktype,jpk,jpi,t3p1_my1,t3p2_my1,zsgn,ptab,nlcj,jn)
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
          imigr=jprecj*jpi*jpk*packsize
C
          CALL shmem_put(t3p1_my1(1,1,1,1,2),t3p1_my1(1,1,1,1,1),imigr,nono)
          CALL shmem_put(t3p2_my1(1,1,1,1,2),t3p2_my1(1,1,1,1,1),imigr,nono)
          CALL barrier()
          CALL shmem_udcflush()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
          imigr=jprecj*jpi*jpk*jpbyt*packsize
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
#endif      
          ENDIF
C
C 4.3 Write north fold conditions
C

C        trcadvparttime = MPI_WTIME()
C

!$omp   parallel default(none) private(jk,ji,mytid)
!$omp&      shared(packsize,npolj,ktype,jpk,nlci,ptab,t3p1_my1,nimpp,nlcj,t3p2_my1,jpjglo,jn)
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

!$omp   parallel default(none) private(jl,jk,jj,mytid)
!$omp&      shared(packsize,npolj,nbondi,iihom,nlci,nreci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab,jn)
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
          imigr=jpreci*jpj*jpk*packsize
C
          IF(nbondi.eq.-1) THEN
              CALL shmem_put(t3we_my1(1,1,1,1,2),t3we_my1(1,1,1,1,1),imigr,noea)
          ELSE IF(nbondi.eq.0) THEN
              CALL shmem_put(t3ew_my1(1,1,1,1,2),t3ew_my1(1,1,1,1,1),imigr,nowe)
              CALL shmem_put(t3we_my1(1,1,1,1,2),t3we_my1(1,1,1,1,1),imigr,noea)
          ELSE IF(nbondi.eq.1) THEN
              CALL shmem_put(t3ew_my1(1,1,1,1,2),t3ew_my1(1,1,1,1,1),imigr,nowe)
          ENDIF
          CALL  barrier()
          CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
          imigr=jpreci*jpj*jpk*jpbyt*packsize
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
#endif
         ENDIF
C
C 5.3 Write Dirichlet lateral conditions
C
C        trcadvparttime = MPI_WTIME()

!$omp   parallel default(none) private(jl,jk,jj,mytid)
!$omp&      shared(packsize,npolj,nbondi,iihom,nlci,jpreci,jpk,jpj,t3ew_my1,t3we_my1,ptab,jn)
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
CC				  at which ptab is defined for 0 
CC                                initialization
CC				  = 1 ,  T- and W-points
CC				  = 2 ,  U-point
CC				  = 3 ,  V-point
CC				  = 4 ,  F-point
CC                                = 11,  T-point only fold treatment
CC                                = 14,  F-point only fold treatment
CC		ksgn		: control of the sign change
CC				  = 0 , the sign is modified following
CC				  the type of b.c. used
CC				  = 1 , the sign of the field is un-
CC				  changed at the boundaries
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

c#include  "parameter.h"
c#include  "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
C
      INTEGER ktype, ksgn
      REAL(8) ptab(jpi,jpj,jpk)
CCC 20/06/2003 Cristiano - Otherwise it does not work
      REAL(8) t3p1(jpi,1,jpk,2)
      REAL(8) t3p2(jpi,1,jpk,2)

#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

      INTEGER ji,jj,jk,jl
      INTEGER imigr,iihom,ijhom,iloc,ijt,iju
      REAL(8) zsgn
      INTEGER reqs1, reqs2, reqr1, reqr2
      INTEGER jn
C
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      imigr=jpreci*jpj*jpk
C
      IF(nbondi.eq.-1) THEN
          CALL shmem_put(t3we(1,1,1,2),t3we(1,1,1,1),imigr,noea)
      ELSE IF(nbondi.eq.0) THEN
          CALL shmem_put(t3ew(1,1,1,2),t3ew(1,1,1,1),imigr,nowe)
          CALL shmem_put(t3we(1,1,1,2),t3we(1,1,1,1),imigr,noea)
      ELSE IF(nbondi.eq.1) THEN
          CALL shmem_put(t3ew(1,1,1,2),t3ew(1,1,1,1),imigr,nowe)
      ENDIF
      CALL  barrier()
      CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
      imigr=jpreci*jpj*jpk*jpbyt
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
#endif
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      imigr=jprecj*jpi*jpk
C
      IF(nbondj.eq.-1) THEN
          CALL shmem_put(t3sn(1,1,1,2),t3sn(1,1,1,1),imigr,nono)
      ELSE IF(nbondj.eq.0) THEN
          CALL shmem_put(t3ns(1,1,1,2),t3ns(1,1,1,1),imigr,noso)
          CALL shmem_put(t3sn(1,1,1,2),t3sn(1,1,1,1),imigr,nono)
      ELSE IF(nbondj.eq.1) THEN
          CALL shmem_put(t3ns(1,1,1,2),t3ns(1,1,1,1),imigr,noso)
      ENDIF
      CALL  barrier()
      CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
      imigr=jprecj*jpi*jpk*jpbyt
C
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
#endif
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
          imigr=jprecj*jpi*jpk
C
          CALL shmem_put(t3p1(1,1,1,2),t3p1(1,1,1,1),imigr,nono)
          CALL shmem_put(t3p2(1,1,1,2),t3p2(1,1,1,1),imigr,nono)
          CALL barrier()
          CALL shmem_udcflush()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
          imigr=jprecj*jpi*jpk*jpbyt
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
#endif      
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
          imigr=jpreci*jpj*jpk
C
          IF(nbondi.eq.-1) THEN
              CALL shmem_put(t3we(1,1,1,2),t3we(1,1,1,1),imigr,noea)
          ELSE IF(nbondi.eq.0) THEN
              CALL shmem_put(t3ew(1,1,1,2),t3ew(1,1,1,1),imigr,nowe)
              CALL shmem_put(t3we(1,1,1,2),t3we(1,1,1,1),imigr,noea)
          ELSE IF(nbondi.eq.1) THEN
              CALL shmem_put(t3ew(1,1,1,2),t3ew(1,1,1,1),imigr,nowe)
          ENDIF
          CALL  barrier()
          CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
          imigr=jpreci*jpj*jpk*jpbyt
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
#endif
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
CC				  at which ptab is defined for 0 
CC                                initialization
CC                                initialization
CC				  = 1 ,  T- and W-points
CC				  = 2 ,  U-point
CC				  = 3 ,  V-point
CC				  = 4 ,  F-point
CC                                = 11,  T-point only fold treatment
CC                                = 14,  F-point only fold treatment
CC		ksgn		: control of the sign change
CC				  = 0 , the sign is modified following
CC				  the type of b.c. used
CC				  = 1 , the sign of the field is un-
CC				  changed at the boundaries
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

c#include  "parameter.h"
c#include  "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
C
      INTEGER ktype, ksgn
      REAL(8) ptab(jpi,jpj)
CCC 20/06/2003 Cristiano - Otherwise it does not work
      REAL(8) t2p1(jpi,1,2)
      REAL(8) t2p2(jpi,1,2)
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      imigr=jpreci*jpj
C
      IF(nbondi.eq.-1) THEN
          CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)
      ELSE IF(nbondi.eq.0) THEN
          CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)
          CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)
      ELSE IF(nbondi.eq.1) THEN
          CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)
      ENDIF
      CALL  barrier()
      CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
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
# else
C
C PVM VERSION
C
      imigr=jpreci*jpj*jpbyt
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
#endif
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
C
      imigr=jprecj*jpi
C
      IF(nbondj.eq.-1) THEN
          CALL shmem_put(t2sn(1,1,2),t2sn(1,1,1),imigr,nono)
      ELSE IF(nbondj.eq.0) THEN
          CALL shmem_put(t2ns(1,1,2),t2ns(1,1,1),imigr,noso)
          CALL shmem_put(t2sn(1,1,2),t2sn(1,1,1),imigr,nono)
      ELSE IF(nbondj.eq.1) THEN
          CALL shmem_put(t2ns(1,1,2),t2ns(1,1,1),imigr,noso)
      ENDIF
      CALL  barrier()
      CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
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
#else
C
C PVM VERSION
C
      imigr=jprecj*jpi*jpbyt
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
#endif
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
          imigr=jprecj*jpi
C
          CALL shmem_put(t2p1(1,1,2),t2p1(1,1,1),imigr,nono)
          CALL shmem_put(t2p2(1,1,2),t2p2(1,1,1),imigr,nono)
          CALL barrier()
          CALL shmem_udcflush()
C
#  elif defined key_mpp_mpi
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
#  else
C
C PVM VERSION
C
          imigr=jprecj*jpi*jpbyt
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
#endif      
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
          imigr=jpreci*jpj
C
          IF(nbondi.eq.-1) THEN
              CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)
          ELSE IF(nbondi.eq.0) THEN
              CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)
              CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)
          ELSE IF(nbondi.eq.1) THEN
              CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)
          ENDIF
          CALL  barrier()
          CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
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
# else
C
C PVM VERSION
C
          imigr=jpreci*jpj*jpbyt
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
#endif
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

c#include  "parameter.h"
c#include  "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
C
      REAL(8) ptab(jpi,jpj)
CCC 20/06/2003 Cristiano - Otherwise it does not work
      REAL(8) t2p1(jpi,1,2)
      REAL(8) t2p2(jpi,1,2)

#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

      INTEGER ji,jj,jl
      INTEGER imigr,iihom,ijhom
      INTEGER reqs1, reqs2, reqr1, reqr2
C
CCC---------------------------------------------------------------------
CCC  OPA8, LODYC (15/11/96)
CCC---------------------------------------------------------------------
C
C
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
          imigr=jprecj*jpi
C
          CALL shmem_put(t2p1(1,1,2),t2p1(1,1,1),imigr,nono)
          CALL barrier()
          CALL shmem_udcflush()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
          imigr=jprecj*jpi
C
          CALL mppsend(3,t2p1(1,1,1),imigr,nono,0,reqs1)
          CALL mpprecv(3,t2p1(1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
C
#  else
C
C PVM VERSION
C
          imigr=jprecj*jpi*jpbyt
C
          CALL mppsend(3,t2p1(1,1,1),imigr,nono,0,reqs1)
          CALL mpprecv(3,t2p1(1,1,2),imigr,reqr1)
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
C
#endif      
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      imigr=jpreci*jpj
C
      IF(nbondi.eq.-1) THEN
          CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)
        ELSE IF(nbondi.eq.0) THEN
          CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)
          CALL shmem_put(t2we(1,1,2),t2we(1,1,1),imigr,noea)
        ELSE IF(nbondi.eq.1) THEN
          CALL shmem_put(t2ew(1,1,2),t2ew(1,1,1),imigr,nowe)
      ENDIF
      CALL  barrier()
      CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
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
#  else
C
C PVM VERSION
C
      imigr=jpreci*jpj*jpbyt
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
#endif
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
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      imigr=jprecj*jpi
C
      IF(nbondj.eq.-1) THEN
          CALL shmem_put(t2sn(1,1,2),t2sn(1,1,1),imigr,nono)
        ELSE IF(nbondj.eq.0) THEN
          CALL shmem_put(t2ns(1,1,2),t2ns(1,1,1),imigr,noso)
          CALL shmem_put(t2sn(1,1,2),t2sn(1,1,1),imigr,nono)
        ELSE IF(nbondj.eq.1) THEN
          CALL shmem_put(t2ns(1,1,2),t2ns(1,1,1),imigr,noso)
      ENDIF
      CALL  barrier()
      CALL  shmem_udcflush()
C
#  elif defined key_mpp_mpi
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
#  else
C
C PVM VERSION
C
      imigr=jprecj*jpi*jpbyt
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
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
      REAL(8) pmess(*)
      INTEGER kbytes,kdest,ktyp,kid, ireqsend
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
C routine not used
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
      INTEGER iflag
      INTEGER itid_dest,info
C
C      write(*,*)  "#### kbytes = ", kbytes
      CALL mpi_isend(pmess,kbytes,mpi_real8,kdest,ktyp,
     $    mpi_comm_world,ireqsend,iflag)
C
#  else
C
C PVM VERSION
C
      INTEGER iflag
      INTEGER itid_dest,info
C		
      itid_dest = npvm_tids(kdest)
      IF(mppsend_print .NE. 0) THEN
          WRITE(nummpp,*) 'mytid=',npvm_mytid,
     $        ' ========== mppsend.f ========== '
          WRITE(nummpp,*) 'mytid=',npvm_mytid,
     $        ' kbytes=',kbytes,' kdest=',kdest
     $        ,' ktyp=',ktyp,' iflag=',iflag
      END IF
      CALL pvmfinitsend(pvmdataraw, info)
      CALL pvmfpack( byte1, pmess, kbytes, 1, info )
      IF(info .LT. 0) STOP ' mppsend.f : problem in pvmfpack '
      CALL pvmfsend  ( itid_dest , ktyp , iflag )
      IF(iflag .LT. 0) STOP ' mppsend.f : problem in pvmfsend '
      IF(mppsend_print .NE. 0) THEN
          WRITE(nummpp,*) 'mytid=',npvm_mytid,' after:itid_dest='
     $        ,itid_dest,' iflag=',iflag
      END IF
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
      REAL(8) pmess(*)
      INTEGER   kbytes,ktyp, ireqrecv
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
C routine not used
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
      INTEGER istatus(mpi_status_size)
      INTEGER iflag
C
      CALL mpi_irecv(pmess,kbytes,mpi_real8,mpi_any_source,ktyp
     $    ,mpi_comm_world,ireqrecv,iflag)
c
#  else
C
C PVM VERSION
C
      INTEGER itid_kexp,info,ibufid
C
      itid_kexp = -1
      IF(mpprecv_print .NE. 0) THEN
          WRITE(nummpp,*) 'mytid=',npvm_mytid,
     $        ' ============= mpprecv.f ============'
          WRITE(nummpp,*) 'mytid=',npvm_mytid,
     $        ' mpprecv, pvmfrecv, itid_kexp=',itid_kexp
      END IF
      CALL pvmfrecv  ( itid_kexp , ktyp , ibufid )
      IF(mpprecv_print .NE. 0) THEN
          WRITE(nummpp,*) 'mytid=',npvm_mytid,' mpprecv,END pvmfrecv'
     $        ,'ibufid=',ibufid,' npvm_me=',npvm_me
      END IF
      CALL pvmfunpack( byte1, pmess, kbytes, 1, info )
      IF(info .LT. 0) kbytes = info
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"
#endif

      integer req
      INTEGER istatus(mpi_status_size), ierr

#ifdef key_mpp
#ifdef key_mpp_mpi
      call MPI_WAIT(req, istatus, ierr)  
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CCC--------------------------------------------------------------------
CCC  OPA8, IDRIS-LODYC (15/11/96)
CCC--------------------------------------------------------------------
C
      INTEGER kk,kp
      REAL(8) ptab(jpi,jpj),pio(jpi,jpj,jpnij)
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      CALL barrier()
      CALL shmem_put(pio(1,1,npvm_me+1),ptab,jpi*jpj,kp)
      CALL barrier()
C
#elif defined key_mpp_mpi
C
C MPI VERSION
C
      INTEGER itaille,ierror
C
      itaille=jpi*jpj
      CALL mpi_gather(ptab,itaille,mpi_real8,pio,itaille
     $    ,mpi_real8,kp,mpi_comm_world,ierror) 
#else
C
C PVM VERSION
C
      INTEGER imess,ic
      INTEGER ji,jj
      INTEGER ii
      INTEGER reqs1, reqs2, reqr1, reqr2
C
      IF(jpnij .EQ. 1) THEN
          DO jj=1,jpj
            DO ji=1,jpi
              pio(ji,jj,1)=ptab(ji,jj)
            END DO 
          END DO 
          RETURN
      ENDIF
      CALL mppsync
      IF (npvm_me .NE. kp) THEN 
c
c     send data to the root member
c
      imess=kk+ 100000*npvm_me
      CALL mppsend(imess,ptab,jpi*jpj*jpbyt,kp,0,reqs1)
      ELSE 
c
c     receive message form other member
c     of the group
c
          DO ji=0,npvm_nproc-1
            IF  (ji .EQ. npvm_me) THEN 
                pio(:,:,ji+1)=ptab(:,:)
            ELSE 
                imess=kk+ 100000*ji
                CALL mpprecv(imess,pio(1,1,ji+1),jpi*jpj*jpbyt,reqr1)
                CALL mppwait(reqr1)
            ENDIF 
          END DO 
      ENDIF
      CALL mppwait(reqs1)
      CALL  mppsync
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CCC--------------------------------------------------------------------
CCC  OPA8, IDRIS-LODYC (15/11/96)
CCC--------------------------------------------------------------------
C
      INTEGER kk,kp
      REAL(8) ptab(jpi,jpj),pio(jpi,jpj,jpnij)
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      CALL barrier()
      CALL shmem_get(ptab,pio(1,1,npvm_me+1),jpi*jpj,kp)
      CALl barrier()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
      INTEGER itaille, ierror
C
      itaille=jpi*jpj

      CALL mpi_scatter(pio,itaille,mpi_real8,ptab,itaille,
     $     mpi_real8,kp,mpi_comm_world,ierror)
C
#  else
C
C PVM VERSION
C
      INTEGER imess,ic
      INTEGER ji,jj
C
      IF (jpnij .EQ. 1) THEN
          DO jj=1,jpj
            DO ji=1,jpi
              ptab(ji,jj)=pio(ji,jj,1)
            END DO 
          END DO 
          RETURN
      ENDIF
      CALL  mppsync
      imess=kk
      CALL pvmfscatter(ptab,pio,jpi*jpj,jpvmreal,imess,opaall,kp,ic)
      IF(ic .NE. 0) THEN
          WRITE(nummpp,*) "problem pvmfscatter"
     $        ," kk=",kk
     $        ," kp=",kp
      END IF
      CALL  mppsync
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
      INTEGER kdim
      INTEGER ktab(kdim),kwork(kdim)
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      INTEGER ji
      INTEGER, SAVE :: IBOOL=0
C
      IF( kdim .GT. jpmppsum ) THEN
          WRITE(numout,*) 'mppiisl routine : kdim is too big'
          WRITE(numout,*) 'change jpmppsum dimension in mpp.h'
          STOP 'mppiisl'
      ENDIF
C
      DO ji=1,kdim
        niitab_shmem(ji)=ktab(ji)
      END DO
      CALL  barrier()
      IF (IBOOL == 0) THEN 
          CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,kdim,0
     $        ,0,N$PES,ni11wrk_shmem,ni11sync_shmem)
          CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,kdim,0
     $        ,0,N$PES,ni12wrk_shmem,ni12sync_shmem)
        ELSE
          CALL shmem_int8_min_to_all (ni11tab_shmem,niitab_shmem,kdim,0
     $        ,0,N$PES,ni21wrk_shmem,ni21sync_shmem)
          CALL shmem_int8_max_to_all (ni12tab_shmem,niitab_shmem,kdim,0
     $        ,0,N$PES,ni22wrk_shmem,ni22sync_shmem)
      ENDIF
      CALL  barrier()
      IBOOL=IBOOL+1
      IBOOL=mod(IBOOL,2)
      DO ji=1,kdim
        IF (ni11tab_shmem(ji).ne.0.) THEN
            ktab(ji)=ni11tab_shmem(ji)
          ELSE
            ktab(ji)=ni12tab_shmem(ji)
        ENDIF
      END DO      
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
      INTEGER ityd
      INTEGER info,itype,ibuf,iroot
      EXTERNAL PvmIsl2
C
      itype= 100
      iroot=0
      ityd=npvm_tids(npvm_me)
      IF(jpnij .EQ. 1) RETURN
      IF(mpprisl_print .NE. 0) THEN
          WRITE(nummpp,*) 'mppiisl me=',npvm_me
     $        ,' ityd=',ityd
      END IF
      CALL pvmfreduce(PvmIsl2, ktab, kdim, jpvmint,
     $    itype, opaall, iroot,info)
      IF (iroot .EQ. npvm_me) THEN
          CALL pvmfinitsend(pvmdataraw, ibuf )
          CALL pvmfpack(jpvmint,ktab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfpack problem'
          STOP 'mppiisl.f'
          END IF
          CALL pvmfbcast(opaall,itype+1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfbcast problem'
              STOP 'mppiisl.f'
          END IF
        ELSE
          CALL pvmfrecv(iroot,itype+1,ibuf)
          CALL pvmfunpack(jpvmint,ktab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfunpack problem'
              STOP 'mppiisl.f'
          END IF
      ENDIF
      CALL pvmfbarrier(opaall,npvm_nproc,info)
#endif
#endif
      RETURN
      END
      SUBROUTINE PvmIsl2(kdtatyp,kx,ky,kdim,knfo)
      INTEGER kx(kdim),ky(kdim)
      INTEGER kdim,knfo,kdtatyp,ji
      DO ji=1,kdim
        IF (ky(ji) .NE. 0) kx(ji)=ky(ji)
      END DO
      RETURN
      END

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
CC                               (Never user by opa!)
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

CC----------------------------------------------------------------------
CC local declarations
CC ==================
      INTEGER kdim
      INTEGER ktab(kdim),kwork(kdim)
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      INTEGER ji
      INTEGER, SAVE :: IBOOL=0
C
      IF( kdim .GT. jpmppsum ) THEN
          WRITE(numout,*) 'mppimin routine : kdim is too big'
          WRITE(numout,*) 'change jpmppsum dimension in mpp.h'
          STOP 'mppimin'
      ENDIF
C
      DO ji=1,kdim
        niltab_shmem(ji)=ktab(ji)
      END DO
      CALL  barrier()
      IF (IBOOL == 0) THEN 
          CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem,kdim,0,0
     $        ,N$PES,nil1wrk_shmem,nil1sync_shmem )
        ELSE
          CALL shmem_int8_min_to_all (niltab_shmem,niltab_shmem,kdim,0,0
     $        ,N$PES,nil2wrk_shmem,nil2sync_shmem )
      ENDIF
      CALL  barrier()
      IBOOL=IBOOL+1
      IBOOL=mod(IBOOL,2)
      DO ji=1,kdim
        ktab(ji)=niltab_shmem(ji)
      END DO
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
      INTEGER ierror
C
      CALL mpi_allreduce(ktab,kwork,kdim,mpi_integer
     $    ,mpi_min,mpi_comm_world,ierror)

      ktab=kwork
C
#  else
C
C PVM VERSION
C
      INTEGER ityd
      INTEGER info,itype,ibuf,iroot
      EXTERNAL PvmMin
C
      itype= 100
      iroot=0
      ityd=npvm_tids(npvm_me)
      IF(jpnij .EQ. 1) RETURN
      IF(mppimin_print .NE. 0) THEN
          WRITE(nummpp,*) 'mppimin me=',npvm_me
     $        ,' ityd=',ityd
      END IF
      CALL pvmfreduce(PvmMin,ktab, kdim, jpvmint,
     $    itype, opaall, iroot, info)
      IF (iroot .EQ. npvm_me) THEN
          CALL pvmfinitsend(pvmdataraw, ibuf )
          CALL pvmfpack(jpvmint,ktab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfpack problem'
          STOP 'mppimin.f'
          END IF
          CALL pvmfbcast(opaall,itype+1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfbcast problem'
              STOP 'mppimin.f'
          END IF
        ELSE
          CALL pvmfrecv(iroot,itype+1,ibuf)
          CALL pvmfunpack(jpvmint,ktab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfunpack problem'
              STOP 'mppimin.f'
          END IF
      ENDIF
      CALL pvmfbarrier(opaall,npvm_nproc,info)
#endif
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
CC                               (Never user by opa!)
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
      INTEGER kdim
      INTEGER ktab(kdim),kwork(kdim)
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      INTEGER ji
      INTEGER, SAVE :: IBOOL=0
C
      IF( kdim .GT. jpmppsum ) THEN
          WRITE(numout,*) 'mppisum routine : kdim is too big'
          WRITE(numout,*) 'change jpmppsum dimension in mpp.h'
          STOP 'mppisum'
      ENDIF
C
      DO ji=1,kdim
        nistab_shmem(ji)=ktab(ji)
      END DO
      CALL  barrier()
      IF (IBOOL == 0) THEN 
          CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem,kdim,0,0,
     $        N$PES,nis1wrk_shmem,nis1sync_shmem)
        ELSE
          CALL shmem_int8_sum_to_all(nistab_shmem,nistab_shmem,kdim,0,0,
     $        N$PES,nis2wrk_shmem,nis2sync_shmem)
      ENDIF
      CALL  barrier()
      IBOOL=IBOOL+1
      IBOOL=mod(IBOOL,2)
      DO ji=1,kdim
        ktab(ji)=nistab_shmem(ji)
      END DO
C 
#  elif defined key_mpp_mpi
C
C MPI VERSION
      INTEGER ierror
C
      CALL mpi_allreduce(ktab,kwork,kdim,mpi_integer
     $    ,mpi_sum,mpi_comm_world,ierror)

      ktab=kwork
C
#  else
C
C PVM VERSION
C
      INTEGER ityd
      INTEGER info,itype,ibuf,iroot
      EXTERNAL PvmSum
C
      itype= 100
      iroot=0
      ityd=npvm_tids(npvm_me)
      IF(jpnij .EQ. 1) RETURN
      IF(mpprsum_print.ne.0) THEN
          WRITE(nummpp,*) 'mppisum me=',npvm_me
     $        ,' ityd=',ityd
      END IF
      CALL pvmfreduce(PvmSum, ktab, kdim, jpvmint,
     $    itype, opaall, iroot, info)
      IF (iroot .EQ. npvm_me) THEN
          CALL pvmfinitsend(pvmdataraw, ibuf )
          CALL pvmfpack(jpvmint,ktab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfpack problem'
          STOP 'mppisum.f'
          END IF
          CALL pvmfbcast(opaall,itype+1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfbcast problem'
              STOP 'mppisum.f'
          END IF
        ELSE
          CALL pvmfrecv(iroot,itype+1,ibuf)
          CALL pvmfunpack(jpvmint,ktab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfunpack problem'
              STOP 'mppisum.f'
          END IF
      END IF
      CALL pvmfbarrier(opaall,npvm_nproc,info)
C
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
      INTEGER kdim
      REAL(8) ptab(kdim),pwork(kdim)
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      INTEGER ji
      INTEGER, SAVE :: IBOOL=0
C
      IF( kdim .GT. jpmppsum ) THEN
          WRITE(numout,*) 'mpprisl routine : kdim is too big'
          WRITE(numout,*) 'change jpmppsum dimension in mpp.h'
          STOP 'mpprisl'
      ENDIF
C
      DO ji=1,kdim
        wiltab_shmem(ji)=ptab(ji)
      END DO
      CALL  barrier()
      IF (IBOOL == 0) THEN 
          CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem,kdim,0
     $        ,0,N$PES,wi11wrk_shmem,ni11sync_shmem)
          CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem,kdim,0
     $        ,0,N$PES,wi12wrk_shmem,ni12sync_shmem)
        ELSE
          CALL shmem_real8_min_to_all (wi1tab_shmem,wiltab_shmem,kdim,0
     $        ,0,N$PES,wi21wrk_shmem,ni21sync_shmem)
          CALL shmem_real8_max_to_all (wi2tab_shmem,wiltab_shmem,kdim,0
     $        ,0,N$PES,wi22wrk_shmem,ni22sync_shmem)
      ENDIF
      CALL  barrier()
      IBOOL=IBOOL+1
      IBOOL=mod(IBOOL,2)
      DO ji=1,kdim
        IF (wi1tab_shmem(ji).ne.0.) THEN
            ptab(ji)=wi1tab_shmem(ji)
          ELSE
            ptab(ji)=wi2tab_shmem(ji)
        ENDIF
      END DO      
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
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
#  else
C
C PVM VERSION
C
      INTEGER ityd
      INTEGER info,itype,ibuf,iroot
      EXTERNAL PvmIsl
C
      itype= 100
      iroot=0
      ityd=npvm_tids(npvm_me)
      IF(jpnij .EQ. 1) RETURN
      IF(mpprisl_print .NE. 0) THEN
          WRITE(nummpp,*) 'mpprisl me=',npvm_me
     $        ,' ityd=',ityd
      END IF
      CALL pvmfreduce(PvmIsl, ptab, kdim, jpvmreal,
     $    itype, opaall, iroot,info)
      IF (iroot .EQ. npvm_me) THEN
          CALL pvmfinitsend(pvmdataraw, ibuf )
          CALL pvmfpack(jpvmreal,ptab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfpack problem'
          STOP 'mpprisl.f'
          END IF
          CALL pvmfbcast(opaall,itype+1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfbcast problem'
              STOP 'mpprisl.f'
          END IF
        ELSE
          CALL pvmfrecv(iroot,itype+1,ibuf)
          CALL pvmfunpack(jpvmreal,ptab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfunpack problem'
              STOP 'mpprisl.f'
          END IF
      ENDIF
      CALL pvmfbarrier(opaall,npvm_nproc,info)
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
      INTEGER kdim
      REAL(8) ptab(kdim),pwork(kdim)
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      INTEGER ji
      INTEGER, SAVE :: IBOOL=0
C
      IF( kdim .GT. jpmppsum ) THEN
          WRITE(numout,*) 'mpprmax routine : kdim is too big'
          WRITE(numout,*) 'change jpmppsum dimension in mpp.h'
          STOP 'mpprmax'
      ENDIF
C
      DO ji=1,kdim
        wintab_shmem(ji)=ptab(ji)
      END DO
      CALL  barrier()
      IF (IBOOL == 0) THEN 
          CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem,kdim,0
     $        ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
        ELSE
          CALL shmem_real8_max_to_all (wintab_shmem,wintab_shmem,kdim,0
     $        ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
      ENDIF
      CALL  barrier()
      IBOOL=IBOOL+1
      IBOOL=mod(IBOOL,2)
      DO ji=1,kdim
        ptab(ji)=wintab_shmem(ji)
      END DO      
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
      INTEGER ierror
C
      CALL mpi_allreduce(ptab,pwork,kdim,mpi_real8
     $    ,mpi_max,mpi_comm_world,ierror)
      ptab=pwork
C
#  else
C
C PVM VERSION
C
      INTEGER ityd
      INTEGER info,itype,ibuf,iroot
      EXTERNAL PvmMax
C
      itype= 100
      iroot=0
      ityd=npvm_tids(npvm_me)
      IF(jpnij .EQ. 1) RETURN
      IF(mpprmax_print .NE. 0) THEN
          WRITE(nummpp,*) 'mpprmax me=',npvm_me
     $        ,' ityd=',ityd
      END IF
      CALL pvmfreduce(PvmMax, ptab, kdim, jpvmreal,
     $    itype, opaall, iroot, info)
      IF (iroot .EQ. npvm_me) THEN
          CALL pvmfinitsend(pvmdataraw, ibuf )
          CALL pvmfpack(jpvmreal,ptab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfpack problem'
          STOP 'mpprmax.f'
          END IF
          CALL pvmfbcast(opaall,itype+1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfbcast problem'
              STOP 'mpprmax.f'
          END IF
        ELSE
          CALL pvmfrecv(iroot,itype+1,ibuf)
          CALL pvmfunpack(jpvmreal,ptab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfunpack problem'
              STOP 'mpprmax.f'
          END IF
      ENDIF
      CALL pvmfbarrier(opaall,npvm_nproc,info)
C
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
      INTEGER kdim
      REAL(8) ptab(kdim),pwork(kdim)
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      INTEGER ji
      INTEGER, SAVE :: IBOOL=0
C
      IF( kdim .GT. jpmppsum ) THEN
          WRITE(numout,*) 'mpprmin routine : kdim is too big'
          WRITE(numout,*) 'change jpmppsum dimension in mpp.h'
          STOP 'mpprmin'
      ENDIF
C
      DO ji=1,kdim
        wintab_shmem(ji)=ptab(ji)
      END DO
      CALL  barrier()
      IF (IBOOL == 0) THEN 
          CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem,kdim,0
     $        ,0,N$PES,wi1wrk_shmem,ni1sync_shmem)
        ELSE
          CALL shmem_real8_min_to_all (wintab_shmem,wintab_shmem,kdim,0
     $        ,0,N$PES,wi2wrk_shmem,ni2sync_shmem)
      ENDIF
      CALL  barrier()
      IBOOL=IBOOL+1
      IBOOL=mod(IBOOL,2)
      DO ji=1,kdim
        ptab(ji)=wintab_shmem(ji)
      END DO      
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
      INTEGER ierror
C
      CALL mpi_allreduce(ptab,pwork,kdim,mpi_real8
     $    ,mpi_min,mpi_comm_world,ierror)
      ptab=pwork
C
#  else
C
C PVM VERSION
C
      INTEGER ityd
      INTEGER info,itype,ibuf,iroot
      EXTERNAL PvmMin
C
      itype= 100
      iroot=0
      ityd=npvm_tids(npvm_me)
      IF(jpnij .EQ. 1) RETURN
      IF(mpprmin_print .NE. 0) THEN
          WRITE(nummpp,*) 'mpprmin me=',npvm_me
     $        ,' ityd=',ityd
      END IF
      CALL pvmfreduce(PvmMin, ptab, kdim, jpvmreal,
     $    itype, opaall, iroot, info)
      IF (iroot .EQ. npvm_me) THEN
          CALL pvmfinitsend(pvmdataraw, ibuf )
          CALL pvmfpack(jpvmreal,ptab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfpack problem'
          STOP 'mpprmin.f'
          END IF
          CALL pvmfbcast(opaall,itype+1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfbcast problem'
              STOP 'mpprmin.f'
          END IF
        ELSE
          CALL pvmfrecv(iroot,itype+1,ibuf)
          CALL pvmfunpack(jpvmreal,ptab,kdim,1,info)
          IF (info.ne.0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfunpack problem'
              STOP 'mpprmin.f'
          END IF
      ENDIF
      CALL pvmfbarrier(opaall,npvm_nproc,info)
C
#endif
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
CC                               (Never user by opa!)
CC
CC   Modifications:
CC   --------------
CC       original  : 93-09 (M. Imbard)
CC       additions : 96-05 (j. Escobar)
CC       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
CC                          SHMEM and MPI versions
C-----------------------------------------------------------------------

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
      INTEGER kdim
      REAL(8) ptab(kdim),pwork(kdim)
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
      INTEGER ji
      INTEGER, SAVE :: IBOOL=0
C
      IF( kdim .GT. jpmppsum ) THEN
          WRITE(numout,*) 'mpprsum routine : kdim is too big'
          WRITE(numout,*) 'change jpmppsum dimension in mpp.h'
          STOP 'mpprsum'
      ENDIF
C
      DO ji=1,kdim
        wrstab_shmem(ji)=ptab(ji)
      END DO
      CALL  barrier()
      IF (IBOOL == 0) THEN 
          CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem,kdim,0
     $        ,0,N$PES,wrs1wrk_shmem,nrs1sync_shmem )
        ELSE
          CALL shmem_real8_sum_to_all (wrstab_shmem,wrstab_shmem,kdim,0
     $        ,0,N$PES,wrs2wrk_shmem,nrs2sync_shmem )
      ENDIF
      CALL  barrier()
      IBOOL=IBOOL+1
      IBOOL=mod(IBOOL,2)
      DO ji=1,kdim
        ptab(ji)=wrstab_shmem(ji)
      END DO
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
      INTEGER ierror
C
      CALL mpi_allreduce(ptab,pwork,kdim,mpi_real8
     $    ,mpi_sum,mpi_comm_world,ierror)
      ptab=pwork
C
#  else
C
C PVM VERSION
C
      INTEGER ityd
      INTEGER info,itype,ibuf,iroot
      EXTERNAL PvmSum
C
      itype= 100
      iroot=0
      ityd=npvm_tids(npvm_me)
      IF(jpnij .EQ. 1) RETURN
      IF(mpprsum_print .NE. 0) THEN
          WRITE(nummpp,*) 'mpprsum me=',npvm_me
     $        ,' ityd=',ityd
      END IF
      CALL pvmfreduce(PvmSum, ptab, kdim, jpvmreal,
     $    itype, opaall, iroot, info)
      IF (iroot .EQ. npvm_me) THEN
          CALL pvmfinitsend(pvmdataraw, ibuf )
          CALL pvmfpack(jpvmreal,ptab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfpack problem'
          STOP 'mpprsum.f'
          END IF
          CALL pvmfbcast(opaall,itype+1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfbcast problem'
              STOP 'mpprsum.f'
          END IF
        ELSE
          CALL pvmfrecv(iroot,itype+1,ibuf)
          CALL pvmfunpack(jpvmreal,ptab,kdim,1,info)
          IF (info .NE. 0) THEN
              WRITE(nummpp,*) 'me=',npvm_me,' pvmfunpack problem'
              STOP 'mpprsum.f'
          END IF
      ENDIF
      CALL pvmfbarrier(opaall,npvm_nproc,info)
C
#endif
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

c#include "parameter.h"
c#include "common.h"
      USE myalloc
      USE myalloc_mpp

CC+CC  Implicit typing is never allowed
        IMPLICIT NONE
CC+CC  Implicit typing is never allowed

C-----------------------------------------------------------------------
C
#ifdef key_mpp

c#include "parameter.mpp.h"
c#include "common.mpp.h"

C
#ifdef key_mpp_shmem
C
C SHMEM VERSION
C
         call barrier()
C
#  elif defined key_mpp_mpi
C
C MPI VERSION
C
      INTEGER ierror
C
      CALL mpi_barrier(mpi_comm_world,ierror)
C
#  else
C
C PVM VERSION
C
      INTEGER info
C
      IF(jpnij .EQ. 1) RETURN
      IF(mppsync_print .NE. 0) THEN
          WRITE(nummpp,*) 'mppsync me=',npvm_me
      END IF
      CALL pvmfbarrier(opaall,npvm_nproc,info)
      IF(info .NE. 0) THEN
          WRITE(nummpp,*) 'me=',npvm_me,' barrier problem'
          STOP 'mppsync.f'
      END IF
#endif
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
#ifdef key_mpp

      CALL mppsync
#endif
C
C 1. Unit close
C -------------
C

      CLOSE( numnam )


C
      IF(lrstar) CLOSE( numrst )
C
#if defined key_islands
      CLOSE( numisl )
      IF(lwp) CLOSE( numisp )
#endif
C
      IF(lwp) THEN
CMAF          IF( ntopo.EQ.1 ) CLOSE( numbat )

          IF( nmsh .EQ.1 ) CLOSE( nummsh )
      ENDIF
C
#if defined key_temdta
      CLOSE( numte1 )
#endif
C
#if defined key_saldta
      CLOSE( numsa1 )
#endif
C
#if defined key_icedta
      IF(lwp) CLOSE( numice )
#endif
C
#if defined key_tau && ! defined key_coupled
      CLOSE( numtau )
#endif
C
#if defined key_flx && ! defined key_coupled
      CLOSE( numfl1 )
#endif
CMAF      CLOSE( numicl )



CMAF      IF(lwp) CLOSE( numsol )
C
#if defined key_diatrdtra || defined key_diatrddyn
      CLOSE( numtrd )
#endif
C
#if defined key_coupled
      CLOSE( numlhf )
      CLOSE( numlts )
#endif
C
#if defined key_sst
      CLOSE( numsst )
#endif
C
C 2. Mpp synchroneus
C ------------------
C
#ifdef key_mpp

      CALL mppsync
#    ifdef key_mpp_mpi
      CALL mpi_finalize(info)
#      else
      CALL pvmfexit( info )
#    endif
#endif
C
      RETURN
      END
