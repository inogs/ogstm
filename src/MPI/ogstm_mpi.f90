
MODULE ogstm_mpi_module

#ifdef ExecDA
#include <petscversion.h>
#if PETSC_VERSION_GE(3,8,0)
#include "petsc/finclude/petscvec.h"
#else
#include <petsc/finclude/petscvecdef.h>
#endif
#endif

USE myalloc
USE mpi

#ifdef ExecDA
use DA_mem
use mpi_str, only: Var3DCommunicator
use petscvec, only: PETSC_COMM_WORLD, PETSC_NULL_CHARACTER
#endif

implicit NONE



contains
!! mpp routines
!!
!! mynode
!! mpplnk_my
!! mpprecv
!! mppsend
!! mppwait
!! mppsync
!! mppstop
!!
!!!---------------------------------------------------------------------
!!!
!!!                       routine mynode
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Massively parallel processors
!!!     Find processor unit
!!!
!!   Input :
!!   -----
!!      argument                :
!!
!!   Modifications:
!!   --------------
!!       original  : 93-09 (M. Imbard)
!!       additions : 96-05 (j. Escobar)
!!       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
!!                          SHMEM and MPI versions
!!----------------------------------------------------------------------

SUBROUTINE mynode
      INTEGER :: ierr

!#ifdef ExecDA
!      PetscErrorCode :: stat
!#endif

#ifdef key_mpp_mpi

      CALL mpi_comm_rank(mpi_comm_world,myrank,ierr)
      CALL mpi_comm_size(mpi_comm_world,mpi_glcomm_size,ierr)
      !call parlec ! in order to read DA_Nprocs

!#ifdef ExecDA
!      if(V3D_VAR_PARALLEL) then
!        call MPI_Comm_split(MPI_COMM_WORLD, nodes * TREd_procs_per_node, myrank, Var3DCommunicator, ierr)
        
!        PETSC_COMM_WORLD = Var3DCommunicator
!        call PetscInitialize(PETSC_NULL_CHARACTER,stat)
!        CHKERRQ(stat)
!      else
!        call MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, myrank, Var3DCommunicator, ierr)
!      endif
!#endif

#else
      mpi_glcomm_size = 1
      myrank = 0
#endif

END SUBROUTINE

     
!!!---------------------------------------------------------------------
!!!
!!!                       routine mpplnk_my
!!!                     ******************
!!!
!!!  Purpose :
!!!  ---------
!!!      Message passing manadgement
!!!
!!   Method :
!!   -------
!!       Use mppsend and mpprecv function for passing mask between
!!       processors following neighboring subdomains.
!!
!!   Input :
!!   -----
!!      argument
!!              ptab            : variable array
!!
!!      common
!!            /COMDOM/ : domain parameters
!!                    nbondi : mark for "east-west local boundary"
!!                    nbondj : mark for "north-south local boundary"
!!                    noea   : number for local neighboring processors
!!                    nowe   : number for local neighboring processors
!!                    noso   : number for local neighboring processors
!!                    nono   : number for local neighboring processors
!!
!!----------------------------------------------------------------------
 SUBROUTINE mpplnk_my(ptab)

      double precision ptab(jpk,jpj,jpi)


#ifdef key_mpp_mpi

      INTEGER jk,jj,ji
      INTEGER reqs1, reqs2, reqr1, reqr2
      INTEGER reqs3, reqs4, reqr3, reqr4
      INTEGER jw, packsize

!!     trcadvparttime = MPI_WTIME()

!!
!!2. East and west directions exchange
!!------------------------------------



!!
!!2.2 Migrations
!!
!!
!               3     4
!               |     ^
!               |     |
!               v     |
!           ________________
!          |                |
!     1<-- |                | 1 <--
!     2--> |                | 2 -->
!          |________________|
!               3     4
!               |     ^
!               |     |
!               v     |

    packsize=jpk*jpj

      IF(nbondi.eq.-1) THEN ! We are at the west side of the domain

          CALL mppsend(2,ptab(:,:,jpi-1),packsize,noea,0,reqs1)
          CALL mpprecv(1,ptab(:,:,  jpi),packsize,reqr1)

      ELSE IF(nbondi.eq.0) THEN
          CALL mppsend(1, ptab(:,:    ,2),packsize,nowe,0,reqs1)
          CALL mppsend(2, ptab(:,:,jpi-1),packsize,noea,0,reqs2)

          CALL mpprecv(1,ptab(:,:,jpi),packsize,reqr1)
          CALL mpprecv(2,ptab(:,:,  1),packsize,reqr2)

      ELSE IF(nbondi.eq.1) THEN ! We are at the east side of the domain

          CALL mppsend(1,ptab(:,:,2), packsize, nowe,0, reqs1)
          CALL mpprecv(2,ptab(:,:,1), packsize, reqr1)


      ENDIF


!!
!!
!!3. North and south directions
!!-----------------------------
!!
!!3.1 Read Dirichlet lateral conditions
!!


      IF(nbondj.eq.0.or.nbondj.eq.-1) THEN
         DO jw=1,NORTH_count_send
              ji = NORTHpoints_send(1,jw)
              jk = NORTHpoints_send(2,jw)
              tn_send(jw) = ptab(jk,jpj-1,ji)
         ENDDO
     ENDIF
     IF(nbondj.eq.0.or.nbondj.eq.1) THEN
         DO jw=1,SOUTH_count_send
             ji = SOUTHpoints_send(1,jw)
             jk = SOUTHpoints_send(2,jw)
             ts_send(jw) = ptab(jk,2,ji)
         ENDDO


      ENDIF!    PACK_LOOP4


!!
!!2.2 Migrations
!!
!!

      IF(nbondj.eq.-1) THEN ! We are at the south side of the domain
          CALL mppsend(4,tn_send,NORTH_count_send,nono,0,reqs4)
          CALL mpprecv(3,tn_recv,NORTH_count_recv,reqr3)
          CALL mppwait(reqs4)
          CALL mppwait(reqr3)
      ELSE IF(nbondj.eq.0) THEN
          CALL mppsend(4, tn_send,NORTH_count_send,nono,0,reqs4)
          CALL mppsend(3, ts_send,SOUTH_count_send,noso,0,reqs3)
          CALL mpprecv(3,tn_recv,NORTH_count_recv,reqr3)
          CALL mpprecv(4,ts_recv,SOUTH_count_recv,reqr4)

          CALL mppwait(reqs4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr3)
          CALL mppwait(reqr4)
      ELSE IF(nbondj.eq.1) THEN ! We are at the north side of the domain
          CALL mppsend(3,ts_send, SOUTH_count_send, noso,0, reqs3)
          CALL mpprecv(4,ts_recv, SOUTH_count_recv, reqr4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr4)
      ENDIF



!!
!!2.3 Write Dirichlet lateral conditions
!!

      IF(nbondj.eq.0.or.nbondj.eq.1) THEN ! All but south boundary, we received from south

         DO jw=1,SOUTH_count_recv
              ji = SOUTHpoints_recv(1,jw)
              jk = SOUTHpoints_recv(2,jw)
              ptab(jk,1,ji)= ts_recv(jw)
         ENDDO

      ENDIF

      IF(nbondj.eq.-1.or.nbondj.eq.0) THEN ! All but north boundary, we received from north

        DO jw=1,NORTH_count_recv
              ji = NORTHpoints_recv(1,jw)
              jk = NORTHpoints_recv(2,jw)
             ptab(jk,jpj,ji)= tn_recv(jw)
        ENDDO

      ENDIF ! PACK_LOOP5


!!!  East - West waits

      IF(nbondi.eq.-1) THEN ! We are at the west side of the domain
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondi.eq.0) THEN
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondi.eq.1) THEN ! We are at the east side of the domain
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ENDIF

#endif

END SUBROUTINE

SUBROUTINE mpplnk_my_openacc(ptab)

      double precision ptab(jpk,jpj,jpi)


#ifdef key_mpp_mpi

      INTEGER jk,jj,ji
      INTEGER reqs1, reqs2, reqr1, reqr2
      INTEGER reqs3, reqs4, reqr3, reqr4
      INTEGER jw, packsize

!!     trcadvparttime = MPI_WTIME()

!!
!!2. East and west directions exchange
!!------------------------------------



!!
!!2.2 Migrations
!!
!!
!               3     4
!               |     ^
!               |     |
!               v     |
!           ________________
!          |                |
!     1<-- |                | 1 <--
!     2--> |                | 2 -->
!          |________________|
!               3     4
!               |     ^
!               |     |
!               v     |

    packsize=jpk*jpj


!$acc host_data use_device(ptab)
    
      IF(nbondi.eq.-1) THEN ! We are at the west side of the domain

          CALL mppsend(2,ptab(:,:,jpi-1),packsize,noea,0,reqs1)
          CALL mpprecv(1,ptab(:,:,  jpi),packsize,reqr1)

      ELSE IF(nbondi.eq.0) THEN
          CALL mppsend(1, ptab(:,:    ,2),packsize,nowe,0,reqs1)
          CALL mppsend(2, ptab(:,:,jpi-1),packsize,noea,0,reqs2)

          CALL mpprecv(1,ptab(:,:,jpi),packsize,reqr1)
          CALL mpprecv(2,ptab(:,:,  1),packsize,reqr2)

      ELSE IF(nbondi.eq.1) THEN ! We are at the east side of the domain

          CALL mppsend(1,ptab(:,:,2), packsize, nowe,0, reqs1)
          CALL mpprecv(2,ptab(:,:,1), packsize, reqr1)


      ENDIF
!$acc end host_data 

!!
!!
!!3. North and south directions
!!-----------------------------
!!
!!3.1 Read Dirichlet lateral conditions
!!

!$acc kernels default(present)
      IF(nbondj.eq.0.or.nbondj.eq.-1) THEN
         DO jw=1,NORTH_count_send
              ji = NORTHpoints_send(1,jw)
              jk = NORTHpoints_send(2,jw)
              tn_send(jw) = ptab(jk,jpj-1,ji)
         ENDDO
     ENDIF
     IF(nbondj.eq.0.or.nbondj.eq.1) THEN
         DO jw=1,SOUTH_count_send
             ji = SOUTHpoints_send(1,jw)
             jk = SOUTHpoints_send(2,jw)
             ts_send(jw) = ptab(jk,2,ji)
         ENDDO


      ENDIF!    PACK_LOOP4
!$acc end kernels

!!
!!2.2 Migrations
!!
!!

      
      !$acc host_data use_device(tn_send,tn_recv,ts_send,ts_recv)
      
      IF(nbondj.eq.-1) THEN ! We are at the south side of the domain
          CALL mppsend(4,tn_send,NORTH_count_send,nono,0,reqs4)
          CALL mpprecv(3,tn_recv,NORTH_count_recv,reqr3)
          CALL mppwait(reqs4)
          CALL mppwait(reqr3)
      ELSE IF(nbondj.eq.0) THEN
          CALL mppsend(4, tn_send,NORTH_count_send,nono,0,reqs4)
          CALL mppsend(3, ts_send,SOUTH_count_send,noso,0,reqs3)
          CALL mpprecv(3,tn_recv,NORTH_count_recv,reqr3)
          CALL mpprecv(4,ts_recv,SOUTH_count_recv,reqr4)

          CALL mppwait(reqs4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr3)
          CALL mppwait(reqr4)
      ELSE IF(nbondj.eq.1) THEN ! We are at the north side of the domain
          CALL mppsend(3,ts_send, SOUTH_count_send, noso,0, reqs3)
          CALL mpprecv(4,ts_recv, SOUTH_count_recv, reqr4)
          CALL mppwait(reqs3)
          CALL mppwait(reqr4)
      ENDIF
      !$acc end host_data


!!
!!2.3 Write Dirichlet lateral conditions
!!

      !$acc kernels default(present)
      IF(nbondj.eq.0.or.nbondj.eq.1) THEN ! All but south boundary, we received from south

         DO jw=1,SOUTH_count_recv
              ji = SOUTHpoints_recv(1,jw)
              jk = SOUTHpoints_recv(2,jw)
              ptab(jk,1,ji)= ts_recv(jw)
         ENDDO

      ENDIF

      IF(nbondj.eq.-1.or.nbondj.eq.0) THEN ! All but north boundary, we received from north

        DO jw=1,NORTH_count_recv
              ji = NORTHpoints_recv(1,jw)
              jk = NORTHpoints_recv(2,jw)
             ptab(jk,jpj,ji)= tn_recv(jw)
        ENDDO

      ENDIF ! PACK_LOOP5
      !$acc end kernels

!!!  East - West waits

      IF(nbondi.eq.-1) THEN ! We are at the west side of the domain
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ELSE IF(nbondi.eq.0) THEN
          CALL mppwait(reqs1)
          CALL mppwait(reqs2)
          CALL mppwait(reqr1)
          CALL mppwait(reqr2)
      ELSE IF(nbondi.eq.1) THEN ! We are at the east side of the domain
          CALL mppwait(reqs1)
          CALL mppwait(reqr1)
      ENDIF

#endif

END SUBROUTINE


!!!---------------------------------------------------------------------
!!!
!!!                       routine mppsend
!!!                     *******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Send messag passing array
!!
!!   Input :
!!   -----
!!      argument                :
!!                   ktyp   -> Tag of the message
!!                   pmess  -> array of double precision to send
!!                   kbytes -> size of pmess in double precision
!!                   kdest  -> receive process number
!!                   kid    _> ? (note used)
!!
!!   Modifications:
!!   --------------
!!       original  : 93-09 (M. Imbard)
!!       additions : 96-05 (j. Escobar)
!!       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
!!                          SHMEM and MPI versions
!!----------------------------------------------------------------------

SUBROUTINE mppsend(ktyp,pmess,kbytes,kdest,kid,ireqsend)


      double precision pmess(*)
      INTEGER kbytes,kdest,ktyp,kid, ireqsend

#ifdef key_mpp_mpi





      INTEGER iflag
      CALL mpi_isend(pmess,kbytes,mpi_real8,kdest,ktyp, &
     &    mpi_comm_world,ireqsend,iflag)


#endif
      RETURN

END SUBROUTINE

!!!---------------------------------------------------------------------
!!!
!!!                       routine mpprecv
!!!                     *******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Receive messag passing array
!!
!!   Input :
!!   -----
!!      argument
!!                   ktyp    -> Tag of the recevied message
!!                   pmess   -> array of double precision
!!                   kbytes  -> suze of the array pmess


!!
!!   Modifications:
!!   --------------
!!       original  : 93-09 (M. Imbard)
!!       additions : 96-05 (j. Escobar)
!!       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
!!                          SHMEM and MPI versions
!!----------------------------------------------------------------------

 SUBROUTINE mpprecv(ktyp,pmess,kbytes,ireqrecv)

      double precision pmess(*)
      INTEGER   kbytes,ktyp, ireqrecv

#ifdef key_mpp_mpi

      INTEGER iflag

      CALL mpi_irecv(pmess,kbytes,mpi_real8,mpi_any_source,ktyp,mpi_comm_world,ireqrecv,iflag)
#endif

      RETURN

END SUBROUTINE

!!!---------------------------------------------------------------------
!!!
!!!                       routine mppwait
!!!                     *******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Wait message passing isend/irecv
!!
!!   Input :
!!   -----
!!      argument
!!----------------------------------------------------------------------
SUBROUTINE mppwait(req)
      INTEGER istatus(mpi_status_size), ierr
      integer req
      


#ifdef key_mpp_mpi
      call MPI_WAIT(req, istatus, ierr)
#endif
      RETURN
END SUBROUTINE

!!!---------------------------------------------------------------------
!!!
!!!                       routine mppsync
!!!                     *******************
!!!
!!!  Purpose :
!!!  ---------
!!!     Massively parallel processors, synchroneous
!!
!!   Modifications:
!!   --------------
!!       original  : 93-09 (M. Imbard)
!!       additions : 96-05 (j. Escobar)
!!       additions : 98-05 (M. Imbard, J. Escobar, L. Colombet )
!!                          SHMEM and MPI versions
!!----------------------------------------------------------------------
SUBROUTINE mppsync()

!!----------------------------------------------------------------------

#ifdef key_mpp_mpi

      INTEGER ierror

      CALL mpi_barrier(mpi_comm_world,ierror)

#endif
      RETURN
END SUBROUTINE

SUBROUTINE mppstop
!!!---------------------------------------------------------------------
!!!
!!!                       routine mppstop
!!!                     *******************
!!!
!!!  purpose :
!!!  --------
!!!     Stop massilively parallel processors method
!!


      INTEGER info

#ifdef key_mpp_mpi
      CALL mppsync
#endif

      RETURN
END SUBROUTINE

END MODULE ogstm_mpi_module
