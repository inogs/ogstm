      SUBROUTINE trcdmp
!---------------------------------------------------------------------
!
!                       ROUTINE trcdmp
!                     ******************
!
!  Purpose :
!  --------
!    Compute (if asked) the passive tracer trend due to a newtonian
!    damping of the tracer field towards given data field and add it
!    to the general tracer trends.
!
!   Method :
!   --------
!    Default key          : empty routine, no damping trend
!    'key_trc_dmp' defined :
!       Newtonian damping towards tdta and sdta computed and add to
!       the general tracer trends:
!                     trc = ta + restotrc * (trcdta - trcb)
!       The trend is computed either throughout the water column
!       (nlmdmptrc=0) or in area of weak vertical mixing (nlmdmptrc=1) or
!       below the well mixed layer (nlmdmptrc=2)



       USE myalloc
       USE myalloc_mpp
       USE BC_mem
       IMPLICIT NONE

!----------------------------------------------------------------------
! local declarations
! ==================

      INTEGER ji, jj, jk, jn
      INTEGER tra_idx
      INTEGER(4) jv
      INTEGER shift
      REAL(8)  ztra
      INTEGER :: mytid, ntids
#ifdef __OPENMP
      INTEGER ::  omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
      EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

#ifdef __OPENMP
      ntids = omp_get_max_threads() ! take the number of threads
      mytid = -1000000
#else
      ntids = 1
      mytid = 0
#endif


      trcdmpparttime = MPI_WTIME() ! cronometer-start


! Boundary conditions fo Gib area
       IF (Gsize .NE. 0) THEN
         DO jn=1, jn_gib, ntids
!$omp   parallel default(none) private(mytid,ji,jj,jk,ztra,jv,tra_idx)
!$omp&                         shared(jn,tra,trb,gib, Gsize, gib_ridxt,restotr,tra_matrix_gib,jn_gib)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
      IF( mytid + jn <= jn_gib ) THEN

          tra_idx=tra_matrix_gib(jn+mytid)
          DO jv=1, Gsize
             ji = gib_ridxt(2,jv)
             jj = gib_ridxt(3,jv)
             jk = gib_ridxt(4,jv)

             ztra = restotr(ji,jj,jk,tra_idx) * ( gib(jv,jn+mytid)-trb(ji,jj,jk,tra_idx) )

             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + ztra
          ENDDO
      ENDIF
!$omp end parallel

         ENDDO 
 
       ENDIF

! Boundary conditions fo Gib area
       IF (GsizeALK .NE. 0) THEN
         DO jn=5,6
          tra_idx=tra_matrix_gib(jn)
          DO jv=1, GsizeALK
             ji = gib_ridxtALK(2,jv)
             jj = gib_ridxtALK(3,jv)
             jk = gib_ridxtALK(4,jv)
             ztra = restotr(ji,jj,jk,tra_idx) * ( gibALK(jv,jn)-trb(ji,jj,jk,tra_idx) )
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + ztra
          ENDDO

         ENDDO

       ENDIF

! Boundary conditions fo Gib area
       IF (Gsize .NE. 0) THEN
         jn=7
          tra_idx=tra_matrix_gib(jn)
          DO jv=1, Gsize
             ji = gib_ridxt(2,jv)
             jj = gib_ridxt(3,jv)
             jk = gib_ridxt(4,jv)
             ztra = restotr(ji,jj,jk,tra_idx) * ( gib(jv,jn)-trb(ji,jj,jk,tra_idx) )
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + ztra

         ENDDO

       ENDIF

! Boundary conditions for rivers
       IF (Rsize .NE. 0) THEN
         DO jn=1, jn_riv,ntids
!$omp   parallel default(none) private(mytid,ji,jj,jk,jv,tra_idx,shift)
!$omp&                         shared(jn,tra, riv,Rsize, riv_ridxt,tra_matrix_riv,jn_riv,tra_DIA)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
        IF( mytid +jn <= jn_riv ) THEN
          tra_idx=tra_matrix_riv(jn+mytid)

          DO jv=1, Rsize

             ji = riv_ridxt(2,jv)
             jj = riv_ridxt(3,jv)
             jk = riv_ridxt(4,jv)
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + riv(jv,jn+mytid)
          ENDDO


        ENDIF
!$omp end parallel


         ENDDO


       ENDIF


! Boundary conditions for Atmosphere
       IF (Asize .NE. 0) THEN
         DO jn=1, jn_atm, ntids
!$omp   parallel default(none) private(mytid,ji,jj,jk,jv,tra_idx,shift)
!$omp&                         shared(jn,tra,atm,Asize, atm_ridxt,tra_matrix_atm,jn_atm,tra_DIA)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif

      IF( mytid + jn <= jn_atm ) THEN
          tra_idx=tra_matrix_atm(jn+mytid)
          DO jv=1, Asize

             ji = atm_ridxt(2,jv)
             jj = atm_ridxt(3,jv)
             jk = atm_ridxt(4,jv)
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + atm(jv,jn+mytid)
         ENDDO


      ENDIF
!$omp end parallel


         ENDDO

       ENDIF




       trcdmpparttime = MPI_WTIME() - trcdmpparttime ! cronometer-stop
       trcdmptottime  = trcdmptottime + trcdmpparttime



      END SUBROUTINE trcdmp
