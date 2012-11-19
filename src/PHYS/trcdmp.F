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
      INTEGER nmldmp_gib,nmldmp_riv,nmldmp_atm
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


       trcdmpparttime = MPI_WTIME() ! F79 cronometer-start



!      tra_matrix(1) = ppO2o
!      tra_matrix(2) = ppN1p
!      tra_matrix(3) = ppN3n
!      tra_matrix(4) = ppN5s

      nmldmp_gib = 0
      nmldmp_riv = 1
      nmldmp_atm = 1


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
             ji = gib_ridxt(jv,2)
             jj = gib_ridxt(jv,3)
             jk = gib_ridxt(jv,4)
c           ztra = 0 ************* DEFAULT CASE
c           IF ( nmldmp_gib.EQ.0 ) THEN
!             write(*,*) 'tra_idx  -->', tra_idx
!             write(*,*) 'restotr --> ',restotr(ji,jj,jk,tra_idx)
!             write(*,*) 'gib(jv,jn) -->', gib(jv,jn)
!             write(*,*) 'trb(ji,jj,jk,tra_idx) -->', trb(ji,jj,jk,tra_idx)
             ztra = restotr(ji,jj,jk,tra_idx) * ( gib(jv,jn+mytid)-trb(ji,jj,jk,tra_idx) )
c           ENDIF
c           IF ( nmldmp_gib.EQ.1 ) THEN
c             ztra = gib(jv,jn)
c           ENDIF
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + ztra
          ENDDO
      ENDIF
!$omp end parallel

         ENDDO 
!        write(103+rank,*) 'NEW - > restotr', restotr(20,62,1,tra_matrix(1))
!        write(103+rank,*) 'NEW - > gib    ', gib(1,1)
!        write(103+rank,*) 'NEW - > trb    ', trb(20,62,1,tra_matrix(1))
!        write(103+rank,*) 'NEW - > ztra   ' , restotr(20,62,1,tra_matrix(1))* ( gib(1,1)-trb(20,62,1,tra_matrix(1)) )
!        CALL flush_(103+rank)
 
       ENDIF


! Boundary conditions for rivers
       IF (Rsize .NE. 0) THEN
         DO jn=1, jn_riv,ntids
!$omp   parallel default(none) private(mytid,ji,jj,jk,jv,tra_idx,shift)
!$omp&                         shared(jn,tra, riv,Rsize, riv_ridxt,tra_matrix_riv,jn_riv,tra_pp)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif
        IF( mytid +jn <= jn_riv ) THEN
          tra_idx=tra_matrix_riv(jn+mytid)

          DO jv=1, Rsize

             ji = riv_ridxt(jv,2)
             jj = riv_ridxt(jv,3)
             jk = riv_ridxt(jv,4)
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + riv(jv,jn+mytid) ! since nmldmp_riv=1
          ENDDO

           !output diagnostic on jptra_dia
             if (tra_idx .EQ. ppN1p) shift = 2
             if (tra_idx .EQ. ppN3n) shift = 4
             if (tra_idx .EQ. ppN5s) shift = 6
             if (tra_idx .EQ. ppO3c) shift = 8
             if (tra_idx .EQ. ppO3h) shift = 10

             DO jv=1, Rsize

                ji = riv_ridxt(jv,2)
                jj = riv_ridxt(jv,3)
                jk = riv_ridxt(jv,4) + shift

                tra_pp(ji,jj,jk,jptra_dia) = riv(jv,jn+mytid) ! since nmldmp_riv=1

             ENDDO

        ENDIF
!$omp end parallel

c           ztra = 0
c           IF ( nmldmp_riv.EQ.0 ) THEN
c             ztra = restotr(ji,jj,jk,tra_idx)* ( riv(jv,jn)-trb(ji,jj,jk,tra_idx) )
c           ENDIF
c           IF ( nmldmp_riv.EQ.1 ) THEN
c             ztra = riv(jv,jn)
c           ENDIF
c             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + ztra

         ENDDO
!        write(103+rank,*) 'OLD - >', restotr(20,62,1,tra_matrix(1)) * ( trdta(20,62,1,tra_matrix(1))-trb(20,62,1,tra_matrix(1)) )
!        write(103+rank,*) 'NEW - >', restotr(20,62,1,tra_matrix(1))* ( riv(1,1)-trb(20,62,1,tra_matrix(1)) )
!        CALL flush_(103+rank)

       ENDIF


! Boundary conditions for Atmosphere
       IF (Asize .NE. 0) THEN
         DO jn=1, jn_atm, ntids
!$omp   parallel default(none) private(mytid,ji,jj,jk,jv,tra_idx,shift)
!$omp&                         shared(jn,tra,atm,Asize, atm_ridxt,tra_matrix_atm,jn_atm,tra_pp)
#ifdef __OPENMP
        mytid = omp_get_thread_num()  ! take the thread ID
#endif

      IF( mytid + jn <= jn_atm ) THEN
          tra_idx=tra_matrix_atm(jn+mytid)
          DO jv=1, Asize

             ji = atm_ridxt(jv,2)
             jj = atm_ridxt(jv,3)
             jk = atm_ridxt(jv,4)
             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + atm(jv,jn+mytid) ! since nmldmp_atm=1
         ENDDO

                 !output diagnostic on jptra_dia
             if (tra_idx .EQ. ppN1p) shift = 12
             if (tra_idx .EQ. ppN3n) shift = 13


             DO jv=1, Asize

                ji = atm_ridxt(jv,2)
                jj = atm_ridxt(jv,3)
                jk = atm_ridxt(jv,4) + shift

                tra_pp(ji,jj,jk,jptra_dia) = atm(jv,jn+mytid) ! since nmldmp_riv=1

             ENDDO
      ENDIF
!$omp end parallel

!            ztra = 0 ************ DEFAULT CASE
!           IF ( nmldmp_atm.EQ.0 ) THEN
!             ztra = restotr(ji,jj,jk,tra_idx)* ( atm(jv,jn)-trb(ji,jj,jk,tra_idx) )
!           ENDIF
!           IF ( nmldmp_atm.EQ.1 ) THEN
!             ztra = atm(jv,jn)
!           ENDIF
!             tra(ji,jj,jk,tra_idx) = tra(ji,jj,jk,tra_idx) + ztra

         ENDDO

       ENDIF




       trcdmpparttime = MPI_WTIME() - trcdmpparttime ! F79 cronometer-stop
       trcdmptottime  = trcdmptottime + trcdmpparttime



      END SUBROUTINE trcdmp
