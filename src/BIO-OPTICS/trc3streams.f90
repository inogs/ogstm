
      SUBROUTINE trc3streams()
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trcopt
!!!                     ******************

       USE myalloc
       USE mpi
       USE OPT_mem
       IMPLICIT NONE


!!! local declarations
!!! ==================

#if defined key_trc_nnpzddom || defined key_trc_npzd || key_trc_bfm

      INTEGER :: jk,jj,ji,jl,bottom
      INTEGER :: MODE ! 0-exact, 1-approx
      double precision :: rmud
      double precision :: Edz(jpk,nlt),Esz(jpk,nlt),Euz(jpk,nlt)
      double precision :: CHLz(jpk,nchl),CDOMz(jpk),NAPz(jpk)
      double precision :: Eu_0m(nlt)


      trcoptparttime = MPI_WTIME() ! cronometer-start

      do ji=1,jpi
         do jj=1,jpj
       
             zenith_angle = 0.
             call getrmud(zenith_angle,rmud)

             if (bfmmask(1,jj,ji) == 0) CYCLE

             bottom = mbathy(jj,ji)
             bottom = min(bottom,37) ! Stop at approx 500 mt
             
             if ( gdept(bottom) > 200.) then 
!               MODE = 1 ! TEST
                MODE = 0 ! exact solution no problem on sea floor reflectance
             else
                MODE = 1 ! approximate solution
             endif

             CHLz(1:bottom,1) =  trn(1:bottom,jj,ji,ppP1l)
             CHLz(1:bottom,2) =  trn(1:bottom,jj,ji,ppP2l)
             CHLz(1:bottom,3) =  trn(1:bottom,jj,ji,ppP3l)
             CHLz(1:bottom,4) =  trn(1:bottom,jj,ji,ppP4l)
            
             CDOMz(1:bottom)  = trn(1:bottom,jj,ji,ppR1c) !!!! WARNING TO BE FIXED

             NAPz(1:bottom)  = trn(1:bottom,jj,ji,ppR6c) 

             call edeseu(MODE,bottom,e3t(:,jj,ji),Ed_0m(:,jj,ji),Es_0m(:,jj,ji),CHLz,CDOMz,NAPz,rmud,Edz,Esz,Euz,Eu_0m)

             Ed(1,jj,ji,:) = Ed_0m(:,jj,ji)
             Es(1,jj,ji,:) = Es_0m(:,jj,ji)
             Eu(1,jj,ji,:) = Eu_0m(:)

             do jk =2, jpk-1
                do jl=1, nlt
                    Ed(jk,jj,ji,jl) = Edz(jk-1,jl)
                    Es(jk,jj,ji,jl) = Esz(jk-1,jl)
                    Eu(jk,jj,ji,jl) = Euz(jk-1,jl)
                enddo
             enddo
            
         enddo
      enddo


      trcoptparttime = MPI_WTIME() - trcoptparttime ! cronometer-stop
      trcopttottime = trcopttottime + trcoptparttime

#else

!!    No optical model

#endif

      END SUBROUTINE trc3streams
