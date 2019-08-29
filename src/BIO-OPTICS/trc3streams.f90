
      SUBROUTINE trc3streams(datestring)
!!!---------------------------------------------------------------------
!!!
!!!                       ROUTINE trc3streams
!!!                     ******************

       USE myalloc
       USE mpi
       USE OPT_mem
       USE TIME_MANAGER

       IMPLICIT NONE

       character(LEN=17), INTENT(IN) ::  datestring

!!! local declarations
!!! ==================

#if defined key_trc_nnpzddom || defined key_trc_npzd || key_trc_bfm

      INTEGER :: jk,jj,ji,jl,bottom
      INTEGER :: MODE ! 0-exact, 1-approx
      INTEGER :: year, month, day
      INTEGER :: it_actual
      double precision :: rmud
      double precision :: Edz(jpk,nlt),Esz(jpk,nlt),Euz(jpk,nlt)
      double precision :: CHLz(jpk,nchl),CDOMz(jpk),NAPz(jpk)
      double precision :: Eu_0m(nlt)
      double precision :: sec


      trcoptparttime = MPI_WTIME() ! cronometer-start

! Compute RT every hour

      if (datestring .eq. DATESTART) then
         call read_date_string(datestring, year, month, day, sec)
         it_check = int(sec/3600.d0) 
      endif

      call read_date_string(datestring, year, month, day, sec)
      it_actual = int(sec/3600.d0)
      
      if (lwp) write(*,*) "it_actual", it_actual
      if (lwp) write(*,*) "it_check", it_check

      if ( (it_actual .eq. it_check) .and. (datestring.ne.DATESTART)) then 

         if (lwp) write(*,*) 'Skip computing RT at : ', datestring

         trcoptparttime = MPI_WTIME() - trcoptparttime ! cronometer-stop
         trcopttottime = trcopttottime + trcoptparttime

         return ! no need to compute RT
      else

         it_check = it_actual

      endif

      
! Start computing  RT when and where needed

      do ji=1,jpi
         do jj=1,jpj

! Controls to avoid  calc where no biology
             if (bfmmask(1,jj,ji) == 0) CYCLE

!to be completed       
             zenith_angle = 0.0d0
!to be completed       
             call getrmud(zenith_angle,rmud)

             bottom = mbathy(jj,ji)
             bottom = min(bottom,37) ! Stop at approx 500 mt

             if ( gdept(bottom) > 200.0d0) then 
!               MODE = 1 ! TEST
                MODE = 0 ! exact solution no problem on sea floor reflectance
             else
!               MODE = 0 ! TEST
                MODE = 1 ! approximate solution
             endif


             if ((maxval(Ed_0m(:,jj,ji)) < 0.0001d0) .AND. (maxval(Es_0m(:,jj,ji))< 0.0001d0)) then
  
                 Ed(1,jj,ji,:) = Ed_0m(:,jj,ji)
                 Es(1,jj,ji,:) = Es_0m(:,jj,ji)

                 if (MODE .EQ. 0)  Eu(1,jj,ji,:) = 0.0d0
                 if (MODE .EQ. 1)  Eu(1,jj,ji,:) = -1.0d0
  
                 Ed(2:bottom,jj,ji,:) = 0.0001d0
                 Es(2:bottom,jj,ji,:) = 0.0001d0
                 Eu(2:bottom,jj,ji,:) = 0.0001d0

             else
             
                 CHLz(1:bottom,1) =  trn(1:bottom,jj,ji,ppP1l)
                 CHLz(1:bottom,2) =  trn(1:bottom,jj,ji,ppP2l)
                 CHLz(1:bottom,3) =  trn(1:bottom,jj,ji,ppP3l)
                 CHLz(1:bottom,4) =  trn(1:bottom,jj,ji,ppP4l)
                
                 CDOMz(1:bottom)  =  trn(1:bottom,jj,ji,ppR1c) !!!! WARNING TO BE FIXED
    
                 NAPz(1:bottom)   =  trn(1:bottom,jj,ji,ppR6c) 
    
                 call edeseu(MODE,bottom,e3t(:,jj,ji),Ed_0m(:,jj,ji),Es_0m(:,jj,ji),CHLz,CDOMz,NAPz,rmud,Edz,Esz,Euz,Eu_0m)
    
                 Ed(1,jj,ji,:) = Ed_0m(:,jj,ji)
                 Es(1,jj,ji,:) = Es_0m(:,jj,ji)
                 Eu(1,jj,ji,:) = Eu_0m(:)


                 do jl=1, nlt
                    do jk =2, bottom
                        Ed(jk,jj,ji,jl) = Edz(jk-1,jl)
                        Es(jk,jj,ji,jl) = Esz(jk-1,jl)
                        Eu(jk,jj,ji,jl) = Euz(jk-1,jl)
                    enddo
                 enddo

             endif
            
         enddo
      enddo


      trcoptparttime = MPI_WTIME() - trcoptparttime ! cronometer-stop
      trcopttottime = trcopttottime + trcoptparttime

#else

!!    No optical model

#endif

      END SUBROUTINE trc3streams
