
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

      INTEGER :: jk,jj,ji,jl,bottom, phys_bottom
      INTEGER :: day_of_year
      INTEGER :: MODE ! 0-exact, 1-approx
      INTEGER :: year, month, day, ihr
      INTEGER :: it_actual
      CHARACTER(LEN=20) :: V_POSITION     
      double precision :: solz(jpj,jpi), rmud(jpj,jpi)
      double precision :: Edz(jpk,nlt),Esz(jpk,nlt),Euz(jpk,nlt),SWRz(jpk)
      double precision, allocatable::  E(:,:,:) !(3,jpk+1,nlt)

      double precision, allocatable  :: PARz(:,:) !, SWRz(:)   !(jpk,nchl+1)
      double precision :: CHLz(jpk,nchl),PCz(jpk,nchl),CDOMz(jpk,3),POCz(jpk)
      double precision :: Eu_0m(nlt)
      double precision :: zgrid(jpk+1)
      double precision :: sec
      integer j3




! Where to compute irradiance options are
!- BOTTOM
!- MID
!- AVERAGE

      V_POSITION = "AVERAGE"


      call tau2julianday(datestringToTAU(datestring), deltaT, day_of_year)
     
      if (lwp) write(*,*) 'day_of_year', day_of_year

      call sfcsolz(year, day_of_year, ihr, solz) 

      call getrmud(solz,rmud)



! Start computing  RT when and where needed

      do ji=1,jpi
         do jj=1,jpj

! Controls to avoid  calc where no biology
             if (bfmmask(1,jj,ji) == 0) CYCLE
         

             phys_bottom = mbathy(jj,ji)
             bottom = min(phys_bottom,jpk_opt) ! Stop at approx 500 mt

             allocate( E(3,bottom+1,nlt))
             allocate( PARz(bottom,nchl+1))
             E(:,:,:)  = 0.0001d0
             PARz(:,:) = 0.0001d0

!            MODE = 0 ! exact solution
!            MODE = 1 ! approximate solution
             MODE = 2 ! library solution



             RMU(jj,ji) = rmud(jj,ji)

       
             if ((maxval(Ed_0m(:,jj,ji)) < 0.0001d0) .AND. (maxval(Es_0m(:,jj,ji))< 0.0001d0)) then
  
  !               Ed(1,jj,ji,:) = Ed_0m(:,jj,ji)
  !               Es(1,jj,ji,:) = Es_0m(:,jj,ji)

                 if (MODE .EQ. 0)  Eu(1,jj,ji,:) = 1.0E-08
                 if (MODE .EQ. 1)  Eu(1,jj,ji,:) = -1.0d0
                 if (MODE .EQ. 2)  Eu(1,jj,ji,:) = 1.0E-08
  
                 Ed(1:bottom,jj,ji,:) = 0.0001d0
                 Es(1:bottom,jj,ji,:) = 0.0001d0
                 Eu(2:bottom,jj,ji,:) = 1.0E-08
                 PAR(1:bottom,jj,ji,:) = 0.0001d0
                 SWR_RT(1:bottom,jj,ji) = 0.0001d0

             else
             
                 CHLz(1:bottom,1) = trn(1:bottom,jj,ji,ppP1l)
                 CHLz(1:bottom,2) = trn(1:bottom,jj,ji,ppP2l)
                 CHLz(1:bottom,3) = trn(1:bottom,jj,ji,ppP3l)
                 CHLz(1:bottom,4) = trn(1:bottom,jj,ji,ppP4l)

                 PCz(1:bottom,1) = trn(1:bottom,jj,ji,ppP1c)
                 PCz(1:bottom,2) = trn(1:bottom,jj,ji,ppP2c)
                 PCz(1:bottom,3) = trn(1:bottom,jj,ji,ppP3c)
                 PCz(1:bottom,4) = trn(1:bottom,jj,ji,ppP4c)
                 
                 CDOMz(1:bottom,1)  = trn(1:bottom,jj,ji,ppR1l)
                 CDOMz(1:bottom,2)  = trn(1:bottom,jj,ji,ppR2l)
                 CDOMz(1:bottom,3)  = trn(1:bottom,jj,ji,ppR3l)

                 POCz(1:bottom)   = trn(1:bottom,jj,ji,ppR6c) 
    
                 IF ( (MODE .EQ. 0) .OR. (MODE .EQ. 1)) then
                     call edeseu(MODE,V_POSITION,bottom,e3t(:,jj,ji),Ed_0m(:,jj,ji),Es_0m(:,jj,ji),CHLz,PCz,CDOMz,POCz,rmud(jj,ji),Edz,Esz,Euz,Eu_0m,PARz,SWRz)
    
                 Ed(1,jj,ji,:) = Ed_0m(:,jj,ji)
                 Es(1,jj,ji,:) = Es_0m(:,jj,ji)
                 Eu(1,jj,ji,:) = Eu_0m(:)



                 do jl=1, nlt
                    do jk =2, bottom
                        Ed(jk,jj,ji,jl) = Edz(jk-1,jl)
                        Es(jk,jj,ji,jl) = Esz(jk-1,jl)
                        Eu(jk,jj,ji,jl) = Euz(jk-1,jl)
!                       write(*,*) "Ed", jl,jk,"=", Ed(jk,jj,ji,jl)

                    enddo
                 enddo


                 do jl=1, nchl+1
                    do jk =1, bottom
                        PAR(jk,jj,ji,jl) = PARz(jk,jl)
!                       write(*,*) "PAR", jl,jk,"=", PAR(jk,jj,ji,jl)
                    enddo
                 enddo

                do jk =1, bottom
                    SWR_RT(jk,jj,ji) = SWRz(jk)
                enddo


                 ENDIF

                 IF (MODE .EQ. 2) then
                     zgrid(1)=0.0D0
                     do jk =1,jpk
                         zgrid(jk+1) = zgrid(jk) + e3t(jk,jj,ji)
                     enddo
                    




                     call edeseu2(MODE,V_POSITION,bottom,zgrid,Ed_0m(:,jj,ji),Es_0m(:,jj,ji), &
                                  CHLz,PCz,CDOMz,POCz,rmud(jj,ji),E,PARz)
                
                 do jl=1, nlt
                    do jk =1, bottom+1 ! Defined on w faces (cell's interfaces)
                        Ed(jk,jj,ji,jl) = E(1,jk,jl)
                        Es(jk,jj,ji,jl) = E(2,jk,jl)
                        Eu(jk,jj,ji,jl) = E(3,jk,jl)

                    enddo

                    if ((bottom+1).lt.phys_bottom) then

                       do jk=bottom+1, phys_bottom+1 ! propagation 500m --> phys_bottom of last value
                          Ed(jk,jj,ji,jl) = Ed(bottom+1, jj,ji, jl)
                          Es(jk,jj,ji,jl) = Es(bottom+1, jj,ji, jl)
                          Eu(jk,jj,ji,jl) = Eu(bottom+1, jj,ji, jl)
                       enddo
                    endif
                 enddo

                 do jl=1, nchl+1
                    do jk =1, bottom
                        PAR(jk,jj,ji,jl) = PARz(jk,jl)
                    enddo
                    if (bottom.lt.phys_bottom) then
                        do jk=bottom+1, phys_bottom ! propagation of last value
                          PAR(jk,jj,ji,jl) = PARz(bottom,jl)
                        enddo
                    endif

                 enddo 


                 ENDIF

                 

             endif
            
            deallocate(E)
            deallocate(PARz)
         enddo
      enddo



#else

!!    No optical model

#endif

      END SUBROUTINE trc3streams
