
      ! * contains Times of:
      ! * Forcings
      ! * Restarts
      ! * Ave
      ! * Start and finish of simulation
      MODULE TIME_MANAGER
      use calendar
      IMPLICIT NONE
      integer, parameter :: Time_Max = 5000

      TYPE TIME_CONTAINER
          CHARACTER(LEN=50)          :: NAME
          INTEGER                    :: N
          CHARACTER(LEN=1024)        :: Filename
          CHARACTER(LEN=17), DIMENSION(Time_Max) :: TimeStrings
          CHARACTER(LEN=17), DIMENSION(Time_Max) :: TimeStringsExtended
          real(8),           DIMENSION(Time_Max) :: Times
          INTEGER                                :: Before
          INTEGER                                :: After
          LOGICAL                                :: Periodic
      END TYPE TIME_CONTAINER

      TYPE DUMP_CONTAINER
          CHARACTER(LEN=50)                      :: NAME
          INTEGER                                :: N
          CHARACTER(LEN=1024)                    :: Filename
          CHARACTER(LEN=17), DIMENSION(Time_Max) :: TimeStrings
      END TYPE DUMP_CONTAINER



      PUBLIC :: Load_Timestrings, CheckStartEnd
      !PRIVATE :: CheckInForcings

      TYPE (TIME_CONTAINER) :: TC_FOR
      TYPE (TIME_CONTAINER) :: TC_TIN
      TYPE (TIME_CONTAINER) :: TC_ATM
      TYPE (TIME_CONTAINER) :: TC_GIB
      TYPE (TIME_CONTAINER) :: TC_LEX

      TYPE (DUMP_CONTAINER) :: RESTARTS
      TYPE (DUMP_CONTAINER) :: AVE


      CHARACTER(LEN=17) :: DATESTART
      CHARACTER(LEN=17) :: DATE__END
      INTEGER           :: TimeStepStart
      INTEGER           :: TimeStep__End

      REAL(8)           :: DELTAT    = 1800.0
      REAL(8)           :: TIME_0


      INTEGER           :: TauAVEfrom
      CHARACTER(LEN=17) :: NOW_datestring
      INTEGER           :: NOW_timestep
      REAL(8)           :: NOW_sec
      CHARACTER(LEN=17) :: BKPdatefrom
      LOGICAL           :: IsStartBackup = .false.




      CONTAINS
      ! *************************************************
      ! ** true if the datestring is in the restart list
      ! *************************************************
      LOGICAL FUNCTION IsaRestart(datestring)
      IMPLICIT NONE
      CHARACTER(LEN=17), INTENT(IN) :: datestring
          ! LOCAL
          INTEGER I

          IsaRestart = .false.




          DO I=1, RESTARTS%N
        if (datestring.eq.RESTARTS%TimeStrings(I)) then
            IsaRestart = .true.
            CYCLE
        endif
          ENDDO

         if (datestring.eq.DATESTART) IsaRestart = .false.





      END FUNCTION IsaRestart

      ! *************************************************
      ! ** true if the datestring is in the ave list
      ! *************************************************
      LOGICAL FUNCTION IsAnAveDump(datestring)
      CHARACTER(LEN=17), INTENT(IN) :: datestring

      ! LOCAL
      INTEGER I

      IsAnAveDump = .false.

      DO I=1, AVE%N
          if (datestring.eq.AVE%TimeStrings(I)) then
        IsAnAveDump = .true.
        CYCLE
          endif
      ENDDO

      if (datestring.eq.DATESTART) IsAnAveDump = .false.
      END FUNCTION IsAnAveDump

      ! *************************************************
      ! ** true if there is consistency between
      ! - datestart and date_end
      ! - datestart and forcings
      ! - dateend   and forcings
      ! *************************************************
      LOGICAL FUNCTION CheckStartEnd()
      use calendar
      IMPLICIT NONE

      ! LOCAL
      real(8) t1,t2
      LOGICAL B


      t1= datestring2sec(DATESTART)
      t2= datestring2sec(DATE__END)

      B = .true.

      if (t1.gt.t2) then
        write(*,*) 'ERROR. Datestart follows dateend'
        B = .false.
      endif



      if (.not.CheckDatelist(DATESTART,TC_FOR)) then ; B = .false.; endif
      if (.not.CheckDatelist(DATESTART,TC_TIN)) then ; B = .false.; endif
      if (.not.CheckDatelist(DATESTART,TC_ATM)) then ; B = .false.; endif
      if (.not.CheckDatelist(DATESTART,TC_GIB)) then ; B = .false.; endif
      if (.not.CheckDatelist(DATESTART,TC_LEX)) then ; B = .false.; endif


      if (.not.CheckDatelist(DATE__END,TC_FOR)) then ; B = .false.; endif
      if (.not.CheckDatelist(DATE__END,TC_TIN)) then ; B = .false.; endif
      if (.not.CheckDatelist(DATE__END,TC_ATM)) then ; B = .false.; endif
      if (.not.CheckDatelist(DATE__END,TC_GIB)) then ; B = .false.; endif
      if (.not.CheckDatelist(DATE__END,TC_LEX)) then ; B = .false.; endif

      CheckStartEnd = B
      END FUNCTION CheckStartEnd






      ! *************************************************
      LOGICAL FUNCTION CheckDatelist(datestring, STRUCT)
      USE calendar
      IMPLICIT NONE
      CHARACTER(LEN=17),     INTENT(IN) :: datestring
      TYPE (TIME_CONTAINER), INTENT(IN) :: STRUCT

      ! LOCAL
      real(8) t

      IF (STRUCT%Periodic) THEN
         CheckDatelist = .true.
      ELSE
          t = datestring2sec(datestring)
          CheckDatelist = .true.

          if (t.lt.STRUCT%Times(1)) then
        CheckDatelist=.false.
        write(*,*) 'DateStart not in ', STRUCT%NAME, 'dates'
          else
        if (t.gt.STRUCT%Times( STRUCT%N )) then
           CheckDatelist=.false.
           write(*,*) 'DateEnd not in ', STRUCT%NAME, 'dates'
        endif

          endif

      ENDIF

      END FUNCTION CheckDatelist

      SUBROUTINE Load_Dump_container(STRUCT)
          integer I,N
          TYPE (DUMP_CONTAINER) STRUCT
          integer TheUnit

          TheUnit = 326
          call countline(STRUCT%Filename,N)
          if ((N+2).gt.Time_Max) then
             write(*,*) 'Increase Time_Max, please'
             STOP
          endif
          STRUCT%N = N

!          ALLOCATE(STRUCT%Timestrings(N))

          OPEN(UNIT=TheUnit,file=STRUCT%Filename,status='old')
          DO I=1,N
            read(TheUnit,'(A)') STRUCT%TimeStrings(I)
          ENDDO
          CLOSE(TheUnit)

      END SUBROUTINE Load_Dump_container


      ! *************************************************
      SUBROUTINE Load_Time_container(STRUCT)
      integer I,N
      TYPE (TIME_CONTAINER) STRUCT
      integer TheUnit

      TheUnit = 326
      call countline(STRUCT%Filename,N)

      if ((N+2).gt.Time_Max) then
         write(*,*) 'Increase Time_Max, please'
         STOP
      endif

      STRUCT%N = N
!      ALLOCATE(STRUCT%Timestrings(N))


      OPEN(UNIT=TheUnit,file=STRUCT%Filename,status='old')
      DO I=1,N
          read(TheUnit,'(A)') STRUCT%TimeStrings(I)
      ENDDO
      CLOSE(TheUnit)

       ! Periodicità ( mensile, stagionale, ... )
      if (STRUCT%TimeStrings(1)(1:4).eq.'yyyy') then
          STRUCT%Periodic = .true.
!          ALLOCATE(STRUCT%Times        (1:N+2))
!          ALLOCATE(STRUCT%TimeStringsExtended(1:N+2))
      else
          STRUCT%Periodic = .false.
!          ALLOCATE(STRUCT%Times(1:N))
          DO I=1,N
              STRUCT%Times(I)= datestring2sec(STRUCT%TimeStrings(I))
          ENDDO
      endif

      END SUBROUTINE Load_Time_container



      ! ****************************************************
      !** Lettura dei tempi di :
      ! * Inizio e fine simulazione,
      ! * Forzanti, ave e restarts
      !** E caricamento nel modulo
      ! ****************************************************
      SUBROUTINE Load_Timestrings
      USE calendar
      IMPLICIT NONE

      character(LEN=1024) FileName
      integer TheUnit

      TheUnit = 326

      TC_FOR%Filename    = 'forcingsTimes'
      TC_FOR%Name        = 'Forcings'

      TC_TIN%Filename    = 'RiversTimes'
      TC_TIN%Name        = 'Terrestrial Inputs'

      TC_ATM%Filename    = 'AtmTimes'
      TC_ATM%Name        = 'Atmospherical'

      TC_GIB%Filename    = 'GibTimes'
      TC_GIB%Name        = 'Gibraltar'

      RESTARTS%Filename  = 'restartTimes'
      RESTARTS%Name      = 'BFM files'

      AVE%Filename       = 'aveTimes'
      TC_LEX%Filename    = 'kextTimes'

      FileName = 'Start_End_Times'
      OPEN(UNIT=TheUnit,file=FileName,status='old')
      read(TheUnit,'(A)') DATESTART
      read(TheUnit,'(A)') DATE__END
      CLOSE(TheUnit)

      call Load_Time_container(TC_FOR)
      call Load_Time_container(TC_TIN)
      call Load_Time_container(TC_ATM)
      call Load_Time_container(TC_GIB)
      call Load_Time_container(TC_LEX)

      ! DUMP
      call Load_Dump_container(RESTARTS)
      call Load_Dump_container(AVE)


      END SUBROUTINE Load_Timestrings


      ! potrei fare
      ! tau2datestring e poi datestring2sec
      ! ma evito il passaggio alle stringhe
      REAL(8) FUNCTION  tau2sec(TAU)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: TAU
      tau2sec = TAU*deltaT + TIME_0
      END FUNCTION tau2sec


      SUBROUTINE tau2datestring(TAU, datestring)
      use calendar
      IMPLICIT NONE

      INTEGER, INTENT(IN):: TAU
      CHARACTER(LEN=17), INTENT(OUT) :: datestring

      ! LOCAL
      integer :: year, month,day
      real(8) :: sec

         call itau2ymds(TAU,deltaT, year, month,day,sec)
         call write_date_string(datestring, year, month, day, sec)

      END SUBROUTINE tau2datestring





      SUBROUTINE TimeExtension(datestring, STRUCT)
      IMPLICIT NONE
      CHARACTER(LEN=17),INTENT( IN) :: datestring
      TYPE (TIME_CONTAINER) STRUCT
      character(LEN=4) yyyy
      integer year,I, LAST,N

      N   = STRUCT%N
      LAST= STRUCT%N+2

      if (STRUCT%Periodic) then

          read(datestring,'(I4)') year
          write(yyyy,'(I4)') year-1 ;
          STRUCT%TimeStringsExtended(1   ) = yyyy//STRUCT%TimeStrings(N)(5:17)

          DO I=2, N+1
        STRUCT%TimeStringsExtended(I) = datestring(1:4)//STRUCT%TimeStrings(I-1)(5:17)
          ENDDO

          write(yyyy,'(I4)') year+1 ;
          STRUCT%TimeStringsExtended(LAST) = yyyy//STRUCT%TimeStrings(1)(5:17)

          ! now, time in sec
          DO I = 1,LAST
        STRUCT%Times(I) = datestring2sec(STRUCT%TimeStringsExtended(I))
          ENDDO

      endif

      END SUBROUTINE TimeExtension



      ! ****************************************
      SUBROUTINE TimeInterpolation(sec, STRUCT,Before, After, t_interp)
      IMPLICIT NONE
      REAL(8),        INTENT(IN ) :: sec
      TYPE (TIME_CONTAINER),INTENT(IN ) :: STRUCT
      INTEGER,        INTENT(OUT) :: Before, After
      REAL(8),        INTENT(OUT) :: t_interp


      ! LOCAL
      integer I, N

      N=STRUCT%N

         I =1 ;
         DO WHILE (STRUCT%Times(I).lt.sec)
           I = I + 1;
         ENDDO
         BEFORE = I - 1;

      AFTER = BEFORE +1
      t_interp = (sec - STRUCT%Times(BEFORE))/(STRUCT%Times(AFTER) - STRUCT%Times(BEFORE))


      IF (STRUCT%Periodic) then
          ! ritorno ai 12
          BEFORE = BEFORE -1 ;
          AFTER  = BEFORE +1
          if (BEFORE.eq.0   ) then ; BEFORE = N; endif
          if (AFTER.eq.(N+1)) then ; AFTER  = 1; endif
      ENDIF



      END SUBROUTINE TimeInterpolation

      ! *******************************************************
      SUBROUTINE YEARLY(datestring)
      IMPLICIT NONE
      CHARACTER(LEN=17),INTENT( IN) :: datestring

      if (datestring(5:17).eq.'0101-00:00:00') then
         call TimeExtension(datestring,TC_FOR)
         call TimeExtension(datestring,TC_TIN)
         call TimeExtension(datestring,TC_ATM)
         call TimeExtension(datestring,TC_GIB)
         call TimeExtension(datestring,TC_LEX)
      endif
      END SUBROUTINE YEARLY

! *******************************************************

      SUBROUTINE DAILY(datestring)
       USE myalloc
       USE CALENDAR
      IMPLICIT NONE

      CHARACTER(LEN=17),INTENT( IN) :: datestring
      ! LOCAL
      integer julianday, ji,jj

      if (datestring(10:17).eq.'00:00:00') then

         call tau2julianday(NOW_timestep, deltaT, julianday)
         do jj =1, jpj
            do ji=1, jpi
               DAY_LENGTH(ji,jj) = photoperiod(julianday, gphit(ji,jj))
            enddo
         enddo


      endif
      END SUBROUTINE DAILY


! **************************************************************

      SUBROUTINE COUNTLINE(FILENAME,LINES)
          implicit none
          character FILENAME*(*)
          integer lines
          integer TheUnit

          TheUnit = 326

          lines=0
          OPEN(UNIT=TheUnit,file=FILENAME,status='old')
          DO WHILE (.true.)
       read(TheUnit, *, END=21)
       lines = lines+1
          ENDDO

21        CLOSE(TheUnit)

      END SUBROUTINE COUNTLINE


      ! ** gets initial and final timestep, integers to
      ! ** implement the main cycle
      SUBROUTINE getTimesteps(timestep1, timestep2)
          use calendar
          IMPLICIT NONE
          integer , INTENT(OUT) :: timestep1, timestep2

          ! LOCAL
          integer year1,  year2
          integer month1, month2
          integer day1,   day2
          real(8) sec1, sec2

          real(8) sec_diff, sec0, julian0
          integer year0, month0, day0


          call read_date_string(DATESTART, year1, month1, day1, sec1)
          call read_date_string(DATE__END, year2, month2, day2, sec2)

          ! ***** Tempo zero = 1 gennaio (anno di startdate) *********
          year0=year1;
          month0=1
          day0=1
          sec0=0.
          call ioconf_startdate(year0, month0, day0, sec0)
          call ymds2ju(year0, month0, day0, sec0, julian0)
          TIME_0 = julian0*86400.0
          ! *********************************************************

          call time_diff(year0, month0, day0, sec0, year1, month1,day1,sec1, sec_diff)
          timestep1 = NINT(sec_diff/deltaT);
          call time_diff(year0, month0, day0, sec0, year2, month2,day2,sec2, sec_diff)
          timestep2 = NINT(sec_diff/deltaT);


      END SUBROUTINE getTimesteps

       INTEGER FUNCTION datestringToTAU(datestring)
       use calendar
       IMPLICIT NONE
       character(LEN=17), INTENT(IN) :: datestring

       !local
       integer year, year0
       integer month, month0
       integer day, day0
       real(8) sec, sec0, sec_diff

       datestringToTAU = -500

       call read_date_string(datestring, year, month, day, sec)
       call ju2ymds(TIME_0/86400.0, year0, month0, day0, sec0)
       call time_diff(year0, month0, day0, sec0, year, month,day,sec, sec_diff)


       datestringToTAU = NINT(sec_diff/deltaT);


       END FUNCTION datestringToTAU


      SUBROUTINE MIDDLEDATE(Tau1, Tau2, datestring)
        USE calendar
        IMPLICIT NONE
        integer,           INTENT(IN ) :: Tau1, Tau2
        CHARACTER(LEN=17), INTENT(OUT) :: datestring

        ! local

        real(8) sec,  seconds
        integer year, month, day

      sec = (tau2sec(tau1) + tau2sec(tau2)) /2
      call ju2ymds(sec/86400., year, month, day, seconds)
      seconds = real(nint(seconds),8)
      call write_date_string(datestring, year, month, day, seconds)

      END SUBROUTINE MIDDLEDATE


      END MODULE TIME_MANAGER
