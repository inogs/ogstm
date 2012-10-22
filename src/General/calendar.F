!
      SUBROUTINE histerr (plev, pcname, pstr1, pstr2, pstr3)
!---------------------------------------------------------------------
!- INPUT
!- plev   : Category of message to be reported to the user
!-          1 = Note to the user
!-          2 = Warning to the user
!-          3 = Fatal error
!- pcname : Name of subroutine which has called histerr
!- pstr1
!- pstr2  : String containing the explanations to the user
!- pstr3
!---------------------------------------------------------------------
         IMPLICIT NONE
!-
         INTEGER :: plev
         CHARACTER(LEN=*) :: pcname,pstr1,pstr2,pstr3
!-
         CHARACTER(LEN=30),DIMENSION(3) :: pemsg = (/ "NOTE TO THE USER FROM ROUTINE ",
     &     "WARNING FROM ROUTINE          ", 
     &     "FATAL ERROR FROM ROUTINE      " /)
!---------------------------------------------------------------------
         IF ( (plev >= 1).AND.(plev <= 3) ) THEN
           WRITE(*,'("     ")')
           WRITE(*,'(A," ",A)') TRIM(pemsg(plev)),TRIM(pcname)
           WRITE(*,'(" --> ",a)') pstr1
           WRITE(*,'(" --> ",a)') pstr2
           WRITE(*,'(" --> ",a)') pstr3
         ENDIF
         IF (plev == 3) THEN
           STOP 'Fatal error from CALENDAR. See stdout for more details'
         ENDIF
!------------------------
         END SUBROUTINE histerr


      MODULE calendar
        !
        !$Header: /cvsroot/opatm-bfm/opa_model/OPA/calendar.F,v 1.2 2009-09-11 09:20:55 cvsogs01 Exp $
        !
        !
        !   This is the calendar which going to be used to do all calculations
        !   on time. Three types of calendars are possible :
        !   - gregorian : The normal calendar. The time origin for the julian day in this case is
        !     24 Nov -4713
        !   - nolap : A 365 day year without leap years. The origin for the julian days
        !     is in this case 1 Jan 0
        !   - xxxd : Year of xxx days with month of equal length. The origin for the julian
        !     days is then also 1 Jan 0
        !   As one can see it is difficult to go from one calendar to the other. All operations
        !   involving julian days will be wrong. This calendar will lock as soon as possible 
        !   the length of the year and  forbid any further modification.
        !
        !   For the non leap-year calendar the method is still brute force. We need to find
        !   an Integer series which takes care of the length of the various month. (Jan)
        !
        !      un_jour : one day in seconds
        !      un_an   : one year in days
        !
        USE stringop,  ONLY : strlowercase
        ! 
        PRIVATE
        PUBLIC :: ymds2ju, ju2ymds, tlen2itau, isittime,
     &   ioconf_calendar, ioget_calendar, itau2date, ioget_timestamp,
     & ioconf_startdate, itau2ymds, time_add, time_diff,
     & write_date_string, read_date_string, datestring2sec,photoperiod,
     & tau2julianday
        !
        INTERFACE ioget_calendar
          MODULE PROCEDURE ioget_calendar_real1
          MODULE PROCEDURE ioget_calendar_real2
          MODULE PROCEDURE ioget_calendar_str
        END INTERFACE
      

        INTERFACE ioconf_startdate
           MODULE PROCEDURE ioconf_startdate_simple
           MODULE PROCEDURE ioconf_startdate_internal
           MODULE PROCEDURE ioconf_startdate_ymds
        END INTERFACE
        !
        REAL(8), PARAMETER :: un_jour = 86400.0
        REAL(8), SAVE :: un_an = 365.2425
        LOGICAL, SAVE :: lock_unan = .FALSE.
        LOGICAL, SAVE :: lock_startdate = .FALSE.
        !
        CHARACTER(LEN=30), SAVE :: time_stamp='XXXXXXXXXXXXXXXX'
        !
        INTEGER, PARAMETER :: mon_len(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)

        CHARACTER(LEN=3),PARAMETER::cal(12)=(/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)
        !
        REAL(8), SAVE :: start_day, start_sec
        !
      CONTAINS
        !
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        !
        SUBROUTINE ymds2ju(year, month, day, sec, julian)
          !
          IMPLICIT NONE
          !
          !     INPUT
          !
          INTEGER, INTENT(IN) :: year, month, day
          REAL(8), INTENT(IN)    :: sec
          !
          !     OUTPUT
          !     
          REAL(8), INTENT(OUT) :: julian
          !
          !     LOCAL
          !
          INTEGER :: julian_day
          REAL(8)    :: julian_sec
          !
        CALL ymds2ju_internal(year, month, day, sec,
     & julian_day, julian_sec)
          !
          julian = julian_day + julian_sec / un_jour
          !
        END SUBROUTINE ymds2ju
        !
        !-----------------------------------------------------------------------
        !
        SUBROUTINE ymds2ju_internal(year, month, day, sec, julian_day, julian_sec)
          !
          IMPLICIT NONE
          !
          !     Converts year, month, day and seconds into a julian day
          !
          !     In 1968 in a letter to the editor of Communications of the ACM
          !     (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
          !     and Thomas C. Van Flandern presented such an algorithm.
          !
          !     See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
          !
          !     In the case of the Gregorian calendar we have chosen to use the Lilian 
          !     day numbers. This is the day counter which starts on the 15th October 1582.
          !     This is the day at which Pope Gregory XIII introduced the Gregorian calendar.
          !     Compared to the true Julian calendar, which starts some 7980 years ago, the Lilian
          !     days are smaler and are dealt with easily on 32 bit machines. With the true Julian 
          !     days you can only the fraction of the day in the real(8) part to a precision of a 1/4 of 
          !     a day with 32 bits.
          !
          !     INPUT
          !
          INTEGER, INTENT(IN) :: year, month, day
          REAL(8), INTENT(IN)    :: sec
          !
          !     OUTPUT
          !     
          INTEGER, INTENT(OUT) :: julian_day
          REAL(8), INTENT(OUT)    :: julian_sec
          !
          !     LOCAL
          !
          INTEGER :: jd, m, y, d, ml
          !
          lock_unan = .TRUE.
          !   
          m = month
          y = year
          d = day
          !
          IF ( un_an .GT. 365.0 .AND. un_an .LT. 366.0) THEN
              !
              jd = ( 1461 * ( y + 4800 + INT(( m - 14 ) / 12 )) ) / 4 + 
     &        ( 367 * ( m - 2 - 12 * ( INT(( m - 14 ) / 12 )) ) ) / 12 - 
     &        ( 3 * ( ( y + 4900 + INT(( m - 14 ) / 12 )) / 100 ) ) / 4 + 
     &        d - 32075
              jd = jd - 2299160
              !          
          ELSE IF ( ABS(un_an - 365.0) .LE. EPSILON(un_an) ) THEN
              !
              ml = SUM(mon_len(1:m-1))
              jd = y*INT(un_an) + ml + (d - 1)
              !
          ELSE
              !
              ml = INT(un_an)/12
              jd = y*INT(un_an) + (m-1)*ml + (d - 1)
              !
          ENDIF
          !
          julian_day = jd - 1 ! -1 by GB
          julian_sec = sec
          !
          RETURN
          !
        END SUBROUTINE ymds2ju_internal
        !
        !=======================================================================
        !
        SUBROUTINE ju2ymds(julian, year, month, day, sec)
          !
          IMPLICIT NONE
          !
          !     INPUT
          !     
          REAL(8), INTENT(IN) :: julian
          !
          !     OUTPUT
          !
          INTEGER, INTENT(OUT) :: year, month, day
          REAL(8), INTENT(OUT)    :: sec
          !
          !     LOCAL
          !
          INTEGER :: julian_day
          REAL(8)    :: julian_sec
          !
          julian_day = INT(julian)
          julian_sec = (julian - julian_day)*un_jour
          !
      CALL ju2ymds_internal(julian_day, julian_sec, year, month,
     & day, sec)
          !
        END SUBROUTINE ju2ymds
        !
        !=======================================================================
        !
        SUBROUTINE ju2ymds_internal(julian_day, julian_sec, year, month, day, sec)
          !
          IMPLICIT NONE
          !
          !     This subroutine computes from the julian day the year, month, day and
          !     seconds
          !
          !     In 1968 in a letter to the editor of Communications of the ACM
          !     (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
          !     and Thomas C. Van Flandern presented such an algorithm.
          !
          !     See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
          !
          !     In the case of the Gregorian calendar we have chosen to use the Lilian 
          !     day numbers. This is the day counter which starts on the 15th October 1582.
          !     This is the day at which Pope Gregory XIII introduced the Gregorian calendar.
          !     Compared to the true Julian calendar, which starts some 7980 years ago, the Lilian
          !     days are smaler and are dealt with easily on 32 bit machines. With the true Julian 
          !     days you can only the fraction of the day in the real(8) part to a precision of a 1/4 of 
          !     a day with 32 bits.
          !
          !     INPUT
          !     
          INTEGER, INTENT(IN) :: julian_day
          REAL(8), INTENT(IN)    :: julian_sec
          !
          !     OUTPUT
          !
          INTEGER, INTENT(OUT) :: year, month, day
          REAL(8), INTENT(OUT)    :: sec
          !
          !     LOCAL
          !
          INTEGER :: l, n, i, jd, j, d, m, y, ml
          INTEGER :: add_day
          !      
          lock_unan = .TRUE.
          !
          jd = julian_day + 1 ! + 1 by GB
          sec = julian_sec
          IF ( sec > un_jour ) THEN
             add_day = INT(sec / un_jour)
             sec = sec - add_day*un_jour
             jd = jd + add_day
          ENDIF
          !
          IF ( un_an .GT. 365.0 .AND. un_an .LT. 366.0) THEN
              !
              jd = jd + 2299160
              !
              l = jd + 68569
              n = ( 4 * l ) / 146097
              l = l - ( 146097 * n + 3 ) / 4
              i = ( 4000 * ( l + 1 ) ) / 1461001 
              l = l - ( 1461 * i ) / 4 + 31
              j = ( 80 * l ) / 2447
              d = l - ( 2447 * j ) / 80
              l = j / 11
              m = j + 2 - ( 12 * l )
              y = 100 * ( n - 49 ) + i + l   
              !
          ELSE IF ( ABS(un_an - 365.0) .LE. EPSILON(un_an) ) THEN
              !
              y = jd/INT(un_an)
              l = jd - y*INT(un_an)
              m = 1
              ml = 0
              DO WHILE (ml + mon_len(m) .LE. l)
                ml = ml + mon_len(m)
                m = m + 1
              ENDDO
              d = l - ml + 1
              !          
          ELSE
              ! 
              ml = INT(un_an)/12
              y = jd/INT(un_an)
              l = jd - y*INT(un_an)
              m = (l/ml)+1
              d = l - (m-1)*ml + 1
              !        
          ENDIF
          !
          day = d
          month = m
          year = y
      
          RETURN
      
        END SUBROUTINE ju2ymds_internal
        !
        !=======================================================================
        !
        SUBROUTINE tlen2itau(input_str, dt, date, itau)
          !
          IMPLICIT NONE
          !
          !     This subroutine transforms a string containing a time length
          !     into a number of time steps.
          !     To do this operation the date (in julian days) is needed as the
          !     length of the month varies.
          !     The following convention is used :
          !
          !         n   : n time steps
          !         nS  : n seconds is transformed into itaus
          !         nH  : n hours
          !         nD  : n days 
          !         nM  : n month  
          !         nY  : n years  
          !     Combinations are also possible
          !         nYmD : nyears plus m days !
          !
          !     INPUT
          !
          CHARACTER(LEN=*), INTENT(IN) :: input_str
          REAL(8), INTENT(IN)             :: dt, date
          !     
          !     OUTPUT
          !
          INTEGER, INTENT(OUT)          :: itau
          !
          !     LOCAL
          !
          INTEGER           :: y_pos, m_pos, d_pos, h_pos, s_pos
          INTEGER           :: read_time
          CHARACTER(LEN=13) :: fmt
          CHARACTER(LEN=80) :: tmp_str
          !
          INTEGER           :: year, month, day
          REAL(8)              :: sec, date_new, dd, ss
          !
          !
          !
          itau = 0
          CALL ju2ymds(date, year, month, day, sec)
          !
          y_pos = MAX(INDEX(input_str,'y'),INDEX(input_str,'Y'))
          m_pos = MAX(INDEX(input_str,'m'),INDEX(input_str,'M'))
          d_pos = MAX(INDEX(input_str,'d'),INDEX(input_str,'D'))
          h_pos = MAX(INDEX(input_str,'h'),INDEX(input_str,'H'))
          s_pos = MAX(INDEX(input_str,'s'),INDEX(input_str,'S'))
          !
          IF ( MAX(y_pos,m_pos,d_pos,s_pos) .GT. 0) THEN
              !
              tmp_str = input_str
              !
              DO WHILE ( MAX(y_pos,m_pos,d_pos,s_pos) .GT. 0)
      
                !            WRITE(*,*) tmp_str
                !            WRITE(*,*) y_pos,m_pos,d_pos,s_pos
      
                IF (y_pos .GT. 0) THEN
                    WRITE(fmt,'("(I",I10.10,")")') y_pos-1
                    READ(tmp_str(1:y_pos-1),fmt) read_time
                CALL ymds2ju(year+read_time, month, day, sec, date_new)
                    dd = date_new-date 
                    ss = INT(dd)*un_jour+ dd-INT(dd)
                    itau = itau + NINT(ss/dt)
                    tmp_str = tmp_str(y_pos+1:LEN_TRIM(tmp_str))
                ELSE IF (m_pos .GT. 0) THEN
                    WRITE(fmt,'("(I",I10.10,")")') m_pos-1
                    READ(tmp_str(1:m_pos-1),fmt) read_time
                CALL ymds2ju(year, month+read_time, day, sec, date_new)
                    dd = date_new-date 
                    ss = INT(dd)*un_jour+ dd-INT(dd)
                    itau = itau + NINT(ss/dt)
                    tmp_str = tmp_str(m_pos+1:LEN_TRIM(tmp_str))
                ELSE IF (d_pos .GT. 0) THEN
                    WRITE(fmt,'("(I",I10.10,")")') d_pos-1
                    READ(tmp_str(1:d_pos-1),fmt) read_time
                    itau = itau + NINT(read_time*un_jour/dt)
                    tmp_str = tmp_str(d_pos+1:LEN_TRIM(tmp_str))
                ELSE IF (h_pos .GT. 0) THEN
                    WRITE(fmt,'("(I",I10.10,")")') h_pos-1
                    READ(tmp_str(1:h_pos-1),fmt) read_time
                    itau = itau + NINT(read_time*60.*60./dt)
                    tmp_str = tmp_str(d_pos+1:LEN_TRIM(tmp_str))
                ELSE IF  (s_pos .GT. 0) THEN
                    WRITE(fmt,'("(I",I10.10,")")') s_pos-1
                    READ(tmp_str(1:s_pos-1),fmt) read_time
                    itau = itau + NINT(read_time/dt)
                    tmp_str = tmp_str(s_pos+1:LEN_TRIM(tmp_str))
                ENDIF
                !
                y_pos = MAX(INDEX(tmp_str,'y'),INDEX(tmp_str,'Y'))
                m_pos = MAX(INDEX(tmp_str,'m'),INDEX(tmp_str,'M'))
                d_pos = MAX(INDEX(tmp_str,'d'),INDEX(tmp_str,'D'))
                h_pos = MAX(INDEX(tmp_str,'h'),INDEX(tmp_str,'H'))
                s_pos = MAX(INDEX(tmp_str,'s'),INDEX(tmp_str,'S'))
                !
      
                !
              ENDDO
              !
          ELSE
              WRITE(fmt,'("(I",I10.10,")")') LEN_TRIM(input_str)
              READ(input_str(1:LEN_TRIM(input_str)),fmt) itau
          ENDIF
          !
      
      
          RETURN
      
        END SUBROUTINE tlen2itau
      
        !
        !=======================================================================
        !
        !
        REAL(8) FUNCTION itau2date(itau, date0, deltat)
          !
          IMPLICIT NONE
          !
          !     This function transforms itau into a date. The date whith which
          !     the time axis is going to be labeled
          !
          !     INPUT 
          !          itau   : current time step
          !          date0  : Date at which itau was equal to 0
          !          deltat : time step between itau s
          !
          !     OUTPUT
          !          itau2date : Date for the given itau
          !
          !
          !
          INTEGER  :: itau
          REAL(8)     :: date0, deltat
          !
          itau2date = float(itau)*deltat/un_jour + date0
          !
        END FUNCTION itau2date
        !
        !=======================================================================
        !
        !
        SUBROUTINE itau2ymds(itau, deltat, year, month, date, sec)
          !
          IMPLICIT NONE
          !
          !     This function transforms itau into a date. The date whith which
          !     the time axis is going to be labeled
          !
          !     INPUT 
          !          itau   : current time step
          !          deltat : time step between itau s
          !
          !     OUTPUT
          !          year : year
          !          month : month
          !          date : date
          !          sec  : seconds since midnight
          !
          ! INPUT
          !
          INTEGER, INTENT(IN) :: itau
          REAL(8), INTENT(IN)    :: deltat
          !
          ! OUTPUT
          !
          INTEGER, INTENT(OUT) :: year, month, date
          REAL(8), INTENT(OUT)    :: sec
          !
          ! Local
          !
          INTEGER :: julian_day
          REAL(8)    :: julian_sec
          !
          julian_day = start_day
          julian_sec = start_sec + float(itau)*deltat
          !
        CALL ju2ymds_internal(julian_day, julian_sec, year, month,
     &   date, sec)
          !
        END SUBROUTINE itau2ymds
        !
        !=======================================================================
        !
        !
        REAL(8) FUNCTION dtchdate(itau, date0, old_dt, new_dt)
          !
          IMPLICIT NONE
          !
          !     This function changes the date so that the simulation can 
          !     continue with the same itau but a different dt.
          !
          !
          !     INPUT 
          !          itau   : current time step
          !          date0  : Date at which itau was equal to 0
          !          old_dt : Old time step between itaus
          !          new_dt : New time step between itaus
          !
          !     OUTPUT
          !          dtchdate : Date for the given itau
          !
          !
          !
          INTEGER, INTENT(IN) :: itau
          REAL(8), INTENT(IN)    :: date0, old_dt, new_dt
          !
          REAL(8) :: rtime
          !
          rtime = itau2date(itau, date0, old_dt)
          dtchdate = rtime - float(itau)*new_dt/un_jour
          !
        END FUNCTION dtchdate
        !
        !=======================================================================
        !
      
        SUBROUTINE isittime(itau, date0, dt, freq, last_action, last_check, do_action)
          !
          !
          !     This function checks the time as come for a given action. This
          !     is computed from the current time-step(itau). Thus we need to have the
          !     time delta (dt), the frequency of the action (freq) and the last time 
          !     it was done (last_action in units of itau). In order to
          !     extrapolate when will be the next check we need the time step of the 
          !     last call (last_check).
          !
          !     The test is done on the following condition : The distance from the current time to the 
          !     time for the next action is smaller than the one from the next expected check to the next action. 
          !     When the test is done on the time steps simplifactions make it more difficult to read in the code.
          !     For the real(8) time case it is easier to understand !
          !
          IMPLICIT NONE
          !
          !     INPUT
          !
          INTEGER,INTENT(IN) :: itau
          REAL(8),INTENT(IN)    :: dt, freq
          INTEGER,INTENT(IN) :: last_action, last_check
          REAL(8),INTENT(IN)    :: date0
          !
          !     OUTPUT
          !
          LOGICAL,INTENT(OUT)  :: do_action
          !
          !     LOCAL
          !
          REAL(8) :: dt_action, dt_check
          REAL(8) :: date_last_act, date_next_check, date_next_act, date_now, date_mp1, date_mpf
          INTEGER :: year, month, monthp1, day, next_check_itau, next_act_itau
          INTEGER :: yearp, dayp
          REAL(8) :: sec, secp
          LOGICAL :: check
          !
          !
          check = .FALSE.
          !
          IF (check) WRITE(*,*) "isittime 1.0 ", itau, date0, dt, freq,
     &      last_action, last_check
          !
          IF (last_check >= 0) THEN
              dt_action = (itau - last_action)*dt
              dt_check = (itau - last_check)*dt
              next_check_itau = itau + (itau - last_check)
              !
              ! We are dealing with frequencies in seconds and thus operation can be done
              ! on the time steps.
              !
              IF ( freq .GT. 0) THEN
                  IF ( ABS(dt_action - freq) .LE. ABS(dt_action +
     &             dt_check - freq) ) THEN
                      do_action = .TRUE.
                  ELSE
                      do_action = .FALSE.
                  ENDIF
                  !
                  ! Here we deal with frequencies in month and work on julian days.
                  !
              ELSE
                  date_now = itau2date(itau, date0, dt)
                  date_last_act = itau2date(last_action, date0, dt)
                  CALL ju2ymds(date_last_act, year, month, day, sec)
                  monthp1 = month - freq
                  yearp = year
                  !
                  ! Here we compute what logically should be the next month
                  !            
                  IF ( month .GE. 13 ) THEN
                    yearp = year + 1
                    monthp1 = monthp1 - 12
                  ENDIF
                  CALL ymds2ju(year, monthp1, day, sec, date_mpf)
                  !
                  ! But it could be that because of a shorter month or a bad starting date
                  ! that we end up further than we should be. Thus we compute the first day
                  ! of the next month. We can not be beyond this date and if we are close 
                  ! then we will take it as it is better.
                  ! 
                  monthp1 = month + abs(freq)
                  yearp=year
                  IF ( monthp1 .GE. 13 ) THEN
                     yearp = year  + 1
                     monthp1 = monthp1 -12
                  ENDIF
                  dayp = 1
                  secp = 0.0
                  CALL ymds2ju(yearp, monthp1, dayp, secp, date_mp1)
                  !
                  ! If date_mp1 is smaller than date_mpf or only less than 4 days
                  ! larger then we take it. This needed to ensure that short month
                  ! like February do not mess up the thing !
                  !
                  IF ( date_mp1 - date_mpf .LT. 4. ) THEN
                     date_next_act = date_mp1
                  ELSE
                     date_next_act = date_mpf
                  ENDIF
              date_next_check  = itau2date(next_check_itau, date0, dt)
                  ! 
      ! Transform the dates into time-steps for the needed precisions.
                  !
                   next_act_itau = last_action +
     &              INT((date_next_act - date_last_act)*(un_jour/dt))
                  !
                  IF ( ABS(itau - next_act_itau)  .LE.
     &             ABS( next_check_itau - next_act_itau) ) THEN
                      do_action = .TRUE.
                      IF ( check ) THEN
       WRITE(*,*) 'ACT-TIME : itau, next_act_itau, next_check_itau : ',
     &               itau, next_act_itau, next_check_itau
                         CALL ju2ymds(date_now, year, month, day, sec)
       WRITE(*,*) 'ACT-TIME : y, m, d, s : ', year, month, day, sec
       WRITE(*,*) 'ACT-TIME : date_mp1, date_mpf : ', date_mp1, date_mpf
                      ENDIF
                  ELSE
                      do_action = .FALSE.
                  ENDIF
              ENDIF
              !
              IF (check) WRITE(*,*) "isittime 2.0 ",
     &         date_next_check, date_next_act, ABS(dt_action - freq),
     & ABS(dt_action + dt_check - freq), dt_action,
     & dt_check, next_check_itau, do_action
              !
          ELSE 
              do_action=.FALSE.
          ENDIF
          !
        END SUBROUTINE isittime
        ! 
        !===========================================================
        !
        SUBROUTINE ioconf_calendar(str)
  !
  !     This routine allows to configure the calendar to be used.
  !     This operation is only allowed once and the first call to
  !     ymds2ju or ju2ymsd will lock the current configuration.
  !     the argument to ioconf_calendar can be any of the following :
  !     - gregorian : This is the gregorian calendar (default here)
  !     - noleap    : A calendar without leap years = 365 days
  !     - xxxd      : A calendar of xxx days (has to be a modulo of 12)
  !                   with 12 month of equal length
          !
          !
          IMPLICIT NONE
          !
          !     INPUT
          !
          CHARACTER(LEN=*), INTENT(IN) :: str
          !
          !     LOCAL
          !
          INTEGER :: leng, ipos
          CHARACTER(LEN=10) :: str10
          !
          !     1.0 Clean up the sring !
          !    
          CALL strlowercase(str)
          !      
          IF (.NOT. lock_unan) THEN
             !
             lock_unan=.TRUE.
             !
             SELECT CASE(str)
                !
             CASE('gregorian')
                un_an = 365.2425
             CASE('noleap')
                un_an = 365.0
             CASE DEFAULT
                ipos = INDEX(str,'d')
                IF ( ipos .EQ. 4) THEN
                   READ(str(1:3),'(I3)') leng
                   IF ( MOD(leng,12) .EQ. 0 .AND. leng .GT. 1) THEN
                      un_an = leng
                   ELSE
                      CALL histerr(3,'ioconf_calendar',
     &             'The length of the year as to be a modulo of 12',
     &'so that it can be divided into 12 month of equal length', str)
                   ENDIF
                ELSE
      CALL histerr(3,'ioconf_calendar',
     & 'Unrecognized input, please ceck the man pages.',
     &             str, ' ')
                ENDIF
             END SELECT
          ELSE
             WRITE(str10,'(f10.4)') un_an
             CALL histerr(2,'ioconf_calendar',
     &'The calendar was already used or configured. You are',
     &'not allowed to change it again.The followingLengthOfYearIsUsed :'
     &,str10)
          ENDIF
          !
          RETURN
          !
        END SUBROUTINE ioconf_calendar
        !
        !=======================================================================
        !
        SUBROUTINE ioconf_startdate_simple(julian)
          !
          IMPLICIT NONE
          !
          !  INPUT
          !
          REAL(8), INTENT(IN) :: julian
          !
          ! LOCAL
          !
          INTEGER :: julian_day
          REAL(8)    :: julian_sec
      
          julian_day = INT(julian)
          julian_sec = (julian - julian_day)*un_jour
      
          CALL ioconf_startdate_internal(julian_day, julian_sec)
      
        END SUBROUTINE ioconf_startdate_simple
        !
        !=======================================================================
        !
        SUBROUTINE ioconf_startdate_ymds(year, month, day, sec)
          !
          IMPLICIT NONE
          !
          !  INPUT
          !
          INTEGER, INTENT(IN) :: year, month, day
          REAL(8), INTENT(IN)    :: sec
          !
          !  LOCAL
          !
          INTEGER :: julian_day
          REAL(8)    :: julian_sec
      
          CALL ymds2ju_internal(year, month, day, sec, julian_day,
     &     julian_sec)
      
          CALL ioconf_startdate_internal(julian_day, julian_sec)
      
        END SUBROUTINE ioconf_startdate_ymds
        !
        !=======================================================================
        !
        SUBROUTINE ioconf_startdate_internal(julian_day, julian_sec)
          !
          ! This subroutine allows to set the startdate for later
          ! use. It allows the applications to access the date directly from
          ! the timestep. In order to avoid any problems the start date will
          ! be locked and can not be changed once set.
          !
          !
          INTEGER, INTENT(IN)  :: julian_day
          REAL(8), INTENT(IN)     :: julian_sec
          !
          CHARACTER(len=70) :: str70a, str70b
          !
          IF ( .NOT. lock_startdate ) THEN
             !
             lock_startdate = .TRUE.
             start_day = julian_day
             start_sec = julian_sec
             !
          ELSE
             !
      WRITE(str70a,'("The date you tried to set : ",f10.4)')
     & julian_day, julian_sec/un_jour
      WRITE(str70b,'("The dateWhichWasAlreadyinTheCalendar: ",f10.4)')
     &   start_day + start_sec/un_jour
             CALL histerr(2,'ioconf_startdate The start date has ',
     &       'already been set and you tried to change it', 
     &  str70a, str70b)
             !
          ENDIF
          !
          lock_startdate = .TRUE.
          !
        END SUBROUTINE ioconf_startdate_internal
        !
        !=======================================================================
        !
        SUBROUTINE ioget_calendar_str(str)
          !
          !     This subroutine returns the name of the calendar used here. Three
          !     options exist :
          !     - gregorian : This is the gregorian calendar (default here)
          !     - noleap    : A calendar without leap years = 365 days
          !     - xxxd      : A calendar of xxx days (has to be a modulo of 12)
          !                   with 12 month of equal length
          !
          IMPLICIT NONE
          !
          !     INPUT
          !
          CHARACTER(LEN=*), INTENT(OUT) :: str
          !
          !     This routine will lock the calendar. You do not want it to change after your inquiry.
          !
          lock_unan = .TRUE. 
          !
          IF ( un_an .GT. 365.0 .AND. un_an .LT. 366.0) THEN
              str = 'gregorian'
          ELSE IF ( ABS(un_an - 365.0) .LE. EPSILON(un_an) ) THEN
              str = 'noleap'
          ELSE
              WRITE(str,'(I3.3,"d")') INT(un_an)
          ENDIF
          !
          RETURN
          !
        END SUBROUTINE ioget_calendar_str
        !
        !=======================================================================
        !
        SUBROUTINE ioget_calendar_real1(long_an)
          !
          !     This subroutine returns the name of the calendar used here. Three
          !     options exist :
          !     - gregorian : This is the gregorian calendar (default here)
          !     - noleap    : A calendar without leap years = 365 days
          !     - xxxd      : A calendar of xxx days (has to be a modulo of 12)
          !                   with 12 month of equal length
          !
          IMPLICIT NONE
          !
          !     INPUT
          !
          REAL(8), INTENT(OUT) :: long_an
          !
          !     This routine will lock the calendar. You do not want it to change after your inquiry.
          !
          lock_unan = .TRUE. 
          !
          long_an = un_an 
          !
          RETURN
          !
        END SUBROUTINE ioget_calendar_real1
        !
        !
        !=======================================================================
        !
        SUBROUTINE ioget_calendar_real2(long_an, long_jour)
          !
          !     This subroutine returns the name of the calendar used here. Three
          !     options exist :
          !     - gregorian : This is the gregorian calendar (default here)
          !     - noleap    : A calendar without leap years = 365 days
          !     - xxxd      : A calendar of xxx days (has to be a modulo of 12)
          !                   with 12 month of equal length
          !
          IMPLICIT NONE
          !
          !     INPUT
          !
          REAL(8), INTENT(OUT) :: long_an, long_jour
          !
          !     This routine will lock the calendar. You do not want it to change after your inquiry.
          !
          lock_unan = .TRUE. 
          !
          long_an = un_an 
          long_jour = un_jour
          !
          RETURN
          !
        END SUBROUTINE ioget_calendar_real2
        !
        !
        !=======================================================================
        !
        SUBROUTINE ioget_timestamp(string)
          !
          IMPLICIT NONE
          !
          !     INPUT
          !
          CHARACTER(LEN=30), INTENT(OUT) :: string
          !
          !     LOCAL
          !
          INTEGER :: date_time(8)
          CHARACTER(LEN=10) :: bigben(3)
          !
          IF ( INDEX(time_stamp,'XXXXXX') .GT. 0) THEN
              !
        CALL DATE_AND_TIME(bigben(1), bigben(2), bigben(3), date_time)
              !
              WRITE(time_stamp,"(I4.4,'-',A3,'-',I2.2,
     &         ' ',I2.2,':',I2.2,':',I2.2,' GMT',a5)") date_time(1), 
     &         cal(date_time(2)),date_time(3),date_time(5),
     & date_time(6), date_time(7), bigben(3)
              !
          ENDIF
          !
          string = time_stamp
          !
          RETURN
          !
        END SUBROUTINE ioget_timestamp
        !
        !
        !=======================================================================
        !
        SUBROUTINE time_add(year_s, month_s, day_s, sec_s, sec_increment, year_e, month_e, day_e, sec_e)
          !
          IMPLICIT NONE
          !
          ! This subroutine allows to increment a date by a number of seconds.
          !
          ! INPUT
          !
          INTEGER, INTENT(IN) :: year_s, month_s, day_s
          REAL(8), INTENT(IN)    :: sec_s
          REAL(8), INTENT(IN)    :: sec_increment            !! Time in seconds to be added to the date
          !
          ! OUTUT
          !
          INTEGER, INTENT(OUT) :: year_e, month_e, day_e
          REAL(8), INTENT(OUT)    :: sec_e
          !
          ! LOCAL
          !
          INTEGER :: julian_day
          REAL(8)    :: julian_sec
          !
          CALL ymds2ju_internal(year_s, month_s, day_s, sec_s,
     &     julian_day, julian_sec)
          !
          julian_sec = julian_sec + sec_increment
          !
          CALL ju2ymds_internal(julian_day, julian_sec, year_e,
     &     month_e, day_e, sec_e)
          !
        END SUBROUTINE time_add
        !
        !
        !=======================================================================
        !
        SUBROUTINE time_diff(year_s, month_s, day_s, sec_s, year_e, month_e, day_e, sec_e, sec_diff)
          !
          IMPLICIT NONE
          !
          ! This subroutine allows to dtermine the number of seconds between two dates.
          !
          ! INPUT
          !
          INTEGER, INTENT(IN) :: year_s, month_s, day_s
          REAL(8), INTENT(IN)    :: sec_s
          INTEGER, INTENT(IN) :: year_e, month_e, day_e
          REAL(8), INTENT(IN)    :: sec_e
          !
          ! OUTPUT
          !
          REAL(8), INTENT(OUT)    :: sec_diff           !! Time in seconds between the two dates
          !
          ! LOCAL
          !
          INTEGER :: julian_day_s, julian_day_e, day_diff
          REAL(8)    :: julian_sec_s, julian_sec_e
          !
          CALL ymds2ju_internal(year_s, month_s, day_s, sec_s,
     &     julian_day_s, julian_sec_s)
          CALL ymds2ju_internal(year_e, month_e, day_e, sec_e,
     &     julian_day_e, julian_sec_e)
          !
          day_diff = julian_day_e - julian_day_s
          sec_diff = julian_sec_e - julian_sec_s
          !
          sec_diff = sec_diff + day_diff*un_jour
          !
        END SUBROUTINE time_diff
        !
       !=======================================================================

      SUBROUTINE read_date_string(datestring,year, month, day, seconds)
      IMPLICIT NONE

      CHARACTER(LEN=17), INTENT(IN ) :: datestring
      integer,      INTENT(OUT) :: year, month, day
      real(8),      INTENT(OUT) ::  seconds
      integer hr, minutes,sec

      read(datestring,'(I4,I2,I2)') year, month, day
      read(datestring(10:11),'(I2)') hr
      read(datestring(13:14),'(I2)') minutes
      read(datestring(16:17),'(I2)') sec

      seconds = hr*3600.0 + minutes*60. + sec;
      END SUBROUTINE read_date_string

       !=======================================================================

      SUBROUTINE write_date_string(datestring,year,month,day, seconds)
      IMPLICIT NONE
!      INPUTS, not declared because of calling of ju2ymds
      INTEGER            :: year, month, day
      REAL(8)            :: seconds
      CHARACTER(LEN=17), INTENT(OUT) :: datestring

      INTEGER           :: hr, minutes,sec
      CHARACTER(LEN=8)  :: thedate, ora
!        ! ********* scamuffo per evitare di avere 20100102-24:00:00 ***
!        ! ********* e avere piuttosto             20100103-00:00:00 ***
       REAL(8) endofday, julian
       endofday=3600.*24 - 1.0
       IF (seconds.gt.endofday) THEN
          CALL ymds2ju(year,  month, day, seconds,julian)
          CALL ju2ymds(julian, year, month,day,seconds)
       ENDIF
       WRITE(thedate,'(I4,I2.2,I2.2)') year, month, day

       hr      = INT(seconds/3600.)
       minutes = INT((seconds - hr*3600)/60.)
       sec     = INT(seconds -hr*3600. -minutes*60.)

       WRITE(ora,"(I2.2,':',I2.2,':',I2.2)") hr, minutes, sec
       WRITE(datestring,"(A8, '-', A8 )")    thedate, ora

      END SUBROUTINE write_date_string

       !=======================================================================

       REAL(8) FUNCTION datestring2sec(datestring)
       IMPLICIT NONE
       CHARACTER(LEN=17), INTENT(IN) :: datestring
       !REAL(8)          , INTENT(OUT):: seconds

       ! LOCAL
       INTEGER  :: year, month, day
       real(8)  :: sec, julian

       call read_date_string(datestring, year, month, day, sec)
       CALL ymds2ju(year,month,day,sec,julian)
       datestring2sec = julian*86400.0

       END FUNCTION datestring2sec

!==============================================================================
      SUBROUTINE tau2julianday(TAU, deltaT,julian)

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: TAU
      REAL(8), INTENT(IN)  :: deltaT
      INTEGER, INTENT(OUT) :: julian

      ! LOCAL
      integer :: year, month,day
      real(8) :: sec, julianFrom1582, julianOnYear


         call itau2ymds(TAU,deltaT, year, month,day,sec)
         call ymds2ju(year, month, day, sec       , julianFrom1582)
         call ymds2ju(year,     1,  0,  real(0.0,8) , julianOnYear)

         julian = INT(julianFrom1582-julianOnYear)


      END SUBROUTINE tau2julianday


!=======================================================================

      REAL(8) FUNCTION photoperiod(jday, latit)
!     jd          = julian day
!     latit       = latitude in degrees
!     photoperiod = day length in hours


      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: jday
      REAL(8), INTENT(IN)    :: latit

      ! LOCAL
      REAL(8)    :: Pi
      REAL(8)    :: declin,refract,sundia,f0,x1


      Pi = 3.1415927

!      declin = 180./Pi * pho_delta(jday)
      declin  = p_delta(jday)
      refract = Pi * 17./30./180.
      sundia  = Pi * 4./15./180.
      f0 = 0

      x1 = tan(declin + 0.5 * sundia + refract) * tan(Pi * latit / 180.)

      if (x1 .GE. 1) then
         f0 = Pi
      endif

      if ( (x1 .GT. -1) .AND. (x1 .LT. 1) ) then
         f0 = Pi / 2. + atan(x1 / sqrt(1. - x1 * x1))
      endif

      photoperiod  = f0/Pi*24.

      END FUNCTION photoperiod

!=======================================================================
      REAL(8) FUNCTION p_delta(jday)

      implicit none

      INTEGER, INTENT(IN)  :: jday
      ! LOCAL
      REAL(8)  :: sdelta


      sdelta  = sin(p_oblique(jday)) * sin(p_lambda(jday))
      p_delta = atan(sdelta / sqrt(1. - sdelta * sdelta))

      END FUNCTION p_delta

!=======================================================================
      REAL(8) FUNCTION p_oblique(jday)

      implicit none

      INTEGER, INTENT(IN)  :: jday
      !LOCAL

      REAL(8)  :: rads

      rads = 3.1415927 / 180.
      p_oblique = 23.439 * rads - 4.e-07 * rads * REAL(jday,8)

      END FUNCTION p_oblique
!=======================================================================

      REAL(8) FUNCTION p_lambda(jday)
      implicit none

      INTEGER, INTENT(IN)  :: jday
      ! LOCAL
      REAL(8)  :: rads, L, g

      rads = 3.1415927 / 180.

      L        = p_range(280.461 * rads + 0.9856474 * rads * jday)
      g        = p_range(357.528 * rads + 0.9856003 * rads * jday)
      p_lambda = p_range(L + 1.915 * rads * sin(g) + 0.02 * rads * sin(2 * g))

      end function p_lambda
!=======================================================================
      REAL(8) FUNCTION p_range(x)
      implicit none


      REAL(8), INTENT(IN)  :: x
      ! LOCAL
      REAL(8)  :: tpi,b,a

      tpi = 2 * 3.1415927

      b = x / tpi
      a = tpi * (b - floor(b))

      if (a .LT. 0) then
         a = tpi + a
      endif

      p_range = a

      END FUNCTION p_range
!=======================================================================
      END MODULE calendar



