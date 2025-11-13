MODULE calendar
!-
!$Id: calendar.f90 2459 2010-12-07 11:17:48Z smasson $
!-
! This software is governed by the CeCILL license
! See IOIPSL/IOIPSL_License_CeCILL.txt
!---------------------------------------------------------------------
!- This is the calendar which going to be used to do all
!- calculations on time. Three types of calendars are possible :
!-
!-  - gregorian :
!-      The normal calendar. The time origin for the
!-      julian day in this case is 24 Nov -4713
!-      (other names : 'standard','proleptic_gregorian')
!-  - noleap :
!-      A 365 day year without leap years.
!-      The origin for the julian days is in this case 1 Jan 0
!-      (other names : '365_day','365d')
!-  - all_leap :
!-      A 366 day year with leap years.
!-      The origin for the julian days is in this case ????
!-      (other names : '366_day','366d'
!-  - julian :
!-      same as gregorian, but with all leap century years
!-  - xxxd :
!-      Year of xxx days with month of equal length.
!-      The origin for the julian days is then also 1 Jan 0
!-
!- As one can see it is difficult to go from one calendar to the other.
!- All operations involving julian days will be wrong.
!- This calendar will lock as soon as possible
!- the length of the year and forbid any further modification.
!-
!- For the non leap-year calendar the method is still brute force.
!- We need to find an Integer series which takes care of the length
!- of the various month. (Jan)
!-
!-   one_day  : one day in seconds
!-   one_year : one year in days
!---------------------------------------------------------------------
  USE stringop,ONLY  : strlowercase
  USE errioipsl,ONLY : ipslerr
!-
  PRIVATE
  PUBLIC :: ymds2ju,ju2ymds,tlen2itau,isittime,ioconf_calendar, &
 &          ioget_calendar,ioget_mon_len,ioget_year_len,itau2date, &
 &          ioget_timestamp,ioconf_startdate,itau2ymds, &
 &          time_diff,time_add,lock_calendar
!-
  INTERFACE ioget_calendar
    MODULE PROCEDURE &
 &    ioget_calendar_real1,ioget_calendar_real2,ioget_calendar_str
  END INTERFACE
!-
  INTERFACE ioconf_startdate
     MODULE PROCEDURE &
 &    ioconf_startdate_simple,ioconf_startdate_internal, &
 &    ioconf_startdate_ymds
  END INTERFACE
!-
  REAL,PARAMETER :: one_day = 86400.0
  LOGICAL,SAVE :: lock_startdate = .FALSE.
!-
  CHARACTER(LEN=30),SAVE :: time_stamp='XXXXXXXXXXXXXXXX'
!-
!- Description of calendar
!-
  CHARACTER(LEN=20),SAVE :: calendar_used="gregorian"
  LOGICAL,SAVE :: lock_one_year = .FALSE.
  REAL,SAVE :: one_year = 365.2425
  INTEGER,SAVE :: mon_len(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
!-
  CHARACTER(LEN=3),PARAMETER :: &
 &  cal(12) = (/'JAN','FEB','MAR','APR','MAY','JUN', &
 &              'JUL','AUG','SEP','OCT','NOV','DEC'/)
!-
  REAL,SAVE :: start_day,start_sec
!-
CONTAINS
!-
!===
!-
SUBROUTINE lock_calendar (new_status,old_status)
!!--------------------------------------------------------------------
!! The "lock_calendar" routine
!! allows to lock or unlock the calendar,
!! and to know the current status of the calendar.
!! Be careful !
!!
!! SUBROUTINE lock_calendar (new_status,old_status)
!!
!! Optional INPUT argument
!!
!! (L) new_status : new status of the calendar
!!
!! Optional OUTPUT argument
!!
!! (L) old_status : current status of the calendar
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  LOGICAL,OPTIONAL,INTENT(IN)  :: new_status
  LOGICAL,OPTIONAL,INTENT(OUT) :: old_status
!---------------------------------------------------------------------
  IF (PRESENT(old_status)) THEN
    old_status = lock_one_year
  ENDIF
  IF (PRESENT(new_status)) THEN
    lock_one_year = new_status
  ENDIF
!---------------------------
END SUBROUTINE lock_calendar
!-
!===
!-
SUBROUTINE ymds2ju (year,month,day,sec,julian)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: year,month,day
  REAL,INTENT(IN)    :: sec
!-
  REAL,INTENT(OUT) :: julian
!-
  INTEGER :: julian_day
  REAL    :: julian_sec
!---------------------------------------------------------------------
  CALL ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)
!-
  julian = julian_day+julian_sec/one_day
!---------------------
END SUBROUTINE ymds2ju
!-
!===
!-
SUBROUTINE ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)
!---------------------------------------------------------------------
!- Converts year, month, day and seconds into a julian day
!-
!- In 1968 in a letter to the editor of Communications of the ACM
!- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
!- and Thomas C. Van Flandern presented such an algorithm.
!-
!- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
!-
!- In the case of the Gregorian calendar we have chosen to use
!- the Lilian day numbers. This is the day counter which starts
!- on the 15th October 1582.
!- This is the day at which Pope Gregory XIII introduced the
!- Gregorian calendar.
!- Compared to the true Julian calendar, which starts some
!- 7980 years ago, the Lilian days are smaler and are dealt with
!- easily on 32 bit machines. With the true Julian days you can only
!- the fraction of the day in the real part to a precision of
!- a 1/4 of a day with 32 bits.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: year,month,day
  REAL,INTENT(IN)    :: sec
!-
  INTEGER,INTENT(OUT) :: julian_day
  REAL,INTENT(OUT)    :: julian_sec
!-
  INTEGER :: jd,m,y,d,ml
!---------------------------------------------------------------------
  lock_one_year = .TRUE.
!-
  m = month
  y = year
  d = day
!-
!- We deduce the calendar from the length of the year as it
!- is faster than an INDEX on the calendar variable.
!-
  IF ( (one_year > 365.0).AND.(one_year < 366.0) ) THEN
!-- "Gregorian"
    jd = (1461*(y+4800+INT((m-14)/12)))/4 &
 &      +(367*(m-2-12*(INT((m-14)/12))))/12 &
 &      -(3*((y+4900+INT((m-14)/12))/100))/4 &
 &      +d-32075
    jd = jd-2299160
  ELSE IF (    (ABS(one_year-365.0) <= EPSILON(one_year))  &
 &         .OR.(ABS(one_year-366.0) <= EPSILON(one_year)) ) THEN
!-- "No leap" or "All leap"
    ml = SUM(mon_len(1:m-1))
    jd = y*NINT(one_year)+ml+(d-1)
  ELSE
!-- Calendar with regular month
    ml = NINT(one_year/12.)
    jd = y*NINT(one_year)+(m-1)*ml+(d-1)
  ENDIF
!-
  julian_day = jd
  julian_sec = sec
!------------------------------
END SUBROUTINE ymds2ju_internal
!-
!===
!-
SUBROUTINE ju2ymds (julian,year,month,day,sec)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  REAL,INTENT(IN) :: julian
!-
  INTEGER,INTENT(OUT) :: year,month,day
  REAL,INTENT(OUT)    :: sec
!-
  INTEGER :: julian_day
  REAL    :: julian_sec
!---------------------------------------------------------------------
  julian_day = INT(julian)
  julian_sec = (julian-julian_day)*one_day
!-
  CALL ju2ymds_internal(julian_day,julian_sec,year,month,day,sec)
!---------------------
END SUBROUTINE ju2ymds
!-
!===
!-
SUBROUTINE ju2ymds_internal (julian_day,julian_sec,year,month,day,sec)
!---------------------------------------------------------------------
!- This subroutine computes from the julian day the year,
!- month, day and seconds
!-
!- In 1968 in a letter to the editor of Communications of the ACM
!- (CACM, volume 11, number 10, October 1968, p.657) Henry F. Fliegel
!- and Thomas C. Van Flandern presented such an algorithm.
!-
!- See also : http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm
!-
!- In the case of the Gregorian calendar we have chosen to use
!- the Lilian day numbers. This is the day counter which starts
!- on the 15th October 1582. This is the day at which Pope
!- Gregory XIII introduced the Gregorian calendar.
!- Compared to the true Julian calendar, which starts some 7980
!- years ago, the Lilian days are smaler and are dealt with easily
!- on 32 bit machines. With the true Julian days you can only the
!- fraction of the day in the real part to a precision of a 1/4 of
!- a day with 32 bits.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: julian_day
  REAL,INTENT(IN)    :: julian_sec
!-
  INTEGER,INTENT(OUT) :: year,month,day
  REAL,INTENT(OUT)    :: sec
!-
  INTEGER :: l,n,i,jd,j,d,m,y,ml
  INTEGER :: add_day
  REAL :: eps_day
!---------------------------------------------------------------------
  eps_day = SPACING(one_day)
  lock_one_year = .TRUE.
!-
  jd = julian_day
  sec = julian_sec
  IF (sec > (one_day-eps_day)) THEN
    add_day = INT(sec/one_day)
    sec = sec-add_day*one_day
    jd = jd+add_day
  ENDIF
  IF (sec < -eps_day) THEN
    sec = sec+one_day
    jd = jd-1
  ENDIF
!-
  IF ( (one_year > 365.0).AND.(one_year < 366.0) ) THEN
!-- Gregorian
    jd = jd+2299160
!-
    l = jd+68569
    n = (4*l)/146097
    l = l-(146097*n+3)/4
    i = (4000*(l+1))/1461001
    l = l-(1461*i)/4+31
    j = (80*l)/2447
    d = l-(2447*j)/80
    l = j/11
    m = j+2-(12*l)
    y = 100*(n-49)+i+l
  ELSE IF (    (ABS(one_year-365.0) <= EPSILON(one_year)) &
 &         .OR.(ABS(one_year-366.0) <= EPSILON(one_year)) ) THEN
!-- No leap or All leap
    y = jd/NINT(one_year)
    l = jd-y*NINT(one_year)
    m = 1
    ml = 0
    DO WHILE (ml+mon_len(m) <= l)
      ml = ml+mon_len(m)
      m = m+1
    ENDDO
    d = l-ml+1
  ELSE
!-- others
    ml = NINT(one_year/12.)
    y = jd/NINT(one_year)
    l = jd-y*NINT(one_year)
    m = (l/ml)+1
    d = l-(m-1)*ml+1
  ENDIF
!-
  day = d
  month = m
  year = y
!------------------------------
END SUBROUTINE ju2ymds_internal
!-
!===
!-
SUBROUTINE tlen2itau (input_str,dt,date,itau)
!---------------------------------------------------------------------
!- This subroutine transforms a string containing a time length
!- into a number of time steps.
!- To do this operation the date (in julian days is needed as the
!- length of the month varies.
!- The following convention is used :
!-   n   : n time steps
!-   nS  : n seconds is transformed into itaus
!-   nH  : n hours
!-   nD  : n days
!-   nM  : n month
!-   nY  : n years
!- Combinations are also possible
!-   nYmD : nyears plus m days !
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: input_str
  REAL,INTENT(IN)             :: dt,date
!-
  INTEGER,INTENT(OUT)         :: itau
!-
  INTEGER           :: y_pos,m_pos,d_pos,h_pos,s_pos
  INTEGER           :: read_time
  CHARACTER(LEN=13) :: fmt
  CHARACTER(LEN=80) :: tmp_str
!-
  INTEGER :: year,month,day
  REAL    :: sec,date_new,dd,ss
!---------------------------------------------------------------------
  itau = 0
  CALL ju2ymds (date,year,month,day,sec)
!-
  y_pos = MAX(INDEX(input_str,'y'),INDEX(input_str,'Y'))
  m_pos = MAX(INDEX(input_str,'m'),INDEX(input_str,'M'))
  d_pos = MAX(INDEX(input_str,'d'),INDEX(input_str,'D'))
  h_pos = MAX(INDEX(input_str,'h'),INDEX(input_str,'H'))
  s_pos = MAX(INDEX(input_str,'s'),INDEX(input_str,'S'))
!-
  IF (MAX(y_pos,m_pos,d_pos,s_pos) > 0) THEN
    tmp_str = input_str
    DO WHILE ( MAX(y_pos,m_pos,d_pos,s_pos) > 0)
!---- WRITE(*,*) tmp_str
!---- WRITE(*,*) y_pos,m_pos,d_pos,s_pos
      IF (y_pos > 0) THEN
        WRITE(fmt,'("(I",I10.10,")")') y_pos-1
        READ(tmp_str(1:y_pos-1),fmt) read_time
        CALL ymds2ju (year+read_time,month,day,sec,date_new)
        dd = date_new-date
        ss = INT(dd)*one_day+dd-INT(dd)
        itau = itau+NINT(ss/dt)
        tmp_str = tmp_str(y_pos+1:LEN_TRIM(tmp_str))
      ELSE IF (m_pos > 0) THEN
        WRITE(fmt,'("(I",I10.10,")")') m_pos-1
        READ(tmp_str(1:m_pos-1),fmt) read_time
        CALL ymds2ju (year,month+read_time,day,sec,date_new)
        dd = date_new-date
        ss = INT(dd)*one_day+dd-INT(dd)
        itau = itau+NINT(ss/dt)
        tmp_str = tmp_str(m_pos+1:LEN_TRIM(tmp_str))
      ELSE IF (d_pos > 0) THEN
        WRITE(fmt,'("(I",I10.10,")")') d_pos-1
        READ(tmp_str(1:d_pos-1),fmt) read_time
        itau = itau+NINT(read_time*one_day/dt)
        tmp_str = tmp_str(d_pos+1:LEN_TRIM(tmp_str))
      ELSE IF (h_pos > 0) THEN
        WRITE(fmt,'("(I",I10.10,")")') h_pos-1
        READ(tmp_str(1:h_pos-1),fmt) read_time
        itau = itau+NINT(read_time*60.*60./dt)
        tmp_str = tmp_str(d_pos+1:LEN_TRIM(tmp_str))
      ELSE IF  (s_pos > 0) THEN
        WRITE(fmt,'("(I",I10.10,")")') s_pos-1
        READ(tmp_str(1:s_pos-1),fmt) read_time
        itau = itau+NINT(read_time/dt)
        tmp_str = tmp_str(s_pos+1:LEN_TRIM(tmp_str))
      ENDIF
!-
      y_pos = MAX(INDEX(tmp_str,'y'),INDEX(tmp_str,'Y'))
      m_pos = MAX(INDEX(tmp_str,'m'),INDEX(tmp_str,'M'))
      d_pos = MAX(INDEX(tmp_str,'d'),INDEX(tmp_str,'D'))
      h_pos = MAX(INDEX(tmp_str,'h'),INDEX(tmp_str,'H'))
      s_pos = MAX(INDEX(tmp_str,'s'),INDEX(tmp_str,'S'))
    ENDDO
  ELSE
    WRITE(fmt,'("(I",I10.10,")")') LEN_TRIM(input_str)
    READ(input_str(1:LEN_TRIM(input_str)),fmt) itau
  ENDIF
!-----------------------
END SUBROUTINE tlen2itau
!-
!===
!-
REAL FUNCTION itau2date (itau,date0,deltat)
!---------------------------------------------------------------------
!- This function transforms itau into a date. The date with which
!- the time axis is going to be labeled
!-
!- INPUT
!-   itau   : current time step
!-   date0  : Date at which itau was equal to 0
!-   deltat : time step between itau s
!-
!- OUTPUT
!-   itau2date : Date for the given itau
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER  :: itau
  REAL     :: date0,deltat
!---------------------------------------------------------------------
  itau2date = REAL(itau)*deltat/one_day+date0
!---------------------
END FUNCTION itau2date
!-
!===
!-
SUBROUTINE itau2ymds (itau,deltat,year,month,day,sec)
!---------------------------------------------------------------------
!- This subroutine transforms itau into a date. The date with which
!- the time axis is going to be labeled
!-
!- INPUT
!-   itau   : current time step
!-   deltat : time step between itau s
!-
!- OUTPUT
!-   year  : year
!-   month : month
!-   day   : day
!-   sec   : seconds since midnight
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: itau
  REAL,INTENT(IN)    :: deltat
!-
  INTEGER,INTENT(OUT) :: year,month,day
  REAL,INTENT(OUT)    :: sec
!-
  INTEGER :: julian_day
  REAL    :: julian_sec
!---------------------------------------------------------------------
  IF (.NOT.lock_startdate) THEN
    CALL ipslerr (2,'itau2ymds', &
 &   'You try to call this function, itau2ymds, but you didn''t', &
 &   ' call ioconf_startdate to initialize date0 in calendar.', &
 &   ' Please call ioconf_startdate before itau2ymds.')
  ENDIF
  julian_day = start_day
  julian_sec = start_sec+REAL(itau)*deltat
  CALL ju2ymds_internal (julian_day,julian_sec,year,month,day,sec)
!-----------------------
END SUBROUTINE itau2ymds
!-
!===
!-
REAL FUNCTION dtchdate (itau,date0,old_dt,new_dt)
!---------------------------------------------------------------------
!- This function changes the date so that the simulation can
!- continue with the same itau but a different dt.
!-
!- INPUT
!-   itau   : current time step
!-   date0  : Date at which itau was equal to 0
!-   old_dt : Old time step between itaus
!-   new_dt : New time step between itaus
!-
!- OUTPUT
!-   dtchdate : Date for the given itau
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: itau
  REAL,INTENT(IN)    :: date0,old_dt,new_dt
!-
  REAL :: rtime
!---------------------------------------------------------------------
  rtime = itau2date (itau,date0,old_dt)
  dtchdate = rtime-REAL(itau)*new_dt/one_day
!--------------------
END FUNCTION dtchdate
!-
!===
!-
SUBROUTINE isittime &
 &  (itau,date0,dt,freq,last_action,last_check,do_action)
!---------------------------------------------------------------------
!- This subroutine checks the time as come for a given action.
!- This is computed from the current time-step(itau).
!- Thus we need to have the time delta (dt), the frequency
!- of the action (freq) and the last time it was done
!- (last_action in units of itau).
!- In order to extrapolate when will be the next check we need
!- the time step of the last call (last_check).
!-
!- The test is done on the following condition :
!- the distance from the current time to the time for the next
!- action is smaller than the one from the next expected
!- check to the next action.
!- When the test is done on the time steps simplifications make
!- it more difficult to read in the code.
!- For the real time case it is easier to understand !
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: itau
  REAL,INTENT(IN)    :: dt,freq
  INTEGER,INTENT(IN) :: last_action,last_check
  REAL,INTENT(IN)    :: date0
!-
  LOGICAL,INTENT(OUT)  :: do_action
!-
  REAL :: dt_action,dt_check
  REAL :: date_last_act,date_next_check,date_next_act, &
 &        date_now,date_mp1,date_mpf
  INTEGER :: year,month,monthp1,day,next_check_itau,next_act_itau
  INTEGER :: yearp,dayp
  REAL :: sec,secp
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  IF (check) THEN
    WRITE(*,*) &
 &    "isittime 1.0 ",itau,date0,dt,freq,last_action,last_check
  ENDIF
!-
  IF (last_check >= 0) THEN
    dt_action = (itau-last_action)*dt
    dt_check = (itau-last_check)*dt
    next_check_itau = itau+(itau-last_check)
!-
!-- We are dealing with frequencies in seconds and thus operation
!-- can be done on the time steps.
!-
    IF (freq > 0) THEN
      IF (ABS(dt_action-freq) <= ABS(dt_action+dt_check-freq)) THEN
        do_action = .TRUE.
      ELSE
        do_action = .FALSE.
      ENDIF
!-
!---- Here we deal with frequencies in month and work on julian days.
!-
    ELSE
      date_now = itau2date (itau,date0,dt)
      date_last_act = itau2date (last_action,date0,dt)
      CALL ju2ymds (date_last_act,year,month,day,sec)
      monthp1 = month-freq
      yearp = year
!-
!---- Here we compute what logically should be the next month
!-
      DO WHILE (monthp1 >= 13)
        yearp = yearp+1
        monthp1 = monthp1-12
      END DO
      CALL ymds2ju (yearp,monthp1,day,sec,date_mpf)
!-
!---- But it could be that because of a shorter month or a bad
!---- starting date that we end up further than we should be.
!---- Thus we compute the first day of the next month.
!---- We can not be beyond this date and if we are close
!---- then we will take it as it is better.
!-
      monthp1 = month+ABS(freq)
      yearp=year
      DO WHILE (monthp1 >= 13)
        yearp = yearp+1
        monthp1 = monthp1-12
      END DO
      dayp = 1
      secp = 0.0
      CALL ymds2ju (yearp,monthp1,dayp,secp,date_mp1)
!-
!---- If date_mp1 is smaller than date_mpf or only less than 4 days
!---- larger then we take it. This needed to ensure that short month
!---- like February do not mess up the thing !
!-
      IF (date_mp1-date_mpf < 4.) THEN
        date_next_act = date_mp1
      ELSE
        date_next_act = date_mpf
      ENDIF
      date_next_check = itau2date (next_check_itau,date0,dt)
!-
!---- Transform the dates into time-steps for the needed precisions.
!-
      next_act_itau = &
 &      last_action+INT((date_next_act-date_last_act)*(one_day/dt))
!-----
      IF (   ABS(itau-next_act_itau) &
 &        <= ABS( next_check_itau-next_act_itau)) THEN
        do_action = .TRUE.
        IF (check) THEN
          WRITE(*,*) &
 &         'ACT-TIME : itau, next_act_itau, next_check_itau : ', &
 &         itau,next_act_itau,next_check_itau
          CALL ju2ymds (date_now,year,month,day,sec)
          WRITE(*,*) 'ACT-TIME : y, m, d, s : ',year,month,day,sec
          WRITE(*,*) &
 &         'ACT-TIME : date_mp1, date_mpf : ',date_mp1,date_mpf
        ENDIF
      ELSE
        do_action = .FALSE.
      ENDIF
    ENDIF
!-
    IF (check) THEN
      WRITE(*,*) "isittime 2.0 ", &
 &     date_next_check,date_next_act,ABS(dt_action-freq), &
 &     ABS(dt_action+dt_check-freq),dt_action,dt_check, &
 &     next_check_itau,do_action
    ENDIF
  ELSE
    do_action=.FALSE.
  ENDIF
!----------------------
END SUBROUTINE isittime
!-
!===
!-
SUBROUTINE ioconf_calendar (str)
!---------------------------------------------------------------------
!- This routine allows to configure the calendar to be used.
!- This operation is only allowed once and the first call to
!- ymds2ju or ju2ymsd will lock the current configuration.
!- the argument to ioconf_calendar can be any of the following :
!-  - gregorian : This is the gregorian calendar (default here)
!-  - noleap    : A calendar without leap years = 365 days
!-  - xxxd      : A calendar of xxx days (has to be a modulo of 12)
!-                with 12 month of equal length
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: str
!-
  INTEGER :: leng,ipos
  CHARACTER(LEN=20) :: str_w
!---------------------------------------------------------------------
!-
! Clean up the string !
!-
  str_w = str
  CALL strlowercase (str_w)
!-
  IF (.NOT.lock_one_year) THEN
!---
    lock_one_year=.TRUE.
!---
    SELECT CASE(TRIM(str_w))
    CASE('gregorian','standard','proleptic_gregorian')
      calendar_used = 'gregorian'
      one_year = 365.2425
      mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    CASE('noleap','365_day','365d')
      calendar_used = 'noleap'
      one_year = 365.0
      mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    CASE('all_leap','366_day','366d')
      calendar_used = 'all_leap'
      one_year = 366.0
      mon_len(:)=(/31,29,31,30,31,30,31,31,30,31,30,31/)
    CASE('360_day','360d')
      calendar_used = '360d'
      one_year = 360.0
      mon_len(:)=(/30,30,30,30,30,30,30,30,30,30,30,30/)
    CASE('julian')
      calendar_used = 'julian'
      one_year = 365.25
      mon_len(:)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    CASE DEFAULT
      ipos = INDEX(str_w,'d')
      IF (ipos == 4) THEN
        READ(str_w(1:3),'(I3)') leng
        IF ( (MOD(leng,12) == 0).AND.(leng > 1) ) THEN
          calendar_used = str_w
          one_year = leng
          mon_len(:) = leng/12
        ELSE
          CALL ipslerr (3,'ioconf_calendar', &
 &         'The length of the year as to be a modulo of 12', &
 &         'so that it can be divided into 12 month of equal length', &
 &         TRIM(str_w))
        ENDIF
      ELSE
        CALL ipslerr (3,'ioconf_calendar', &
 &       'Unrecognized input, please check the man pages.', &
 &       TRIM(str_w),' ')
      ENDIF
    END SELECT
  ELSE IF (TRIM(str_w) /= TRIM(calendar_used)) THEN
    WRITE(str_w,'(f10.4)') one_year
    CALL ipslerr (2,'ioconf_calendar', &
 &   'The calendar was already used or configured to : '// &
 &    TRIM(calendar_used)//'.', &
 &   'You are not allowed to change it to : '//TRIM(str)//'.', &
 &   'The following length of year is used : '//TRIM(ADJUSTL(str_w)))
  ENDIF
!-----------------------------
END SUBROUTINE ioconf_calendar
!-
!===
!-
SUBROUTINE ioconf_startdate_simple (julian)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  REAL,INTENT(IN) :: julian
!-
  INTEGER :: julian_day
  REAL    :: julian_sec
!---------------------------------------------------------------------
  julian_day = INT(julian)
  julian_sec = (julian-julian_day)*one_day
!-
  CALL ioconf_startdate_internal (julian_day,julian_sec)
!-------------------------------------
END SUBROUTINE ioconf_startdate_simple
!-
!===
!-
SUBROUTINE ioconf_startdate_ymds (year,month,day,sec)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: year,month,day
  REAL,INTENT(IN)    :: sec
!-
  INTEGER :: julian_day
  REAL    :: julian_sec
!---------------------------------------------------------------------
  CALL ymds2ju_internal (year,month,day,sec,julian_day,julian_sec)
!-
  CALL ioconf_startdate_internal (julian_day,julian_sec)
!-----------------------------------
END SUBROUTINE ioconf_startdate_ymds
!-
!===
!-
SUBROUTINE ioconf_startdate_internal (julian_day,julian_sec)
!---------------------------------------------------------------------
! This subroutine allows to set the startdate for later
! use. It allows the applications to access the date directly from
! the timestep. In order to avoid any problems the start date will
! be locked and can not be changed once set.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN)  :: julian_day
  REAL,INTENT(IN)     :: julian_sec
!-
  CHARACTER(len=70) :: str70a,str70b
!---------------------------------------------------------------------
  IF (.NOT.lock_startdate) THEN
    start_day = julian_day
    start_sec = julian_sec
    lock_startdate = .TRUE.
  ELSE
    WRITE(str70a,'("The date you tried to set : ",f10.4)') &
 &   julian_day,julian_sec/one_day
    WRITE(str70b, &
 &   '("The date which was already set in the calendar : ",f10.4)') &
 &   start_day+start_sec/one_day
    CALL ipslerr (2,'ioconf_startdate', &
 &   'The start date has already been set and you tried to change it', &
 &   str70a,str70b)
  ENDIF
!---------------------------------------
END SUBROUTINE ioconf_startdate_internal
!-
!===
!-
SUBROUTINE ioget_calendar_str (str)
!---------------------------------------------------------------------
!- This subroutine returns the name of the calendar used here.
!- Three options exist :
!-  - gregorian : This is the gregorian calendar (default here)
!-  - noleap    : A calendar without leap years = 365 days
!-  - xxxd      : A calendar of xxx days (has to be a modulo of 12)
!-                with 12 month of equal length
!-
!- This routine will lock the calendar.
!- You do not want it to change after your inquiry.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(OUT) :: str
!---------------------------------------------------------------------
  lock_one_year = .TRUE.
!-
  str = calendar_used
!--------------------------------
END SUBROUTINE ioget_calendar_str
!-
!===
!-
SUBROUTINE ioget_calendar_real1 (long_year)
!---------------------------------------------------------------------
!- This subroutine returns the name of the calendar used here.
!- Three options exist :
!-  - gregorian : This is the gregorian calendar (default here)
!-  - noleap    : A calendar without leap years = 365 days
!-  - xxxd      : A calendar of xxx days (has to be a modulo of 12)
!-                with 12 month of equal length
!-
!- This routine will lock the calendar.
!- You do not want it to change after your inquiry.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  REAL,INTENT(OUT) :: long_year
!---------------------------------------------------------------------
  long_year = one_year
  lock_one_year = .TRUE.
!----------------------------------
END SUBROUTINE ioget_calendar_real1
!-
!===
!-
SUBROUTINE ioget_calendar_real2 (long_year,long_day)
!---------------------------------------------------------------------
!- This subroutine returns the name of the calendar used here.
!- Three options exist :
!-  - gregorian : This is the gregorian calendar (default here)
!-  - noleap    : A calendar without leap years = 365 days
!-  - xxxd      : A calendar of xxx days (has to be a modulo of 12)
!-                with 12 month of equal length
!-
!- This routine will lock the calendar.
!- You do not want it to change after your inquiry.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  REAL,INTENT(OUT) :: long_year,long_day
!---------------------------------------------------------------------
  long_year = one_year
  long_day  = one_day
  lock_one_year = .TRUE.
!----------------------------------
END SUBROUTINE ioget_calendar_real2
!-
!===
!-
INTEGER FUNCTION ioget_mon_len (year,month)
!!--------------------------------------------------------------------
!! The "ioget_mon_len" function returns
!! the number of days in a "month" of a "year",
!! in the current calendar.
!!
!! INTEGER FUNCTION ioget_mon_len (year,month)
!!
!! INPUT
!!
!! (I) year  : year
!! (I) month : month in the year (1 --> 12)
!!
!! OUTPUT
!!
!! (I) ioget_mon_len : number of days in the month
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: year,month
!-
  INTEGER :: ml
!---------------------------------------------------------------------
  IF ( (month >= 1).AND.(month <= 12) ) THEN
    IF ( (one_year > 365.0).AND.(one_year < 366.0) ) THEN
!---- "Gregorian" or "Julian"
      ml = mon_len(month)
      IF (month == 2) THEN
        IF (ABS(one_year-365.2425) <= EPSILON(one_year) ) THEN
!-------- "Gregorian"
          IF (    ((MOD(year,4) == 0).AND.(MOD(year,100) /= 0)) &
              .OR.(MOD(year,400) == 0) ) THEN
            ml = ml+1
          ENDIF
        ELSE
!-------- "Julian"
          IF (MOD(year,4) == 0) THEN
            ml = ml+1
          ENDIF
        ENDIF
      ENDIF
      ioget_mon_len = ml
    ELSE
!---- "No leap" or "All leap" or "Calendar with regular month"
      ioget_mon_len = mon_len(month)
    ENDIF
  ELSE
    CALL ipslerr (3,'ioget_mon_len', &
 &    'The number of the month','must be between','1 and 12')
  ENDIF
!-------------------------
END FUNCTION ioget_mon_len
!-
!===
!-
INTEGER FUNCTION ioget_year_len (year)
!!--------------------------------------------------------------------
!! The "ioget_year_len" function returns
!! the number of days in "year", in the current calendar.
!!
!! INTEGER FUNCTION ioget_year_len (year)
!!
!! INPUT
!!
!! (I) year  : year
!!
!! OUTPUT
!!
!! (I) ioget_year_len : number of days in the year
!!--------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: year
!-
  INTEGER :: yl
!---------------------------------------------------------------------
  SELECT CASE(TRIM(calendar_used))
  CASE('gregorian')
    yl = 365
    IF (    ((MOD(year,4) == 0).AND.(MOD(year,100) /= 0)) &
        .OR.(MOD(year,400) == 0) ) THEN
      yl = yl+1
    ENDIF
  CASE('julian')
    yl = 365
    IF (MOD(year,4) == 0) THEN
      yl = yl+1
    ENDIF
  CASE DEFAULT
    yl = NINT(one_year)
  END SELECT
  ioget_year_len = yl
!--------------------------
END FUNCTION ioget_year_len
!-
!===
!-
SUBROUTINE ioget_timestamp (string)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=30),INTENT(OUT) :: string
!-
  INTEGER :: date_time(8)
  CHARACTER(LEN=10) :: bigben(3)
!---------------------------------------------------------------------
  IF (INDEX(time_stamp,'XXXXXX') > 0) THEN
    CALL DATE_AND_TIME (bigben(1),bigben(2),bigben(3),date_time)
!---
    WRITE(time_stamp, &
 &   "(I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2,' GMT',a5)") &
 &   date_time(1),cal(date_time(2)),date_time(3),date_time(5), &
 &   date_time(6),date_time(7),bigben(3)
  ENDIF
!-
  string = time_stamp
!-----------------------------
END SUBROUTINE ioget_timestamp
!-
!===
!-
SUBROUTINE time_add &
 &  (year_s,month_s,day_s,sec_s,sec_increment, &
 &   year_e,month_e,day_e,sec_e)
!---------------------------------------------------------------------
!- This subroutine allows to increment a date by a number of seconds.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: year_s,month_s,day_s
  REAL,INTENT(IN)    :: sec_s
!-
! Time in seconds to be added to the date
!-
  REAL,INTENT(IN)    :: sec_increment
!-
  INTEGER,INTENT(OUT) :: year_e,month_e,day_e
  REAL,INTENT(OUT)    :: sec_e
!-
  INTEGER :: julian_day
  REAL    :: julian_sec
!---------------------------------------------------------------------
  CALL ymds2ju_internal &
 &  (year_s,month_s,day_s,sec_s,julian_day,julian_sec)
!-
  julian_sec = julian_sec+sec_increment
!-
  CALL ju2ymds_internal &
 &  (julian_day,julian_sec,year_e,month_e,day_e,sec_e)
!----------------------
END SUBROUTINE time_add
!-
!===
!-
SUBROUTINE time_diff &
 &  (year_s,month_s,day_s,sec_s,year_e,month_e,day_e,sec_e,sec_diff)
!---------------------------------------------------------------------
!- This subroutine allows to determine the number of seconds
!- between two dates.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: year_s,month_s,day_s
  REAL,INTENT(IN)    :: sec_s
  INTEGER,INTENT(IN) :: year_e,month_e,day_e
  REAL,INTENT(IN)    :: sec_e
!-
! Time in seconds between the two dates
!-
  REAL,INTENT(OUT)    :: sec_diff
!-
  INTEGER :: julian_day_s,julian_day_e,day_diff
  REAL    :: julian_sec_s,julian_sec_e
!---------------------------------------------------------------------
  CALL ymds2ju_internal &
 &  (year_s,month_s,day_s,sec_s,julian_day_s,julian_sec_s)
  CALL ymds2ju_internal &
 &  (year_e,month_e,day_e,sec_e,julian_day_e,julian_sec_e)
!-
  day_diff = julian_day_e-julian_day_s
  sec_diff = julian_sec_e-julian_sec_s
!-
  sec_diff = sec_diff+day_diff*one_day
!-----------------------
END SUBROUTINE time_diff
!-
!===
!-
END MODULE calendar
