
  program test

  use esmf

  implicit none

  type(ESMF_Time) :: time1,time2,time3,time4,time5,time6,time7,time8
  type(ESMF_TimeInterval) :: timeint1,timeint2,timeint3,timeint4,timeint5
  type(ESMF_Calkind_Flag) :: calkindflag

  integer :: year,month,day,hour,min,sec,jday
  integer :: year1,month1,day1,hour1,min1,sec1,jday1
  integer :: year2,month2,day2,hour2,min2,sec2,jday2
  integer :: iyear,imonth,iday,ihour,imin,isec
  integer :: dyear,dmonth,dday,dhour,dmin,dsec
  integer :: icyear,icmonth,icday,ichour,icmin,icsec
  integer :: ical,i1,i2,delta
  integer :: errcnt, totcnt
  logical :: errfound
  character(len=8) :: dstr,calstr
  character(len=32) :: estr1,estr2

  INTEGER, PARAMETER :: mday(12)   &
                      = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER, PARAMETER :: mdayleap(12) &
                      = (/31,29,31,30,31,30,31,31,30,31,30,31/)

  character(len=*),parameter :: F01 = "(2x,a,1x,a6,i6,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2,1x,a8,i12,a8,i6)"
  character(len=*),parameter :: F02 = "(a,1x,a6,2x,i6,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,':',i2.2,1x,a8,i12)"
  character(len=*),parameter :: F03 = "(a,1x,i6,'-',i2.2,'-',i2.2,1x,a8,a8,i12)"

  call ESMF_Initialize()

  totcnt = 0
  errcnt = 0

  do icyear = 1,8
  do icmonth = 1,12
  do icday = 1,4
  do ichour = 2,2
  do icmin = 30,30
  do icsec = 10,10
  do ical = 1,2

     write(6,*) ' '
     write(estr1,'(i2.2,i2.2,i2.2,i2.2,i2.2,i2.2,i2.2)') icyear,icmonth,icday,ichour,icmin,icsec,ical

     if (icyear == 1) iyear = 0
     if (icyear == 2) iyear = 1
     if (icyear == 3) iyear = 1900
     if (icyear == 4) iyear = 1995
     if (icyear == 5) iyear = 1996
     if (icyear == 6) iyear = 2000
     if (icyear == 7) iyear = 9900
     if (icyear == 8) iyear = 9999

     imonth = icmonth

     if (icday == 1) iday = 1
     if (icday == 2) iday = 20
     if (icday == 3) iday = mday(imonth)-1
     if (icday == 4) iday = mday(imonth)

     ihour = ichour

     imin = icmin

     isec = icsec

     if (ical == 1) then
       calstr = 'noleap'
       calkindflag = ESMF_CALKIND_NOLEAP
     endif
     if (ical == 2) then
       calstr = 'gregor'
       calkindflag = ESMF_CALKIND_GREGORIAN
     endif

     write(6,F02) trim(estr1),'jd0 ',iyear,imonth,iday,ihour,imin,isec,trim(calstr)

     call ESMF_TimeSet(time1,yy=iyear,mm=imonth,dd=iday,h=ihour,m=imin,s=isec,calkindflag=calkindflag)

     call ESMF_TimeGet(time1,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
     call ESMF_TimeGet(time1,d=jday)
     write(6,F02) trim(estr1),'jd1 ',year,month,day,hour,min,sec,trim(calstr),jday
     call checkdate(year,month,day,hour,min,sec,trim(calstr))

     call ESMF_TimeSet(time2,d=jday,calkindflag=calkindflag)
     call ESMF_TimeGet(time2,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
     call ESMF_TimeGet(time2,d=jday)
     write(6,F02) trim(estr1),'jd2 ',year,month,day,hour,min,sec,trim(calstr),jday
     call checkdate(year,month,day,hour,min,sec,trim(calstr))

     if (year /= iyear .or. month /= imonth .or. day /= iday) then
         call wrf_error_fatal('ERROR: jday conversion')
     endif

     do i1 = 1,7
     do i2 = 1,4
        write(6,*) ' '
        write(estr2,'(a,i2.2,i2.2)') trim(estr1),i1,i2

        if (i2 == 1) delta = 1
        if (i2 == 2) delta = -1
        if (i2 == 3) delta = 150
        if (i2 == 4) delta = -150

        dyear = 0
        dmonth = 0
        dday = 0
        dhour =0
        dmin = 0
        dsec = 0

        if (i1 == 1) then
           dstr = 'year'
           dyear = delta
        endif
        if (i1 == 2) then
           dstr = 'month'
           dmonth = delta
        endif
        if (i1 == 3) then
           dstr = 'day'
           dday = delta
        endif
        if (i1 == 4) then
           dstr = 'hour'
           dhour = delta
        endif
        if (i1 == 5) then
           dstr = 'min'
           dmin = delta
        endif
        if (i1 == 6) then
           dstr = 'sec'
           dsec = delta
        endif
        if (i1 == 7) then
           dstr = 'all'
           dyear = delta
           dmonth = delta
           dday = delta
           dhour = delta
           dmin = delta
           dsec = delta
        endif

        call ESMF_TimeIntervalSet(timeint1,yy= dyear,mm= dmonth,d= dday,h= dhour,m= dmin,s= dsec)
        call ESMF_TimeIntervalSet(timeint2,yy=2*dyear,mm=2*dmonth,d=2*dday,h=2*dhour,m=2*dmin,s=2*dsec)
        call ESMF_TimeIntervalSet(timeint3,yy=-dyear,mm=-dmonth,d=-dday,h=-dhour,m=-dmin,s=-dsec)

       !time1 =                   ! zero
        time2 = time1 + timeint1  ! + delta
        timeint4 = time2 - time1  ! this should be same as timeint1 but only for time2-time1
        time3 = time2 - timeint4  ! zero
        time4 = time3 + timeint2  ! + 2*delta
        time5 = time4 - timeint1  ! + delta
        time6 = time5 + timeint3  ! zero
        time7 = time6 + timeint3  ! - delta
        time8 = time7 - timeint3  ! zero

        call ESMF_TimeGet(time1,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
        call ESMF_TimeGet(time1,d=jday)
        write(6,F01) trim(estr2),'ti1 ',year,month,day,hour,min,sec,trim(calstr),jday,trim(dstr),delta
        call checkdate(year,month,day,hour,min,sec,trim(calstr))

        call ESMF_TimeGet(time2,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
        call ESMF_TimeGet(time2,d=jday)
        write(6,F01) trim(estr2),'ti2 ',year,month,day,hour,min,sec,trim(calstr),jday,trim(dstr),delta
        call checkdate(year,month,day,hour,min,sec,trim(calstr))

        call ESMF_TimeGet(time3,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
        call ESMF_TimeGet(time3,d=jday)
        write(6,F01) trim(estr2),'ti3 ',year,month,day,hour,min,sec,trim(calstr),jday,trim(dstr),delta
        call checkdate(year,month,day,hour,min,sec,trim(calstr))

        call ESMF_TimeGet(time4,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
        call ESMF_TimeGet(time4,d=jday)
        write(6,F01) trim(estr2),'ti4 ',year,month,day,hour,min,sec,trim(calstr),jday,trim(dstr),delta
        call checkdate(year,month,day,hour,min,sec,trim(calstr))

        call ESMF_TimeGet(time5,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
        call ESMF_TimeGet(time5,d=jday)
        write(6,F01) trim(estr2),'ti5 ',year,month,day,hour,min,sec,trim(calstr),jday,trim(dstr),delta
        call checkdate(year,month,day,hour,min,sec,trim(calstr))

        call ESMF_TimeGet(time6,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
        call ESMF_TimeGet(time6,d=jday)
        write(6,F01) trim(estr2),'ti6 ',year,month,day,hour,min,sec,trim(calstr),jday,trim(dstr),delta
        call checkdate(year,month,day,hour,min,sec,trim(calstr))

        call ESMF_TimeGet(time7,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
        call ESMF_TimeGet(time7,d=jday)
        write(6,F01) trim(estr2),'ti7 ',year,month,day,hour,min,sec,trim(calstr),jday,trim(dstr),delta
        call checkdate(year,month,day,hour,min,sec,trim(calstr))

        call ESMF_TimeGet(time8,yy=year,mm=month,dd=day,h=hour,m=min,s=sec)
        call ESMF_TimeGet(time8,d=jday)
        write(6,F01) trim(estr2),'ti8 ',year,month,day,hour,min,sec,trim(calstr),jday,trim(dstr),delta
        call checkdate(year,month,day,hour,min,sec,trim(calstr))

        call ESMF_TimeGet(time1,yy=year1,mm=month1,dd=day1,h=hour1,m=min1,s=sec1)
        call ESMF_TimeGet(time1,d=jday1)
        call ESMF_TimeGet(time8,yy=year2,mm=month2,dd=day2,h=hour2,m=min2,s=sec2)
        call ESMF_TimeGet(time8,d=jday2)

        totcnt = totcnt + 1
        errfound = .false.

        if (time1 /= time3) then
           if (trim(dstr) == 'month' .or. trim(dstr) == 'all') then
              write(6,F03) 'ERROR: timediff non fatal',year1,month1,day1,trim(calstr),trim(dstr),delta
              if (.not. errfound) errcnt = errcnt + 1
              errfound = .true.
           else
              call wrf_error_fatal('ERROR: timeinc time')
           endif
        endif

        if (time3 /= time6) then
           if (trim(dstr) == 'month' .or. trim(dstr) == 'all') then
              write(6,F03) 'ERROR: time2x   non fatal',year1,month1,day1,trim(calstr),trim(dstr),delta
              if (.not. errfound) errcnt = errcnt + 1
              errfound = .true.
           else
              call wrf_error_fatal('ERROR: timeinc time')
           endif
        endif

        if (time6 /= time8) then
           if (trim(dstr) == 'month' .or. trim(dstr) == 'all') then
              write(6,F03) 'ERROR: timeneg  non fatal',year1,month1,day1,trim(calstr),trim(dstr),delta
              if (.not. errfound) errcnt = errcnt + 1
              errfound = .true.
           else
              call wrf_error_fatal('ERROR: timeinc time')
           endif
        endif

        if (time2 /= time5) then
           if (trim(dstr) == 'month' .or. trim(dstr) == 'all') then
              write(6,F03) 'ERROR: timecomp non fatal',year1,month1,day1,trim(calstr),trim(dstr),delta
              if (.not. errfound) errcnt = errcnt + 1
              errfound = .true.
           else
              call wrf_error_fatal('ERROR: timeinc time')
           endif
        endif

        if (year1 /= year2 .or. month1 /= month2 .or. day1 /= day2 .or. &
            hour1 /= hour2 .or. min1 /= min2 .or. sec1 /= sec2 .or. jday1 /= jday2) then
           if (trim(dstr) == 'month' .or. trim(dstr) == 'all') then
              write(6,F03) 'ERROR: ymdhms   non fatal',year1,month1,day1,trim(calstr),trim(dstr),delta
              if (.not. errfound) errcnt = errcnt + 1
              errfound = .true.
           else
              call wrf_error_fatal('ERROR: timeinc ymdhms')
           endif
        endif

     enddo
     enddo

  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  write(6,*) ' '
  write(6,*) 'tests run = ',totcnt,' error tests = ',errcnt
  write(6,*) 'esmf_wrf_timemgr test program completed successfully '
  write(6,*) ' '

  end program test


  subroutine checkdate(year,month,day,hour,min,sec,calstr)

  implicit none
  integer, intent(in) :: year,month,day,hour,min,sec
  character(len=*),intent(in) :: calstr
  INTEGER, PARAMETER :: mday(12)   &
                      = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  INTEGER, PARAMETER :: mdayleap(12) &
                      = (/31,29,31,30,31,30,31,31,30,31,30,31/)
  logical :: error

  error = .false.

  if (month < 1 .or. month > 12) error = .true.
  if (trim(calstr) == 'noleap') then
    if (day < 1 .or. day > mday(month)) error = .true.
  elseif (trim(calstr) == 'gregor') then
    if (day < 1 .or. day > mdayleap(month)) error = .true.
  else
    error = .true.
  endif
  if (hour < 0 .or. hour > 23) error = .true.
  if (min < 0 .or. min > 59) error = .true.
  if (sec < 0 .or. sec > 59) error = .true.

  if (error) then
     write(6,*) 'ERROR checkdate ',year,month,day,hour,min,sec,trim(calstr)
     call wrf_error_fatal('ERROR: checkdate')
  endif

  end subroutine checkdate
