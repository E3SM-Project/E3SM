
      subroutine TIMCON ( buff, &
                          time, &
                          lout )

      use IO, only : buffh

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)           ::    lout          ! output unit number
      character (len=*), intent(in) ::    buff          ! input time character string
      real, intent(out)             ::    time          ! converted time in seconds


!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer ::   retcod, l, i, j, k, slen
      integer ::   mnth, tokcnt
      integer ::   toklen(6)
      real    ::   units(6)
      real    ::   days
      real    ::   ndym(12) = (/ 31., 28., 31., 30., 31., 30., &
                                 31., 31., 30., 31., 30., 31. /)
      character(len=16) ::  tokens(6)
      character(len=16) ::  number
      character(len=3)  ::  timsym(6) = &
                    (/ 'Y  ', 'MON', 'D  ', 'H  ', 'MIN', 'S  ' /)

      units(:) = 0.
      slen = LEN_TRIM( buff )
      buffh = buff(:slen)
      call UPCASE( buffh )
      call GETTOKENS( buffh, slen, ':', 16, tokens, toklen, 6, tokcnt )
      if( tokcnt == 0 ) then
	 write(*,*) ' TIMCON : Improper input time string'
	 stop
      end if
      do i = 1,tokcnt
         l = VERIFY( tokens(i)(:toklen(i)), 'YMONDHIS', back = .true. )
	 if( l == 0 ) then
	    write(*,*) ' TIMCON : Improper input time number'
	    stop
	 end if
	 if( l == toklen(i) ) then
	    j = 6
	 else
	    l = l + 1
	    do j = 1,6
	       if( tokens(i)(l:toklen(i)) == timsym(j)(:LEN_TRIM(timsym(j))) ) then
	 	  l = l - 1
		  exit
	       end if
	    end do
	 end if

	 number = tokens(i)(:l)
         if( j /= 2 ) then
            call RELCON( number, l, units(j), retcod )
            if( retcod /= 0 ) then
               call ERRMES( 'number format error in time input #@', &
                           lout, &
                           number(:l), &
                           l, &
                           buff )
            end if
         else
            call INTCON( number, l, mnth, retcod )
            if( retcod /= 0 ) then
               call ERRMES( 'number format error in time input #@', &
                            lout, &
                            number(:l), &
                            l, &
                            buff )
            end if
            units(2) = SUM( ndym(:mnth-1) )
         end if
      end do

      time = 365.*units(1) + SUM( units(2:3) )
      time = time * 8.64e4
      time = time + 60.*(60.*units(4) + units(5)) + units(6)

      end subroutine TIMCON

      subroutine TIMCON_D ( buff, &
                            days0, &
                            secs )

      use IO, only : lout, buffh

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy arguments
!-----------------------------------------------------------------------
      character (len=*), intent(in) ::    buff          ! input time character string
      real, intent(out)             ::    days0         ! elapsed days since 0/0/0
      real, intent(out)             ::    secs          ! elapsed secs

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer ::   retcod, l, i, j, k, slen
      integer ::   mnth, tokcnt
      integer ::   toklen(6)
      real    ::   units(6)
      real    ::   days
      real    ::   ndym(12) = (/ 31., 28., 31., 30., 31., 30., &
                                 31., 31., 30., 31., 30., 31. /)
      character(len=16) ::  tokens(6)
      character(len=16) ::  number
      character(len=3)  ::  timsym(6) = &
                    (/ 'Y  ', 'MON', 'D  ', 'H  ', 'MIN', 'S  ' /)

      units(:) = 0.
      slen = LEN_TRIM(buff)
      buffh = buff(:slen)
      call UPCASE( buffh )
      call GETTOKENS( buffh, slen, ':', 16, tokens, toklen, 6, tokcnt )
      if( tokcnt == 0 ) then
	 write(*,*) ' TIMCON_D : Improper input time string'
	 stop
      end if
      do i = 1,tokcnt
         l = VERIFY( tokens(i)(:toklen(i)), 'YMONDHIS', back = .true. )
	 if( l == 0 ) then
	    write(*,*) ' TIMCON_D : Improper input time number'
	    stop
	 end if
	 if( l == toklen(i) ) then
	    j = 6
	 else
	    l = l + 1
	    do j = 1,6
	       if( tokens(i)(l:toklen(i)) == timsym(j)(:LEN_TRIM(timsym(j))) ) then
	 	  l = l - 1
		  exit
	       end if
	    end do
	 end if

	 number = tokens(i)(:l)
         if( j /= 2 ) then
            call RELCON( number, l, units(j), retcod )
            if( retcod /= 0 ) then
               call ERRMES( 'number format error in time input #@', &
                           lout, &
                           number(:l), &
                           l, &
                           buff )
            end if
         else
            call INTCON( number, l, mnth, retcod )
            if( retcod /= 0 ) then
               call ERRMES( 'number format error in time input #@', &
                            lout, &
                            number(:l), &
                            l, &
                            buff )
            end if
            units(2) = SUM( ndym(:mnth-1) )
         end if
      end do

      days0 = 365.*units(1) + SUM( units(2:3) )
      secs  = 60.*(60.*units(4) + units(5)) + units(6)

      end subroutine TIMCON_D

      subroutine CARDIN( lin, card, chars )
!-----------------------------------------------------------------------
!       ... Cardin reads on logical input unit lin an 80 character
!           card image right filled with blanks.
!           The image has all imbedded whitespace removed and
!           non whitespace character count returned in chars
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)  ::   lin
      integer, intent(out) ::   chars
      character(len=*), intent(out) ::  card

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer, parameter ::  ht = 9
      integer      ::  i, ios, slen
      character(len=6) ::  format
      logical      ::  compress


      format = ' '
      slen = LEN( card ) 
      if( slen < 10 ) then
         write(format,'(''(a'',i1,'')'')') slen
      else if( slen < 100 ) then
         write(format,'(''(a'',i2,'')'')') slen
      else if( slen < 1000 ) then
         write(format,'(''(a'',i3,'')'')') slen
      end if
      do
         compress = .true.
	 card = ' '
         read(lin,format,iostat=ios) card
	 if( ios /= 0 ) then
	    write(*,*) ' CARDIN : Read error = ',ios 
            call ERRMES( 'Read Error in driver file@', 6, card, 1, card )
	 end if
!-----------------------------------------------------------------------
!       ... Remove blanks and horizontal tabs
!-----------------------------------------------------------------------
	 slen = LEN_TRIM( card )
	 if( slen == 0 ) then
	    cycle
	 end if
         chars = VERIFY( card(:slen), ' ' )
	 if( card(chars:chars) == '*' ) then
	    cycle
	 end if
         chars = 0
         do i = 1,slen
            if( compress ) then
               if( card(i:i) /= ' ' .and.  ICHAR( card(i:i) ) /= ht ) then
                  if( card(i:i) == '"' ) then
                     compress = .false.
                     cycle
                  end if
                  chars = chars + 1
                  card(chars:chars) = card(i:i)
               end if
            else
               if( card(i:i) == '"' ) then
                  compress = .true.
                  cycle
               end if
               chars = chars + 1
               card(chars:chars) = card(i:i)
            end if
         end do
	 if( card(chars:chars) == ',' ) then
	    card(chars:chars) = ' '
	    chars = chars - 1
	    if( chars == 0 ) then
               cycle
	    end if
	 end if
!-----------------------------------------------------------------------
!       ... Ignore "blank" or comment card
!-----------------------------------------------------------------------
         if( card == ' ' ) then
            cycle
         else
            card(chars+1:) = ' '
            exit
         end if
      end do

      end subroutine CARDIN

      subroutine ERRMES( string, &
                         lout, &
                         instng, &
                         count, &
                         card )
!-----------------------------------------------------------------------
!     	... Prints the input string error message;
!           stops the pre-processor
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::  count
      integer, intent(in) ::  lout
      character(len=*), intent(in) ::  string, instng, card

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::  ls, i
      character(len=320) :: copy

      copy = '0 *** '
      ls  = 6

      do i = 1,MIN( LEN(string),120 )
         if( string(i:i) == '@' ) then
            exit
         else if( string(i:i) == '#' ) then
            copy(ls+1:ls+count) = instng(:count)
            ls = ls + count
         else
            ls = ls + 1
            copy(ls:ls) = string(i:i)
         end if
      end do

      write(lout,'(a)') copy(:ls)
      write(   *,'(a)') copy(:ls)
      if( card /= ' ' ) then
         write(lout,'('' Input line:'')')
         write(   *,'('' Input line:'')')
         write(lout,'(1x,a80)') card
         write(   *,'(1x,a80)') card
      end if

      stop 'abort'

      end subroutine ERRMES

      subroutine ALTCON( string, &
                         ls, &
                         alt, &
                         retcod )
!-----------------------------------------------------------------------
!        altcon converts the input character string altcon
!        to a real number returned in alt.  The input string
!        must have length ls and can be of the following forms:
!                  1.  %km
!                  2.  %
!        where % is a generalized e format
!        Successful conversion returns zero in retcod and -12
!        otherwise
!
!     inputs:
!
!       string  =  character string to convert
!       ls      =  length of string (max value = 80)
!
!     outputs:
!
!       nout    =  integer value if retcod = 0
!       retcod  =  error flag; = 0 => proper format for string
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)  ::  ls
      integer, intent(out) ::  retcod
      real, intent(out)    ::  alt
      character(len=*), intent(in) ::  string

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::  len

      if( string(ls:ls) == 'm' .or.  string(ls:ls) == 'M' ) then
        len = ls - 1
      else
        len = ls
      end if

      if( len /= 0 ) then
         call RELCON( string, &
                      len, &
                      alt, &
                      retcod )
         if( retcod == 0 ) then
            if( len == ls ) then
               alt = alt*1.e5
            else
               alt = alt*1.e2
            end if
            retcod = 0
            return
         end if
      end if

      retcod = -12

      end subroutine ALTCON

      subroutine INTCON( string, &
                         ls, &
                         nout, &
                         retcod )
!-----------------------------------------------------------------------
!     intcon converts a character string of length ls to an integer
!     format errors are trapped and retcod is set to -12
!
!     inputs:
!
!       string  =  character string to convert
!       ls      =  length of string (max value = 80)
!
!     outputs:
!
!       nout    =  integer value if retcod = 0
!       retcod  =  error flag; = 0 => proper format for string
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)  :: ls
      integer, intent(out) :: nout, retcod
      character(len=*), intent(in) :: string

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::  ios

      read (string(:ls),*,iostat=ios) nout
      if( ios == 0 ) then
         retcod = 0
      else
         retcod = -12
      end if

      end subroutine INTCON

      subroutine RELCON( string, &
                         ls, &
                         flpout, &
                         retcod )
!-----------------------------------------------------------------------
!     relcon converts a character string of length ls to a real number
!     format errors are trapped and retcod is set to -12
!
!     inputs:
!
!       string  =  character string to convert
!       ls      =  length of string (max value = 80)
!
!     outputs:
!
!       flpout  =  real value if retcod = 0
!       retcod  =  error flag; = 0 => proper format for string
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)  ::  ls
      integer, intent(out) ::  retcod
      real, intent(out)    ::  flpout
      character(len=*), intent(in) ::  string

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::  ios, l

      l = LEN_TRIM( string(:ls) )
      if( l == 0 ) then
	 retcod = -12
      else
         read(string(:l),*,iostat=ios) flpout
         if( ios == 0 ) then
	    retcod = 0
         else
	    retcod = -12
         end if
      end if

      end subroutine RELCON

      subroutine NUMCON( string, &
                         num, &
                         jus )
!-----------------------------------------------------------------------
!     numcon converts a real number to a generalized
!     e format character string
!
!     inputs:
!
!      num      =  real number to convert 
!      jus      =  character code for output string justification
!                  'l' = left justify
!                  'c' = center
!                  'r' = right justify
!
!     outputs:
!
!       string  =  converted character string 
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      real, intent(in) ::  num
      character(len=*), intent(out) :: string
      character(len=1), intent(in)  :: jus

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::  wpart, il, i, iu
      real     ::  frac
      character(len=16) :: mask
      character(len=8)  :: copy(2)

      equivalence (mask,copy)

      if( num == 0. ) then
         string  =  '0'
      else
         wpart  =  INT( num )
         if( wpart == 0 ) then
            il  =  9
         else
            write(copy(1),'(i8)') wpart
            do i  =  1,8
              if( copy(1)(i:i) /= ' ') then
                go to 12
              end if
            end do
12          il  =  i
         end if
         frac  =  num - INT( num )
         if( frac == 0. ) then
            iu  =  8
         else
            if( wpart == 0 .and. frac < 0.e0 ) then
               write(copy(2),'(f8.6)') frac
            else
               write(copy(2),'(f8.7)') ABS( frac )
            end if
            do i = 8,2,-1
               if( copy(2)(i:i) /= '0') then
                 go to 22
               end if
            end do
22          iu = i + 8
         end if
         string = mask(il:iu)
      end if
      if( jus == 'c' ) then
         call CENTER( string, 16 )
      end if

      end subroutine NUMCON

      subroutine CENTER( string, ls )
!-----------------------------------------------------------------------
!     	... Center a character string of length ls
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!     	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::  ls
      character(len=*), intent(inout) ::    string

!-----------------------------------------------------------------------
!     	... Local variables
!-----------------------------------------------------------------------
      integer  ::  offset, i, il, iu, len, j
      character(len=320) :: copy

      len = MIN( ls,120 )
      copy(:len) = string(:len)
      il = VERIFY( string(:len), ' ' )
      if( il /= 0) then
	 iu = VERIFY( string(:len), ' ', back = .true. )
         len = iu - il + 1
         offset = MAX( 0,(ls - len)/2 + 1)
         string(:ls) = ' '
         string(offset:offset+len) = copy(il:iu)
      end if

      end subroutine CENTER

      integer function LENOF( ls, string )
!-----------------------------------------------------------------------
!    	... Returns the length of a string by finding
!           the first non-blank character from the right side
!           of the input string.
!           This function will not scan beyond ls characters
!           in the input string.  If the string consists of only
!           blanks then an "error" value of zero is returned.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::  ls
      character(len=*), intent(in) ::  string

      LENOF = LEN_TRIM( string(:ls) )

      end function LENOF
      
      integer function ALTCHK( alt, sptgrd, np )
      
      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::  np
      real, intent(in)    ::  alt
      real, intent(in)    ::  sptgrd(np)
      
!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: i
      
      ALTCHK = 0
      do i = 1,np
         if( sptgrd(i) == alt ) then
            ALTCHK = 1
            exit
         end if
      end do
      
      end function ALTCHK

      subroutine UPCASE( lstring )
!----------------------------------------------------------------------
!       ... Convert character string lstring to upper case
!----------------------------------------------------------------------
      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      character(len=*), intent(inout) ::  lstring

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: i

      do i = 1,LEN_TRIM( lstring )
         if( ICHAR(lstring(i:i)) >= 97 .and.  ICHAR(lstring(i:i)) <= 122 ) then
            lstring(i:i) = CHAR(ICHAR(lstring(i:i)) - 32)
         end if
      end do

      end subroutine UPCASE

      integer function STRLEN ( string )
!-----------------------------------------------------------------------
!  	... Returns the length of a string by finding
!           the first non-blank character from the right side
!           of the input string.
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      character(len=*), intent(in) ::  string

      STRLEN = LEN_TRIM( string )

      end function STRLEN

      subroutine PARSE_FLPTH( fullpath, filename, filepath )

      implicit none

!-------------------------------------------------------
!       ... Dummy args
!-------------------------------------------------------
      character(len=*), intent(in) ::    fullpath   ! incoming full pathname
      character(len=*), intent(out) ::   filename   ! the file name
      character(len=*), intent(out) ::   filepath   ! the file path

!-------------------------------------------------------
!       ... Local variables
!-------------------------------------------------------
      integer  ::  i

      i = INDEX( fullpath(:LEN_TRIM(fullpath)), '/', back = .true. )
      filename = fullpath(i+1:)
      if( i /= 0 ) then
         filepath = fullpath(:i) 
      else
         filepath = ' '
      end if

      end subroutine PARSE_FLPTH

      logical function ISNUM( char )

      implicit none

      character(len=1), intent(in) ::  char

      if( char <= '9' .and. char >= '0' ) then
         ISNUM = .true.
      else
         ISNUM = .false.
      end if

      end function ISNUM

      subroutine MKDATE( time, date )

      implicit none

!-------------------------------------------------------
!       ... Dummy args
!-------------------------------------------------------
      real, intent(in) ::     time             ! time to convert in days
      character(len=6), intent(out) :: date    ! date in form yymmdd

!-------------------------------------------------------
!       ... Local variables
!-------------------------------------------------------
      integer  :: years, days
      integer  :: mnth
      integer  :: mdys(0:12) = (/ 0, 31, 59, 90, 120, 151, &
				  181, 212, 243, 273, 304, 334, 365 /)

      years = INT( time/365. )
      if( years > 1900 ) then
         years = years - 1900
      end if
      write(date(1:2),'(i2)') years

      days = INT( time - REAL(years)*365. )
      do mnth = 1,12
         if( days <= mdys(mnth) ) then
            exit
         end if
      end do

      mnth = MAX( 1,MIN( 12, mnth ) )
      if( mnth < 10 ) then
         write(date(3:4),'(''0'',i1)') mnth
      else
         write(date(3:4),'(i2)') mnth
      end if

      days = days - mdys(mnth-1)
      if( days < 10 ) then
         write(date(5:6),'(''0'',i1)') days
      else
         write(date(5:6),'(i2)') days
      end if
      
      end subroutine MKDATE

      integer function INCLIST( target, list, cnt )
!-------------------------------------------------------
!       ... Check for match in character list
!-------------------------------------------------------

      implicit none

!-------------------------------------------------------
!       ... Input arguments
!-------------------------------------------------------
      integer, intent(in) ::    cnt      ! no elements in list
      character(len=*), intent(in) ::  target   ! match string
      character(len=*), intent(in) ::  list(*)  ! list to search

!-------------------------------------------------------
!       ... Local variables
!-------------------------------------------------------
      integer :: i

      INCLIST = 0
      do i = 1,cnt
         if( target == list(i) ) then
            INCLIST = i
            exit
         end if
      end do

      end function INCLIST

      integer function INILIST( target, list, cnt )
!-------------------------------------------------------
!       ... Check for match in integer list
!-------------------------------------------------------

      implicit none

!-------------------------------------------------------
!       ... Input arguments
!-------------------------------------------------------
      integer, intent(in) ::   cnt      ! no elements in list
      integer, intent(in) ::   target   ! match integer
      integer, intent(in) ::   list(*)  ! list to search

!-------------------------------------------------------
!       ... Local variables
!-------------------------------------------------------
      integer :: i

      INILIST = 0
      do i = 1,cnt
         if( target == list(i) ) then
            INILIST = i
            exit
         end if
      end do

      end function INILIST

      subroutine r2c( string, num, jus )
!-----------------------------------------------------------------------
!     r2c converts a real number to a generalized e format character string
!
!     inputs:
!
!      num      =  real number to convert 
!      jus      =  character code for output string justification
!                  'l' = left justify
!                  'c' = center
!                  'r' = right justify
!
!     outputs:
!
!       char    =  converted character string 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      real, intent(in) ::  num
      character(len=*), intent(out) :: string
      character(len=1) ::  jus

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      real, parameter ::  epsilon = .0000005 

      integer  ::  wpart, il, i, iu
      integer  ::  power
      real     ::  frac
      real     ::  wrk_num
      character(len=24) :: mask
      character(len=8)  :: copy(3)

      equivalence (mask,copy)

      if( num == 0. ) then
         string  =  '0'
      else
         wrk_num = num
         power = 0
         if( abs(wrk_num) < 1.e-4 ) then
            power = 1
            wrk_num = wrk_num
            do while( power < 40  )
               wrk_num = 10.*wrk_num
               if( abs(wrk_num) >= 1. ) then
                  go to 100
               end if
               power = power + 1
            end do
            string = '0'
            return
         end if

100      wrk_num = wrk_num + epsilon
         wpart   = int( wrk_num )
         if( wpart == 0 ) then
            il  =  9
         else
            write(copy(1),'(i8)') wpart
            do i = 1,8
               if( copy(1)(i:i) /= ' ') then
                  exit
               end if
            end do
            il = i
         end if
         frac  =  wrk_num - int(wrk_num)
         if( frac == 0. ) then
            mask(9:9) = '.'
            iu  =  9
         else
            if( wpart == 0 .and. frac < 0.e0 ) then
               write(copy(2),'(f8.6)') frac
            else
               write(copy(2),'(f8.7)') abs(frac)
            end if
            do i = 8,2,-1
               if( copy(2)(i:i) /= '0') then
                  exit
               end if
            end do
            iu = i + 8
         end if
         if( frac /= 0. ) then
            if( mask(iu-3:iu-1) == '000') then
               iu = iu - 4
            else if( mask(iu-4:iu-2) == '000' ) then
               iu = iu - 5
            end if
         end if
         if( power /= 0 ) then
            if( power < 10 ) then
               write(mask(iu+1:),'(''e'',i2)') -power
               iu = iu + 3
            else
               write(mask(iu+1:),'(''e'',i3)') -power
               iu = iu + 4
            end if
         end if
         if( num > 0. ) then
            string  =  mask(il:iu)
         else
            string  =  '(' // mask(il:iu) // ')'
         end if
      end if

      if( jus == 'c' ) then
         call center( string, 16 )
      end if

      end subroutine r2c

      integer function XLATE( match )
!------------------------------------------------------------------------
!	... Translate between overall indexing and method indexing
!------------------------------------------------------------------------
     
      use VAR_MOD, only : var_lim, clsmap

      implicit none

!------------------------------------------------------------------------
!	... Dummy args
!------------------------------------------------------------------------
      integer, intent(inout) ::     match
      
!------------------------------------------------------------------------
!	... Local variables
!------------------------------------------------------------------------
      integer  :: class
      
      do class = 1,5
        if( clsmap(match,class,1) /= 0 ) then
           match = clsmap(match,class,1)
           XLATE = class
           return
        end if
      end do
      
      XLATE = 0
      
      end function XLATE
