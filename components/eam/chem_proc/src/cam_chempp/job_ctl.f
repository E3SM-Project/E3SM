      subroutine JOB_CTL( lin, &
                          lout, &
                          jobctl )

      implicit none

!-----------------------------------------------------------------------
!	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)            ::  lin,      lout
      character(len=16), intent(out) ::  jobctl(8)          ! job control variables

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::  kpar, nchar, k
      integer  ::  parsw(8)

      real     ::  time

      character(len=80) :: buff
      character(len=80) :: buffh
      character(len=20) :: parkey(8),     keywrd

      logical  :: found

      integer  ::  LENOF

      parkey(1) = 'SIMULATIONTIMESTEP'
      parkey(2) = 'CRAYTIMELIMIT'
      parkey(3) = 'SIMULATIONLENGTH'
      parkey(4) = 'CRAYMEMORY'
      parkey(5) = 'ACCOUNT'
      parkey(6) = 'CASE'
      parkey(7) = 'RESTART'
      parkey(8) = 'CRAYQUE'

      parsw  = 0

!-----------------------------------------------------------------------
!   	... Scan for valid option keyword
!-----------------------------------------------------------------------
      do
         call CARDIN( lin, buff, nchar )
	 buffh = buff
	 call UPCASE ( buffh )
         if( buffh == 'ENDJOBCONTROL' ) then
            exit
	 end if
	 k = INDEX( buffh(:nchar), '=' )
         if( k /= 0 ) then
	    keywrd = buffh(:k-1)
	    found = .false.
            do kpar = 1,8
               if( keywrd == parkey(kpar) ) then
		  found = .true.
		  exit
	       end if
	    end do
	 else
!-----------------------------------------------------------------------
!  	... Invalid parameter keyword; terminate the program
!-----------------------------------------------------------------------
            call ERRMES ( ' job ctl specification has no = operator@', lout, buff, 1, buff )
         end if
	 if( .not. found) then
            call ERRMES ( ' # is an invalid job control parameter keyword@', &
                          lout, &
                          keywrd, &
                          LENOF(20,keywrd), &
                          buffh )
         end if

!-----------------------------------------------------------------------
!     	... Valid parameter keyword; now check for duplicate keyword
!-----------------------------------------------------------------------
         if( parsw(kpar) /= 0 ) then
            call ERRMES( '0 *** # has already been specified@', lout, parkey(kpar), k, ' ' )
         end if

!-----------------------------------------------------------------------
!     	... Set individual options
!-----------------------------------------------------------------------
	 if( kpar <= 3 ) then
	    if( kpar == 3 ) then
	       if( buffh(nchar-4:nchar) == 'STEPS' ) then
	          jobctl(3) = buff(k+1:nchar-5)
	       else
	          call TIMCON( buff(k+1:nchar), time, lout )
	          jobctl(3) = buff(k+1:nchar)
	       end if
	    else
	       call TIMCON( buff(k+1:nchar), time, lout )
	       jobctl(kpar) = buff(k+1:nchar)
	    end if
	 else
	    jobctl(kpar) = buff(k+1:nchar)
	 end if
	 parsw(kpar) = 1
      end do

      end subroutine JOB_CTL
