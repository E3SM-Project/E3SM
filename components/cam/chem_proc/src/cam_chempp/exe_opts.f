      subroutine EXE_OPTS( options, &
                           lin, &
                           lout )
!-----------------------------------------------------------------------
!	... Set the execution options
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)  ::  lin
      integer, intent(in)  ::  lout
      logical, intent(out) ::  options(3)

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::  kpar, nchar, k
      integer  ::  parsw(3)

      character(len=80) :: buff
      character(len=20) :: parkey(3), keywrd
      logical :: found

      integer  ::  LENOF

      parkey(1) = 'QSUBFILE'
      parkey(2) = 'SUBMIT'
      parkey(3) = 'FIXER'

      parsw = 0

!-----------------------------------------------------------------------
!   	... Scan for valid option keyword
!-----------------------------------------------------------------------
      do
         call CARDIN( lin, buff, nchar )
	 call UPCASE ( buff )
         if( buff == 'ENDEXECUTIONOPTIONS' ) then
            exit
	 end if
	 k = INDEX( buff(:nchar), '=' )
         if( k /= 0 ) then
	    keywrd = buff(:k-1)
	    found = .false.
            do kpar = 1,6
               if( keywrd == parkey(kpar) ) then
		  found = .true.
	          exit
	       end if
	    end do
	    if( .not. found ) then
               call ERRMES ( ' # is an invalid options' &
                             // ' parameter keyword@', lout, keywrd, &
                             LENOF(20,keywrd), buff )
	    end if
	 else
!-----------------------------------------------------------------------
!  	... Invalid parameter keyword; terminate the program
!-----------------------------------------------------------------------
            call ERRMES( ' option specification has no = operator@', &
                         lout, buff, 1, buff )
         end if

!-----------------------------------------------------------------------
!     	... Valid parameter keyword; now check for duplicate keyword
!-----------------------------------------------------------------------
         if( parsw(kpar) /= 0 ) then
            call ERRMES( '0 *** # has already been specified@', &
                          lout, parkey(kpar), k, ' ' )
         end if

!-----------------------------------------------------------------------
!     	... Set individual options
!-----------------------------------------------------------------------
	 if( buff(k+1:nchar) == 'ON' .or. &
             buff(k+1:nchar) == 'YES' ) then
	    options(kpar) = .true.
	 else
	    options(kpar) = .false.
	 end if
	 parsw(kpar) = 1
      end do

      end subroutine EXE_OPTS
