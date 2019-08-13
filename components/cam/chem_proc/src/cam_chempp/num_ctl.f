
      subroutine NUM_CTL( iter_counts )

      use IO

      implicit none

!-----------------------------------------------------------------------
!	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(inout) ::      iter_counts(4)          ! iteration counts

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: max_parm = 4
      integer  ::  kpar, nchar, k
      integer  ::  parsw(max_parm)
      integer  ::  retcod
      character(len=20) :: keywrd
      character(len=20) :: parkey(max_parm)
      logical  :: found

      parkey(1) = 'HOVITERATIONS'
      parkey(2) = 'IMPLICITITERATIONS'
      parkey(3) = 'JACOBIANITERATIONS'
      parkey(4) = 'EBIITERATIONS'

      parsw = 0

!-----------------------------------------------------------------------
!   	... Scan for valid option keyword
!-----------------------------------------------------------------------
      do
         call CARDIN( lin, buff, nchar )
	 buffh = buff
	 call UPCASE ( buffh )
         if( buffh == 'ENDNUMERICALCONTROL' ) then
            exit
	 end if
	 k = INDEX( buffh(:nchar), '=' )
         if( k /= 0 ) then
	    found = .false.
	    keywrd = buffh(:k-1)
            do kpar = 1,max_parm
               if( keywrd == parkey(kpar) ) then
	          found = .true.
		  exit
	       end if
	    end do
	 else
!-----------------------------------------------------------------------
!  	... Invalid parameter keyword; terminate the program
!-----------------------------------------------------------------------
            call ERRMES ( ' Num ctl specification has no = operator@', lout, buff, 1, buff )
         end if
!-----------------------------------------------------------------------
!  	... Invalid parameter keyword; terminate the program
!-----------------------------------------------------------------------
	 if( .not. found ) then
            call ERRMES ( ' # is an invalid Num control parameter keyword@', &
                          lout, &
                          keywrd, &
                          LEN_TRIM(keywrd), &
                          buffh )
!-----------------------------------------------------------------------
!     	... Valid parameter keyword; now check for duplicate keyword
!-----------------------------------------------------------------------
         else if( parsw(kpar) /= 0 ) then
            call ERRMES( '0 *** # has already been specified@', lout, parkey(kpar), k, ' ' )
         end if

!-----------------------------------------------------------------------
!     	... Set individual iteration counts
!-----------------------------------------------------------------------
	 call INTCON( buff(k+1:nchar), nchar - k, iter_counts(kpar), retcod )
!-----------------------------------------------------------------------
!     	... Check itertion limit for validity
!-----------------------------------------------------------------------
         if( retcod /= 0 ) then
            call ERRMES ( ' # is an invalid iteration count@', &
                          lout, &
                          buff(k+1:nchar), &
                          nchar - k, &
                          buffh )
         else if( iter_counts(kpar) <= 0 ) then
            call ERRMES ( ' # is an invalid iteration count@', &
                          lout, &
                          buff(k+1:nchar), &
                          nchar - k, &
                          buffh )
	 end if
	 parsw(kpar) = 1
      end do

      end subroutine NUM_CTL
