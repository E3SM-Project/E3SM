
      module MO_SPAT_DIMS

      CONTAINS

      subroutine SPAT_DIMS( buff, dimensions )
!-----------------------------------------------------------------------
!   	... Set the simulation spatial dimensions
!-----------------------------------------------------------------------

      use IO, only : lin, lout

      implicit none

!-----------------------------------------------------------------------
!   	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(inout)           ::  dimensions(6)
      character(len=80), intent(inout) ::  buff

!-----------------------------------------------------------------------
!   	... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: maxparms = 6
      integer  ::  kpar, nchar, retcod, i, k
      character(len=20) :: parkey(maxparms),     keywrd
      logical  ::  found
      logical  ::  processed(maxparms)

!-----------------------------------------------------------------------
!   	... Function declarations
!-----------------------------------------------------------------------
      integer  ::  LENOF

      parkey(1) = 'LONGITUDEPOINTS'
      parkey(2) = 'LATITUDEPOINTS'
      parkey(3) = 'VERTICALPOINTS'
      parkey(4) = 'NXPT'
      parkey(5) = 'JINTMX'
      parkey(6) = 'PLONL'

      processed = .false.
!-----------------------------------------------------------------------
!   	... Scan for valid numerical control parameter keyword
!-----------------------------------------------------------------------
      do
         call CARDIN( lin, buff, nchar )
	 call UPCASE( buff )
         if( buff == 'ENDSPATIALDIMENSIONS' ) then
            exit
	 end if
	 k = INDEX( buff(:nchar), '=' )
         if( k /= 0 ) then
	    keywrd = buff(:k-1)
            found = .false.
            do kpar = 1,maxparms
               if( keywrd == parkey(kpar) ) then
		  found = .true.
		  exit
	       end if
	    end do
	    if( .not. found ) then
               call ERRMES ( ' # is an invalid numerical control' &
                          // ' parameter keyword@', lout, keywrd, &
                             LENOF(20,keywrd), buff )
	    end if
         else
!-----------------------------------------------------------------------
!  	... Invalid parameter keyword; terminate the program
!-----------------------------------------------------------------------
            call ERRMES ( ' numerical specification has no = operator@', &
                          lout, buff, 1, buff )
         end if

!-----------------------------------------------------------------------
!     	... Valid parameter keyword; now check for duplicate keyword
!-----------------------------------------------------------------------
         if( processed(kpar) ) then
            call ERRMES( '0 *** # has already been specified@', &
                          lout, parkey(kpar), k, ' ' )
         end if
         call INTCON ( buff(k+1:), &
                       nchar-k, &
                       dimensions(kpar), &
                       retcod )
!-----------------------------------------------------------------------
!     	... Check for numeric parameter syntax error
!-----------------------------------------------------------------------
         if( retcod /= 0 ) then
            call ERRMES ( ' # is an invalid real or integer in ' &
                        // 'numeric controls@', lout, buff(k+1:), &
                        LENOF( nchar-k, buff(k+1:)), buff )
         end if
         processed(kpar) = .true.
      end do

      end subroutine SPAT_DIMS

      end module MO_SPAT_DIMS
