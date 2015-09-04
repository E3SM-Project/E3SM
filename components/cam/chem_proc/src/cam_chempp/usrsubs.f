
      subroutine USRSUBS ( sub_names, sub_cnt )

      use IO

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(out) ::  sub_cnt                    ! count of user subroutines
      character(len=128), intent(out) :: sub_names(*)      ! user filenames

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer, parameter ::  symlen = 64

      integer  ::  sublim = 100
      integer  ::  nchar, pos
      integer  ::  toklen(20)
      integer  ::  j, count
      integer  ::  no_tokens

      character(len=64) :: filepath
      character(len=128) :: filespec
      character(len=64) :: tokens(20)

      logical  ::  lexist

      sub_cnt = 0
      count   = 0
      filepath = ' '
      filespec = ' '
!-----------------------------------------------------------------------
!       ... Read the subroutine pathnames
!-----------------------------------------------------------------------
      do
         call CARDIN ( lin, buff, nchar )
         buffh = buff
         call UPCASE( buffh )
         if( buffh /= 'ENDUSERSUBROUTINES' ) then
            call GETTOKENS(  buff,     nchar,    ',',      symlen, &
                             tokens,   toklen,   20,       no_tokens )
            if( no_tokens == 0 ) then
               call ERRMES( ' Files input line in error@', lout, buff, 1, ' ' )
!-----------------------------------------------------------------------
!       ... Check for filepath setting
!-----------------------------------------------------------------------
	    else if( no_tokens == 1 .and. tokens(1)(9:9) == '='  ) then
	       filepath = tokens(1)(:8)
	       call UPCASE( filepath )
	       if( filepath == 'FILEPATH' ) then
		  filepath = tokens(1)(10:)
		  filepath = TRIM( filepath )
	       else
                  call ERRMES( ' # is not FILEPATH keyword@', lout, tokens(1)(:8), 8, buff )
	       end if
	       cycle
	    end if
!-----------------------------------------------------------------------
!       ... Process the user subroutine filespec
!-----------------------------------------------------------------------
            do j = 1,no_tokens
               count = count + 1
               if( count > sublim ) then
                  call ERRMES( ' Files count exceeds limit@', lout, buff, 1, buff )
               end if
	       if( tokens(j)(1:1) == '/' ) then
	          filespec = tokens(j)(:toklen(j))
	       else
	          if( filepath /= ' ' ) then
	             pos = LEN_TRIM(filepath)
	             if( filepath(pos:pos) /= '/' ) then
		        filepath(pos+1:pos+1) = '/'
		     end if
		     filespec = TRIM( filepath) // tokens(j)(:toklen(j))
	          else
	             filespec = tokens(j)(:toklen(j))
	          end if
	       end if
!-----------------------------------------------------------------------
!       ... Check for file existence
!-----------------------------------------------------------------------
	       INQUIRE( file  = TRIM( filespec ), exist = lexist )
	       if( .not. lexist ) then
                  call ERRMES( ' File # does NOT exist@', lout, filespec, LEN_TRIM(filespec), buff )
	       end if
	       sub_names(count) = filespec
            end do
            cycle
         else
            sub_cnt = count
	    exit
         end if
      end do

      end subroutine USRSUBS
