      subroutine GETTOKENS( string, &
                            ls, &
                            delim, &
                            maxlen, &
                            tokens, &
                            toklen, &
                            maxtok, &
                            tokcnt )

      implicit none

!-----------------------------------------------------------------------     
!     Input arguments:
!        string - character string to crack into tokens
!        ls     - length of string
!        delim  - token delimiter character
!        maxlen - maximum length of any single token
!        maxtok - maximum number of tokens
!     Output arguments:
!        tokcnt - number of actual tokens
!                 < 0 => hit maxtok before end of string
!                 = 0 => error in input string
!        toklen - array containing length of each token
!        tokens - character array of tokens
!-----------------------------------------------------------------------     

      integer, intent(in)  ::  ls, maxlen, maxtok
      integer, intent(out) ::  tokcnt
      integer, intent(out) ::  toklen(*)
      
      character(len=*), intent(in)  :: string
      character(len=*), intent(out) :: tokens(*)
      character(len=1), intent(in)  :: delim
      
!-----------------------------------------------------------------------     
!	... Local variables
!-----------------------------------------------------------------------     
      integer  ::   marker, i, length, endpos

      tokcnt = 0
      marker = 1
      do i = 1,ls
         if( string(i:i) == delim .or. i == ls ) then
            if( i == ls ) then
               if( string(i:i) == delim ) then
                  tokcnt = 0
                  exit
               end if
               length = i - marker + 1
               endpos = i
            else
               length = i - marker
               endpos = i - 1
            end if
            if( length < 1 .or. length > maxlen ) then
               tokcnt = 0
               exit
            end if
            tokcnt = tokcnt + 1
            if( tokcnt > maxtok ) then
               tokcnt = -tokcnt
               exit
            end if
            tokens(tokcnt) = ' '
            tokens(tokcnt)(:length) = string(marker:endpos)
            toklen(tokcnt) = length
            marker = i + 1
         end if
      end do
      
      end subroutine GETTOKENS
