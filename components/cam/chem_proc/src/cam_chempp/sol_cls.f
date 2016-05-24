
      subroutine SOL_CLS( iout )
!-----------------------------------------------------------------------
!	... Map solution species to solution method groups
!-----------------------------------------------------------------------

      use IO
      use VAR_MOD, only : spccnt => new_nq, spcsym => new_solsym, clscnt, clsmap

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      character(len=80), intent(inout) :: iout(*)

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer, parameter ::  symlen = 16

      integer  ::  kpar, i, parsw(5), nchar
      integer  ::  toklen(20)
      integer  ::  j, l
      integer  ::  no_tokens
      integer  ::  class
      character(len=16) :: tokens(20)
      character(len=10) :: clshdr(5) =  (/ 'EXPLICIT  ', 'EBI       ', &
					   'HOV       ', 'IMPLICIT  ', &
					   'RODAS     ' /)
      character(len=11) :: clsend(5) =  (/ 'ENDEXPLICIT', 'ENDEBI     ', &
					   'ENDHOV     ', 'ENDIMPLICIT', &
					   'ENDRODAS   ' /)
      character(len=1)  :: char
      logical  :: found

      integer :: INILIST

!-----------------------------------------------------------------------
!	... Initialization
!-----------------------------------------------------------------------
      parsw(:)  = 0 ; clscnt(:) = 0 ; clsmap(:,:,:) = 0

      call CARDIN( lin, buff, nchar )
      buffh = buff
      call UPCASE( buffh )
      if( buffh /= 'SOLUTIONCLASSES' ) then
         call ERRMES( '"Solution classes" card missing; run terminated@', &
                      lout, char, 1, buff )
      end if

      do
         call CARDIN(lin, buff, nchar )
         buffh = buff
         call UPCASE( buffh )
         if( buffh == 'ENDSOLUTIONCLASSES' ) then
!-----------------------------------------------------------------------
!       ... Check for all species in class
!-----------------------------------------------------------------------
	    if( sum( clscnt(:) ) /= spccnt ) then
	       write(lout,*) ' '
	       write(lout,*) 'Following species not in a class'
	       write(lout,*) ' '
	       do l = 1,spccnt
		  found = .false.
	          do class = 1,5
		     if( clscnt(class) /= 0 ) then
		        j = INILIST( l, clsmap(1,class,2), clscnt(class) )
		        if( j /= 0 ) then
			   found = .true.
			   exit
		        end if
		     end if
	          end do
		  if( .not. found ) then
	             write(lout,*) trim(spcsym(l))
	          end if
	       end do
	       stop 'abort'
	    end if
            exit
         end if

	 found = .false.
         do kpar = 1,5
           if( buffh == clshdr(kpar) ) then
	      found = .true.
	      exit
	   end if
         end do
	 if( .not. found ) then
            call ERRMES( '# is an invalid class header@',  &
                         lout, &
                         buff(:8), &
                         8, &
                         buff )
         else if( parsw(kpar) /= 0 ) then
            call ERRMES( '# solution class already declared@', &
                      lout, &
                      clshdr(kpar), &
                      LEN_TRIM(clshdr(kpar)), &
                      buff )
         else
            parsw(kpar) = 1
         end if

!-----------------------------------------------------------------------
!       ... Read the solution class members
!-----------------------------------------------------------------------
Methods : &
         do
            call CARDIN(lin, buff, nchar)
            buffh = buff
            call UPCASE( buffh )
            if( buffh /= clsend(kpar) ) then
               if( buffh(:nchar) == 'ALL' ) then
		  clscnt(:5) = 0
                  clscnt(kpar) = spccnt
		  clsmap(:,:,:) = 0
                  do j = 1,spccnt
                     clsmap(j,kpar,1) = j
                     clsmap(j,kpar,2) = j
                  end do
                  cycle
	       else if( buffh(:nchar) == 'ALLOTHERS' ) then
		  clscnt(kpar) = 0
		  clsmap(:,kpar,:) = 0
                  do j = 1,spccnt
		     if( SUM( clsmap(j,:5,1) ) == 0 ) then
			clscnt(kpar) = clscnt(kpar) + 1
                        clsmap(j,kpar,1) = clscnt(kpar)
                        clsmap(clscnt(kpar),kpar,2) = j
		     end if
                  end do
                  cycle
               end if
               call GETTOKENS( buff, &
                               nchar, &
                               ',', &
                               symlen, &
                               tokens, &
                               toklen, &
                               20, &
                               no_tokens )
               if( no_tokens == 0 ) then
                  call ERRMES( ' Species input line in error@', lout, buff, 1, ' ' )
               end if

Tok_loop:      do j = 1,no_tokens
                  do l = 1,spccnt
                     if( trim(tokens(j)) == trim(spcsym(l)) ) then
                        clscnt(kpar) = clscnt(kpar) + 1
                        if( clscnt(kpar) > spccnt ) then
                           call ERRMES( ' Species count exceeds limit@', &
                                        lout, &
                                        buff, 1, buff )
                        end if
			if( SUM( clsmap(l,:5,1) ) /= 0 ) then
                           call ERRMES( ' # in two or more classes@', &
                                        lout, &
                                        tokens(j), &
                                        toklen(j), &
                                        buff )
		        end if
                        clsmap(l,kpar,1)            = clscnt(kpar)
                        clsmap(clscnt(kpar),kpar,2) = l
                        cycle tok_loop
                     end if
                  end do
                  call ERRMES( ' Class member # not in solution list@', &
                               lout, &
                               tokens(j), &
                               toklen(j), &
                               buff )
               end do Tok_loop
	    else
	       exit
            end if
	 end do Methods
      end do

      end subroutine SOL_CLS
