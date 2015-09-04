      subroutine HIST_INP( lin, &
                           lout, &
                           histinp, &
			   dyn_hst_fld_cnt )

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)            :: lin
      integer, intent(in)            :: lout
      integer, intent(out)           :: dyn_hst_fld_cnt(2)
      character(len=64), intent(out) :: histinp(4)          ! hist tape inputs

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer  ::  kpar, nchar, k
      integer  ::  retcod, slen
      integer  ::  parsw(6)
      real     ::  time
      character(len=80) :: buff
      character(len=80) :: buffh
      character(len=20) :: parkey(6), keywrd
      logical  ::  found

      integer  ::  LENOF

      parkey(1) = 'DYNAMICSMSSFILE'
      parkey(2) = 'STARTTIME'
      parkey(3) = 'ICMSSFILE'
      parkey(4) = 'DYNHISTTAPE'
      parkey(5) = 'MULTILEVELFIELDS'
      parkey(6) = 'SINGLELEVELFIELDS'

      parsw  = 0
      dyn_hst_fld_cnt = -1

!-----------------------------------------------------------------------
!   	... Scan for valid option keyword
!-----------------------------------------------------------------------
      do
         call CARDIN( lin, buff, nchar )
	 buffh = buff
	 call UPCASE ( buffh )
         if( buffh == 'ENDINPUTS' ) then
	    if( dyn_hst_fld_cnt(1) == -1 .and. dyn_hst_fld_cnt(2) == -1 ) then
	       if( histinp(4) == 'LONG' ) then
		  dyn_hst_fld_cnt(1) = 57
		  dyn_hst_fld_cnt(2) = 44
	       else if( histinp(4) == 'SHORT' ) then
		  dyn_hst_fld_cnt(1) = 10
		  dyn_hst_fld_cnt(2) = 4
	       end if
	    end if
            exit
	 end if
	 k = INDEX( buffh(:nchar), '=' )
         if( k /= 0 ) then
	    keywrd = buffh(:k-1)
	    found = .false.
            do kpar = 1,6
               if( keywrd == parkey(kpar) ) then
		  found = .true.
		  exit
	       end if
	    end do
	    if( .not. found ) then
               call ERRMES ( ' # is an invalid job control' &
                          // ' parameter keyword@', &
                             lout, &
                             keywrd, &
                             LENOF(20,keywrd), &
                             buffh )
            end if
	 else
!-----------------------------------------------------------------------
!  	... Invalid parameter keyword; terminate the program
!-----------------------------------------------------------------------
            call ERRMES ( ' Job ctl specification has no = operator@', &
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
	 if( kpar == 2 ) then
	    call TIMCON( buff(k+1:nchar), time, lout )
	    histinp(2) = buff(k+1:nchar)
	 else if( kpar == 4 ) then
	    histinp(4) = buffh(k+1:nchar)
	 else if( kpar == 5 ) then
	    slen = LEN_TRIM( buff(k+1:nchar) )
	    call INTCON( buff(k+1:nchar), slen, dyn_hst_fld_cnt(1), retcod )
	    if( retcod /= 0 .or. dyn_hst_fld_cnt(1) < 0 ) then
               call ERRMES ( ' # is an invalid Dyn hst tape field count@', &
                             lout, &
                             buff(k+1:nchar), &
                             slen, &
                             buffh )
	    end if
	 else if( kpar == 6 ) then
	    slen = LEN_TRIM( buff(k+1:nchar) )
	    call INTCON( buff(k+1:nchar), slen, dyn_hst_fld_cnt(2), retcod )
	    if( retcod /= 0 .or. dyn_hst_fld_cnt(2) < 0 ) then
               call ERRMES ( ' Dyn hst tape has invalid field count@', &
                             lout, &
                             buff(k+1:nchar), &
                             slen, &
                             buffh )
	    end if
	 else
	    histinp(kpar) = buff(k+1:nchar)
	 end if
	 parsw(kpar) = 1
      end do

      end subroutine HIST_INP
