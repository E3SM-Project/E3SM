
      module mo_hist_out

      private
      public :: hist_out

      contains

      subroutine hist_out( histout, longnames, hst_file_cnt )
!-----------------------------------------------------------------------
!	... Process all output history file controls
!-----------------------------------------------------------------------

      use io
      use var_mod, only : indexh2o, srf_flx_cnt, dvel_cnt, spccnt, &
                          histout_cnt, histout_map, user_hst_names, &
                          spcsym, class_prod_cnt, class_loss_cnt, &
                          clscnt, clsmap, hst_map_lim
      use rxt_mod, only : hetcnt, hetmap, usrcnt, usrmap, gascnt, phtcnt

      implicit none

!-----------------------------------------------------------------------
!	... Dummy args
!-----------------------------------------------------------------------
      integer, intent(out)              :: hst_file_cnt        ! count of history files
      character (len=64), intent(inout) :: histout(6)          ! hist tape outputs
      logical, intent(inout)            :: longnames

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: maxparms = 26
      integer, parameter :: inst = 1, avgr = 2
      integer, parameter :: singl = 1, multi = 2
      integer, parameter :: symlen = 8

      integer ::   class, kpar, nchar, k, kkk, j, khold, m
      integer ::   time_ind, level_ind
      integer ::   fileno = 1
      integer ::   retcod
      integer ::   kindex
      integer ::   tokcnt, cnt
      integer ::   parsw(maxparms,2)
      integer ::   toklen(20)
      integer ::   wrk(hst_map_lim)
      real    ::   time
      character (len=23) :: parkeyend(maxparms)
      character (len=20) :: parkey(maxparms), keywrd
      character (len=16)  :: tokens(20)
      character (len=16)  :: temp
      logical ::   found, doing_file

!-----------------------------------------------------------------------
!	... Function declarations
!-----------------------------------------------------------------------
      integer  ::  lenof
      integer  ::  inclist
      integer  ::  inilist
      integer  ::  xlate

      parkey(1) = 'RETENTIONTIME'
      parkey(2) = 'WRITEFREQUENCY'
      parkey(3) = 'STARTFILENUMBER'
      parkey(4) = 'DENSITY'
      parkey(5) = 'PASSWORD'

      parkey(6)  = 'TRANSPORTEDSPECIES'
      parkey(7)  = 'GROUPMEMBERS'
      parkey(8)  = 'SURFACEFLUX'
      parkey(9)  = 'DEPOSITIONVELOCITY'
      parkey(10) = 'TEMPERATURE'
      parkey(11) = 'WATERVAPOR'
      parkey(12) = 'SURFACEPRESSURE'
      parkey(13) = 'PHOTORATES'
      parkey(14) = 'REACTIONRATES'
      parkey(15) = 'WASHOUTRATES'
      parkey(16) = 'EXTERNALFORCING'

      parkey(17) = 'DEFAULTOUTPUTS'
      parkey(18) = 'USERDEFINED'
      parkey(19) = 'PRODUCTION'
      parkey(20) = 'LOSS'
      parkey(21) = 'MASSDIAGS'
      parkey(22) = 'DEPOSITIONFLUX'
      parkey(23) = 'WASHOUTFLUX'
      parkey(24) = 'LIFETIME'
      parkey(25) = 'PRINTFREQUENCY'
      parkey(26) = 'LONGNAMES'

      parkeyend(6)  = 'ENDTRANSPORTEDSPECIES'
      parkeyend(7)  = 'ENDGROUPMEMBERS'
      parkeyend(8)  = 'ENDSURFACEFLUX'
      parkeyend(9)  = 'ENDDEPOSITIONVELOCITY'
      parkeyend(10) = 'ENDTEMPERATURE'
      parkeyend(11) = 'ENDWATERVAPOR'
      parkeyend(12) = 'ENDSURFACEPRESSURE'
      parkeyend(13) = 'ENDPHOTORATES'
      parkeyend(14) = 'ENDREACTIONRATES'
      parkeyend(15) = 'ENDWASHOUTRATES'
      parkeyend(16) = 'ENDEXTERNALFORCING'

      parkeyend(17) = 'ENDDEFAULTOUTPUTS'
      parkeyend(18) = 'ENDUSERDEFINED'
      parkeyend(19) = 'ENDPRODUCTION'
      parkeyend(20) = 'ENDLOSS'
      parkeyend(21) = 'ENDMASSDIAGS'
      parkeyend(22) = 'ENDDEPOSITIONFLUX'
      parkeyend(23) = 'ENDWASHOUTFLUX'
      parkeyend(24) = 'ENDLIFETIME'

      parsw(:,:)         = 0
      histout_cnt(:,:,:) = 0

!-----------------------------------------------------------------------
!   	... Check for intial file header
!-----------------------------------------------------------------------
      call cardin( lin, buff, nchar )
      buffh = buff
      call upcase ( buffh )
      if( buffh /= 'FILE' ) then
         call errmes( ' FILE keyword must follow OUTPUTS keyword@', &
                       lout, &
                       buff, &
                       nchar, &
                       buff )
      end if
!-----------------------------------------------------------------------
!   	... Scan for valid option keyword
!-----------------------------------------------------------------------
      doing_file = .true.
read_loop : &
      do
         call cardin( lin, buff, nchar )
         buffh = buff
         call upcase ( buffh )
         if( buffh == 'FILE' ) then
            if( doing_file ) then
               call errmes( ' ENDFILE must follow FILE keyword@', &
                             lout, buff, nchar, buff )
            end if
            if( fileno >= 10 ) then
               call errmes( ' Only 10 output files allowed@', &
                             lout, buff, nchar, buff )
            end if
            fileno = fileno + 1
            doing_file = .true.
            parsw(:,:) = 0
            cycle read_loop
         else if( buffh == 'ENDFILE' ) then
!-----------------------------------------------------------------------
!   	... Recast mass diags into individual deltas and fluxes
!           Note : only allow for 25 diagnostics
!-----------------------------------------------------------------------
            if( .not. doing_file ) then
               call errmes( ' ENDFILE not preceeded by FILE keyword@', &
                             lout, buff, nchar, buff )
            end if
            do time_ind = inst,avgr
               if( histout_cnt(16,time_ind,fileno) /= 0 ) then
                  k = min( histout_cnt(16,time_ind,fileno),25 )
                  j = 8*k
                  histout_cnt(16,time_ind,fileno) = j
                  wrk(:k) = histout_map(:k,16,time_ind,fileno)
                  do j = 1,k
                     histout_map(8*(j-1)+1:8*j,16,time_ind,fileno) = wrk(j)
                  end do
               end if
            end do
            doing_file = .false.
            cycle read_loop
         else if( buffh == 'ENDOUTPUTS' ) then
            if( doing_file ) then
               call errmes( ' ENDFILE must follow FILE keyword@', &
                             lout, buff, nchar, buff )
            end if
            hst_file_cnt = fileno
            exit
         end if
         k = index( buffh(:nchar), '=' )
         if( k /= 0 ) then
            keywrd = buffh(:k-1)
         else
            keywrd = buffh(:nchar)
         end if
         found = .false.
         do kpar = 1,maxparms
            if( keywrd == parkey(kpar) ) then
               found = .true.
               exit
            end if
         end do
         if( .not. found ) then
!-----------------------------------------------------------------------
!  	... Invalid parameter keyword; terminate the program
!-----------------------------------------------------------------------
            call errmes ( ' # is an invalid job control parameter keyword@', &
                          lout, keywrd, lenof(20,keywrd), buffh )
         end if
!-----------------------------------------------------------------------
!     	... Check for instantaneous or averaged qualifier
!-----------------------------------------------------------------------
output_type_qualifier : &
         if( kpar > 5 .and. kpar < maxparms-1 ) then
assignment_test : &
	    if( k /= 0 ) then
	       if( kpar == 18 ) then
	          j = index( buffh(:nchar), ',' )
		  if( j == 0 ) then
		     j = nchar
		  else
		     j = j - 1
		  end if
	       else
		 j = nchar
	       end if
	       if( buffh(k+1:j) == 'INST' ) then
	          time_ind = inst
	       else if( buffh(k+1:j) == 'AVRG' ) then
	          time_ind = avgr
	       else
                  call errmes( '0 *** # is invalid time qualifier@', &
                               lout, buff(k+1:j), j-k, buff )
	       end if
	       if( kpar == 18 ) then
	          k = index( buffh(:nchar), ',' )
	          if( k /= 0 ) then
	             if( buffh(k+1:nchar) == 'SINGLE' ) then
	                level_ind = singl
	             else if( buffh(k+1:nchar) == 'MULTI' ) then
	                level_ind = multi
	             else
                        call errmes( '0 *** # is invalid level qualifier@', &
                                     lout, buff(k+1:), nchar-k, buff )
                     end if
                  else
                     level_ind = multi
                  end if
               end if
            else
               if( kpar == 18 ) then
                  call errmes( ' User defined type must have time qualifier@', &
                               lout, buff, nchar, buff )
	       end if
	       time_ind = inst
	    end if assignment_test
	 else
	    time_ind = inst
	 end if output_type_qualifier
!-----------------------------------------------------------------------
!     	... Valid parameter keyword; now check for duplicate keyword
!-----------------------------------------------------------------------
	 if( kpar == 18 ) then
            if( parsw(kpar,time_ind) == level_ind ) then
               call errmes( '0 *** # has already been specified@', &
                             lout, parkey(kpar), k, ' ' )
            end if
	    parsw(kpar,time_ind) = level_ind
	 else
            if( parsw(kpar,time_ind) /= 0 ) then
               call errmes( '0 *** # has already been specified@', &
                             lout, parkey(kpar), k, ' ' )
            end if
	    parsw(kpar,time_ind) = 1
         end if

!-----------------------------------------------------------------------
!     	... Set individual options
!-----------------------------------------------------------------------
	 if( kpar >= 6 ) then
output_variables : &
	    select case( kpar )
	    case( 17 )
!-----------------------------------------------------------------------
!     	... The "default" option
!-----------------------------------------------------------------------
	       histout_cnt(1:20,time_ind,fileno) = 0
	       histout_cnt(1,time_ind,fileno) = spccnt(6)
	       histout_cnt(2,time_ind,fileno) = spccnt(7)
	       histout_cnt(3,time_ind,fileno) = spccnt(6)
	       histout_cnt(4,time_ind,fileno) = spccnt(6)
	       if( time_ind == inst ) then
	          histout_cnt(5,inst,fileno) = 1
	          if( indexh2o /= 0 ) then
	             histout_cnt(6,inst,fileno) = 1
	          end if
	          histout_cnt(7,inst,fileno) = 1
	       end if
	       do k = 1,spccnt(6)
		  histout_map(k,1,time_ind,fileno) = k
		  histout_map(k,3,time_ind,fileno) = k
		  histout_map(k,4,time_ind,fileno) = k
	       end do
	       do k = 1,spccnt(7)
		  histout_map(k,2,time_ind,fileno) = k
	       end do
	    case(maxparms-1)
!-----------------------------------------------------------------------
!     	... The printout frequency option
!-----------------------------------------------------------------------
	       call timcon( buff(k+1:nchar), &
                            time, &
                            lout )
	       histout(6) = buff(k+1:nchar)
	    case( maxparms )
!-----------------------------------------------------------------------
!     	... The longnames flag
!-----------------------------------------------------------------------
	       longnames = .true.
	    case( 10:12 )
!-----------------------------------------------------------------------
!     	... The temp, water vapor, and surf press options
!-----------------------------------------------------------------------
	       if( time_ind == inst ) then
	          histout_cnt(kpar-5,inst,fileno) = 1
               end if
	    case default 
!-----------------------------------------------------------------------
!     	... All other options
!-----------------------------------------------------------------------
	       call cardin( lin, buff, nchar )
	       buffh = buff
	       call upcase( buffh )
	       khold = kpar
	       do
	          if( buffh == parkeyend(khold) ) then
		     exit
		  end if
		  kpar = khold
		  call GETTOKENS( buff, &
                                  nchar, &
                                  ',', &
                                  symlen, &
                                  tokens, &
                                  toklen, &
                                  20, &
                                  tokcnt )
		  if( tokcnt == 0 ) then
		     call errmes( ' Hist tape output list in error@', &
                                  lout, &
                                  buff, &
                                  1, &
                                  buff )
		  end if
		  if( kpar == 18 ) then
		     if( histout_cnt(11+level_ind,time_ind,fileno) + tokcnt > hst_map_lim ) then
		        call errmes( ' Hist tape output list > hst_map_lim elements@', &
                                     lout, &
                                     buff, &
                                     1, &
                                     buff )
		     end if
		  else if( histout_cnt(kpar-5,time_ind,fileno) + tokcnt > hst_map_lim ) then
		     call errmes( ' Hist tape output list > hst_map_lim elements@', &
                                  lout, &
                                  buff, &
                                  1, &
                                  buff )
		  end if
	          if( kpar == 18 ) then
	             do j = 1,tokcnt
	                cnt = histout_cnt(11+level_ind,time_ind,fileno) + 1
	                histout_cnt(11+level_ind,time_ind,fileno) = cnt
	                user_hst_names(cnt,2*(time_ind-1)+level_ind) = tokens(j)
	             end do
		  else if( kpar <= 9 .or. kpar >= 13 ) then
		     temp = tokens(1)
	             if( kpar /= 7 ) then
			kindex = 6
	             else
			kindex = 7
		     end if
	             kpar = kpar - 5
		     call upcase( temp )
!-----------------------------------------------------------------------
!     	... Handle the "all" list specifier
!-----------------------------------------------------------------------
		     if( tokcnt == 1 .and. temp == 'ALL' ) then
			if( kpar <= 4 ) then
			   do j = 1,spccnt(kindex)
			      histout_map(j,kpar,time_ind,fileno) = j
			   end do
			   histout_cnt(kpar,time_ind,fileno) = spccnt(kindex)
			else if( kpar >= 8 ) then
			   select case( kpar )
			      case( 8 )
			         do j = 1,phtcnt
			            histout_map(j,kpar,time_ind,fileno) = j
			         end do
			         histout_cnt(kpar,time_ind,fileno) = phtcnt
			      case( 9 )
			         do j = 1,gascnt
			            histout_map(j,kpar,time_ind,fileno) = j
			         end do
			         histout_cnt(kpar,time_ind,fileno) = gascnt
			      case( 10 )
			         do j = 1,hetcnt
			            histout_map(j,kpar,time_ind,fileno) = j
			         end do
			         histout_cnt(kpar,time_ind,fileno) = hetcnt
			      case( 11 )
			         do j = 1,usrcnt
			            histout_map(j,kpar,time_ind,fileno) = j
			         end do
			         histout_cnt(kpar,time_ind,fileno) = usrcnt
			   end select
			end if
!-----------------------------------------------------------------------
!     	... Handle individual list elements
!-----------------------------------------------------------------------
		     else
			do j = 1,tokcnt
			   if( kpar == 8 .or. kpar == 9 ) then
	                      call intcon( tokens(j), &      ! input string to convert
                                           toklen(j), &      ! length of input string
                                           k, &              ! surrogate for converted number
                                           retcod )          ! return code
	                      if( retcod /= 0 ) then
	                         call errmes( ' # is not a valid integer@', &
                                              lout, tokens(j), toklen(j), buff )
	                      end if
			      if( kpar == 8 .and. k > phtcnt) then
	                         call errmes( ' # out of photolysis rate numbering@', &
                                              lout, tokens(j), toklen(j), buff )
                              else if( k > gascnt ) then
                                 call errmes( ' # out of reaction rate numbering@', &
                                              lout, tokens(j), toklen(j), buff )
			      end if
			      histout_cnt(kpar,time_ind,fileno) = histout_cnt(kpar,time_ind,fileno) + 1
			      histout_map(histout_cnt(kpar,time_ind,fileno),kpar,time_ind,fileno) = k
			   else
			      class = 0
			      k = inclist( tokens(j), &
                                           spcsym(1,kindex), &
                                           spccnt(kindex) )
			      if( k == 0 ) then
			         call errmes( '# not in list@', &
                                              lout, &
                                              tokens(j), &
                                              toklen(j), &
                                              buff )
			      end if
			      if( kpar >= 10 ) then
			         if( kpar == 10 ) then
			            k = inilist( k, hetmap, hetcnt )
			         else if( kpar == 11 ) then
			            k = inilist( k, usrmap, usrcnt )
			         else if( kpar == 18 ) then
			            class = inilist( k, hetmap, hetcnt )
			         end if
			         if( k == 0 ) then
			            call errmes( '# not in list@', &
                                                 lout, &
                                                 tokens(j), &
                                                 toklen(j), &
                                                 buff )
			         end if
			         if( kpar == 14 .or. kpar == 15 .or. kpar == 19 ) then
				    kkk   = k
				    class = xlate( kkk )
                                    do m = 1,clscnt(class)
                                       if( clsmap(m,class,2) == k ) then
                                          exit
                                       end if
                                    end do
                                    k = m
                                    if( class /= 0 ) then
                                       if( kpar == 14 ) then
                                          class_prod_cnt(class,time_ind) = class_prod_cnt(class,time_ind) + 1
                                       else if( kpar == 15 ) then
                                          class_loss_cnt(class,time_ind) = class_loss_cnt(class,time_ind) + 1
                                       end if
                                    end if
                                 end if
                              end if
                              histout_cnt(kpar,time_ind,fileno) = histout_cnt(kpar,time_ind,fileno) + 1
                              histout_map(histout_cnt(kpar,time_ind,fileno),kpar,time_ind,fileno) = 1000*class + k
                           end if
                        end do
                     end if
                  end if
                  call cardin( lin, buff, nchar )
                  buffh = buff
                  call upcase( buffh )
               end do
            end select output_variables
         else if( kpar <= 2 ) then
            call timcon( buff(k+1:nchar), time, lout )
            histout(kpar) = buff(k+1:nchar)
         else if( kpar == 3 ) then
            call intcon( buff(k+1:nchar), & ! input string to convert
                         nchar - k,       & ! length of input string
                         toklen(1),       & ! surrogate for converted number
                         toklen(2) )        ! surrogate for error code
            if( toklen(2) /= 0 ) then
               call errmes( ' # is not a valid integer@', &
                            lout, buff(k+1:nchar), nchar - k, buff )
            end if
            histout(kpar) = buff(k+1:nchar)
         else
            histout(kpar) = buff(k+1:nchar)
         end if
      end do read_loop

      end subroutine hist_out

      end module mo_hist_out
