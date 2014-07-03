
      module MASS_DIAGS
!--------------------------------------------------------------------
!	... General purpose mozart2 diagnostic module
!--------------------------------------------------------------------

      implicit none

      integer, parameter :: max_diags = 100

      type CONSERVATION
	 integer, dimension(2) :: long_ind, lat_ind, lev_ind
	 real, dimension(2) :: longitudes, latitudes, levels
         character(len=16) :: species
      end type CONSERVATION

      integer :: ndiags = 0
      type(CONSERVATION) :: mdiags(max_diags)
      real :: bigneg

      CONTAINS

      subroutine INIDIAGS( name, plon, plev, plat )
!--------------------------------------------------------------------
!	... Initialize the diagnostic variables
!--------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      integer, intent(in) :: plon, plev, plat
      character(len=16), intent(in) :: name

      if( ndiags >= max_diags ) then
         write(*,*) ' INIDIAGS: Exceeded diagnostic limit'
	 stop
      end if
      ndiags = ndiags + 1
      mdiags(ndiags)%species = name
      mdiags(ndiags)%long_ind = (/ 1,plon /)
      mdiags(ndiags)%lat_ind = (/ 1,plat /)
      mdiags(ndiags)%lev_ind = (/ 1,plev /)
      bigneg = -HUGE( bigneg )
      mdiags(ndiags)%longitudes = bigneg
      mdiags(ndiags)%latitudes  = bigneg
      mdiags(ndiags)%levels     = bigneg

      end subroutine INIDIAGS

      subroutine MASS_DIAGNOSTICS( spcsym, spccnt, plon, plev, plat )
!--------------------------------------------------------------------
!	... Process the mass diagnostics
!--------------------------------------------------------------------

      use IO

      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      integer, intent(in) :: plon, plev, plat
      integer, intent(in) :: spccnt(:)
      character(len=16), intent(in) :: spcsym(:,:)

!--------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------
      integer, parameter :: symlen = 8
      integer, parameter :: numlen = 16
      integer :: tail, head
      integer :: retcod, nchar
      integer :: tokcnt, m, sep
      integer :: toklen(20)
      integer :: iind(2)
      real    :: rind(2)
      character(len=numlen) :: tokens(20)
      character(len=symlen) :: tstring
      logical :: do_tokens, do_grid

      do_tokens = .true.
      do_grid   = .false.
      tokcnt    = 0
      head      = 0
      do
         call CARDIN( lin, buff, nchar )
	 buffh = buff
	 call UPCASE( buffh )
	 if( buffh == 'ENDMASS_DIAGNOSTICS' ) then
	    return
	 end if
	 if( do_tokens ) then
	    call GETTOKENS( buff, &
	                    nchar, &
			    ',', &
			    symlen, &
			    tokens, &
			    toklen, &
			    20, &
			    tokcnt )
	    if( tokcnt == 0 ) then
	       call ERRMES( 'Diagnostic spieces list error@', &
	                    lout, &
			    buff, &
			    1, &
			    buff )
	    end if
	    tstring = tokens(1)
	    call UPCASE( tstring )
	    if( tokcnt == 1 .and. tstring == 'ALL' ) then
	       do m = 1,spccnt(6)
	          call INIDIAGS( spcsym(m,6), plon, plev, plat )
	       end do
	    else
	       do m = 1,tokcnt
	          call INIDIAGS( tokens(m), plon, plev, plat )
	       end do
	    end if
	    do_grid = .true.
	    do_tokens = .false.
	    tail = head + 1
            head = ndiags
	    cycle
	 else if( do_grid ) then
	    if( buffh /= 'GRID' ) then
	       call ERRMES( '"GRID" keyword missing@', &
	                    lout, &
			    buff, &
			    1, &
			    buff )
	    end if
	    do
               call CARDIN( lin, buff, nchar )
	       buffh = buff
	       call UPCASE( buffh )
	       if( buffh == 'ENDLST' ) then
	          do_grid = .false.
		  do_tokens = .true.
		  exit
	       end if
	       sep = INDEX( buffh, '=' )
	       if( sep == 0 .or. sep == nchar ) then
	          call ERRMES( 'Grid keyword missing "=" separator@', &
	                       lout, &
			       buff, &
			       1, &
			       buff )
	       end if
	       call GETTOKENS( buff(sep+1:), &
	                       nchar-sep, &
			       ',', &
			       numlen, &
			       tokens, &
			       toklen, &
			       20, &
			       tokcnt )
	       if( tokcnt == 0 ) then
	          call ERRMES( 'Grid value format invalid@', &
	                       lout, &
			       buff, &
			       1, &
			       buff )
	       end if
               tokcnt = MIN(tokcnt,2)
	       select case( buffh(:sep-1) )
	          case( 'LONGITUDES', 'LATITUDES', 'LEVELS' )
		     do m = 1,tokcnt
                        call RELCON( tokens(m), &
			             toklen(m), &
				     rind(m), &
				     retcod )
                        if( retcod /= 0 ) then
	                   call ERRMES( '# is an invalid real number@', &
			                lout, &
                                        tokens(m), &
					toklen(m), &
					buff )
	                end if
		     end do
	             select case( buffh(:sep-1) )
	                case( 'LONGITUDES' )
			   do m = tail,head
			      mdiags(m)%longitudes(:tokcnt) = rind(:tokcnt)
			   end do
	                case( 'LATITUDES' )
			   do m = tail,head
			      mdiags(m)%latitudes(:tokcnt) = rind(:tokcnt)
			   end do
	                case( 'LEVELS' )
			   do m = tail,head
			      mdiags(m)%levels(:tokcnt) = rind(:tokcnt)
			   end do
		     end select
	          case( 'LONG_INDEX', 'LAT_INDEX', 'LEV_INDEX' )
		     do m = 1,tokcnt
                        call INTCON( tokens(m), &
			             toklen(m), &
				     iind(m), &
				     retcod )
	                if( retcod /= 0 ) then
                           call ERRMES( ' # is an invalid integer@', &
                                        lout, &
                                        tokens(m), &
				        toklen(m), &
					buff )
	                end if
		     end do
	             select case( buffh(:sep-1) )
	                case( 'LONG_INDEX' )
			   do m = tail,head
			      mdiags(m)%long_ind(:tokcnt) = iind(:tokcnt)
			   end do
	                case( 'LAT_INDEX' )
			   do m = tail,head
			      mdiags(m)%lat_ind(:tokcnt) = iind(:tokcnt)
			   end do
	                case( 'LEV_INDEX' )
			   do m = tail,head
			      mdiags(m)%lev_ind(:tokcnt) = iind(:tokcnt)
			   end do
		     end select
		  case default
	             call ERRMES( 'Grid keyword invalid@', &
	                          lout, &
			          buff, &
			          1, &
			          buff )
	       end select
	    end do
         end if
      end do
      call CHECK_RANGE()

      end subroutine MASS_DIAGNOSTICS

      subroutine CHECK_RANGE( )
!--------------------------------------------------------------------
!	... Check spatial domain ranges
!--------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------
      integer :: k, m
      real    :: value
      real, parameter :: eps = 1.e-5

      do m = 1,ndiags
         do k = 1,2
	    value = mdiags(m)%longitudes(k)
	    if( value /= bigneg ) then
	       if( value < 0. ) then
	          value = MOD( value,360. ) + 360.
	       else
	          value = MOD( value+eps,360. )
	       end if
	       mdiags(m)%longitudes(k) = value
            end if
	    value = mdiags(m)%latitudes(k)
	    if( value /= bigneg ) then
	       if( ABS( value) > 90. ) then
	          value = SIGN( 90.-eps,value )
	       end if
	       mdiags(m)%latitudes(k) = value
            end if
	 end do
      end do

      end subroutine CHECK_RANGE

      subroutine MASS_DIAGS_SERIALIZE( unit )
!--------------------------------------------------------------------
!	... Write conservation info to output file
!--------------------------------------------------------------------

      implicit none

!--------------------------------------------------------------------
!	... Dummy args
!--------------------------------------------------------------------
      integer, intent(in) :: unit

!--------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------
      integer :: m

      write(unit,'(i3)') ndiags
      do m = 1,ndiags
         write(unit,*) mdiags(m)
      end do

      end subroutine MASS_DIAGS_SERIALIZE

      end module MASS_DIAGS
