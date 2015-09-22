
       subroutine srfflx( lin, &
                          lout, &
                          new_nq, &
                          new_solsym, &
                          srf_flx_map, &
                          srf_flx_cnt, &
			  tag )

      use var_mod, only : var_lim

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)  ::  lin, &       ! input unit number
                               lout, &      ! output unit number
                               new_nq, &    ! species count
                               tag          ! emission or deposition tag ( 1,2 )
      integer, intent(inout) ::  srf_flx_cnt  ! count of species with srf flux
      integer, intent(out) ::  srf_flx_map(*)  ! srf flux "map"
      character(len=16), intent(in) ::  new_solsym(*) ! species names

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer  ::  nchar
      integer  ::  toklen(20)
      integer  ::  j, k
      integer  ::  no_tokens

      character(len=320) ::  buff
      character(len=320) ::  buffh
      character(len=16)  ::  tokens(20)

      logical  ::  found

      integer, parameter ::  symlen = 8

!-----------------------------------------------------------------------
!       ... Read the surface flux species
!-----------------------------------------------------------------------
      do
         call cardin( lin, buff, nchar )
         buffh = buff
         call upcase( buffh )
         if( tag == 1 .and. buffh == 'ENDSURFACEFLUX' ) then
	    exit
         else if( tag == 2 .and. buffh == 'ENDSURFACEDEPOSITION' ) then
	    exit
	 end if
         call gettokens(  buff,     nchar,    ',',      symlen, &
                          tokens,   toklen,   20,       no_tokens )
         if( no_tokens == 0 ) then
            call errmes( ' SRFFLX: Species input line in error@', lout, buff, 1, buff )
         end if
         do j = 1,no_tokens
            srf_flx_cnt = srf_flx_cnt + 1
            if( srf_flx_cnt > var_lim ) then
               call errmes( ' SRFFLX: Species count exceeds limit@', lout, buff, 1, buff )
            end if
	    found = .false.
	    do k = 1,new_nq
	       if( tokens(j) == new_solsym(k) ) then
		  if( srf_flx_cnt > 1 ) then
		     if( any( srf_flx_map(:srf_flx_cnt) == k ) ) then
			if( tag == 1 ) then
                           call errmes( '# is already in srf emis list@', lout, tokens(j), toklen(j), buff )
			else if( tag == 2 ) then
                           call errmes( '# is already in dry dep list@', lout, tokens(j), toklen(j), buff )
			end if
		     end if
		  end if
	          srf_flx_map(srf_flx_cnt) = k
	          found = .true.
		  exit
	       end if
	    end do
	    if( .not. found ) then
               call errmes( '# is not in solution species list@', lout, tokens(j), toklen(j), buff )
            end if
         end do
      end do
         
      end subroutine srfflx
