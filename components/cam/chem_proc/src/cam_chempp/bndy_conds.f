
      module mo_bndy_conds

      private
      public :: bndy_conds

      contains

      subroutine bndy_conds( lin, lout, new_nq, new_solsym, bc_is_fixed, bc_cnt )

      use var_mod, only : var_lim

      implicit none

!-----------------------------------------------------------------------
!       ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)    ::  lin, &              ! input unit number
                                 lout, &             ! output unit number
                                 new_nq              ! species count
      integer, intent(inout) ::  bc_cnt(:)           ! count of species with fixed bc
      logical, intent(inout) ::  bc_is_fixed(:,:)    ! fixed bndy condition matrix
      character(len=16), intent(in) ::  new_solsym(:) ! species names

!-----------------------------------------------------------------------
!       ... Local variables
!-----------------------------------------------------------------------
      integer  ::  nchar
      integer  ::  toklen(20)
      integer  ::  j, k
      integer  ::  no_tokens
      integer  ::  parsw(2) = 0
      integer  ::  bndy

      character(len=320) ::  buff
      character(len=320) ::  buffh
      character(len=16)  ::  tokens(20)

      logical  ::  found

      integer, parameter ::  symlen = 8

!-----------------------------------------------------------------------
!       ... Read the species boundary conditions
!-----------------------------------------------------------------------
section_loop : &
      do
         call cardin( lin, buff, nchar )
         buffh = buff
         call upcase( buffh )
         if( buffh == 'ENDBNDYCONDS' ) then
	    exit
	 end if
         if( buffh == 'FIXEDLOWERBC' ) then
            bndy = 1
            if( parsw(bndy) /= 0 ) then
               call errmes( ' BNDY_COND: Fixed Lower BC already specified@', lout, buff, 1, buff )
            end if
         else if( buffh == 'FIXEDUPPERBC' ) then
            bndy  = 2
            if( parsw(bndy) /= 0 ) then
               call errmes( ' BNDY_COND: Fixed Upper BC already specified@', lout, buff, 1, buff )
            end if
	 else
            call errmes( ' BNDY_COND: # is an invalid keyword @', lout, buff, 1, buff )
	 end if
         parsw(bndy) = 1

bndy_loop : &
         do
            call cardin( lin, buff, nchar )
            buffh = buff
            call upcase( buffh )
            if( buffh == 'ENDFIXEDLOWERBC' ) then
               if( bndy /= 1 ) then
                  call errmes( ' BNDY_COND: In Fixed Upper BC @', lout, buff, 1, buff )
               else if( parsw(bndy) /= 1 ) then
                  call errmes( ' BNDY_COND: Fixed Lower BC not entered@', lout, buff, 1, buff )
               end if
               exit
            else if( buffh == 'ENDFIXEDUPPERBC' ) then
               if( bndy /= 2 ) then
                  call errmes( ' BNDY_COND: In Fixed Lower BC @', lout, buff, 1, buff )
               else if( parsw(bndy) /= 1 ) then
                  call errmes( ' BNDY_COND: Fixed Upper BC not entered@', lout, buff, 1, buff )
               end if
               exit
            end if
            call gettokens(  buff,     nchar,    ',',      symlen, &
                             tokens,   toklen,   20,       no_tokens )
            if( no_tokens == 0 ) then
               call errmes( ' BNDY_COND: Species input line in error@', lout, buff, 1, buff )
            end if
token_loop : &
            do j = 1,no_tokens
               bc_cnt(bndy) = bc_cnt(bndy) + 1
               if( bc_cnt(bndy) > var_lim ) then
                  call errmes( ' BNDY_COND: Species count exceeds limit@', lout, buff, 1, buff )
               end if
	       found = .false.
	       do k = 1,new_nq
	          if( tokens(j) == new_solsym(k) ) then
		     if( bc_is_fixed(k,bndy) ) then
                        call errmes( '# is already specified @', lout, tokens(j), toklen(j), buff )
		     end if
	             bc_is_fixed(k,bndy) = .true.
	             found = .true.
		     exit
	          end if
	       end do
	       if( .not. found ) then
                  call errmes( '# is not in solution species list@', lout, tokens(j), toklen(j), buff )
               end if
            end do token_loop
         end do bndy_loop
      end do section_loop
         
      end subroutine bndy_conds

      end module mo_bndy_conds
