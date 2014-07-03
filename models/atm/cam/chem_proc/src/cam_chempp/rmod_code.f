
      subroutine make_rmod( rxt2rel_pntr, rel_rxt_map, rxt2grp_pntr, &
                            grp_rxt_map, hetmap, hetcnt, rxntot, model, march )
!-----------------------------------------------------------------------
!        ... Make the group ratios reaction rate adjustment code
!-----------------------------------------------------------------------

      use rxt_mod, only : rxt_lim
      use io,      only : temp_path

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::    rxt2rel_pntr(rxt_lim,2)
      integer, intent(in) ::    rel_rxt_map(rxt_lim,3,2)
      integer, intent(in) ::    rxt2grp_pntr(rxt_lim,2)
      integer, intent(in) ::    grp_rxt_map(rxt_lim,3,2)
      integer, intent(in) ::    hetmap(rxt_lim)
      integer, intent(in) ::    rxntot
      integer, intent(in) ::    hetcnt
      character(len=*), intent(in) ::    model
      character(len=*), intent(in) ::    march

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: max_len = 90
      integer  ::   k, l, rxno, row, index
      character(len=max_len) :: line
      logical  ::  first
      logical  ::  found
      logical  ::  lexist
      
      integer  ::  strlen
      
      inquire( file = trim( temp_path ) // 'mo_rxtmod.F', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'mo_rxt_mod.F' )
      end if
      open( unit = 30, file = trim( temp_path ) // 'mo_rxt_mod.F' )

      line = ' '
      write(30,100) trim(line)
      line(7:) = 'module mo_rxt_mod'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'private'
      write(30,100) trim(line)
      line(7:) = 'public :: rxt_mod'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'contains'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      select case( model )
         case( 'MOZART' )
            line(7:) = 'subroutine rxt_mod( rate, het_rates, grp_ratios, plnplv )'
         case( 'CAM' )
            line(7:) = 'subroutine rxt_mod( rate, het_rates, grp_ratios, chnkpnts )'
         case( 'WRF' )
            line(7:) = 'subroutine rxt_mod( rate, grp_ratios )'
      end select
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim( line )
      if( model /= 'WRF' ) then
         line(7:) = 'use chem_mods, only : rxntot, hetcnt, grpcnt'
         write(30,100) trim( line )
      end if
      if( model == 'CAM' ) then
         line(7:) = 'use shr_kind_mod, only : r8 => shr_kind_r8'
         write(30,100) trim( line )
      end if
      line = ' '
      write(30,100) trim( line )
      line(7:) = 'implicit none '
      write(30,100) trim( line )
      line = ' '
      write(30,100) trim( line )
      line = '!---------------------------------------------------------------------------'
      write(30,100) trim(line)
      line = '!       ... dummy arguments'
      write(30,100) trim(line)
      line = '!---------------------------------------------------------------------------'
      write(30,100) trim(line)
      if( model == 'CAM' .and. march == 'VECTOR' ) then
         line = '      integer, intent(in) ::  chnkpnts'
         write(30,100) trim(line)
      else if( model /= 'WRF' ) then
         line = '      integer, intent(in) ::  plnplv'
         write(30,100) trim(line)
      end if
      if( model == 'MOZART' ) then
         line = '      real, intent(inout) ::  rate(plnplv,rxntot)'
         write(30,100) trim(line)
         line = '      real, intent(inout) ::  het_rates(plnplv,hetcnt)'
         write(30,100) trim(line)
         line = '      real, intent(in)    ::  grp_ratios(plnplv,grpcnt)'
      else if( model == 'CAM' ) then
         if( march /= 'VECTOR' ) then
            line = '      real(r8), intent(inout) ::  rate(:,:)'
            write(30,100) trim(line)
            line = '      real(r8), intent(inout) ::  het_rates(:,:)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  grp_ratios(:,:)'
         else
            line = '      real(r8), intent(inout) ::  rate(chnkpnts,max(1,rxntot))'
            write(30,100) trim(line)
            line = '      real(r8), intent(inout) ::  het_rates(chnkpnts,max(1,hetcnt))'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  grp_ratios(chnkpnts,max(1,grpcnt))'
         end if
      else if( model == 'WRF' ) then
         line = '      real, intent(in)    ::  grp_ratios(:,:)'
         write(30,100) trim(line)
         line = '      real, intent(inout) ::  rate(:,:)'
      end if
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      
      first = .true.
      do k = 1,rxntot
!-----------------------------------------------------------------------
!        ... Scan the group map
!-----------------------------------------------------------------------
	 found = .false.
         index = rxt2grp_pntr(k,1)
         row   = rxt2grp_pntr(k,2)
         do l = 1,index
            found = .true.
            if( first ) then
               line = ' '
               first = .false.
            end if
            rxno = grp_rxt_map(row,1,index)
            if( l == 1 ) then
	       line = ' '
               line(7:) = 'rate(:,   ) = rate(:,   )'
               write(line(14:16),'(i3)') rxno
               write(line(28:30),'(i3)') rxno
               line(strlen(line)+2:) = ' * grp_ratios(:,'
               write(line(strlen(line)+1:),'(i2)') grp_rxt_map(row,l+1,index)
               line(strlen(line)+1:) = ')'
	    else
	       line(len_trim(line)+1:) = ' &'
	       write(30,100) trim(line)
	       line(6:) = ' '
               line(33:) = ' * grp_ratios(:,'
               write(line(strlen(line)+1:),'(i2)') grp_rxt_map(row,l+1,index)
               line(strlen(line)+1:) = ')'
            end if
         end do
         if( found ) then
            write(30,100) trim(line)
         end if
      end do

      do k = 1,hetcnt
	 if( hetmap(k) /= 0 ) then
            line = ' '
            if( first ) then
               first = .false.
            end if
            line(7:) = '    het_rates(j,   ) = het_rates(j,   )'
            write(line(19:21),'(i3)') k
            write(line(38:40),'(i3)') k
            line(strlen(line)+2:) = ' * grp_ratios(:,'
            write(line(strlen(line)+1:),'(i2)') hetmap(k)
            line(strlen(line)+1:) = ')'
            write(30,100) trim(line)
         end if
      end do

      line = ' '
      write(30,100) trim(line)
      line = '      end subroutine rxt_mod'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end module mo_rxt_mod'
      write(30,100) trim(line)
      
      close(30)
      
100   format(a)      
      
      end subroutine make_rmod
