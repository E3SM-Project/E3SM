      
      subroutine mak_grp_vmr( grp_mem_cnt, mem2grp_map, model, march )
!-------------------------------------------------------------------
!	... Write the group volume mixing ratios code
!-------------------------------------------------------------------

      use io, only : temp_path

      implicit none

!-------------------------------------------------------------------
!	... Dummy args
!-------------------------------------------------------------------
      integer, intent(in) ::  grp_mem_cnt
      integer, intent(in) ::  mem2grp_map(*)
      character(len=*), intent(in) ::  model
      character(len=*), intent(in) ::  march

!-------------------------------------------------------------------
!	... Local variables
!-------------------------------------------------------------------
      integer, parameter :: max_len= 90
      integer  ::  m
      character(len=max_len) :: line
      logical  ::  lexist

      inquire( file = trim( temp_path ) // 'mo_make_grp_vmr.F', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'mo_make_grp_vmr.F' )
      end if
      open( unit = 30, file = trim( temp_path ) // 'mo_make_grp_vmr.F' )

      line = ' '
      write(30,100) trim(line)
      line(7:) = 'module mo_make_grp_vmr'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'private'
      write(30,100) trim(line)
      line(7:) = 'public :: mak_grp_vmr'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'contains'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      select case( model )
         case( 'MOZART' )
            line(7:) = 'subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, plonl )'
         case( 'CAM' )
            if( march /= 'VECTOR' ) then
               line(7:) = 'subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, plonl )'
            else
               line(7:) = 'subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs, chnkpnts )'
            end if
         case( 'WRF' )
            line(7:) = 'subroutine mak_grp_vmr( vmr, group_ratios, group_vmrs )'
      end select
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( model == 'MOZART' ) then
         line(7:) = 'use mo_grid,   only : plev, pcnstm1'
         write(30,100) trim(line)
         line(7:) = 'use chem_mods, only : grpcnt'
      else if( model == 'CAM' ) then
         line(7:) = 'use chem_mods, only : grpcnt, gas_pcnst'
         write(30,100) trim(line)
         if( march /= 'VECTOR' ) then
            line(7:) = 'use ppgrid, only : pver'
            write(30,100) trim(line)
         end if
         line(7:) = 'use shr_kind_mod, only : r8 => shr_kind_r8'
      else if( model == 'WRF' ) then
         line(7:) = ' '
      end if
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'implicit none '
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '!----------------------------------------------------------------------------'
      write(30,100) trim(line)
      line = '!        ... dummy arguments'
      write(30,100) trim(line)
      line = '!----------------------------------------------------------------------------'
      write(30,100) trim(line)
      select case( model )
         case( 'MOZART' )
            line = '      integer, intent(in) :: plonl'
            write(30,100) trim(line)
         case( 'CAM' )
            if( march /= 'VECTOR' ) then
               line = '      integer, intent(in) :: plonl'
            else
               line = '      integer, intent(in) :: chnkpnts'
            end if
            write(30,100) trim(line)
      end select
      if( model == 'MOZART' ) then
         line = '      real, intent(in)    :: vmr(plonl,plev,pcnstm1)'
         write(30,100) trim(line)
         line = '      real, intent(in)    :: group_ratios(plonl,plev,grpcnt)'
         write(30,100) trim(line)
         line = '      real, intent(out)   :: group_vmrs(plonl,plev,grpcnt)'
      else if( model == 'CAM' ) then
         if( march /= 'VECTOR' ) then
            line = '      real(r8), intent(in)    :: vmr(:,:,:)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    :: group_ratios(:,:,:)'
            write(30,100) trim(line)
            line = '      real(r8), intent(out)   :: group_vmrs(:,:,:)'
         else
            line = '      real(r8), intent(in)    :: vmr(chnkpnts,max(1,gas_pcnst))'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    :: group_ratios(chnkpnts,max(1,grpcnt))'
            write(30,100) trim(line)
            line = '      real(r8), intent(out)   :: group_vmrs(chnkpnts,max(1,grpcnt))'
         end if
      else if( model == 'WRF' ) then
         line = '      real, intent(in)    :: vmr(:,:)'
         write(30,100) trim(line)
         line = '      real, intent(in)    :: group_ratios(:,:)'
         write(30,100) trim(line)
         line = '      real, intent(out)   :: group_vmrs(:,:)'
      end if
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( model /= 'WRF' ) then
         line = '!----------------------------------------------------------------------------'
         write(30,100) trim(line)
         line = '!        ... local variables'
         write(30,100) trim(line)
         line = '!----------------------------------------------------------------------------'
         write(30,100) trim(line)
         line = '      integer ::  k'
         write(30,100) trim(line)
      end if
      if( grp_mem_cnt > 0 ) then
         line = ' '
         write(30,100) trim(line)
         select case( model )
            case( 'MOZART' )
               line(7:) = 'do k = 1,plev'
               write(30,100) trim(line)
            case( 'CAM' )
               if( march /= 'VECTOR' ) then
                  line(7:) = 'do k = 1,plev'
               else
                  line(7:) = 'do k = 1,chnkpnts'
               end if
               write(30,100) trim(line)
         end select
         do m = 1,grp_mem_cnt
	    line = ' '
            select case( model )
               case( 'MOZART' )
	          line(10:) = 'group_vmrs(:,k,  ) = group_ratios(:,k,  )'
	          write(line(25:26),'(i2)') m
	          write(line(48:49),'(i2)') m
	          line(len_trim(line)+1:) = ' * vmr(:,k,'
               case( 'CAM' )
                  if( march /= 'VECTOR' ) then
	             line(10:) = 'group_vmrs(:,k,  ) = group_ratios(:,k,  )'
	             write(line(25:26),'(i2)') m
	             write(line(48:49),'(i2)') m
	             line(len_trim(line)+1:) = ' * vmr(:,k,'
                  else
	             line(10:) = 'group_vmrs(k,  ) = group_ratios(k,  )'
                     write(line(10:),'(''group_vmrs(k,'',i2,'') = group_ratios(k,'',i2,'')'')') m, m
                     line(len_trim(line)+1:) = ' * vmr(:,k,'
                  end if
               case( 'WRF' )
                  line(7:) = 'group_vmrs(:,'
                  write(line(len_trim(line)+1:),*) m
                  line(len_trim(line)+1:) = ') = group_ratios(:,'
                  write(line(len_trim(line)+1:),*) m
                  line(len_trim(line)+1:) = ') * vmr(:,'
            end select
            write(line(len_trim(line)+1:),'(i2,'')'')') mem2grp_map(m)
            write(30,100) trim(line)
         end do
         if( model /= 'WRF' ) then
            line = ' '
            line(7:) = 'end do'
            write(30,100) trim(line)
         end if
      end if
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end subroutine mak_grp_vmr'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end module mo_make_grp_vmr'
      write(30,100) trim(line)

100   format(a)

      end subroutine mak_grp_vmr
