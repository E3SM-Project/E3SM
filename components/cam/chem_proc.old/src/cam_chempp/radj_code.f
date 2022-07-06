
      subroutine make_radj( fixmap, fixcnt, rxmap, rxmcnt, phtcnt, &
                            model, march )
!-----------------------------------------------------------------------
!        ... Make the reaction rate "adjustment" code
!-----------------------------------------------------------------------

      use var_mod, only : var_lim
      use rxt_mod, only : rxt_lim
      use io,      only : temp_path

      implicit none

!-----------------------------------------------------------------------
!        ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::  phtcnt
      integer, intent(in) ::  rxmcnt
      integer, intent(in) ::  fixcnt(2)
      integer, intent(in) ::  fixmap(var_lim,3,2)
      integer, intent(in) ::  rxmap(rxt_lim)
      character(len=*), intent(in) ::  model
      character(len=*), intent(in) ::  march

!-----------------------------------------------------------------------
!        ... local variables
!-----------------------------------------------------------------------
      integer  ::   j, k, l, rxno
      character(len=128):: line
      character(len=4)  :: dec_suffix
      logical  ::  first
      logical  ::  divide
      logical  ::  doloop
      logical  ::  lexist
      
      integer  ::  strlen
      
      inquire( file = trim( temp_path ) // 'mo_adjrxt.F', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'mo_adjrxt.F' )
      end if
      open( unit = 30, file = trim( temp_path ) // 'mo_adjrxt.F' )

      if( model == 'CAM' ) then
         dec_suffix = '(r8)'
      else
         dec_suffix = ' '
      end if

      line = ' '
      write(30,100) trim(line)
      line(7:) = 'module mo_adjrxt'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'private'
      write(30,100) trim(line)
      line(7:) = 'public :: adjrxt'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'contains'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( model == 'MOZART' ) then
         line(7:) = 'subroutine adjrxt( rate, inv, m, plnplv )'
      else if( trim(model) == 'CAM' ) then
         if( march /= 'VECTOR' ) then
            line(7:) = 'subroutine adjrxt( rate, inv, m, ncol )'
         else
            line(7:) = 'subroutine adjrxt( rate, inv, m, chnkpnts )'
         end if
      else if( trim(model) == 'WRF' ) then
         line(7:) = 'subroutine adjrxt( rate, inv, m )'
      end if
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( trim(model) == 'CAM' ) then
         if( trim(march) /= 'VECTOR' ) then
            line(7:) = 'use ppgrid, only : pver'
            write(30,100) trim(line)
         end if
         line(7:) = 'use shr_kind_mod, only : r8 => shr_kind_r8'
         write(30,100) trim(line)
         line(7:) = 'use chem_mods, only : nfs, rxntot'
         write(30,100) trim(line)
      end if
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'implicit none '
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      line = '!       ... dummy arguments'
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      if( trim(model) == 'MOZART' ) then
         line = '      integer, intent(in) :: plnplv'
         write(30,100) trim(line)
         line = '      real, intent(in)    :: inv(plnplv,nfs)'
         write(30,100) trim(line)
         line = '      real, intent(in)    :: m(plnplv)'
         write(30,100) trim(line)
         line = '      real, intent(inout) :: rate(plnplv,rxntot)'
      else if( trim(model) == 'CAM' ) then
         if( march /= 'VECTOR' ) then
            line = '      integer, intent(in) :: ncol'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    :: inv(ncol,pver,nfs)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    :: m(ncol,pver)'
            write(30,100) trim(line)
            line = '      real(r8), intent(inout) :: rate(ncol,pver,rxntot)'
         else
            line = '      integer, intent(in)     :: chnkpnts'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    :: inv(chnkpnts,max(1,nfs))'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    :: m(chnkpnts)'
            write(30,100) trim(line)
            line = '      real(r8), intent(inout) :: rate(chnkpnts,max(1,rxntot))'
         end if
      else if( trim(model) == 'WRF' ) then
         line = '      real, intent(in)    :: inv(:,:)'
         write(30,100) trim(line)
         line = '      real, intent(in)    :: m(:)'
         write(30,100) trim(line)
         line = '      real, intent(inout) :: rate(:,:)'
      end if
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      line = '!       ... local variables'
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      if( trim(model) == 'MOZART' ) then
         line = '      real    ::  im(plnplv)'
      else if( trim(model) == 'CAM' ) then
         if( march /= 'VECTOR' ) then
            line = '      real(r8) :: im(ncol,pver)'
         else
            line = '      real(r8) :: im(chnkpnts)'
         end if
      end if
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      write(30,100) trim(line)
      
      first  = .true.
      divide = .false.
      doloop = .false.
!---------------------------------------------------------
!	... First check reactions with invariants for 
!	    potential modification
!---------------------------------------------------------
      do j = 1,2
         do k = 1,fixcnt(j)
            rxno = abs( fixmap(k,1,j) )
            if( j == 2 .or. rxno > phtcnt ) then
               if( first ) then
                  first  = .false.
                  doloop = .true.
                  if( trim(model) == 'CAM' ) then
                     line = ' '
                  end if
               end if
               if( .not. divide .and. fixmap(k,1,j) < 0 ) then
                  if( trim(model) /= 'CAM' ) then
                     line(10:) = 'im(:) = 1. / m(:)'
                  else if( trim(march) /= 'VECTOR' ) then
                     line(10:) = 'im(:,:) = 1._r8 / m(:,:)'
                  else
                     line(10:) = 'im(k) = 1._r8 / m(k)'
                  end if
                  write(30,100) trim(line)
                  divide = .true.
               end if
               if( trim(model) /= 'CAM' ) then
                  line(10:) = 'rate(:,'
                  write(line(len_trim(line)+1:),'(i3)') rxno 
                  line(len_trim(line)+1:) = ') = rate(:,'
                  write(line(len_trim(line)+1:),'(i3)') rxno 
                  line(len_trim(line)+1:) = ')'
                  do l = 2,j+1
                     line(strlen(line)+1:) = ' * inv(:,'
                     write(line(strlen(line)+1:),'(i2)') fixmap(k,l,j)
                     line(strlen(line)+1:) = ')'
                  end do
               else
                  if( trim(march) /= 'VECTOR' ) then
                     line(10:) = 'rate(:,:,'
                     write(line(len_trim(line)+1:),'(i3)') rxno 
                     line(len_trim(line)+1:) = ') = rate(:,:,'
                     write(line(len_trim(line)+1:),'(i3)') rxno 
                     line(len_trim(line)+1:) = ')'
                     do l = 2,j+1
                        line(strlen(line)+1:) = ' * inv(:,:,'
                        write(line(strlen(line)+1:),'(i2)') fixmap(k,l,j)
                        line(strlen(line)+1:) = ')'
                     end do
                  else
                     line(10:) = 'rate(k,'
                     write(line(len_trim(line)+1:),'(i3)') rxno 
                     line(len_trim(line)+1:) = ') = rate(k,'
                     write(line(len_trim(line)+1:),'(i3)') rxno 
                     line(len_trim(line)+1:) = ')'
                     do l = 2,j+1
                        line(strlen(line)+1:) = ' * inv(k,'
                        write(line(strlen(line)+1:),'(i2)') fixmap(k,l,j)
                        line(strlen(line)+1:) = ')'
                     end do
                  end if
               end if
               if( fixmap(k,1,j) < 0 ) then
                  if( trim(model) /= 'CAM' ) then
                     line(strlen(line)+1:) = ' * im(:)'
                  else
                     if( trim(march) /= 'VECTOR' ) then
                        line(strlen(line)+1:) = ' * im(:,:)'
                     else
                        line(strlen(line)+1:) = ' * im(k)'
                     end if
                  end if
               end if
               write(30,100) trim(line)
            end if
         end do
      end do

!---------------------------------------------------------
!	... Now do all nonlinear reactions
!---------------------------------------------------------
      first = .true.
      do k = 1,rxmcnt
         rxno = rxmap(k)
         line = ' '
         if( model /= 'CAM' ) then
            line(10:) = 'rate(:,'
            write(line(len_trim(line)+1:),'(i3)') rxno 
            line(len_trim(line)+1:) = ') = rate(:,'
            write(line(len_trim(line)+1:),'(i3)') rxno 
            line(len_trim(line)+1:) = ') * m(:)'
         else
            if( march /= 'VECTOR' ) then
               line(10:) = 'rate(:,:,'
               write(line(len_trim(line)+1:),'(i3)') rxno 
               line(len_trim(line)+1:) = ') = rate(:,:,'
               write(line(len_trim(line)+1:),'(i3)') rxno 
               line(len_trim(line)+1:) = ') * m(:,:)'
            else
               line(10:) = 'rate(k,'
               write(line(len_trim(line)+1:),'(i3)') rxno 
               line(len_trim(line)+1:) = ') = rate(k,'
               write(line(len_trim(line)+1:),'(i3)') rxno 
               line(len_trim(line)+1:) = ') * m(k)'
            end if
         end if
         write(30,100) trim(line)
      end do

      if( model == 'CAM' ) then
         if ( rxmcnt > 0 .or. fixcnt(1) > 0 .or. fixcnt(2) > 0 ) then
         endif
      end if
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end subroutine adjrxt'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end module mo_adjrxt'
      write(30,100) trim(line)
      
      close(30)
      
100   format(a)      
      
      end subroutine make_radj
