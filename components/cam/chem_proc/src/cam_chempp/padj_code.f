
      subroutine make_padj( fixmap, fixcnt, phtcnt, model, march )
!-----------------------------------------------------------------------
!        ... Write the photorate adjustment code
!-----------------------------------------------------------------------

      use var_mod, only : var_lim
      use io,      only : temp_path

      implicit none

!-----------------------------------------------------------------------
!        ... The arguments
!-----------------------------------------------------------------------
      integer, intent(in) ::  fixcnt
      integer, intent(in) ::  phtcnt
      integer, intent(in) ::  fixmap(var_lim,2)
      character(len=*), intent(in) ::  model
      character(len=*), intent(in) ::  march

!-----------------------------------------------------------------------
!        ... The local variables
!-----------------------------------------------------------------------
      integer  ::   k, rxno
      character(len=72) :: line
      logical  ::  first
      logical  ::  lexist

      
      inquire( file = trim( temp_path ) // 'mo_phtadj.F', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'mo_phtadj.F' )
      end if
      open( unit = 30, file = trim( temp_path ) // 'mo_phtadj.F' )

      line = ' '
      write(30,100) trim(line)
      line(7:) = 'module mo_phtadj'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'private'
      write(30,100) trim(line)
      line(7:) = 'public :: phtadj'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'contains'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      select case( model )
         case( 'MOZART' )
            line(7:) = 'subroutine phtadj( p_rate, inv, m, plnplv )'
         case ( 'CAM' )
            if( march /= 'VECTOR' ) then
               line(7:) = 'subroutine phtadj( p_rate, inv, m, ncol )'
            else
               line(7:) = 'subroutine phtadj( p_rate, inv, m, chnkpnts )'
            end if
         case ( 'WRF' )
            line(7:) = 'subroutine phtadj( p_rate, inv, m, n )'
      end select
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( model /= 'WRF' ) then
         line(7:) = 'use chem_mods,    only : nfs, phtcnt'
         write(30,100) trim(line)
      end if
      if( model == 'CAM' ) then
         line(7:) = 'use shr_kind_mod, only : r8 => shr_kind_r8'
         write(30,100) trim(line)
         if( march /= 'VECTOR' ) then
            line(7:) = 'use ppgrid,       only : pver'
            write(30,100) trim(line)
         end if
      end if
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'implicit none'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      line = '!       ... dummy arguments'
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      select case( model )
         case( 'MOZART' )
            line = '      integer, intent(in) :: plnplv'
         case ( 'CAM' )
            if( march /= 'VECTOR' ) then
               line = '      integer, intent(in)     :: ncol'
            else
               line = '      integer, intent(in) :: chnkpnts'
            end if
         case ( 'WRF' )
            line = '      integer, intent(in) :: n'
      end select
      write(30,100) trim(line)
      select case( model )
      case( 'MOZART' )
         line = '      real, intent(in)    :: inv(plnplv,nfs)'
         write(30,100) trim(line)
         line = '      real, intent(in)    :: m(plnplv)'
         write(30,100) trim(line)
         line = '      real, intent(inout) :: p_rate(plnplv,phtcnt)'
      case( 'CAM' )
         if( march /= 'VECTOR' ) then
            line = '      real(r8), intent(in)    :: inv(:,:,:)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    :: m(:,:)'
            write(30,100) trim(line)
            line = '      real(r8), intent(inout) :: p_rate(:,:,:)'
         else
            line = '      real(r8), intent(in)    :: inv(chnkpnts,max(1,nfs))'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    :: m(chnkpnts)'
            write(30,100) trim(line)
            line = '      real(r8), intent(inout) :: p_rate(chnkpnts,max(1,phtcnt))'
         end if
      case( 'WRF' )
         line = '      real, intent(in)    :: inv(:,:)'
         write(30,100) trim(line)
         line = '      real, intent(in)    :: m(:)'
         write(30,100) trim(line)
         line = '      real, intent(inout) :: p_rate(:,:)'
      end select
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      line = '!       ... local variables'
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      select case( model )
      case( 'MOZART' )
         line = '      real    ::  im(plnplv)'
      case( 'CAM' )
         if( march /= 'VECTOR' ) then
            line = '      integer  ::  k'
            write(30,100) trim(line)
            line = '      real(r8) ::  im(ncol)'
         else
            line = '      real(r8) ::  im(chnkpnts)'
         end if
      case( 'WRF' )
         line = '      real ::  im(n)'
      end select
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)

      if( model == 'CAM' .and. march /= 'VECTOR' ) then
         line = '      do k = 1,pver'
         write(30,100) trim(line)
      end if
      
      first = .true.
      do k = 1,fixcnt
         rxno = abs( fixmap(k,1) )
         if( fixmap(k,1) < 0 .and. rxno <= phtcnt ) then
            if( first ) then
               select case( model )
               case( 'CAM' )
                  if( march /= 'VECTOR' ) then
                     line(7:) = '   im(:ncol) = 1._r8 / m(:ncol,k)'
                  else
                     line(7:) = 'im(:) = 1._r8 / m(:)'
                  end if
               case default
                  line(7:) = 'im(:) = 1. / m(:)'
               end select
               write(30,100) trim(line)
               line = ' '
               first = .false.
            end if
            select case( model )
            case( 'CAM' )
               if( march /= 'VECTOR' ) then
                  write(line(7:),'(''   p_rate(:,k,'',i3,'') = p_rate(:,k,'',i3,'')'')') rxno,rxno
                  line(len_trim(line)+2:) = ' * inv(:,k,'
               else
                  line(7:) = 'p_rate(:,   ) = p_rate(:,   )'
                  write(line(16:18),'(i3)') rxno
                  write(line(32:34),'(i3)') rxno
                  line(len_trim(line)+2:) = ' * inv(:,'
               end if
               write(line(len_trim(line)+1:),'(i2)') fixmap(k,2)
               line(len_trim(line)+1:) = ') * im(:)'
            case default
               line(7:) = 'p_rate(:,   ) = p_rate(:,   )'
               write(line(16:18),'(i3)') rxno
               write(line(32:34),'(i3)') rxno
               line(len_trim(line)+2:) = ' * inv(:,'
               write(line(len_trim(line)+1:),'(i2)') fixmap(k,2)
               line(len_trim(line)+1:) = ') * im(:)'
            end select
            write(30,100) trim(line)
         end if
      end do

      if( model == 'CAM' .and. march /= 'VECTOR' ) then
         line = '      end do'
         write(30,100) trim(line)
      end if

      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end subroutine phtadj'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line(7:) = 'end module mo_phtadj'
      write(30,100) trim(line)
      
      close(30)
      
100   format(a)      
      
      end subroutine make_padj
