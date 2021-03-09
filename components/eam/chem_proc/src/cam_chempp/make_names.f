
      subroutine make_name_mod
!--------------------------------------------------------------------------------
!	... Makes a module of parameter species names
!--------------------------------------------------------------------------------

      use var_mod, only : spc_cnt => new_nq, spc_names => new_solsym, &
			  grp_mem_cnt, grp_mem_names => grp_mem_sym
      use io, only : temp_path

      implicit none

!--------------------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------------------
      integer :: i, j
      integer :: beg, end
      character(len=80) :: buff
      character(len=63) :: legal = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' &
				// 'abcdefghijklmnopqrstuvwxyz' &
				// '0123456789_'
      character(len=16)  :: name
      logical :: lexist

      inquire( file = trim( temp_path ) // 'spc_names.mod', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'spc_names.mod' )
      end if
      open( unit = 30, &
	    file = trim( temp_path ) // 'spc_names.mod' )
             
      buff = ' '
      write(30,*) ' '
      buff(7:) = 'module m_spc_id'
      write(30,'(a)') buff
      write(30,*) ' '
      buff(7:) = 'implicit none'
      write(30,'(a)') buff
      buff = ' '
      write(30,*) ' '

      do i = 1,spc_cnt
	 name = spc_names(i)
	 end = len_trim(name)
	 beg = 1
	 do
	    j = VERIFY( name(beg:end), legal )
	    if( j == 0 ) then
	       exit
	    end if
	    j = j + beg - 1
	    if( j == end ) then
	       end = end - 1
	       exit
	    end if
	    name(j:j) = '_'
	    if( j >= end ) then
	       exit
	    end if
	    beg = j + 1
	 end do
	 write(buff(7:),'(''integer, parameter :: id_'',a,1x,''='',1x,i3)') &
	             name(:end), i
         write(30,'(a)') buff(:len_trim(buff))
      end do
      write(30,*) ' '
      do i = 1,grp_mem_cnt
	 name = grp_mem_names(i)
	 end = len_trim(name)
	 beg = 1
	 do
	    j = VERIFY( name(beg:end), legal )
	    if( j == 0 ) then
	       exit
	    end if
	    j = j + beg - 1
	    if( j == end ) then
	       end = end - 1
	       exit
	    end if
	    name(j:j) = '_'
	    if( j >= end ) then
	       exit
	    end if
	    beg = j + 1
	 end do
	 write(buff(7:),'(''integer, parameter :: id_'',a,1x,''='',1x,i3)') &
	             name(:end), i
         write(30,'(a)') buff(:len_trim(buff))
      end do
      buff = ' '
      write(30,*) ' '
      buff(7:) = 'end module m_spc_id'
      write(30,'(a)') buff
      CLOSE(30)

      end subroutine MAKE_NAME_MOD
