
      subroutine MAKE_HET_NAME_MOD
!--------------------------------------------------------------------------------
!	... Makes a module of parameter reaction names
!--------------------------------------------------------------------------------

      use RXT_MOD, only : hetcnt, hetmap
      use VAR_MOD, only : spc_cnt => new_nq, spc_names => new_solsym
      use IO,      only : temp_path

      implicit none

!--------------------------------------------------------------------------------
!	... Local variables
!--------------------------------------------------------------------------------
      integer :: i, m
      character(len=80) :: buff
      character(len=5)  :: num
      logical :: lexist

!--------------------------------------------------------------------------------
!	... Check mod file existence; remove if found
!--------------------------------------------------------------------------------
      INQUIRE( file = TRIM( temp_path ) // 'het_names.mod', exist = lexist )
      if( lexist ) then
         call SYSTEM( 'rm ' // TRIM( temp_path ) // 'het_names.mod' )
      end if
      OPEN( unit = 30, file = TRIM( temp_path ) // 'het_names.mod' )
             
      buff = ''
      write(30,'(a)') buff
      buff(7:) = 'module m_het_id'
      write(30,'(a)') buff
      buff = ''
      write(30,'(a)') buff
      buff(7:) = 'implicit none'
      write(30,'(a)') buff
      buff = ''
      write(30,'(a)') buff

      do i = 1,hetcnt
	 m = hetmap(i,1)
	 write(buff(7:),'(''integer, parameter :: hid_'',a,1x,''='',1x,i4)') &
	             spc_names(m)(:LEN_TRIM(spc_names(m))), i
         write(30,'(a)') buff
      end do

      buff = ''
      write(30,'(a)') buff
      buff(7:) = 'end module m_het_id'
      write(30,'(a)') buff
      CLOSE(30)

      end subroutine MAKE_HET_NAME_MOD
