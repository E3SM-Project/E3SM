
      module lu_factor

      use io, only : temp_path

      implicit none

      character(len=4) :: hdr, up_hdr
      character(len=4) :: num_suffix 
      character(len=4) :: dec_suffix

      contains

      subroutine make_lu_fac( n, class, lu_sp_pat, mat_sp_pat, sp_map, &
                              march, model )
!-----------------------------------------------------------------------
!        ... Write the fortran code for the sparse matrix decomposition
!-----------------------------------------------------------------------

      implicit none
     
!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: n                       ! species in class count
      integer, intent(in) :: class                   ! class number
      integer, intent(in) :: sp_map(n,n)             ! sparsity matrix map
      character(len=16), intent(in) :: march          ! target architecture
      character(len=16), intent(in) :: model          ! target model
      logical, intent(in), dimension(n,n) :: lu_sp_pat, mat_sp_pat
      
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: max_lines = 50
      integer           :: i, ip1, j, k, l, row, col, sub_cnt
      integer           :: indx, pos, line_cnt
      character(len=90) :: code
      character(len=72) :: comment, blank, buff
      character(len= 6) :: mat_piece
      character(len= 4) :: num
      logical           :: lexist
      logical           :: sp_pat(n,n)
      
!-----------------------------------------------------------------------
!        ... Create and open code file; if it exists remove first
!-----------------------------------------------------------------------
      line_cnt = 0
      if( class == 4 ) then
         inquire( file = trim( temp_path ) // 'lu_fac.F', exist = lexist )
         if( lexist ) then
	    call system( 'rm ' // trim( temp_path ) // 'lu_fac.F' )
         end if
         open( unit = 30, file = trim( temp_path ) // 'lu_fac.F' )
         if( model /= 'CAM' ) then
	    hdr    = 'imp_'
	    up_hdr = 'imp_'
         else
	    hdr    = ' '
	    up_hdr = ' '
         end if
      else
         open( unit = 30, file = trim( temp_path ) // 'lu_fac.F', position='append' )
	 hdr    = 'rod_'
	 up_hdr = 'rod_'
      end if

      if( model == 'CAM' ) then
         num_suffix = '_r8'
         dec_suffix = '(r8)'
      else
         num_suffix = ' '
         dec_suffix = ' '
      end if

      if( n == 0 ) then
	 sub_cnt = 0
      else
	 sub_cnt = 1
      end if
      code = ' '
      write(30,100) trim(code)
      if( model == 'MOZART' ) then
         if( class == 4 ) then
            code = '      module mo_imp_factor'
         else if( class == 5 ) then
            code = '      module mo_rod_factor'
         end if
      else if( model == 'CAM' ) then
         if( class == 4 ) then
            code = '      module mo_lu_factor'
         end if
      else if( model == 'WRF' ) then
         if( class == 4 ) then
            code = '      module mo_imp_factor'
         end if
      end if
      write(30,100) trim(code)
      code = ' '
      write(30,100) trim(code)
      code = '      private'
      write(30,100) trim(code)
      code = '      public :: lu_fac'
      write(30,100) trim(code)
      code = ' '
      write(30,100) trim(code)
      code = '      contains'
      write(30,100) trim(code)
      call make_lu_fac_hdr( sub_cnt, march, model )

      code = ' ' ; blank = ' '
      if( n > 0 ) then
         sp_pat = mat_sp_pat
         comment = '!------------------------------------------------------------------------'
         select case( march )
            case( 'SCALAR' )
	       mat_piece = 'lu('
	    case default
	       mat_piece = 'lu(k,'
         end select
      end if

Column_loop : &
      do i = 1,n
!-----------------------------------------------------------------------
!        ... Form diagonal inverse
!-----------------------------------------------------------------------
         indx = sp_map(i,i)
         write(num,'(i4)') indx
	 num = adjustl( num )
	 l   = len_trim( num )
         code(10:) = trim( mat_piece ) // num(:l) // ') = 1.' // trim(num_suffix) // ' / ' // trim( mat_piece ) // num(:l) // ')'
         write(30,100) trim(code)
	 line_cnt = line_cnt + 1
         buff = ' * ' // trim( mat_piece ) // num(:l) // ')'
	 ip1  = i + 1
!-----------------------------------------------------------------------
!        ... Multiply column below diagonal
!-----------------------------------------------------------------------
	 do row = ip1,n
	    if( sp_pat(row,i) ) then
               indx = sp_map(row,i)
               write(num,'(i4)') indx
	       num = adjustl( num )
	       l   = len_trim( num )
               code(10:) = trim( mat_piece ) // num(:l) // ') = ' // trim( mat_piece ) // num(:l) // ')' &
				   // buff(:len_trim(buff))
               write(30,100) trim(code)
	       line_cnt = line_cnt + 1
	    end if
	 end do
!-----------------------------------------------------------------------
!        ... Modify sub-matrix
!-----------------------------------------------------------------------
	 do col = ip1,n
	    if( sp_pat(i,col) ) then
               indx = sp_map(i,col)
               write(num,'(i4)') indx
	       num = adjustl( num )
	       l = len_trim( num )
               buff = ' * ' // trim( mat_piece ) // num(:l) // ')'
	       do row = ip1,n
	          if( sp_pat(row,i) ) then
                     indx = sp_map(row,col)
                     write(num,'(i4)') indx
	             num = adjustl( num )
	             l = len_trim( num )
	             if( sp_pat(row,col) ) then
                        code(10:) = trim( mat_piece ) // num(:l) // ') = ' // trim( mat_piece ) // num(:l) // ')'
                        indx = sp_map(row,i)
                        write(num,'(i4)') indx
	                num = adjustl( num )
	                l = len_trim( num )
                        code(len_trim(code)+2:) = '- ' // trim( mat_piece ) // num(:l) // ')' // buff(:len_trim(buff))
                        write(30,100) trim(code)
		        code(6:) = ' '
		     else
			sp_pat(row,col) = .true.
                        code(10:) = trim( mat_piece ) // num(:l) // ') = '
                        indx = sp_map(row,i)
                        write(num,'(i4)') indx
	                num = adjustl( num )
	                l = len_trim( num )
		        pos = index( code,'=' ) + 2
                        code(pos:) = '- ' // trim( mat_piece ) // num(:l) // ')' // buff(:len_trim(buff))
                        write(30,100) trim(code)
	                line_cnt = line_cnt + 1
	             end if
	          end if
	       end do
	    end if
	 end do
         write(30,100) blank
	 if( line_cnt > max_lines ) then
            if( march /= 'SCALAR' ) then
               code = '      end do'
               write(30,100) trim(code)
            end if
            write(30,100) blank
            write(num,'(i3)') 100+sub_cnt
            write(code,'(''      end subroutine '',a,''lu_fac'',a)') trim(up_hdr),num(2:3)
            write(30,100) trim(code)
	    line_cnt = 0
	    if( i /= n ) then
	       sub_cnt  = sub_cnt + 1
               call make_lu_fac_hdr( sub_cnt, march, model )
	    end if
            code = ' '
	 end if
      end do Column_loop

      if( line_cnt /= 0 ) then
         if( march /= 'SCALAR' ) then
            code(7:) = 'end do'
            write(30,100) trim(code)
         end if
         write(30,100) blank
         write(num,'(i3)') 100+sub_cnt
         write(code,'(''      end subroutine '',a,''lu_fac'',a)') trim(up_hdr),num(2:3)
         write(30,100) trim(code)
      end if

      if( n > 0 ) then
         call make_lu_fac_hdr( 0, march, model )
      end if
      do k = 1,sub_cnt
         write(num,'(i3)') 100+k
	 select case( march )
	    case( 'SCALAR' )
               write(code,'(''      call '',a,''lu_fac'',a,''( lu )'')') trim(up_hdr),num(2:3)
	    case( 'VECTOR' )
               write(code,'(''      call '',a,''lu_fac'',a,''( ofl, ofu, lu )'')') trim(up_hdr),num(2:3)
	    case default
               if( model /= 'CAM' ) then
                  write(code,'(''      call '',a,''lu_fac'',a,''( lu )'')') trim(up_hdr),num(2:3)
               else
                  write(code,'(''      call '',a,''lu_fac'',a,''( lu, cols )'')') trim(up_hdr),num(2:3)
               end if
	 end select
         write(30,100) trim(code)
      end do
      write(30,100) blank
      code(7:) = 'end subroutine ' // trim(up_hdr) // 'lu_fac'
      write(30,100) trim(code)
      write(30,100) blank

      if( model == 'MOZART' ) then
         if( class == 4 ) then
            code = '      end module mo_imp_factor'
         else if( class == 5 ) then
            code = '      end module mo_rod_factor'
         end if
      else if( model == 'CAM' ) then
         if( class == 4 ) then
            code = '      end module mo_lu_factor'
         end if
      else if( model == 'WRF' ) then
         if( class == 4 ) then
            code = '      end module mo_imp_factor'
         end if
      end if
      write(30,100) trim(code)

      close( 30 )

100   format(a)

      end subroutine make_lu_fac

      subroutine make_lu_fac_hdr( sub_cnt, march, model )
!-----------------------------------------------------------------------
!        ... Write the fortran header code for the sparse matrix decomposition
!-----------------------------------------------------------------------

      implicit none
     
!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)          :: sub_cnt
      character(len=16), intent(in) :: march
      character(len=16), intent(in) :: model
      
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer  ::   i, ip1, j, k, l, row, col
      integer  ::   indx, pos
      character(len=72) :: code, comment, blank, buff
      character(len=3)  :: num

      code = ' ' ; blank = ' '
      comment = '!------------------------------------------------------------------------'

      write(30,100) blank
      write(num,'(i3)') 100+sub_cnt
      select case( march )
         case( 'SCALAR' )
            if( sub_cnt /= 0 ) then
               write(code,'(''      subroutine '',a,''lu_fac'',a,''( lu )'')') trim(up_hdr),num(2:3)
            else
               write(code,'(''      subroutine '',a,''lu_fac( lu )'')') trim(up_hdr)
            end if
         case( 'VECTOR' )
            if( sub_cnt /= 0 ) then
               write(code,'(''      subroutine '',a,''lu_fac'',a,''( ofl, ofu, lu )'')') trim(up_hdr),num(2:3)
            else
               write(code,'(''      subroutine '',a,''lu_fac( ofl, ofu, lu )'')') trim(up_hdr)
            end if
         case default
            if( sub_cnt /= 0 ) then
               if( model /= 'CAM' ) then
                  write(code,'(''      subroutine '',a,''lu_fac'',a,''( lu )'')') trim(up_hdr),num(2:3)
               else
                  write(code,'(''      subroutine '',a,''lu_fac'',a,''( lu, cols )'')') trim(up_hdr),num(2:3)
               end if
            else
               if( model /= 'CAM' ) then
                  write(code,'(''      subroutine '',a,''lu_fac( lu )'')') trim(up_hdr)
               else
                  write(code,'(''      subroutine '',a,''lu_fac( lu, cols )'')') trim(up_hdr)
               end if
            end if
      end select
      write(30,100) trim(code)
      write(30,100) blank
      if( march == 'SCALAR' ) then
         code(:) = ' '
      else if( march == 'VECTOR' ) then
         if( model /= 'CAM' ) then
            code(7:) = 'use mo_grid,   only : plnplv'
            write(30,100) trim(code)
            code(7:) = 'use chem_mods, only : ' // hdr // 'nzcnt'
         else
            code(:) = ' '
         end if
      else
         if( model /= 'WRF' ) then
            code(7:) = 'use chem_mods, only : ' // hdr // 'nzcnt, clsze'
         else
            code(:) = ' '
         end if
      end if
      write(30,100) trim(code)
      if( model == 'CAM' ) then
         code(7:) = 'use shr_kind_mod, only : r8 => shr_kind_r8'
         write(30,100) trim(code)
      end if
      write(30,100) blank
      code(7:) = 'implicit none '
      write(30,100) trim(code)
      write(30,100) blank
      write(30,100) comment
      code = '!       ... dummy args'
      write(30,100) trim(code)
      write(30,100) comment
      code = ' '
      if( model == 'CAM' .and. march == 'CACHE' ) then
         code(7:) = 'integer, intent(in) :: cols'
         write(30,100) trim(code)
      end if
      select case( march )
         case( 'SCALAR' )
            code(7:) = 'real' // trim(dec_suffix) // ', intent(inout) ::   lu(:)'
         case( 'VECTOR' )
            code(7:) = 'integer, intent(in) ::   ofl'
            write(30,100) trim(code)
            code(7:) = 'integer, intent(in) ::   ofu'
            write(30,100) trim(code)
            if( model /= 'CAM' ) then
               code(7:) = 'real' // trim(dec_suffix) // ', intent(inout) ::   lu(plnplv,' // hdr // 'nzcnt)'
            else
               code(7:) = 'real' // trim(dec_suffix) // ', intent(inout) ::   lu(:,:)'
            end if
         case default
            code(7:) = 'real' // trim(dec_suffix) // ', intent(inout) ::   lu(clsze,' // hdr // 'nzcnt)'
      end select
      write(30,100) trim(code)
      write(30,100) blank
      if( sub_cnt /= 0 ) then
         if( march /= 'SCALAR' ) then
            write(30,100) comment
            code = '!       ... local variables'
            write(30,100) trim(code)
            write(30,100) comment
            code = ' '
            code(7:) = 'integer :: k'
            write(30,100) trim(code)
            write(30,100) blank
            code = ' '
            if( march == 'VECTOR' ) then
               code(7:) = 'do k = ofl,ofu'
            else if( march == 'CACHE' ) then
               if( model == 'MOZART' ) then
                  code(7:) = 'do k = 1,clsze'
               else if( model == 'CAM' ) then
                  code(7:) = 'do k = 1,cols'
               end if
            end if
            write(30,100) trim(code)
         end if
      end if

100   format(a)

      end subroutine make_lu_fac_hdr

      end module lu_factor
