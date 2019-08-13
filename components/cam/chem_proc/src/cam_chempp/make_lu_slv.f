
      module lu_solve

      use io, only : temp_path

      implicit none

      character(len=4) :: hdr, up_hdr
      character(len=4) :: num_suffix 
      character(len=4) :: dec_suffix

      contains

      subroutine make_lu_slv( n, class, lu_sp_pat, march, model )
!-----------------------------------------------------------------------
!        ... Write the fortran code for the sparse matrix solver
!-----------------------------------------------------------------------

      implicit none
     
!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) :: n                            ! count of species in class
      integer, intent(in) :: class                        ! class number
      character(len=16), intent(in) :: march               ! target architecture
      character(len=16), intent(in) :: model               ! target model
      logical, intent(in), dimension(n,n) :: lu_sp_pat
      
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer           :: i, ip1, j, k, l, row, col, sub_cnt
      integer           :: indx, pos, line_cnt
      integer           :: sp_map(n,n)
      character(len=72) :: code, comment, blank, buff
      character(len= 6) :: mat_piece
      character(len= 4) :: b_piece
      character(len= 4) :: num
      logical           :: lexist
      
!-----------------------------------------------------------------------
!        ... Create and open code file; if it exists remove first
!-----------------------------------------------------------------------
      line_cnt = 0
      if( class == 4 ) then
         inquire( file = trim( temp_path ) // 'lu_slv.F', exist = lexist )
         if( lexist ) then
	    call system( 'rm ' // trim( temp_path ) // 'lu_slv.F' )
         end if
         open( unit = 30, file = trim( temp_path ) // 'lu_slv.F' )
         if( model /= 'CAM' ) then
	    hdr    = 'imp_'
	    up_hdr = 'imp_'
         else
	    hdr    = ' '
	    up_hdr = ' '
         end if
      else
         open( unit = 30, file = trim( temp_path ) // 'lu_slv.F', position='append' )
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
            code = '      module mo_imp_solve'
         else if( class == 5 ) then
            code = '      module mo_rod_solve'
         end if
      else if( model == 'CAM' ) then
         if( class == 4 ) then
            code = '      module mo_lu_solve'
         end if
      else if( model == 'WRF' ) then
         if( class == 4 ) then
            code = '      module mo_imp_solve'
         end if
      end if
      write(30,100) trim(code)
      code = ' '
      write(30,100) trim(code)
      code = '      private'
      write(30,100) trim(code)
      code = '      public :: lu_slv'
      write(30,100) trim(code)
      code = ' '
      write(30,100) trim(code)
      code = '      contains'
      write(30,100) trim(code)
      call make_lu_slv_hdr( n, class, sub_cnt, march, model )

      code = ' ' ; blank = ' '
      if( n > 0 ) then
!-----------------------------------------------------------------------
!        ... Form the lu matrix map
!-----------------------------------------------------------------------
         k = 0 ; sp_map = 0
         do i = 1,n
            do j = 1,n
	       if( lu_sp_pat(j,i) ) then
	          k = k + 1
	          sp_map(j,i) = k
	       end if
            end do
         end do

         code = ' ' ; blank = ' '
         comment = '!------------------------------------------------------------------------'

         if( march == 'SCALAR' ) then
	    mat_piece = 'lu('
	    b_piece = 'b('
         else
	    mat_piece = 'lu(k,'
	    b_piece = 'b(k,'
         end if
      end if

!-----------------------------------------------------------------------
!        ... Solve L * y = b
!-----------------------------------------------------------------------
Forward_loop : &
      do col = 1,n-1
         write(num,'(i4)') col
	 num = adjustl( num )
	 l = len_trim( num )
         buff = ' * ' // trim( b_piece ) // num(:l) // ')'
	 do row = col+1,n
	    if( lu_sp_pat(row,col) ) then
               write(num,'(i4)') row
	       num = adjustl( num )
	       l = len_trim( num )
               code(10:) = trim( b_piece ) // num(:l) // ') = ' // trim( b_piece ) // num(:l) // ')'
               indx = sp_map(row,col)
               write(num,'(i4)') indx
	       num = adjustl( num )
	       l = len_trim( num )
               code(len_trim(code)+2:) = '- ' // trim( mat_piece ) // num(:l) // ')' // buff(:len_trim(buff))
               write(30,100) trim(code)
	       line_cnt = line_cnt + 1
	    end if
	 end do
         write(30,100) blank
	 if( line_cnt > 200 ) then
	    if( march /= 'SCALAR' ) then
	       code = '      end do'
               write(30,100) trim(code)
	    end if
            write(30,100) blank
            write(num,'(i3)') 100+sub_cnt
            write(code,'(''      end subroutine '',a,''lu_slv'',a)') trim(up_hdr),num(2:3)
            write(30,100) trim(code)
	    line_cnt = 0
	    sub_cnt = sub_cnt + 1
	    call make_lu_slv_hdr( n, class, sub_cnt, march, model )
	    code = ' '
	 end if
      end do Forward_loop

      if( line_cnt /= 0 ) then
	 if( march /= 'SCALAR' ) then
	    code = '      end do'
            write(30,100) trim(code)
	 end if
         write(30,100) blank
         write(num,'(i3)') 100+sub_cnt
         write(code,'(''      end subroutine '',a,''lu_slv'',a)') trim(up_hdr),num(2:3)
         write(30,100) trim(code)
	 line_cnt = 0
	 sub_cnt = sub_cnt + 1
	 call make_lu_slv_hdr( n,  class, sub_cnt, march, model )
	 code = ' '
      end if

      if( n > 0 ) then
         write(30,100) blank
         write(30,100) comment
         code = '!       ... Solve U * x = y'
         write(30,100) trim(code)
         write(30,100) comment
         code = ' '
      end if

!-----------------------------------------------------------------------
!        ... Solve U * x = y
!-----------------------------------------------------------------------
Backward_loop : &
      do col = n,1,-1
         write(num,'(i4)') col
	 num = adjustl( num )
	 l = len_trim( num )
         code(10:) = trim( b_piece) // num(:l) // ') = ' // trim( b_piece ) // num(:l) // ')'
         buff = ' * ' // trim( b_piece ) // num(:l) // ')'
         write(num,'(i4)') sp_map(col,col)
	 num = adjustl( num )
	 l = len_trim( num )
         code(len_trim(code)+2:) = '* ' // trim( mat_piece ) // num(:l) // ')'
	 write(30,100) trim(code)
	 line_cnt = line_cnt + 1
	 do row = col-1,1,-1
	    if( lu_sp_pat(row,col) ) then
               write(num,'(i4)') row
	       num = adjustl( num )
	       l = len_trim( num )
               code(10:) = trim( b_piece ) // num(:l) // ') = ' // trim( b_piece ) // num(:l) // ')'
               indx = sp_map(row,col)
               write(num,'(i4)') indx
	       num = adjustl( num )
	       l = len_trim( num )
               code(len_trim(code)+2:) = '- ' // trim( mat_piece ) // num(:l) // ')' // buff(:len_trim(buff))
               write(30,100) trim(code)
	       line_cnt = line_cnt + 1
	    end if
	 end do
         write(30,100) blank
	 if( line_cnt > 200 ) then
	    if( march /= 'SCALAR' ) then
	       code = '      end do'
               write(30,100) trim(code)
	    end if
            write(30,100) blank
            write(num,'(i3)') 100+sub_cnt
            write(code,'(''      end subroutine '',a,''lu_slv'',a)') trim(up_hdr),num(2:3)
            write(30,100) trim(code)
	    line_cnt = 0
	    if( col /= 1 ) then
	       sub_cnt = sub_cnt + 1
	       call make_lu_slv_hdr( n,  class, sub_cnt, march, model )
	    end if
	    code = ' '
	 end if
      end do Backward_loop

      if( line_cnt /= 0 ) then
	 if( march /= 'SCALAR' ) then
	    code = '      end do'
            write(30,100) trim(code)
	 end if
         write(30,100) blank
         write(num,'(i3)') 100+sub_cnt
         write(code,'(''      end subroutine '',a,''lu_slv'',a)') trim(up_hdr),num(2:3)
         write(30,100) trim(code)
      end if

      if( n > 0 ) then
         call make_lu_slv_hdr( n,  class, 0, march, model )
      end if

      do k = 1,sub_cnt
         write(num,'(i3)') 100+k
	 select case( march )
	    case( 'SCALAR' )
               write(code,'(''      call '',a,''lu_slv'',a,''( lu, b )'')') trim(up_hdr),num(2:3)
	    case( 'VECTOR' )
	       write(code,'(''      call '',a,''lu_slv'',a,''( ofl, ofu, lu, b )'')') trim(up_hdr),num(2:3)
	    case default
               if( model /= 'CAM' ) then
                  write(code,'(''      call '',a,''lu_slv'',a,''( lu, b )'')') trim(up_hdr),num(2:3)
               else
                  write(code,'(''      call '',a,''lu_slv'',a,''( lu, b, cols )'')') trim(up_hdr),num(2:3)
               end if
	 end select
         write(30,100) trim(code)
      end do

      write(30,100) blank
      code(7:) = 'end subroutine ' // trim(up_hdr) // 'lu_slv'
      write(30,100) trim(code)
      write(30,100) blank
      if( model == 'MOZART' ) then
         if( class == 4 ) then
            code = '      end module mo_imp_solve'
         else if( class == 5 ) then
            code = '      end module mo_rod_solve'
         end if
      else if( model == 'CAM' ) then
         if( class == 4 ) then
            code = '      end module mo_lu_solve'
         end if
      else if( model == 'WRF' ) then
         if( class == 4 ) then
            code = '      end module mo_imp_solve'
         end if
      end if
      write(30,100) trim(code)

      close( 30 )

100   format(a)

      end subroutine make_lu_slv

      subroutine make_lu_slv_hdr( n, class, sub_cnt, march, model )
!-----------------------------------------------------------------------
!        ... Write the fortran header code for the sparse matrix solver
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!        ... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: n, class
      integer, intent(in) :: sub_cnt
      character(len=16), intent(in) :: march
      character(len=16), intent(in) :: model

!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      character(len=3)  :: num
      character(len=72) :: code, comment, blank


      code = ' ' ; blank = ' '
      comment = '!------------------------------------------------------------------------'

      write(30,100) blank
      write(num,'(i3)') 100+sub_cnt
      select case( march )
         case( 'SCALAR' )
            if( sub_cnt /= 0 ) then
               write(code,'(''      subroutine '',a,''lu_slv'',a,''( lu, b )'')') trim(up_hdr),num(2:3)
            else
               write(code,'(''      subroutine '',a,''lu_slv( lu, b )'')') trim(up_hdr)
            end if
         case( 'VECTOR' )
            if( sub_cnt /= 0 ) then
               write(code,'(''      subroutine '',a,''lu_slv'',a,''( ofl, ofu, lu, b )'')') trim(up_hdr),num(2:3)
            else
               write(code,'(''      subroutine '',a,''lu_slv( ofl, ofu, lu, b )'')') trim(up_hdr)
            end if
         case default
            if( sub_cnt /= 0 ) then
               if( model /= 'CAM' ) then
                  write(code,'(''      subroutine '',a,''lu_slv'',a,''( lu, b )'')') trim(up_hdr),num(2:3)
               else
                  write(code,'(''      subroutine '',a,''lu_slv'',a,''( lu, b, cols )'')') trim(up_hdr),num(2:3)
               end if
            else
               if( model /= 'CAM' ) then
                  write(code,'(''      subroutine '',a,''lu_slv( lu, b )'')') trim(up_hdr)
               else
                  write(code,'(''      subroutine '',a,''lu_slv( lu, b, cols )'')') trim(up_hdr)
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
            write(code(7:),'(''use chem_mods, only : ' // hdr // 'nzcnt, clscnt'',i1)') class
         else
            code(:) = ' '
         end if
      else
         if( model /= 'WRF' ) then
            write(code(7:),'(''use chem_mods, only : ' // hdr // 'nzcnt, clsze, clscnt'',i1)') class
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
      code = '!       ... Dummy args'
      write(30,100) trim(code)
      write(30,100) comment
      code = ' '
      if( model == 'CAM' .and. march == 'CACHE' ) then
         code(7:) = 'integer, intent(in) :: cols'
         write(30,100) trim(code)
      end if
      code = ' '
      if( march == 'SCALAR' ) then
         if( model /= 'CAM' ) then
            code(7:) = 'real, intent(in)    ::   lu(:)'
         else
            code(7:) = 'real(r8), intent(in) ::   lu(:)'
         end if
      else
         if( march == 'VECTOR' ) then
            code(7:) = 'integer, intent(in) ::   ofl'
            write(30,100) trim(code)
            code(7:) = 'integer, intent(in) ::   ofu'
            write(30,100) trim(code)
            if( model /= 'CAM' ) then
               code(7:) = 'real, intent(in)    ::   lu(plnplv,' // hdr // 'nzcnt)'
            else
               code(7:) = 'real(r8), intent(in)    ::   lu(:,:)'
            end if
         else
            code(7:) = 'real, intent(in)    ::   lu(clsze,' // hdr // 'nzcnt)'
         end if
      end if
      write(30,100) trim(code)
      write(num,'(i3)') n
      num = adjustl( num )
      if( march == 'SCALAR' ) then
         code(7:) = 'real' // trim(dec_suffix) // ', intent(inout) ::   b(:)'
      else if( march == 'VECTOR' ) then
         if( model /= 'CAM' ) then
            write(code(7:),'(''real' // trim(dec_suffix) // ', intent(inout) ::   b(plnplv,clscnt'',i1,'')'')') class
         else
            write(code(7:),'(''real' // trim(dec_suffix) // ', intent(inout) ::   b(:,:)'')'')')
         end if
      else
         write(code(7:),'(''real' // trim(dec_suffix) // ', intent(inout) ::   b(clsze,clscnt'',i1,'')'')') class
      end if
      write(30,100) trim(code)
      write(30,100) blank
      if( sub_cnt /= 0 ) then
         write(30,100) comment
         code = '!       ... Local variables'
         write(30,100) trim(code)
         write(30,100) comment
         code = ' '
         if( march /= 'SCALAR' ) then
            code(7:) = 'integer :: k'
            write(30,100) trim(code)
         end if
         write(30,100) blank
         write(30,100) comment
         code = '!       ... solve L * y = b'
         write(30,100) trim(code)
         write(30,100) comment
         if( march /= 'SCALAR' ) then
            code = ' '
            if( march == 'VECTOR' ) then
               code(7:) = 'do k = ofl,ofu'
            else
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

      end subroutine make_lu_slv_hdr

      end module lu_solve
