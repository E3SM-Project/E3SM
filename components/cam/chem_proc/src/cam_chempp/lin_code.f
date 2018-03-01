
      module lin_matrix

      use io, only : temp_path

      implicit none

      character(len=4) :: hdr, up_hdr
      character(len=4) :: num_suffix 
      character(len=4) :: dec_suffix

      contains

      subroutine make_lin( clscnt, clsmap, cls_rxt_cnt, cls_rxt_map, pcoeff_ind, &
                           pcoeff, permute, mat_map, class, &
                           lin_mat_pat, march, model )
!-----------------------------------------------------------------------
!        ... Write the fortran code for the linear components
!	     of the Jacobian matrix
!-----------------------------------------------------------------------

      use var_mod, only : var_lim
      use rxt_mod, only : rxt_lim, prd_lim

      implicit none
     
!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::      clscnt, &
                                  class, &
                                  clsmap(var_lim,5,2), &
                                  cls_rxt_map(rxt_lim,prd_lim+3), &
                                  cls_rxt_cnt(4)
      integer, intent(in) ::      mat_map(max(1,clscnt),max(1,clscnt))
      integer, intent(in) ::      permute(max(1,clscnt))
      integer, intent(in) ::      pcoeff_ind(*)
      real, intent(in)    ::      pcoeff(prd_lim,*)
      character(len=16), intent(in) ::  model                  ! target model
      character(len=16), intent(in) ::  march                  ! target architecture
      logical, intent(out)::      lin_mat_pat(:)
      
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer, parameter ::   max_len= 90
      integer  ::   i, j, k, l, m
      integer  ::   length, index
      integer  ::   row, col, sub_cnt
      integer  ::   line_pos, rxno, target, line_cnt
      integer  ::   base
      integer  ::   species
      integer  ::   mat_ind
      integer  ::   match_cnt
      integer  ::   list_cnt
      integer  ::   other_ind
      integer  ::   match_ind(var_lim)
      integer  ::   scan(var_lim,3)
      real     ::   rate
      character(len=max_len+2) :: line
      character(len=72) :: buff
      character(len=12) :: het_piece
      character(len= 6) :: mat_piece, rxt_piece
      character(len= 4) :: sol_piece, num
      logical  ::  beg_line, flush
      logical  ::  lexist
      
      if( class == 4 ) then
         inquire( file = trim( temp_path ) // 'linmat.F', exist = lexist )
         if( lexist ) then
            call system( 'rm ' // trim( temp_path ) // 'linmat.F' )
         end if
         open( unit = 30, file = trim( temp_path ) // 'linmat.F' )
         up_hdr = 'imp_'
         hdr    = 'imp_'
         if( model == 'CAM' ) then
            up_hdr = ' '
            hdr    = ' '
         end if
      else
         open( unit = 30, file = trim( temp_path ) // 'linmat.F', position='append' )
         up_hdr = 'rod_'
         hdr    = 'rod_'
      end if

      if( model == 'CAM' ) then
         num_suffix = '_r8'
         dec_suffix = '(r8)'
      else
         num_suffix = ' '
         dec_suffix = ' '
      end if

      line_cnt = 0
      line = ' '
      write(30,100) trim(line)
      line = '      module mo_' // trim(up_hdr) // 'lin_matrix'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '      private'
      write(30,100) trim(line)
      line = '      public :: linmat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '      contains'
      write(30,100) trim(line)
      if( clscnt == 0 .or. (cls_rxt_cnt(2)+cls_rxt_cnt(4)) == 0 ) then
         sub_cnt = 0
      else
         sub_cnt = 1
      end if
      call make_lin_hdr( clscnt, sub_cnt, march, model )

      select case( march )
         case( 'SCALAR' )
            mat_piece = 'mat('
            rxt_piece = 'rxt('
            sol_piece = 'y('
            het_piece = 'het_rates('
         case ( 'CACHE','VECTOR' )
            mat_piece = 'mat(k,'
            rxt_piece = 'rxt(k,'
            sol_piece = 'y(k,'
            het_piece = 'het_rates(k,'
         case default
            mat_piece = 'mat(k,'
            rxt_piece = 'rxt(k,'
            sol_piece = 'y(k,'
            het_piece = 'het_rates(k,'
      end select

      lin_mat_pat(:) = .false.
Species_loop : &
      do species = 1,clscnt
         target = clsmap(species,class,2)
         flush  = .false.
!-----------------------------------------------------------------------
!       ...Write code for linear loss entries
!-----------------------------------------------------------------------
         row     = permute(species)
         mat_ind = mat_map(row,row)
         write(num,'(i4)') mat_map(row,row)
         num       = adjustl( num )
         l         = len_trim( num )
         line      = ' '
         line(10:) = trim( mat_piece ) // num(:l) // ') = -('
         line_pos  = len_trim( line ) + 2
         base      = cls_rxt_cnt(1)
         beg_line  = .true.
         do k = base+1,base+cls_rxt_cnt(2)
            if( cls_rxt_map(k,2) == target ) then
               lin_mat_pat(mat_ind) = .true.
               flush                = .true.
               write(num,'(i4)') cls_rxt_map(k,1)
               num  = adjustl( num )
               l    = len_trim( num )
               buff = trim( rxt_piece ) // num(:l) // ')'
               if( cls_rxt_map(k,3) > 0 ) then
                  write(num,'(i4)') cls_rxt_map(k,3)
                  num = adjustl( num )
                  l   = len_trim( num )
                  buff(len_trim(buff)+1:) = '*' // trim( sol_piece ) // num(:l) // ')'
               end if
               length = len_trim(buff)
               if( (line_pos + length) <= max_len-2 ) then
                  if( beg_line ) then
                     line(line_pos:) = buff(:length)
                     beg_line        = .false.
                  else
                     line(line_pos:) = ' + ' // buff(:length)
                  end if
               else
                  line(len_trim(line)+1:) = ' &'
                  write(30,100) trim(line)
                  line_cnt  = line_cnt + 1
                  line      = ' '
                  line(23:) = '+ ' // buff(:length)
               end if
               line_pos = len_trim( line ) + 1
            end if
         end do
         base = base + cls_rxt_cnt(2) + cls_rxt_cnt(3)
         do k = base+1,base+cls_rxt_cnt(4)
            if( cls_rxt_map(k,2) == species ) then
               lin_mat_pat(mat_ind) = .true.
               flush = .true.
               write(num,'(i4)') cls_rxt_map(k,1)
               num = adjustl( num )
               l = len_trim( num )
               buff = trim( het_piece ) // num(:l) // ')'
               length = len_trim(buff)
               if( (line_pos + length) <= max_len-2 ) then
                  if( beg_line ) then
                     line(line_pos:) = buff(:length)
                     beg_line = .false.
                  else
                     line(line_pos:) = ' + ' // buff(:length)
                  end if
               else
                  line(len_trim(line)+1:) = ' &'
                  write(30,100) trim(line)
                  line_cnt = line_cnt + 1
                  line = ' '
                  line(18:) = '+ ' // buff(:length)
               end if
               line_pos = len_trim( line ) + 1
            end if
         end do
         if( flush ) then
            if( line_pos <= max_len-2 ) then
               line(line_pos+1:) = ')'
            else
               line(len_trim(line)+1:) = ' &'
               write(30,100) trim(line)
               line_cnt = line_cnt + 1
               line = '       )'
            end if
            write(30,100) trim(line)
            line_cnt = line_cnt + 1
         end if
               
!-----------------------------------------------------------------------
!       ... Scan for production matches
!-----------------------------------------------------------------------
         match_cnt = 0
         base = cls_rxt_cnt(1)
         do k = base+1,base+cls_rxt_cnt(2)
            other_ind = 0
            do l = 4,prd_lim+3
               if( cls_rxt_map(k,l) == species ) then
                  if( other_ind == 0 ) then
                     match_cnt = match_cnt + 1
                     scan(match_cnt,1) = k
                     scan(match_cnt,2) = ABS(cls_rxt_map(k,2))
                  end if
                  other_ind = other_ind + 1
               end if
            end do
            if( other_ind /= 0 ) then
               scan(match_cnt,3) = other_ind
            end if
         end do
         list_cnt = match_cnt
         do while( list_cnt > 0 )
            do j = 1,match_cnt
               if( scan(j,2) /= 0 ) then
                  index = scan(j,2)
                  exit
               end if
            end do
            m = 0
            do j = 1,match_cnt
               if( scan(j,2) == index ) then
                  m = m + 1
                  match_ind(m) = j
                  scan(j,2)    = 0
                  list_cnt = list_cnt - 1
               end if
            end do
            row = permute(species)
            col = permute(clsmap(index,class,1))
            mat_ind = mat_map(row,col)
            lin_mat_pat(mat_ind) = .true.
            write(num,'(i4)') mat_map(row,col)
            num = adjustl( num )
            l = len_trim( num )
            line = ' '
            line(10:) = trim( mat_piece ) // num(:l) // ') ='
            line_pos = len_trim( line )
            if( clsmap(index,class,1) == species ) then
               line(len_trim(line)+2:) = trim( mat_piece ) // num(:l) // ') +'
            end if
            line_pos = len_trim( line ) + 2
            if( m > 0 ) then
               beg_line = .true.
            else
               line(line_pos:) = '0.'
            end if
            do j = 1,m
               l     = match_ind(j)
               rxno  = cls_rxt_map(scan(l,1),1)
               index = pcoeff_ind(rxno)
               rate = 0.
               if( index /= 0 ) then
                  do i = 4,prd_lim+3
                     if( cls_rxt_map(scan(l,1),i) == species ) then
                        rate = rate + pcoeff(i-3,index)
                     end if
                  end do
               else if( scan(l,3) /= 1 ) then
                  rate = REAL(scan(l,3))
               end if
               buff = ' '
               if( rate /= 0. .and. rate /= 1. ) then
                  call r2c( buff, rate, 'l' )
                  buff(len_trim(buff)+1:) = trim(num_suffix) // '*'
               end if
               write(num,'(i4)') rxno
               num = adjustl( num )
               buff(len_trim(buff)+1:) = trim( rxt_piece ) // num(:len_trim(num)) // ')'
               if( cls_rxt_map(scan(l,1),3) > 0 ) then
                  write(num,'(i4)') cls_rxt_map(scan(l,1),3)
                  num = adjustl( num )
                  buff(len_trim(buff)+1:) = '*' // trim( sol_piece )  // num(:len_trim(num)) // ')'
               end if
               length = len_trim(buff)
               if( (line_pos + length) <= max_len-2 ) then
                  if( beg_line ) then
                     line(line_pos:) = buff(:length)
                     beg_line = .false.
                  else
                     line(line_pos:) = ' + ' // buff(:length)
                  end if
               else
                  line(len_trim(line)+1:) = ' &'
                  write(30,100) trim(line)
                  line_cnt = line_cnt + 1
                  line = ' '
                  line(23:) = '+ ' // buff(:length)
               end if
               line_pos = len_trim( line ) + 1
            end do
            write(30,100) trim(line)
            line_cnt = line_cnt + 1
         end do       
         line = ' '
         write(30,100) trim(line)
         if( line_cnt > 200 ) then
            if( march /= 'SCALAR' ) then
               line = '      end do'
               write(30,100) trim(line)
            end if
            line = ' '
            write(30,100) trim(line)
            write(num,'(i4)') 1000+sub_cnt
            write(line,'(''      end subroutine '',a,''linmat'',a)') trim(up_hdr),num(3:4)
            write(30,100) trim(line)
            line_cnt = 0
            if( species /= clscnt ) then
               sub_cnt = sub_cnt + 1
               call make_lin_hdr( clscnt, sub_cnt, march, model )
            end if
         end if
      end do Species_loop

      if( line_cnt /= 0  ) then
         if( march /= 'SCALAR' ) then
            line = '      end do'
            write(30,100) trim(line)
         end if
         line = ' '
         write(30,100) trim(line)
         write(num,'(i4)') 1000+sub_cnt
         write(line,'(''      end subroutine '',a,''linmat'',a)') trim(up_hdr),num(3:4)
         write(30,100) trim(line)
      end if

      if( clscnt > 0 .and. (cls_rxt_cnt(2)+cls_rxt_cnt(4)) > 0 ) then
         call make_lin_hdr( clscnt, 0, march, model )
      end if
      do m = 1,sub_cnt
         write(num,'(i4)') 1000+m
         select case( march )
            case ( 'VECTOR' )
               write(line,'(''      call '',a,''linmat'',a,''( ofl, ofu, mat, y, rxt, het_rates )'')') trim(up_hdr),num(3:4)
            case ( 'SCALAR' )
               write(line,'(''      call '',a,''linmat'',a,''( mat, y, rxt, het_rates )'')') trim(up_hdr),num(3:4)
            case default
               if( model == 'MOZART' ) then
                  write(line,'(''      call '',a,''linmat'',a,''( mat, y, rxt, het_rates )'')') trim(up_hdr),num(3:4)
               else if( model == 'CAM' ) then
                  write(line,'(''      call '',a,''linmat'',a,''( mat, y, rxt, het_rates, cols )'')') trim(up_hdr),num(3:4)
               else if( model == 'WRF' ) then
                  write(line,'(''      call '',a,''linmat'',a,''( mat, y, rxt )'')') trim(up_hdr),num(3:4)
               end if
         end select
         write(30,100) trim(line)
      end do
      line = ' '
      write(30,100) trim(line)
      line = '      end subroutine ' // trim(up_hdr) // 'linmat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '      end module mo_' // trim(up_hdr) // 'lin_matrix'
      write(30,100) trim(line)

      close( 30 )

100   format(a)

      end subroutine make_lin

      subroutine make_lin_hdr( clscnt, sub_cnt, march, model )
!-----------------------------------------------------------------------
!        ... Write the fortran header code for the linear components
!	     of the Jacobian matrix
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!	... Dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: clscnt, sub_cnt
      character(len=16), intent(in) :: march
      character(len=16), intent(in) :: model

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: length
      character(len=72) :: line
      character(len=3)  :: num

      line = ' '
      write(30,100) trim(line)
      write(num,'(i3)') 100+sub_cnt
      select case( march )
         case ( 'SCALAR' )
            if( sub_cnt /= 0 ) then
               if( model == 'MOZART' .or. model == 'CAM' ) then
                  write(line,'(''      subroutine '',a,''linmat'',a,''( mat, y, rxt, het_rates )'')') trim(up_hdr),num(2:3)
               else if( model == 'WRF' ) then
                  write(line,'(''      subroutine '',a,''linmat'',a,''( mat, y, rxt )'')') trim(up_hdr),num(2:3)
               end if
            else
               if( model == 'MOZART' .or. model == 'CAM' ) then
                  write(line,'(''      subroutine '',a,''linmat( mat, y, rxt, het_rates )'')') trim(up_hdr)
               else if( model == 'WRF' ) then
                  write(line,'(''      subroutine '',a,''linmat( mat, y, rxt )'')') trim(up_hdr)
               end if
            end if
         case ( 'VECTOR' )
            if( sub_cnt /= 0 ) then
               write(line,'(''      subroutine '',a,''linmat'',a,''( ofl, ofu, mat, y, rxt, het_rates )'')') trim(up_hdr),num(2:3)
            else
               write(line,'(''      subroutine '',a,''linmat( ofl, ofu, mat, y, rxt, het_rates )'')') trim(up_hdr)
            end if
         case default
            if( sub_cnt /= 0 ) then
               if( model == 'MOZART' ) then
                  write(line,'(''      subroutine '',a,''linmat'',a,''( mat, y, rxt, het_rates )'')') trim(up_hdr),num(2:3)
               else if( model == 'CAM' ) then
                  write(line,'(''      subroutine '',a,''linmat'',a,''( mat, y, rxt, het_rates, cols )'')') trim(up_hdr),num(2:3)
               else if( model == 'WRF' ) then
                  write(line,'(''      subroutine '',a,''linmat'',a,''( mat, y, rxt )'')') trim(up_hdr),num(2:3)
               end if
            else
               if( model == 'MOZART' ) then
                  write(line,'(''      subroutine '',a,''linmat( mat, y, rxt, het_rates )'')') trim(up_hdr)
               else if( model == 'CAM' ) then
                  write(line,'(''      subroutine '',a,''linmat( mat, y, rxt, het_rates, cols )'')') trim(up_hdr)
               else if( model == 'WRF' ) then
                  write(line,'(''      subroutine '',a,''linmat( mat, y, rxt )'')') trim(up_hdr)
               end if
            end if
      end select
      write(30,100) trim(line)
      line = '!----------------------------------------------'
      write(30,100) trim(line)
      line = '!       ... linear matrix entries for'
      length = len_trim( line ) + 2
      line(length:) = 'implicit species'
      write(30,100) trim(line)
      line = '!----------------------------------------------'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( model == 'MOZART' ) then
         select case( march )
            case( 'SCALAR' )
               line = '      use mo_grid,   only : pcnstm1'
               write(30,100) trim(line)
               line = '      use chem_mods, only : rxntot, ' // hdr // 'nzcnt'
            case ( 'VECTOR' )
               line = '      use mo_grid,   only : pcnstm1, plnplv'
               write(30,100) trim(line)
               line = '      use chem_mods, only : rxntot, ' // hdr // 'nzcnt'
            case default
               line = '      use mo_grid,   only : pcnstm1'
               write(30,100) trim(line)
               line = '      use chem_mods, only : rxntot, ' // hdr // 'nzcnt, clsze'
         end select
      else if( model == 'CAM' ) then
         select case( march )
            case( 'SCALAR' )
               line = '      use chem_mods, only : gas_pcnst, rxntot, nzcnt'
            case ( 'VECTOR' )
               if( model /= 'CAM' ) then
                  line = '      use chem_mods, only : gas_pcnst, rxntot, nzcnt'
               else
                  line = ' '
               end if
            case default
               line = '      use chem_mods, only : gas_pcnst, rxntot, nzcnt, clsze'
         end select
         write(30,100) trim(line)
         line = '      use shr_kind_mod, only : r8 => shr_kind_r8'
      else if( model == 'WRF' ) then
         line = ''
      end if
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '      implicit none '
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '!----------------------------------------------'
      write(30,100) trim(line)
      line = '!       ... dummy arguments'
      write(30,100) trim(line)
      line = '!----------------------------------------------'
      write(30,100) trim(line)
      if( model /= 'CAM' ) then
         select case( march )
         case( 'SCALAR' )
            if( model /= 'WRF' ) then
               line = '      real, intent(in)    ::  y(pcnstm1)'
               write(30,100) trim(line)
               line = '      real, intent(in)    ::  rxt(rxntot)'
               write(30,100) trim(line)
               line = '      real, intent(in)    ::  het_rates(gas_pcnst)'
               write(30,100) trim(line)
               line = '      real, intent(inout) ::  mat(' // hdr // 'nzcnt)'
            else
               line = '      real, intent(in)    ::  y(:)'
               write(30,100) trim(line)
               line = '      real, intent(in)    ::  rxt(:)'
               write(30,100) trim(line)
               line = '      real, intent(inout) ::  mat(:)'
            end if
            write(30,100) trim(line)
	 case ( 'VECTOR' )
            line = '      integer, intent(in) ::  ofl'
            write(30,100) trim(line)
            line = '      integer, intent(in) ::  ofu'
            write(30,100) trim(line)
            line = '      real, intent(in)    ::  y(:,:)'
            write(30,100) trim(line)
            line = '      real, intent(in)    ::  rxt(:,:)'
            write(30,100) trim(line)
            if( model /= 'WRF' ) then
               line = '      real, intent(in)    ::  het_rates(:,:)'
               write(30,100) trim(line)
            end if
            line = '      real, intent(inout) ::  mat(:,:)'
            write(30,100) trim(line)
            if( sub_cnt /= 0 ) then
               line = ' '
               line = '!----------------------------------------------'
               write(30,100) trim(line)
               line = '!       ... Local variables'
               write(30,100) trim(line)
               line = '!----------------------------------------------'
               write(30,100) trim(line)
               line = '      integer :: k'
               write(30,100) trim(line)
               line = ' '
               write(30,100) trim(line)
               line = ' '
               write(30,100) trim(line)
               line(7:) = 'do k = ofl,ofu'
               write(30,100) trim(line)
            end if
         case default
            line = '      real, intent(in)    ::  y(clsze,pcnstm1)'
            write(30,100) trim(line)
            line = '      real, intent(in)    ::  rxt(clsze,rxntot)'
            write(30,100) trim(line)
            if( model /= 'WRF' ) then
               line = '      real, intent(in)    ::  het_rates(clsze,gas_pcnst)'
               write(30,100) trim(line)
	    end if
            line = '      real, intent(inout) ::  mat(clsze,' // hdr // 'nzcnt)'
            write(30,100) trim(line)
	    if( sub_cnt /= 0 ) then
               line = ' '
               line = '!----------------------------------------------'
               write(30,100) trim(line)
               line = '!       ... local variables'
               write(30,100) trim(line)
               line = '!----------------------------------------------'
               write(30,100) trim(line)
               line = '      integer :: k'
               write(30,100) trim(line)
               line = ' '
               write(30,100) trim(line)
               line = ' '
               write(30,100) trim(line)
               if( clscnt /= 0 ) then
                  line(7:) = 'do k = 1,clsze'
                  write(30,100) trim(line)
               end if
            end if
         end select
      else
         select case( march )
         case ( 'SCALAR' )
            line = '      real(r8), intent(in)    ::  y(gas_pcnst)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  rxt(rxntot)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  het_rates(max(1,gas_pcnst))'
            write(30,100) trim(line)
            line = '      real(r8), intent(inout) ::  mat(nzcnt)'
            write(30,100) trim(line)
         case ( 'VECTOR' )
            line = '      integer,  intent(in)    ::  ofl, ofu'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  y(:,:)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  rxt(:,:)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  het_rates(:,:)'
            write(30,100) trim(line)
            line = '      real(r8), intent(inout) ::  mat(:,:)'
            write(30,100) trim(line)
	    if( sub_cnt /= 0 ) then
               line = ' '
               line = '!----------------------------------------------'
               write(30,100) trim(line)
               line = '!       ... local variables'
               write(30,100) trim(line)
               line = '!----------------------------------------------'
               write(30,100) trim(line)
               line = '      integer :: k'
               write(30,100) trim(line)
               line = ' '
               write(30,100) trim(line)
               line = ' '
               write(30,100) trim(line)
               if( clscnt /= 0 ) then
                  if( model == 'CAM' ) then
                     line(7:) = 'do k = ofl,ofu'
                  end if
                  write(30,100) trim(line)
               end if
            end if
         case default
            line = '      integer, intent(in)     ::  cols'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  y(clsze,gas_pcnst)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  rxt(clsze,rxntot)'
            write(30,100) trim(line)
            line = '      real(r8), intent(in)    ::  het_rates(clsze,max(1,gas_pcnst))'
            write(30,100) trim(line)
            line = '      real(r8), intent(inout) ::  mat(clsze,nzcnt)'
            write(30,100) trim(line)
	    if( sub_cnt /= 0 ) then
               line = ' '
               line = '!----------------------------------------------'
               write(30,100) trim(line)
               line = '!       ... local variables'
               write(30,100) trim(line)
               line = '!----------------------------------------------'
               write(30,100) trim(line)
               line = '      integer :: k'
               write(30,100) trim(line)
               line = ' '
               write(30,100) trim(line)
               line = ' '
               write(30,100) trim(line)
               if( clscnt /= 0 ) then
                  if( model == 'MOZART' ) then
                     line(7:) = 'do k = 1,clsze'
                  else if( model == 'CAM' ) then
                     line(7:) = 'do k = 1,cols'
                  end if
                  write(30,100) trim(line)
               end if
            end if
         end select
      end if
      line = ' '
      write(30,100) trim(line)

100   format(a)

      end subroutine make_lin_hdr

      end module lin_matrix
