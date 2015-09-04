
      module nln_matrix

      use io, only : temp_path

      implicit none

      character(len=9) :: spc_cnt
      character(len=4) :: hdr, up_hdr
      character(len=4) :: num_suffix 
      character(len=4) :: dec_suffix

      contains

      subroutine make_nln( clscnt, clsmap, cls_rxt_cnt, cls_rxt_map, pcoeff_ind, &
                           pcoeff, permute, mat_map, class, &
                           lin_mat_pat, nzcnt, diag_map, march, model )
!-----------------------------------------------------------------------
!        ... Write the fortran code for the non-linear components
!	     of the Jacobian matrix
!-----------------------------------------------------------------------
     
      use var_mod, only : var_lim
      use rxt_mod, only : rxt_lim, prd_lim

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::      clscnt                   ! count of class members
      integer, intent(in) ::      class                    ! class index
      integer, intent(in) ::      nzcnt                    ! matrix non-zero count
      integer, intent(in) ::      clsmap(var_lim,5,2)   
      integer, intent(in) ::      cls_rxt_map(rxt_lim,prd_lim+3)   
      integer, intent(in) ::      cls_rxt_cnt(4)           ! class rxtns count
      integer, intent(in) ::      permute(clscnt)
      integer, intent(in) ::      mat_map(clscnt,clscnt)
      integer, intent(in) ::      diag_map(:)
      integer, intent(in) ::      pcoeff_ind(*)            ! map for nonunity prod
      real, intent(in)    ::      pcoeff(prd_lim,*)
      character(len=16), intent(in) ::  march               ! target architecture  
      character(len=16), intent(in) ::  model               ! target model  
      logical, intent(in) ::      lin_mat_pat(:)
      
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: max_len   = 90
      integer, parameter :: max_lines = 200
      integer  ::   i, j, k, l, m, m2, n, r1, r2
      integer  ::   length, index, pindx, mat_ind
      integer  ::   row, col, sub_cnt
      integer  ::   line_pos, buf_pos, rxno, target, line_cnt
      integer  ::   base
      integer  ::   species
      integer  ::   match_cnt
      integer  ::   list_cnt
      integer  ::   rxtnt_cnt, rxtnt1, rxtnt2
      integer  ::   other_ind
      integer  ::   match_ind(rxt_lim)
      integer  ::   rxt_match_ind(rxt_lim)
      integer  ::   scan(rxt_lim,4)
      integer  ::   rxtnt(2)
      real     ::   rate
      character(len=max_len+10) :: line
      character(len=max_len) :: buff
      character(len= 6) :: mat_piece, rxt_piece
      character(len= 4) :: sol_piece, num, num1
      logical  ::  beg_line
      logical  ::  lexist
      logical, allocatable  ::  nln_mat_pat(:)
      logical :: hdr_made = .false.

      allocate( nln_mat_pat(nzcnt),stat=pindx )
      if( pindx /= 0 ) then
	 stop
      end if
      nln_mat_pat(:) = .false.
      
      if( class == 4 ) then
         inquire( file = trim( temp_path ) // 'nlnmat.F', exist = lexist )
         if( lexist ) then
	    call system( 'rm ' // trim( temp_path ) // 'nlnmat.F' )
         end if
         open( unit = 30, file = trim( temp_path ) // 'nlnmat.F' )
         if( model /= 'CAM' ) then
	    up_hdr = 'imp_'
	    hdr    = 'imp_'
         else
	    up_hdr = ' '
	    hdr    = ' '
         end if
      else
         open( unit = 30, file = trim( temp_path ) // 'nlnmat.F', position='append' )
	 up_hdr = 'rod_'
	 hdr    = 'rod_'
      end if

      if( model == 'CAM' ) then
         num_suffix = '_r8'
         dec_suffix = '(r8)'
         spc_cnt    = 'gas_pcnst'
      else
         num_suffix = ' '
         dec_suffix = ' '
         spc_cnt    = 'pcnstm1'
      end if

      line_cnt = 0
      line = ' '
      write(30,100) trim(line)
      line = '      module mo_' // trim(up_hdr) // 'nln_matrix'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( model == 'CAM' ) then
         line = '      use shr_kind_mod, only : r8 => shr_kind_r8'
         write(30,100) trim(line)
         line = ' '
         write(30,100) trim(line)
      end if
      line = '      private'
      write(30,100) trim(line)
      line = '      public :: nlnmat'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '      contains'
      write(30,100) trim(line)
      if( clscnt == 0 .or. cls_rxt_cnt(3) == 0 ) then
	 sub_cnt = 0
      else
	 sub_cnt = 1
      end if
      call make_nln_hdr( sub_cnt, march, model )
      if (sub_cnt>0) hdr_made = .true.

      select case ( march )
         case( 'SCALAR' )
	    mat_piece = 'mat('
	    rxt_piece = 'rxt('
	    sol_piece = 'y('
	 case default
	    mat_piece = 'mat(k,'
	    rxt_piece = 'rxt(k,'
	    sol_piece = 'y(k,'
      end select

      base = sum( cls_rxt_cnt(:2) )
Species_loop : &
      do species = 1,clscnt
        target = clsmap(species,class,2)
	 line   = ' '
!-----------------------------------------------------------------------
!       ... Write code for nonlinear loss entries
!-----------------------------------------------------------------------
         match_cnt = 0
         do k = base+1,base+cls_rxt_cnt(3)
!-----------------------------------------------------------------------
!       ... Find all reactions with target reactant
!-----------------------------------------------------------------------
            other_ind = 0
            do l = 2,3
               if( cls_rxt_map(k,l) == target ) then
                  if( other_ind == 0 ) then
                     match_cnt = match_cnt + 1
                     scan(match_cnt,1) = k
                     if( l == 2 ) then
                        scan(match_cnt,2) = abs(cls_rxt_map(k,3))
                     else
                        scan(match_cnt,2) = abs(cls_rxt_map(k,2))
                     end if
                     scan(match_cnt,4) = l
                  end if
                  other_ind = other_ind + 1
               end if
            end do
         end do
!-----------------------------------------------------------------------
!       ... Write the diagonal loss entry
!-----------------------------------------------------------------------
         if( match_cnt > 0 ) then
            scan(:match_cnt,3) = scan(:match_cnt,2)
	    pindx = permute(species)
	    mat_ind = mat_map(pindx,pindx)
            write(num,'(i4)') mat_map(pindx,pindx)
	    num = adjustl( num )
	    n = len_trim( num )
            line = ' '
            line(10:) = trim( mat_piece ) // num(:n) // ') = -('
            line_pos = len_trim( line ) + 1
            beg_line = .true.
         end if
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
            do j = 1,m
               l = match_ind(j)
               rxno = cls_rxt_map(scan(l,1),1)
               buff = ' '
	       buf_pos = 1
               if( j == 1 .and. m > 1 ) then
                  if( scan(l,3) == target ) then
                     buff(buf_pos:) =  '(4.' // trim(num_suffix) // '*'
                  else
                     buff(buf_pos:) =  '('
                  end if
               else if( scan(l,3) == target ) then
                  buff(buf_pos:) =  '4.' // trim(num_suffix) // '*'
               end if
               write(num,'(i4)') rxno
	       num = adjustl( num )
	       n = len_trim( num )
               buff(len_trim(buff)+1:) = trim( rxt_piece ) // num(:n) // ')'
               length = len_trim(buff)
               if( (line_pos + length) <= max_len-3 ) then
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
            write(num,'(i4)') scan(l,3)
	    num = adjustl( num )
            if( m > 1 ) then
               buff = ') * ' // trim( sol_piece ) // num(:len_trim(num)) // ')'
            else
               buff = '*' // trim( sol_piece ) // num(:len_trim(num)) // ')'
            end if
            length = len_trim(buff)
            if( (line_pos + length) <= max_len-3 ) then
               line(line_pos:) = buff(:length)
            else
	       line(len_trim(line)+1:) = ' &'
               write(30,100) trim(line)
	       line_cnt = line_cnt + 1
               line = ' '
               line(23:) = buff(:length)
            end if
            line_pos = len_trim( line ) + 1
	    nln_mat_pat(mat_ind) = .true.
         end do
	 if( match_cnt /= 0 ) then
	    line(len_trim(line)+1:) = ')'
	 end if
	 if( line /= ' ' ) then
            write(30,100) trim(line)
            line_cnt = line_cnt + 1
	 end if
         
!-----------------------------------------------------------------------
!       ... Write nondiagonal loss entries
!-----------------------------------------------------------------------
         list_cnt = match_cnt
         do j = 1,match_cnt
            if( scan(j,3) == target ) then
               scan(j,2) = 0
               list_cnt = list_cnt - 1
            else
               scan(j,2) = scan(j,3)
            end if
         end do
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
	    pindx = permute(clsmap(index,class,1))
            mat_ind = mat_map(permute(species),pindx)
            write(num,'(i4)') mat_map(permute(species),pindx)
	    num = adjustl( num )
	    n = len_trim( num )
            line = ' '
            line(10:) = trim( mat_piece ) // num(:n) // ') = -'
            line_pos = len_trim( line ) + 1
            if( m > 0 ) then
               beg_line = .true.
	       nln_mat_pat(mat_ind) = .true.
            else
               if( model /= 'CAM' ) then
                  line(line_pos:) = '0.'
               else
                  line(line_pos:) = '0._r8'
               end if
            end if
            do j = 1,m
               l = match_ind(j)
               rxno = cls_rxt_map(scan(l,1),1)
               buff = ' '
               if( j == 1 .and. m > 1 ) then
                  buff =  '('
                  buf_pos = 2
               else
                  buf_pos = 1
               end if
               write(num,'(i4)') rxno
	       num = adjustl( num )
	       n = len_trim( num )
               buff(buf_pos:) = trim( rxt_piece ) // num(:n) // ')'
               if( j == 1 ) then
                  if ( scan(l,4) == 2 ) then
                     index = 3
                  else
                     index = 2
                  end if
               end if
               length = len_trim(buff)
               if( (line_pos + length) <= max_len-3 ) then
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
            write(num,'(i4)') target
	    num = adjustl( num )
            if( m > 1 ) then
               buff = ') * '// trim( sol_piece ) // num(:len_trim(num)) // ')'
            else
               buff = '*' // trim( sol_piece ) // num(:len_trim(num)) // ')'
            end if
            length = len_trim(buff)
            if( (line_pos + length) <= max_len-3 ) then
               line(line_pos:) = buff(:length)
            else
	       line(len_trim(line)+1:) = ' &'
               write(30,100) trim(line)
	       line_cnt = line_cnt + 1
               line = ' '
               line(23:) = buff(:length)
            end if
            write(30,100) trim(line)
            line_cnt = line_cnt + 1
         end do               
         line = ' '
         write(30,100) trim(line)

!-----------------------------------------------------------------------
!       ... Scan for production matches
!-----------------------------------------------------------------------
         match_cnt = 0
Product_match : &
         do k = base+1,base+cls_rxt_cnt(3)
            other_ind = 0
            do l = 4,prd_lim+3
               if( cls_rxt_map(k,l) == species ) then
                  if( other_ind == 0 ) then
                     match_cnt = match_cnt + 1
                     scan(match_cnt,1) = k
                     scan(match_cnt,2) = abs(cls_rxt_map(k,2))
                     scan(match_cnt,4) = abs(cls_rxt_map(k,3))
                  end if
                  other_ind = other_ind + 1
               end if
            end do
            if( other_ind /= 0 ) then
               scan(match_cnt,3) = other_ind
            end if
         end do Product_match
	 if( match_cnt == 0 ) then
	    cycle
	 end if
!-----------------------------------------------------------------------
!       ... "Order" the match list reactants
!-----------------------------------------------------------------------
         do j = 1,match_cnt
            if( scan(j,2) > scan(j,4) ) then
               l = scan(j,2)
               scan(j,2) = scan(j,4)
               scan(j,4) = l
            end if
         end do
!-----------------------------------------------------------------------
!       ... Search matching reactions for reactant match
!-----------------------------------------------------------------------
Reactant_match : &
         do r1 = 1,clscnt
	    m = 0
	    rxtnt1 = clsmap(r1,class,2)
            do j = 1,match_cnt
               if( scan(j,2) == rxtnt1 .or. scan(j,4) == rxtnt1 ) then
                  m = m + 1
                  match_ind(m) = j
               end if
	    end do
	    if( m == 0 ) then
	       cycle
	    end if
	    pindx   = permute(clsmap(rxtnt1,class,1))
            mat_ind = mat_map(permute(species),pindx)
            write(num,'(i4)') mat_ind
	    num = adjustl( num )
            line      = ' '
            line(10:) = trim( mat_piece ) // num(:len_trim(num)) // ') ='
            beg_line  = .true.
	    if( nln_mat_pat(mat_ind) ) then
	       line_pos        = len_trim(line) + 2
               line(line_pos:) = trim( mat_piece ) // num(:len_trim(num)) // ')'
               beg_line        = .false.
	    end if
            nln_mat_pat(mat_ind) = .true.
Second_reactant_match : &
            do r2 = 1,clscnt
	       m2 = 0
	       rxtnt2 = clsmap(r2,class,2)
               do n = 1,m
	          j = match_ind(n)
		  if( rxtnt2 /= rxtnt1 ) then
                     if( scan(j,2) == rxtnt2 .or. scan(j,4) == rxtnt2 ) then
                        m2 = m2 + 1
                        rxt_match_ind(m2) = j
                     end if
		  else if( scan(j,2) == rxtnt2 .and. scan(j,4) == rxtnt2 ) then
                     m2 = m2 + 1
                     rxt_match_ind(m2) = j
                  end if
	       end do
	       if( m2 == 0 ) then
	          cycle
	       end if
	       if( .not. beg_line ) then
		  line(len_trim(line)+2:) = '+'
	       else
		  beg_line = .false.
	       end if
	       if( m2 > 1 ) then
	          line(len_trim(line)+2:) = '('
                  line_pos = len_trim( line ) + 1
	       else
                  line_pos = len_trim( line ) + 2
	       end if
Rates_loop : &
               do n = 1,m2
!-----------------------------------------------------------------------
!       ... The reaction rate
!-----------------------------------------------------------------------
                  l        = rxt_match_ind(n)
                  rxno     = cls_rxt_map(scan(l,1),1)
                  index    = pcoeff_ind(rxno)
                  rate     = 0.
                  if( index /= 0 ) then
                     do i = 4,prd_lim+3
                        if( cls_rxt_map(scan(l,1),i) == species ) then
                           rate = rate + pcoeff(i-3,index)
                        end if
                     end do
                  else if( scan(l,3) /= 1 ) then
                     rate = REAL(scan(l,3))
                  end if
                  if( rxtnt1 == rxtnt2 ) then
                     if( rate == 0. ) then
                        rate = 2.
                     else
                        rate = 2.*rate
                     end if
                  end if
                  buff = ' '
		  if( n > 1 .and. m2 > 1 ) then
                     buff = '+'
		  end if
                  if( rate /= 0. .and. rate /= 1. ) then
                     call r2c( buff(len_trim(buff)+1:), rate, 'l' )
                     buff(len_trim(buff)+1:) = trim(num_suffix) // '*'
                  end if
                  write(num,'(i4)') rxno
	          num = adjustl( num )
                  buff(len_trim(buff)+1:) = trim( rxt_piece ) // num(:len_trim(num)) // ')'
                  length = len_trim(buff)
                  if( (line_pos + length) <= max_len-3 ) then
                     line(line_pos:) = buff(:length)
                  else
	             if( line(len_trim(line):len_trim(line)) /= '+' ) then
	                line(len_trim(line)+1:) = ' &'
		     else
	                line(len_trim(line):) = ' &'
		     end if
                     write(30,100) trim(line)
		     line_cnt  = line_cnt + 1
                     line      = ' '
		     if( buff(1:1) /= '+' ) then
                        line(23:) = '+ ' // buff(:length)
		     else
                        line(23:) = ' ' // buff(:length)
		     end if
!                    line(23:) = '+ ' // buff(:length)
                  end if
                  line_pos = len_trim( line ) + 1
               end do Rates_loop
!-----------------------------------------------------------------------
!       ... The reactant
!-----------------------------------------------------------------------
	       if( m2 > 1 ) then
	          line(len_trim(line)+1:) = ')'
                  line_pos = len_trim( line ) + 1
	       end if
               write(num,'(i4)') rxtnt2
	       num    = adjustl( num )
               buff   = '*' // trim( sol_piece ) // num(:len_trim(num)) // ')'
               length = len_trim(buff)
               if( (line_pos + length) <= max_len-3 ) then
                  line(line_pos:) = buff(:length)
               else
	          line(len_trim(line)+1:) = ' &'
                  write(30,100) trim(line)
	          line_cnt  = line_cnt + 1
                  line      = ' '
                  line(23:) = buff(:length)
               end if
            end do Second_Reactant_match
            write(30,100) trim(line)
	    line_cnt  = line_cnt + 1
         end do Reactant_match
         line = ' '
         write(30,100) trim(line)
	 if( line_cnt > max_lines ) then
            if( march /= 'SCALAR' ) then
               line = '      end do'
               write(30,100) trim(line)
            end if
            line = ' '
            write(30,100) trim(line)
            write(num,'(i3)') 100+sub_cnt
            write(line,'(''      end subroutine '',a,''nlnmat'',a)') up_hdr,num(2:3)
            write(30,100) trim(line)
            hdr_made = .false.
	    line_cnt = 0
	    if( species /= clscnt ) then
	       sub_cnt  = sub_cnt + 1
               call make_nln_hdr( sub_cnt, march, model )
               hdr_made = .true.
	    end if
	 end if
      end do Species_loop
     if ( hdr_made ) then
         if( march /= 'SCALAR' ) then
            line = '      end do'
            write(30,100) trim(line)
         end if

         line = ' '
         write(30,100) trim(line)
         write(num,'(i3)') 100+sub_cnt
         write(line,'(''      end subroutine '',a,''nlnmat'',a)') up_hdr,num(2:3)
         write(30,100) trim(line)
      end if
!-----------------------------------------------------------------------
!	... Make the inclusion routine
!-----------------------------------------------------------------------
      if( clscnt > 0 ) then
	 if( cls_rxt_cnt(3) == 0 ) then
            if( model == 'MOZART' ) then
	       select case( march )
	       case ( 'VECTOR' )
                  write(line,'(''      call '',a,''nlnmat_finit( ofl, ofu, mat, lmat, dti )'')') up_hdr
	       case default
                  write(line,'(''      call '',a,''nlnmat_finit( mat, lmat, dti )'')') up_hdr
	       end select
            else
               write(line,'(''      call '',a,''nlnmat_finit( mat, lmat, dti )'')') up_hdr
	    end if
            write(30,100) trim(line)
            line = ' '
            write(30,100) trim(line)
            line = '      end subroutine ' // trim(up_hdr) // 'nlnmat'
            write(30,100) trim(line)
	 end if
         call make_nln_hdr( -1, march, model )
         line = ' '
	 do n = 1,size(lin_mat_pat)
	    if( lin_mat_pat(n) ) then
               write(num,'(i4)') n
	       m = len_trim( num )
	       if( nln_mat_pat(n) ) then
                  line(10:) = trim( mat_piece ) // num(:m) // ') = ' // trim(mat_piece) // num(:m) // ') + l' &
			      // trim(mat_piece) // num(:m) // ')'
               else
                  line(10:) = trim( mat_piece ) // num(:m) // ') = l' // trim(mat_piece) // num(:m) // ')'
	       end if
               write(30,100) trim(line)
	    end if
	 end do
         line = ' '
	 do n = 1,size(lin_mat_pat)
	    if( .not. lin_mat_pat(n) .and. .not. nln_mat_pat(n) ) then
               write(num,'(i4)') n
	       m = len_trim( num )
               if( model /= 'CAM' ) then
                  line(10:) = trim( mat_piece ) // num(:m) // ') = 0.'
               else
                  line(10:) = trim( mat_piece ) // num(:m) // ') = 0._r8'
               end if
               write(30,100) trim(line)
	    end if
	 end do
	 do n = 1,size(diag_map)
	    l = diag_map(n)
	    if( lin_mat_pat(l) .or. nln_mat_pat(l) ) then
               write(num,'(i4)') l
	       m = len_trim( num )
               line(10:) = trim( mat_piece ) // num(:m) // ') = ' // trim(mat_piece) // num(:m) // ') - dti'
	    else
               write(num,'(i4)') l
	       m = len_trim( num )
               line(10:) = trim( mat_piece ) // num(:m) // ') = -dti'
	    end if
            write(30,100) trim(line)
	 end do
	 if( march /= 'SCALAR' ) then
            line = '      end do'
            write(30,100) trim(line)
	 end if
         line = ' '
         write(30,100) trim(line)
         line = '      end subroutine ' // trim(up_hdr) // 'nlnmat_finit'
         write(30,100) trim(line)
      end if
!-----------------------------------------------------------------------
!	... Now make the driver routine
!-----------------------------------------------------------------------
      if( clscnt > 0 .and. cls_rxt_cnt(3) > 0 ) then
         call make_nln_hdr( 0, march, model )
      end if
      do n = 1,sub_cnt
         write(num,'(i3)') 100+n
	 select case( march )
	    case ( 'SCALAR' )
               write(line,'(''      call '',a,''nlnmat'',a,''( mat, y, rxt )'')') up_hdr,num(2:3)
	    case ( 'VECTOR' )
               write(line,'(''      call '',a,''nlnmat'',a,''( ofl, ofu, mat, y, rxt )'')') up_hdr,num(2:3)
	    case default
               if( model /= 'CAM' ) then
                  write(line,'(''      call '',a,''nlnmat'',a,''( mat, y, rxt )'')') up_hdr,num(2:3)
               else
                  write(line,'(''      call '',a,''nlnmat'',a,''( mat, y, rxt, cols )'')') up_hdr,num(2:3)
               end if
	 end select
         write(30,100) trim(line)
      end do
      if( clscnt > 0 .and. cls_rxt_cnt(3) > 0 ) then
	 select case( march )
	    case ( 'SCALAR' )
               write(line,'(''      call '',a,''nlnmat_finit( mat, lmat, dti )'')') up_hdr
	    case ( 'VECTOR' )
               write(line,'(''      call '',a,''nlnmat_finit( ofl, ofu, mat, lmat, dti )'')') up_hdr
	    case default
               if( model /= 'CAM' ) then
                  write(line,'(''      call '',a,''nlnmat_finit( mat, lmat, dti )'')') up_hdr
               else
                  write(line,'(''      call '',a,''nlnmat_finit( mat, lmat, dti, cols )'')') up_hdr
               end if
	 end select
         write(30,100) trim(line)
         line = ' '
         write(30,100) trim(line)
         line = '      end subroutine ' // trim(up_hdr) // 'nlnmat'
         write(30,100) trim(line)
      end if
      if( clscnt == 0 ) then
         line = ' '
         write(30,100) trim(line)
         line = '      end subroutine ' // trim(up_hdr) // 'nlnmat'
         write(30,100) trim(line)
      end if
      line = ' '
      write(30,100) trim(line)
      line = '      end module mo_' // trim(up_hdr) // 'nln_matrix'
      write(30,100) trim(line)

      if( allocated( nln_mat_pat ) ) then
         deallocate( nln_mat_pat )
      end if

      close( 30 )

100   format(a)

      end subroutine make_nln

      subroutine make_nln_hdr( sub_cnt, march, model )
!-----------------------------------------------------------------------
!        ... Write the fortran header code for the non-linear components
!	     of the Jacobian matrix
!-----------------------------------------------------------------------
     
      implicit none

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in)          :: sub_cnt              ! subroutine counter
      character(len=16), intent(in) :: march                ! targe  architecture
      character(len=16), intent(in) :: model                ! target model
      
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer           :: length
      character(len=72) :: line
      character(len=3)  :: num
      
      line = ' '
      write(30,100) trim(line)
      write(num,'(i3)') 100+sub_cnt
      select case( march )
	 case ( 'SCALAR' )
            if( sub_cnt > 0 ) then
               write(line,'(''      subroutine '',a,''nlnmat'',a,''( mat, y, rxt )'')') up_hdr,num(2:3)
            else if( sub_cnt < 0 ) then
               write(line,'(''      subroutine '',a,''nlnmat_finit( mat, lmat, dti )'')') up_hdr
            else
               write(line,'(''      subroutine '',a,''nlnmat( mat, y, rxt, lmat, dti )'')') up_hdr
            end if
	 case ( 'VECTOR' )
            if( sub_cnt > 0 ) then
               write(line,'(''      subroutine '',a,''nlnmat'',a,''( ofl, ofu, mat, y, rxt )'')') up_hdr,num(2:3)
            else if( sub_cnt < 0 ) then
               write(line,'(''      subroutine '',a,''nlnmat_finit( ofl, ofu, mat, lmat, dti )'')') up_hdr
            else
               write(line,'(''      subroutine '',a,''nlnmat( ofl, ofu, mat, y, rxt, lmat, dti )'')') up_hdr
            end if
	 case default
            if( model /= 'CAM' ) then
               if( sub_cnt > 0 ) then
                  write(line,'(''      subroutine '',a,''nlnmat'',a,''( mat, y, rxt )'')') up_hdr,num(2:3)
               else if( sub_cnt < 0 ) then
                  write(line,'(''      subroutine '',a,''nlnmat_finit( mat, lmat, dti )'')') up_hdr
               else
                  write(line,'(''      subroutine '',a,''nlnmat( mat, y, rxt, lmat, dti )'')') up_hdr
               end if
            else
               if( sub_cnt > 0 ) then
                  write(line,'(''      subroutine '',a,''nlnmat'',a,''( mat, y, rxt, cols )'')') up_hdr,num(2:3)
               else if( sub_cnt < 0 ) then
                  write(line,'(''      subroutine '',a,''nlnmat_finit( mat, lmat, dti, cols )'')') up_hdr
               else
                  write(line,'(''      subroutine '',a,''nlnmat( mat, y, rxt, lmat, dti, cols )'')') up_hdr
               end if
            end if
      end select
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      select case( march )
	 case( 'SCALAR' )
            if( model /= 'WRF' ) then
               if( model == 'MOZART' ) then
                  line = '      use mo_grid,   only : ' // trim(spc_cnt)
                  write(30,100) trim(line)
                  line = '      use chem_mods, only : rxntot, ' // hdr // 'nzcnt'
               else
                  line = '      use chem_mods, only : gas_pcnst, rxntot, ' // hdr // 'nzcnt'
               end if
            else
               line = ' '
            end if
	 case ( 'VECTOR' )
            if( model == 'MOZART' ) then
               line = '      use mo_grid,   only : plnplv, ' // trim(spc_cnt)
               write(30,100) trim(line)
               line = '      use chem_mods, only : rxntot, ' // hdr // 'nzcnt'
            end if
	 case default
            if( model == 'MOZART' ) then
               line = '      use mo_grid,   only : ' // trim(spc_cnt)
               write(30,100) trim(line)
               line = '      use chem_mods, only : rxntot, ' // hdr // 'nzcnt, clsze'
            else if( model == 'CAM' ) then
               line = '      use chem_mods, only : gas_pcnst, rxntot, nzcnt, clsze'
               write(30,100) trim(line)
               line = ' '
            end if
      end select
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
      if( model == 'CAM' .and. march == 'CACHE' ) then
         line = '      integer, intent(in)    ::  cols'
         write(30,100) trim(line)
      end if
      select case( march )
	 case( 'SCALAR' )
	    if( sub_cnt <= 0 ) then
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  dti'
               write(30,100) trim(line)
	       if( model /= 'WRF' ) then
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  lmat(' // trim(hdr) // 'nzcnt)'
                  write(30,100) trim(line)
	       else
                  line = '      real, intent(in)    ::  lmat(:)'
                  write(30,100) trim(line)
	       end if
	    end if
	    if( sub_cnt >= 0 ) then
	       if( model /= 'WRF' ) then
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(' // trim(spc_cnt) // ')'
                  write(30,100) trim(line)
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(rxntot)'
                  write(30,100) trim(line)
	       else
                  line = '      real, intent(in)    ::  y(:)'
                  write(30,100) trim(line)
                  line = '      real, intent(in)    ::  rxt(:)'
                  write(30,100) trim(line)
	       end if
	    end if
	    if( model /= 'WRF' ) then
               line = '      real' // trim(dec_suffix) // ', intent(inout) ::  mat(' // trim(hdr) // 'nzcnt)'
	    else
               line = '      real, intent(inout) ::  mat(:)'
	    end if
	 case ( 'VECTOR' )
            line = '      integer, intent(in) ::  ofl'
            write(30,100) trim(line)
            line = '      integer, intent(in) ::  ofu'
            write(30,100) trim(line)
	    if( sub_cnt <= 0 ) then
               line = '      real, intent(in)    ::  dti'
               write(30,100) trim(line)
               if( model /= 'CAM' ) then
                  line = '      real, intent(in)    ::  lmat(plnplv,' // trim(hdr) // 'nzcnt)'
               else
                  line = '      real, intent(in)    ::  lmat(:,:)'
               end if
               write(30,100) trim(line)
	    end if
	    if( sub_cnt >= 0 ) then
               if( model /= 'CAM' ) then
                  line = '      real, intent(in)    ::  y(plnplv,' // trim(spc_cnt) // ')'
                  write(30,100) trim(line)
                  line = '      real, intent(in)    ::  rxt(plnplv,rxntot)'
               else
                  line = '      real, intent(in)    ::  y(:,:)'
                  write(30,100) trim(line)
                  line = '      real, intent(in)    ::  rxt(:,:)'
               end if
               write(30,100) trim(line)
	    end if
            if( model /= 'CAM' ) then
               line = '      real, intent(inout) ::  mat(plnplv,' // hdr // 'nzcnt)'
            else
               line = '      real, intent(inout) ::  mat(:,:)'
            end if
	 case ( 'CACHE' )
	    if( sub_cnt <= 0 ) then
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  dti'
               write(30,100) trim(line)
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  lmat(clsze,' // trim(hdr) // 'nzcnt)'
               write(30,100) trim(line)
	    end if
	    if( sub_cnt >= 0 ) then
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(clsze,' // trim(spc_cnt) // ')'
               write(30,100) trim(line)
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(clsze,rxntot)'
               write(30,100) trim(line)
	    end if
            line = '      real' // trim(dec_suffix) // ', intent(inout) ::  mat(clsze,' // trim(hdr) // 'nzcnt)'
	 case default
	    if( sub_cnt <= 0 ) then
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  dti'
               write(30,100) trim(line)
               line = '      real' // trim(dec_suffix) // ', intent(in) ::  lmat(plnplv,' // trim(hdr) // 'nzcnt)'
               write(30,100) trim(line)
	    end if
	    if( sub_cnt >= 0 ) then
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(plnplv,' // trim(spc_cnt) // ')'
               write(30,100) trim(line)
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(plnplv,rxntot)'
               write(30,100) trim(line)
	    end if
            line = '      real' // trim(dec_suffix) // ', intent(inout) ::  mat(plnplv,' // trim(hdr) // 'nzcnt)'
      end select
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( sub_cnt /= 0 ) then
         line = ' '
         write(30,100) trim(line)
         line = '!----------------------------------------------'
         write(30,100) trim(line)
         line = '!       ... local variables'
         write(30,100) trim(line)
         line = '!----------------------------------------------'
         write(30,100) trim(line)
         if( march /= 'SCALAR' ) then
            line = '      integer :: k'
            write(30,100) trim(line)
         end if
         line = ' '
         write(30,100) trim(line)

         line = '!----------------------------------------------'
         write(30,100) trim(line)
         line = '!       ... complete matrix entries'
         length = len_trim( line ) + 2
         line(length:) = 'implicit species'
         write(30,100) trim(line)
         line = '!----------------------------------------------'
         write(30,100) trim(line)
         line = ' '
         write(30,100) trim(line)
	 if( march == 'VECTOR' ) then
            line(7:) = 'do k = ofl,ofu'
         else if( march == 'CACHE' ) then
            if( model == 'MOZART' ) then
               line(7:) = 'do k = 1,clsze'
            else if( model == 'CAM' ) then
               line(7:) = 'do k = 1,cols'
	    end if
	 end if
         write(30,100) trim(line)
      end if

100   format(a)

      end subroutine make_nln_hdr
      
      end module nln_matrix
