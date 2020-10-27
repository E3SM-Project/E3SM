
      module prod_loss

      use io, only : temp_path

      contains

      subroutine pl_code( spccnt, clscnt, clsmap, cls_rxt_cnt, cls_rxt_map, &
                          pcoeff_ind, pcoeff, permute, march, model )
!-----------------------------------------------------------------------
!	... Write the fortran production and loss code
!-----------------------------------------------------------------------
     
      use var_mod, only : var_lim
      use rxt_mod, only : rxt_lim, prd_lim

      implicit none

!-----------------------------------------------------------------------
!        ... The arguments
!
!            The columns of the cls_rxt_cnt represent the reaction count
!	     for each class with the following row conontation:
!		(1) - independent reactions
!		(2) - linear reactions
!		(3) - nonlinear reactions
!		(4) - heterogeneous processes
!-----------------------------------------------------------------------
      integer, intent(in) ::  spccnt
      integer, intent(in) ::  clscnt(5), &
                              clsmap(var_lim,5,2), &
                              cls_rxt_map(rxt_lim,prd_lim+3,5), &
                              cls_rxt_cnt(4,5)
      integer, intent(in) ::  pcoeff_ind(rxt_lim)
      integer, intent(in) ::  permute(var_lim,5)
      real, intent(in)    ::  pcoeff(prd_lim,rxt_lim)
      character(len=16), intent(in) :: march
      character(len=16), intent(in) :: model
      
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: max_len = 90
      integer, parameter :: explicit = 1
      integer, parameter :: ebi      = 2
      integer, parameter :: hov      = 3
      integer, parameter :: implicit = 4
      integer, parameter :: rodas    = 5
      integer  ::   i, k, kl, ku, l, m, ml
      integer  ::   length, index, cnt
      integer  ::   line_pos, target
      integer  ::   class
      integer  ::   base
      integer  ::   species
      integer  ::   match_cnt
      integer  ::   other_ind
      integer  ::   match_ind(4)
      integer  ::   max_loc(1)
      integer  ::   freq(spccnt)
      integer, allocatable  ::   indexer(:)
      real     ::   rate
      character(len=max_len) :: line
      character(len=72) :: buff
      character(len=14) :: het_piece
      character(len= 9) :: l_piece
      character(len= 9) :: p_piece
      character(len= 8) :: rxt_piece
      character(len= 6) :: sol_piece
      character(len=4)  :: num_suffix 
      character(len=4)  :: dec_suffix
      character(len= 3) :: num
      logical, allocatable :: match_mask(:,:)
      logical, allocatable :: pmask(:,:)
      logical  ::  beg_line
      logical  ::  lexist, first
      
      inquire( file = trim( temp_path ) // 'prd_loss.F', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'prd_loss.F' )
      end if
      open( unit = 30, file = trim( temp_path ) // 'prd_loss.F' )

      if( model == 'CAM' ) then
         num_suffix = '_r8'
         dec_suffix = '(r8)'
      else
         num_suffix = ' '
         dec_suffix = ' '
      end if

      line = '      module mo_prod_loss'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( model == 'CAM' ) then
         line = '      use shr_kind_mod, only : r8 => shr_kind_r8'
         write(30,100) trim(line)
         line = ' '
         write(30,100) trim(line)
         line = '      private'
         write(30,100) trim(line)
      end if
      do class = 1,5
         select case( class )
            case( 1 )
               line = '      public :: exp_prod_loss'
               write(30,100) trim(line)
            case( 4 )
               line = '      public :: imp_prod_loss'
               write(30,100) trim(line)
         end select
      end do
      line = ' '
      write(30,100) trim(line)
      line = '      contains'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)

Class_loop : &
      do class = 1,5
         if( class == ebi .or. class == hov ) then
	    cycle
         else if( class == rodas .and. model /= 'MOZART' ) then
	    cycle
	 end if
         if( march == 'SCALAR' ) then
	    select case( class )
	       case( explicit )
                  line = '      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )'
	       case( ebi )
                  line = '      subroutine ebi_prod_loss( prod, loss, y, rxt, het_rates )'
	       case( hov )
                  line = '      subroutine hov_prod_loss( prod, loss, y, rxt, het_rates )'
	       case( implicit )
                  line = '      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )'
	       case( rodas )
                  line = '      subroutine rodas_prod_loss( prod, loss, y, rxt, het_rates )'
	    end select
         else if( march == 'VECTOR' ) then
	    select case( class )
	       case( explicit )
                  if( model /= 'CAM' ) then
                     line = '      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )'
                  else
                     line = '      subroutine exp_prod_loss( ofl, ofu, prod, loss, y, &'
                     write(30,100) trim(line)
                     line = '                                rxt, het_rates )'
                  end if
	       case( ebi )
                  line = '      subroutine ebi_prod_loss( prod, loss, y, rxt, het_rates )'
	       case( hov )
                  line = '      subroutine hov_prod_loss( prod, loss, y, rxt, het_rates )'
	       case( implicit )
                  line = '      subroutine imp_prod_loss( ofl, ofu, prod, loss, y, &'
                  write(30,100) trim(line)
                  line = '                                rxt, het_rates )'
	       case( rodas )
                  line = '      subroutine rodas_prod_loss( ofl, ofu, prod, loss, y, &'
                  write(30,100) trim(line)
                  line = '                                  rxt, het_rates )'
	    end select
	 else
            if( model == 'MOZART' ) then
	       select case( class )
	          case( explicit )
                     line = '      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )'
	          case( ebi )
                     line = '      subroutine ebi_prod_loss( prod, loss, y, rxt, het_rates )'
	          case( hov )
                     line = '      subroutine hov_prod_loss( prod, loss, y, rxt, het_rates )'
	          case( implicit )
                     line = '      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates )'
	          case( rodas )
                     line = '      subroutine rodas_prod_loss( prod, loss, y, rxt, het_rates )'
	       end select
            else if( model == 'CAM' ) then
	       select case( class )
	          case( explicit )
                     line = '      subroutine exp_prod_loss( prod, loss, y, rxt, het_rates )'
	          case( ebi )
                     line = '      subroutine ebi_prod_loss( prod, loss, y, rxt, het_rates, cols )'
	          case( hov )
                     line = '      subroutine hov_prod_loss( prod, loss, y, rxt, het_rates, cols )'
	          case( implicit )
                     line = '      subroutine imp_prod_loss( prod, loss, y, rxt, het_rates, cols )'
	          case( rodas )
                     line = '      subroutine rodas_prod_loss( prod, loss, y, rxt, het_rates, cols )'
	       end select
            else if( model == 'WRF' ) then
	       select case( class )
	          case( explicit )
                     line = '      subroutine exp_prod_loss( prod, loss, y, rxt )'
	          case( ebi )
                     line = '      subroutine ebi_prod_loss( prod, loss, y, rxt )'
	          case( hov )
                     line = '      subroutine hov_prod_loss( prod, loss, y, rxt )'
	          case( implicit )
                     line = '      subroutine imp_prod_loss( prod, loss, y, rxt )'
	          case( rodas )
                     line = '      subroutine rodas_prod_loss( prod, loss, y, rxt )'
	       end select
	    end if
	 end if
         write(30,100) trim(line)
         line = ' '
         write(30,100) trim(line)
         if( model == 'CAM' ) then
            if( march /= 'VECTOR' ) then
               line = '      use ppgrid,       only : pver'
               write(30,100) trim(line)
            end if
         end if
         line = ' '
         write(30,100) trim(line)
         line = '      implicit none '
         write(30,100) trim(line)
         line = ' '
         write(30,100) trim(line)
         line = '!--------------------------------------------------------------------'
         write(30,100) trim(line)
         line = '!     ... dummy args'
         write(30,100) line
         line = '!--------------------------------------------------------------------'
         write(30,100) trim(line)
         if( model == 'CAM'  ) then
            if( class == implicit ) then
               if( march == 'CACHE' ) then
                  line = '      integer :: cols'
                  write(30,100) trim(line)
               end if
               line = '      real' // trim(dec_suffix) // ', dimension(:,:), intent(out) :: &'
            else
               if( march /= 'VECTOR' ) then
                  line = '      real' // trim(dec_suffix) // ', dimension(:,:,:), intent(out) :: &'
               else
                  line = '      real' // trim(dec_suffix) // ', dimension(:,:), intent(out) :: &'
               end if
            end if
         else
            line = '      real' // trim(dec_suffix) // ', dimension(:,:), intent(out) :: &'
         end if
	 if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
	    if( march == 'SCALAR' .and. class /= explicit ) then
                  line = '      real' // trim(dec_suffix) // ', dimension(:), intent(out) :: &'
		  p_piece = 'prod('
		  l_piece = 'loss('
		  rxt_piece = 'rxt('
		  het_piece = 'het_rates('
		  sol_piece = 'y('
	    else if( march == 'VECTOR' ) then
		  p_piece = 'prod(k,'
		  l_piece = 'loss(k,'
		  rxt_piece = 'rxt(k,'
		  het_piece = 'het_rates(k,'
		  sol_piece = 'y(k,'
	    else
               if( class == explicit ) then
		     p_piece = 'prod(:,:,'
		     l_piece = 'loss(:,:,'
		     rxt_piece = 'rxt(:,:,'
		     het_piece = 'het_rates(:,:,'
		     sol_piece = 'y(:,:,'
               else
		     p_piece = 'prod(k,'
		     l_piece = 'loss(k,'
		     rxt_piece = 'rxt(k,'
		     het_piece = 'het_rates(k,'
		     sol_piece = 'y(k,'
               end if
	    end if
	 end if
         write(30,100) trim(line)
         line = '            prod, &'
         write(30,100) trim(line)
         line = '            loss'
         write(30,100) trim(line)
	 if( class == explicit ) then
	    if( model /= 'CAM' ) then
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(:,:)'
               write(30,100) trim(line)
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(:,:)'
               write(30,100) trim(line)
               if( model /= 'WRF' ) then
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  het_rates(:,:)'
                  write(30,100) trim(line)
               end if
	    else
               if( march /= 'VECTOR' ) then
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(:,:,:)'
                  write(30,100) trim(line)
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(:,:,:)'
                  write(30,100) trim(line)
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  het_rates(:,:,:)'
               else
                  line = '      integer, intent(in) :: ofl, ofu'
                  write(30,100) trim(line)
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(:,:)'
                  write(30,100) trim(line)
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(:,:)'
                  write(30,100) trim(line)
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  het_rates(:,:)'
               end if
               write(30,100) trim(line)
	    end if
	 else
            if( march == 'SCALAR' ) then
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(:)'
               write(30,100) trim(line)
               line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(:)'
               write(30,100) trim(line)
               if( model /= 'WRF' ) then
                  line = '      real' // trim(dec_suffix) // ', intent(in)    ::  het_rates(:)'
               end if
	    else if( march == 'VECTOR' ) then
                     line = '      integer, intent(in)    ::  ofl'
                     write(30,100) trim(line)
                     line = '      integer, intent(in)    ::  ofu'
                     write(30,100) trim(line)
                     line = '      real' // trim(dec_suffix) // ', intent(in)       ::  y(:,:)'
                     write(30,100) trim(line)
                     line = '      real' // trim(dec_suffix) // ', intent(in)       ::  rxt(:,:)'
                     write(30,100) trim(line)
                     if( model /= 'WRF' ) then
                        line = '      real' // trim(dec_suffix) // ', intent(in)       ::  het_rates(:,:)'
                     end if
            else
	       if( model /= 'CAM' ) then
                        line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(:,:)'
                        write(30,100) trim(line)
                        line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(:,:)'
                        write(30,100) trim(line)
                        if( model /= 'WRF' ) then
                           line = '      real' // trim(dec_suffix) // ', intent(in)    ::  het_rates(:,:)'
                        end if
	       else
                  select case( class )
                     case( explicit )
                           line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(:,:,:)'
                           write(30,100) trim(line)
                           line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(:,:,:)'
                           write(30,100) trim(line)
                           line = '      real' // trim(dec_suffix) // ', intent(in)    ::  het_rates(:,:,:)'
                     case( implicit )
                           line = '      real' // trim(dec_suffix) // ', intent(in)    ::  y(:,:)'
                           write(30,100) trim(line)
                           line = '      real' // trim(dec_suffix) // ', intent(in)    ::  rxt(:,:)'
                           write(30,100) trim(line)
                           line = '      real' // trim(dec_suffix) // ', intent(in)    ::  het_rates(:,:)'
                  end select
               end if
            end if
            if( model /= 'WRF' ) then
               write(30,100) trim(line)
	    end if
	 end if
         line = ' '
         write(30,100) trim(line)
	 if( clscnt(class) /= 0 ) then
            line = ' '
            write(30,100) trim(line)
            buff = ' '
            if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
!              if( march /= 'SCALAR' .and. class /= explicit ) then
               if( march == 'VECTOR' ) then
                  line = '!--------------------------------------------------------------------'
                  write(30,100) trim(line)
                  line = '!     ... local variables'
                  write(30,100) line
                  line = '!--------------------------------------------------------------------'
                  write(30,100) trim(line)
                  line = ' '
                  line(7:) = 'integer :: k'
                  write(30,100) trim(line)
	       end if
               line = ' '
               write(30,100) trim(line)
               line = ' '
	    end if
	    if( allocated( match_mask ) ) then
	       deallocate( match_mask )
	    end if
	    if( allocated( pmask ) ) then
	       deallocate( pmask )
	    end if
	    if( allocated( indexer ) ) then
	       deallocate( indexer )
	    end if
	    k = sum( cls_rxt_cnt(:,class) )
	    if( k == 0 ) then
		  if( class == implicit .or. class == rodas ) then
	             line(10:) = trim( l_piece ) // ':) = 0.' // trim(num_suffix)
                  else
                     if( model == 'CAM' ) then
                        if( march /= 'VECTOR' ) then
                           line(10:) = 'loss(:,:,:) = 0.' // trim(num_suffix)
                        end if
                     else
                        line(7:) = 'loss(:,:) = 0.' // trim(num_suffix)
                     end if
		  end if
	          write(30,100) trim(line)

		  if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
	             line(10:) = trim( p_piece ) // ':) = 0.' // trim(num_suffix)
	          else
	             line(7:) = 'prod(:,:) = 0.' // trim(num_suffix)
	          end if
		  write(30,100) trim(line)

               call terminate_subroutine
	       cycle Class_loop
	    end if
	    allocate( match_mask(k,3) )
	    allocate( indexer(k) )
	    if( sum( cls_rxt_cnt(2:3,class) ) /= 0 ) then
	       allocate( pmask(k,prd_lim) )
	    end if
            line = '!--------------------------------------------------------------------'
	    write(30,100) trim(line)
            line = '!       ... loss and production for'
	    length = len_trim( line ) + 2
	    select case( class )
	       case( 1 )
	          line(length:) = 'Explicit method'
	       case( 2 )
	          line(length:) = 'Ebi-gs method'
	       case( 3 )
	          line(length:) = 'Hov-gs method'
	       case( 4 )
	          line(length:) = 'Implicit method'
	       case( 5 )
	          line(length:) = 'Rodas3 method'
	    end select
	    write(30,100) trim(line)
            line = '!--------------------------------------------------------------------'
	    write(30,100) trim(line)
	    line = ' '
	    write(30,100) trim(line)
100   format(a)
	    if( class == ebi ) then
	       line = ' '
	       line(10:) = 'select case( index )'
	       write(30,100) trim(line)
	       line = ' '
	    else if( model == 'CAM' .or. class == 4 .or. class == 5 ) then
	       line = ' '
	       if( march /= 'SCALAR' ) then
	          if( march == 'VECTOR' ) then
	             line(7:) = 'do k = ofl,ofu'
		  else
                     if( model == 'MOZART' ) then
	                line(7:) = 'do k = 1,clsze'
                     else if( model == 'CAM' .and. class /= explicit ) then
	                line(7:) = 'do k = 1,cols'
		     end if
		  end if
	       end if
	       write(30,100) trim(line)
	       line = ' '
	    end if
Species_loop : &
            do species = 1,clscnt(class)
	       if( class == ebi .or. class == hov ) then
		  write(num,'(i3)') permute(species,2)
	          line = ' '
	          line(13:) = 'case( ' // num(:len_trim(num)) // ' )'
	          write(30,100) trim(line)
	          line = ' '
	       end if
!-----------------------------------------------------------------------
!   	...Write code for loss processes; linear, nonlinear, and heterogeneous
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!   	... Setup indicies and check whether species is in any loss reactions
!-----------------------------------------------------------------------
	       target = clsmap(species,class,2)
	       kl = cls_rxt_cnt(1,class) + 1
	       ku = sum( cls_rxt_cnt(:3,class) )
	       do i = 1,2
	          match_mask(kl:ku,i) = cls_rxt_map(kl:ku,i+1,class) == target
	          where( match_mask(kl:ku,i) )
	             indexer(kl:ku) = 6/(i+1)
	          endwhere
	       end do
	       match_mask(kl:ku,1) = match_mask(kl:ku,1) .or. match_mask(kl:ku,2)
	       kl = ku + 1 ; ku = sum(cls_rxt_cnt(:,class))
	       match_mask(kl:ku,1) = cls_rxt_map(kl:ku,2,class) == species
	       kl = cls_rxt_cnt(1,class) + 1
	       if( class == ebi .or. class == hov ) then
		  num = '1'
	       else
	          write(num,'(i3)') permute(species,class)
	          num =  adjustl( num )
	       end if
	       if( count( match_mask(kl:ku,1) ) == 0 ) then
!-----------------------------------------------------------------------
!   	... no loss for this species
!-----------------------------------------------------------------------
		  if( class == implicit .or. class == rodas ) then
	             line(10:) = trim( l_piece ) // num(:len_trim(num)) // ') = 0.' // trim(num_suffix)
                  else
                     if( model == 'CAM' ) then
                        if( march /= 'VECTOR' ) then
                           line(10:) = 'loss(:,:,' // num(:len_trim(num)) // ') = 0.' // trim(num_suffix)
                        end if
                     else
                        line(7:) = 'loss(:,' // num(:len_trim(num)) // ') = 0.' // trim(num_suffix)
                     end if
		  end if
	          write(30,100) trim(line)
	       else
	          line = ' '
	          if( class == explicit .or. class >= implicit ) then
		     if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
	                line(10:) = trim( l_piece ) // num(:len_trim(num)) // ') = ('
                     else
                        line(7:) = 'loss(:,' // num(:len_trim(num)) // ') = ('
		     end if
	             line_pos = len_trim( line ) + 1
	          else
	             line(7:) = 'loss(:,' // num(:len_trim(num)) // ') ='
	             line_pos = len_trim( line ) + 2
	          end if
!-----------------------------------------------------------------------
!	... Scan loss reactions for common terms
!-----------------------------------------------------------------------
	          ku = ku - cls_rxt_cnt(4,class)
		  first = .true.
	          do m = 1,spccnt
		     match_mask(kl:ku,2) = .false.
		     do k = kl,ku
			if( match_mask(k,1) ) then
			   if( abs( cls_rxt_map(k,indexer(k),class) ) == m ) then
			      match_mask(k,2) = .true.
			   end if
			end if
		     end do
		     cnt = count( match_mask(kl:ku,2) )
		     if( cnt == 0 ) then
		        cycle
		     end if
		     if( m == target ) then
		        if( cnt > 1 ) then
			   if( first ) then
		              buff = '2.' // trim(num_suffix) // '*('
		           else
		              buff = ' + 2.' // trim(num_suffix) // '*('
			   end if
		        else
			   if( first ) then
		              buff = '2.' // trim(num_suffix) // '*'
		           else
		              buff = ' + 2.' // trim(num_suffix) // '*'
			   end if
		        end if
		     else
		        if( cnt > 1 ) then
			   if( first ) then
		              buff = '('
		           else
		              buff = ' + ('
			   end if
		        else
			   if( first ) then
		              buff = ' '
		           else
		              buff = ' + '
			   end if
		        end if
		     end if
		     if( first ) then
			first = .false.
		     end if
		     l = 0
		     do k = kl,ku
		        if( match_mask(k,2) ) then
			   l = l + 1
	                   write(num,'(i3)') cls_rxt_map(k,1,class)
	                   num =  adjustl( num )
			   if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
		              buff(len_trim(buff)+1:) = trim( rxt_piece ) // num(:len_trim(num)) // ') +'
		           else
		              buff(len_trim(buff)+1:) = 'rxt(:,' // num(:len_trim(num)) // ') +'
		           end if
			   if( l == cnt ) then
		              if( cnt > 1 ) then
			         buff(len_trim(buff)-1:) = ')'
			      else
			         buff(len_trim(buff)-1:) = ' '
			      end if
			      call put_in_line
			      write(num,'(i3)') m
	                      num =  adjustl( num )
			      if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
		                 buff(len_trim(buff)+1:) = '* ' // trim( sol_piece ) // num(:len_trim(num)) // ')'
			      else
		                 buff(len_trim(buff)+1:) = '* y(:,' // num(:len_trim(num)) // ')'
			      end if
			   end if
			   call put_in_line
		        end if
		     end do
		     where( match_mask(kl:ku,2) )
		        match_mask(kl:ku,1) = .false.
		     endwhere
	          end do
!-----------------------------------------------------------------------
!	... Strictly unimolecular losses
!-----------------------------------------------------------------------
	          ku = sum( cls_rxt_cnt(:,class) )
	          cnt = count( match_mask(kl:ku,1) )
		  if( cnt > 0 ) then
	             do k = kl,ku
		        if( match_mask(k,1) ) then
		           cnt = cnt - 1
	                   write(num,'(i3)') cls_rxt_map(k,1,class)
	                   num =  adjustl( num )
			   if( k <= sum(cls_rxt_cnt(1:3,class)) ) then
			      if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
		                 buff(len_trim(buff)+1:) = ' + ' // trim( rxt_piece ) // num(:len_trim(num)) // ')'
		              else
		                 buff(len_trim(buff)+1:) = ' + rxt(:,' // num(:len_trim(num)) // ')'
		              end if
		           else
			      if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
		                 buff(len_trim(buff)+1:) = ' + ' // trim( het_piece ) // num(:len_trim(num)) // ')'
		              else
		                 buff(len_trim(buff)+1:) = ' + het_rates(:,' // num(:len_trim(num)) // ')'
		              end if
		           end if
	                   if( cnt == 0 .and. (class == explicit .or. class >= implicit) ) then
		              buff(len_trim(buff)+1:) = ')'
		           end if
		           call put_in_line
		        end if
	             end do
		  else if( class == explicit .or. class >= implicit ) then
		     buff(len_trim(buff)+1:) = ')'
		     call put_in_line
		  end if
		  if( class == explicit .or. class >= implicit ) then
	             write(num,'(i3)') target
	             num =  adjustl( num )
		     if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
		        buff(len_trim(buff)+1:) = '* ' // trim( sol_piece ) // num(:len_trim(num)) // ')'
	             else
		        buff(len_trim(buff)+1:) = '* y(:,' // num(:len_trim(num)) // ')'
		     end if
		     call put_in_line
		  end if
	          if( line(7:) /= ' ' ) then
		     write(30,100) trim(line)
	          end if
	       end if
!-----------------------------------------------------------------------
!   	...Write code for production from linear and nonlinear reactions
!-----------------------------------------------------------------------
	       ku = sum( cls_rxt_cnt(:3,class) )
	       do k = kl,ku
		  pmask(k,:) = cls_rxt_map(k,4:prd_lim+3,class) == species
	          match_mask(k,1) = ANY( pmask(k,:) )
	       end do
	       if( class == ebi .or. class == hov ) then
		  num = '1'
	       else
	          write(num,'(i3)') permute(species,class)
	          num =  adjustl( num )
	       end if
	       line = ' '
!-----------------------------------------------------------------------
!	... No species products
!-----------------------------------------------------------------------
	       if( count( match_mask(kl:ku,1) ) == 0 ) then
		  if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
	             line(10:) = trim( p_piece ) // num(:len_trim(num)) // ') = 0.' // trim(num_suffix)
	          else
	             line(7:) = 'prod(:,' // num(:len_trim(num)) // ') = 0.' // trim(num_suffix)
	          end if
		  write(30,100) trim(line)
		  cycle Species_loop
	       else
		  if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
	             line(10:) = trim( p_piece ) // num(:len_trim(num)) // ') = '
	          else
	             line(7:) = 'prod(:,' // num(:len_trim(num)) // ') = '
	          end if
	       end if
	       first = .true.
	       do
	          do m = 1,spccnt
		     match_mask(kl:ku,3) = (abs(cls_rxt_map(kl:ku,2,class)) == m .or. &
		                           abs(cls_rxt_map(kl:ku,3,class)) == m) .and.&
					   match_mask(kl:ku,1)
		     freq(m) = count( match_mask(kl:ku,3) )
	          end do
		  max_loc = MAXLOC( freq(:spccnt) )
		  cnt = MAXVAL( freq(:spccnt) )
		  match_mask(kl:ku,3) = (abs(cls_rxt_map(kl:ku,2,class)) == max_loc(1) .or. &
		                        abs(cls_rxt_map(kl:ku,3,class)) == max_loc(1)) .and. &
					match_mask(kl:ku,1)
		  do k = kl,ku
		     if( match_mask(k,3) ) then
			if( abs( cls_rxt_map(k,2,class) ) == max_loc(1) ) then
			   indexer(k) = 3
			else
			   indexer(k) = 2
			end if
		     end if
		  end do
		  if( cnt > 1 ) then
		     if( first ) then
		        buff = ' ('
		     else
		        buff = ' + ('
		     end if
		  else if( first ) then
		     buff = ' '
		  else
		     buff = ' +'
		  end if
		  if( first ) then
		     first = .false.
		  end if
		  m = cnt
		  do k = kl,ku
		     if( match_mask(k,3) ) then
		        index = pcoeff_ind(cls_rxt_map(k,1,class))
			if( index /= 0 ) then
			   rate = 0.
			   do ml = 1,prd_lim
			      if( pmask(k,ml) ) then
			         rate = rate + pcoeff(ml,index)
			      end if
			   end do
			else
			   rate = REAL( count( abs( cls_rxt_map(k,4:prd_lim+3,class) ) == species ) )
		        end if
		        if( rate /= 0. .and. rate /= 1. ) then
			   call r2c( buff(len_trim(buff)+1:), rate, 'l' )
			   buff(len_trim( buff )+1:) = trim(num_suffix) // '*'
		        end if
	                write(num,'(i3)') cls_rxt_map(k,1,class)
	                num =  adjustl( num )
			if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
		           buff(len_trim(buff)+1:) = trim( rxt_piece ) // num(:len_trim(num)) // ')'
			else
		           buff(len_trim(buff)+1:) = 'rxt(:,' // num(:len_trim(num)) // ')'
			end if
			if( abs( cls_rxt_map(k,indexer(k),class) ) /= 0 ) then
	                   write(num,'(i3)') abs( cls_rxt_map(k,indexer(k),class) )
	                   num =  adjustl( num )
			   if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
			      if( m > 1 ) then
		                 buff(len_trim(buff)+1:) = '*' // trim( sol_piece ) // num(:len_trim(num)) // ') +'
			      else if( cnt > 1 ) then
		                 buff(len_trim(buff)+1:) = '*' // trim( sol_piece ) // num(:len_trim(num)) // '))'
			      else
		                 buff(len_trim(buff)+1:) = '*' // trim( sol_piece ) // num(:len_trim(num)) // ')'
			      end if
			   else
			      if( m > 1 ) then
		                 buff(len_trim(buff)+1:) = '*y(:,' // num(:len_trim(num)) // ') +'
			      else if( cnt > 1 ) then
		                 buff(len_trim(buff)+1:) = '*y(:,' // num(:len_trim(num)) // '))'
			      else
		                 buff(len_trim(buff)+1:) = '*y(:,' // num(:len_trim(num)) // ')'
			      end if
			   end if
			else
			   if( m > 1 ) then
		              buff(len_trim(buff)+1:) = ' +'
			   else if( cnt > 1 ) then
		              buff(len_trim(buff)+1:) = ')'
			   end if
			end if
			call put_in_line
			if( m == 1 ) then
	                   write(num,'(i3)') max_loc(1)
	                   num =  adjustl( num )
			   if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
		              buff = '*' // trim( sol_piece ) // num(:len_trim(num)) // ')'
			   else
		              buff = '*y(:,' // num(:len_trim(num)) // ')'
			   end if
			   call put_in_line
			   exit
			end if
			m = m - 1
		     end if
		  end do
		  where( match_mask(kl:ku,3) )
		     match_mask(kl:ku,1) = .false.
		  endwhere
		  if( count( match_mask(kl:ku,1) ) == 0 ) then
		     exit
		  end if
	       end do
	       if( line(7:) /= ' ' ) then
		  write(30,100) trim(line)
		  line = ' '
	       end if
	    end do Species_loop
	    if( class == ebi ) then
	       line = ' '
	       line(10:) = 'end select'
	       write(30,100) trim(line)
	       line = ' '
	    end if
	    if( model == 'CAM' .or. class == implicit .or. class == rodas ) then
!              if( march /= 'SCALAR' .and. class /= explicit ) then
               if( march == 'VECTOR' ) then
	          line = '      end do'
	          write(30,100) trim(line)
	       end if
	    end if
	 end if
!-----------------------------------------------------------------------
!	... Terminate the subroutine
!-----------------------------------------------------------------------
         call terminate_subroutine
      end do Class_loop

      line = '      end module mo_prod_loss'
      write(30,100) trim(line)


      CLOSE( 30 )
      if( allocated( match_mask ) ) then
	 deallocate( match_mask )
      end if
      if( allocated( pmask ) ) then
	 deallocate( pmask )
      end if
      if( allocated( indexer ) ) then
	 deallocate( indexer )
      end if

      contains

      subroutine put_in_line
!-----------------------------------------------------------------------
!	... Put line piece in buff into the line
!-----------------------------------------------------------------------

      implicit none

!-----------------------------------------------------------------------
!	... Local variables
!-----------------------------------------------------------------------
      integer :: blen, llen

      blen = len_trim( buff )
      llen = len_trim( line ) + 1
      if( blen + llen < max_len-2 ) then
	 line(llen:) = buff(:blen)
      else
	 line(len_trim(line)+1:) = ' &'
	 write(30,'(a)') trim(line)
	 line = ' '
	 line(18:) = buff(:length)
      end if
      buff = ' '

      end subroutine put_in_line

      subroutine terminate_subroutine
!-----------------------------------------------------------------------
!	... Terminate the subroutine
!-----------------------------------------------------------------------

      implicit none

      line = ' '
      write(30,100) trim(line)
      select case( class )
	 case( 1 )
            line = '      end subroutine exp_prod_loss'
	 case( 2 )
            line = '      end subroutine ebi_prod_loss'
	 case( 3 )
            line = '      end subroutine hov_prod_loss'
	 case( 4 )
            line = '      end subroutine imp_prod_loss'
	 case( 5 )
            line = '      end subroutine rodas_prod_loss'
      end select
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)

100   format(a)

      end subroutine terminate_subroutine

      end subroutine pl_code
       
      end module prod_loss
