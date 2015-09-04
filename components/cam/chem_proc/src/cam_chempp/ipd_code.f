
      module ind_prod

      use io, only : temp_path

      contains

      subroutine ipd_code( spccnt, &
                           clscnt, &
                           clsmap, &
                           cls_rxt_cnt, &
                           extcnt, &
                           cls_rxt_map, &
                           pcoeff_ind, &
                           pcoeff, &
                           permute, &
                           model, &
                           march )
     
      use var_mod, only : var_lim
      use rxt_mod, only : rxt_lim, prd_lim

      implicit none

!-----------------------------------------------------------------------
!        ... Dummy args
!-----------------------------------------------------------------------
      integer, intent(in) ::      spccnt
      integer, intent(in) ::      clscnt(5), &
                                  extcnt(5), &
                                  clsmap(var_lim,5,2), &
                                  cls_rxt_map(rxt_lim,prd_lim+3,5), &
                                  cls_rxt_cnt(4,5)
      integer, intent(in) ::      permute(var_lim,5)
      integer, intent(in) ::      pcoeff_ind(*)
      real, intent(in)    ::      pcoeff(prd_lim,*)
      character(len=*), intent(in) :: model
      character(len=*), intent(in) :: march
      
!-----------------------------------------------------------------------
!        ... Local variables
!-----------------------------------------------------------------------
      integer, parameter :: max_len = 90
      integer  ::   i, k, kl, ku, l, m, n, prdndx
      integer  ::   length, index
      integer  ::   line_pos, cnt
      integer  ::   class
      integer  ::   base
      integer  ::   species
      integer  ::   match_cnt
      integer  ::   max_loc(1)
      integer  ::   match_ind(4)
      integer  ::   freq(spccnt)
      integer, allocatable :: indexer(:)
      real     ::   rate
      character(len=max_len) :: line
      character(len=max_len) :: buff
      character(len=4) :: num_suffix 
      character(len=4) :: dec_suffix
      character(len=3) ::  num
      logical  ::  beg_line
      logical  ::  lexist, first, indprds
      logical  ::  first_class = .true.
      logical  ::  is_vector
      logical, allocatable, dimension(:,:) :: match_mask, pmask
      
      inquire( file = trim( temp_path ) // 'indprd.F', exist = lexist )
      if( lexist ) then
         call system( 'rm ' // trim( temp_path ) // 'indprd.F' )
      end if
      open( unit = 30, file = trim( temp_path ) // 'indprd.F' )

      if( model == 'CAM' ) then
         num_suffix = '_r8'
         dec_suffix = '(r8)'
      else
         num_suffix = ' '
         dec_suffix = ' '
      end if

      is_vector = march == 'VECTOR'

      line = ' '
      write(30,100) trim(line)
      line = '      module mo_indprd'
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
      line = '      public :: indprd'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '      contains'
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      if( model /= 'CAM' ) then
         line = '      subroutine indprd( class, prod, y, extfrc, rxt )'
      else
         if( .not. is_vector ) then
            line = '      subroutine indprd( class, prod, nprod, y, extfrc, rxt, ncol )'
         else
            line = '      subroutine indprd( class, prod, y, extfrc, rxt )'
         end if
         write(30,100) trim(line)
         line = ' '
         write(30,100) trim(line)
         line = '      use chem_mods, only : gas_pcnst, extcnt, rxntot'
         if( .not. is_vector ) then
            write(30,100) trim(line)
            line = '      use ppgrid,    only : pver'
         end if
      end if
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '      implicit none '
      write(30,100) trim(line)
      line = ' '
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      line = '!       ... dummy arguments'
      write(30,100) trim(line)
      line = '!--------------------------------------------------------------------'
      write(30,100) trim(line)
      line = '      integer, intent(in) :: class'
      write(30,100) trim(line)
      if( model /= 'CAM' ) then
         line = '      real' // trim(dec_suffix) // ', intent(in)    :: y(:,:)'
         write(30,100) trim(line)
         line = '      real' // trim(dec_suffix) // ', intent(in)    :: rxt(:,:)'
         write(30,100) trim(line)
         line = '      real' // trim(dec_suffix) // ', intent(in)    :: extfrc(:,:)'
         write(30,100) trim(line)
         line = '      real' // trim(dec_suffix) // ', intent(inout) :: prod(:,:)'
         write(30,100) trim(line)
      else
         if( .not. is_vector ) then
            line = '      integer, intent(in) :: ncol'
            write(30,100) trim(line)
            line = '      integer, intent(in) :: nprod'
            write(30,100) trim(line)
            line = '      real' // trim(dec_suffix) // ', intent(in)    :: y(ncol,pver,gas_pcnst)'
            write(30,100) trim(line)
            line = '      real' // trim(dec_suffix) // ', intent(in)    :: rxt(ncol,pver,rxntot)'
            write(30,100) trim(line)
            line = '      real' // trim(dec_suffix) // ', intent(in)    :: extfrc(ncol,pver,extcnt)'
            write(30,100) trim(line)
            line = '      real' // trim(dec_suffix) // ', intent(inout) :: prod(ncol,pver,nprod)'
            write(30,100) trim(line)
         else
            line = '      real' // trim(dec_suffix) // ', intent(in)    :: y(:,:)'
            write(30,100) trim(line)
            line = '      real' // trim(dec_suffix) // ', intent(in)    :: rxt(:,:)'
            write(30,100) trim(line)
            line = '      real' // trim(dec_suffix) // ', intent(in)    :: extfrc(:,:)'
            write(30,100) trim(line)
            line = '      real' // trim(dec_suffix) // ', intent(inout) :: prod(:,:)'
            write(30,100) trim(line)
         end if
      end if
      line = ' '
      write(30,100) trim(line)
      line = ' '

Class_loop : &
      do class = 1,5
         if( clscnt(class) /= 0 ) then
            if( allocated( match_mask ) ) then
               deallocate( match_mask )
            end if
            if( allocated( pmask ) ) then
               deallocate( pmask )
            end if
            if( allocated( indexer ) ) then
               deallocate( indexer )
            end if
            line = '!--------------------------------------------------------------------'
            write(30,100) trim(line)
            line = '!       ... "independent" production for'
            length = len_trim( line ) + 2
            if( class == 1 ) then
               line(length:) = 'Explicit species'
            else if( class == 2 ) then
               line(length:) = 'Ebi species'
            else if( class == 3 ) then
               line(length:) = 'Hov species'
             else if( class == 4 ) then
               line(length:) = 'Implicit species'
             else if( class == 5 ) then
               line(length:) = 'Rodas species'
             end if
             write(30,100) trim(line)
            line = '!--------------------------------------------------------------------'
	    write(30,100) trim(line)
	    if( first_class ) then
	       line = '      if( class =='
	       first_class = .false.
	    else
	       line = '      else if( class =='
	    end if
	    write(line(len_trim(line)+2:),'(i1)') class
	    line(len_trim(line)+2:) = ') then'
	    write(30,100) trim(line)
100   format(a)
	    ku = MAX( cls_rxt_cnt(1,class),extcnt(class) )
	    if( ku == 0 ) then
               do species = 1,clscnt(class)
	          write(num,'(i3)') species
	          num =  adjustl( num )
	          line = ' '
                  if( model /= 'CAM' ) then
	             line(10:) = 'prod(:,' // num(:len_trim(num)) // ') = 0.'
                  else
                     if( .not. is_vector ) then
	                line(10:) = 'prod(:,:,' // num(:len_trim(num)) // ') = 0._r8'
                     else
                        line(10:) = 'prod(:,' // num(:len_trim(num)) // ') = 0._r8'
                     end if
                  end if
                  write(30,100) trim(line)
               end do
               cycle Class_loop
            end if
            kl = 1
            allocate( match_mask(ku,3) )
            allocate( indexer(ku) )
            allocate( pmask(ku,prd_lim) )

Species_loop : &
            do species = 1,clscnt(class)
               line = ' '
	       write(num,'(i3)') permute(species,class)
	       num =  adjustl( num )
               if( model /= 'CAM' ) then
	          line(10:) = 'prod(:,' // num(:len_trim(num)) // ') = '
               else
                  if( .not. is_vector ) then
	             line(10:) = 'prod(:,:,' // num(:len_trim(num)) // ') = '
                  else
	             line(10:) = 'prod(:,' // num(:len_trim(num)) // ') = '
                  end if
               end if
               ku = cls_rxt_cnt(1,class)
!-----------------------------------------------------------------------
!   	...Write code for "independent" production processes
!-----------------------------------------------------------------------
	       do k = kl,ku
		  pmask(k,:) = cls_rxt_map(k,4:prd_lim+3,class) == species
	          match_mask(k,1) = any( pmask(k,:) )
	       end do
!-----------------------------------------------------------------------
!	... No species products
!-----------------------------------------------------------------------
	       if( count( match_mask(kl:ku,1) ) /= 0 ) then
		  indprds = .true.
	          first = .true.
	          do
	             do m = 1,spccnt
		        match_mask(kl:ku,3) = (abs(cls_rxt_map(kl:ku,2,class)) == m .or. &
		                              abs(cls_rxt_map(kl:ku,3,class)) == m) .and.&
                                                  match_mask(kl:ku,1)
		        freq(m) = count( match_mask(kl:ku,3) )
	             end do
		     max_loc = maxloc( freq(:spccnt) )
		     cnt     = maxvaL( freq(:spccnt) )
		     if( cnt /= 0 ) then
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
		     else
		        match_mask(kl:ku,3) = match_mask(kl:ku,1)
			indexer(kl:ku) = 0
			cnt = count( match_mask(kl:ku,3) )
		     end if
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
                           if( index == 0 ) then
                              rate = 1.
                           else
                              rate = 0.
                              do prdndx = 1,prd_lim
                                 if( pmask(k,prdndx) ) then
                                    rate = rate + pcoeff(prdndx,index)
                                  end if
                              end do
                           end if
                           if( rate /= 0. .and. rate /= 1. ) then
                              call r2c( buff(len_trim(buff)+1:), rate, 'l' )
                              buff(len_trim( buff )+1:) = trim(num_suffix) // '*'
                           end if
                           write(num,'(i3)') cls_rxt_map(k,1,class)
	                   num =  adjustl( num )
		           if( model /= 'CAM' ) then
		              buff(len_trim(buff)+1:) = 'rxt(:,' // num(:len_trim(num)) // ')'
		           else
		              if( .not. is_vector ) then
		                 buff(len_trim(buff)+1:) = 'rxt(:,:,' // num(:len_trim(num)) // ')'
		              else
		                 buff(len_trim(buff)+1:) = 'rxt(:,' // num(:len_trim(num)) // ')'
		              end if
		           end if
			   if( indexer(k) /= 0 ) then
                              if( abs( cls_rxt_map(k,indexer(k),class) ) /= 0 ) then
                                 write(num,'(i3)') abs( cls_rxt_map(k,indexer(k),class) )
                                 num =  adjustl( num )
		                 if( model /= 'CAM' ) then
                                    if( m > 1 ) then
                                       buff(len_trim(buff)+1:) = '*y(:,' // num(:len_trim(num)) // ') +'
                                    else if( cnt > 1 ) then
                                       buff(len_trim(buff)+1:) = '*y(:,' // num(:len_trim(num)) // '))'
                                    else
                                       buff(len_trim(buff)+1:) = '*y(:,' // num(:len_trim(num)) // ')'
                                    end if
		                 else
                                    if( .not. is_vector ) then
                                       if( m > 1 ) then
                                          buff(len_trim(buff)+1:) = '*y(:,:,' // num(:len_trim(num)) // ') +'
                                       else if( cnt > 1 ) then
                                          buff(len_trim(buff)+1:) = '*y(:,:,' // num(:len_trim(num)) // '))'
                                       else
                                          buff(len_trim(buff)+1:) = '*y(:,:,' // num(:len_trim(num)) // ')'
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
                                 end if
                              else
                                 if( m > 1 ) then
                                    buff(len_trim(buff)+1:) = ' +'
                                 else if( cnt > 1 ) then
                                    buff(len_trim(buff)+1:) = ')'
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
			      if( indexer(k) /= 0 ) then
	                         write(num,'(i3)') max_loc(1)
	                         num =  adjustl( num )
		                 if( model /= 'CAM' ) then
		                    buff = '*y(:,' // num(:len_trim(num)) // ')'
		                 else
		                    if( .not. is_vector ) then
		                       buff = '*y(:,:,' // num(:len_trim(num)) // ')'
		                    else
		                       buff = '*y(:,' // num(:len_trim(num)) // ')'
                                    end if
			         end if
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
	       else
	          indprds = .false.
	       end if
!-----------------------------------------------------------------------
!   	... Write code for "extraneous" production processes
!-----------------------------------------------------------------------
	       base = sum( cls_rxt_cnt(1:4,class) )
	       match_mask(:,2) = .false.
	       match_mask(:extcnt(class),2) = cls_rxt_map(base+1:base+extcnt(class),2,class) == species
	       if( count( match_mask(:extcnt(class),2) ) /= 0 ) then
	          do k = base+1,base+extcnt(class)
		     if( cls_rxt_map(k,2,class) == species ) then
                        write(num,'(i3)') cls_rxt_map(k,1,class)
		        num = adjustl( num )
		        n = len_trim( num )
		        if( model /= 'CAM' ) then
		           buff = ' + extfrc(:,' // num(:n) // ')'
		        else
		           if( .not. is_vector ) then
                              buff = ' + extfrc(:,:,' // num(:n) // ')'
		           else
                              buff = ' + extfrc(:,' // num(:n) // ')'
		           end if
		        end if
		        call put_in_line
		     end if
	          end do
	       else if( .not. indprds ) then
		  buff = ' 0._r8'
		  call put_in_line
	       end if
	       if( line /= ' ' ) then
		  write(30,100) trim(line)
	       end if
	       line = ' '
	       write(30,100) line
	    end do Species_loop
	 end if
      end do Class_loop
      
      if ( .not. first_class ) then
         line = '      end if'
         write(30,100) line
      endif

      line = ' '
      write(30,100) line
      line = '      end subroutine indprd'
      write(30,100) line
      line = ' '
      write(30,100) line
      line = '      end module mo_indprd'
      write(30,100) line

      if( allocated( match_mask ) ) then
         deallocate( match_mask )
      end if
      if( allocated( pmask ) ) then
	 deallocate( pmask )
      end if
      if( allocated( indexer ) ) then
	 deallocate( indexer )
      end if
      close( 30 )

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
	 line(18:) = buff(:blen)
      end if
      buff = ' '

      end subroutine put_in_line

      end subroutine ipd_code

      end module ind_prod
