module set_rxt_rates

  use rxt_mod, only : rxntot, cph_flg
  use io,      only : temp_path

  private
  public :: make_rate

contains

  subroutine make_rate( sym_rates, rxptab, rxpcnt, machine, vec_ftns, &
       model, march )
    !-----------------------------------------------------------------------
    !        ... write fortran "internal" reaction rates
    !-----------------------------------------------------------------------

    use rxt_mod, only : troecnt, troetab, troe_sym_rates, rxparm

    implicit none

    !-----------------------------------------------------------------------
    !        ... dummy arguments
    !-----------------------------------------------------------------------
    integer, intent(in)           ::    rxpcnt
    integer, intent(in)           ::    rxptab(*)
    character(len=16), intent(in) :: sym_rates(2,*)
    character(len=16), intent(in)  :: machine
    character(len=16), intent(in)  :: march
    character(len=16), intent(in)  :: model
    logical, intent(in)           :: vec_ftns

    !-----------------------------------------------------------------------
    !        ... local variables
    !-----------------------------------------------------------------------
    integer  ::   i, cnt, indp, inde, l, m, m1, pos
    integer  ::   ftn_lim
    integer  ::   ftn_mode
    integer  ::   cph_cnt
    integer  ::   match_cnt
    integer  ::   subs_lim
    integer  ::   subs
    integer  ::   match_ind(rxpcnt)
    real     ::   rate
    character(len=132) :: line
    character(len=32)  :: wrk, buff
    character(len=16)   :: vec_dim
    character(len=4)   :: num_suffix 
    character(len=4)   :: dec_suffix
    character(len=3)   :: num
    character(len=3)   :: numa
    logical  ::  lexist
    logical  ::  vftns
    logical  ::  do_tmp
    logical  ::  do_rxt
    logical  ::  t_dependent(rxpcnt)



    inquire( file = trim( temp_path ) // 'mo_setrxt.F', exist = lexist )
    if( lexist ) then
       call system( 'rm ' // trim( temp_path ) // 'mo_setrxt.F' )
    end if
    open( unit = 30, file = trim( temp_path ) // 'mo_setrxt.F' )

    vec_dim = 'plnplv'
    if( model == 'CAM' ) then
       num_suffix = '_r8'
       dec_suffix = '(r8)'
       if( march == 'VECTOR' ) then
          vec_dim = 'chnkpnts'
       end if
    else
       num_suffix = ' '
       dec_suffix = ' '
    end if

    if( model == 'CAM' ) then
       subs_lim = 2
    else
       subs_lim = 1
    end if
    cph_cnt = count( cph_flg(:) )

    line = ' '
    write(30,100) trim(line)
    line(7:) = 'module mo_setrxt'
    write(30,100) trim(line)
    line = ' '
    write(30,100) trim(line)
    line = '      use shr_kind_mod, only : r8 => shr_kind_r8'
    write(30,100) trim(line)
    line = ' '
    write(30,100) trim(line)
    line(7:) = 'private'
    write(30,100) trim(line)
    line(7:) = 'public :: setrxt'
    write(30,100) trim(line)
    line(7:) = 'public :: setrxt_hrates'
    write(30,100) trim(line)
    line = ' '
    write(30,100) trim(line)
    line(7:) = 'contains'
    write(30,100) trim(line)
    line = ' '
    subs_loop : do subs = 1,subs_lim
       write(30,100) trim(line)
       if( march == 'VECTOR' ) then
          if( subs == 1 ) then
             line(7:) = 'subroutine setrxt( rate, temp, m, ' // trim(vec_dim) // ' )'
          else
             line(7:) = 'subroutine setrxt_hrates( rate, temp, m, ' // trim(vec_dim) // ' )'
          end if
          write(30,100) trim(line)
          line = ' '
          write(30,100) trim(line)
       else
          select case( model )
          case( 'MOZART' )
             line(7:) = 'subroutine setrxt( rate, temp, m, plonl )'
             write(30,100) trim(line)
             line = ' '
             write(30,100) trim(line)
             line = ' '
             line(7:) = 'use mo_grid,   only : plev, plnplv'
          case( 'CAM' )
             if( subs == 1 ) then
                line(7:) = 'subroutine setrxt( rate, temp, m, ncol )'
             else
                line(7:) = 'subroutine setrxt_hrates( rate, temp, m, ncol, kbot )'
             end if
             write(30,100) trim(line)
             line = ' '
             write(30,100) trim(line)
             line = ' '
             line(7:) = 'use ppgrid,       only : pver, pcols'
             write(30,100) trim(line)
             line(7:) = 'use shr_kind_mod, only : r8 => shr_kind_r8'
          case( 'WRF' )
             line(7:) = 'subroutine setrxt( rate, temp, m, n )'
             write(30,100) trim(line)
             line = ' '
             write(30,100) trim(line)
             line(7:) = 'use mo_jpl, only : jpl'
          end select
       end if
       write(30,100) trim(line)
       if( model /= 'WRF' ) then
          line(7:) = 'use chem_mods, only : rxntot'
          write(30,100) trim(line)
          line(7:) = 'use mo_jpl,    only : jpl'
          write(30,100) trim(line)
       end if
       line = ' '
       write(30,100) trim(line)
       line(7:) = 'implicit none '
       write(30,100) trim(line)
       line = ' '
       write(30,100) trim(line)
       line = '!-------------------------------------------------------'
       write(30,100) trim(line)
       line = '!       ... dummy arguments'
       write(30,100) trim(line)
       line = '!-------------------------------------------------------'
       write(30,100) trim(line)
       if( machine == 'NEC' .or. machine == 'FUJITSU' ) then
          line = '      integer, intent(in) :: ' // trim(vec_dim)
          write(30,100) trim(line)
          line = '      real' // trim(dec_suffix) // ', intent(in)    :: temp(' // trim(vec_dim) // ')'
          write(30,100) trim(line)
          line = '      real' // trim(dec_suffix) // ', intent(in)    :: m(' // trim(vec_dim) // ')'
          write(30,100) trim(line)
          line = '      real' // trim(dec_suffix) // ', intent(inout) :: rate(' // trim(vec_dim) // ',max(1,rxntot))'
       else
          select case( model )
          case( 'MOZART' )
             line = '      integer, intent(in) :: plonl'
             write(30,100) trim(line)
             line = '      real' // trim(dec_suffix) // ', intent(in)    :: temp(plonl,plev)'
             write(30,100) trim(line)
             line = '      real' // trim(dec_suffix) // ', intent(in)    :: m(plonl,plev)'
             write(30,100) trim(line)
             line = '      real' // trim(dec_suffix) // ', intent(inout) :: rate(plonl,plev,rxntot)'
          case( 'CAM' )
             line = '      integer, intent(in) :: ncol'
             write(30,100) trim(line)
             if( subs == 2 ) then
                line = '      integer, intent(in) :: kbot'
                write(30,100) trim(line)
             end if
             line = '      real' // trim(dec_suffix) // ', intent(in)    :: temp(pcols,pver)'
             write(30,100) trim(line)
             line = '      real' // trim(dec_suffix) // ', intent(in)    :: m(ncol,pver)'
             write(30,100) trim(line)
             line = '      real' // trim(dec_suffix) // ', intent(inout) :: rate(ncol,pver,rxntot)'
          case( 'WRF' )
             line = '      integer, intent(in)  :: n'
             write(30,100) trim(line)
             line = '      real' // trim(dec_suffix) // ', intent(in)    :: temp(:)'
             write(30,100) trim(line)
             line = '      real' // trim(dec_suffix) // ', intent(in)    :: m(:)'
             write(30,100) trim(line)
             line = '      real' // trim(dec_suffix) // ', intent(inout) :: rate(:,:)'
          end select
       end if
       write(30,100) trim(line)
       line = ' '
       write(30,100) trim(line)

       if( rxpcnt == 0 ) then
          if( subs == 1 ) then
             line(7:) = 'end subroutine setrxt'
          else
             line(7:) = 'end subroutine setrxt_hrates'
          end if
          write(30,100) trim(line)
          line = ' '
          write(30,100) trim(line)
          if( subs == 2 ) then
             line(7:) = 'end module mo_setrxt'
             write(30,100) trim(line)
          end if
          cnt = 0
          cycle subs_loop
       else
          t_dependent(:rxpcnt) = sym_rates(2,:rxpcnt) /= ' '
          cnt = count( t_dependent(:rxpcnt) )
       end if


       !-----------------------------------------------------------------------
       !        ... check for temp dependent rates
       !-----------------------------------------------------------------------
       temp_dep_rxts : if( cnt /= 0 .or. troecnt /= 0 ) then
          line = '!-------------------------------------------------------'
          write(30,100) trim(line)
          line = '!       ... local variables'
          write(30,100) trim(line)
          line = '!-------------------------------------------------------'
          write(30,100) trim(line)
          line = ' '
          if( model /= 'WRF' ) then
             line(7:) = 'integer  ::  n'
             write(30,100) trim(line)
          end if
          if( machine == 'NEC' .or. machine == 'FUJITSU' ) then
             line(7:) = 'real' // trim(dec_suffix) // '  ::  itemp(' // trim(vec_dim) // ')'
             write(30,100) trim(line)
             line(7:) = 'real' // trim(dec_suffix) // '  ::  exp_fac(' // trim(vec_dim) // ')'
          else
             select case( model )
             case( 'MOZART' )
                line(7:) = 'real' // trim(dec_suffix) // '  ::  itemp(plonl,plev)'
                write(30,100) trim(line)
                line(7:) = 'real' // trim(dec_suffix) // '  ::  exp_fac(plonl,plev)'
             case( 'CAM' )
                if( subs == 1 ) then
                   line(7:) = 'real' // trim(dec_suffix) // '  ::  itemp(ncol,pver)'
                   write(30,100) trim(line)
                   line(7:) = 'real' // trim(dec_suffix) // '  ::  exp_fac(ncol,pver)'
                else
                   line(7:) = 'real' // trim(dec_suffix) // '  ::  itemp(ncol,kbot)'
                   write(30,100) trim(line)
                   line(7:) = 'real' // trim(dec_suffix) // '  ::  exp_fac(ncol,kbot)'
                end if
             case( 'WRF' )
                line(7:) = 'real' // trim(dec_suffix) // '  ::  itemp(n)'
                write(30,100) trim(line)
                line(7:) = 'real' // trim(dec_suffix) // '  ::  exp_fac(n)'
             end select
          end if
          write(30,100) trim(line)
       end if temp_dep_rxts
       troe_rxts : if( troecnt /= 0 ) then
          if( machine == 'NEC' .or. machine == 'FUJITSU' ) then
             line(7:) = 'real' // trim(dec_suffix) // '  :: ko(' // trim(vec_dim) // ')'
             write(30,100) trim(line)
             line(7:) = 'real' // trim(dec_suffix) // '  :: kinf(' // trim(vec_dim) // ')'
          else
             select case( model )
             case( 'MOZART' )
                line(7:) = 'real' // trim(dec_suffix) // '  :: ko(plonl,plev)'
                write(30,100) trim(line)
                line(7:) = 'real' // trim(dec_suffix) // '  :: kinf(plonl,plev)'
             case( 'CAM' )
                if( subs == 1 ) then
                   line(7:) = 'real' // trim(dec_suffix) // '  :: ko(ncol,pver)'
                   write(30,100) trim(line)
                   line(7:) = 'real' // trim(dec_suffix) // '  :: kinf(ncol,pver)'
                else
                   line(7:) = 'real' // trim(dec_suffix) // '  :: ko(ncol,kbot)'
                   write(30,100) trim(line)
                   line(7:) = 'real' // trim(dec_suffix) // '  :: kinf(ncol,kbot)'
                   write(30,100) trim(line)
                   line(7:) = 'real' // trim(dec_suffix) // '  :: wrk(ncol,kbot)'
                end if
             case( 'WRF' )
                line(7:) = 'real' // trim(dec_suffix) // '  :: ko(n)'
                write(30,100) trim(line)
                line(7:) = 'real' // trim(dec_suffix) // '  :: kinf(n)'
             end select
          end if
          write(30,100) trim(line)
       end if troe_rxts
       line = ' '
       write(30,100) trim(line)

       !-----------------------------------------------------------------------
       !        ... first do all temperature independent rates
       !-----------------------------------------------------------------------
       const_rxts : if( rxpcnt > 0 .and. cnt /= rxpcnt ) then
          line = ' '
          do i = 1,rxpcnt
             if( sym_rates(2,i) == ' ' ) then
                if( subs == 2 ) then
                   if( .not. cph_flg(rxptab(i)) ) then
                      cycle
                   end if
                end if
                if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                   line(7:) = 'rate(:,'
                else
                   if( subs == 1 ) then
                      line(7:) = 'rate(:,:,'
                   else
                      line(7:) = 'rate(:,:kbot,'
                   end if
                end if
                write(num,'(i3)') rxptab(i)
                num = adjustl( num )
                l   = len_trim( sym_rates(1,i) )
                wrk = sym_rates(1,i)(:l)
                indp = scan( wrk(:l), '.' )
                inde = scan( wrk(:l), 'eE' )
                if( indp == 0 .and. inde == 0 ) then
                   l = l + 1
                   wrk(l:l) = '.'
                end if
                line(len_trim(line)+1:) = num(:len_trim(num)) // ') = ' // wrk(:l) // trim(num_suffix)
                write(30,100) trim(line)
             end if
          end do
       end if const_rxts

       !-----------------------------------------------------------------------
       !        ... now do temp dependent rxts
       !-----------------------------------------------------------------------
       do_tmp = cnt /= 0 .and. (subs == 1 .or. cph_cnt > 0)
       any_temp_dep :  if( do_tmp ) then
          if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
             line(7:) = 'itemp(:) = 1.' // trim(num_suffix)// ' / temp(:)'
          else if( model == 'MOZART' ) then
             line(7:) = 'itemp(:,:) = 1.' // trim(num_suffix)// ' / temp(:,:)'
          else if( model == 'CAM' ) then
             if( subs == 1 ) then
                line(7:) = 'itemp(:ncol,:) = 1.' // trim(num_suffix)// ' / temp(:ncol,:)'
             else
                line(7:) = 'itemp(:ncol,:kbot) = 1.' // trim(num_suffix)// ' / temp(:ncol,:kbot)'
             end if
          end if
          write(30,100) trim(line)
          if( vec_ftns ) then
             if( model == 'MOZART' ) then
                line(7:) = 'n = plonl*plev'
             else if( model == 'CAM' ) then
                if( march /= 'VECTOR' ) then
                   if( subs == 1 ) then
                      line(7:) = 'n = ncol*pver'
                   else
                      line(7:) = 'n = ncol*kbot'
                   end if
                else
                   line(7:) = 'n = ' // trim(vec_dim)
                end if
             end if
             write(30,100) trim(line)
          else if( model == 'CAM' ) then
             if( march /= 'VECTOR' ) then
                if( subs == 1 ) then
                   line(7:) = 'n = ncol*pver'
                else
                   line(7:) = 'n = ncol*kbot'
                end if
             else
                line(7:) = 'n = ' // trim(vec_dim)
             end if
             write(30,100) trim(line)
          end if
          line = ' '
          if( model == 'CAM' .and. machine == 'IBM' ) then
             ftn_lim = 2
          else
             ftn_lim = 1
          end if
          ftn_mode_loop : do ftn_mode = 1,ftn_lim
             if( model == 'CAM' .and. machine == 'IBM' ) then
                vftns = ftn_mode == 2
                if( ftn_mode == 1 ) then
                   !                 line = '#ifndef AIX'
                   line = '#if ( !defined AIX || defined NO_VEXP )'
                else if( ftn_mode == 2 ) then
                   line = '#else'
                end if
                write(30,100) trim(line)
                line = ' '
             else
                vftns = vec_ftns
             end if
             if( ftn_mode == 2 ) then
                t_dependent(:rxpcnt) = sym_rates(2,:rxpcnt) /= ' '
             end if
             rxt_loop : do i = 1,rxpcnt
                do_rxt = t_dependent(i) 
                if( subs == 2 ) then
                   if( do_rxt ) then
                      do_rxt = cph_flg(rxptab(i))
                   end if
                end if
                is_temp_dep : if( do_rxt ) then
                   match_cnt = 0
                   do m = i,rxpcnt
                      if( rxparm(2,i) == rxparm(2,m) ) then
                         if( subs == 1 ) then
                            match_cnt = match_cnt + 1
                            match_ind(match_cnt) = m
                            t_dependent(m) = .false.
                         else if( cph_flg(rxptab(m)) ) then
                            match_cnt = match_cnt + 1
                            match_ind(match_cnt) = m
                            t_dependent(m) = .false.
                         end if
                      end if
                   end do
                   multiple_matches : if( match_cnt > 1 ) then
                      if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                         line(7:) = 'exp_fac(:) = '
                      else if( .not. vftns ) then
                         line(7:) = 'exp_fac(:,:) = '
                      end if
                      l    = len_trim( sym_rates(2,i) )
                      wrk  = sym_rates(2,i)(:l)
                      indp = scan( wrk(:l), '.' )
                      inde = scan( wrk(:l), 'eE' )
                      if( indp == 0 .and. inde == 0 ) then
                         l = l + 1
                         wrk(l:l) = '.'
                      end if
                      pos = len_trim( line )
                      if( .not. vftns ) then
                         line(pos+1:) = ' exp( ' // wrk(:l) // trim(num_suffix)
                      else
                         line(7:) = 'call vexp( exp_fac, ' // wrk(:l) // trim(num_suffix)
                      end if
                      if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                         line(len_trim(line)+1:) = ' * itemp(:) )'
                      else if( .not. vftns ) then
                         line(len_trim(line)+1:) = ' * itemp(:,:) )'
                      else if( vftns ) then
                         line(len_trim(line)+1:) = ' * itemp, n )'
                      end if
                      write(30,100) trim(line)
                      do m = 1,match_cnt
                         if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                            line(7:) = 'rate(:,'
                         else
                            if( subs == 1 ) then
                               line(7:) = 'rate(:,:,'
                            else
                               line(7:) = 'rate(:,:kbot,'
                            end if
                         end if
                         m1 = match_ind(m)
                         write(num,'(i3)') rxptab(m1)
                         num = adjustl( num )
                         l = len_trim( sym_rates(1,m1) )
                         wrk = sym_rates(1,m1)(:l)
                         indp = scan( wrk(:l), '.' )
                         inde = scan( wrk(:l), 'eE' )
                         if( indp == 0 .and. inde == 0 ) then
                            l = l + 1
                            wrk(l:l) = '.'
                         end if
                         if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                            line(len_trim(line)+1:) = num(:len_trim(num)) // ') = ' // wrk(:l) // &
                                                      trim(num_suffix) // ' * exp_fac(:)'
                         else
                            line(len_trim(line)+1:) = num(:len_trim(num)) // ') = ' // wrk(:l) // &
                                                      trim(num_suffix) // ' * exp_fac(:,:)'
                         end if
                         write(30,100) trim(line)
                      end do
                   else multiple_matches
                      if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                         line(7:) = 'rate(:,'
                      else if( .not. vftns ) then
                         if( subs == 1 ) then
                            line(7:) = 'rate(:,:,'
                         else
                            line(7:) = 'rate(:,:kbot,'
                         end if
                      end if
                      write(num,'(i3)') rxptab(i)
                      num = adjustl( num )
                      l = len_trim( sym_rates(1,i) )
                      wrk  = sym_rates(1,i)(:l)
                      indp = scan( wrk(:l), '.' )
                      inde = scan( wrk(:l), 'eE' )
                      if( indp == 0 .and. inde == 0 ) then
                         l = l + 1
                         wrk(l:l) = '.'
                      end if
                      if( .not. vftns ) then
                         line(len_trim(line)+1:) = num(:len_trim(num)) // ') = ' // wrk(:l) // trim(num_suffix)
                      else
                         buff = wrk
                      end if
                      l = len_trim( sym_rates(2,i) )
                      wrk  = sym_rates(2,i)(:l)
                      indp = scan( wrk(:l), '.' )
                      inde = scan( wrk(:l), 'eE' )
                      if( indp == 0 .and. inde == 0 ) then
                         l = l + 1
                         wrk(l:l) = '.'
                      end if
                      if( vftns ) then
                         line(7:) = 'call vexp( exp_fac, ' // wrk(:l) // trim(num_suffix) // '*itemp, n )'
                         write(30,100) trim(line)
                         if( subs == 1 ) then
                            line(7:) = 'rate(:,:,' // trim(num) // ') = ' // trim(buff) // trim(num_suffix) // ' * exp_fac(:,:)'
                         else
                            line(7:) = 'rate(:,:kbot,' // trim(num) // ') = ' // trim(buff) // trim(num_suffix) // ' * exp_fac(:,:)'
                         end if
                      end if
                      if( .not. vftns ) then
                         pos = len_trim( line )
                         line(pos+1:) = ' * exp( ' // wrk(:l) // trim(num_suffix)
                         if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                            line(len_trim(line)+1:) = ' * itemp(:) )'
                         else
                            line(len_trim(line)+1:) = ' * itemp(:,:) )'
                         end if
                      end if
                      write(30,100) trim(line)
                   end if multiple_matches
                end if is_temp_dep
             end do rxt_loop
             if( model == 'CAM' .and. machine == 'IBM' .and. ftn_mode == 2 ) then
                line = '#endif'
                write(30,100) trim(line)
                line = ' '
             end if
          end do ftn_mode_loop
       end if any_temp_dep

       !-----------------------------------------------------------------------
       !        ... troe rates
       !-----------------------------------------------------------------------
       do_tmp = troecnt /= 0 .and. (subs == 1 .or. cph_cnt > 0)
       any_troe_rxts : if( do_tmp ) then
          line = ' '
          write(30,100) trim(line)
          if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
             line(7:) = 'itemp(:) = 300.' // trim(num_suffix) // ' * itemp(:)'
          else
             line(7:) = 'itemp(:,:) = 300.' // trim(num_suffix) // ' * itemp(:,:)'
          end if
          write(30,100) trim(line)
          troe_rxt_loop : do i = 1,troecnt
             if( subs == 2 ) then
                m1 = troetab(i)
                do_rxt = cph_flg(m1)
             else
                do_rxt = .true.
             end if
             line = ' '
             write(30,100) trim(line)
             do_troe_rate : if( do_rxt ) then
                if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                   line(7:) = 'ko(:)'
                else
                   line(7:) = 'ko(:,:)'
                end if
                l    = len_trim( troe_sym_rates(1,i) )
                wrk  = troe_sym_rates(1,i)(:l)
                indp = scan( wrk(:l), '.' )
                inde = scan( wrk(:l), 'eE' )
                if( indp == 0 .and. inde == 0 ) then
                   l = l + 1
                   wrk(l:l) = '.'
                end if
                line(len_trim(line)+1:) =  ' = ' // wrk(:l) // trim(num_suffix)
                if( troe_sym_rates(2,i) /= ' ' ) then
                   l = len_trim( troe_sym_rates(2,i) )
                   read(troe_sym_rates(2,i)(:l),*) rate
                   wrk  = troe_sym_rates(2,i)(:l)
                   indp = scan( wrk(:l), '.' )
                   inde = scan( wrk(:l), 'eE' )
                   if( indp == 0 .and. inde == 0 ) then
                      l = l + 1
                      wrk(l:l) = '.'
                   end if
                   if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                      if( rate /= 0. ) then
                         if ( rate < 0. ) then
                            line(len_trim(line)+1:) = ' * itemp(:)**(' // wrk(:l) // trim(num_suffix)//')'
                         else
                            line(len_trim(line)+1:) = ' * itemp(:)**' // wrk(:l) // trim(num_suffix)
                         endif
                      end if
                   else
                      if( rate /= 0. ) then
                         if ( rate < 0. ) then
                            line(len_trim(line)+1:) = ' * itemp(:,:)**(' // wrk(:l) // trim(num_suffix)//')'
                         else
                            line(len_trim(line)+1:) = ' * itemp(:,:)**' // wrk(:l) // trim(num_suffix)
                         endif
                      end if
                   end if
                end if
                write(30,100) trim(line)
                if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                   line(7:) = 'kinf(:)'
                else
                   line(7:) = 'kinf(:,:)'
                end if
                l = len_trim( troe_sym_rates(3,i) )
                wrk = troe_sym_rates(3,i)(:l)
                indp = scan( wrk(:l), '.' )
                inde = scan( wrk(:l), 'eE' )
                if( indp == 0 .and. inde == 0 ) then
                   l = l + 1
                   wrk(l:l) = '.'
                end if
                line(len_trim(line)+1:) =  ' = ' // wrk(:l) // trim(num_suffix)
                if( troe_sym_rates(4,i) /= ' ' ) then
                   l = len_trim( troe_sym_rates(4,i) )
                   read(troe_sym_rates(4,i)(:l),*) rate
                   wrk = troe_sym_rates(4,i)(:l)
                   indp = scan( wrk(:l), '.' )
                   inde = scan( wrk(:l), 'eE' )
                   if( indp == 0 .and. inde == 0 ) then
                      l = l + 1
                      wrk(l:l) = '.'
                   end if
                   if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                      if( rate /= 0. ) then
                         if( rate /= 1. ) then
                            if ( rate < 0. ) then
                               line(len_trim(line)+1:) = ' * itemp(:)**(' // wrk(:l) // trim(num_suffix)//')'
                            else
                               line(len_trim(line)+1:) = ' * itemp(:)**' // wrk(:l) // trim(num_suffix)
                            endif
                         else
                            line(len_trim(line)+1:) = ' * itemp(:)'
                         end if
                      end if
                   else
                      if( rate /= 0. ) then
                         if( rate /= 1. ) then
                            if ( rate < 0. ) then
                               line(len_trim(line)+1:) = ' * itemp(:,:)**(' // wrk(:l) // trim(num_suffix)//')'
                            else
                               line(len_trim(line)+1:) = ' * itemp(:,:)**' // wrk(:l) // trim(num_suffix)
                            endif
                         else
                            line(len_trim(line)+1:) = ' * itemp(:,:)'
                         end if
                      end if
                   end if
                end if
                write(30,100) trim(line)
                if( model == 'WRF' .or. machine == 'NEC' .or. machine == 'FUJITSU' ) then
                   line(7:) = 'call jpl( rate(1,'
                else
                   if( subs == 1 ) then
                      line(7:) = 'call jpl( rate(1,1,'
                   else
                      line(7:) = 'call jpl( wrk'
                   end if
                end if
                write(numa,'(i3)') troetab(i)
                if( subs == 1 ) then
                   num = numa
                else
                   num = ' '
                end if
                num = adjustl( num )
                l = len_trim( troe_sym_rates(5,i) )
                wrk = troe_sym_rates(5,i)(:l)
                indp = scan( wrk(:l), '.' )
                inde = scan( wrk(:l), 'eE' )
                if( indp == 0 .and. inde == 0 ) then
                   l = l + 1
                   wrk(l:l) = '.'
                end if
                select case( model )
                case( 'MOZART' )
                   line(len_trim(line)+1:) = num(:len_trim(num)) // '), m, ' // wrk(:l) // &
                                             trim(num_suffix) // ', ko, kinf, plnplv )'
                case( 'CAM' )
                   if( subs == 1 ) then
                      line(len_trim(line)+1:) = num(:len_trim(num)) // '), m, ' // wrk(:l) // &
                                                trim(num_suffix) // ', ko, kinf, n )'
                   else
                      line(len_trim(line)+1:) = num(:len_trim(num)) // ', m, ' // wrk(:l) // &
                                                trim(num_suffix) // ', ko, kinf, n )'
                      write(30,100) trim(line)
                      line = '      rate(:,:kbot,' // numa(:len_trim(numa)) // ') = wrk(:,:)'
                   end if
                case( 'WRF' )
                   line(len_trim(line)+1:) = num(:len_trim(num)) // '), m, ' // wrk(:l) // &
                                             trim(num_suffix) // ', ko, kinf, n )'
                end select
                write(30,100) trim(line)
             end if do_troe_rate
          end do troe_rxt_loop
       end if any_troe_rxts

       line = ' '
       write(30,100) trim(line)
       if( subs == 1 ) then
          line(7:) = 'end subroutine setrxt'
       else
          line(7:) = 'end subroutine setrxt_hrates'
       end if
       write(30,100) trim(line)
       line = ' '
       write(30,100) trim(line)
       if( subs == 2 ) then
          line(7:) = 'end module mo_setrxt'
          write(30,100) trim(line)
       end if
    end do subs_loop

    close(30)

100 format(a)

  end subroutine make_rate

end module set_rxt_rates
