

module mo_imp_sol
  use shr_kind_mod, only : r8 => shr_kind_r8
  use chem_mods,    only : clscnt4, gas_pcnst, clsmap
  use cam_logfile,  only : iulog
  implicit none
  private
  public :: imp_slv_inti, imp_sol
  save
  real(r8), parameter :: rel_err = 1.e-3_r8
  real(r8), parameter :: high_rel_err = 1.e-4_r8
  !-----------------------------------------------------------------------
  ! Newton-Raphson iteration limits
  !-----------------------------------------------------------------------
  integer, parameter :: itermax = 11
  integer, parameter :: cut_limit = 5
  real(r8) :: small
  real(r8) :: epsilon(clscnt4)
  logical :: factor(itermax)
  integer :: ox_ndx
  integer :: o1d_ndx = -1
  integer :: h2o_ndx = -1
  integer :: oh_ndx, ho2_ndx, ch3o2_ndx, po2_ndx, ch3co3_ndx
  integer :: c2h5o2_ndx, isopo2_ndx, macro2_ndx, mco3_ndx, c3h7o2_ndx
  integer :: ro2_ndx, xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx
  integer :: c2h4_ndx, c3h6_ndx, isop_ndx, mvk_ndx, c10h16_ndx
  integer :: ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, ox_p4_ndx, ox_p5_ndx
  integer :: ox_p6_ndx, ox_p7_ndx, ox_p8_ndx, ox_p9_ndx, ox_p10_ndx
  integer :: ox_p11_ndx
  integer :: ox_l1_ndx, ox_l2_ndx, ox_l3_ndx, ox_l4_ndx, ox_l5_ndx
  integer :: ox_l6_ndx, ox_l7_ndx, ox_l8_ndx, ox_l9_ndx, usr4_ndx
  integer :: usr16_ndx, usr17_ndx, r63_ndx,c2o3_ndx,ole_ndx
  integer :: tolo2_ndx, terpo2_ndx, alko2_ndx, eneo2_ndx, eo2_ndx, meko2_ndx
  integer :: ox_p17_ndx,ox_p12_ndx,ox_p13_ndx,ox_p14_ndx,ox_p15_ndx,ox_p16_ndx
  integer :: lt_cnt
  logical :: full_ozone_chem = .false.
  logical :: reduced_ozone_chem = .false.
  ! for xnox ozone chemistry diagnostics
  integer :: o3a_ndx, xno2_ndx, no2xno3_ndx, xno2no3_ndx, xno3_ndx, o1da_ndx, xno_ndx
  integer :: usr4a_ndx, usr16a_ndx, usr16b_ndx, usr17b_ndx
contains
  subroutine imp_slv_inti
    !-----------------------------------------------------------------------
    ! ... Initialize the implict solver
    !-----------------------------------------------------------------------
    use mo_chem_utls, only : get_spc_ndx, get_rxt_ndx
    use cam_abortutils, only : endrun
    use cam_history, only : addfld, add_default
    use ppgrid, only : pver
    use mo_tracname, only : solsym
    implicit none
    !-----------------------------------------------------------------------
    ! ... Local variables
    !-----------------------------------------------------------------------
    integer :: m
    real(r8) :: eps(gas_pcnst)
    integer :: wrk(27)
    integer :: i,j
    ! small = 1.e6_r8 * tiny( small )
    small = 1.e-40_r8
    factor(:) = .true.
    eps(:) = rel_err
    ox_ndx = get_spc_ndx( 'OX' )
    h2o_ndx = get_spc_ndx( 'H2O' )
    if( ox_ndx < 1 ) then
       ox_ndx = get_spc_ndx( 'O3' )
       o1d_ndx = get_spc_ndx( 'O1D' )
       o1da_ndx = get_spc_ndx( 'O1DA' )
    end if
    if( ox_ndx > 0 ) then
       eps(ox_ndx) = high_rel_err
    end if
    m = get_spc_ndx( 'NO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'N2O5' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'OH' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'HO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    o3a_ndx = get_spc_ndx( 'O3A' )
    if( o3a_ndx > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XHO2NO2' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'XNO2NO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    m = get_spc_ndx( 'NO2XNO3' )
    if( m > 0 ) then
       eps(m) = high_rel_err
    end if
    !do m = 1,max(1,clscnt4)
    do m = 1,clscnt4
       epsilon(m) = eps(clsmap(m,4))
    end do
    has_o3_chem: if( ox_ndx > 0 ) then
       ox_p1_ndx = get_rxt_ndx( 'ox_p1' )
       ox_p2_ndx = get_rxt_ndx( 'ox_p2' )
       ox_p3_ndx = get_rxt_ndx( 'ox_p3' )
       ox_p4_ndx = get_rxt_ndx( 'ox_p4' )
       ox_p5_ndx = get_rxt_ndx( 'ox_p5' )
       ox_p6_ndx = get_rxt_ndx( 'ox_p6' )
       ox_p7_ndx = get_rxt_ndx( 'ox_p7' )
       ox_p8_ndx = get_rxt_ndx( 'ox_p8' )
       ox_p9_ndx = get_rxt_ndx( 'ox_p9' )
       ox_p10_ndx = get_rxt_ndx( 'ox_p10' )
       ox_p11_ndx = get_rxt_ndx( 'ox_p11' )
       ox_p12_ndx = get_rxt_ndx( 'ox_p12' )
       ox_p13_ndx = get_rxt_ndx( 'ox_p13' )
       ox_p14_ndx = get_rxt_ndx( 'ox_p14' )
       ox_p15_ndx = get_rxt_ndx( 'ox_p15' )
       ox_p16_ndx = get_rxt_ndx( 'ox_p16' )
       ox_p17_ndx = get_rxt_ndx( 'ox_p17' )
       wrk(1:17) = (/ ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, ox_p4_ndx, ox_p5_ndx, &
            ox_p6_ndx, ox_p7_ndx, ox_p8_ndx, ox_p9_ndx, ox_p10_ndx, ox_p11_ndx, &
            ox_p12_ndx, ox_p13_ndx, ox_p14_ndx, ox_p15_ndx, ox_p16_ndx, ox_p17_ndx /)
       if( all( wrk(1:17) > 0 ) ) then
          full_ozone_chem = .true.
       end if
       if ( .not. full_ozone_chem ) then
          r63_ndx = get_rxt_ndx( 'r63' )
          wrk(1:4) = (/ ox_p1_ndx, ox_p2_ndx, ox_p3_ndx, r63_ndx/)
          if( all( wrk(1:4) > 0 ) ) then
             reduced_ozone_chem = .true.
          end if
       endif
       if( full_ozone_chem .or. reduced_ozone_chem ) then
          ox_l1_ndx = get_rxt_ndx( 'ox_l1' )
          ox_l2_ndx = get_rxt_ndx( 'ox_l2' )
          ox_l3_ndx = get_rxt_ndx( 'ox_l3' )
          ox_l4_ndx = get_rxt_ndx( 'ox_l4' )
          ox_l5_ndx = get_rxt_ndx( 'ox_l5' )
          ox_l6_ndx = get_rxt_ndx( 'ox_l6' )
          ox_l7_ndx = get_rxt_ndx( 'ox_l7' )
          ox_l8_ndx = get_rxt_ndx( 'ox_l8' )
          ox_l9_ndx = get_rxt_ndx( 'ox_l9' )
          if( ox_l9_ndx < 1 ) then
             ox_l9_ndx = get_rxt_ndx( 'soa1' )
          end if
          if( ox_l9_ndx < 1 ) then
             ox_l9_ndx = get_rxt_ndx( 'C10H16_O3' )
          end if
          usr4_ndx = get_rxt_ndx( 'usr4' )
          if (usr4_ndx < 1) then
             usr4_ndx = get_rxt_ndx( 'tag_NO2_OH' )
          endif
          usr16_ndx = get_rxt_ndx( 'usr16' )
          if (usr16_ndx < 1) then
            usr16_ndx = get_rxt_ndx( 'usr_N2O5_aer' )
          endif
          usr17_ndx = get_rxt_ndx( 'usr17' )
          if (usr17_ndx < 1) then
            usr17_ndx = get_rxt_ndx( 'usr_NO3_aer' )
          endif
          usr4a_ndx = get_rxt_ndx( 'usr4a' )
          if (usr4a_ndx < 1) then
            usr4a_ndx = get_rxt_ndx( 'tag_XNO2_OH' )
          endif
          usr16b_ndx = get_rxt_ndx( 'usr16b' )
          if (usr16b_ndx < 1) then
            usr16b_ndx = get_rxt_ndx( 'usr_NO2XNO3_aer' )
          endif
          usr16a_ndx = get_rxt_ndx( 'usr16a' )
          if (usr16a_ndx < 1) then
            usr16a_ndx = get_rxt_ndx( 'usr_XNO2NO3_aer' )
          endif
          usr17b_ndx = get_rxt_ndx( 'usr17b' )
          if (usr17b_ndx < 1) then
            usr17b_ndx = get_rxt_ndx( 'usr_NO2_aer' )
          endif
          if ( full_ozone_chem ) then
             wrk(1:12) = (/ ox_l1_ndx, ox_l2_ndx, ox_l3_ndx, ox_l4_ndx, ox_l5_ndx, &
                  ox_l6_ndx, ox_l7_ndx, ox_l8_ndx, ox_l9_ndx, usr4_ndx, &
                  usr16_ndx, usr17_ndx /)
             if( any( wrk(1:12) < 1 ) ) then
                full_ozone_chem = .false.
             endif
          endif
          if ( reduced_ozone_chem ) then
             wrk(1:9) = (/ ox_l1_ndx, ox_l2_ndx, ox_l3_ndx, ox_l5_ndx, &
                  ox_l6_ndx, ox_l7_ndx, usr4_ndx, &
                  usr16_ndx, usr17_ndx /)
             if( any( wrk(1:9) < 1 ) ) then
                reduced_ozone_chem = .false.
             end if
          endif
       end if
       if( full_ozone_chem .or. reduced_ozone_chem ) then
          oh_ndx = get_spc_ndx( 'OH' )
          ho2_ndx = get_spc_ndx( 'HO2' )
          ch3o2_ndx = get_spc_ndx( 'CH3O2' )
          po2_ndx = get_spc_ndx( 'PO2' )
          ch3co3_ndx = get_spc_ndx( 'CH3CO3' )
          c2h5o2_ndx = get_spc_ndx( 'C2H5O2' )
          macro2_ndx = get_spc_ndx( 'MACRO2' )
          mco3_ndx = get_spc_ndx( 'MCO3' )
          c3h7o2_ndx = get_spc_ndx( 'C3H7O2' )
          ro2_ndx = get_spc_ndx( 'RO2' )
          xo2_ndx = get_spc_ndx( 'XO2' )
          no_ndx = get_spc_ndx( 'NO' )
          xno_ndx = get_spc_ndx( 'XNO' )
          no2_ndx = get_spc_ndx( 'NO2' )
          xno2_ndx = get_spc_ndx( 'XNO2' )
          no3_ndx = get_spc_ndx( 'NO3' )
          xno3_ndx = get_spc_ndx( 'XNO3' )
          n2o5_ndx = get_spc_ndx( 'N2O5' )
          xno2no3_ndx = get_spc_ndx( 'XNO2NO3' )
          no2xno3_ndx = get_spc_ndx( 'NO2XNO3' )
          c2h4_ndx = get_spc_ndx( 'C2H4' )
          c3h6_ndx = get_spc_ndx( 'C3H6' )
          isop_ndx = get_spc_ndx( 'ISOP' )
          isopo2_ndx = get_spc_ndx( 'ISOPO2' )
          mvk_ndx = get_spc_ndx( 'MVK' )
          c10h16_ndx = get_spc_ndx( 'C10H16' )
          tolo2_ndx = get_spc_ndx('TOLO2')
          terpo2_ndx =get_spc_ndx('TERPO2')
          alko2_ndx = get_spc_ndx('ALKO2')
          eneo2_ndx = get_spc_ndx('ENEO2')
          eo2_ndx = get_spc_ndx('EO2')
          meko2_ndx = get_spc_ndx('MEKO2')
          if ( full_ozone_chem ) then
             wrk(1:27) = (/ oh_ndx, ho2_ndx, ch3o2_ndx, po2_ndx, ch3co3_ndx, &
                  c2h5o2_ndx, macro2_ndx, mco3_ndx, c3h7o2_ndx, ro2_ndx, &
                  xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx, &
                  c2h4_ndx, c3h6_ndx, isop_ndx, isopo2_ndx, mvk_ndx, c10h16_ndx, &
                  tolo2_ndx, terpo2_ndx, alko2_ndx, eneo2_ndx, eo2_ndx, meko2_ndx /)
             if( any( wrk(1:27) < 1 ) ) then
                full_ozone_chem = .false.
             end if
          endif
          if ( reduced_ozone_chem ) then
             c2o3_ndx = get_spc_ndx( 'C2O3' )
             ole_ndx = get_spc_ndx( 'OLE' )
             wrk(1:12) = (/ oh_ndx, ho2_ndx, ch3o2_ndx, &
                  ole_ndx,c2o3_ndx, &
                  xo2_ndx, no_ndx, no2_ndx, no3_ndx, n2o5_ndx, &
                  c2h4_ndx, isop_ndx /)
             if( any( wrk(1:12) < 1 ) ) then
                reduced_ozone_chem = .false.
             end if
          endif
       end if
    else
       reduced_ozone_chem = .false.
       full_ozone_chem = .false.
    end if has_o3_chem
    if ( reduced_ozone_chem .and. full_ozone_chem ) then
       write(iulog,*) 'can not have both full_ozone_chem and reduced_ozone_chem'
       call endrun
    endif
    do i = 1,clscnt4
       j = clsmap(i,4)
       call addfld( trim(solsym(j))//'_CHMP', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate' )
       call addfld( trim(solsym(j))//'_CHML', (/ 'lev' /), 'I', '/cm3/s', 'chemical loss rate' )
    enddo
  end subroutine imp_slv_inti
  subroutine imp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, &
       xhnm, ncol, lchnk, ltrop )
    !-----------------------------------------------------------------------
    ! ... imp_sol advances the volumetric mixing ratio
    ! forward one time step via the fully implicit euler scheme.
    ! this source is meant for small l1 cache machines such as
    ! the intel pentium and itanium cpus
    !-----------------------------------------------------------------------
    use chem_mods, only : rxntot, extcnt, nzcnt, permute, cls_rxt_cnt
    use mo_tracname, only : solsym
    use ppgrid, only : pver
    use mo_lin_matrix, only : linmat
    use mo_nln_matrix, only : nlnmat
    use mo_lu_factor, only : lu_fac
    use mo_lu_solve, only : lu_slv
    use mo_prod_loss, only : imp_prod_loss
    use mo_indprd, only : indprd
    use time_manager, only : get_nstep
    use cam_history, only : outfld
    implicit none
    !-----------------------------------------------------------------------
    ! ... dummy args
    !-----------------------------------------------------------------------
    integer, intent(in) :: ncol ! columns in chunck
    integer, intent(in) :: lchnk ! chunk id
    real(r8), intent(in) :: delt ! time step (s)
    real(r8), intent(in) :: reaction_rates(ncol,pver,max(1,rxntot)), & ! rxt rates (1/cm^3/s)
         extfrc(ncol,pver,max(1,extcnt)), & ! external in-situ forcing (1/cm^3/s)
         het_rates(ncol,pver,max(1,gas_pcnst)) ! washout rates (1/s)
    real(r8), intent(inout) :: base_sol(ncol,pver,gas_pcnst) ! species mixing ratios (vmr)
    real(r8), intent(in) :: xhnm(ncol,pver)
    integer, intent(in) :: ltrop(ncol) ! chemistry troposphere boundary (index)
    !-----------------------------------------------------------------------
    ! ... local variables
    !-----------------------------------------------------------------------
    integer :: nr_iter, &
         lev, &
         i, &
         j, &
         k, l, &
         m
    integer :: fail_cnt, cut_cnt, stp_con_cnt
    integer :: nstep
    real(r8) :: interval_done, dt, dti, wrk
    real(r8) :: max_delta(max(1,clscnt4))
    real(r8) :: sys_jac(max(1,nzcnt))
    real(r8) :: lin_jac(max(1,nzcnt))
    real(r8), dimension(max(1,clscnt4)) :: &
         solution, &
         forcing, &
         iter_invariant, &
         prod, &
         loss
    real(r8) :: lrxt(max(1,rxntot))
    real(r8) :: lsol(max(1,gas_pcnst))
    real(r8) :: lhet(max(1,gas_pcnst))
    real(r8), dimension(ncol,pver,max(1,clscnt4)) :: &
         ind_prd
    logical :: convergence
    logical :: frc_mask, iter_conv
    logical :: converged(max(1,clscnt4))
    real(r8), dimension(ncol,pver,max(1,clscnt4)) :: prod_out, loss_out
    prod_out(:,:,:) = 0._r8
    loss_out(:,:,:) = 0._r8
    solution(:) = 0._r8
    !-----------------------------------------------------------------------
    ! ... class independent forcing
    !-----------------------------------------------------------------------
    if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
       call indprd( 4, ind_prd, clscnt4, base_sol, extfrc, &
            reaction_rates, ncol )
    else
       do m = 1,max(1,clscnt4)
          ind_prd(:,:,m) = 0._r8
       end do
    end if
    level_loop : do lev = 1,pver
       column_loop : do i = 1,ncol
          IF (lev <= ltrop(i)) CYCLE column_loop
          !-----------------------------------------------------------------------
          ! ... transfer from base to local work arrays
          !-----------------------------------------------------------------------
          do m = 1,rxntot
             lrxt(m) = reaction_rates(i,lev,m)
          end do
          if( gas_pcnst > 0 ) then
             do m = 1,gas_pcnst
                lhet(m) = het_rates(i,lev,m)
             end do
          end if
          !-----------------------------------------------------------------------
          ! ... time step loop
          !-----------------------------------------------------------------------
          dt = delt
          cut_cnt = 0
          fail_cnt = 0
          stp_con_cnt = 0
          interval_done = 0._r8
          time_step_loop : do
             dti = 1._r8 / dt
             !-----------------------------------------------------------------------
             ! ... transfer from base to local work arrays
             !-----------------------------------------------------------------------
             do m = 1,gas_pcnst
                lsol(m) = base_sol(i,lev,m)
             end do
             !-----------------------------------------------------------------------
             ! ... transfer from base to class array
             !-----------------------------------------------------------------------
             do k = 1,clscnt4
                j = clsmap(k,4)
                m = permute(k,4)
                solution(m) = lsol(j)
             end do
             !-----------------------------------------------------------------------
             ! ... set the iteration invariant part of the function f(y)
             !-----------------------------------------------------------------------
             if( cls_rxt_cnt(1,4) > 0 .or. extcnt > 0 ) then
                do m = 1,clscnt4
                   iter_invariant(m) = dti * solution(m) + ind_prd(i,lev,m)
                end do
             else
                do m = 1,clscnt4
                   iter_invariant(m) = dti * solution(m)
                end do
             end if
             !-----------------------------------------------------------------------
             ! ... the linear component
             !-----------------------------------------------------------------------
             !if( cls_rxt_cnt(2,4) > 0 ) then
                call linmat( lin_jac, lsol, lrxt, lhet )
             !end if
             !=======================================================================
             ! the newton-raphson iteration for f(y) = 0
             !=======================================================================
             iter_loop : do nr_iter = 1,itermax
                !-----------------------------------------------------------------------
                ! ... the non-linear component
                !-----------------------------------------------------------------------
                if( factor(nr_iter) ) then
                   call nlnmat( sys_jac, lsol, lrxt, lin_jac, dti )
                   !-----------------------------------------------------------------------
                   ! ... factor the "system" matrix
                   !-----------------------------------------------------------------------
                   call lu_fac( sys_jac )
                end if
                !-----------------------------------------------------------------------
                ! ... form f(y)
                !-----------------------------------------------------------------------
                call imp_prod_loss( prod, loss, lsol, lrxt, lhet )
                do m = 1,clscnt4
                   forcing(m) = solution(m)*dti - (iter_invariant(m) + prod(m) - loss(m))
                end do
                !-----------------------------------------------------------------------
                ! ... solve for the mixing ratio at t(n+1)
                !-----------------------------------------------------------------------
                call lu_slv( sys_jac, forcing )
                do m = 1,clscnt4
                   solution(m) = solution(m) + forcing(m)
                end do
                !-----------------------------------------------------------------------
                ! ... convergence measures
                !-----------------------------------------------------------------------
                if( nr_iter > 1 ) then
                   do k = 1,clscnt4
                      m = permute(k,4)
                      if( abs(solution(m)) > 1.e-20_r8 ) then
                         max_delta(k) = abs( forcing(m)/solution(m) )
                      else
                         max_delta(k) = 0._r8
                      end if
                   end do
                end if
                !-----------------------------------------------------------------------
                ! ... limit iterate
                !-----------------------------------------------------------------------
                where( solution(:) < 0._r8 )
                   solution(:) = 0._r8
                endwhere
                !-----------------------------------------------------------------------
                ! ... transfer latest solution back to work array
                !-----------------------------------------------------------------------
                do k = 1,clscnt4
                   j = clsmap(k,4)
                   m = permute(k,4)
                   lsol(j) = solution(m)
                end do
                !-----------------------------------------------------------------------
                ! ... check for convergence
                !-----------------------------------------------------------------------
                converged(:) = .true.
                if( nr_iter > 1 ) then
                   do k = 1,clscnt4
                      m = permute(k,4)
                      frc_mask = abs( forcing(m) ) > small
                      if( frc_mask ) then
                         converged(k) = abs(forcing(m)) <= epsilon(k)*abs(solution(m))
                      else
                         converged(k) = .true.
                      end if
                   end do
                   convergence = all( converged(:) )
                   if( convergence ) then
                      exit
                   end if
                end if
             end do iter_loop
             !-----------------------------------------------------------------------
             ! ... check for newton-raphson convergence
             !-----------------------------------------------------------------------
             if( .not. convergence ) then
                !-----------------------------------------------------------------------
                ! ... non-convergence
                !-----------------------------------------------------------------------
                fail_cnt = fail_cnt + 1
                nstep = get_nstep()
                write(iulog,'('' imp_sol: Time step '',1p,e21.13,'' failed to converge @ (lchnk,lev,col,nstep) = '',4i6)') &
                     dt,lchnk,lev,i,nstep
                stp_con_cnt = 0
                if( cut_cnt < cut_limit ) then
                   cut_cnt = cut_cnt + 1
                   if( cut_cnt < cut_limit ) then
                      dt = .5_r8 * dt
                   else
                      dt = .1_r8 * dt
                   end if
                   cycle time_step_loop
                else
                   write(iulog,'('' imp_sol: Failed to converge @ (lchnk,lev,col,nstep,dt,time) = '',4i6,1p,2e21.13)') &
                        lchnk,lev,i,nstep,dt,interval_done+dt
                   do m = 1,clscnt4
                      if( .not. converged(m) ) then
                         write(iulog,'(1x,a8,1x,1pe10.3)') solsym(clsmap(m,4)), max_delta(m)
                      end if
                   end do
                end if
             end if
             !-----------------------------------------------------------------------
             ! ... check for interval done
             !-----------------------------------------------------------------------
             interval_done = interval_done + dt
             if( abs( delt - interval_done ) <= .0001_r8 ) then
                if( fail_cnt > 0 ) then
                   write(iulog,*) 'imp_sol : @ (lchnk,lev,col) = ',lchnk,lev,i,' failed ',fail_cnt,' times'
                end if
                exit time_step_loop
             else
                !-----------------------------------------------------------------------
                ! ... transfer latest solution back to base array
                !-----------------------------------------------------------------------
                if( convergence ) then
                   stp_con_cnt = stp_con_cnt + 1
                end if
                do m = 1,gas_pcnst
                   base_sol(i,lev,m) = lsol(m)
                end do
                if( stp_con_cnt >= 2 ) then
                   dt = 2._r8*dt
                   stp_con_cnt = 0
                end if
                dt = min( dt,delt-interval_done )
                ! write(iulog,'('' imp_sol: New time step '',1p,e21.13)') dt
             end if
          end do time_step_loop
          !-----------------------------------------------------------------------
          ! ... Transfer latest solution back to base array
          !-----------------------------------------------------------------------
          cls_loop: do k = 1,clscnt4
             j = clsmap(k,4)
             m = permute(k,4)
             base_sol(i,lev,j) = solution(m)
          end do cls_loop
          !-----------------------------------------------------------------------
          ! ... Prod/Loss history buffers...
          !-----------------------------------------------------------------------
          cls_loop2: do k = 1,clscnt4
             j = clsmap(k,4)
             m = permute(k,4)
             has_o3_chem: if( ( full_ozone_chem .or. reduced_ozone_chem ) .and. (j == ox_ndx .or. j == o3a_ndx )) then
                if( o1d_ndx < 1 ) then
                   loss_out(i,lev,k) = reaction_rates(i,lev,ox_l1_ndx)
                else
                   if (j == ox_ndx) &
                      loss_out(i,lev,k) = reaction_rates(i,lev,ox_l1_ndx) * base_sol(i,lev,o1d_ndx)/base_sol(i,lev,ox_ndx)
                   if (j == o3a_ndx) &
                      loss_out(i,lev,k) = reaction_rates(i,lev,ox_l1_ndx) * base_sol(i,lev,o1da_ndx)/base_sol(i,lev,o3a_ndx)
                   if ( h2o_ndx > 0 ) &
                      loss_out(i,lev,k) = loss_out(i,lev,k) * base_sol(i,lev,h2o_ndx)
                end if
                if ( full_ozone_chem ) then
                   prod_out(i,lev,k) = reaction_rates(i,lev,ox_p1_ndx) * base_sol(i,lev,ho2_ndx) &
                        + reaction_rates(i,lev,ox_p2_ndx) * base_sol(i,lev,ch3o2_ndx) &
                        + reaction_rates(i,lev,ox_p3_ndx) * base_sol(i,lev,po2_ndx) &
                        + reaction_rates(i,lev,ox_p4_ndx) * base_sol(i,lev,ch3co3_ndx) &
                        + reaction_rates(i,lev,ox_p5_ndx) * base_sol(i,lev,c2h5o2_ndx) &
                        + .92_r8* reaction_rates(i,lev,ox_p6_ndx) * base_sol(i,lev,isopo2_ndx) &
                        + reaction_rates(i,lev,ox_p7_ndx) * base_sol(i,lev,macro2_ndx) &
                        + reaction_rates(i,lev,ox_p8_ndx) * base_sol(i,lev,mco3_ndx) &
                        + reaction_rates(i,lev,ox_p9_ndx) * base_sol(i,lev,c3h7o2_ndx) &
                        + reaction_rates(i,lev,ox_p10_ndx)* base_sol(i,lev,ro2_ndx) &
                        + reaction_rates(i,lev,ox_p11_ndx)* base_sol(i,lev,xo2_ndx) &
                        + .9_r8*reaction_rates(i,lev,ox_p12_ndx)*base_sol(i,lev,tolo2_ndx) &
                        + reaction_rates(i,lev,ox_p13_ndx)*base_sol(i,lev,terpo2_ndx)&
                        + .9_r8*reaction_rates(i,lev,ox_p14_ndx)*base_sol(i,lev,alko2_ndx) &
                        + reaction_rates(i,lev,ox_p15_ndx)*base_sol(i,lev,eneo2_ndx) &
                        + reaction_rates(i,lev,ox_p16_ndx)*base_sol(i,lev,eo2_ndx) &
                        + reaction_rates(i,lev,ox_p17_ndx)*base_sol(i,lev,meko2_ndx)
                   loss_out(i,lev,k) = loss_out(i,lev,k) &
                        + reaction_rates(i,lev,ox_l2_ndx) * base_sol(i,lev,oh_ndx) &
                        + reaction_rates(i,lev,ox_l3_ndx) * base_sol(i,lev,ho2_ndx) &
                        + reaction_rates(i,lev,ox_l6_ndx) * base_sol(i,lev,c2h4_ndx) &
                        + reaction_rates(i,lev,ox_l4_ndx) * base_sol(i,lev,c3h6_ndx) &
                        + .9_r8* reaction_rates(i,lev,ox_l5_ndx) * base_sol(i,lev,isop_ndx) &
                        + .8_r8*( reaction_rates(i,lev,ox_l7_ndx) * base_sol(i,lev,mvk_ndx) &
                        + reaction_rates(i,lev,ox_l8_ndx) * base_sol(i,lev,macro2_ndx)) &
                        + .235_r8*reaction_rates(i,lev,ox_l9_ndx) * base_sol(i,lev,c10h16_ndx)
                else if ( reduced_ozone_chem ) then
                   prod_out(i,lev,k) = reaction_rates(i,lev,ox_p1_ndx ) * base_sol(i,lev,ho2_ndx) &
                        + reaction_rates(i,lev,ox_p2_ndx ) * base_sol(i,lev,ch3o2_ndx) &
                        + reaction_rates(i,lev,ox_p3_ndx ) * base_sol(i,lev,c2o3_ndx) &
                        + reaction_rates(i,lev,r63_ndx ) * base_sol(i,lev,xo2_ndx)
                   loss_out(i,lev,k) = loss_out(i,lev,k) &
                        + reaction_rates(i,lev,ox_l2_ndx) * base_sol(i,lev,oh_ndx) &
                        + reaction_rates(i,lev,ox_l3_ndx) * base_sol(i,lev,ho2_ndx) &
                        + .9_r8* reaction_rates(i,lev,ox_l5_ndx) * base_sol(i,lev,isop_ndx) &
                        + reaction_rates(i,lev,ox_l6_ndx) * base_sol(i,lev,c2h4_ndx) &
                        + reaction_rates(i,lev,ox_l7_ndx) * base_sol(i,lev,ole_ndx)
                endif
                if (j == ox_ndx) then
                   loss_out(i,lev,k) = loss_out(i,lev,k) &
                        + ( reaction_rates(i,lev,usr4_ndx) * base_sol(i,lev,no2_ndx) * base_sol(i,lev,oh_ndx) &
                        + 3._r8 * reaction_rates(i,lev,usr16_ndx) * base_sol(i,lev,n2o5_ndx) &
                        + 2._r8 * reaction_rates(i,lev,usr17_ndx) * base_sol(i,lev,no3_ndx) ) &
                        / max( base_sol(i,lev,ox_ndx),1.e-20_r8 )
                   loss_out(i,lev,k) = loss_out(i,lev,k) * base_sol(i,lev,ox_ndx)
                   prod_out(i,lev,k) = prod_out(i,lev,k) * base_sol(i,lev,no_ndx)
                else if (j == o3a_ndx) then
                   loss_out(i,lev,k) = loss_out(i,lev,k) &
                        + ( reaction_rates(i,lev,usr4a_ndx) * base_sol(i,lev,xno2_ndx) * base_sol(i,lev,oh_ndx) &
                        + 1._r8 * reaction_rates(i,lev,usr16a_ndx) * base_sol(i,lev,xno2no3_ndx) &
                        + 2._r8 * reaction_rates(i,lev,usr16b_ndx) * base_sol(i,lev,no2xno3_ndx) &
                        + 2._r8 * reaction_rates(i,lev,usr17b_ndx) * base_sol(i,lev,xno3_ndx) ) &
                        / max( base_sol(i,lev,o3a_ndx),1.e-20_r8 )
                   loss_out(i,lev,k) = loss_out(i,lev,k) * base_sol(i,lev,o3a_ndx)
                   prod_out(i,lev,k) = prod_out(i,lev,k) * base_sol(i,lev,xno_ndx)
                endif
             else
                prod_out(i,lev,k) = prod(m) + ind_prd(i,lev,m)
                loss_out(i,lev,k) = loss(m)
             endif has_o3_chem
          end do cls_loop2
       end do column_loop
    end do level_loop
    do i = 1,clscnt4
       j = clsmap(i,4)
       prod_out(:,:,i) = prod_out(:,:,i)*xhnm
       loss_out(:,:,i) = loss_out(:,:,i)*xhnm
       call outfld( trim(solsym(j))//'_CHMP', prod_out(:,:,i), ncol, lchnk )
       call outfld( trim(solsym(j))//'_CHML', loss_out(:,:,i), ncol, lchnk )
    enddo
  end subroutine imp_sol
end module mo_imp_sol
