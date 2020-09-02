
module mo_exp_sol

  private
  public :: exp_sol
  public :: exp_sol_inti

  save

  integer :: uci_CH2O_ndx,uci_CH3O2_ndx,uci_CH3OOH_ndx,uci_PAN_ndx,uci_CO_ndx, &
             uci_C2H6_ndx,uci_C3H8_ndx,uci_C2H4_ndx,uci_ROHO2_ndx, &
             uci_CH3COCH3_ndx,uci_C2H5O2_ndx,uci_C2H5OOH_ndx,uci_CH3CHO_ndx, &
             uci_CH3CO3_ndx,uci_ISOP_ndx,uci_ISOPO2_ndx,uci_MVKMACR_ndx, &
             uci_MVKO2_ndx

contains

  subroutine exp_sol_inti

    use mo_tracname, only : solsym
    use chem_mods,   only : clscnt1, clsmap
    use ppgrid,      only : pver
    use cam_history, only : addfld, add_default
    use mo_chem_utls, only : get_rxt_ndx
    use cam_logfile,  only : iulog
    use spmd_utils,   only : masterproc

    implicit none

    integer :: i,j

    do i = 1,clscnt1

       j = clsmap(i,1)
       call addfld( trim(solsym(j))//'_CHMP', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate' )
       call addfld( trim(solsym(j))//'_CHML', (/ 'lev' /), 'I', '/cm3/s', 'chemical loss rate' )

    enddo

    !-----------------------------------------------------------------------      
    !    	... Reaction rate indices
    !-----------------------------------------------------------------------      
    uci_CH2O_ndx     = get_rxt_ndx( 'uci_CH2O' )
    uci_CH3O2_ndx    = get_rxt_ndx( 'uci_CH3O2' )
    uci_CH3OOH_ndx   = get_rxt_ndx( 'uci_CH3OOH' )
    uci_PAN_ndx      = get_rxt_ndx( 'uci_PAN' )
    uci_CO_ndx       = get_rxt_ndx( 'uci_CO' )
    uci_C2H6_ndx     = get_rxt_ndx( 'uci_C2H6' )
    uci_C3H8_ndx     = get_rxt_ndx( 'uci_C3H8' )
    uci_C2H4_ndx     = get_rxt_ndx( 'uci_C2H4' )
    uci_ROHO2_ndx    = get_rxt_ndx( 'uci_ROHO2' )
    uci_CH3COCH3_ndx = get_rxt_ndx( 'uci_CH3COCH3' )
    uci_C2H5O2_ndx   = get_rxt_ndx( 'uci_C2H5O2' )
    uci_C2H5OOH_ndx  = get_rxt_ndx( 'uci_C2H5OOH' )
    uci_CH3CHO_ndx   = get_rxt_ndx( 'uci_CH3CHO' )
    uci_CH3CO3_ndx   = get_rxt_ndx( 'uci_CH3CO3' )
    uci_ISOP_ndx     = get_rxt_ndx( 'uci_ISOP' )
    uci_ISOPO2_ndx   = get_rxt_ndx( 'uci_ISOPO2' )
    uci_MVKMACR_ndx  = get_rxt_ndx( 'uci_MVKMACR' )
    uci_MVKO2_ndx    = get_rxt_ndx( 'uci_MVKO2' )

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'exp_sol_inti: diagnostics uci rates indices'
       write(iulog,'(10i5)') uci_CH2O_ndx,uci_CH3O2_ndx,uci_CH3OOH_ndx,uci_PAN_ndx,uci_CO_ndx, &
             uci_C2H6_ndx,uci_C3H8_ndx,uci_C2H4_ndx,uci_ROHO2_ndx, &
             uci_CH3COCH3_ndx,uci_C2H5O2_ndx,uci_C2H5OOH_ndx,uci_CH3CHO_ndx, &
             uci_CH3CO3_ndx,uci_ISOP_ndx,uci_ISOPO2_ndx,uci_MVKMACR_ndx, &
             uci_MVKO2_ndx
    end if

  end subroutine exp_sol_inti


  subroutine exp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, xhnm, ncol, lchnk, ltrop )
    !-----------------------------------------------------------------------
    !      	... Exp_sol advances the volumetric mixing ratio
    !           forward one time step via the fully explicit
    !           Euler scheme
    !-----------------------------------------------------------------------

    use chem_mods,     only : clscnt1, extcnt, gas_pcnst, clsmap, rxntot
    use ppgrid,        only : pcols, pver
    use mo_prod_loss,  only : exp_prod_loss
    use mo_indprd,     only : indprd
    use shr_kind_mod,  only : r8 => shr_kind_r8
    use cam_history,   only : outfld
    use mo_tracname,   only : solsym
    use cam_logfile,  only : iulog
    use spmd_utils,   only : masterproc

    implicit none
    !-----------------------------------------------------------------------
    !     	... Dummy arguments
    !-----------------------------------------------------------------------
    integer,  intent(in)    ::  ncol                                ! columns in chunck
    integer,  intent(in)    ::  lchnk                               ! chunk id
    real(r8), intent(in)    ::  delt                                ! time step (s)
    real(r8), intent(in)    ::  het_rates(ncol,pver,max(1,gas_pcnst))  ! het rates (1/cm^3/s)
    real(r8), intent(in)    ::  reaction_rates(ncol,pver,rxntot)    ! rxt rates (1/cm^3/s)
    real(r8), intent(in)    ::  extfrc(ncol,pver,extcnt)            ! "external insitu forcing" (1/cm^3/s)
    real(r8), intent(in)    ::  xhnm(ncol,pver)
    integer,  intent(in)    ::  ltrop(pcols)                        ! chemistry troposphere boundary (index)
    real(r8), intent(inout) ::  base_sol(ncol,pver,gas_pcnst)       ! working mixing ratios (vmr)

    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    integer  ::  i, k, l, m
    real(r8), dimension(ncol,pver,clscnt1) :: &
         prod, &
         loss, &
         ind_prd

    real(r8), dimension(ncol,pver) :: wrk

    !-----------------------------------------------------------------------      
    !        ... Put "independent" production in the forcing
    !-----------------------------------------------------------------------      
    call indprd( 1, ind_prd, clscnt1, base_sol, extfrc, &
         reaction_rates, ncol )

    !-----------------------------------------------------------------------      
    !      	... Form F(y)
    !-----------------------------------------------------------------------      
    call exp_prod_loss( prod, loss, base_sol, reaction_rates, het_rates )

    !-----------------------------------------------------------------------      
    !    	... Solve for the mixing ratio at t(n+1)
    !-----------------------------------------------------------------------      
    do m = 1,clscnt1
       l = clsmap(m,1)
       ! apply E90 loss in all levels, including stratosphere
       if (trim(solsym(l)) == 'E90') then
          do i = 1,ncol
             do k = 1,pver
                base_sol(i,k,l)  = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
          end do
       else
          do i = 1,ncol
             !
             ! troposphere
             !
             do k = ltrop(i)+1,pver
                base_sol(i,k,l)  = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
             !
             ! stratospheric decay for chemUCI
             !
             do k = 1,ltrop(i)
                if (uci_CH2O_ndx > 0 .and. trim(solsym(l)) == 'CH2O') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_CH2O_ndx)
                end if
                if (uci_CH3O2_ndx > 0 .and. trim(solsym(l)) == 'CH3O2') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_CH3O2_ndx)
                end if
                if (uci_CH3OOH_ndx > 0 .and. trim(solsym(l)) == 'CH3OOH') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_CH3OOH_ndx)
                end if
                if (uci_PAN_ndx > 0 .and. trim(solsym(l)) == 'PAN') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_PAN_ndx)
                end if
                if (uci_CO_ndx > 0 .and. trim(solsym(l)) == 'CO') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_CO_ndx)
                end if
                if (uci_C2H6_ndx > 0 .and. trim(solsym(l)) == 'C2H6') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_C2H6_ndx)
                end if
                if (uci_C3H8_ndx > 0 .and. trim(solsym(l)) == 'C3H8') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_C3H8_ndx)
                end if
                if (uci_C2H4_ndx > 0 .and. trim(solsym(l)) == 'C2H4') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_C2H4_ndx)
                end if
                if (uci_ROHO2_ndx > 0 .and. trim(solsym(l)) == 'ROHO2') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_ROHO2_ndx)
                end if
                if (uci_CH3COCH3_ndx > 0 .and. trim(solsym(l)) == 'CH3COCH3') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_CH3COCH3_ndx)
                end if
                if (uci_C2H5O2_ndx > 0 .and. trim(solsym(l)) == 'C2H5O2') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_C2H5O2_ndx)
                end if
                if (uci_C2H5OOH_ndx > 0 .and. trim(solsym(l)) == 'C2H5OOH') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_C2H5OOH_ndx)
                end if
                if (uci_CH3CHO_ndx > 0 .and. trim(solsym(l)) == 'CH3CHO') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_CH3CHO_ndx)
                end if
                if (uci_CH3CO3_ndx > 0 .and. trim(solsym(l)) == 'CH3CO3') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_CH3CO3_ndx)
                end if
                if (uci_ISOP_ndx > 0 .and. trim(solsym(l)) == 'ISOP') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_ISOP_ndx)
                end if
                if (uci_ISOPO2_ndx > 0 .and. trim(solsym(l)) == 'ISOPO2') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_ISOPO2_ndx)
                end if
                if (uci_MVKMACR_ndx > 0 .and. trim(solsym(l)) == 'MVKMACR') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_MVKMACR_ndx)
                end if
                if (uci_MVKO2_ndx > 0 .and. trim(solsym(l)) == 'MVKO2') then
                   base_sol(i,k,l)  = base_sol(i,k,l) - delt * reaction_rates(i,k,uci_MVKO2_ndx)
                end if
             end do
          end do
       end if

       wrk(:,:) = (prod(:,:,m) + ind_prd(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHMP', wrk(:,:), ncol, lchnk )
       wrk(:,:) = (loss(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHML', wrk(:,:), ncol, lchnk )
       
    end do

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'exp_sol: diagnostics uci rates'
       write(iulog,*) reaction_rates(1,1,uci_CH2O_ndx),reaction_rates(1,1,uci_CH3O2_ndx),&
         reaction_rates(1,1,uci_CH3OOH_ndx),reaction_rates(1,1,uci_PAN_ndx),reaction_rates(1,1,uci_CO_ndx), &
         reaction_rates(1,1,uci_C2H6_ndx),reaction_rates(1,1,uci_C3H8_ndx),reaction_rates(1,1,uci_C2H4_ndx),&
         reaction_rates(1,1,uci_ROHO2_ndx),reaction_rates(1,1,uci_CH3COCH3_ndx),&
         reaction_rates(1,1,uci_C2H5O2_ndx),reaction_rates(1,1,uci_C2H5OOH_ndx),reaction_rates(1,1,uci_CH3CHO_ndx),&
         reaction_rates(1,1,uci_CH3CO3_ndx),reaction_rates(1,1,uci_ISOP_ndx),reaction_rates(1,1,uci_ISOPO2_ndx),&
         reaction_rates(1,1,uci_MVKMACR_ndx),reaction_rates(1,1,uci_MVKO2_ndx)
    end if
  end subroutine exp_sol

end module mo_exp_sol
