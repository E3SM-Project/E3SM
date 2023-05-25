
module mo_exp_sol

  private
  public :: exp_sol
  public :: exp_sol_inti

contains

  subroutine exp_sol_inti

    use mo_tracname, only : solsym
    use chem_mods,   only : clscnt1, clsmap
    use ppgrid,      only : pver
    use cam_history, only : addfld, add_default

    implicit none

    integer :: i,j

    do i = 1,clscnt1

       j = clsmap(i,1)
       call addfld( trim(solsym(j))//'_CHMP', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate' )
       call addfld( trim(solsym(j))//'_CHML', (/ 'lev' /), 'I', '/cm3/s', 'chemical loss rate' )

    enddo
  end subroutine exp_sol_inti


  subroutine exp_sol( base_sol, reaction_rates, het_rates, extfrc, delt, xhnm, ncol, lchnk, ltrop, diags_reaction_rates, chem_prod, chem_loss, chemmp_prod, chemmp_loss )
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
    use mo_chem_utls,      only :  get_inv_ndx
    use cam_logfile,      only : iulog
    use chem_mods,     only : nfs 
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
    real(r8), intent(in)    ::  diags_reaction_rates(ncol,pver,rxntot)    ! rxt rates (1/cm^3/s)
    real(r8), intent(out)   ::  chem_prod(ncol,pver,gas_pcnst)      ! production rate (vmr/delt)
    real(r8), intent(out)   ::  chem_loss(ncol,pver,gas_pcnst)      ! loss rate (vmr/delt)
    real(r8), intent(out)   ::  chemmp_prod(ncol,pver,gas_pcnst)      ! production rate (vmr/delt)
    real(r8), intent(out)   ::  chemmp_loss(ncol,pver,gas_pcnst)      ! loss rate (vmr/delt)

    !-----------------------------------------------------------------------
    !     	... Local variables
    !-----------------------------------------------------------------------
    integer  ::  i, k, l, m
    real(r8), dimension(ncol,pver,clscnt1) :: &
         prod, &
         loss, &
         ind_prd

    real(r8), dimension(ncol,pver) :: wrk
    real(r8), dimension(ncol,pver,gas_pcnst) :: base_sol_reset

    chem_prod(:,:,:) = 0._r8
    chem_loss(:,:,:) = 0._r8
    chemmp_prod(:,:,:) = 0._r8
    chemmp_loss(:,:,:) = 0._r8
    !-----------------------------------------------------------------------      
    !        ... Put "independent" production in the forcing
    !-----------------------------------------------------------------------      
    base_sol_reset = base_sol
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
       if (trim(solsym(l)) == 'E90' ) then
          do i = 1,ncol
             do k = 1,pver
                ! change the old equation
                ! base_sol(i,k,l) = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
                ! to use the e-folding decay to handle the loss part to ensure non-negative solutions
                ! loss/base_sol is the loss frequency (dx/x).
                ! base_sol(i,k,l)*exp(-delt*loss(i,k,m)/base_sol(i,k,l)) is what left after decay
                ! and delt*(prod(i,k,m)+ind_prd(i,k,m)) is the production.
                base_sol(i,k,l) = base_sol(i,k,l)*exp(-delt*loss(i,k,m)/base_sol(i,k,l)) + delt*(prod(i,k,m)+ind_prd(i,k,m))
             end do
          end do
       ! SO2 and H2SO4 can be dead zeros due to aerosol processes
       ! use a different equation for them to avoid debug built issues
#if (defined MODAL_AERO_5MODE)  
       ! for MAM5, H2SO4, SO2, and DMS needs to be solved in the stratosphere     
       elseif (trim(solsym(l)) == 'H2SO4' .or. trim(solsym(l)) == 'SO2') then
         ! V2-like explicit equation is used to solve H2SO4 and SO2 due to dead zero values
         do i = 1,ncol
             do k = 1,pver
                chem_loss(i,k,l) = -loss(i,k,m)
                chem_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                base_sol(i,k,l) = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
         end do
       elseif (trim(solsym(l)) == 'DMS') then
         ! DMS doesn't have dead zero value issue
         do i = 1,ncol
             do k = 1,pver
                chem_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                chem_loss(i,k,l) = (base_sol(i,k,l)*exp(-delt*loss(i,k,m)/base_sol(i,k,l)) - base_sol(i,k,l))/delt
                base_sol(i,k,l) = base_sol(i,k,l)*exp(-delt*loss(i,k,m)/base_sol(i,k,l)) + delt*(prod(i,k,m)+ind_prd(i,k,m))
             end do
          end do
#else
       elseif (trim(solsym(l)) == 'H2SO4' .or. trim(solsym(l)) == 'SO2') then
          do i = 1,ncol
              do k = ltrop(i)+1,pver
                 chem_loss(i,k,l) = -loss(i,k,m)
                 chem_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m) 
                 base_sol(i,k,l) = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
              end do
          end do
#endif         
       else
          do i = 1,ncol
             do k = ltrop(i)+1,pver
                chem_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                chem_loss(i,k,l) = (base_sol(i,k,l)*exp(-delt*loss(i,k,m)/base_sol(i,k,l)) - base_sol(i,k,l))/delt
                base_sol(i,k,l) = base_sol(i,k,l)*exp(-delt*loss(i,k,m)/base_sol(i,k,l)) + delt*(prod(i,k,m)+ind_prd(i,k,m))
             end do
          end do
       end if

       wrk(:,:) = (prod(:,:,m) + ind_prd(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHMP', wrk(:,:), ncol, lchnk )
       wrk(:,:) = (loss(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHML', wrk(:,:), ncol, lchnk )
       
    end do
!----------- for UCIchem diagnostics ------------
    ! Note: base_sol_rest and diags_reaction_rates are passed to indprd and exp_prod_loss
    !       so the reset tendency for tropopause folds (i.e., stratospheric boxes below the tropopause)
    !       and the explicit solver tendency are included in chemmp_prod and chemmp_loss below.
    call indprd( 1, ind_prd, clscnt1, base_sol_reset, extfrc, &
         diags_reaction_rates, ncol )
    call exp_prod_loss( prod, loss, base_sol_reset, diags_reaction_rates, het_rates )
    do m = 1,clscnt1
       l = clsmap(m,1)
       if (trim(solsym(l)) == 'H2SO4' .or. trim(solsym(l)) == 'SO2') then
           do i = 1,ncol
              do k = ltrop(i)+1,pver
                 chemmp_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                 chemmp_loss(i,k,l) = -loss(i,k,m)
              end do
           end do
       else
           do i = 1,ncol
              do k = ltrop(i)+1,pver
                 chemmp_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                 chemmp_loss(i,k,l) = (base_sol_reset(i,k,l)*exp(-delt*loss(i,k,m)/base_sol_reset(i,k,l)) - base_sol_reset(i,k,l))/delt
              end do
           end do
       endif
    end do

  end subroutine exp_sol

end module mo_exp_sol
