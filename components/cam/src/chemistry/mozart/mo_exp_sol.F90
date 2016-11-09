
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
       do i = 1,ncol
          do k = ltrop(i)+1,pver
             base_sol(i,k,l)  = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
          end do
       end do

       wrk(:,:) = (prod(:,:,m) + ind_prd(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHMP', wrk(:,:), ncol, lchnk )
       wrk(:,:) = (loss(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHML', wrk(:,:), ncol, lchnk )
       
    end do

  end subroutine exp_sol

end module mo_exp_sol
