
module mo_exp_sol

  use shr_kind_mod, only : r8 => shr_kind_r8

  private
  public :: exp_sol
  public :: exp_sol_inti

  save

  real(r8), parameter :: r30days = 3.858e-7_r8 ! e-folding decay rate of 30 days
  integer :: uci1_ndx

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
    uci1_ndx     = get_rxt_ndx( 'uci1' )

    if (masterproc) then
       write(iulog,*) ' '
       write(iulog,*) 'exp_sol_inti: diagnostics uci rates indices'
       write(iulog,'(10i5)') uci1_ndx
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
          !
          ! tropospheric solver for non-E90 species
          !
          do i = 1,ncol
             do k = ltrop(i)+1,pver
                base_sol(i,k,l)  = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
          end do
       end if
       !
       ! stratospheric decay for chemUCI
       !
       if (uci1_ndx > 0) then
          if (trim(solsym(l))=='CH2O' .or. trim(solsym(l))=='CH3O2' &
            .or. trim(solsym(l))=='CH3OOH' .or. trim(solsym(l))=='PAN' &
            .or. trim(solsym(l))=='CO' .or. trim(solsym(l))=='C2H6' &
            .or. trim(solsym(l))=='C3H8' .or. trim(solsym(l))=='C2H4' &
            .or. trim(solsym(l))=='ROHO2'.or.trim(solsym(l))=='CH3COCH3' &
            .or. trim(solsym(l))=='C2H5O2'.or.trim(solsym(l))=='C2H5OOH' &
            .or. trim(solsym(l))=='CH3CHO'.or.trim(solsym(l))=='CH3CO3' &
            .or. trim(solsym(l))=='ISOP' .or. trim(solsym(l))=='ISOPO2' &
            .or. trim(solsym(l))=='MVKMACR'.or.trim(solsym(l))=='MVKO2' &
            ) then
             do i = 1,ncol
                do k = 1,ltrop(i)
                   base_sol(i,k,l)  = base_sol(i,k,l) * (1 - delt * r30days)
                end do
             end do
          end if
       end if

       wrk(:,:) = (prod(:,:,m) + ind_prd(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHMP', wrk(:,:), ncol, lchnk )
       wrk(:,:) = (loss(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHML', wrk(:,:), ncol, lchnk )
       
    end do

  end subroutine exp_sol

end module mo_exp_sol
