
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
       if (trim(solsym(l)) == 'E90') then
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
       ! apply brc -> pom decay to troposphere
       ! apply ext emission to whole atmosphere in case of large fires
       elseif (trim(solsym(l)) == 'brc_a4' .or. trim(solsym(l)) == 'brc_a1' .or. trim(solsym(l)) == 'brc_a3') then
               ! V2-like explicit equation 
          do i = 1,ncol
             do k = 1,ltrop(i)
             ! above tropopause, allow BrC emission, but do not consider BrC decay for now because of low radical concentration.
                chem_loss(i,k,l) = 0.0_r8
                chem_prod(i,k,l) = ind_prd(i,k,m) 
                base_sol(i,k,l) = base_sol(i,k,l) + delt * ind_prd(i,k,m)  
             end do               
             do k = ltrop(i)+1,pver ! brown carbon oxidation to form pom  only in troposphere
                chem_loss(i,k,l) = -loss(i,k,m)
                chem_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                base_sol(i,k,l) = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
          end do
       elseif (trim(solsym(l)) == 'pom_a4' .or. trim(solsym(l)) == 'pom_a1' .or. trim(solsym(l)) == 'pom_a3') then
               ! V2-like explicit equation 
          do i = 1,ncol
             do k = 1,ltrop(i)
             ! above tropopause, emission only
                chem_loss(i,k,l) = 0.0_r8
                chem_prod(i,k,l) = ind_prd(i,k,m)
                base_sol(i,k,l) = base_sol(i,k,l) + delt * ind_prd(i,k,m)
             end do
             do k = ltrop(i)+1,pver ! brown carbon oxidation to form pom only in troposphere 
                chem_loss(i,k,l) = -loss(i,k,m)
                chem_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                base_sol(i,k,l) = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
          end do
       ! no reaction for now but in case emission in stratosphere 
       elseif (trim(solsym(l)) == 'bc_a4' .or. trim(solsym(l)) == 'bc_a1' .or. trim(solsym(l)) == 'bc_a3') then
               ! V2-like explicit equation
          do i = 1,ncol
             do k = 1,ltrop(i)
             ! above tropopause, emission only
                chem_loss(i,k,l) = 0.0_r8
                chem_prod(i,k,l) = ind_prd(i,k,m)
                base_sol(i,k,l) = base_sol(i,k,l) + delt * ind_prd(i,k,m)
             end do 
             do k = ltrop(i)+1,pver
                chem_loss(i,k,l) = -loss(i,k,m)
                chem_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                base_sol(i,k,l) = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
          end do
       ! no reaction for now but in case emission in stratosphere
       elseif (trim(solsym(l)) == 'num_a4' .or. trim(solsym(l)) == 'num_a1' .or. trim(solsym(l)) == 'num_a3') then
               ! V2-like explicit equation
          do i = 1,ncol
             do k = 1,ltrop(i)
             ! above tropopause, emission only
                chem_loss(i,k,l) = 0.0_r8
                chem_prod(i,k,l) = ind_prd(i,k,m)
                base_sol(i,k,l) = base_sol(i,k,l) + delt * ind_prd(i,k,m)
             end do
             do k = ltrop(i)+1,pver
                chem_loss(i,k,l) = -loss(i,k,m)
                chem_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                base_sol(i,k,l) = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
          end do
       elseif (trim(solsym(l)) == 'H2SO4' .or. trim(solsym(l)) == 'SO2' .or. trim(solsym(l)) == 'DMS') then
               ! V2-like explicit equation is used to solve H2SO4 and SO2 due to dead zero values
          do i = 1,ncol
             do k = 1,pver
                chem_loss(i,k,l) = -loss(i,k,m)
                chem_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                base_sol(i,k,l) = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
          end do
       else
          do i = 1,ncol
             do k = 1,ltrop(i) ! for emissions in the stratosphere: NO2 as high as 15 km
                chem_prod(i,k,l) = ind_prd(i,k,m)
                chem_loss(i,k,l) = 0.0_r8
                base_sol(i,k,l) = base_sol(i,k,l) + delt*ind_prd(i,k,m)
             end do
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
    call indprd( 1, ind_prd, clscnt1, base_sol_reset, extfrc, &
         diags_reaction_rates, ncol )
    call exp_prod_loss( prod, loss, base_sol_reset, diags_reaction_rates, het_rates )
    do m = 1,clscnt1
       l = clsmap(m,1)
       ! v2 like treatment
       if (trim(solsym(l)) == 'H2SO4' .or. trim(solsym(l)) == 'SO2' .or. trim(solsym(l)) == 'DMS') then
           do i = 1,ncol
              do k = ltrop(i)+1,pver
                 chemmp_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                 chemmp_loss(i,k,l) = -loss(i,k,m)
               end do
           end do
       elseif (trim(solsym(l)) == 'brc_a1' .or. trim(solsym(l)) == 'brc_a4' .or. trim(solsym(l)) == 'brc_a3') then
           do i = 1,ncol
              do k = ltrop(i)+1,pver
                 chemmp_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                 chemmp_loss(i,k,l) = -loss(i,k,m)
               end do
           end do
       elseif (trim(solsym(l)) == 'pom_a1' .or. trim(solsym(l)) == 'pom_a4' .or. trim(solsym(l)) == 'pom_a3') then 
           do i = 1,ncol
              do k = ltrop(i)+1,pver
                 chemmp_prod(i,k,l) = prod(i,k,m)+ind_prd(i,k,m)
                 chemmp_loss(i,k,l) = -loss(i,k,m)
               end do
           end do
       elseif (trim(solsym(l)) == 'bc_a1' .or. trim(solsym(l)) == 'bc_a4' .or. trim(solsym(l)) == 'bc_a3') then
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
       end if
    end do

  end subroutine exp_sol

end module mo_exp_sol
