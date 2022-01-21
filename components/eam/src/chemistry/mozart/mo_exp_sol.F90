
module mo_exp_sol

!LXu@09/2021+++
  use spmd_utils,       only : masterproc
  use cam_logfile,      only : iulog
!LXu@09/2021---

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

       if (trim(solsym(j)) == 'CO') then
       call addfld( trim(solsym(j))//'_CHMP_jch2o_a', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate by photolysis(CH2O+hv->CO+2HO2)' )
       call addfld( trim(solsym(j))//'_CHMP_jch2o_b', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate by photolysis(CH2O+hv->CO+H2)' )
       call addfld( trim(solsym(j))//'_CHMP_jch3cho', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate by photolysis(CH3CHO+hv)' )
       call addfld( trim(solsym(j))//'_CHMP_jacet', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate by photolysis(CH3COCH3+hv)' )
       call addfld( trim(solsym(j))//'_CHMP_jmvk', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate by photolysis(MVKMACR+hv)' )
       call addfld( trim(solsym(j))//'_CHMP_c2h4_o3', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate from C2H4+O3' )
       call addfld( trim(solsym(j))//'_CHMP_ch2o_oh', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate from CH2O+OH' )
       call addfld( trim(solsym(j))//'_CHMP_ch3cho_oh', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate by CH3CHO+OH' )
       call addfld( trim(solsym(j))//'_CHMP_isopo2_ho2', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate from ISOPO2+HO2' )
       call addfld( trim(solsym(j))//'_CHMP_mvko2_ho2', (/ 'lev' /), 'I', '/cm3/s', 'chemical production rate from MVKO2+HO2' )
       call addfld( trim(solsym(j))//'_CHML_co_oh', (/ 'lev' /), 'I', '/cm3/s', 'chemical loss rate (CO+OH)' )
       call addfld( trim(solsym(j))//'_CHML_co_oh_m', (/ 'lev' /), 'I', '/cm3/s', 'chemical loss rate (CO+OH+M)' )
       call addfld( trim(solsym(j))//'_CHML_rainout', (/ 'lev' /), 'I', '/cm3/s', 'chemical loss rate (raintout)' )
       end if
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
    real(r8), intent(in)    ::  reaction_rates(ncol,pver,rxntot)    ! reaction_rates rates (1/cm^3/s)
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

!added by LXu@09/2021+++
    real(r8), dimension(ncol,pver,clscnt1) :: &
         ind_prd01, &
         ind_prd02, &
         ind_prd03, &
         ind_prd04, &
         ind_prd05, &
         ind_prd06, &
         ind_prd07, &
         ind_prd08, &
         ind_prd09, &
         ind_prd10, &
         loss01, &
         loss02, &
         loss03
	 
    real(r8), dimension(ncol,pver) :: wrk1,wrk2,wrk3,wrk4,wrk5,wrk6,wrk7,wrk8,wrk9,wrk10
!added by LXu@09/2021---
!    if (masterproc) write(iulog,*) 'mo_exp_sol write ',clscnt1, extcnt, gas_pcnst, clsmap, rxntot
    !-----------------------------------------------------------------------      
    !        ... Put "independent" production in the forcing
    !-----------------------------------------------------------------------      
    call indprd( 1, ind_prd, clscnt1, base_sol, extfrc, &
         reaction_rates, ncol )

!added by LXu@09/2021+++
!chmp_photo
!         ind_prd01(:,:,1) = (reaction_rates(:,:,4) +reaction_rates(:,:,5))*base_sol(:,:,5) &
!;	                  +  reaction_rates(:,:,16)*base_sol(:,:,23) &
!			  +  0.330_r8*reaction_rates(:,:,18)*base_sol(:,:,20) &
! 			  +  1.500_r8*reaction_rates(:,:,19)*base_sol(:,:,27) 
         ind_prd01(:,:,1) =          reaction_rates(:,:,4) *base_sol(:,:,5) 
         ind_prd02(:,:,1) =          reaction_rates(:,:,5) *base_sol(:,:,5) 
         ind_prd03(:,:,1) =          reaction_rates(:,:,16)*base_sol(:,:,23) 
         ind_prd04(:,:,1) = 0.330_r8*reaction_rates(:,:,18)*base_sol(:,:,20)
         ind_prd05(:,:,1) = 1.500_r8*reaction_rates(:,:,19)*base_sol(:,:,27)
!chmp_c2h4_o3
         ind_prd06(:,:,1) =          reaction_rates(:,:,30)*base_sol(:,:,18)*base_sol(:,:,1) 
!chmp_ch2o
         ind_prd07(:,:,1) =          reaction_rates(:,:,33)*base_sol(:,:, 5)*base_sol(:,:,2) 
         ind_prd08(:,:,1) = 0.500_r8*reaction_rates(:,:,67)*base_sol(:,:,23)*base_sol(:,:,2) 
         ind_prd09(:,:,1) = 2.000_r8*reaction_rates(:,:,78)*base_sol(:,:,26)*base_sol(:,:,3) 
         ind_prd10(:,:,1) =          reaction_rates(:,:,83)*base_sol(:,:,28)*base_sol(:,:,3) 
!added by LXu@09/2021---
    !-----------------------------------------------------------------------      
    !      	... Form F(y)
    !-----------------------------------------------------------------------      
    call exp_prod_loss( prod, loss, base_sol, reaction_rates, het_rates )

!added by LXu@09/2021+++
!         loss(:,:,1) = ((rxt(:,:,23) +rxt(:,:,24))* y(:,:,2) + het_rates(:,:,15)) &
!                 * y(:,:,15)
!         prod(:,:,1) =.330_r8*rxt(:,:,18)*y(:,:,20)
         loss01(:,:,1) = reaction_rates(:,:,23)* base_sol(:,:,15) * base_sol(:,:,2)
         loss02(:,:,1) = reaction_rates(:,:,24)* base_sol(:,:,15) * base_sol(:,:,2)
         loss03(:,:,1) = het_rates(:,:,15)* base_sol(:,:,15) 
!added by LXu@09/2021---
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
             do k = ltrop(i)+1,pver
                base_sol(i,k,l)  = base_sol(i,k,l) + delt * (prod(i,k,m) + ind_prd(i,k,m) - loss(i,k,m))
             end do
          end do
       end if

       wrk(:,:) = (prod(:,:,m) + ind_prd(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHMP', wrk(:,:), ncol, lchnk )
       wrk(:,:) = (loss(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHML', wrk(:,:), ncol, lchnk )

!added by LXu@09/2021+++
       
       if (trim(solsym(l)) == 'CO') then
       wrk1(:,:)  = (ind_prd01(:,:,m))*xhnm
       wrk2(:,:)  = (ind_prd02(:,:,m))*xhnm
       wrk3(:,:)  = (ind_prd03(:,:,m))*xhnm
       wrk4(:,:)  = (ind_prd04(:,:,m))*xhnm
       wrk5(:,:)  = (ind_prd05(:,:,m))*xhnm
       wrk6(:,:)  = (ind_prd06(:,:,m))*xhnm
       wrk7(:,:)  = (ind_prd07(:,:,m))*xhnm
       wrk8(:,:)  = (ind_prd08(:,:,m))*xhnm
       wrk9(:,:)  = (ind_prd09(:,:,m))*xhnm
       wrk10(:,:) = (ind_prd10(:,:,m))*xhnm
!       call outfld( trim(solsym(l))//'_CHMP_photo', wrk1(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_jch2o_a',    wrk1(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_jch2o_b',    wrk2(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_jch3cho',    wrk3(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_jacet',      wrk4(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_jmvk',       wrk5(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_c2h4_o3',    wrk6(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_ch2o_oh',    wrk7(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_ch3cho_oh',  wrk8(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_isopo2_ho2', wrk9(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHMP_mvko2_ho2', wrk10(:,:), ncol, lchnk )

       wrk1(:,:) = (loss01(:,:,m))*xhnm
       wrk2(:,:) = (loss02(:,:,m))*xhnm
       wrk3(:,:) = (loss03(:,:,m))*xhnm
       call outfld( trim(solsym(l))//'_CHML_co_oh',   wrk1(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHML_co_oh_m', wrk2(:,:), ncol, lchnk )
       call outfld( trim(solsym(l))//'_CHML_rainout', wrk3(:,:), ncol, lchnk )
       end if
!added by LXu@09/2021---
    end do

  end subroutine exp_sol

end module mo_exp_sol
