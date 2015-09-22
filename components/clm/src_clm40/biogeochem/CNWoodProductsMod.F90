module CNWoodProductsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNWoodProductsMod
!
! !DESCRIPTION:
! Calculate loss fluxes from wood products pools, and update product pool state variables
!
! !USES:
    use decompMod   , only : get_proc_bounds
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varcon  , only: istsoil
    use spmdMod     , only: masterproc
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CNWoodProducts
!
! !REVISION HISTORY:
! 5/20/2009: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNWoodProducts
!
! !INTERFACE:
subroutine CNWoodProducts(num_soilc, filter_soilc)
!
! !DESCRIPTION:
! Update all loss fluxes from wood product pools, and update product pool state variables
! for both loss and gain terms.  Gain terms are calculated in pftdyn_cnbal() for gains associated
! with changes in landcover, and in CNHarvest(), for gains associated with wood harvest.
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varctl, only : use_c13
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 5/21/09: Created by Peter Thornton
!
! !LOCAL VARIABLES:

   integer :: fc        ! lake filter indices
   integer :: c         ! indices
   real(r8):: dt        ! time step (seconds)
   type(column_type),   pointer :: cptr         ! pointer to column derived subtype
   real(r8) :: kprod10       ! decay constant for 10-year product pool
   real(r8) :: kprod100      ! decay constant for 100-year product pool

!EOP
!-----------------------------------------------------------------------

   cptr => col
	
   ! calculate column-level losses from product pools
	! the following (1/s) rate constants result in ~90% loss of initial state over 10 and 100 years,
	! respectively, using a discrete-time fractional decay algorithm.
	kprod10 = 7.2e-9
	kprod100 = 7.2e-10

   do fc = 1,num_soilc
      c = filter_soilc(fc)

		! calculate fluxes (1/sec)
		ccf%prod10c_loss(c)    = ccs%prod10c(c)    * kprod10
		ccf%prod100c_loss(c)   = ccs%prod100c(c)   * kprod100
                if (use_c13) then
                   cc13f%prod10c_loss(c)  = cc13s%prod10c(c)  * kprod10
                   cc13f%prod100c_loss(c) = cc13s%prod100c(c) * kprod100
                end if
		cnf%prod10n_loss(c)    = cns%prod10n(c)    * kprod10
		cnf%prod100n_loss(c)   = cns%prod100n(c)   * kprod100
	end do

   ! set time steps
   dt = real( get_step_size(), r8 )

   ! update wood product state variables
   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! column-level fluxes

      ! fluxes into wood product pools, from landcover change
      ccs%prod10c(c)    = ccs%prod10c(c)    + ccf%dwt_prod10c_gain(c)*dt
      ccs%prod100c(c)   = ccs%prod100c(c)   + ccf%dwt_prod100c_gain(c)*dt
      if (use_c13) then
         cc13s%prod10c(c)  = cc13s%prod10c(c)  + cc13f%dwt_prod10c_gain(c)*dt
         cc13s%prod100c(c) = cc13s%prod100c(c) + cc13f%dwt_prod100c_gain(c)*dt
      end if
      cns%prod10n(c)    = cns%prod10n(c)    + cnf%dwt_prod10n_gain(c)*dt
      cns%prod100n(c)   = cns%prod100n(c)   + cnf%dwt_prod100n_gain(c)*dt

      ! fluxes into wood product pools, from harvest
      ccs%prod10c(c)    = ccs%prod10c(c)    + ccf%hrv_deadstemc_to_prod10c(c)*dt
      ccs%prod100c(c)   = ccs%prod100c(c)   + ccf%hrv_deadstemc_to_prod100c(c)*dt
      if (use_c13) then
         cc13s%prod10c(c)  = cc13s%prod10c(c)  + cc13f%hrv_deadstemc_to_prod10c(c)*dt
         cc13s%prod100c(c) = cc13s%prod100c(c) + cc13f%hrv_deadstemc_to_prod100c(c)*dt
      end if
      cns%prod10n(c)    = cns%prod10n(c)    + cnf%hrv_deadstemn_to_prod10n(c)*dt
      cns%prod100n(c)   = cns%prod100n(c)   + cnf%hrv_deadstemn_to_prod100n(c)*dt
     
      ! fluxes out of wood product pools, from decomposition
      ccs%prod10c(c)    = ccs%prod10c(c)    - ccf%prod10c_loss(c)*dt
      ccs%prod100c(c)   = ccs%prod100c(c)   - ccf%prod100c_loss(c)*dt
      if (use_c13) then
         cc13s%prod10c(c)  = cc13s%prod10c(c)  - cc13f%prod10c_loss(c)*dt
         cc13s%prod100c(c) = cc13s%prod100c(c) - cc13f%prod100c_loss(c)*dt
      end if
      cns%prod10n(c)    = cns%prod10n(c)    - cnf%prod10n_loss(c)*dt
      cns%prod100n(c)   = cns%prod100n(c)   - cnf%prod100n_loss(c)*dt
 
   end do ! end of column loop

end subroutine CNWoodProducts
!-----------------------------------------------------------------------

end module CNWoodProductsMod
