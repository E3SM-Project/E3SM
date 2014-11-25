module CNWoodProductsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate loss fluxes from wood products pools, and update product pool state variables
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use decompMod           , only : get_proc_bounds
  use spmdMod             , only : masterproc
  use landunit_varcon     , only : istsoil
  use clm_time_manager    , only : get_step_size
  use clm_varctl          , only : use_c13, use_c14
  use CNCarbonStateType   , only : carbonstate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNWoodProducts
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNWoodProducts(num_soilc, filter_soilc, &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
       carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, nitrogenflux_vars)
    !
    ! !DESCRIPTION:
    ! Update all loss fluxes from wood product pools, and update product pool state variables
    ! for both loss and gain terms.  Gain terms are calculated in pftdyn_cnbal() for gains associated
    ! with changes in landcover, and in CNHarvest(), for gains associated with wood harvest.
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(carbonstate_type)   , intent(in)    :: c13_carbonstate_vars
    type(carbonstate_type)   , intent(in)    :: c14_carbonstate_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: fc        ! lake filter indices
    integer :: c         ! indices
    real(r8):: dt        ! time step (seconds)
    real(r8) :: kprod10       ! decay constant for 10-year product pool
    real(r8) :: kprod100      ! decay constant for 100-year product pool
    !-----------------------------------------------------------------------


    ! calculate column-level losses from product pools
    ! the following (1/s) rate constants result in ~90% loss of initial state over 10 and 100 years,
    ! respectively, using a discrete-time fractional decay algorithm.
    kprod10 = 7.2e-9
    kprod100 = 7.2e-10

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! calculate fluxes (1/sec)
       carbonflux_vars%prod10c_loss_col(c)    = carbonstate_vars%prod10c_col(c)    * kprod10
       carbonflux_vars%prod100c_loss_col(c)   = carbonstate_vars%prod100c_col(c)   * kprod100

       if ( use_c13 ) then
          c13_carbonflux_vars%prod10c_loss_col(c)  = c13_carbonstate_vars%prod10c_col(c)  * kprod10
          c13_carbonflux_vars%prod100c_loss_col(c) = c13_carbonstate_vars%prod100c_col(c) * kprod100
       endif

       if ( use_c14 ) then
          c14_carbonflux_vars%prod10c_loss_col(c)  = c14_carbonstate_vars%prod10c_col(c)  * kprod10
          c14_carbonflux_vars%prod100c_loss_col(c) = c14_carbonstate_vars%prod100c_col(c) * kprod100
       endif

       nitrogenflux_vars%prod10n_loss_col(c)    = nitrogenstate_vars%prod10n_col(c)    * kprod10
       nitrogenflux_vars%prod100n_loss_col(c)   = nitrogenstate_vars%prod100n_col(c)   * kprod100
    end do

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! update wood product state variables
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! fluxes into wood product pools, from landcover change
       carbonstate_vars%prod10c_col(c)    = carbonstate_vars%prod10c_col(c)    + &
            carbonflux_vars%dwt_prod10c_gain_col(c)*dt
       carbonstate_vars%prod100c_col(c)   = carbonstate_vars%prod100c_col(c)   + &
            carbonflux_vars%dwt_prod100c_gain_col(c)*dt

       if ( use_c13 ) then
          c13_carbonstate_vars%prod10c_col(c)  = c13_carbonstate_vars%prod10c_col(c)  + &
               c13_carbonflux_vars%dwt_prod10c_gain_col(c) *dt
          c13_carbonstate_vars%prod100c_col(c) = c13_carbonstate_vars%prod100c_col(c) + &
               c13_carbonflux_vars%dwt_prod100c_gain_col(c)*dt
       endif

       if ( use_c14 ) then
          c14_carbonstate_vars%prod10c_col(c)  = c14_carbonstate_vars%prod10c_col(c)  + &
               c14_carbonflux_vars%dwt_prod10c_gain_col(c) *dt
          c14_carbonstate_vars%prod100c_col(c) = c14_carbonstate_vars%prod100c_col(c) + &
               c14_carbonflux_vars%dwt_prod100c_gain_col(c)*dt
       endif

       nitrogenstate_vars%prod10n_col(c)    = nitrogenstate_vars%prod10n_col(c)    + &
            nitrogenflux_vars%dwt_prod10n_gain_col(c)*dt
       nitrogenstate_vars%prod100n_col(c)   = nitrogenstate_vars%prod100n_col(c)   + &
            nitrogenflux_vars%dwt_prod100n_gain_col(c)*dt

       ! fluxes into wood product pools, from harvest
       carbonstate_vars%prod10c_col(c)    = carbonstate_vars%prod10c_col(c)    + &
            carbonflux_vars%hrv_deadstemc_to_prod10c_col(c)*dt
       carbonstate_vars%prod100c_col(c)   = carbonstate_vars%prod100c_col(c)   + &
            carbonflux_vars%hrv_deadstemc_to_prod100c_col(c)*dt

       if ( use_c13 ) then
          c13_carbonstate_vars%prod10c_col(c)  = c13_carbonstate_vars%prod10c_col(c)  + &
               c13_carbonflux_vars%hrv_deadstemc_to_prod10c_col(c)*dt
          c13_carbonstate_vars%prod100c_col(c) = c13_carbonstate_vars%prod100c_col(c) + &
               c13_carbonflux_vars%hrv_deadstemc_to_prod100c_col(c)*dt
       endif

       if ( use_c14 ) then
          c14_carbonstate_vars%prod10c_col(c)  = c14_carbonstate_vars%prod10c_col(c)  + &
               c14_carbonflux_vars%hrv_deadstemc_to_prod10c_col(c)*dt
          c14_carbonstate_vars%prod100c_col(c) = c14_carbonstate_vars%prod100c_col(c) + &
               c14_carbonflux_vars%hrv_deadstemc_to_prod100c_col(c)*dt
       endif

       nitrogenstate_vars%prod10n_col(c)    = nitrogenstate_vars%prod10n_col(c)    + &
            nitrogenflux_vars%hrv_deadstemn_to_prod10n_col(c)*dt
       nitrogenstate_vars%prod100n_col(c)   = nitrogenstate_vars%prod100n_col(c)   + &
            nitrogenflux_vars%hrv_deadstemn_to_prod100n_col(c)*dt

       ! fluxes out of wood product pools, from decomposition
       carbonstate_vars%prod10c_col(c)    = carbonstate_vars%prod10c_col(c)    - &
            carbonflux_vars%prod10c_loss_col(c)*dt
       carbonstate_vars%prod100c_col(c)   = carbonstate_vars%prod100c_col(c)   - &
            carbonflux_vars%prod100c_loss_col(c)*dt

       if ( use_c13 ) then
          c13_carbonstate_vars%prod10c_col(c)  = c13_carbonstate_vars%prod10c_col(c)  - &
               c13_carbonflux_vars%prod10c_loss_col(c)*dt
          c13_carbonstate_vars%prod100c_col(c) = c13_carbonstate_vars%prod100c_col(c) - &
               c13_carbonflux_vars%prod100c_loss_col(c)*dt
       endif

       if ( use_c14 ) then
          c14_carbonstate_vars%prod10c_col(c)  = c14_carbonstate_vars%prod10c_col(c)  - &
               c14_carbonflux_vars%prod10c_loss_col(c)*dt
          c14_carbonstate_vars%prod100c_col(c) = c14_carbonstate_vars%prod100c_col(c) - &
               c14_carbonflux_vars%prod100c_loss_col(c)*dt
       endif

       nitrogenstate_vars%prod10n_col(c)    = nitrogenstate_vars%prod10n_col(c)    - &
            nitrogenflux_vars%prod10n_loss_col(c)*dt
       nitrogenstate_vars%prod100n_col(c)   = nitrogenstate_vars%prod100n_col(c)   - &
            nitrogenflux_vars%prod100n_loss_col(c)*dt

    end do ! end of column loop

  end subroutine CNWoodProducts

end module CNWoodProductsMod
