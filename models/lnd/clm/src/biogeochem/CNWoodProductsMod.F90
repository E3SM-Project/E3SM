module CNWoodProductsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate loss fluxes from wood products pools, and update product pool state variables
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use decompMod               , only : get_proc_bounds
  use spmdMod                 , only : masterproc
  use landunit_varcon         , only : istsoil
  use clm_time_manager        , only : get_step_size
  use clm_varctl              , only : use_c13, use_c14
  use CNVegCarbonStateType   , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType    , only : cnveg_carbonflux_type
  use CNVegnitrogenstateType , only : cnveg_nitrogenstate_type
  use CNVegnitrogenfluxType  , only : cnveg_nitrogenflux_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNWoodProducts
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNWoodProducts(num_soilc, filter_soilc, &
       cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
       cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst, &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst)
    !
    ! !DESCRIPTION:
    ! Update all loss fluxes from wood product pools, and update product pool state variables
    ! for both loss and gain terms.  Gain terms are calculated in pftdyn_cnbal() for gains associated
    ! with changes in landcover, and in CNHarvest(), for gains associated with wood harvest.
    !
    ! !ARGUMENTS:
    integer                        , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                        , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnveg_carbonstate_type)   , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(in)    :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)   , intent(in)    :: c14_cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)    , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type) , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)  , intent(inout) :: cnveg_nitrogenflux_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fc       ! lake filter indices
    integer  :: c        ! indices
    real(r8) :: dt       ! time step (seconds)
    real(r8) :: kprod10  ! decay constant for 10-year product pool
    real(r8) :: kprod100 ! decay constant for 100-year product pool
    !-----------------------------------------------------------------------


    ! calculate column-level losses from product pools
    ! the following (1/s) rate constants result in ~90% loss of initial state over 10 and 100 years,
    ! respectively, using a discrete-time fractional decay algorithm.
    kprod10 = 7.2e-9
    kprod100 = 7.2e-10

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! calculate fluxes (1/sec)
       cnveg_carbonflux_inst%prod10c_loss_col(c)    = cnveg_carbonstate_inst%prod10c_col(c)    * kprod10
       cnveg_carbonflux_inst%prod100c_loss_col(c)   = cnveg_carbonstate_inst%prod100c_col(c)   * kprod100

       if ( use_c13 ) then
          cnveg_carbonflux_inst%prod10c_loss_col(c)  = cnveg_carbonstate_inst%prod10c_col(c)  * kprod10
          cnveg_carbonflux_inst%prod100c_loss_col(c) = cnveg_carbonstate_inst%prod100c_col(c) * kprod100
       endif

       if ( use_c14 ) then
          cnveg_carbonflux_inst%prod10c_loss_col(c)  = cnveg_carbonstate_inst%prod10c_col(c)  * kprod10
          cnveg_carbonflux_inst%prod100c_loss_col(c) = cnveg_carbonstate_inst%prod100c_col(c) * kprod100
       endif

       cnveg_nitrogenflux_inst%prod10n_loss_col(c)    = cnveg_nitrogenstate_inst%prod10n_col(c)    * kprod10
       cnveg_nitrogenflux_inst%prod100n_loss_col(c)   = cnveg_nitrogenstate_inst%prod100n_col(c)   * kprod100
    end do

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! update wood product state variables
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! fluxes into wood product pools, from landcover change
       cnveg_carbonstate_inst%prod10c_col(c)    = cnveg_carbonstate_inst%prod10c_col(c)    + &
            cnveg_carbonflux_inst%dwt_prod10c_gain_col(c)*dt
       cnveg_carbonstate_inst%prod100c_col(c)   = cnveg_carbonstate_inst%prod100c_col(c)   + &
            cnveg_carbonflux_inst%dwt_prod100c_gain_col(c)*dt

       if ( use_c13 ) then
          cnveg_carbonstate_inst%prod10c_col(c)  = cnveg_carbonstate_inst%prod10c_col(c)  + &
               cnveg_carbonflux_inst%dwt_prod10c_gain_col(c) *dt
          cnveg_carbonstate_inst%prod100c_col(c) = cnveg_carbonstate_inst%prod100c_col(c) + &
               cnveg_carbonflux_inst%dwt_prod100c_gain_col(c)*dt
       endif

       if ( use_c14 ) then
          cnveg_carbonstate_inst%prod10c_col(c)  = cnveg_carbonstate_inst%prod10c_col(c)  + &
               cnveg_carbonflux_inst%dwt_prod10c_gain_col(c) *dt
          cnveg_carbonstate_inst%prod100c_col(c) = cnveg_carbonstate_inst%prod100c_col(c) + &
               cnveg_carbonflux_inst%dwt_prod100c_gain_col(c)*dt
       endif

       cnveg_nitrogenstate_inst%prod10n_col(c)    = cnveg_nitrogenstate_inst%prod10n_col(c)    + &
            cnveg_nitrogenflux_inst%dwt_prod10n_gain_col(c)*dt
       cnveg_nitrogenstate_inst%prod100n_col(c)   = cnveg_nitrogenstate_inst%prod100n_col(c)   + &
            cnveg_nitrogenflux_inst%dwt_prod100n_gain_col(c)*dt

       ! fluxes into wood product pools, from harvest
       cnveg_carbonstate_inst%prod10c_col(c)    = cnveg_carbonstate_inst%prod10c_col(c)    + &
            cnveg_carbonflux_inst%hrv_deadstemc_to_prod10c_col(c)*dt
       cnveg_carbonstate_inst%prod100c_col(c)   = cnveg_carbonstate_inst%prod100c_col(c)   + &
            cnveg_carbonflux_inst%hrv_deadstemc_to_prod100c_col(c)*dt

       if ( use_c13 ) then
          cnveg_carbonstate_inst%prod10c_col(c)  = cnveg_carbonstate_inst%prod10c_col(c)  + &
               cnveg_carbonflux_inst%hrv_deadstemc_to_prod10c_col(c)*dt
          cnveg_carbonstate_inst%prod100c_col(c) = cnveg_carbonstate_inst%prod100c_col(c) + &
               cnveg_carbonflux_inst%hrv_deadstemc_to_prod100c_col(c)*dt
       endif

       if ( use_c14 ) then
          cnveg_carbonstate_inst%prod10c_col(c)  = cnveg_carbonstate_inst%prod10c_col(c)  + &
               cnveg_carbonflux_inst%hrv_deadstemc_to_prod10c_col(c)*dt
          cnveg_carbonstate_inst%prod100c_col(c) = cnveg_carbonstate_inst%prod100c_col(c) + &
               cnveg_carbonflux_inst%hrv_deadstemc_to_prod100c_col(c)*dt
       endif

       cnveg_nitrogenstate_inst%prod10n_col(c)    = cnveg_nitrogenstate_inst%prod10n_col(c)    + &
            cnveg_nitrogenflux_inst%hrv_deadstemn_to_prod10n_col(c)*dt
       cnveg_nitrogenstate_inst%prod100n_col(c)   = cnveg_nitrogenstate_inst%prod100n_col(c)   + &
            cnveg_nitrogenflux_inst%hrv_deadstemn_to_prod100n_col(c)*dt

       ! fluxes out of wood product pools, from decomposition
       cnveg_carbonstate_inst%prod10c_col(c)    = cnveg_carbonstate_inst%prod10c_col(c)    - &
            cnveg_carbonflux_inst%prod10c_loss_col(c)*dt
       cnveg_carbonstate_inst%prod100c_col(c)   = cnveg_carbonstate_inst%prod100c_col(c)   - &
            cnveg_carbonflux_inst%prod100c_loss_col(c)*dt

       if ( use_c13 ) then
          cnveg_carbonstate_inst%prod10c_col(c)  = cnveg_carbonstate_inst%prod10c_col(c)  - &
               cnveg_carbonflux_inst%prod10c_loss_col(c)*dt
          cnveg_carbonstate_inst%prod100c_col(c) = cnveg_carbonstate_inst%prod100c_col(c) - &
               cnveg_carbonflux_inst%prod100c_loss_col(c)*dt
       endif

       if ( use_c14 ) then
          cnveg_carbonstate_inst%prod10c_col(c)  = cnveg_carbonstate_inst%prod10c_col(c)  - &
               cnveg_carbonflux_inst%prod10c_loss_col(c)*dt
          cnveg_carbonstate_inst%prod100c_col(c) = cnveg_carbonstate_inst%prod100c_col(c) - &
               cnveg_carbonflux_inst%prod100c_loss_col(c)*dt
       endif

       cnveg_nitrogenstate_inst%prod10n_col(c)    = cnveg_nitrogenstate_inst%prod10n_col(c)    - &
            cnveg_nitrogenflux_inst%prod10n_loss_col(c)*dt
       cnveg_nitrogenstate_inst%prod100n_col(c)   = cnveg_nitrogenstate_inst%prod100n_col(c)   - &
            cnveg_nitrogenflux_inst%prod100n_loss_col(c)*dt

    end do ! end of column loop

  end subroutine CNWoodProducts

end module CNWoodProductsMod
