module CNCropHarvestPoolsMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate loss fluxes from crop harvest pools, and update product pool state variables
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
  use PhosphorusStateType , only : phosphorusstate_type
  use PhosphorusFluxType  , only : phosphorusflux_type
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: CNCropHarvestPools
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNCropHarvestPools(num_soilc, filter_soilc, &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, phosphorusstate_vars,&
       carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, nitrogenflux_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! Update all loss fluxes from crop harvest pools, and update harvest pool state variables
    ! for both loss and gain terms. GAin terms are calculated in CNCropHarvestPools() for gains
    ! associated with crop harvest. 
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(carbonstate_type)   , intent(in)    :: c13_carbonstate_vars
    type(carbonstate_type)   , intent(in)    :: c14_carbonstate_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(phosphorusstate_type), intent(in)   :: phosphorusstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type), intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: fc        ! lake filter indices
    integer :: c         ! indices
    real(r8):: dt        ! time step (seconds)
    real(r8) :: kprod1       ! decay constant for 1-year product pool
    !-----------------------------------------------------------------------


    ! calculate column-level losses from product pools
    ! the following (1/s) rate constants result in ~90% loss of initial state over 1 year,
    ! using a discrete-time fractional decay algorithm.
    kprod1 = 7.2e-9

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! calculate fluxes (1/sec)
       carbonflux_vars%prod1c_loss_col(c)    = carbonstate_vars%prod1c_col(c)    * kprod1

       if ( use_c13 ) then
          c13_carbonflux_vars%prod1c_loss_col(c)  = c13_carbonstate_vars%prod1c_col(c)  * kprod1
       endif

       if ( use_c14 ) then
          c14_carbonflux_vars%prod1c_loss_col(c)  = c14_carbonstate_vars%prod1c_col(c)  * kprod1
       endif

       nitrogenflux_vars%prod1n_loss_col(c)    = nitrogenstate_vars%prod1n_col(c)    * kprod1
       phosphorusflux_vars%prod1p_loss_col(c)  = phosphorusstate_vars%prod1p_col(c)  * kprod1
    end do

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! update wood product state variables
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! fluxes into wood product pools, from harvest
       carbonstate_vars%prod1c_col(c)    = carbonstate_vars%prod1c_col(c)    + &
            carbonflux_vars%hrv_cropc_to_prod1c_col(c)*dt

       if ( use_c13 ) then
          c13_carbonstate_vars%prod1c_col(c)  = c13_carbonstate_vars%prod1c_col(c)  + &
               c13_carbonflux_vars%hrv_cropc_to_prod1c_col(c)*dt
       endif

       if ( use_c14 ) then
          c14_carbonstate_vars%prod1c_col(c)  = c14_carbonstate_vars%prod1c_col(c)  + &
               c14_carbonflux_vars%hrv_cropc_to_prod1c_col(c)*dt
       endif

       nitrogenstate_vars%prod1n_col(c)    = nitrogenstate_vars%prod1n_col(c)    + &
            nitrogenflux_vars%hrv_cropn_to_prod1n_col(c)*dt
       phosphorusstate_vars%prod1p_col(c)    = phosphorusstate_vars%prod1p_col(c)    + &
            phosphorusflux_vars%hrv_cropp_to_prod1p_col(c)*dt

       ! fluxes out of wood product pools, from decomposition
       carbonstate_vars%prod1c_col(c)    = carbonstate_vars%prod1c_col(c)    - &
            carbonflux_vars%prod1c_loss_col(c)*dt

       if ( use_c13 ) then
          c13_carbonstate_vars%prod1c_col(c)  = c13_carbonstate_vars%prod1c_col(c)  - &
               c13_carbonflux_vars%prod1c_loss_col(c)*dt
       endif

       if ( use_c14 ) then
          c14_carbonstate_vars%prod1c_col(c)  = c14_carbonstate_vars%prod1c_col(c)  - &
               c14_carbonflux_vars%prod1c_loss_col(c)*dt
       endif

       nitrogenstate_vars%prod1n_col(c)    = nitrogenstate_vars%prod1n_col(c)    - &
            nitrogenflux_vars%prod1n_loss_col(c)*dt
       phosphorusstate_vars%prod1p_col(c)  = phosphorusstate_vars%prod1p_col(c)  - &
            phosphorusflux_vars%prod1p_loss_col(c)*dt

    end do ! end of column loop

  end subroutine CNCropHarvestPools

end module CNCropHarvestPoolsMod
