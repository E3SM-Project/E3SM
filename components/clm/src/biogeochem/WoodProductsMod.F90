module WoodProductsMod
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

  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type

  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: WoodProducts
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine WoodProducts(num_soilc, filter_soilc, &
       cs, c13_cs, c14_cs, ns, &
       cf, c13_cf, c14_cf, nf,&
       ps, pf)
    !
    ! !DESCRIPTION:
    ! Update all loss fluxes from wood product pools, and update product pool state variables
    ! for both loss and gain terms.  Gain terms are calculated in pftdyn_cnbal() for gains associated
    ! with changes in landcover, and in CNHarvest(), for gains associated with wood harvest.
    !
    ! !ARGUMENTS:
    integer                    , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                    , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(carbonstate_type)     , intent(inout) :: cs
    type(carbonstate_type)     , intent(inout) :: c13_cs
    type(carbonstate_type)     , intent(inout) :: c14_cs
    type(nitrogenstate_type)   , intent(inout) :: ns
    type(carbonflux_type)      , intent(inout) :: cf
    type(carbonflux_type)      , intent(inout) :: c13_cf
    type(carbonflux_type)      , intent(inout) :: c14_cf
    type(nitrogenflux_type)    , intent(inout) :: nf

    type(phosphorusstate_type) , intent(in)    :: ps
    type(phosphorusflux_type)  , intent(inout) :: pf

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
       cf%prod10c_loss_col(c)    = cs%prod10c_col(c)    * kprod10
       cf%prod100c_loss_col(c)   = cs%prod100c_col(c)   * kprod100

       if ( use_c13 ) then
          c13_cf%prod10c_loss_col(c)  = c13_cs%prod10c_col(c)  * kprod10
          c13_cf%prod100c_loss_col(c) = c13_cs%prod100c_col(c) * kprod100
       endif

       if ( use_c14 ) then
          c14_cf%prod10c_loss_col(c)  = c14_cs%prod10c_col(c)  * kprod10
          c14_cf%prod100c_loss_col(c) = c14_cs%prod100c_col(c) * kprod100
       endif

       nf%prod10n_loss_col(c)    = ns%prod10n_col(c)    * kprod10
       nf%prod100n_loss_col(c)   = ns%prod100n_col(c)   * kprod100

       pf%prod10p_loss_col(c)    = ps%prod10p_col(c)    * kprod10
       pf%prod100p_loss_col(c)   = ps%prod100p_col(c)   * kprod100
    end do

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! update wood product state variables
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       cs%prod10c_col(c)    = cs%prod10c_col(c)             &
            + cf%dwt_prod10c_gain_col(c)*dt                 & ! from landcover change
            + cf%hrv_deadstemc_to_prod10c_col(c)*dt         & ! from harvest
            - cf%prod10c_loss_col(c)*dt                       ! from decomposition

       cs%prod100c_col(c)   = cs%prod100c_col(c)            &
            + cf%dwt_prod100c_gain_col(c)*dt                & ! from landcover change
            + cf%hrv_deadstemc_to_prod100c_col(c)*dt        & ! from harvest
            - cf%prod100c_loss_col(c)*dt                      ! from decomposition


       if ( use_c13 ) then
          c13_cs%prod10c_col(c)  = c13_cs%prod10c_col(c)    &
               + c13_cf%dwt_prod10c_gain_col(c) *dt         & ! from landcover change
               + c13_cf%hrv_deadstemc_to_prod10c_col(c)*dt  & ! from harvest
               - c13_cf%prod10c_loss_col(c)*dt                ! from decomposition

          c13_cs%prod100c_col(c) = c13_cs%prod100c_col(c)   &
               + c13_cf%dwt_prod100c_gain_col(c)*dt         & ! from landcover change
               + c13_cf%hrv_deadstemc_to_prod100c_col(c)*dt & ! from harvest
               - c13_cf%prod100c_loss_col(c)*dt               ! from decomposition
       endif

       if ( use_c14 ) then
          c14_cs%prod10c_col(c)  = c14_cs%prod10c_col(c)    &
               + c14_cf%dwt_prod10c_gain_col(c) *dt         & ! from landcover change
               + c14_cf%hrv_deadstemc_to_prod10c_col(c)*dt  & ! from harvest
               - c14_cf%prod10c_loss_col(c)*dt                ! from decomposition

          c14_cs%prod100c_col(c) = c14_cs%prod100c_col(c)   &
               + c14_cf%dwt_prod100c_gain_col(c)*dt         & ! from landcover change
               + c14_cf%hrv_deadstemc_to_prod100c_col(c)*dt & ! from harvest
               - c14_cf%prod100c_loss_col(c)*dt               ! from decomposition
       endif

       ns%prod10n_col(c)    = ns%prod10n_col(c)             &
            + nf%dwt_prod10n_gain_col(c)*dt                 & ! from landcover change
            + nf%hrv_deadstemn_to_prod10n_col(c)*dt         & ! from harvest
            - nf%prod10n_loss_col(c)*dt                       ! from decomposition

       ns%prod100n_col(c)   = ns%prod100n_col(c)            &
            + nf%dwt_prod100n_gain_col(c)*dt                & ! from landcover change
            + nf%hrv_deadstemn_to_prod100n_col(c)*dt        & ! from harvest
            - nf%prod100n_loss_col(c)*dt                      ! from decomposition

       ps%prod10p_col(c)    = ps%prod10p_col(c)             &
            + pf%dwt_prod10p_gain_col(c)*dt                 & ! from landcover change
            + pf%hrv_deadstemp_to_prod10p_col(c)*dt         & ! from harvest
            - pf%prod10p_loss_col(c)*dt                       ! from decomposition

       ps%prod100p_col(c)   = ps%prod100p_col(c)            &
            + pf%dwt_prod100p_gain_col(c)*dt                & ! from landcover change
            + pf%hrv_deadstemp_to_prod100p_col(c)*dt        & ! from harvest
            - pf%prod100p_loss_col(c)*dt                      ! from decomposition

    end do ! end of column loop

  end subroutine WoodProducts

end module WoodProductsMod
