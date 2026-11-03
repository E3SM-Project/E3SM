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
  use elm_time_manager    , only : get_step_size
  use elm_varctl          , only : use_c13, use_c14

  use ColumnDataType      , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType      , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType      , only : col_ns, col_nf, col_ps, col_pf
  use timeInfoMod
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
  subroutine WoodProducts(num_soilc, filter_soilc )
    !
    ! !DESCRIPTION:
    ! Update all loss fluxes from wood product pools, and update product pool state variables
    ! for both loss and gain terms.  Gain terms are calculated in pftdyn_cnbal() for gains associated
    ! with changes in landcover, and in CNHarvest(), for gains associated with wood harvest.
    !
    ! !ARGUMENTS:
    integer                    , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                    , intent(in)    :: filter_soilc(:) ! filter for soil columns
    real(r8)   :: dt       ! time step (seconds)


    !
    ! !LOCAL VARIABLES:
    integer  :: fc       ! lake filter indices
    integer  :: c        ! indices
    real(r8) :: kprod10  ! decay constant for 10-year product pool
    real(r8) :: kprod100 ! decay constant for 100-year product pool
    !-----------------------------------------------------------------------
    dt = dtime_mod
    ! calculate column-level losses from product pools
    ! the following (1/s) rate constants result in ~90% loss of initial state over 10 and 100 years,
    ! respectively, using a discrete-time fractional decay algorithm.
    kprod10 = 7.2e-9
    kprod100 = 7.2e-10

    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! calculate fluxes (1/sec)
       col_cf%core%prod10_loss(c)    = col_cs%core%prod10(c)    * kprod10
       col_cf%core%prod100_loss(c)   = col_cs%core%prod100(c)   * kprod100

       if ( use_c13 ) then
          c13_col_cf%core%prod10_loss(c)  = c13_col_cs%core%prod10(c)  * kprod10
          c13_col_cf%core%prod100_loss(c) = c13_col_cs%core%prod100(c) * kprod100
       endif

       if ( use_c14 ) then
          c14_col_cf%core%prod10_loss(c)  = c14_col_cs%core%prod10(c)  * kprod10
          c14_col_cf%core%prod100_loss(c) = c14_col_cs%core%prod100(c) * kprod100
       endif

       col_nf%core%prod10_loss(c)    = col_ns%core%prod10(c)    * kprod10
       col_nf%core%prod100_loss(c)   = col_ns%core%prod100(c)   * kprod100

       col_pf%core%prod10_loss(c)    = col_ps%core%prod10(c)    * kprod10
       col_pf%core%prod100_loss(c)   = col_ps%core%prod100(c)   * kprod100
    end do


    ! update wood product state variables
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       col_cs%core%prod10(c)    = col_cs%core%prod10(c)             &
            + col_cf%core%hrv_deadstem_to_prod10(c)*dt         & ! from harvest
            - col_cf%core%prod10_loss(c)*dt                       ! from decomposition

       col_cs%core%prod100(c)   = col_cs%core%prod100(c)            &
            + col_cf%core%hrv_deadstem_to_prod100(c)*dt        & ! from harvest
            - col_cf%core%prod100_loss(c)*dt                      ! from decomposition


       if ( use_c13 ) then
          c13_col_cs%core%prod10(c)  = c13_col_cs%core%prod10(c)    &
               + c13_col_cf%core%hrv_deadstem_to_prod10(c)*dt  & ! from harvest
               - c13_col_cf%core%prod10_loss(c)*dt                ! from decomposition

          c13_col_cs%core%prod100(c) = c13_col_cs%core%prod100(c)   &
               + c13_col_cf%core%hrv_deadstem_to_prod100(c)*dt & ! from harvest
               - c13_col_cf%core%prod100_loss(c)*dt               ! from decomposition
       endif

       if ( use_c14 ) then
          c14_col_cs%core%prod10(c)  = c14_col_cs%core%prod10(c)    &
               + c14_col_cf%core%hrv_deadstem_to_prod10(c)*dt  & ! from harvest
               - c14_col_cf%core%prod10_loss(c)*dt                ! from decomposition

          c14_col_cs%core%prod100(c) = c14_col_cs%core%prod100(c)   &
               + c14_col_cf%core%hrv_deadstem_to_prod100(c)*dt & ! from harvest
               - c14_col_cf%core%prod100_loss(c)*dt               ! from decomposition
       endif

       col_ns%core%prod10(c)    = col_ns%core%prod10(c)             &
            + col_nf%core%hrv_deadstem_to_prod10(c)*dt         & ! from harvest
            - col_nf%core%prod10_loss(c)*dt                       ! from decomposition

       col_ns%core%prod100(c)   = col_ns%core%prod100(c)            &
            + col_nf%core%hrv_deadstem_to_prod100(c)*dt        & ! from harvest
            - col_nf%core%prod100_loss(c)*dt                      ! from decomposition

       col_ps%core%prod10(c)    = col_ps%core%prod10(c)             &
            + col_pf%core%hrv_deadstem_to_prod10(c)*dt         & ! from harvest
            - col_pf%core%prod10_loss(c)*dt                       ! from decomposition

       col_ps%core%prod100(c)   = col_ps%core%prod100(c)            &
            + col_pf%core%hrv_deadstem_to_prod100(c)*dt        & ! from harvest
            - col_pf%core%prod100_loss(c)*dt                      ! from decomposition

    end do ! end of column loop

  end subroutine WoodProducts

end module WoodProductsMod
