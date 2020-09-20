module ExternalModelBETRMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module provides wrapper for BeTR in ALM
  !
  use abortutils                   , only : endrun
  use shr_kind_mod                 , only : r8 => shr_kind_r8
  use shr_log_mod                  , only : errMsg => shr_log_errMsg
  use EMI_DataMod, only : emi_data_list, emi_data
  use decompMod                    , only : bounds_type
  use ExternalModelConstants
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ColumnType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_Filter_Constants
  use EMI_Landunit_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  !
  implicit none
  !

  integer :: index_l2e_state_frac_h2osfc
  integer :: index_l2e_state_finundated
  integer :: index_l2e_state_h2osoi_liq
  integer :: index_l2e_state_h2osoi_ice
  integer :: index_l2e_state_h2osoi_liqvol
  integer :: index_l2e_state_h2osoi_icevol
  integer :: index_l2e_state_h2osoi_vol
  integer :: index_l2e_state_air_vol
  integer :: index_l2e_state_rho_vap
  integer :: index_l2e_state_rhvap_soi
  integer :: index_l2e_state_smp_l

  integer :: index_l2e_filter_nolakec
  integer :: index_l2e_filter_num_nolakec

  public :: EM_BETR_Populate_L2E_List, &
            EM_BETR_Populate_E2L_List, &
            EM_BETR_Solve

contains

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_L2E_List(l2e_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables needed by BeTR from ALM
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list), intent(inout) :: l2e_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages
    integer                  :: id
    integer                  :: index

    number_em_stages = 1
    allocate(em_stages(number_em_stages))

    em_stages(1) = EM_BETR_PRE_DIAG_WATER_FLUX_STAGE

    id                            = L2E_STATE_FRAC_H2OSFC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_frac_h2osfc   = index

    id                            = L2E_STATE_FRAC_INUNDATED
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_finundated    = index

    id                            = L2E_STATE_H2OSOI_LIQ_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_h2osoi_liq    = index

    id                            = L2E_STATE_H2OSOI_ICE_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_h2osoi_ice    = index

    id                            = L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_h2osoi_liqvol = index

    id                            = L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_h2osoi_icevol = index

    id                            = L2E_STATE_H2OSOI_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_h2osoi_vol    = index

    id                            = L2E_STATE_AIR_VOL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_air_vol       = index

    id                            = L2E_STATE_RHO_VAP_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_rho_vap       = index

    id                            = L2E_STATE_RHVAP_SOI_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_rhvap_soi     = index

    id                            = L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_state_smp_l         = index

    id                            = L2E_FILTER_NOLAKEC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_filter_nolakec      = index

    id                            = L2E_FILTER_NUM_NOLAKEC
    call l2e_list%AddDataByID(id, number_em_stages, em_stages, index)
    index_l2e_filter_num_nolakec  = index

    deallocate(em_stages)

  end subroutine EM_BETR_Populate_L2E_List

  !------------------------------------------------------------------------
  subroutine EM_BETR_Populate_E2L_List(e2l_list)
    !
    ! !DESCRIPTION:
    ! Create a list of all variables returned by BeTR to ALM
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list), intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    class(emi_data), pointer :: data
    integer        , pointer :: em_stages(:)
    integer                  :: number_em_stages


  end subroutine EM_BETR_Populate_E2L_List

    !------------------------------------------------------------------------
  subroutine EM_BETR_Solve(em_stage, dt, nstep, bounds, l2e_list, e2l_list, &
       bounds_clump)
    !
    ! !DESCRIPTION:
    ! 
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use clm_varctl                , only : iulog
    !
    implicit none
    !
    ! !ARGUMENTS:
    integer              , intent(in)    :: em_stage
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    type(bounds_type)    , intent(in)    :: bounds
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    type(bounds_type)    , intent (in)   :: bounds_clump

    select case(em_stage)

    case (EM_BETR_BEGIN_MASS_BALANCE_STAGE)
       call EM_BETR_BeginMassBalance_Solve(dt, nstep, bounds, l2e_list, e2l_list)

    case (EM_BETR_PRE_DIAG_WATER_FLUX_STAGE)
       call EM_BETR_PreDiagSoilColWaterFlux_Solve(bounds, l2e_list, e2l_list)

    case default
       write(iulog,*)'EM_BETR_Solve: Unknown em_stage.'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine EM_BETR_Solve

    !------------------------------------------------------------------------
  subroutine EM_BETR_BeginMassBalance_Solve(dt, nstep, bounds_clump, l2e_list, e2l_list)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use clm_varctl                , only : iulog
#ifdef BETR_VIA_EMI
    use elm_instMod               , only : ep_betr
#endif
    !
    implicit none
    !
    ! !ARGUMENTS:
    real(r8)             , intent(in)    :: dt
    integer              , intent(in)    :: nstep
    type(bounds_type)    , intent(in)    :: bounds_clump
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list

#ifdef BETR_VIA_EMI
    call ep_betr%SetClock(dtime = dt, nelapstep = nstep)

    call ep_betr%BeginMassBalanceCheck(bounds_clump)
#else
    call endrun('BeTR is on but code was not compiled with -DBETR_VIA_EMI')
#endif


  end subroutine EM_BETR_BeginMassBalance_Solve

    !------------------------------------------------------------------------
  subroutine EM_BETR_PreDiagSoilColWaterFlux_Solve(bounds, l2e_list, e2l_list)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use abortutils                , only : endrun
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use clm_varctl                , only : iulog
    use clm_varpar                , only : nlevsoi
#ifdef BETR_VIA_EMI
    use elm_instMod               , only : ep_betr
    use BeTR_decompMod            , only : betr_bounds_type
#endif
    !
    implicit none
    !
    ! !ARGUMENTS:
    type(bounds_type)    , intent(in)    :: bounds
    class(emi_data_list) , intent(in)    :: l2e_list
    class(emi_data_list) , intent(inout) :: e2l_list
    !
    ! !LOCAL VARIABLES:
    real(r8), pointer      :: l2e_frac_h2osfc(:)
    real(r8), pointer      :: l2e_finundated(:)
    real(r8), pointer      :: l2e_h2osoi_liq(:,:)
    real(r8), pointer      :: l2e_h2osoi_ice(:,:)
    real(r8), pointer      :: l2e_h2osoi_liqvol(:,:)
    real(r8), pointer      :: l2e_h2osoi_icevol(:,:)
    real(r8), pointer      :: l2e_h2osoi_vol(:,:)
    real(r8), pointer      :: l2e_air_vol(:,:)
    real(r8), pointer      :: l2e_rho_vap(:,:)
    real(r8), pointer      :: l2e_rhvap_soi(:,:)
    real(r8), pointer      :: l2e_smp_l(:,:)
    integer, pointer       :: l2e_filter_nolakec(:)
    integer                :: l2e_num_nolakec
    integer                :: cc, c, fc, lbj, ubj
#ifdef BETR_VIA_EMI
    type(betr_bounds_type) :: betr_bounds

    call l2e_list%GetPointerToReal1D(index_l2e_state_frac_h2osfc  , l2e_frac_h2osfc    )
    call l2e_list%GetPointerToReal1D(index_l2e_state_finundated   , l2e_finundated     )
    call l2e_list%GetPointerToReal2D(index_l2e_state_h2osoi_liq   , l2e_h2osoi_liq     )
    call l2e_list%GetPointerToReal2D(index_l2e_state_h2osoi_ice   , l2e_h2osoi_ice     )
    call l2e_list%GetPointerToReal2D(index_l2e_state_h2osoi_liqvol, l2e_h2osoi_liqvol  )
    call l2e_list%GetPointerToReal2D(index_l2e_state_h2osoi_icevol, l2e_h2osoi_icevol  )
    call l2e_list%GetPointerToReal2D(index_l2e_state_h2osoi_vol   , l2e_h2osoi_vol     )
    call l2e_list%GetPointerToReal2D(index_l2e_state_air_vol      , l2e_air_vol        )
    call l2e_list%GetPointerToReal2D(index_l2e_state_rho_vap      , l2e_rho_vap        )
    call l2e_list%GetPointerToReal2D(index_l2e_state_rhvap_soi    , l2e_rhvap_soi      )
    call l2e_list%GetPointerToReal2D(index_l2e_state_smp_l        , l2e_smp_l          )

    call l2e_list%GetPointerToInt1D(index_l2e_filter_nolakec      , l2e_filter_nolakec )
    call l2e_list%GetIntValue(index_l2e_filter_num_nolakec        , l2e_num_nolakec    )

    cc  = 1
    lbj = 1
    ubj = nlevsoi

    do c = bounds%begc, bounds%endc
      if (.not. ep_betr%active_col(c)) cycle

      !assign waterstate
      ep_betr%biophys_forc(c)%finundated_col    (cc)         = l2e_finundated   (c)
      ep_betr%biophys_forc(c)%frac_h2osfc_col   (cc)         = l2e_frac_h2osfc  (c)
      ep_betr%biophys_forc(c)%h2osoi_liq_col    (cc,lbj:ubj) = l2e_h2osoi_liq   (c,lbj:ubj)
      ep_betr%biophys_forc(c)%h2osoi_ice_col    (cc,lbj:ubj) = l2e_h2osoi_ice   (c,lbj:ubj)
      ep_betr%biophys_forc(c)%h2osoi_liqvol_col (cc,lbj:ubj) = l2e_h2osoi_liqvol(c,lbj:ubj)
      ep_betr%biophys_forc(c)%h2osoi_icevol_col (cc,lbj:ubj) = l2e_h2osoi_icevol(c,lbj:ubj)
      ep_betr%biophys_forc(c)%h2osoi_vol_col    (cc,lbj:ubj) = l2e_h2osoi_vol   (c,lbj:ubj)
      ep_betr%biophys_forc(c)%air_vol_col       (cc,lbj:ubj) = l2e_air_vol      (c,lbj:ubj)
      ep_betr%biophys_forc(c)%rho_vap           (cc,lbj:ubj) = l2e_rho_vap      (c,lbj:ubj)
      ep_betr%biophys_forc(c)%rhvap_soi         (cc,lbj:ubj) = l2e_rhvap_soi    (c,lbj:ubj)
      ep_betr%biophys_forc(c)%smp_l_col         (cc,lbj:ubj) = l2e_smp_l        (c,lbj:ubj)

    enddo

    call ep_betr%BeTRSetBounds(betr_bounds)

    do fc = 1, l2e_num_nolakec
       c = l2e_filter_nolakec(fc)
       if (.not. ep_betr%active_col(c)) cycle

       call ep_betr%betr(c)%pre_diagnose_soilcol_water_flux(betr_bounds, ep_betr%num_soilc, &
            ep_betr%filter_soilc, ep_betr%biophys_forc(c))
   enddo

#else
    call endrun('BeTR is on but code was not compiled with -DBETR_VIA_EMI')
#endif

  end subroutine EM_BETR_PreDiagSoilColWaterFlux_Solve

end module ExternalModelBETRMod
