module EMI_WaterStateType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use clm_varctl                            , only : iulog
  use ExternalModelInterfaceDataMod         , only : emi_data_list, emi_data
  use ExternalModelIntefaceDataDimensionMod , only : emi_data_dimension_list_type
  use WaterStateType                        , only : waterstate_type
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ChemStateType_Constants
  use EMI_EnergyFluxType_Constants
  use EMI_SoilHydrologyType_Constants
  use EMI_SoilStateType_Constants
  use EMI_TemperatureType_Constants
  use EMI_WaterFluxType_Constants
  use EMI_WaterStateType_Constants
  use EMI_Filter_Constants
  use EMI_ColumnType_Constants
  use EMI_Landunit_Constants
  !
  implicit none
  !
  !
  public :: EMI_Pack_WaterStateType_at_Column_Level_for_EM
  public :: EMI_Unpack_WaterStateType_at_Column_Level_from_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_WaterStateType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM waterstate_vars for EM
    !
    ! !USES:
    use clm_varpar             , only : nlevsoi, nlevgrnd, nlevsno
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    integer                , intent(in) :: num_filter
    integer                , intent(in) :: filter(:)
    type(waterstate_type)  , intent(in) :: waterstate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         h2osoi_liq    => waterstate_vars%h2osoi_liq_col    , &
         h2osoi_ice    => waterstate_vars%h2osoi_ice_col    , &
         soilp         => waterstate_vars%soilp_col         , &
         frac_h2osfc   => waterstate_vars%frac_h2osfc_col   , &
         finundated    => waterstate_vars%finundated_col    , &
         h2osoi_liqvol => waterstate_vars%h2osoi_liqvol_col , &
         h2osoi_icevol => waterstate_vars%h2osoi_icevol_col , &
         h2osoi_vol    => waterstate_vars%h2osoi_vol_col    , &
         air_vol       => waterstate_vars%air_vol_col       , &
         rho_vap       => waterstate_vars%rho_vap_col       , &
         rhvap_soi     => waterstate_vars%rhvap_soi_col     , &
         smp_l         => waterstate_vars%smp_l_col         , &
         h2osno        => waterstate_vars%h2osno_col        , &
         h2osfc        => waterstate_vars%h2osfc_col        , &
         frac_sno_eff  => waterstate_vars%frac_sno_eff_col    &
         )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_STATE_H2OSOI_LIQ_NLEVGRND)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = h2osoi_liq(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_ICE_NLEVGRND)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = h2osoi_ice(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_VSFM_PROGNOSTIC_SOILP)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = soilp(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_FRAC_H2OSFC)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = frac_h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_FRAC_INUNDATED)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = finundated(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_liqvol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_icevol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_VOL_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_vol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_AIR_VOL_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = air_vol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_RHO_VAP_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = rho_vap(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_RHVAP_SOI_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = rhvap_soi(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = smp_l(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_LIQ_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_liq(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_ICE_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = h2osoi_ice(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_LIQ_NLEVSNOW)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno + 1, 0
                   cur_data%data_real_2d(c,j) = h2osoi_liq(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSOI_ICE_NLEVSNOW)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno + 1, 0
                   cur_data%data_real_2d(c,j) = h2osoi_ice(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSNOW)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = h2osno(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_H2OSFC)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_FRAC_SNOW_EFFECTIVE)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = frac_sno_eff(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_WaterStateType_at_Column_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Unpack_WaterStateType_at_Column_Level_from_EM(data_list, em_stage, &
        num_filter, filter, waterstate_vars)
    !
    ! !DESCRIPTION:
    ! Unpack data for ALM waterstate_vars from EM
    !
    ! !USES:
    use clm_varpar             , only : nlevsoi, nlevgrnd, nlevsno
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    integer                , intent(in) :: num_filter
    integer                , intent(in) :: filter(:)
    type(waterstate_type)  , intent(in) :: waterstate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         h2osoi_liq => waterstate_vars%h2osoi_liq_col , &
         h2osoi_ice => waterstate_vars%h2osoi_ice_col , &
         soilp      => waterstate_vars%soilp_col        &
         )

    count = 0
    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit
       count = count + 1

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (E2L_STATE_H2OSOI_LIQ)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   h2osoi_liq(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_H2OSOI_ICE)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   h2osoi_ice(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_VSFM_PROGNOSTIC_SOILP)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   soilp(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Unpack_WaterStateType_at_Column_Level_from_EM


end module EMI_WaterStateType_ExchangeMod
