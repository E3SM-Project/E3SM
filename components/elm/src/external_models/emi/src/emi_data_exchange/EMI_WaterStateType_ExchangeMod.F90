module EMI_WaterStateType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod                           , only : emi_data_list, emi_data
  use EMI_DataDimensionMod                  , only : emi_data_dimension_list_type
  use WaterStateType                        , only : waterstate_type
  use ColumnDataType                        , only : col_ws
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
    use elm_varpar             , only : nlevsoi, nlevgrnd, nlevsno
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
         h2osoi_liq    => col_ws%h2osoi_liq    , &
         h2osoi_ice    => col_ws%h2osoi_ice    , &
         soilp         => col_ws%soilp         , &
         frac_h2osfc   => col_ws%frac_h2osfc   , &
         finundated    => col_ws%finundated    , &
         h2osoi_liqvol => col_ws%h2osoi_liqvol , &
         h2osoi_icevol => col_ws%h2osoi_icevol , &
         h2osoi_vol    => col_ws%h2osoi_vol    , &
         air_vol       => col_ws%air_vol       , &
         rho_vap       => waterstate_vars%rho_vap_col       , &
         rhvap_soi     => waterstate_vars%rhvap_soi_col     , &
         smp_l         => col_ws%smp_l         , &
         h2osno        => col_ws%h2osno        , &
         h2osfc        => col_ws%h2osfc        , &
         frac_sno_eff  => col_ws%frac_sno_eff    &
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
    use elm_varpar             , only : nlevsoi, nlevgrnd, nlevsno
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
         h2osoi_liq => col_ws%h2osoi_liq , &
         h2osoi_ice => col_ws%h2osoi_ice , &
         soilp      => col_ws%soilp        &
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
