module EMI_Atm2LndType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use clm_varctl                            , only : iulog
  use EMI_DataMod                           , only : emi_data_list, emi_data
  use EMI_DataDimensionMod                  , only : emi_data_dimension_list_type
  use Atm2LndType                           , only : atm2lnd_type
  use TopounitDataType                      , only : top_as
  use ColumnType                            , only : col_pp
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
  public :: EMI_Pack_Atm2LndType_at_Column_Level_for_EM
  public :: EMI_Pack_Atm2LndType_at_Grid_Level_for_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_Atm2LndType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter, atm2lndtype_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM atm2lndtype_vars for EM
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
    type(atm2lnd_type)     , intent(in) :: atm2lndtype_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,t,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         forc_pbot_downscaled => top_as%pbot , &
         forc_t_downscaled    => top_as%tbot   &
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

          case (L2E_STATE_FORC_PBOT_DOWNSCALED)
             do fc = 1, num_filter
                c = filter(fc)
                t = col_pp%topounit(c)
                cur_data%data_real_1d(c) = forc_pbot_downscaled(t)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_FORC_T_DOWNSCALED)
             do fc = 1, num_filter
                c = filter(fc)
                t = col_pp%topounit(c)
                cur_data%data_real_1d(c) = forc_t_downscaled(t)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_Atm2LndType_at_Column_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Pack_Atm2LndType_at_Grid_Level_for_EM(data_list, em_stage, &
        num_filter, filter, atm2lndtype_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM atm2lndtype_vars for EM
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
    type(atm2lnd_type)     , intent(in) :: atm2lndtype_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fg,g,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         forc_solad => atm2lndtype_vars%forc_solad_grc , &
         forc_solai => atm2lndtype_vars%forc_solai_grc   &
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

          case (L2E_FLUX_SOLAR_DIRECT_RADDIATION)
             do fg = 1, num_filter
                g = filter(fg)
                do j = 1, 2
                   cur_data%data_real_2d(g,j) = forc_solad(g,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_FLUX_SOLAR_DIFFUSE_RADDIATION)
             do fg = 1, num_filter
                g = filter(fg)
                do j = 1, 2
                   cur_data%data_real_2d(g,j) = forc_solai(g,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_Atm2LndType_at_Grid_Level_for_EM


end module EMI_Atm2LndType_ExchangeMod
