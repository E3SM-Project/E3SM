module EMI_CNCarbonStateType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod                           , only : emi_data_list, emi_data
  use EMI_DataDimensionMod                  , only : emi_data_dimension_list_type
  use CNCarbonStateType    , only : carbonstate_type
  use EMI_Atm2LndType_Constants
  use EMI_CanopyStateType_Constants
  use EMI_ChemStateType_Constants
  use EMI_CNCarbonStateType_Constants
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
  public :: EMI_Pack_CNCarbonStateType_at_Column_Level_for_EM
  public :: EMI_Unpack_CNCarbonStateType_at_Column_Level_from_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_CNCarbonStateType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter, carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM carbonstate_vars for EM
    !
    ! !USES:
    use elm_varpar             , only : nlevdecomp_full
    use elm_varpar             , only : ndecomp_pools
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    integer                , intent(in) :: num_filter
    integer                , intent(in) :: filter(:)
    type(carbonstate_type) , intent(in) :: carbonstate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j,k
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         decomp_cpools_vr => carbonstate_vars%decomp_cpools_vr_col   &
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

          case (L2E_STATE_CARBON_POOLS_VERTICALLY_RESOLVED)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      cur_data%data_real_3d(c,j,k) = decomp_cpools_vr(c,j,k)
                   enddo
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_CNCarbonStateType_at_Column_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Unpack_CNCarbonStateType_at_Column_Level_from_EM(data_list, em_stage, &
        num_filter, filter, carbonstate_vars)
    !
    ! !DESCRIPTION:
    ! Unpack data for ALM carbonstate_vars from EM
    !
    ! !USES:
    use elm_varpar             , only : nlevdecomp_full
    use elm_varpar             , only : ndecomp_pools
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(emi_data_list)   , intent(in) :: data_list
    integer                , intent(in) :: em_stage
    integer                , intent(in) :: num_filter
    integer                , intent(in) :: filter(:)
    type(carbonstate_type) , intent(in) :: carbonstate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j,k
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         decomp_cpools_vr => carbonstate_vars%decomp_cpools_vr_col   &
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

          case (E2L_STATE_CARBON_POOLS_VERTICALLY_RESOLVED)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      decomp_cpools_vr(c,j,k) = cur_data%data_real_3d(c,j,k)
                   enddo
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Unpack_CNCarbonStateType_at_Column_Level_from_EM


end module EMI_CNCarbonStateType_ExchangeMod
