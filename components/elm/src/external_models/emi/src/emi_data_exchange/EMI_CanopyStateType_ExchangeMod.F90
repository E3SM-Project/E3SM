module EMI_CanopyStateType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod         , only : emi_data_list, emi_data
  use EMI_DataDimensionMod , only : emi_data_dimension_list_type
  use CanopyStateType                       , only : canopystate_type
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
  public :: EMI_Pack_CanopyStateType_at_Column_Level_for_EM
  public :: EMI_Pack_CanopyStateType_at_Patch_Level_for_EM
  public :: EMI_Unpack_CanopyStateType_at_Patch_Level_from_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_CanopyStateType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter, canopystate_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM canopystate_vars for EM
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
    type(canopystate_type) , intent(in) :: canopystate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         altmax          => canopystate_vars%altmax_col          , &
         altmax_lastyear => canopystate_vars%altmax_lastyear_col   &
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

          case (L2E_STATE_ALTMAX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = altmax(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_ALTMAX_LASTYEAR)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = altmax_lastyear(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_CanopyStateType_at_Column_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Pack_CanopyStateType_at_Patch_Level_for_EM(data_list, em_stage, &
        num_filter, filter, canopystate_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM canopystate_vars for EM
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
    type(canopystate_type) , intent(in) :: canopystate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fp,p,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         lbl_rsc_h2o => canopystate_vars%lbl_rsc_h2o_patch , &
         elai        => canopystate_vars%elai_patch          &
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

          case (L2E_STATE_LBL_RSC_H2O)
             do fp = 1, num_filter
                p = filter(fp)
                cur_data%data_real_1d(p) = lbl_rsc_h2o(p)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_ELAI)
             do fp = 1, num_filter
                p = filter(fp)
                cur_data%data_real_1d(p) = elai(p)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_CanopyStateType_at_Patch_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Unpack_CanopyStateType_at_Patch_Level_from_EM(data_list, em_stage, &
        num_filter, filter, canopystate_vars)
    !
    ! !DESCRIPTION:
    ! Unpack data for ALM canopystate_vars from EM
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
    type(canopystate_type) , intent(in) :: canopystate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fp,p,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         fsun   => canopystate_vars%fsun_patch   , &
         laisun => canopystate_vars%laisun_patch , &
         laisha => canopystate_vars%laisha_patch   &
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

          case (E2L_STATE_FSUN)
             do fp = 1, num_filter
                p = filter(fp)
                fsun(p) = cur_data%data_real_1d(p)
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_LAISUN)
             do fp = 1, num_filter
                p = filter(fp)
                laisun(p) = cur_data%data_real_1d(p)
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_LAISHA)
             do fp = 1, num_filter
                p = filter(fp)
                laisha(p) = cur_data%data_real_1d(p)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Unpack_CanopyStateType_at_Patch_Level_from_EM


end module EMI_CanopyStateType_ExchangeMod
