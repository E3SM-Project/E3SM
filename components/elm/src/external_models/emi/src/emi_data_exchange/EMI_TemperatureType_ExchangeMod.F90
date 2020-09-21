module EMI_TemperatureType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use clm_varctl                            , only : iulog
  use EMI_DataMod                           , only : emi_data_list, emi_data
  use EMI_DataDimensionMod                  , only : emi_data_dimension_list_type
  use TemperatureType                       , only : temperature_type
  use ColumnDataType                        , only : col_es
  use VegetationDataType                    , only : veg_es
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
  public :: EMI_Pack_TemperatureType_at_Column_Level_for_EM
  public :: EMI_Pack_TemperatureType_at_Patch_Level_for_EM
  public :: EMI_Unpack_TemperatureType_at_Column_Level_from_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_TemperatureType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM temperature_vars for EM
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
    type(temperature_type) , intent(in) :: temperature_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         t_soisno  => col_es%t_soisno  , &
         t_h2osfc  => col_es%t_h2osfc  , &
         t_soi10cm => col_es%t_soi10cm   &
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

          case (L2E_STATE_TSOIL_NLEVGRND)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = t_soisno(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_TSNOW)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno + 1, 0
                   cur_data%data_real_2d(c,j) = t_soisno(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_TH2OSFC)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = t_h2osfc(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_TSOI10CM)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_real_1d(c) = t_soi10cm(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_STATE_TSOIL_NLEVSOI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevsoi
                   cur_data%data_real_2d(c,j) = t_soisno(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_TemperatureType_at_Column_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Pack_TemperatureType_at_Patch_Level_for_EM(data_list, em_stage, &
        num_filter, filter, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM temperature_vars for EM
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
    type(temperature_type) , intent(in) :: temperature_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fp,p,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         t_veg => veg_es%t_veg   &
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

          case (L2E_STATE_TVEG)
             do fp = 1, num_filter
                p = filter(fp)
                cur_data%data_real_1d(p) = t_veg(p)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_TemperatureType_at_Patch_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Unpack_TemperatureType_at_Column_Level_from_EM(data_list, em_stage, &
        num_filter, filter, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Unpack data for ALM temperature_vars from EM
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
    type(temperature_type) , intent(in) :: temperature_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         t_soisno => col_es%t_soisno , &
         t_h2osfc => col_es%t_h2osfc   &
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

          case (E2L_STATE_TSOIL_NLEVGRND)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   t_soisno(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_TSNOW_NLEVSNOW)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno + 1, 0
                   t_soisno(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (E2L_STATE_TH2OSFC)
             do fc = 1, num_filter
                c = filter(fc)
                t_h2osfc(c) = cur_data%data_real_1d(c)
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Unpack_TemperatureType_at_Column_Level_from_EM


end module EMI_TemperatureType_ExchangeMod
