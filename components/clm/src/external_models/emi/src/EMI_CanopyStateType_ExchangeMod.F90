module EMI_CanopyStateType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use clm_varctl                            , only : iulog
  use ExternalModelInterfaceDataMod         , only : emi_data_list, emi_data
  use ExternalModelIntefaceDataDimensionMod , only : emi_data_dimension_list_type
  use CanopyStateType                       , only : canopystate_type
  !
  implicit none
  !
  !
  public :: EMI_Unpack_CanopyStateType_at_Patch_Level_from_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Unpack_CanopyStateType_at_Patch_Level_from_EM(data_list, em_stage, &
        num_filter, filter, canopystate_vars)
    !
    ! !DESCRIPTION:
    ! Unack data for ALM canopystate_vars from EM
    !
    ! !USES:
    use ExternalModelConstants , only : E2L_STATE_FSUN
    use ExternalModelConstants , only : E2L_STATE_LAISUN
    use ExternalModelConstants , only : E2L_STATE_LAISHA
    use clm_varpar             , only : nlevsoi, nlevgrnd, nlevsno
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
