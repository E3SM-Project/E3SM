module EMI_Filter_Exchange
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod         , only : emi_data_list, emi_data
  use EMI_DataDimensionMod , only : emi_data_dimension_list_type
  use EMI_Filter_Constants
  !
  implicit none
  !
  !
  public :: EMI_Pack_Filter_for_EM

contains

  !-----------------------------------------------------------------------
  subroutine EMI_Pack_Filter_for_EM(data_list, em_stage, &
        num_filter, filter)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's filter for EM
    !
    ! !USES:
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_filter ! number of column soil points in column filter
    integer              , intent(in) :: filter(:)  ! column filter for soil points
    !
    integer                           :: i
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage

    cur_data => data_list%first
    do
       if (.not.associated(cur_data)) exit

       need_to_pack = .false.
       do istage = 1, cur_data%num_em_stages
          if (cur_data%em_stage_ids(istage) == em_stage) then
             need_to_pack = .true.
             exit
          endif
       enddo

       if (need_to_pack) then

          select case (cur_data%id)

          case (L2E_FILTER_HYDROLOGYC, L2E_FILTER_NOLAKEC, L2E_FILTER_NOLAKEC_AND_NOURBANC)
             do i = 1, num_filter
                cur_data%data_int_1d(i) = filter(i)
             enddo
             cur_data%is_set = .true.

          case (L2E_FILTER_NUM_HYDROLOGYC, L2E_FILTER_NUM_NOLAKEC, L2E_FILTER_NUM_NOLAKEC_AND_NOURBANC)

             cur_data%data_int_1d(1) = num_filter
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

  end subroutine EMI_Pack_Filter_for_EM

end module EMI_Filter_Exchange
