module EMI_Landunit_Exchange
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod         , only : emi_data_list, emi_data
  use EMI_DataDimensionMod , only : emi_data_dimension_list_type
  use EMI_Landunit_Constants
  !
  implicit none
  !
  !
  public :: EMI_Pack_Landunit_for_EM

contains

  !-----------------------------------------------------------------------
  subroutine EMI_Pack_Landunit_for_EM(data_list, em_stage, &
        num_filter, filter)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's landunit type for EM
    !
    ! !USES:
    use LandunitType              , only : lun_pp
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_filter ! number of column soil points in column filter
    integer              , intent(in) :: filter(:)  ! column filter for soil points
    !
    integer                           :: fl,l
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

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

          case (L2E_LANDUNIT_TYPE)
             do fl = 1, num_filter
                l = filter(fl)
                cur_data%data_int_1d(l) = lun_pp%itype(l)
             enddo
             cur_data%is_set = .true.

          case (L2E_LANDUNIT_LAKEPOINT)
             do fl = 1, num_filter
                l = filter(fl)
                if (lun_pp%lakpoi(l)) cur_data%data_int_1d(l) = 1
             enddo
             cur_data%is_set = .true.

          case (L2E_LANDUNIT_URBANPOINT)
             do fl = 1, num_filter
                l = filter(fl)
                if (lun_pp%urbpoi(l)) cur_data%data_int_1d(l) = 1
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

  end subroutine EMI_Pack_Landunit_for_EM

end module EMI_Landunit_Exchange
