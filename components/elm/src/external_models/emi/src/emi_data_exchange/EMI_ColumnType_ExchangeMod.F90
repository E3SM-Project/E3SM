module EMI_ColumnType_Exchange
  
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod         , only : emi_data_list, emi_data
  use EMI_DataDimensionMod , only : emi_data_dimension_list_type
  use EMI_ColumnType_Constants
  
  !
  implicit none
  !
  !
  public :: EMI_Pack_ColumnType_for_EM

contains

  !-----------------------------------------------------------------------
  subroutine EMI_Pack_ColumnType_for_EM(data_list, em_stage, &
        num_filter, filter)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM's column type for EM
    !
    ! !USES:
    use ColumnType                , only : col_pp
    use elm_varpar                , only : nlevgrnd, nlevsno
    !
    implicit none
    !
    class(emi_data_list) , intent(in) :: data_list
    integer              , intent(in) :: em_stage
    integer              , intent(in) :: num_filter ! number of column soil points in column filter
    integer              , intent(in) :: filter(:)  ! column filter for soil points
    !
    integer                           :: fc,c,j
    class(emi_data), pointer          :: cur_data
    logical                           :: need_to_pack
    integer                           :: istage
    integer                           :: count

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

          case (L2E_COLUMN_ACTIVE)
             do fc = 1, num_filter
                c = filter(fc)
                if (col_pp%active(c)) cur_data%data_int_1d(c) = 1
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_TYPE)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_int_1d(c) = col_pp%itype(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_LANDUNIT_INDEX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_int_1d(c) = col_pp%landunit(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_ZI)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 0, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%zi(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_DZ)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%dz(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_Z)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%z(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_AREA)
             do fc = 1, num_filter
                c = filter(fc)
                !cur_data%data_real_1d(c) = col_pp%area(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_GRIDCELL_INDEX)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_int_1d(c) = col_pp%gridcell(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_PATCH_INDEX_BEGIN)
             do fc = 1, num_filter
                c = filter(fc)
#ifndef FATES_VIA_EMI
                cur_data%data_int_1d(c) = col_pp%pfti(c)
#else
                cur_data%data_int_1d(c) = col_pp%patchi(c)
#endif
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_PATCH_INDEX_END)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_int_1d(c) = col_pp%pftf(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_NUM_SNOW_LAYERS)
             do fc = 1, num_filter
                c = filter(fc)
                cur_data%data_int_1d(c) = col_pp%snl(c)
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_ZI_SNOW_AND_SOIL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%zi(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_DZ_SNOW_AND_SOIL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno+1, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%dz(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_COLUMN_Z_SNOW_AND_SOIL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = -nlevsno+1, nlevgrnd
                   cur_data%data_real_2d(c,j) = col_pp%z(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

  end subroutine EMI_Pack_ColumnType_for_EM

end module EMI_ColumnType_Exchange

