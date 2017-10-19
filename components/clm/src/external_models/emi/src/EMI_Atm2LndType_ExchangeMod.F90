module EMI_Atm2LndType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use clm_varctl                            , only : iulog
  use ExternalModelInterfaceDataMod         , only : emi_data_list, emi_data
  use ExternalModelIntefaceDataDimensionMod , only : emi_data_dimension_list_type
  use Atm2LndType                           , only : atm2lnd_type
  !
  implicit none
  !
  !
  public :: EMI_Pack_Atm2LndType_at_Grid_Level_for_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_Atm2LndType_at_Grid_Level_for_EM(data_list, em_stage, &
        num_filter, filter, atm2lndtype_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM atm2lndtype_vars for EM
    !
    ! !USES:
    use ExternalModelConstants , only : L2E_FLUX_SOLAR_DIRECT_RADDIATION
    use ExternalModelConstants , only : L2E_FLUX_SOLAR_DIFFUSE_RADDIATION
    use clm_varpar             , only : nlevsoi, nlevgrnd, nlevsno
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
