module EMI_SoilStateType_ExchangeMod
  !
  use shr_kind_mod                          , only : r8 => shr_kind_r8
  use shr_log_mod                           , only : errMsg => shr_log_errMsg
  use abortutils                            , only : endrun
  use elm_varctl                            , only : iulog
  use EMI_DataMod         , only : emi_data_list, emi_data
  use EMI_DataDimensionMod , only : emi_data_dimension_list_type
  use SoilStateType                         , only : soilstate_type
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
  public :: EMI_Pack_SoilStateType_at_Column_Level_for_EM
  public :: EMI_Pack_SoilStateType_at_Patch_Level_for_EM
  public :: EMI_Unpack_SoilStateType_at_Column_Level_from_EM

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Pack_SoilStateType_at_Column_Level_for_EM(data_list, em_stage, &
        num_filter, filter, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM soilstate_vars for EM
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
    type(soilstate_type)   , intent(in) :: soilstate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         watsat       => soilstate_vars%watsat_col       , &
         hksat        => soilstate_vars%hksat_col        , &
         bsw          => soilstate_vars%bsw_col          , &
         sucsat       => soilstate_vars%sucsat_col       , &
         eff_porosity => soilstate_vars%eff_porosity_col , &
         csol         => soilstate_vars%csol_col         , &
         tkmg         => soilstate_vars%tkmg_col         , &
         tkdry        => soilstate_vars%tkdry_col        , &
         cellorg      => soilstate_vars%cellorg_col      , &
         cellclay     => soilstate_vars%cellclay_col     , &
         cellsand     => soilstate_vars%cellsand_col     , &
         bd           => soilstate_vars%bd_col           , &
         watfc        => soilstate_vars%watfc_col          &
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

          case (L2E_PARAMETER_WATSATC)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = watsat(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_HKSATC)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = hksat(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_BSWC)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = bsw(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_SUCSATC)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = sucsat(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_EFFPOROSITYC)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = eff_porosity(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_CSOL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = csol(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_TKMG)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = tkmg(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_TKDRY)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = tkdry(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_CELLORG)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = cellorg(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_CELLCLAY)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = cellclay(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_CELLSAND)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = cellsand(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_BD)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = bd(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          case (L2E_PARAMETER_WATFC)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(c,j) = watfc(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_SoilStateType_at_Column_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Pack_SoilStateType_at_Patch_Level_for_EM(data_list, em_stage, &
        num_filter, filter, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! Pack data from ALM soilstate_vars for EM
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
    type(soilstate_type)   , intent(in) :: soilstate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fp,p,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         rootfr => soilstate_vars%rootfr_patch   &
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

          case (L2E_PARAMETER_ROOTFR_PATCH)
             do fp = 1, num_filter
                p = filter(fp)
                do j = 1, nlevgrnd
                   cur_data%data_real_2d(p,j) = rootfr(p,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Pack_SoilStateType_at_Patch_Level_for_EM

!-----------------------------------------------------------------------
  subroutine EMI_Unpack_SoilStateType_at_Column_Level_from_EM(data_list, em_stage, &
        num_filter, filter, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! Unpack data for ALM soilstate_vars from EM
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
    type(soilstate_type)   , intent(in) :: soilstate_vars
    !
    ! !LOCAL_VARIABLES:
    integer                             :: fc,c,j
    class(emi_data), pointer            :: cur_data
    logical                             :: need_to_pack
    integer                             :: istage
    integer                             :: count

    associate(& 
         smp_l => soilstate_vars%smp_l_col   &
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

          case (E2L_STATE_SOIL_MATRIC_POTENTIAL)
             do fc = 1, num_filter
                c = filter(fc)
                do j = 1, nlevgrnd
                   smp_l(c,j) = cur_data%data_real_2d(c,j)
                enddo
             enddo
             cur_data%is_set = .true.

          end select

       endif

       cur_data => cur_data%next
    enddo

    end associate

  end subroutine EMI_Unpack_SoilStateType_at_Column_Level_from_EM


end module EMI_SoilStateType_ExchangeMod
