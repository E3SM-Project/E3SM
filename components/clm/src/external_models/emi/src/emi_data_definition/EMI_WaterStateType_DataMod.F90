module EMI_WaterStateType_DataMod
  !
  use EMI_WaterStateType_Constants
  !
  implicit none
  !
  public :: EMI_WaterStateType_DataInfoByID

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_WaterStateType_DataInfoByID(data_id, id_val, name_val, long_name_val,&
        units_val, is_int_type, is_real_type, ndim, &
        dim1_beg_name, dim1_end_name, dim2_beg_name, dim2_end_name, &
        dim3_beg_name, dim3_end_name, dim4_beg_name, dim4_end_name, &
        data_found)
    !
    ! !DESCRIPTION:
    ! Defines information of data exchanged between ELM and EM
    !
    ! !USES: 
    use EMI_DataDimensionMod
    implicit none
    !
    ! !ARGUMENTS:
    integer            , intent(in)  :: data_id
    integer            , intent(out) :: id_val
    character (len=*)  , intent(out) :: name_val
    character (len=128), intent(out) :: long_name_val
    character (len=24) , intent(out) :: units_val
    logical            , intent(out) :: is_int_type
    logical            , intent(out) :: is_real_type
    integer            , intent(out) :: ndim
    character (len=24) , intent(out) :: dim1_beg_name
    character (len=24) , intent(out) :: dim1_end_name
    character (len=24) , intent(out) :: dim2_beg_name
    character (len=24) , intent(out) :: dim2_end_name
    character (len=24) , intent(out) :: dim3_beg_name
    character (len=24) , intent(out) :: dim3_end_name
    character (len=24) , intent(out) :: dim4_beg_name
    character (len=24) , intent(out) :: dim4_end_name
    logical            , intent(out) :: data_found

    is_int_type    = .false.
    is_real_type   = .false.
    dim1_beg_name  = ''
    dim2_beg_name  = ''
    dim3_beg_name  = ''
    dim4_beg_name  = ''
    dim1_end_name  = ''
    dim2_end_name  = ''
    dim3_end_name  = ''
    dim4_end_name  = ''

    select case(data_id)

    case(L2E_STATE_H2OSOI_LIQ_NLEVGRND)
       id_val         =  L2E_STATE_H2OSOI_LIQ_NLEVGRND
       name_val       =  'soil liq water'
       long_name_val  =  'soil liq water: ELM to EM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_ICE_NLEVGRND)
       id_val         =  L2E_STATE_H2OSOI_ICE_NLEVGRND
       name_val       =  'soil ice water'
       long_name_val  =  'soil ice water: ELM to EM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_STATE_VSFM_PROGNOSTIC_SOILP)
       id_val         =  L2E_STATE_VSFM_PROGNOSTIC_SOILP
       name_val       =  'Soil pressure'
       long_name_val  =  'Soil pressure: ELM to EM'
       units_val      =  '[Pa]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_STATE_FRAC_H2OSFC)
       id_val         =  L2E_STATE_FRAC_H2OSFC
       name_val       =  'Fraction of standing water'
       long_name_val  =  'Fraction of standing water: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_FRAC_INUNDATED)
       id_val         =  L2E_STATE_FRAC_INUNDATED
       name_val       =  'Fraction inundated'
       long_name_val  =  'Fraction inundated: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI)
       id_val         =  L2E_STATE_H2OSOI_LIQ_VOL_NLEVSOI
       name_val       =  'Vol sol liq water'
       long_name_val  =  'Vol sol liq water: ELM to EM'
       units_val      =  '[m3/m3]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI)
       id_val         =  L2E_STATE_H2OSOI_ICE_VOL_NLEVSOI
       name_val       =  'Vol oil liq ice'
       long_name_val  =  'Vol oil liq ice: ELM to EM'
       units_val      =  '[m3/m3]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_VOL_NLEVSOI)
       id_val         =  L2E_STATE_H2OSOI_VOL_NLEVSOI
       name_val       =  'Vol. soil liq water'
       long_name_val  =  'Vol. soil liq water: ELM to EM'
       units_val      =  '[m3/m3]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(L2E_STATE_AIR_VOL_NLEVSOI)
       id_val         =  L2E_STATE_AIR_VOL_NLEVSOI
       name_val       =  'Air vol'
       long_name_val  =  'Air vol: ELM to EM'
       units_val      =  '[m3/m3]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(L2E_STATE_RHO_VAP_NLEVSOI)
       id_val         =  L2E_STATE_RHO_VAP_NLEVSOI
       name_val       =  'Water vapor pressure'
       long_name_val  =  'Water vapor pressure: ELM to EM'
       units_val      =  '[Pa]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(L2E_STATE_RHVAP_SOI_NLEVSOI)
       id_val         =  L2E_STATE_RHVAP_SOI_NLEVSOI
       name_val       =  'Relative humidity'
       long_name_val  =  'Relative humidity: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI)
       id_val         =  L2E_STATE_SOIL_MATRIC_POTENTIAL_NLEVSOI
       name_val       =  'Matric potential'
       long_name_val  =  'Matric potential: ELM to EM'
       units_val      =  '[mm]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_LIQ_NLEVSOI)
       id_val         =  L2E_STATE_H2OSOI_LIQ_NLEVSOI
       name_val       =  'Soil liquid water'
       long_name_val  =  'Soil liquid water: ELM to EM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_ICE_NLEVSOI)
       id_val         =  L2E_STATE_H2OSOI_ICE_NLEVSOI
       name_val       =  'Soil ice water'
       long_name_val  =  'Soil ice water: ELM to EM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_LIQ_NLEVSNOW)
       id_val         =  L2E_STATE_H2OSOI_LIQ_NLEVSNOW
       name_val       =  'Soil liquid water'
       long_name_val  =  'Soil liquid water: ELM to EM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_nlevsno_plus_one
       dim2_end_name  =  dimname_zero
       data_found   =  .true.

    case(L2E_STATE_H2OSOI_ICE_NLEVSNOW)
       id_val         =  L2E_STATE_H2OSOI_ICE_NLEVSNOW
       name_val       =  'Soil ice water'
       long_name_val  =  'Soil ice water: ELM to EM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_nlevsno_plus_one
       dim2_end_name  =  dimname_zero
       data_found   =  .true.

    case(L2E_STATE_H2OSNOW)
       id_val         =  L2E_STATE_H2OSNOW
       name_val       =  'Snow water'
       long_name_val  =  'Snow water: ELM to EM'
       units_val      =  '[mm]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_H2OSFC)
       id_val         =  L2E_STATE_H2OSFC
       name_val       =  'Standing surface water'
       long_name_val  =  'Standing surface water: ELM to EM'
       units_val      =  '[mm]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_FRAC_SNOW_EFFECTIVE)
       id_val         =  L2E_STATE_FRAC_SNOW_EFFECTIVE
       name_val       =  'Frac. of grnd covered with snow'
       long_name_val  =  'Frac. of grnd covered with snow: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(E2L_STATE_H2OSOI_LIQ)
       id_val         =  E2L_STATE_H2OSOI_LIQ
       name_val       =  'Soil liquid water'
       long_name_val  =  'Soil liquid water: EM to ELM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(E2L_STATE_H2OSOI_ICE)
       id_val         =  E2L_STATE_H2OSOI_ICE
       name_val       =  'Soil ice water'
       long_name_val  =  'Soil ice water: EM to ELM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(E2L_STATE_VSFM_PROGNOSTIC_SOILP)
       id_val         =  E2L_STATE_VSFM_PROGNOSTIC_SOILP
       name_val       =  'Soil matric pressure'
       long_name_val  =  'Soil matric pressure: EM to ELM'
       units_val      =  '[Pa]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.
    end select
    
  end subroutine EMI_WaterStateType_DataInfoByID
    
end module EMI_WaterStateType_DataMod
