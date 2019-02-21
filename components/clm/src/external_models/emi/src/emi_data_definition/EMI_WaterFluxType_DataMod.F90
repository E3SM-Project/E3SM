module EMI_WaterFluxType_DataMod
  !
  use EMI_WaterFluxType_Constants
  !
  implicit none
  !
  public :: EMI_WaterFluxType_DataInfoByID

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_WaterFluxType_DataInfoByID(data_id, id_val, name_val, long_name_val,&
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

    case(L2E_FLUX_INFIL_MASS_FLUX)
       id_val         =  L2E_FLUX_INFIL_MASS_FLUX
       name_val       =  'Soil infiltration source'
       long_name_val  =  'Soil infiltration source: ELM to EM'
       units_val      =  '[kg/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_VERTICAL_ET_MASS_FLUX)
       id_val         =  L2E_FLUX_VERTICAL_ET_MASS_FLUX
       name_val       =  'Evapotranspiration sink'
       long_name_val  =  'Evapotranspiration sink: ELM to EM'
       units_val      =  '[kg/s]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_FLUX_DEW_MASS_FLUX)
       id_val         =  L2E_FLUX_DEW_MASS_FLUX
       name_val       =  'Dew sink'
       long_name_val  =  'Dew sink: ELM to EM'
       units_val      =  '[kg/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX)
       id_val         =  L2E_FLUX_SNOW_SUBLIMATION_MASS_FLUX
       name_val       =  'Snow layer disappearance sink'
       long_name_val  =  'Snow layer disappearance sink: ELM to EM'
       units_val      =  '[kg/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
       id_val         =  L2E_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
       name_val       =  'Snow layer disappearance sink'
       long_name_val  =  'Snow layer disappearance sink: ELM to EM'
       units_val      =  '[kg/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
       id_val         =  L2E_FLUX_RESTART_SNOW_LYR_DISAPPERANCE_MASS_FLUX
       name_val       =  'Snow layer disappearance sink needed at restart'
       long_name_val  =  'Snow layer disappearance sink needed at restart: ELM to EM'
       units_val      =  '[kg/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_DRAINAGE_MASS_FLUX)
       id_val         =  L2E_FLUX_DRAINAGE_MASS_FLUX
       name_val       =  'Drainage sink'
       long_name_val  =  'Drainage sink: ELM to EM'
       units_val      =  '[kg/s]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_FLUX_INFL)
       id_val         =  L2E_FLUX_INFL
       name_val       =  'Infiltration'
       long_name_val  =  'Infiltration: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_TOTDRAIN)
       id_val         =  L2E_FLUX_TOTDRAIN
       name_val       =  'Drainage sink'
       long_name_val  =  'Drainage sink: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_GROSS_EVAP_SOIL)
       id_val         =  L2E_FLUX_GROSS_EVAP_SOIL
       name_val       =  'Gross evaporation from soil'
       long_name_val  =  'Gross evaporation from soil: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_GROSS_INFL_SOIL)
       id_val         =  L2E_FLUX_GROSS_INFL_SOIL
       name_val       =  'Gross evaporation infiltration'
       long_name_val  =  'Gross evaporation infiltration: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_SURF)
       id_val         =  L2E_FLUX_SURF
       name_val       =  'Surface runoff'
       long_name_val  =  'Surface runoff: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_DEW_GRND)
       id_val         =  L2E_FLUX_DEW_GRND
       name_val       =  'Surface dew added to ground'
       long_name_val  =  'Surface dew added to ground: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_DEW_SNOW)
       id_val         =  L2E_FLUX_DEW_SNOW
       name_val       =  'Surface dew added to snow pacK'
       long_name_val  =  'Surface dew added to snow pacK: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_SUB_SNOW_VOL)
       id_val         =  L2E_FLUX_SUB_SNOW_VOL
       name_val       =  'Sublimation vol from snow pack'
       long_name_val  =  'Sublimation vol from snow pack: ELM to EM'
       units_val      =  '[mm H2O]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_SUB_SNOW)
       id_val         =  L2E_FLUX_SUB_SNOW
       name_val       =  'sublimation rate from snow pack'
       long_name_val  =  'sublimation rate from snow pack: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_H2OSFC2TOPSOI)
       id_val         =  L2E_FLUX_H2OSFC2TOPSOI
       name_val       =  'Rate of surface standing water entering top soil'
       long_name_val  =  'Rate of surface standing water entering top soil: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_SNOW2TOPSOI)
       id_val         =  L2E_FLUX_SNOW2TOPSOI
       name_val       =  'Rate of snow entering top soil'
       long_name_val  =  'Rate of snow entering top soil: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_ROOTSOI)
       id_val         =  L2E_FLUX_ROOTSOI
       name_val       =  'Root and soil water exchange'
       long_name_val  =  'Root and soil water exchange: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_FLUX_ADV)
       id_val         =  L2E_FLUX_ADV
       name_val       =  'Advective flux across different soil layer interfaces'
       long_name_val  =  'Advective flux across different soil layer interfaces: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_zero
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_FLUX_DRAIN_VR)
       id_val         =  L2E_FLUX_DRAIN_VR
       name_val       =  'liquid water losted as drainage'
       long_name_val  =  'liquid water losted as drainage: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_FLUX_TRAN_VEG)
       id_val         =  L2E_FLUX_TRAN_VEG
       name_val       =  'vegetation transpiration'
       long_name_val  =  'vegetation transpiration: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_ROOTSOI_FRAC)
       id_val         =  L2E_FLUX_ROOTSOI_FRAC
       name_val       =  'Root soil fraction'
       long_name_val  =  'Root soil fraction: ELM to EM'
       units_val      =  '[mm H2O/s]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begp
       dim1_end_name  =  dimname_endp
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevsoi
       data_found   =  .true.

    case(E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX)
       id_val         =  E2L_FLUX_SNOW_LYR_DISAPPERANCE_MASS_FLUX
       name_val       =  'Snow layer disappearance sink'
       long_name_val  =  'Snow layer disappearance sink: EM to ELM'
       units_val      =  '[kg/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.
    end select
    
  end subroutine EMI_WaterFluxType_DataInfoByID
    
end module EMI_WaterFluxType_DataMod
