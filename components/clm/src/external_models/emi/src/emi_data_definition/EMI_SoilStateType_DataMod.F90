module EMI_SoilStateType_DataMod
  !
  use EMI_SoilStateType_Constants
  !
  implicit none
  !
  public :: EMI_SoilStateType_DataInfoByID

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_SoilStateType_DataInfoByID(data_id, id_val, name_val, long_name_val,&
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
    character (len=*)  , intent(out) :: long_name_val
    character (len=*)  , intent(out) :: units_val
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

    case(L2E_PARAMETER_WATSATC)
       id_val         =  L2E_PARAMETER_WATSATC
       name_val       =  'Soil porosity'
       long_name_val  =  'Soil porosity: ELM to EM'
       units_val      =  '[m^3/m^3]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_HKSATC)
       id_val         =  L2E_PARAMETER_HKSATC
       name_val       =  'Soil hydraulic conductivity'
       long_name_val  =  'Soil hydraulic conductivity: ELM to EM'
       units_val      =  '[mm/s]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_BSWC)
       id_val         =  L2E_PARAMETER_BSWC
       name_val       =  'Clapp and Hornberger parameter'
       long_name_val  =  'Clapp and Hornberger parameter: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_SUCSATC)
       id_val         =  L2E_PARAMETER_SUCSATC
       name_val       =  'Minimum soil sucsatc'
       long_name_val  =  'Minimum soil sucsatc: ELM to EM'
       units_val      =  '[mm]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_EFFPOROSITYC)
       id_val         =  L2E_PARAMETER_EFFPOROSITYC
       name_val       =  'Effective porosity'
       long_name_val  =  'Effective porosity: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_CSOL)
       id_val         =  L2E_PARAMETER_CSOL
       name_val       =  'Heat capacity'
       long_name_val  =  'Heat capacity: ELM to EM'
       units_val      =  '[W/m/K]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_TKMG)
       id_val         =  L2E_PARAMETER_TKMG
       name_val       =  'Thermal conductivity minearls'
       long_name_val  =  'Thermal conductivity minearls: ELM to EM'
       units_val      =  '[W/m/K]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_TKDRY)
       id_val         =  L2E_PARAMETER_TKDRY
       name_val       =  'Thermal conductivity dry soils'
       long_name_val  =  'Thermal conductivity dry soils: ELM to EM'
       units_val      =  '[W/m/K]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_CELLORG)
       id_val         =  L2E_PARAMETER_CELLORG
       name_val       =  'Organic matter'
       long_name_val  =  'Organic matter: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_CELLCLAY)
       id_val         =  L2E_PARAMETER_CELLCLAY
       name_val       =  'Clay value'
       long_name_val  =  'Clay value: ELM to EM'
       units_val      =  '[%]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_CELLSAND)
       id_val         =  L2E_PARAMETER_CELLSAND
       name_val       =  'Sand value'
       long_name_val  =  'Sand value: ELM to EM'
       units_val      =  '[%]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_BD)
       id_val         =  L2E_PARAMETER_BD
       name_val       =  'Bulk density'
       long_name_val  =  'Bulk density: ELM to EM'
       units_val      =  '[kg/m^3]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_WATFC)
       id_val         =  L2E_PARAMETER_WATFC
       name_val       =  'Volumetric soil water at field capacity'
       long_name_val  =  'Volumetric soil water at field capacity: ELM to EM'
       units_val      =  '[m^3/m^3]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(L2E_PARAMETER_ROOTFR_PATCH)
       id_val         =  L2E_PARAMETER_ROOTFR_PATCH
       name_val       =  'Fraction of roots in each soil layer'
       long_name_val  =  'Fraction of roots in each soil layer: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begp
       dim1_end_name  =  dimname_endp
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(E2L_STATE_SOIL_MATRIC_POTENTIAL)
       id_val         =  E2L_STATE_SOIL_MATRIC_POTENTIAL
       name_val       =  'Soil matric potential'
       long_name_val  =  ': EM to ELM'
       units_val      =  '[Pa]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.
    end select
    
  end subroutine EMI_SoilStateType_DataInfoByID
    
end module EMI_SoilStateType_DataMod
