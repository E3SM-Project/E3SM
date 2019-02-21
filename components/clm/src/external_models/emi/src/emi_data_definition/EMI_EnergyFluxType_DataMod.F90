module EMI_EnergyFluxType_DataMod
  !
  use EMI_EnergyFluxType_Constants
  !
  implicit none
  !
  public :: EMI_EnergyFluxType_DataInfoByID

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_EnergyFluxType_DataInfoByID(data_id, id_val, name_val, long_name_val,&
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

    case(L2E_FLUX_ABSORBED_SOLAR_RADIATION)
       id_val         =  L2E_FLUX_ABSORBED_SOLAR_RADIATION
       name_val       =  'Absorbed solar radiation'
       long_name_val  =  'Absorbed solar radiation: ELM to EM'
       units_val      =  '[Wm^{-2}]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_nlevsno_plus_one
       dim2_end_name  =  dimname_one
       data_found   =  .true.

    case(L2E_FLUX_SOIL_HEAT_FLUX)
       id_val         =  L2E_FLUX_SOIL_HEAT_FLUX
       name_val       =  'Soil heat flux'
       long_name_val  =  'Soil heat flux: ELM to EM'
       units_val      =  '[Wm^{-2}]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_SNOW_HEAT_FLUX)
       id_val         =  L2E_FLUX_SNOW_HEAT_FLUX
       name_val       =  'Heat flux on snow'
       long_name_val  =  'Heat flux on snow: ELM to EM'
       units_val      =  '[Wm^{-2}]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_H2OSFC_HEAT_FLUX)
       id_val         =  L2E_FLUX_H2OSFC_HEAT_FLUX
       name_val       =  'Heat flux on water'
       long_name_val  =  'Heat flux on water: ELM to EM'
       units_val      =  '[Wm^{-2}]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX)
       id_val         =  L2E_FLUX_DERIVATIVE_OF_HEAT_FLUX
       name_val       =  'Derivative of heat flux'
       long_name_val  =  'Derivative of heat flux: ELM to EM'
       units_val      =  '[Wm^{-2}]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.
    end select
    
  end subroutine EMI_EnergyFluxType_DataInfoByID
    
end module EMI_EnergyFluxType_DataMod
