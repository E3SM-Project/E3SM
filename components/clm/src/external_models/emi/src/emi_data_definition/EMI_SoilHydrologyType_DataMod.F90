module EMI_SoilHydrologyType_DataMod
  !
  use EMI_SoilHydrologyType_Constants
  !
  implicit none
  !
  public :: EMI_SoilHydrologyType_DataInfoByID

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_SoilHydrologyType_DataInfoByID(data_id, id_val, name_val, long_name_val,&
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

    case(L2E_STATE_WTD)
       id_val         =  L2E_STATE_WTD
       name_val       =  'Water table depth'
       long_name_val  =  'Water table depth: ELM to EM'
       units_val      =  '[m]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_QCHARGE)
       id_val         =  L2E_STATE_QCHARGE
       name_val       =  'Flux recharge rate'
       long_name_val  =  'Flux recharge rate: ELM to EM'
       units_val      =  '[mm/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_FRACICE)
       id_val         =  L2E_STATE_FRACICE
       name_val       =  'Fractional impermeability'
       long_name_val  =  'Fractional impermeability: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevgrnd
       data_found   =  .true.

    case(E2L_STATE_WTD)
       id_val         =  E2L_STATE_WTD
       name_val       =  'Water table depth'
       long_name_val  =  'Water table depth: EM to ELM'
       units_val      =  '[m]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(E2L_FLUX_AQUIFER_RECHARGE)
       id_val         =  E2L_FLUX_AQUIFER_RECHARGE
       name_val       =  'Aquifer recharge rate'
       long_name_val  =  'Aquifer recharge rate: EM to ELM'
       units_val      =  '[mm/s]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.
    end select
    
  end subroutine EMI_SoilHydrologyType_DataInfoByID
    
end module EMI_SoilHydrologyType_DataMod
