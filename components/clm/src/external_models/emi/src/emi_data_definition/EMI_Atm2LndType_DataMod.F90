module EMI_Atm2LndType_DataMod
  !
  use EMI_Atm2LndType_Constants
  !
  implicit none
  !
  public :: EMI_Atm2LndType_DataInfoByID

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_Atm2LndType_DataInfoByID(data_id, id_val, name_val, long_name_val,&
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

    case(L2E_STATE_FORC_PBOT_DOWNSCALED)
       id_val         =  L2E_STATE_FORC_PBOT_DOWNSCALED
       name_val       =  'Downscaled atm pressure'
       long_name_val  =  'Downscaled atm pressure: ELM to EM'
       units_val      =  '[Pa]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_FORC_T_DOWNSCALED)
       id_val         =  L2E_STATE_FORC_T_DOWNSCALED
       name_val       =  'Downscaled atm temperature'
       long_name_val  =  'Downscaled atm temperature: ELM to EM'
       units_val      =  '[K]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_FLUX_SOLAR_DIRECT_RADDIATION)
       id_val         =  L2E_FLUX_SOLAR_DIRECT_RADDIATION
       name_val       =  'Incident direct solar radiation'
       long_name_val  =  'Incident direct solar radiation: ELM to EM'
       units_val      =  '[W/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begg
       dim1_end_name  =  dimname_endg
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_two
       data_found   =  .true.

    case(L2E_FLUX_SOLAR_DIFFUSE_RADDIATION)
       id_val         =  L2E_FLUX_SOLAR_DIFFUSE_RADDIATION
       name_val       =  'Incident diffuse solar radiation'
       long_name_val  =  'Incident diffuse solar radiation: ELM to EM'
       units_val      =  '[W/m2]'
       is_real_type   =  .true.
       ndim           =  2
       dim1_beg_name  =  dimname_begg
       dim1_end_name  =  dimname_endg
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_two
       data_found   =  .true.
    end select
    
  end subroutine EMI_Atm2LndType_DataInfoByID
    
end module EMI_Atm2LndType_DataMod
