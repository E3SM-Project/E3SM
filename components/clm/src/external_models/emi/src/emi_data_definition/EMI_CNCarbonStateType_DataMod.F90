module EMI_CNCarbonStateType_DataMod
  !
  use EMI_CNCarbonStateType_Constants
  !
  implicit none
  !
  public :: EMI_CNCarbonStateType_DataInfoByID

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_CNCarbonStateType_DataInfoByID(data_id, id_val, name_val, long_name_val,&
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
    character (len=24) , intent(out) :: name_val
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

    case(L2E_STATE_CARBON_POOLS_VERTICALLY_RESOLVED)
       id_val         =  L2E_STATE_CARBON_POOLS_VERTICALLY_RESOLVED
       name_val       =  'decomp cpools vr'
       long_name_val  =  'decomp cpools vr: ELM to EM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  3
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevdecomp_full
       dim3_beg_name  =  dimname_one
       dim3_end_name  =  dimname_ndecomp_pools
       data_found   =  .true.

    case(E2L_STATE_CARBON_POOLS_VERTICALLY_RESOLVED)
       id_val         =  E2L_STATE_CARBON_POOLS_VERTICALLY_RESOLVED
       name_val       =  'decomp cpools vr'
       long_name_val  =  'decomp cpools vr: EM to ELM'
       units_val      =  '[kg/m2]'
       is_real_type   =  .true.
       ndim           =  3
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       dim2_beg_name  =  dimname_one
       dim2_end_name  =  dimname_nlevdecomp_full
       dim3_beg_name  =  dimname_one
       dim3_end_name  =  dimname_ndecomp_pools
       data_found   =  .true.
    end select
    
  end subroutine EMI_CNCarbonStateType_DataInfoByID
    
end module EMI_CNCarbonStateType_DataMod
