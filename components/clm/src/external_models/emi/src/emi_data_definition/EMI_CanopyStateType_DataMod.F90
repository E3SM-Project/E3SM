module EMI_CanopyStateType_DataMod
  !
  use EMI_CanopyStateType_Constants
  !
  implicit none
  !
  public :: EMI_CanopyStateType_DataInfoByID

contains
  
!-----------------------------------------------------------------------
  subroutine EMI_CanopyStateType_DataInfoByID(data_id, id_val, name_val, long_name_val,&
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

    case(L2E_STATE_ALTMAX)
       id_val         =  L2E_STATE_ALTMAX
       name_val       =  'Maximum annual depth of thaw'
       long_name_val  =  'Maximum annual depth of thaw: ELM to EM'
       units_val      =  '[m]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_ALTMAX_LASTYEAR)
       id_val         =  L2E_STATE_ALTMAX_LASTYEAR
       name_val       =  'Prior year maximum annual depth of thaw'
       long_name_val  =  'Prior year maximum annual depth of thaw: ELM to EM'
       units_val      =  '[m]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begc
       dim1_end_name  =  dimname_endc
       data_found   =  .true.

    case(L2E_STATE_LBL_RSC_H2O)
       id_val         =  L2E_STATE_LBL_RSC_H2O
       name_val       =  'Laminar boundary layer resistance for water over dry leaf'
       long_name_val  =  'Laminar boundary layer resistance for water over dry leaf: ELM to EM'
       units_val      =  '[s/m]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begp
       dim1_end_name  =  dimname_endp
       data_found   =  .true.

    case(L2E_STATE_ELAI)
       id_val         =  L2E_STATE_ELAI
       name_val       =  'Canopy one-sided leaf area index with burying by snow'
       long_name_val  =  'Canopy one-sided leaf area index with burying by snow: ELM to EM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begp
       dim1_end_name  =  dimname_endp
       data_found   =  .true.

    case(E2L_STATE_FSUN)
       id_val         =  E2L_STATE_FSUN
       name_val       =  'Sunlit fraction'
       long_name_val  =  'Sunlit fraction: EM to ELM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begp
       dim1_end_name  =  dimname_endp
       data_found   =  .true.

    case(E2L_STATE_LAISUN)
       id_val         =  E2L_STATE_LAISUN
       name_val       =  'Sunlit projected leaf area index'
       long_name_val  =  'Sunlit projected leaf area index: EM to ELM'
       units_val      =  '[-]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begp
       dim1_end_name  =  dimname_endp
       data_found   =  .true.

    case(E2L_STATE_LAISHA)
       id_val         =  E2L_STATE_LAISHA
       name_val       =  'Shaded projected leaf area index'
       long_name_val  =  'Shaded projected leaf area index: EM to ELM'
       units_val      =  '[]'
       is_real_type   =  .true.
       ndim           =  1
       dim1_beg_name  =  dimname_begp
       dim1_end_name  =  dimname_endp
       data_found   =  .true.
    end select
    
  end subroutine EMI_CanopyStateType_DataInfoByID
    
end module EMI_CanopyStateType_DataMod
