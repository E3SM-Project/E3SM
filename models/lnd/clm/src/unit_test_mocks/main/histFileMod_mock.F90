module histFileMod

  ! This is a stub for histFileMod. Currently all it does is provide empty
  ! implementations for hist_addfld calls, to satisfy the interface that is expected
  ! throughout the CLM code.

  use shr_kind_mod, only : r8 => shr_kind_r8

  implicit none
  private
  save

  integer , public , parameter :: no_snow_normal = 1              ! normal treatment, which should be used for most fields (use spval when snow layer not present)

  public :: hist_addfld1d
  public :: hist_addfld2d
  public :: hist_addfld_decomp

contains
  
  subroutine hist_addfld1d (fname, units, avgflag, long_name, type1d_out, &
                        ptr_gcell, ptr_lunit, ptr_col, ptr_patch, ptr_lnd, &
                        ptr_atm, p2c_scale_type, c2l_scale_type, &
                        l2g_scale_type, set_lake, set_nolake, set_urb, set_nourb, &
                        set_noglcmec, set_spec, default)
      character(len=*), intent(in)           :: fname          ! field name
    character(len=*), intent(in)           :: units          ! units of field
    character(len=1), intent(in)           :: avgflag        ! time averaging flag
    character(len=*), intent(in)           :: long_name      ! long name of field
    character(len=*), optional, intent(in) :: type1d_out     ! output type (from clmtype)
    real(r8)        , optional, pointer    :: ptr_gcell(:)   ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_lunit(:)   ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:)     ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_patch(:)     ! pointer to pft array
    real(r8)        , optional, pointer    :: ptr_lnd(:)     ! pointer to lnd array
    real(r8)        , optional, pointer    :: ptr_atm(:)     ! pointer to atm array
    real(r8)        , optional, intent(in) :: set_lake       ! value to set lakes to
    real(r8)        , optional, intent(in) :: set_nolake     ! value to set non-lakes to
    real(r8)        , optional, intent(in) :: set_urb        ! value to set urban to
    real(r8)        , optional, intent(in) :: set_nourb      ! value to set non-urban to
    real(r8)        , optional, intent(in) :: set_noglcmec   ! value to set non-glacier_mec to
    real(r8)        , optional, intent(in) :: set_spec       ! value to set special to
    character(len=*), optional, intent(in) :: p2c_scale_type ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape

    ! Do nothing

  end subroutine hist_addfld1d

  
  subroutine hist_addfld2d (fname, type2d, units, avgflag, long_name, type1d_out, &
                        ptr_gcell, ptr_lunit, ptr_col, ptr_patch, ptr_lnd, ptr_atm, &
                        p2c_scale_type, c2l_scale_type, l2g_scale_type, &
                        set_lake, set_nolake, set_urb, set_nourb, set_spec, &
                        no_snow_behavior, default)

    character(len=*), intent(in) :: fname                      ! field name
    character(len=*), intent(in) :: type2d                     ! 2d output type
    character(len=*), intent(in) :: units                      ! units of field
    character(len=1), intent(in) :: avgflag                    ! time averaging flag
    character(len=*), intent(in) :: long_name                  ! long name of field
    character(len=*), optional, intent(in) :: type1d_out       ! output type (from clmtype)
    real(r8)        , optional, pointer    :: ptr_atm(:,:)     ! pointer to atm array
    real(r8)        , optional, pointer    :: ptr_lnd(:,:)     ! pointer to lnd array
    real(r8)        , optional, pointer    :: ptr_gcell(:,:)   ! pointer to gridcell array
    real(r8)        , optional, pointer    :: ptr_lunit(:,:)   ! pointer to landunit array
    real(r8)        , optional, pointer    :: ptr_col(:,:)     ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_patch(:,:)     ! pointer to pft array
    real(r8)        , optional, intent(in) :: set_lake         ! value to set lakes to
    real(r8)        , optional, intent(in) :: set_nolake       ! value to set non-lakes to
    real(r8)        , optional, intent(in) :: set_urb          ! value to set urban to
    real(r8)        , optional, intent(in) :: set_nourb        ! value to set non-urban to
    real(r8)        , optional, intent(in) :: set_spec         ! value to set special to
    integer         , optional, intent(in) :: no_snow_behavior ! if a multi-layer snow field, behavior to use for absent snow layers (should be one of the public no_snow_* parameters defined above)
    character(len=*), optional, intent(in) :: p2c_scale_type   ! scale type for subgrid averaging of pfts to column
    character(len=*), optional, intent(in) :: c2l_scale_type   ! scale type for subgrid averaging of columns to landunits
    character(len=*), optional, intent(in) :: l2g_scale_type   ! scale type for subgrid averaging of landunits to gridcells
    character(len=*), optional, intent(in) :: default          ! if set to 'inactive, field will not appear on primary tape

    ! Do nothing

  end subroutine hist_addfld2d

  subroutine hist_addfld_decomp (fname, type2d, units, avgflag, long_name, ptr_col, ptr_patch, default)

    !
    ! !USES:
    use clm_varpar  , only : nlevdecomp_full, crop_prog
    use clm_varctl  , only : iulog
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fname                    ! field name
    character(len=*), intent(in) :: type2d                   ! 2d output type
    character(len=*), intent(in) :: units                    ! units of field
    character(len=1), intent(in) :: avgflag                  ! time averaging flag
    character(len=*), intent(in) :: long_name                ! long name of field
    real(r8)        , optional, pointer    :: ptr_col(:,:)   ! pointer to column array
    real(r8)        , optional, pointer    :: ptr_patch(:,:)   ! pointer to pft array
    character(len=*), optional, intent(in) :: default        ! if set to 'inactive, field will not appear on primary tape

    ! Do nothing

  end subroutine hist_addfld_decomp

end module histFileMod
