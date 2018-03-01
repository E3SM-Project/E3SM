#include "dtypes.h"
module restUtilMod
  use shr_kind_mod, only: r8=>shr_kind_r8, r4 => shr_kind_r4, i4=>shr_kind_i4
  use ncdio_pio, only : file_desc_t
implicit none
save
private

interface restartvar
   !TYPE text,int,double
   !DIMS 0,1,2
   module procedure restartvar_0d_text
   !TYPE text,int,double
   !DIMS 0,1,2
   module procedure restartvar_1d_text
   !TYPE text,int,double
   !DIMS 0,1,2
   module procedure restartvar_2d_text
   !TYPE text,int,double
   !DIMS 0,1,2
   module procedure restartvar_0d_int
   !TYPE text,int,double
   !DIMS 0,1,2
   module procedure restartvar_1d_int
   !TYPE text,int,double
   !DIMS 0,1,2
   module procedure restartvar_2d_int
   !TYPE text,int,double
   !DIMS 0,1,2
   module procedure restartvar_0d_double
   !TYPE text,int,double
   !DIMS 0,1,2
   module procedure restartvar_1d_double
   !TYPE text,int,double
   !DIMS 0,1,2
   module procedure restartvar_2d_double
   module procedure restartvar_2d_double_bounds
end interface restartvar
public :: restartvar

private :: is_restart
contains

  subroutine restartvar_0d_text(&
     ncid, flag, varname, xtype, &
     long_name, units, interpinic_flag, data, readvar, &
     comment, flag_meanings, missing_value, fill_value, &
     imissing_value, ifill_value, flag_values, nvalid_range )

  !----------------------------------------------------
  ! Arguments
  type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
  character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
  character(len=*)  , intent(in)           :: varname          ! variable name
  integer           , intent(in)           :: xtype            ! netcdf data type
  character(len=*)  , intent(in)           :: long_name        ! long name for variable
  character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
  character(len=*)           , intent(inout)        :: data
  logical           , intent(out)          :: readvar          ! was var read?
  character(len=*)  , intent(in), optional :: units            ! long name for variable
  character(len=*)  , intent(in), optional :: comment          ! attribute
  character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
  real(r8)          , intent(in), optional :: missing_value    ! attribute for real
  real(r8)          , intent(in), optional :: fill_value       ! attribute for real
  integer           , intent(in), optional :: imissing_value   ! attribute for int
  integer           , intent(in), optional :: ifill_value      ! attribute for int
  integer           , intent(in), optional :: flag_values(:)   ! attribute for int
  integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int


  end subroutine restartvar_0d_text


  subroutine restartvar_0d_int(&
       ncid, flag, varname, xtype, &
       long_name, units, interpinic_flag, data, readvar, &
       comment, flag_meanings, missing_value, fill_value, &
       imissing_value, ifill_value, flag_values, nvalid_range )

    !----------------------------------------------------
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    integer           , intent(in)           :: xtype            ! netcdf data type
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    integer(i4)           , intent(inout)        :: data
    logical           , intent(out)          :: readvar          ! was var read?
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)          , intent(in), optional :: missing_value    ! attribute for real
    real(r8)          , intent(in), optional :: fill_value       ! attribute for real
    integer           , intent(in), optional :: imissing_value   ! attribute for int
    integer           , intent(in), optional :: ifill_value      ! attribute for int
    integer           , intent(in), optional :: flag_values(:)   ! attribute for int
    integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int
    !
    ! Local variables


  end subroutine restartvar_0d_int

  subroutine restartvar_0d_double(&
      ncid, flag, varname, xtype, &
     long_name, units, interpinic_flag, data, readvar, &
     comment, flag_meanings, missing_value, fill_value, &
     imissing_value, ifill_value, flag_values, nvalid_range )

      !----------------------------------------------------
      ! Arguments
      type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
      character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
      character(len=*)  , intent(in)           :: varname          ! variable name
      integer           , intent(in)           :: xtype            ! netcdf data type
      character(len=*)  , intent(in)           :: long_name        ! long name for variable
      character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
      real(r8)           , intent(inout)        :: data
      logical           , intent(out)          :: readvar          ! was var read?
      character(len=*)  , intent(in), optional :: units            ! long name for variable
      character(len=*)  , intent(in), optional :: comment          ! attribute
      character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
      real(r8)          , intent(in), optional :: missing_value    ! attribute for real
      real(r8)          , intent(in), optional :: fill_value       ! attribute for real
      integer           , intent(in), optional :: imissing_value   ! attribute for int
      integer           , intent(in), optional :: ifill_value      ! attribute for int
      integer           , intent(in), optional :: flag_values(:)   ! attribute for int
      integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int
      !
      ! Local variables



  end subroutine restartvar_0d_double




  subroutine restartvar_1d_text(&
     ncid, flag, varname, xtype, dim1name, dim2name, &
     long_name, units, interpinic_flag, data, readvar, &
     comment, flag_meanings, missing_value, fill_value, &
     imissing_value, ifill_value, flag_values, nvalid_range )

      !----------------------------------------------------
      ! Arguments
      type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
      character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
      character(len=*)  , intent(in)           :: varname          ! variable name
      integer           , intent(in)           :: xtype            ! netcdf data type
      character(len=*)  , intent(in)           :: long_name        ! long name for variable
      character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
      character(len=*)           , pointer              :: data(:)
      logical           , intent(inout)        :: readvar          ! was var read?
      character(len=*)  , intent(in), optional :: dim1name         ! dimension name
      character(len=*)  , intent(in), optional :: dim2name         ! dimension name
      character(len=*)  , intent(in), optional :: units            ! long name for variable
      character(len=*)  , intent(in), optional :: comment          ! attribute
      character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
      real(r8)          , intent(in), optional :: missing_value    ! attribute for real
      real(r8)          , intent(in), optional :: fill_value       ! attribute for real
      integer           , intent(in), optional :: imissing_value   ! attribute for int
      integer           , intent(in), optional :: ifill_value      ! attribute for int
      integer           , intent(in), optional :: flag_values(:)   ! attribute for int
      integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int
      !
      ! Local variables


  end subroutine restartvar_1d_text



  subroutine restartvar_2d_text(&
       ncid, flag, varname, xtype, dim1name, dim2name, &
       long_name, units, interpinic_flag, data, readvar, &
       comment, flag_meanings, missing_value, fill_value, &
       imissing_value, ifill_value, flag_values, nvalid_range )

        !----------------------------------------------------
        ! Arguments
        type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
        character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
        character(len=*)  , intent(in)           :: varname          ! variable name
        integer           , intent(in)           :: xtype            ! netcdf data type
        character(len=*)  , intent(in)           :: long_name        ! long name for variable
        character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
        character(len=*)           , pointer              :: data(:,:)
        logical           , intent(inout)        :: readvar          ! was var read?
        character(len=*)  , intent(in), optional :: dim1name         ! dimension name
        character(len=*)  , intent(in), optional :: dim2name         ! dimension name
        character(len=*)  , intent(in), optional :: units            ! long name for variable
        character(len=*)  , intent(in), optional :: comment          ! attribute
        character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
        real(r8)          , intent(in), optional :: missing_value    ! attribute for real
        real(r8)          , intent(in), optional :: fill_value       ! attribute for real
        integer           , intent(in), optional :: imissing_value   ! attribute for int
        integer           , intent(in), optional :: ifill_value      ! attribute for int
        integer           , intent(in), optional :: flag_values(:)   ! attribute for int
        integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int
        !


  end subroutine restartvar_2d_text



  subroutine restartvar_1d_int(&
         ncid, flag, varname, xtype, dim1name, dim2name, &
         long_name, units, interpinic_flag, data, readvar, &
         comment, flag_meanings, missing_value, fill_value, &
         imissing_value, ifill_value, flag_values, nvalid_range )

          !----------------------------------------------------
          ! Arguments
          type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
          character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
          character(len=*)  , intent(in)           :: varname          ! variable name
          integer           , intent(in)           :: xtype            ! netcdf data type
          character(len=*)  , intent(in)           :: long_name        ! long name for variable
          character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
          integer(i4)           , pointer              :: data(:)
          logical           , intent(inout)        :: readvar          ! was var read?
          character(len=*)  , intent(in), optional :: dim1name         ! dimension name
          character(len=*)  , intent(in), optional :: dim2name         ! dimension name
          character(len=*)  , intent(in), optional :: units            ! long name for variable
          character(len=*)  , intent(in), optional :: comment          ! attribute
          character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
          real(r8)          , intent(in), optional :: missing_value    ! attribute for real
          real(r8)          , intent(in), optional :: fill_value       ! attribute for real
          integer           , intent(in), optional :: imissing_value   ! attribute for int
          integer           , intent(in), optional :: ifill_value      ! attribute for int
          integer           , intent(in), optional :: flag_values(:)   ! attribute for int
          integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int
          !
          ! Local variables


  end subroutine restartvar_1d_int


  subroutine restartvar_2d_int(&
       ncid, flag, varname, xtype, dim1name, dim2name, &
       long_name, units, interpinic_flag, data, readvar, &
       comment, flag_meanings, missing_value, fill_value, &
       imissing_value, ifill_value, flag_values, nvalid_range )

    !----------------------------------------------------
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    integer           , intent(in)           :: xtype            ! netcdf data type
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    integer(i4)           , pointer              :: data(:,:)
    logical           , intent(inout)        :: readvar          ! was var read?
    character(len=*)  , intent(in), optional :: dim1name         ! dimension name
    character(len=*)  , intent(in), optional :: dim2name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)          , intent(in), optional :: missing_value    ! attribute for real
    real(r8)          , intent(in), optional :: fill_value       ! attribute for real
    integer           , intent(in), optional :: imissing_value   ! attribute for int
    integer           , intent(in), optional :: ifill_value      ! attribute for int
    integer           , intent(in), optional :: flag_values(:)   ! attribute for int
    integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int
    !



  end subroutine restartvar_2d_int

  subroutine restartvar_1d_double(&
       ncid, flag, varname, xtype, dim1name, dim2name, &
       long_name, units, interpinic_flag, data, readvar, &
       comment, flag_meanings, missing_value, fill_value, &
       imissing_value, ifill_value, flag_values, nvalid_range )

    !----------------------------------------------------
    ! Arguments
    type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
    character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*)  , intent(in)           :: varname          ! variable name
    integer           , intent(in)           :: xtype            ! netcdf data type
    character(len=*)  , intent(in)           :: long_name        ! long name for variable
    character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    real(r8)           , pointer              :: data(:)
    logical           , intent(inout)        :: readvar          ! was var read?
    character(len=*)  , intent(in), optional :: dim1name         ! dimension name
    character(len=*)  , intent(in), optional :: dim2name         ! dimension name
    character(len=*)  , intent(in), optional :: units            ! long name for variable
    character(len=*)  , intent(in), optional :: comment          ! attribute
    character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)          , intent(in), optional :: missing_value    ! attribute for real
    real(r8)          , intent(in), optional :: fill_value       ! attribute for real
    integer           , intent(in), optional :: imissing_value   ! attribute for int
    integer           , intent(in), optional :: ifill_value      ! attribute for int
    integer           , intent(in), optional :: flag_values(:)   ! attribute for int
    integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int


    end subroutine restartvar_1d_double

    subroutine restartvar_2d_double(&
         ncid, flag, varname, xtype, dim1name, dim2name, &
         long_name, units, interpinic_flag, data, readvar, &
         comment, flag_meanings, missing_value, fill_value, &
         imissing_value, ifill_value, flag_values, nvalid_range )

      !----------------------------------------------------
      ! Arguments
      type(file_desc_t) , intent(inout)        :: ncid             ! netcdf file id
      character(len=*)  , intent(in)           :: flag             ! 'read' or 'write'
      character(len=*)  , intent(in)           :: varname          ! variable name
      integer           , intent(in)           :: xtype            ! netcdf data type
      character(len=*)  , intent(in)           :: long_name        ! long name for variable
      character(len=*)  , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
      real(r8)           , pointer              :: data(:,:)
      logical           , intent(inout)        :: readvar          ! was var read?
      character(len=*)  , intent(in), optional :: dim1name         ! dimension name
      character(len=*)  , intent(in), optional :: dim2name         ! dimension name
      character(len=*)  , intent(in), optional :: units            ! long name for variable
      character(len=*)  , intent(in), optional :: comment          ! attribute
      character(len=*)  , intent(in), optional :: flag_meanings(:) ! attribute
      real(r8)          , intent(in), optional :: missing_value    ! attribute for real
      real(r8)          , intent(in), optional :: fill_value       ! attribute for real
      integer           , intent(in), optional :: imissing_value   ! attribute for int
      integer           , intent(in), optional :: ifill_value      ! attribute for int
      integer           , intent(in), optional :: flag_values(:)   ! attribute for int
      integer           , intent(in), optional :: nvalid_range(2)  ! attribute for int

  end subroutine restartvar_2d_double

  subroutine restartvar_2d_double_bounds(ncid, flag, varname, xtype, &
       dim1name, dim2name, switchdim, lowerb2, upperb2, &
       long_name, units, interpinic_flag, data, readvar, &
       comment, flag_meanings, missing_value, fill_value, &
       imissing_value, ifill_value, flag_values, nvalid_range )

    !----------------------------------------------------
    ! Arguments
    type(file_desc_t), intent(inout)        :: ncid             ! netcdf file id
    character(len=*) , intent(in)           :: flag             ! 'read' or 'write'
    character(len=*) , intent(in)           :: varname          ! variable name
    integer          , intent(in)           :: xtype            ! netcdf data type
    character(len=*) , intent(in)           :: dim1name         ! dimension name
    character(len=*) , intent(in)           :: dim2name         ! dimension name
    logical          , intent(in)           :: switchdim
    character(len=*) , intent(in)           :: long_name        ! long name for variable
    character(len=*) , intent(in)           :: interpinic_flag  ! interpolate variable using interpinic
    real(r8)         , pointer              :: data(:,:)        ! raw data
    logical          , intent(out)          :: readvar          ! was var read?
    integer          , intent(in), optional :: lowerb2
    integer          , intent(in), optional :: upperb2
    character(len=*) , intent(in), optional :: units            ! long name for variable
    character(len=*) , intent(in), optional :: comment          ! attribute
    character(len=*) , intent(in), optional :: flag_meanings(:) ! attribute
    real(r8)         , intent(in), optional :: missing_value    ! attribute for real
    real(r8)         , intent(in), optional :: fill_value       ! attribute for real
    integer          , intent(in), optional :: imissing_value   ! attribute for int
    integer          , intent(in), optional :: ifill_value      ! attribute for int
    integer          , intent(in), optional :: flag_values(:)   ! attribute for int
    integer          , intent(in), optional :: nvalid_range(2)  ! attribute for int

  end subroutine restartvar_2d_double_bounds


  logical function is_restart( )
    ! Determine if restart run
    use clm_varctl, only : nsrest, nsrContinue

    if (nsrest == nsrContinue) then
       is_restart = .true.
    else
       is_restart = .false.
    end if
  end function is_restart

end module restUtilMod
