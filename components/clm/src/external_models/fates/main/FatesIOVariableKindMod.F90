module FatesIOVariableKindMod

  use FatesConstantsMod, only : fates_long_string_length
  use FatesGlobals, only : fates_log
  use FatesIODimensionsMod, only : fates_io_dimension_type

  implicit none

  ! FIXME(bja, 2016-10) do these need to be strings, or can they be integer enumerations?
  ! FIXME(rgk, 2016-11) these should probably be moved to varkindmod?
  
  character(*), parameter :: patch_r8 = 'PA_R8'
  character(*), parameter :: patch_ground_r8 = 'PA_GRND_R8'
  character(*), parameter :: patch_size_pft_r8 = 'PA_SCPF_R8'
  character(*), parameter :: site_r8 = 'SI_R8'
  character(*), parameter :: site_int = 'SI_INT'
  character(*), parameter :: site_ground_r8 = 'SI_GRND_R8'
  character(*), parameter :: site_size_pft_r8 = 'SI_SCPF_R8'
  character(*), parameter :: site_size_r8 = 'SI_SCLS_R8'
  character(*), parameter :: patch_int = 'PA_INT'
  character(*), parameter :: cohort_r8 = 'CO_R8'
  character(*), parameter :: cohort_int = 'CO_INT'
  character(*), parameter :: site_pft_r8 = 'SI_PFT_R8'
  character(*), parameter :: site_age_r8 = 'SI_AGE_R8'
  character(*), parameter :: site_fuel_r8 = 'SI_FUEL_R8'
  character(*), parameter :: site_cwdsc_r8 = 'SI_CWDSC_R8'
  character(*), parameter :: site_can_r8 = 'SI_CAN_R8'
  character(*), parameter :: site_cnlf_r8 = 'SI_CNLF_R8'
  character(*), parameter :: site_cnlfpft_r8 = 'SI_CNLFPFT_R8'
  character(*), parameter :: site_scag_r8 = 'SI_SCAG_R8'

  ! NOTE(RGK, 2016) %active is not used yet. Was intended as a check on the HLM->FATES
  ! control parameter passing to ensure all active dimension types received all
  ! dimensioning specifications from the host, but we currently arent using those
  ! passing functions..

  ! This structure is not multi-threaded
  type fates_io_variable_kind_type
     character(len=fates_long_string_length) :: name ! String labelling this IO type
     integer              :: ndims       ! number of dimensions in this IO type
     integer, allocatable :: dimsize(:)  ! The size of each dimension
     logical, private :: active_
     integer :: dim1_index
     integer :: dim2_index

   contains

     procedure, public :: Init
     procedure, public :: set_active
     procedure, public :: is_active

  end type fates_io_variable_kind_type



contains

  ! ===================================================================================
  subroutine Init(this, name, num_dims)

    use FatesConstantsMod, only : fates_unset_int
    
    implicit none

    class(fates_io_variable_kind_type), intent(inout) :: this
    character(*), intent(in) :: name
    integer, intent(in) :: num_dims
    
    this%name = trim(name)
    this%ndims = num_dims
    allocate(this%dimsize(this%ndims))
    this%dimsize(:) = fates_unset_int
    this%active_ = .false.
    this%dim1_index = fates_unset_int
    this%dim2_index = fates_unset_int
    
  end subroutine Init
  
 ! =======================================================================
 subroutine set_active(this)
   implicit none
   class(fates_io_variable_kind_type), intent(inout) :: this
   this%active_ = .true.
 end subroutine set_active

 logical function is_active(this)
   implicit none
   class(fates_io_variable_kind_type), intent(in) :: this
   is_active = this%active_
 end function is_active

  ! ====================================================================================

  function iotype_index(iotype_name, num_dim_kinds, dim_kinds) result(dk_index)

    ! argument
    character(len=*), intent(in) :: iotype_name
    integer, intent(in) :: num_dim_kinds
    type(fates_io_variable_kind_type), intent(in) :: dim_kinds(:)

    ! local
    integer :: dk_index

    do dk_index=1, num_dim_kinds
       if (trim(iotype_name) .eq. trim(dim_kinds(dk_index)%name)) then
          return
       end if
    end do
    write(fates_log(),*) 'An IOTYPE THAT DOESNT EXIST WAS SPECIFIED'
    !end_run

  end function iotype_index
  
end module FatesIOVariableKindMod
