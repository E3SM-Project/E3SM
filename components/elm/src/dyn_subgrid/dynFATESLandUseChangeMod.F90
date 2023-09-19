module dynFATESLandUseChangeMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the land use harmonization (LUH2) dataset

  ! !USES:
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use spmdMod               , only : masterproc
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use abortutils            , only : endrun
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use elm_varcon            , only : grlnd
  use elm_varctl            , only : iulog

  implicit none

  private

  real(r8), allocatable, public :: landuse_transitions(:,:)
  real(r8), allocatable, public :: landuse_states(:,:)
  real(r8), allocatable, public :: landuse_harvest(:,:)

  integer, public, parameter    :: num_landuse_transition_vars = 108
  integer, public, parameter    :: num_landuse_state_vars = 12
  integer, public, parameter    :: num_landuse_harvest_vars = 5

  ! Define the fates landuse namelist mode switch values
  character(len=13), public, parameter    :: fates_harvest_no_logging   = 'no_harvest'
  character(len=13), public, parameter    :: fates_harvest_logging_only = 'event_code'
  character(len=13), public, parameter    :: fates_harvest_clmlanduse   = 'surfdata_file'
  character(len=13), public, parameter    :: fates_harvest_luh_area     = 'luhdata_area'
  character(len=13), public, parameter    :: fates_harvest_luh_mass     = 'luhdata_mass'

  ! Define landuse harvest unit integer representation
  integer, public, parameter    :: landuse_harvest_area_units = 1
  integer, public, parameter    :: landuse_harvest_mass_units = 2
  integer, public               :: landuse_harvest_units

  type(dyn_file_type), target   :: dynFatesLandUse_file

  ! LUH2 raw wood harvest area fraction
  character(len=10), target :: landuse_harvest_area_varnames(num_landuse_harvest_vars) = &
       [character(len=10) :: 'primf_harv', 'primn_harv', 'secmf_harv', 'secyf_harv', 'secnf_harv']

  ! LUH2 raw wood harvest biomass carbon
  character(len=10), target :: landuse_harvest_mass_varnames(num_landuse_harvest_vars) = &
       [character(len=10) :: 'primf_bioh', 'primn_bioh', 'secmf_bioh', 'secyf_bioh', 'secnf_bioh']

  character(len=10), public, pointer :: landuse_harvest_varnames(:) => null()

  ! Land use name arrays
  character(len=5), public, parameter  :: landuse_state_varnames(num_landuse_state_vars) = &
                    [character(len=5)  :: 'primf', &  ! forested primary land
                                          'primn', &  ! non-forested primary land
                                          'secdf', &  ! potentially forested secondary land
                                          'secdn', &  ! potentially non-forested secondary land
                                          'pastr', &  ! managed pasture
                                          'range', &  ! rangeland
                                          'urban', &  ! urban land
                                          'c3ann', &  ! C3 annual crops
                                          'c4ann', &  ! C4 annual crops
                                          'c3per', &  ! C3 perennial crops
                                          'c4per', &  ! C4 perennial crops
                                          'c3nfx']    ! C3 nitrogen-fixing crops

  character(len=14), public, parameter :: landuse_transition_varnames(num_landuse_transition_vars) = &
                    [character(len=14) :: 'primf_to_secdn','primf_to_pastr','primf_to_range','primf_to_urban', &
                                          'primf_to_c3ann','primf_to_c4ann','primf_to_c3per','primf_to_c4per','primf_to_c3nfx', &
                                          'primn_to_secdf','primn_to_pastr','primn_to_range','primn_to_urban', &
                                          'primn_to_c3ann','primn_to_c4ann','primn_to_c3per','primn_to_c4per','primn_to_c3nfx', &
                                          'secdf_to_secdn','secdf_to_pastr','secdf_to_range','secdf_to_urban', &
                                          'secdf_to_c3ann','secdf_to_c4ann','secdf_to_c3per','secdf_to_c4per','secdf_to_c3nfx', &
                                          'secdn_to_secdf','secdn_to_pastr','secdn_to_range','secdn_to_urban', &
                                          'secdn_to_c3ann','secdn_to_c4ann','secdn_to_c3per','secdn_to_c4per','secdn_to_c3nfx', &
                                          'pastr_to_secdf','pastr_to_secdn','pastr_to_range','pastr_to_urban', &
                                          'pastr_to_c3ann','pastr_to_c4ann','pastr_to_c3per','pastr_to_c4per','pastr_to_c3nfx', &
                                          'range_to_secdf','range_to_secdn','range_to_pastr','range_to_urban', &
                                          'range_to_c3ann','range_to_c4ann','range_to_c3per','range_to_c4per','range_to_c3nfx', &
                                          'urban_to_secdf','urban_to_secdn','urban_to_pastr','urban_to_range', &
                                          'urban_to_c3ann','urban_to_c4ann','urban_to_c3per','urban_to_c4per','urban_to_c3nfx', &
                                          'c3ann_to_c4ann','c3ann_to_c3per','c3ann_to_c4per','c3ann_to_c3nfx', &
                                          'c3ann_to_secdf','c3ann_to_secdn','c3ann_to_pastr','c3ann_to_range','c3ann_to_urban', &
                                          'c4ann_to_c3ann','c4ann_to_c3per','c4ann_to_c4per','c4ann_to_c3nfx', &
                                          'c4ann_to_secdf','c4ann_to_secdn','c4ann_to_pastr','c4ann_to_range','c4ann_to_urban', &
                                          'c3per_to_c3ann','c3per_to_c4ann','c3per_to_c4per','c3per_to_c3nfx', &
                                          'c3per_to_secdf','c3per_to_secdn','c3per_to_pastr','c3per_to_range','c3per_to_urban', &
                                          'c4per_to_c3ann','c4per_to_c4ann','c4per_to_c3per','c4per_to_c3nfx', &
                                          'c4per_to_secdf','c4per_to_secdn','c4per_to_pastr','c4per_to_range','c4per_to_urban', &
                                          'c3nfx_to_c3ann','c3nfx_to_c4ann','c3nfx_to_c3per','c3nfx_to_c4per', &
                                          'c3nfx_to_secdf','c3nfx_to_secdn','c3nfx_to_pastr','c3nfx_to_range','c3nfx_to_urban']

  type(dyn_var_time_uninterp_type) :: landuse_transition_vars(num_landuse_transition_vars) ! value of each transitions variable
  type(dyn_var_time_uninterp_type) :: landuse_state_vars(num_landuse_state_vars)           ! value of each state variable
  type(dyn_var_time_uninterp_type) :: landuse_harvest_vars(num_landuse_harvest_vars)       ! value of each harvest variable

  public :: dynFatesLandUseInit
  public :: dynFatesLandUseInterp

contains

  !-----------------------------------------------------------------------
  subroutine dynFatesLandUseInit(bounds, landuse_filename)

    ! !DESCRIPTION:
    ! Initialize data structures for land use information.

    ! !USES:
    use elm_varctl            , only : use_fates_luh, fates_harvest_mode
    use elm_varctl            , only : use_fates_potentialveg
    use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
    use dynTimeInfoMod        , only : YEAR_POSITION_START_OF_TIMESTEP
    use dynTimeInfoMod        , only : YEAR_POSITION_END_OF_TIMESTEP

    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds                ! proc-level bounds
    character(len=*) , intent(in) :: landuse_filename      ! name of file containing landuse timeseries information (fates luh2)

    ! !LOCAL VARIABLES
    integer :: varnum, i      ! counter for harvest variables
    integer :: landuse_shape(1)  ! land use shape
    integer :: num_points ! number of spatial points
    integer :: ier        ! error code
    real(r8), allocatable :: this_data(:) ! data for a single harvest variable
    !
    character(len=*), parameter :: subname = 'dynFatesLandUseInit'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    ! Allocate and initialize the land use arrays
    allocate(landuse_states(num_landuse_state_vars,bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for landuse_states'//errMsg(__FILE__, __LINE__))
    end if
    allocate(landuse_transitions(num_landuse_transition_vars,bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for landuse_transitions'//errMsg(__FILE__, __LINE__))
    end if
    allocate(landuse_harvest(num_landuse_harvest_vars,bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for landuse_harvest'//errMsg(__FILE__, __LINE__))
    end if

    ! Initialize the states, transitions and harvest mapping percentages as zero by defaut
    landuse_states = 0._r8
    landuse_transitions = 0._r8
    landuse_harvest = 0._r8

    ! Avoid initializing the landuse timeseries file if in fates potential vegetation mode
    if (.not. use_fates_potentialveg) then
       if (use_fates_luh) then

          ! Generate the dyn_file_type object
          ! Start calls get_prev_date, whereas end calls get_curr_date
          dynFatesLandUse_file = dyn_file_type(landuse_filename, YEAR_POSITION_END_OF_TIMESTEP)

          ! Get initial land use data from the fates luh2 timeseries dataset
          num_points = (bounds%endg - bounds%begg + 1)
          landuse_shape(1) = num_points ! Does this need an explicit array shape to be passed to the constructor?
          do varnum = 1, num_landuse_transition_vars
             landuse_transition_vars(varnum) = dyn_var_time_uninterp_type( &
                  dyn_file=dynFatesLandUse_file, varname=landuse_transition_varnames(varnum), &
                  dim1name=grlnd, conversion_factor=1.0_r8, &
                  do_check_sums_equal_1=.false., data_shape=landuse_shape)
          end do
          do varnum = 1, num_landuse_state_vars
             landuse_state_vars(varnum) = dyn_var_time_uninterp_type( &
                  dyn_file=dynFatesLandUse_file, varname=landuse_state_varnames(varnum), &
                  dim1name=grlnd, conversion_factor=1.0_r8, &
                  do_check_sums_equal_1=.false., data_shape=landuse_shape)
          end do

          ! Get the harvest rate data from the fates luh2 timeseries dataset if enabled
          if (trim(fates_harvest_mode) .eq. fates_harvest_luh_area .or. &
               trim(fates_harvest_mode) .eq. fates_harvest_luh_mass) then

             ! change the harvest varnames being used depending on the mode selected
             if (trim(fates_harvest_mode) .eq. fates_harvest_luh_area ) then
                landuse_harvest_varnames => landuse_harvest_area_varnames
                landuse_harvest_units = landuse_harvest_area_units
             elseif (trim(fates_harvest_mode) .eq. fates_harvest_luh_mass ) then
                landuse_harvest_varnames => landuse_harvest_mass_varnames
                landuse_harvest_units = landuse_harvest_mass_units
             else
                call endrun(msg=' undefined fates harvest mode selected'//errMsg(__FILE__, __LINE__))
             end if

             do varnum = 1, num_landuse_harvest_vars
                landuse_harvest_vars(varnum) = dyn_var_time_uninterp_type( &
                     dyn_file=dynFatesLandUse_file, varname=landuse_harvest_varnames(varnum), &
                     dim1name=grlnd, conversion_factor=1.0_r8, &
                     do_check_sums_equal_1=.false., data_shape=landuse_shape)
             end do
          end if

       end if

       ! Since fates needs state data during initialization, make sure to call
       ! the interpolation routine at the start
       call dynFatesLandUseInterp(bounds,init_state=.true.)

    end if

  end subroutine dynFatesLandUseInit


  !-----------------------------------------------------------------------
  subroutine dynFatesLandUseInterp(bounds, init_state)

    use dynTimeInfoMod , only : time_info_type
    use elm_varctl     , only : fates_harvest_mode

    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds       ! proc-level bounds
    logical, optional, intent(in) :: init_state   ! fates needs state for initialization

    ! !LOCAL VARIABLES:
    integer                     :: varnum
    integer                     :: i
    logical                     :: init_flag
    real(r8), allocatable       :: this_data(:)
    character(len=*), parameter :: subname = 'dynFatesLandUseInterp'
    !-----------------------------------------------------------------------
    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    init_flag = .false.
    if (present(init_state)) then
       init_flag = init_state
    end if

    ! Get the data for the current year
    call dynFatesLandUse_file%time_info%set_current_year_get_year()

    if (dynFatesLandUse_file%time_info%is_before_time_series() .and. .not.(init_flag)) then
       ! Reset the land use transitions to zero for safety
       landuse_transitions(1:num_landuse_transition_vars,bounds%begg:bounds%endg) = 0._r8
       landuse_states(1:num_landuse_state_vars,bounds%begg:bounds%endg) = 0._r8
       landuse_harvest(1:num_landuse_harvest_vars,bounds%begg:bounds%endg) = 0._r8
    else
       ! Right now we don't account for the topounits
       allocate(this_data(bounds%begg:bounds%endg))
       do varnum = 1, num_landuse_transition_vars
          call landuse_transition_vars(varnum)%get_current_data(this_data)
          landuse_transitions(varnum,bounds%begg:bounds%endg) = this_data(bounds%begg:bounds%endg)
       end do
       do varnum = 1, num_landuse_state_vars
          call landuse_state_vars(varnum)%get_current_data(this_data)
          landuse_states(varnum,bounds%begg:bounds%endg) = this_data(bounds%begg:bounds%endg)
       end do
       if (trim(fates_harvest_mode) .eq. fates_harvest_luh_area .or. &
           trim(fates_harvest_mode) .eq. fates_harvest_luh_mass) then
          do varnum = 1, num_landuse_harvest_vars
             call landuse_harvest_vars(varnum)%get_current_data(this_data)
             landuse_harvest(varnum,bounds%begg:bounds%endg) = this_data(bounds%begg:bounds%endg)
          end do
       end if
       deallocate(this_data)
    end if

  end subroutine dynFatesLandUseInterp


end module dynFATESLandUseChangeMod
