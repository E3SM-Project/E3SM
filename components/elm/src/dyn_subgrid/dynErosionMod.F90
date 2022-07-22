module dynErosionMod 

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Handle reading of the dataset that specifies time-varying ELM-Erosion 
  ! parameters.
  
  ! !USES:
#include "shr_assert.h"
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use elm_varctl            , only : iulog
  use elm_varcon            , only : grlnd
  use spmdMod               , only : masterproc
  use ColumnType            , only : col_pp
  use topounit_varcon       , only : max_topounits
  use GridcellType          , only : grc_pp
  
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save
  public :: dynErosion_init     ! initialize information read from landuse.timeseries dataset
  public :: dynErosion_interp   ! get parameter data for the current time step, if needed
  
  ! ! PRIVATE TYPES
  type(dyn_file_type), target      :: dynErosion_file ! information for the file containing parameter data
  type(dyn_var_time_uninterp_type) :: pfactor      ! rainfall-driven erosion scaling factor
  type(dyn_var_time_uninterp_type) :: qfactor      ! runoff-driven erosion scaling factor  
  type(dyn_var_time_uninterp_type) :: tillage      ! conserved tillage fraction 

  ! Names of variables on file
  character(len=*), parameter :: c1_varname = 'parEro_c1'
  character(len=*), parameter :: c2_varname = 'parEro_c2'
  character(len=*), parameter :: tillage_varname = 'Tillage'

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !---------------------------------------------------------------------------

contains
  !-----------------------------------------------------------------------
  subroutine dynErosion_init(bounds, dynErosion_filename)
    
    ! !DESCRIPTION:
    ! Initialize dataset containing ELM-Erosion parameters info (position it to the right time
    ! samples that bound the initial model date)
    
    ! !USES:
    use elm_varctl     , only : use_erosion 
    use dynTimeInfoMod , only : YEAR_POSITION_START_OF_TIMESTEP
    
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds           ! proc-level bounds
    character(len=*)  , intent(in) :: dynErosion_filename ! name of file containing parameter information
    
    ! !LOCAL VARIABLES:
    integer :: num_points       ! number of spatial points
    integer :: c1_shape(2)      ! Shape of the parEro_c1 data   
    integer :: c2_shape(2)      ! shape of the parEro_c2 data      
    integer :: tillage_shape(2) ! shape of the Tillage data   
    character(len=*), parameter :: subname = 'dynErosion_init'
    !-----------------------------------------------------------------------
    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    if (masterproc) then
       write(iulog,*) 'Attempting to read ELM-Erosion dynamic parameter data .....'
    end if

    dynErosion_file = dyn_file_type(dynErosion_filename, YEAR_POSITION_START_OF_TIMESTEP)
    
    ! Get initial parameter data 
    if (use_erosion) then
       num_points = (bounds%endg - bounds%begg + 1)     
       ! parEro_c1
       c1_shape = [num_points, max_topounits]       
       pfactor = dyn_var_time_uninterp_type( &
            dyn_file = dynErosion_file, varname=c1_varname, &
            dim1name=grlnd, conversion_factor=1._r8, &
            do_check_sums_equal_1=.false., data_shape=c1_shape)         
       ! parEro_c2
       c2_shape = [num_points, max_topounits]                   
       qfactor = dyn_var_time_uninterp_type( &
            dyn_file = dynErosion_file, varname=c2_varname, &
            dim1name=grlnd, conversion_factor=1._r8, &
            do_check_sums_equal_1=.false., data_shape=c2_shape)
       ! Tillage
       tillage_shape = [num_points, max_topounits]                     
       tillage = dyn_var_time_uninterp_type( &
            dyn_file = dynErosion_file, varname=tillage_varname, &
            dim1name=grlnd, conversion_factor=1._r8, &
            do_check_sums_equal_1=.false., data_shape=tillage_shape)
    end if

  end subroutine dynErosion_init

  !-----------------------------------------------------------------------
  subroutine dynErosion_interp(bounds, soilstate_vars, sedflux_vars)
    
    ! !DESCRIPTION:
    ! Get ELM-Erosion parameters for model time, when needed.
    
    ! Note that ELM-Erosion data are stored as values (not weights) and so time interpolation
    ! is not necessary - the parameter value is held constant through the year.  This is
    ! consistent with the treatment of changing land use, where interpolation of the
    ! annual endpoint weights leads to a constant rate of change in PFT weight through the
    ! year, with abrupt changes in the rate at annual boundaries.

    ! !USES:
    use elm_varctl        , only : use_erosion
    use SedFluxType       , only : sedflux_type 
    use SoilStateType     , only : soilstate_type 
    
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    type(soilstate_type), intent(inout) :: soilstate_vars
    type(sedflux_type), intent(inout) :: sedflux_vars 
    
    ! !LOCAL VARIABLES:
    integer               :: c,g,t,ti,topi      ! indices  
    real(r8), allocatable :: pfactor_cur(:,:) ! current pfactor
    real(r8), allocatable :: qfactor_cur(:,:) ! current qfactor 
    real(r8), allocatable :: tillage_cur(:,:) ! current conserved tillage fraction
    character(len=*), parameter :: subname = 'dynErosion_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    call dynErosion_file%time_info%set_current_year()

    ! Set new parameters
    if (use_erosion) then
       allocate(pfactor_cur(bounds%begg:bounds%endg,max_topounits))    
       allocate(qfactor_cur(bounds%begg:bounds%endg,max_topounits))
       allocate(tillage_cur(bounds%begg:bounds%endg,max_topounits))

       call pfactor%get_current_data(pfactor_cur)
       call qfactor%get_current_data(qfactor_cur)
       call tillage%get_current_data(tillage_cur)
       do c = bounds%begc, bounds%endc
          g = col_pp%gridcell(c)
          t = col_pp%topounit(c)
          topi = grc_pp%topi(g)
          ti = t - topi + 1
          
          soilstate_vars%tillage_col(c) = tillage_cur(g,ti)
          sedflux_vars%pfactor_col(c)   = pfactor_cur(g,ti)
          sedflux_vars%qfactor_col(c)   = qfactor_cur(g,ti)
       end do

       deallocate(pfactor_cur)
       deallocate(qfactor_cur)
       deallocate(tillage_cur)
    end if
  end subroutine dynErosion_interp

end module dynErosionMod 
