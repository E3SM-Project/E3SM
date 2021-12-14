!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_main.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM. 
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glad_main

  ! This module provides an interface to GCMs in the case where fields have already been
  ! downscaled to the ice sheet grid (and the GCM does its own upscaling from the ice
  ! sheet grid to the land grid).
  !
  ! This only provides code for the SMB case, not for the PDD case.

  use glimmer_global, only: dp, fname_length
  use glad_type
  use glad_constants
  use glimmer_config
  use glimmer_filenames, only : process_path
  use parallel, only: main_task
  use glad_input_averages, only : get_av_start_time, accumulate_averages, &
       calculate_averages, reset_glad_input_averages
  
  use glimmer_paramets, only: stdout, GLC_DEBUG

  implicit none
  private
  
  ! ------------------------------------------------------------
  ! glad_params derived type definition
  ! This is where default values are set.
  ! ------------------------------------------------------------

  type, public :: glad_params 

     !> Derived type containing parameters relevant to all instances of 
     !> the model - i.e. those parameters which pertain to the global model. 

     ! Ice model instances --------------------------------------

     integer                                   :: ninstances = 1       !> Number of ice model instances
     character(fname_length),pointer,dimension(:) :: config_fnames => null()    ! array of config filenames
     type(glad_instance),pointer,dimension(:) :: instances  => null() !> Array of glimmer\_instances

     ! Global model parameters ----------------------------------

     integer  :: tstep_mbal = 1        !> Mass-balance timestep (hours)
     integer  :: start_time            !> Time of first call to glad (hours)
     integer  :: time_step             !> Calling timestep of global model (hours)

     ! Parameters that can be set by the GCM calling Glad

     logical  :: gcm_restart = .false. !> If true, restart the model from a GCM restart file
     character(fname_length) :: gcm_restart_file   !> Name of restart file
     integer  :: gcm_fileunit = 99     !> Fileunit specified by GCM for reading config files

  end type glad_params

  !---------------------------------------------------------------------------------------
  ! Use of the routines here:
  !
  ! NOTE(wjs, 2015-03-24) I think this is going to need some rework in order to handle
  ! multiple instances the way I'm planning to do it in CESM, with the coupler managing
  ! these multiple instances: I think we're going to want a totally separate glad
  ! instance for each ice sheet instance. Then some of these initialization routines
  ! could be combined.
  !
  ! In model initialization:
  ! - Call glad_initialize once
  ! - Call glad_initialize_instance once per instance
  ! - Call glad_get_grid_size once per instance
  !   (this is needed so that the caller can allocate arrays appropriately)
  ! - Call glad_get_initial_outputs once per instance
  ! - Call glad_initialization_wrapup once
  !
  ! In the model run loop:
  ! - Call glad_gcm once per instance
  !---------------------------------------------------------------------------------------
  
  public :: glad_initialize
  public :: glad_initialize_instance
  public :: glad_get_grid_size
  public :: glad_get_initial_outputs
  public :: glad_initialization_wrapup

  public :: glad_get_grid_indices
  public :: glad_get_lat_lon
  public :: glad_get_areas
  
  public :: glad_gcm

  public :: end_glad
  
  !---------------------------------------------------------------------------------------
  ! Some notes on coupling to the Community Earth System Model (CESM).  These may be applicable
  ! for coupling to other GCMs:
  !
  ! When coupled to CESM, Glad receives two fields from the coupler on the ice sheet grid:
  !   qsmb = surface mass balance (kg/m^2/s)
  !   tsfc = surface ground temperature (deg C)
  ! Both qsmb and tsfc are computed in the CESM land model.
  ! Seven fields are returned to CESM on the ice sheet grid:
  !   ice_covered = whether a grid cell is ice-covered [0,1]
  !   topo = surface elevation (m)
  !   hflx = heat flux from the ice interior to the surface (W/m^2)
  !   rofi = ice runoff (i.e., calving) (kg/m^2/s)
  !   rofl = liquid runoff (i.e., basal melting; the land model handles sfc runoff) (kg/m^2/s)
  !   ice_sheet_grid_mask = mask of ice sheet grid coverage
  !   icemask_coupled_fluxes = mask of ice sheet grid coverage where we are potentially
  !     sending non-zero fluxes
  !
  ! Note about ice_sheet_grid_mask and icemask_coupled_fluxes: ice_sheet_grid_mask is
  ! non-zero wherever CISM is operating - i.e., grid cells with icesheet or bare land (but
  ! not ocean). icemask_coupled_fluxes is similar, but is 0 for icesheet instances that
  ! have zero_gcm_fluxes = .true. Thus, icemask_coupled_fluxes can be used to determine
  ! the regions of the world in which CISM is operating and potentially sending non-zero
  ! fluxes to the climate model.
  !
  ! The land model has the option to update its ice coverage and surface elevation, given
  ! the fields returned from Glad.
  !
  !---------------------------------------------------------------------------------------
  
contains

  subroutine glad_initialize(params, time_step, paramfile, daysinyear, start_time, &
                             gcm_restart, gcm_restart_file, gcm_debug, gcm_fileunit)

    ! Initialize the model for runs coupled to a GCM. This routine initializes variables
    ! shared between instances. See above for documentation of the full initialization
    ! sequence.

    ! Subroutine argument declarations --------------------------------------------------------
    
    type(glad_params),              intent(inout) :: params      !> parameters to be set
    integer,                           intent(in)  :: time_step   !> Timestep of calling model (hours)
    character(*),dimension(:),         intent(in)  :: paramfile   !> array of configuration filenames.
    integer,                  optional,intent(in)  :: daysinyear  !> Number of days in the year
    integer,                  optional,intent(in)  :: start_time  !> Time of first call to glad (hours)
    logical,                  optional,intent(in)  :: gcm_restart ! logical flag to restart from a GCM restart file
    character(*),             optional,intent(in)  :: gcm_restart_file ! restart filename for a GCM restart
                                                                  ! (currently assumed to be CESM)
    logical,                  optional,intent(in)  :: gcm_debug   ! logical flag from GCM to output debug information
    integer,                  optional,intent(in)  :: gcm_fileunit! fileunit for reading config files
    
    ! Internal variables -----------------------------------------------------------------------

    type(ConfigSection), pointer :: global_config
    
    ! Begin subroutine code --------------------------------------------------------------------


    if (present(gcm_debug)) then
       GLC_DEBUG = gcm_debug
    endif

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Initializing glad'
    end if

    ! Initialise start time and calling model time-step (time_step = integer number of hours)
    ! We ignore t=0 by default 

    params%time_step = time_step

    ! Note: start_time = nhour_glad = 0 for an initial run.
    !       Does this create problems given that Glad convention is to ignore t = 0?

    if (present(start_time)) then
       params%start_time = start_time
    else
       params%start_time = time_step
    end if

    params%gcm_restart = .false.
    if (present(gcm_restart)) then
       params%gcm_restart = gcm_restart
    endif

    params%gcm_restart_file = ''
    if (present(gcm_restart_file)) then
       params%gcm_restart_file = gcm_restart_file
    endif

    params%gcm_fileunit = 99
    if (present(gcm_fileunit)) then
       params%gcm_fileunit = gcm_fileunit
    endif

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'time_step     =', params%time_step
       write(stdout,*) 'start_time    =', params%start_time
    end if

    ! Initialise year-length -------------------------------------------------------------------

    if (present(daysinyear)) then
       call glad_set_year_length(daysinyear)
    end if

    ! ---------------------------------------------------------------
    ! Determine how many instances there are, according to what
    ! configuration files we've been provided with
    ! ---------------------------------------------------------------

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Read paramfile'
       write(stdout,*) 'paramfile =', paramfile
    end if

    if (size(paramfile) == 1) then
       ! Load the configuration file into the linked list
       call ConfigRead(process_path(paramfile(1)), global_config, params%gcm_fileunit)    
       ! Parse the list
       call glad_readconfig(global_config, params%ninstances, params%config_fnames, paramfile)
    else
       params%ninstances = size(paramfile)
       allocate(params%config_fnames(params%ninstances))
       params%config_fnames(:) = paramfile(:)
    end if

    allocate(params%instances(params%ninstances))

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Number of instances =', params%ninstances
    end if

  end subroutine glad_initialize

  !===================================================================

  subroutine glad_initialize_instance(params, instance_index)

    ! Initialize one instance in the params structure. See above for documentation of
    ! the full initialization sequence.

    use glad_initialise, only : glad_i_initialise_gcm
    
    ! Subroutine argument declarations --------------------------------------------------------

    type(glad_params),              intent(inout) :: params          !> parameters to be set
    integer,                         intent(in)    :: instance_index  !> index of current ice sheet instance

    ! Internal variables -----------------------------------------------------------------------

    type(ConfigSection), pointer :: instance_config
    
    ! Begin subroutine code --------------------------------------------------------------------

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Read config file and initialize instance #', instance_index
    end if

    call ConfigRead(process_path(params%config_fnames(instance_index)),&
         instance_config, params%gcm_fileunit)

    call glad_i_initialise_gcm(instance_config,     params%instances(instance_index), &
                                params%start_time,   params%time_step,        &
                                params%gcm_restart,  params%gcm_restart_file, &
                                params%gcm_fileunit )

  end subroutine glad_initialize_instance

  !===================================================================

  subroutine glad_get_grid_size(params, instance_index, &
       ewn, nsn, npts, &
       ewn_tot, nsn_tot, npts_tot)

    ! Get the size of a grid corresponding to this instance.
    !
    ! Returns both the size of local arrays (ewn, nsn, npts) and the size of global arrays
    ! (ewn_tot, nsn_tot, npts_tot).
    !
    ! The size is returned withOUT halo cells - note that the other routines here assume
    ! that inputs and outputs do not have halo cells.
    !
    ! The caller can then allocate arrays (inputs to and outputs from glad) with size
    ! (ewn, nsn).

    use parallel, only : own_ewn, own_nsn, global_ewn, global_nsn
    
    type(glad_params), intent(in) :: params
    integer, intent(in) :: instance_index  ! index of current ice sheet instance
    integer, intent(out) :: ewn  ! number of east-west points owned by this proc (first dimension of arrays)
    integer, intent(out) :: nsn  ! number of north-south points owned by this proc (second dimension of arrays)
    integer, intent(out) :: npts ! total number of points owned by this proc
    integer, intent(out) :: ewn_tot ! total number of east-west points in grid
    integer, intent(out) :: nsn_tot ! total number of north-south points in grid
    integer, intent(out) :: npts_tot ! total number of points in grid
    
    ewn = own_ewn
    nsn = own_nsn
    npts = ewn * nsn

    ewn_tot = global_ewn
    nsn_tot = global_nsn
    npts_tot = ewn_tot * nsn_tot

  end subroutine glad_get_grid_size
    
  !===================================================================
  
  subroutine glad_get_initial_outputs(params,         instance_index,        &
                                      ice_covered,    topo,                  &
                                      rofi,           rofl,           hflx,  &
                                      ice_sheet_grid_mask,                   &
                                      icemask_coupled_fluxes,                &
                                      output_flag)

    ! Get initial outputs for one instance. See above for documentation of the full
    ! initialization sequence.
    !
    ! Output arrays are assumed to NOT have halo cells.

    ! Subroutine argument declarations --------------------------------------------------------

    type(glad_params),               intent(in)    :: params
    integer,                         intent(in)    :: instance_index  !> index of current ice sheet instance

    real(dp),dimension(:,:),intent(out) :: ice_covered  ! whether each grid cell is ice-covered [0,1]
    real(dp),dimension(:,:),intent(out) :: topo         ! output surface elevation (m)
    real(dp),dimension(:,:),intent(out) :: hflx         ! output heat flux (W/m^2, positive down)
    real(dp),dimension(:,:),intent(out) :: rofi         ! output ice runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: rofl         ! output liquid runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: ice_sheet_grid_mask !mask of ice sheet grid coverage
    real(dp),dimension(:,:),intent(out) :: icemask_coupled_fluxes !mask of ice sheet grid coverage where we are potentially sending non-zero fluxes
    
    logical,                  optional,intent(out) :: output_flag !> Flag to show output set (provided for consistency)

    ! Begin subroutine code --------------------------------------------------------------------

    call glad_set_output_fields(params%instances(instance_index), &
         ice_covered, topo, rofi, rofl, hflx, &
         ice_sheet_grid_mask, icemask_coupled_fluxes)
    
    if (present(output_flag)) output_flag = .true.

  end subroutine glad_get_initial_outputs
  
  !===================================================================
  
  subroutine glad_initialization_wrapup(params, ice_dt)

    type(glad_params),              intent(inout) :: params      !> parameters to be set
    integer,                  optional,intent(out) :: ice_dt      !> Ice dynamics time-step in hours

    ! Wrapup glad initialization - perform error checks, etc. See above for documentation
    ! of the full initialization sequence

    ! Check that all mass-balance time-steps are the same length and
    ! assign that value to the top-level variable

    params%tstep_mbal = check_mbts(params%instances(:)%mbal_tstep)

    if (present(ice_dt)) then
       ice_dt = check_mbts(params%instances(:)%ice_tstep)
    end if

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'tstep_mbal =', params%tstep_mbal
       write(stdout,*) 'start_time =', params%start_time
       write(stdout,*) 'time_step =',  params%time_step
       if (present(ice_dt)) write(stdout,*) 'ice_dt =', ice_dt
    end if

    ! Check time-steps divide into one another appropriately.

    if (.not. (mod (params%tstep_mbal, params%time_step) == 0)) then
       call write_log('The mass-balance timestep must be an integer multiple of the forcing time-step', &
                       GM_FATAL,__FILE__,__LINE__)
    end if
    

  end subroutine glad_initialization_wrapup

  !===================================================================
  
  subroutine glad_get_grid_indices(params, instance_index, &
                                   global_indices, local_indices)

    ! Get 1-d indices for each grid cell.
    !
    ! The global indices are unique across all tasks (i.e., the global grid). The local
    ! indices go from 1 .. ncells on each task. The global indices increase going from
    ! left to right, and then from bottom to top. So the indices for the bottom
    ! (southernmost) row go 1 .. (# east-west points), etc. The local indices go in the
    ! same order.
    !
    ! The global_indices and local_indices arrays should NOT include halo cells. The
    ! returned indices also ignore halo cells.

    use parallel, only : own_ewn, own_nsn, global_row_offset, global_col_offset, global_ewn
    
    ! Subroutine argument declarations --------------------------------------------------------

    type(glad_params), intent(in) :: params
    integer, intent(in) :: instance_index  ! index of current ice sheet index
    integer, intent(out) :: global_indices(:,:)
    integer, intent(out) :: local_indices(:,:)

    ! Internal variables -----------------------------------------------------------------------

    integer :: own_points  ! number of points this proc is responsible for
    integer, allocatable :: counts(:)  ! count number of times each local index has been set
    integer :: local_row, local_col
    integer :: global_row, global_col
    integer :: local_index, global_index
    character(len=*), parameter :: subname = 'glad_get_grid_indices'

    ! Begin subroutine code --------------------------------------------------------------------
    
    ! Perform error checking on inputs
    
    if (size(global_indices, 1) /= own_ewn .or. size(global_indices, 2) /= own_nsn) then
       call write_log(subname // ' ERROR: Wrong size for global_indices', &
            GM_FATAL, __FILE__, __LINE__)
    end if

    if (size(local_indices, 1) /= own_ewn .or. size(local_indices, 2) /= own_nsn) then
       call write_log(subname // ' ERROR: Wrong size for local_indices', &
            GM_FATAL, __FILE__, __LINE__)
    end if

    ! Set global and local indices

    own_points = own_ewn * own_nsn
    allocate(counts(own_points))
    counts(:) = 0

    do local_row = 1, own_nsn
       do local_col = 1, own_ewn
          local_index = (local_row - 1)*own_ewn + local_col
          if (local_index < 1 .or. local_index > own_points) then
             write(stdout,*) subname//' ERROR: local_index out of bounds: ', &
                  local_index, own_points
             call write_log(subname // ' ERROR: local_index out of bounds', &
                  GM_FATAL, __FILE__, __LINE__)
          end if
          local_indices(local_col,local_row) = local_index
          counts(local_index) = counts(local_index) + 1
          
          global_row = local_row + global_row_offset
          global_col = local_col + global_col_offset
          global_index = (global_row - 1)*global_ewn + global_col
          global_indices(local_col,local_row) = global_index
       end do
    end do
          
    ! Make sure that each local index has been assigned exactly once
    if (any(counts /= 1)) then
       call write_log(subname // ' ERROR: not all local indices have been assigned exactly once', &
            GM_FATAL, __FILE__, __LINE__)
    end if
    
  end subroutine glad_get_grid_indices

  !===================================================================
  
  subroutine glad_get_lat_lon(params, instance_index, &
                              lats, lons)

    ! Get latitude and longitude for each grid cell

    ! Output arrays do NOT have halo cells

    use parallel, only : own_ewn, own_nsn, parallel_convert_haloed_to_nonhaloed
    
    ! Subroutine argument declarations --------------------------------------------------------

    type(glad_params), intent(in) :: params
    integer, intent(in) :: instance_index  ! index of current ice sheet index
    real(dp), intent(out) :: lats(:,:)      ! latitudes (degrees)
    real(dp), intent(out) :: lons(:,:)      ! longitudes (degrees)

    ! Internal variables -----------------------------------------------------------------------
    character(len=*), parameter :: subname = 'glad_get_lat_lon'
    
    ! Begin subroutine code --------------------------------------------------------------------
    
    ! Perform error checking on inputs

    if (size(lats, 1) /= own_ewn .or. size(lats, 2) /= own_nsn) then
       call write_log(subname // ' ERROR: Wrong size for lats', &
            GM_FATAL, __FILE__, __LINE__)
    end if

    if (size(lons, 1) /= own_ewn .or. size(lons, 2) /= own_nsn) then
       call write_log(subname // ' ERROR: Wrong size for lons', &
            GM_FATAL, __FILE__, __LINE__)
    end if

    call parallel_convert_haloed_to_nonhaloed(params%instances(instance_index)%lat, lats)
    call parallel_convert_haloed_to_nonhaloed(params%instances(instance_index)%lon, lons)
    
  end subroutine glad_get_lat_lon

    !===================================================================
  
  subroutine glad_get_areas(params, instance_index, areas)

    ! Get area of each grid cell

    ! Subroutine argument declarations --------------------------------------------------------

    type(glad_params), intent(in) :: params
    integer, intent(in) :: instance_index  ! index of current ice sheet index
    real(dp), intent(out) :: areas(:,:)     ! areas (m^2)

    areas(:,:) = get_dns(params%instances(instance_index)%model) * &
                 get_dew(params%instances(instance_index)%model)
    
  end subroutine glad_get_areas

  
  !===================================================================

  subroutine glad_gcm(params,         instance_index, time,  &
                      qsmb,           tsfc,                  &
                      ice_covered,    topo,                  &
                      rofi,           rofl,           hflx,  &
                      ice_sheet_grid_mask,                   &
                      icemask_coupled_fluxes,                &
                      output_flag,    ice_tstep)

    ! Main Glad subroutine for GCM coupling.
    !
    ! It does all necessary temporal averaging, 
    ! and calls the dynamic ice sheet model when required. 
    !
    ! Input fields should be taken as means over the period since the last call.
    ! See the user documentation for more information.
    !
    ! Input fields are assumed to NOT have halo cells

    use glimmer_utils
    use glad_timestep, only: glad_i_tstep_gcm
    use glimmer_log
    use glimmer_paramets, only: scyr
    use parallel, only : parallel_convert_nonhaloed_to_haloed
    use glide_types, only : get_ewn, get_nsn
    use glad_output_fluxes, only : calculate_average_output_fluxes
    
    implicit none

    ! Subroutine argument declarations -------------------------------------------------------------

    type(glad_params),              intent(inout) :: params          !> parameters for this run
    integer,                         intent(in)    :: instance_index  !> index of current ice sheet instance
    integer,                         intent(in)    :: time            !> Current model time        (hours)

    real(dp),dimension(:,:),intent(in)    :: qsmb          ! input surface mass balance of glacier ice (kg/m^2/s)
    real(dp),dimension(:,:),intent(in)    :: tsfc          ! input surface ground temperature (deg C)

    real(dp),dimension(:,:),intent(inout) :: ice_covered  ! whether each grid cell is ice-covered [0,1]
    real(dp),dimension(:,:),intent(inout) :: topo         ! output surface elevation (m)
    real(dp),dimension(:,:),intent(inout) :: hflx         ! output heat flux (W/m^2, positive down)
    real(dp),dimension(:,:),intent(inout) :: rofi         ! output ice runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(inout) :: rofl         ! output liquid runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(inout) :: ice_sheet_grid_mask !mask of ice sheet grid coverage
    real(dp),dimension(:,:),intent(inout) :: icemask_coupled_fluxes !mask of ice sheet grid coverage where we are potentially sending non-zero fluxes

    logical,optional,intent(out)   :: output_flag     ! Set true if outputs are set
    logical,optional,intent(out)   :: ice_tstep       ! Set when an ice dynamic timestep has been done
                                                      !  and new output is available

    ! Internal variables ----------------------------------------------------------------------------

    integer :: ewn,nsn    ! dimensions of local grid

    ! version of input fields with halo cells
    real(dp),dimension(:,:),allocatable :: qsmb_haloed
    real(dp),dimension(:,:),allocatable :: tsfc_haloed

    logical :: icets
    character(250) :: message

    integer :: av_start_time  ! value of time from the last occasion averaging was restarted (hours)

    ! Begin subroutine code --------------------------------------------------------------------

    ! Reset output flag

    if (present(output_flag)) output_flag = .false.
    if (present(ice_tstep))   ice_tstep = .false.

    ! Accumulate input fields for later averaging

    ewn = get_ewn(params%instances(instance_index)%model)
    nsn = get_nsn(params%instances(instance_index)%model)
    allocate(qsmb_haloed(ewn,nsn))
    allocate(tsfc_haloed(ewn,nsn))
    call parallel_convert_nonhaloed_to_haloed(qsmb, qsmb_haloed)
    call parallel_convert_nonhaloed_to_haloed(tsfc, tsfc_haloed)

    call accumulate_averages(params%instances(instance_index)%glad_inputs, &
         qsmb = qsmb_haloed, tsfc = tsfc_haloed, time = time)

    ! ---------------------------------------------------------
    ! If this is a mass balance timestep, prepare global fields, and do a timestep
    ! for each model instance
    ! ---------------------------------------------------------

    av_start_time = get_av_start_time(params%instances(instance_index)%glad_inputs)
    
    if (mod (time - av_start_time, params%time_step) /= 0) then
       
       write(message,*) 'Unexpected calling of GLAD at time ', time
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    
    else if (time - av_start_time + params%time_step > params%tstep_mbal) then

       write(message,*) &
            'Incomplete forcing of GLAD mass-balance time-step detected at time ', time
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
       
    else if (time - av_start_time + params%time_step == params%tstep_mbal) then

       ! Set output_flag

       ! At present, outputs are done for each mass-balance timestep, since
       ! that involved least change to the code. However, it might be good
       ! to change the output to occur with user-specified frequency.

       if (present(output_flag)) output_flag = .true.

       ! Do a timestep for this instance

       if (time == params%instances(instance_index)%next_time) then

          params%instances(instance_index)%next_time = &
               params%instances(instance_index)%next_time + &
               params%instances(instance_index)%mbal_tstep

          ! Calculate averages by dividing by number of steps elapsed
          ! since last model timestep.

          call calculate_averages(params%instances(instance_index)%glad_inputs, &
               qsmb = params%instances(instance_index)%acab, &
               tsfc = params%instances(instance_index)%artm)

          ! Calculate total surface mass balance - multiply by time since last model timestep
          ! Note on units: We want acab to have units of meters w.e. (accumulated over mass balance time step)
          ! Initial units are kg m-2 s-1 = mm s-1
          ! Divide by 1000 to convert from mm to m
          ! Multiply by hours2seconds = 3600 to convert from 1/s to 1/hr.  (tstep_mbal has units of hours)

          !TODO - Modify code so that qsmb and acab are always in kg m-2 s-1 water equivalent?
          params%instances(instance_index)%acab(:,:) = &
               params%instances(instance_index)%acab(:,:) * &
               params%tstep_mbal * hours2seconds / 1000.d0

          if (GLC_DEBUG .and. main_task) write(stdout,*) 'Take a glad time step, instance', instance_index
          call glad_i_tstep_gcm(time,                  &
               params%instances(instance_index),   &
               icets)

          call calculate_average_output_fluxes( &
               params%instances(instance_index)%glad_output_fluxes, &
               rofi_tavg = params%instances(instance_index)%rofi_tavg, &
               rofl_tavg = params%instances(instance_index)%rofl_tavg, &
               hflx_tavg = params%instances(instance_index)%hflx_tavg)

          call glad_set_output_fields(params%instances(instance_index), &
               ice_covered, topo, rofi, rofl, hflx, &
               ice_sheet_grid_mask, icemask_coupled_fluxes)
          

          ! Set flag
          if (present(ice_tstep)) then
             ice_tstep = (ice_tstep .or. icets)
          end if

       endif   ! time = next_time

       ! ---------------------------------------------------------
       ! Reset averaging fields, flags and counters
       ! ---------------------------------------------------------

       call reset_glad_input_averages(params%instances(instance_index)%glad_inputs, &
            next_av_start = time + params%time_step)
       
       if (GLC_DEBUG .and. main_task) then
          write(stdout,*) 'Done in glad_gcm'
       endif

   endif    ! time - av_start_time + params%time_step > params%tstep_mbal

  end subroutine glad_gcm

  !===================================================================

  subroutine end_glad(params,close_logfile)

    !> tidy-up operations for Glad
    use glad_initialise
    use glimmer_log
    implicit none

    type(glad_params),intent(inout) :: params          ! parameters for this run
    logical, intent(in), optional    :: close_logfile   ! if true, then close the log file
                                                        ! (GCM may do this elsewhere)                                  
    integer :: i

    ! end individual instances

    do i = 1, params%ninstances
       call glad_i_end(params%instances(i))
    enddo

    if (present(close_logfile)) then
       if (close_logfile) call close_log
    else
       call close_log
    endif

    deallocate(params%config_fnames)
    deallocate(params%instances)

  end subroutine end_glad

  !----------------------------------------------------------------------
  ! PRIVATE INTERNAL GLIMMER SUBROUTINES FOLLOW.............
  !----------------------------------------------------------------------

  subroutine glad_set_output_fields(instance,        &
                                    ice_covered,    topo,                  &
                                    rofi,           rofl,           hflx,  &
                                    ice_sheet_grid_mask,                   &
                                    icemask_coupled_fluxes)

    ! Sets output fields for this instance.
    !
    ! Arguments are assumed to NOT have halo cells. This routine handles the removal of
    ! the halo cells.

    use glad_output_states, only : set_output_states
    use parallel, only : parallel_convert_haloed_to_nonhaloed
    use glide_types, only : get_ewn, get_nsn

    ! Subroutine argument declarations --------------------------------------------------------

    type(glad_instance), intent(in) :: instance

    real(dp),dimension(:,:),intent(out) :: ice_covered  ! whether each grid cell is ice-covered [0,1]
    real(dp),dimension(:,:),intent(out) :: topo         ! output surface elevation (m)
    real(dp),dimension(:,:),intent(out) :: hflx         ! output heat flux (W/m^2, positive down)
    real(dp),dimension(:,:),intent(out) :: rofi         ! output ice runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: rofl         ! output liquid runoff (kg/m^2/s = mm H2O/s)
    real(dp),dimension(:,:),intent(out) :: ice_sheet_grid_mask !mask of ice sheet grid coverage
    real(dp),dimension(:,:),intent(out) :: icemask_coupled_fluxes !mask of ice sheet grid coverage where we are potentially sending non-zero fluxes
    
    ! Internal variables -----------------------------------------------------------------------

    integer :: ewn,nsn    ! dimensions of local grid
    
    ! temporary versions of output fields with halo cells
    real(dp),dimension(:,:),allocatable :: ice_covered_haloed
    real(dp),dimension(:,:),allocatable :: topo_haloed
    real(dp),dimension(:,:),allocatable :: hflx_haloed
    real(dp),dimension(:,:),allocatable :: rofi_haloed
    real(dp),dimension(:,:),allocatable :: rofl_haloed
    real(dp),dimension(:,:),allocatable :: ice_sheet_grid_mask_haloed
    real(dp),dimension(:,:),allocatable :: icemask_coupled_fluxes_haloed

    ! Begin subroutine code --------------------------------------------------------------------

    ewn = get_ewn(instance%model)
    nsn = get_nsn(instance%model)

    allocate(ice_covered_haloed(ewn,nsn))
    allocate(topo_haloed(ewn,nsn))
    allocate(hflx_haloed(ewn,nsn))
    allocate(rofi_haloed(ewn,nsn))
    allocate(rofl_haloed(ewn,nsn))
    allocate(ice_sheet_grid_mask_haloed(ewn,nsn))
    allocate(icemask_coupled_fluxes_haloed(ewn,nsn))
    
    call set_output_states(instance, &
         ice_covered_haloed, topo_haloed, ice_sheet_grid_mask_haloed)

    if (instance%zero_gcm_fluxes == ZERO_GCM_FLUXES_TRUE) then
       icemask_coupled_fluxes_haloed(:,:) = 0.d0
       hflx_haloed(:,:) = 0.d0
       rofi_haloed(:,:) = 0.d0
       rofl_haloed(:,:) = 0.d0
    else
       icemask_coupled_fluxes_haloed(:,:) = ice_sheet_grid_mask_haloed(:,:)
       hflx_haloed(:,:) = instance%hflx_tavg(:,:)
       rofi_haloed(:,:) = instance%rofi_tavg(:,:)
       rofl_haloed(:,:) = instance%rofl_tavg(:,:)
    end if

    call parallel_convert_haloed_to_nonhaloed(ice_covered_haloed, ice_covered)
    call parallel_convert_haloed_to_nonhaloed(topo_haloed, topo)
    call parallel_convert_haloed_to_nonhaloed(hflx_haloed, hflx)
    call parallel_convert_haloed_to_nonhaloed(rofi_haloed, rofi)
    call parallel_convert_haloed_to_nonhaloed(rofl_haloed, rofl)
    call parallel_convert_haloed_to_nonhaloed(ice_sheet_grid_mask_haloed, ice_sheet_grid_mask)
    call parallel_convert_haloed_to_nonhaloed(icemask_coupled_fluxes_haloed, icemask_coupled_fluxes)

  end subroutine glad_set_output_fields
  
  !TODO - Move subroutine glad_readconfig to a glad_setup module, in analogy to glide_setup?

  subroutine glad_readconfig(config, ninstances, fnames, infnames)

    !> Determine whether a given config file is a
    !> top-level glad config file, and return parameters
    !> accordingly.

    use glimmer_config
    use glimmer_log
    implicit none

    ! Arguments -------------------------------------------

    type(ConfigSection),      pointer :: config !> structure holding sections of configuration file
    integer,              intent(out) :: ninstances !> Number of instances to create
    character(fname_length),dimension(:),pointer :: fnames !> list of filenames (output)
    character(fname_length),dimension(:) :: infnames !> list of filenames (input)

    ! Internal variables ----------------------------------

    type(ConfigSection), pointer :: section
    character(len=100) :: message
    integer :: i

    if (associated(fnames)) nullify(fnames)

    call GetSection(config,section,'GLAD')
    if (associated(section)) then
       call GetValue(section,'n_instance',ninstances)
       allocate(fnames(ninstances))
       do i=1,ninstances
          call GetSection(section%next,section,'GLAD instance')
          if (.not.associated(section)) then
             write(message,*) 'Must specify ',ninstances,' instance config files'
             call write_log(message,GM_FATAL,__FILE__,__LINE__)
          end if
          call GetValue(section,'name',fnames(i))
       end do
    else
       ninstances=1
       allocate(fnames(1))
       fnames=infnames
    end if

    ! Print some configuration information

!!$    call write_log('GLAD global')
!!$    call write_log('------------')
!!$    write(message,*) 'number of instances :',params%ninstances
!!$    call write_log(message)
!!$    call write_log('')

  end subroutine glad_readconfig


  !========================================================

  integer function check_mbts(timesteps)

    !> Checks to see that all mass-balance time-steps are
    !> the same. Flags a fatal error if not, else assigns that
    !> value to the output

    use glimmer_log

    implicit none

    integer,dimension(:) :: timesteps !> Array of mass-balance timsteps

    integer :: n,i

    n = size(timesteps)
    if (n==0) then
       check_mbts = 0
       return
    endif

    check_mbts = timesteps(1)

    do i = 2,n
       if (timesteps(i) /= check_mbts) then
          call write_log('All instances must have the same mass-balance and ice timesteps', &
               GM_FATAL,__FILE__,__LINE__)
       endif
    enddo

  end function check_mbts

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glad_main

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
