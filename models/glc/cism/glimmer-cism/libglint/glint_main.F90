!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_main.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glint_main

  !*FD  This is the main glimmer module, which contains the top-level 
  !*FD  subroutines and derived types comprising the glimmer ice model.

  use glimmer_global
  use glint_type
  use glint_global_grid
  use glint_constants
  use glimmer_anomcouple

  use glimmer_paramets, only: stdout, GLC_DEBUG

  implicit none

  ! ------------------------------------------------------------
  ! GLIMMER_PARAMS derived type definition
  ! This is where default values are set.
  ! ------------------------------------------------------------

  !TODO - Move this definition to glint_type?
  !       Create a simpler type (glint_params_gcm) for GCM coupling?

  type glint_params 

     !*FD Derived type containing parameters relevant to all instances of 
     !*FD the model - i.e. those parameters which pertain to the global model. 

     ! Global grids used ----------------------------------------

     type(global_grid) :: g_grid      !*FD The main global grid, used for 
                                      !*FD input and most outputs
     type(global_grid) :: g_grid_orog !*FD Global grid used for orography output.

     ! Ice model instances --------------------------------------

     integer                                   :: ninstances = 1       !*FD Number of ice model instances
     type(glint_instance),pointer,dimension(:) :: instances  => null() !*FD Array of glimmer\_instances

     ! Global model parameters ----------------------------------

     integer  :: tstep_mbal = 1        !*FD Mass-balance timestep (hours)
     integer  :: start_time            !*FD Time of first call to glint (hours)
     integer  :: time_step             !*FD Calling timestep of global model (hours)

     ! Parameters that can be set by the GCM calling Glint

     logical  :: gcm_smb = .false.     !*FD If true, receive surface mass balance from the GCM 
     logical  :: gcm_restart = .false. !*FD If true, restart the model from a GCM restart file
     character(fname_length) :: gcm_restart_file   !*FD Name of restart file
     integer  :: gcm_fileunit = 99     !*FD Fileunit specified by GCM for reading config files
   
     ! Averaging parameters -------------------------------------

     integer  :: av_start_time = 0   !*FD Holds the value of time from 
                                     !*FD the last occasion averaging was restarted (hours)
     integer  :: av_steps      = 0   !*FD Holds the number of times glimmer has 
                                     !*FD been called in current round of averaging.
     integer  :: next_av_start = 0   !*FD Time when we expect next averaging to start
     logical  :: new_av     = .true. !*FD Set to true if the next correct call starts a new averaging round

     ! Averaging arrays -----------------------------------------

     real(rk),pointer,dimension(:,:) :: g_av_precip  => null()  !*FD globally averaged precip
     real(rk),pointer,dimension(:,:) :: g_av_temp    => null()  !*FD globally averaged temperature 
     real(rk),pointer,dimension(:,:) :: g_max_temp   => null()  !*FD global maximum temperature
     real(rk),pointer,dimension(:,:) :: g_min_temp   => null()  !*FD global minimum temperature
     real(rk),pointer,dimension(:,:) :: g_temp_range => null()  !*FD global temperature range
     real(rk),pointer,dimension(:,:) :: g_av_zonwind => null()  !*FD globally averaged zonal wind 
     real(rk),pointer,dimension(:,:) :: g_av_merwind => null()  !*FD globally averaged meridional wind 
     real(rk),pointer,dimension(:,:) :: g_av_humid   => null()  !*FD globally averaged humidity (%)
     real(rk),pointer,dimension(:,:) :: g_av_lwdown  => null()  !*FD globally averaged downwelling longwave (W/m$^2$)
     real(rk),pointer,dimension(:,:) :: g_av_swdown  => null()  !*FD globally averaged downwelling shortwave (W/m$^2$)
     real(rk),pointer,dimension(:,:) :: g_av_airpress => null() !*FD globally averaged surface air pressure (Pa)
     real(rk),pointer,dimension(:,:,:) :: g_av_qsmb => null()   ! globally averaged surface mass balance (kg m-2 s-1)
     real(rk),pointer,dimension(:,:,:) :: g_av_tsfc => null()   ! globally averaged surface temperature (deg C)
     real(rk),pointer,dimension(:,:,:) :: g_av_topo => null()   ! globally averaged surface elevation   (m)

     ! Fractional coverage information --------------------------
     ! Note: these are only valid on the main task
     real(rk),pointer,dimension(:,:) :: total_coverage  => null()     !*FD Fractional coverage by 
                                                                      !*FD all ice model instances.
     real(rk),pointer,dimension(:,:) :: cov_normalise   => null()     !*FD Normalisation values 
                                                                      !*FD for coverage calculation.
     real(rk),pointer,dimension(:,:) :: total_cov_orog  => null()     !*FD Fractional coverage by 
                                                                      !*FD all ice model instances (orog).
     real(rk),pointer,dimension(:,:) :: cov_norm_orog   => null()     !*FD Normalisation values 
                                                                      !*FD for coverage calculation (orog).
     logical                         :: coverage_calculated = .false. !*FD Have we calculated the
                                                                      !*FD coverage map yet?
     ! File information -----------------------------------------

     character(fname_length) :: paramfile      !*FD Name of global parameter file

     ! Accumulation/averaging flags -----------------------------

     logical :: need_winds=.false. !*FD Set if we need the winds to be accumulated/downscaled
     logical :: enmabal=.false.    !*FD Set if we're using the energy balance mass balance model anywhere

     ! Anomaly coupling for global climate ------------------------------------------

     type(anomaly_coupling) :: anomaly_params !*FD Parameters for anomaly coupling

  end type glint_params

  ! Private names -----------------------------------------------

  private glint_allocate_arrays
  private glint_readconfig, calc_bounds, check_init_args

    !---------------------------------------------------------------------------------------
    ! Some notes on coupling to the Community Earth System Model (CESM).  These may be applicable
    ! for coupling to other GCMs:
    !
    ! When coupled to CESM, Glint receives three fields from the coupler on a global grid
    ! in each of several elevation classes:
    !   qsmb = surface mass balance (kg/m^2/s)
    !   tsfc = surface ground temperature (deg C)
    !   topo = surface elevation (m)
    ! Both qsmb and tsfc are computed in the CESM land model.
    ! Five fields are returned to CESM on the global grid for each elevation class:
    !   gfrac = fractional ice coverage
    !   gtopo = surface elevation (m)
    !   grofi = ice runoff (i.e., calving) (kg/m^2/s)
    !   grofl = liquid runoff (i.e., basal melting; the land model handles sfc runoff) (kg/m^2/s)
    !   ghflx = heat flux from the ice interior to the surface (W/m^2)
    ! The land model has the option to update its ice coverage and surface elevation, given
    ! the fields returned from Glint.
    !
    ! There are two driver subroutines in this module for CESM coupling: 
    !  initialise_glint_gcm (for initialization) and glint_gcm (for timestepping).
    ! These drivers loop over the ice sheet model instances (just Greenland for now,
    !  but will simulate Antarctica and other ice sheets later).
    !  
    ! The other driver subroutines, based on the original Glint code in Glimmer version 1,
    !  are initialise_glint and glint.
    ! These subroutines are usually run with temp (= 2-m air temperature) and precip as input.
    !  The surface mass balance is computed using a daily or annual PDD scheme.
    ! It would be possible to call these subroutines from CESM and use the PDD scheme,
    !  but this option has not been tested.
    !
    !---------------------------------------------------------------------------------------

contains

  !TODO - Try calling this subroutine from CESM to estimate SMB in PDD mode?
  !       We would have only one elevation class per grid cell and would not return upscaled fields.

  subroutine initialise_glint(params,                         &
                              lats,         longs,            &
                              time_step,    paramfile,        &
                              latb,         lonb,             &
                              orog,         albedo,           &
                              ice_frac,     veg_frac,         &
                              snowice_frac, snowveg_frac,     &
                              snow_depth,                     &
                              orog_lats,    orog_longs,       &
                              orog_latb,    orog_lonb,        &
                              output_flag,  daysinyear,       &
                              snow_model,   ice_dt,           &
                              extraconfigs, start_time,       &
                              gmask,                          &
                              gcm_restart,  gcm_restart_file, &
                              gcm_debug,    gcm_fileunit)

    !*FD Initialises the model
    !*FD For a multi-processor run, the main task should specify lats & longs spanning
    !*FD the full global domain; the other tasks should give 0-size lats & longs arrays
    !*FD Output arrays on the global grid are only valid on the main task

    use glimmer_config
    use glint_initialise
    use glimmer_log
    use glimmer_filenames
    use glint_timestep, only: get_i_upscaled_fields
    use parallel, only: main_task
    implicit none

    ! Subroutine argument declarations --------------------------------------------------------

    type(glint_params),              intent(inout) :: params      !*FD parameters to be set
    real(rk),dimension(:),           intent(in)    :: lats,longs  !*FD location of gridpoints 
                                                                  !*FD in global data.
    integer,                         intent(in)    :: time_step   !*FD Timestep of calling model (hours)
    character(*),dimension(:),       intent(in)    :: paramfile   !*FD array of configuration filenames.
    real(rk),dimension(:),  optional,intent(in)    :: latb        !*FD Locations of the latitudinal 
                                                                  !*FD boundaries of the grid-boxes.
    real(rk),dimension(:),  optional,intent(in)    :: lonb        !*FD Locations of the longitudinal
                                                                  !*FD boundaries of the grid-boxes.
    real(rk),dimension(:,:),optional,intent(out)   :: orog        !*FD Initial global orography
    real(rk),dimension(:,:),optional,intent(out)   :: albedo      !*FD Initial albedo
    real(rk),dimension(:,:),optional,intent(out)   :: ice_frac    !*FD Initial ice fraction 
    real(rk),dimension(:,:),optional,intent(out)   :: veg_frac    !*FD Initial veg fraction
    real(rk),dimension(:,:),optional,intent(out)   :: snowice_frac !*FD Initial snow-covered ice fraction
    real(rk),dimension(:,:),optional,intent(out)   :: snowveg_frac !*FD Initial snow-covered veg fraction
    real(rk),dimension(:,:),optional,intent(out)   :: snow_depth  !*FD Initial snow depth 
    real(rk),dimension(:),  optional,intent(in)    :: orog_lats   !*FD Latitudinal location of gridpoints 
                                                                  !*FD for global orography output.
    real(rk),dimension(:),  optional,intent(in)    :: orog_longs  !*FD Longitudinal location of gridpoints 
                                                                  !*FD for global orography output.
    real(rk),dimension(:),  optional,intent(in)    :: orog_latb   !*FD Locations of the latitudinal 
                                                                  !*FD boundaries of the grid-boxes (orography).
    real(rk),dimension(:),  optional,intent(in)    :: orog_lonb   !*FD Locations of the longitudinal
                                                                  !*FD boundaries of the grid-boxes (orography).
    logical,                optional,intent(out)   :: output_flag !*FD Flag to show output set (provided for
                                                                  !*FD consistency)
    integer,                optional,intent(in)    :: daysinyear  !*FD Number of days in the year
    logical,                optional,intent(out)   :: snow_model  !*FD Set if the mass-balance scheme has a snow-depth model
    integer,                optional,intent(out)   :: ice_dt      !*FD Ice dynamics time-step in hours
    type(ConfigData),dimension(:),optional ::  extraconfigs !*FD Additional configuration information - overwrites
                                                                  !*FD config data read from files
    integer,                optional,intent(in)    :: start_time  !*FD Time of first call to glint (hours)
    integer, dimension(:,:),  optional,intent(in)  :: gmask       !*FD mask = 1 where global data are valid
    logical,                  optional,intent(in)  :: gcm_restart ! logical flag to restart from a GCM restart file
    character(*),             optional, intent(in) :: gcm_restart_file ! restart filename for a GCM restart
                                                                  ! (currently assumed to be CESM)
    logical,                  optional,intent(in)  :: gcm_debug   ! logical flag from GCM to output debug information
    integer,                  optional,intent(in)  :: gcm_fileunit! fileunit for reading config files

    ! Internal variables -----------------------------------------------------------------------

    type(ConfigSection), pointer :: global_config, instance_config, section  ! configuration stuff
    character(len=100) :: message                                            ! For log-writing
    character(fname_length),dimension(:),pointer :: config_fnames=>null()    ! array of config filenames
    type(ConfigSection), pointer :: econf
    integer :: i, j, n
    real(rk),dimension(:,:),allocatable :: orog_temp, if_temp, vf_temp, sif_temp,  &
                                           svf_temp,  sd_temp, alb_temp      ! Temporary output arrays
    integer,dimension(:),allocatable :: mbts,idts ! Array of mass-balance and ice dynamics timesteps
    logical :: anomaly_check ! Set if we've already initialised anomaly coupling

    if (present(gcm_debug)) then
       GLC_DEBUG = gcm_debug
    endif

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Starting initialise_glint'
    end if

    ! Initialise start time and calling model time-step ----------------------------------------
    ! We ignore t=0 by default 

    params%time_step = time_step

    if (present(start_time)) then
       params%start_time = start_time
    else
       params%start_time = time_step
    end if

    params%next_av_start = params%start_time

    ! Initialisation for runs coupled to a GCM
 
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

!    nec = 1
!    if (present(gcm_nec)) then
!       nec = gcm_nec
!    endif

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'time_step =', params%time_step
       write(stdout,*) 'start_time =', params%start_time
       write(stdout,*) 'next_av_start =', params%next_av_start
    end if

    ! Initialise year-length -------------------------------------------------------------------

    if (present(daysinyear)) then
       call glint_set_year_length(daysinyear)
    end if

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Initialize global grid'
       write(stdout,*) 'present =', present(gmask)
    end if

    ! Initialise main global grid --------------------------------------------------------------

    if (present(gmask)) then
       call new_global_grid(params%g_grid, longs, lats, lonb=lonb, latb=latb, mask=gmask)
    else
       call new_global_grid(params%g_grid, longs, lats, lonb=lonb, latb=latb)
    endif

    if (GLC_DEBUG .and. main_task) then
       write (stdout,*) ' ' 
       write (stdout,*) 'time_step (hr)  =', params%time_step
       write (stdout,*) 'start_time (hr) =', params%start_time
       write (stdout,*) 'Called new_global_grid '
       write (stdout,*) 'g_grid%nx =', params%g_grid%nx
       write (stdout,*) 'g_grid%ny =', params%g_grid%ny
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lons =', params%g_grid%lons
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lats =', params%g_grid%lats
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lon_bound =', params%g_grid%lon_bound
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lat_bound =', params%g_grid%lat_bound
       do j = 5, 10
          write (stdout,*)
          write (stdout,*) 'j, g_grid%mask =', j, params%g_grid%mask(:,j)
       enddo
    end if

    ! Initialise orography grid ------------------------------------

    call check_init_args(orog_lats, orog_longs, orog_latb, orog_lonb)

    if (present(orog_lats) .and. present(orog_longs)) then
       call new_global_grid(params%g_grid_orog, orog_longs, orog_lats,  &
                            lonb=orog_lonb, latb=orog_latb)
    else
       call copy_global_grid(params%g_grid, params%g_grid_orog)
    end if

    ! Allocate arrays -----------------------------------------------

    call glint_allocate_arrays(params)

    ! Initialise arrays ---------------------------------------------

    params%g_av_precip  = 0.0
    params%g_av_temp    = 0.0
    params%g_max_temp   = -1000.0
    params%g_min_temp   = 1000.0
    params%g_temp_range = 0.0
    params%g_av_zonwind = 0.0
    params%g_av_merwind = 0.0
    params%g_av_humid   = 0.0
    params%g_av_lwdown  = 0.0
    params%g_av_swdown  = 0.0
    params%g_av_airpress = 0.0

    ! ---------------------------------------------------------------
    ! Zero coverage maps and normalisation fields for main grid and
    ! orography grid
    ! ---------------------------------------------------------------

    params%total_coverage = 0.0
    params%total_cov_orog = 0.0

    params%cov_normalise = 0.0
    params%cov_norm_orog = 0.0

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Read paramfile'
       write(stdout,*) 'paramfile =', paramfile
    end if

    ! ---------------------------------------------------------------
    ! Determine how many instances there are, according to what
    ! configuration files we've been provided with
    ! ---------------------------------------------------------------

    if (size(paramfile) == 1) then
       ! Load the configuration file into the linked list
       call ConfigRead(process_path(paramfile(1)), global_config, params%gcm_fileunit)    
       call glint_readconfig(global_config, params%ninstances, config_fnames, paramfile) ! Parse the list
    else
       params%ninstances = size(paramfile)
       allocate(config_fnames(params%ninstances))
       config_fnames = paramfile
    end if

    allocate(params%instances(params%ninstances))
    allocate(mbts(params%ninstances), idts(params%ninstances))

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Number of instances =', params%ninstances
       write(stdout,*) 'Read config files and initialize each instance'
    end if

    ! ---------------------------------------------------------------
    ! Read config files, and initialise instances accordingly
    ! ---------------------------------------------------------------

    call write_log('Reading instance configurations')
    call write_log('-------------------------------')

    anomaly_check = .false.

    do i = 1, params%ninstances

       call ConfigRead(process_path(config_fnames(i)),instance_config, params%gcm_fileunit)
       if (present(extraconfigs)) then
          if (size(extraconfigs)>=i) then
             call ConfigCombine(instance_config,extraconfigs(i))
          end if
       end if

       call glint_i_initialise(instance_config,    params%instances(i),     &
                               params%g_grid,      params%g_grid_orog,      &
                               mbts(i),            idts(i),                 &
                               params%need_winds,  params%enmabal,          &
                               params%start_time,  params%time_step,        &
                               params%gcm_restart, params%gcm_restart_file, &
                               params%gcm_fileunit )

       params%total_coverage = params%total_coverage + params%instances(i)%frac_coverage
       params%total_cov_orog = params%total_cov_orog + params%instances(i)%frac_cov_orog

       where (params%total_coverage > 0.0) params%cov_normalise = params%cov_normalise + 1.0
       where (params%total_cov_orog > 0.0) params%cov_norm_orog = params%cov_norm_orog + 1.0

       ! Initialise anomaly coupling
       if (.not.anomaly_check) then 
          call anomaly_init(params%anomaly_params, instance_config)
          if (params%anomaly_params%enabled .and. &
               (params%anomaly_params%nx/=params%g_grid%nx .or. &
                params%anomaly_params%ny/=params%g_grid%ny) ) then
             call write_log("Anomaly coupling grids have different "// &
                  "sizes to GLINT coupling grids",GM_FATAL,__FILE__,__LINE__)
          end if
          if (params%anomaly_params%enabled) anomaly_check=.true.
       end if

    end do

    ! Check that all mass-balance time-steps are the same length and 
    ! assign that value to the top-level variable

    params%tstep_mbal = check_mbts(mbts)
    if (present(ice_dt)) then
       ice_dt = check_mbts(idts)
    end if

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'tstep_mbal =', params%tstep_mbal
       write(stdout,*) 'start_time =', params%start_time
       write(stdout,*) 'time_step =',  params%time_step
       if (present(ice_dt)) write(stdout,*) 'ice_dt =', ice_dt
    end if

    ! Check time-steps divide into one another appropriately.

    if (.not.(mod(params%tstep_mbal,params%time_step)==0)) then
       print*,params%tstep_mbal,params%time_step
       call write_log('The mass-balance timestep must be an integer multiple of the forcing time-step', &
            GM_FATAL,__FILE__,__LINE__)
    end if

    ! Check we don't have coverage greater than one at any point.

    where (params%total_coverage > 1.0) params%total_coverage = 1.0
    where (params%total_cov_orog > 1.0) params%total_cov_orog = 1.0
    params%coverage_calculated=.true.

    ! Zero optional outputs, if present

    if (present(orog)) orog = 0.0
    if (present(albedo)) albedo = 0.0
    if (present(ice_frac)) ice_frac = 0.0
    if (present(veg_frac)) veg_frac = 0.0
    if (present(snowice_frac)) snowice_frac = 0.0
    if (present(snowveg_frac)) snowveg_frac = 0.0
    if (present(snow_depth)) snow_depth = 0.0

    ! Allocate arrays

    allocate(orog_temp(params%g_grid_orog%nx, params%g_grid_orog%ny))
    allocate(alb_temp (params%g_grid%nx, params%g_grid%ny))
    allocate(if_temp  (params%g_grid%nx, params%g_grid%ny))
    allocate(vf_temp  (params%g_grid%nx, params%g_grid%ny))
    allocate(sif_temp (params%g_grid%nx, params%g_grid%ny))
    allocate(svf_temp (params%g_grid%nx, params%g_grid%ny))
    allocate(sd_temp  (params%g_grid%nx, params%g_grid%ny))

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Upscale and splice the initial fields'
    end if

    ! Get initial fields from instances, splice together and return

    do i=1,params%ninstances

       call get_i_upscaled_fields(params%instances(i),   &
                                  orog_temp,   alb_temp, &
                                  if_temp,     vf_temp,  &
                                  sif_temp,    svf_temp, &
                                  sd_temp)

       ! Add this contribution to the global output
       ! Only the main task has valid values for the global output fields
       if (main_task) then
          
          if (present(orog)) &
               orog = splice_field(orog, orog_temp, params%instances(i)%frac_cov_orog, &
               params%cov_norm_orog)

          if (present(albedo)) &
               albedo = splice_field(albedo, alb_temp, params%instances(i)%frac_coverage, &
               params%cov_normalise)

          if (present(ice_frac)) &
               ice_frac = splice_field(ice_frac, if_temp, params%instances(i)%frac_coverage, &
               params%cov_normalise)

          if (present(veg_frac)) &
               veg_frac = splice_field(veg_frac, vf_temp, params%instances(i)%frac_coverage, &
               params%cov_normalise)

          if (present(snowice_frac)) &
               snowice_frac = splice_field(snowice_frac,sif_temp,params%instances(i)%frac_coverage, &
               params%cov_normalise)

          if (present(snowveg_frac)) &
               snowveg_frac = splice_field(snowveg_frac,svf_temp,params%instances(i)%frac_coverage, &
               params%cov_normalise)

          if (present(snow_depth)) &
               snow_depth = splice_field(snow_depth,sd_temp,params%instances(i)%frac_coverage, &
               params%cov_normalise)

       end if

    end do  ! ninstances

    ! Deallocate

    deallocate(orog_temp, alb_temp, if_temp, vf_temp, sif_temp, svf_temp,sd_temp)

    ! Sort out snow_model flag

    if (present(snow_model)) then
       snow_model = .false.
       do i=1, params%ninstances
          snow_model = (snow_model .or. glint_has_snow_model(params%instances(i)))
       end do
    end if

    ! Set output flag

    if (present(output_flag)) output_flag = .true.

  end subroutine initialise_glint

  !================================================================================

  subroutine initialise_glint_gcm(params,                         &
                                  lats,         longs,            &
                                  time_step,    paramfile,        &
                                  daysinyear,   start_time,       &
                                  ice_dt,       output_flag,      &
                                  glc_nec,                        &
                                  gfrac,        gtopo,            &
                                  grofi,        grofl,            &
                                  ghflx,        gmask,            &
                                  gcm_restart,  gcm_restart_file, &
                                  gcm_debug,    gcm_fileunit)

    ! Initialise the model
    ! For a multi-processor run, the main task should specify lats & longs spanning
    !  the full global domain; the other tasks should give 0-size lats & longs arrays
    ! Output arrays on the global grid are only valid on the main task

    use glimmer_config
    use glint_initialise
    use glimmer_log
    use glimmer_filenames
    use glint_timestep, only: get_i_upscaled_fields_gcm
    use parallel, only: main_task

    implicit none

    ! Subroutine argument declarations --------------------------------------------------------

    type(glint_params),              intent(inout) :: params      !*FD parameters to be set
    real(rk),dimension(:),           intent(in)    :: lats,longs  !*FD location of gridpoints 
                                                                  !*FD in global data.
    integer,                         intent(in)    :: time_step   !*FD Timestep of calling model (hours)
    character(*),dimension(:),       intent(in)    :: paramfile   !*FD array of configuration filenames.
    integer,                optional,intent(in)    :: daysinyear  !*FD Number of days in the year
    integer,                optional,intent(in)    :: start_time  !*FD Time of first call to glint (hours)
    integer,                optional,intent(out)   :: ice_dt      !*FD Ice dynamics time-step in hours
    logical,                optional,intent(out)   :: output_flag !*FD Flag to show output set (provided for consistency)
    integer,                  optional,intent(in)  :: glc_nec     !*FD number of elevation classes for GCM input
    real(rk),dimension(:,:,:),optional,intent(out) :: gfrac       !*FD ice fractional area [0,1]
    real(rk),dimension(:,:,:),optional,intent(out) :: gtopo       !*FD surface elevation (m)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofi       !*FD ice runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofl       !*FD liquid runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: ghflx       !*FD heat flux (W/m^2, positive down)
    integer, dimension(:,:),  optional,intent(in)  :: gmask       !*FD mask = 1 where global data are valid
    logical,                  optional,intent(in)  :: gcm_restart ! logical flag to restart from a GCM restart file
    character(*),             optional, intent(in) :: gcm_restart_file ! restart filename for a GCM restart
                                                                  ! (currently assumed to be CESM)
    logical,                  optional,intent(in)  :: gcm_debug   ! logical flag from GCM to output debug information
    integer,                  optional,intent(in)  :: gcm_fileunit! fileunit for reading config files

    ! Internal variables -----------------------------------------------------------------------

    type(ConfigSection), pointer :: global_config, instance_config, section  ! configuration stuff

    character(len=100) :: message                                            ! For log-writing
    character(fname_length),dimension(:),pointer :: config_fnames=>null()    ! array of config filenames

    type(ConfigSection), pointer :: econf

    integer :: i

    integer,dimension(:),allocatable :: mbts,idts ! Array of mass-balance and ice dynamics timesteps

    real(rk),dimension(:,:,:),allocatable ::   &
               gfrac_temp, gtopo_temp, grofi_temp, grofl_temp, ghflx_temp    ! Temporary output arrays

    integer :: n
    integer :: nec       ! number of elevation classes
    integer :: j, ii, jj

    if (present(gcm_debug)) then
       GLC_DEBUG = gcm_debug
    endif

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Initializing glint'
    end if

    ! Initialise start time and calling model time-step (time_step = integer number of hours)
    ! We ignore t=0 by default 

    params%time_step = time_step

    ! Note: start_time = nhour_glint = 0 for an initial run.
    !       Does this create problems given that Glint convention is to ignore t = 0?

    if (present(start_time)) then
       params%start_time = start_time
    else
       params%start_time = time_step
    end if

    params%next_av_start = params%start_time

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

    nec = 1
    if (present(glc_nec)) then
       nec = glc_nec
    endif

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'time_step     =', params%time_step
       write(stdout,*) 'start_time    =', params%start_time
       write(stdout,*) 'next_av_start =', params%next_av_start
    end if

    ! Initialise year-length -------------------------------------------------------------------

    if (present(daysinyear)) then
       call glint_set_year_length(daysinyear)
    end if

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Initialize global grid: present(gmask) =', present(gmask)
    end if

    ! Initialise main global grid --------------------------------------------------------------

    !TODO - Will gmask always be present for GCM runs?

    if (present(gmask)) then
       call new_global_grid(params%g_grid, longs, lats, nec=nec, mask=gmask)
    else
       call new_global_grid(params%g_grid, longs, lats, nec=nec)
    endif

    if (GLC_DEBUG .and. main_task) then
       write (stdout,*) ' ' 
       write (stdout,*) 'time_step (hr)  =', params%time_step
       write (stdout,*) 'start_time (hr) =', params%start_time
       write (stdout,*) 'Called new_global_grid '
       write (stdout,*) 'g_grid%nx =', params%g_grid%nx
       write (stdout,*) 'g_grid%ny =', params%g_grid%ny
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lons =', params%g_grid%lons
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lats =', params%g_grid%lats
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lon_bound =', params%g_grid%lon_bound
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lat_bound =', params%g_grid%lat_bound
       do j = 5, 10
          write (stdout,*)
          write (stdout,*) 'j, g_grid%mask =', j, params%g_grid%mask(:,j)
       enddo
    end if

    ! Allocate arrays -----------------------------------------------

    allocate(params%total_coverage(params%g_grid%nx, params%g_grid%ny))
    allocate(params%cov_normalise (params%g_grid%nx, params%g_grid%ny))

    allocate(params%g_av_qsmb (params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
    allocate(params%g_av_tsfc (params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
    allocate(params%g_av_topo (params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))

    ! Initialise arrays ---------------------------------------------

    !TODO - Make these arrays dp
    params%g_av_qsmb    = 0.0
    params%g_av_tsfc    = 0.0
    params%g_av_topo    = 0.0

    ! ---------------------------------------------------------------
    ! Zero coverage maps and normalisation fields for main grid
    ! (Not using orography grid for GCM coupling)
    ! ---------------------------------------------------------------

    params%total_coverage = 0.0
    params%cov_normalise = 0.0

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Read paramfile'
       write(stdout,*) 'paramfile =', paramfile
    end if

    ! ---------------------------------------------------------------
    ! Determine how many instances there are, according to what
    ! configuration files we've been provided with
    ! ---------------------------------------------------------------

    if (size(paramfile) == 1) then
       ! Load the configuration file into the linked list
       call ConfigRead(process_path(paramfile(1)), global_config, params%gcm_fileunit)    
       call glint_readconfig(global_config, params%ninstances, config_fnames, paramfile) ! Parse the list
    else
       params%ninstances = size(paramfile)
       allocate(config_fnames(params%ninstances))
       config_fnames = paramfile
    end if

    allocate(params%instances(params%ninstances))

    allocate(mbts(params%ninstances), idts(params%ninstances))

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Number of instances =', params%ninstances
       write(stdout,*) 'Read config files and initialize each instance'
    end if

    ! ---------------------------------------------------------------
    ! Read config files, and initialise instances accordingly
    ! ---------------------------------------------------------------

    call write_log('Reading instance configurations')
    call write_log('-------------------------------')

    do i=1,params%ninstances

       call ConfigRead(process_path(config_fnames(i)),instance_config, params%gcm_fileunit)

         !WHL - I don't think this will be needed; commented out for now
!!       if (present(extraconfigs)) then
!!          if (size(extraconfigs)>=i) then
!!             call ConfigCombine(instance_config,extraconfigs(i))
!!          end if
!!       end if
 
       call glint_i_initialise_gcm(instance_config,     params%instances(i),     &
                                   params%g_grid,                                &
                                   mbts(i),             idts(i),                 &
                                   params%start_time,   params%time_step,        &
                                   params%gcm_restart,  params%gcm_restart_file, &
                                   params%gcm_fileunit )

       params%total_coverage = params%total_coverage + params%instances(i)%frac_coverage
       where (params%total_coverage > 0.0) params%cov_normalise = params%cov_normalise + 1.0

    end do    ! ninstances

    ! Check that all mass-balance time-steps are the same length and 
    ! assign that value to the top-level variable

    params%tstep_mbal = check_mbts(mbts)

    if (present(ice_dt)) then
       ice_dt = check_mbts(idts)
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

    ! Make sure we don't have coverage greater than one at any point.

    where (params%total_coverage > 1.0) params%total_coverage = 1.0
    params%coverage_calculated = .true.

    ! Zero optional outputs, if present

    !TODO - Change to dp
    if (present(gfrac)) gfrac = 0.0
    if (present(gtopo)) gtopo = 0.0
    if (present(grofi)) grofi = 0.0
    if (present(grofl)) grofl = 0.0
    if (present(ghflx)) ghflx = 0.0

    ! Allocate arrays

    allocate(gfrac_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
    allocate(gtopo_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
    allocate(grofi_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
    allocate(grofl_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
    allocate(ghflx_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Upscale and splice the initial fields'
    end if

    ! Get initial fields from instances, splice together and return

    do i = 1, params%ninstances

       !TODO - These *_temp arrays are not yet upscaled correctly

       ! Upscale the output fields for this instance

       call get_i_upscaled_fields_gcm(params%instances(i), params%g_grid%nec,  &
                                      params%instances(i)%lgrid%size%pt(1),    &
                                      params%instances(i)%lgrid%size%pt(2),    &
                                      params%g_grid%nx,    params%g_grid%ny,   &
                                      gfrac_temp,          gtopo_temp,         &
                                      grofi_temp,          grofl_temp,         &
                                      ghflx_temp)

       ! Splice together with the global output

       call splice_fields_gcm(gfrac_temp, gtopo_temp,    &
                              grofi_temp, grofl_temp,    &
                              ghflx_temp,                &
                              gfrac,      gtopo,         &
                              grofi,      grofl,         &
                              ghflx,                     &
                              params%g_grid%nec,         &
                              params%instances(i)%frac_coverage, &
                              params%cov_normalise)

    end do       ! ninstances

    ! Deallocate

    deallocate(gfrac_temp, gtopo_temp, grofi_temp, grofl_temp, ghflx_temp)

    ! Set output flag       !TODO - Is this ever used?
    if (present(output_flag)) output_flag = .true.

  end subroutine initialise_glint_gcm

  !================================================================================

  subroutine glint(params,         time,            &
                   rawtemp,        rawprecip,       &
                   orog,                            &
                   zonwind,        merwind,         &
                   humid,          lwdown,          &
                   swdown,         airpress,        &
                   output_flag,                     &
                   orog_out,       albedo,          &
                   ice_frac,       veg_frac,        &
                   snowice_frac,   snowveg_frac,    &
                   snow_depth,                      &
                   water_in,       water_out,       &
                   total_water_in, total_water_out, &
                   ice_volume,     ice_tstep)

    !*FD Main Glint subroutine.
    !*FD
    !*FD This should be called daily or hourly, depending on
    !*FD the mass-balance scheme being used. It does all necessary 
    !*FD spatial and temporal averaging, and calls the dynamics 
    !*FD part of the model when required. 
    !*FD
    !*FD Input fields should be taken as means over the period since the last call.
    !*FD See the user documentation for more information.
    !*FD
    !*FD Global output fields are only valid on the main task. Fields that are integrated
    !*FD over the whole domain (total_water_in, total_water_out, ice_volume) are only
    !*FD valid in single-task runs; trying to compute these in multi-task runs will generate a
    !*FD fatal error.
    !*FD
    !*FD Note that the total ice volume returned is the total at the end of the time-step;
    !*FD the water fluxes are valid over the duration of the timestep. Thus the difference
    !*FD between \texttt{total\_water\_in} and \texttt{total\_water\_out} should be equal
    !*FD to the change in \texttt{ice\_volume}, after conversion between m$^3$ and kg.

    use glimmer_utils
    use glint_interp
    use glint_timestep
    use glimmer_log
    use glimmer_paramets, only: scyr
    use parallel, only: main_task, tasks
    implicit none

    ! Subroutine argument declarations -------------------------------------------------------------

    type(glint_params),              intent(inout) :: params          !*FD parameters for this run
    integer,                         intent(in)    :: time            !*FD Current model time        (hours)
    real(rk),dimension(:,:),target,  intent(in)    :: rawtemp         !*FD Surface temperature field (deg C)
    real(rk),dimension(:,:),target,  intent(in)    :: rawprecip       !*FD Precipitation rate        (mm/s)
    real(rk),dimension(:,:),         intent(in)    :: orog            !*FD The large-scale orography (m)
    real(rk),dimension(:,:),optional,intent(in)    :: zonwind,merwind !*FD Zonal and meridional components 
                                                                      !*FD of the wind field         (m/s)
    real(rk),dimension(:,:),optional,intent(in)    :: humid           !*FD Surface humidity (%)
    real(rk),dimension(:,:),optional,intent(in)    :: lwdown          !*FD Downwelling longwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: swdown          !*FD Downwelling shortwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: airpress        !*FD surface air pressure (Pa)
    logical,                optional,intent(out)   :: output_flag     !*FD Set true if outputs set
    real(rk),dimension(:,:),optional,intent(inout) :: orog_out        !*FD The fed-back, output orography (m)
    real(rk),dimension(:,:),optional,intent(inout) :: albedo          !*FD surface albedo
    real(rk),dimension(:,:),optional,intent(inout) :: ice_frac        !*FD grid-box ice-fraction
    real(rk),dimension(:,:),optional,intent(inout) :: veg_frac        !*FD grid-box veg-fraction
    real(rk),dimension(:,:),optional,intent(inout) :: snowice_frac    !*FD grid-box snow-covered ice fraction
    real(rk),dimension(:,:),optional,intent(inout) :: snowveg_frac    !*FD grid-box snow-covered veg fraction
    real(rk),dimension(:,:),optional,intent(inout) :: snow_depth      !*FD grid-box mean snow depth (m water equivalent)
    real(rk),dimension(:,:),optional,intent(inout) :: water_in        !*FD Input water flux          (mm)
    real(rk),dimension(:,:),optional,intent(inout) :: water_out       !*FD Output water flux         (mm)
    real(rk),               optional,intent(inout) :: total_water_in  !*FD Area-integrated water flux in (kg)
    real(rk),               optional,intent(inout) :: total_water_out !*FD Area-integrated water flux out (kg)
    real(rk),               optional,intent(inout) :: ice_volume      !*FD Total ice volume (m$^3$)
    logical,                optional,intent(out)   :: ice_tstep       !*FD Set when an ice-timestep has been done, and
                                                                      !*FD water balance information is available

    ! Internal variables ----------------------------------------------------------------------------

    integer :: i, n
    real(rk),dimension(:,:),allocatable :: albedo_temp, if_temp, vf_temp, sif_temp, svf_temp,  &
                                           sd_temp, wout_temp, orog_out_temp, win_temp
    real(rk) :: twin_temp,twout_temp,icevol_temp
    type(output_flags) :: out_f
    logical :: icets
    character(250) :: message
    real(rk),dimension(size(rawprecip,1),size(rawprecip,2)),target :: anomprecip
    real(rk),dimension(size(rawtemp,1),  size(rawtemp,2)),  target :: anomtemp
    real(rk),dimension(:,:),pointer :: precip
    real(rk),dimension(:,:),pointer :: temp
    real(rk) :: yearfrac
    integer :: j, ig, jg

    if (GLC_DEBUG .and. main_task) then
!       write (stdout,*) 'In subroutine glint, current time (hr) =', time
!       write (stdout,*) 'av_start_time =', params%av_start_time
!       write (stdout,*) 'next_av_start =', params%next_av_start
!       write (stdout,*) 'new_av =', params%new_av
!       write (stdout,*) 'tstep_mbal =', params%tstep_mbal
    end if

    ! Check we're expecting a call now --------------------------------------------------------------

    if (params%new_av) then
       if (time == params%next_av_start) then
          params%av_start_time = time
          params%new_av = .false.
       else
          write(message,*) 'Unexpected calling of GLINT at time ', time
          call write_log(message,GM_FATAL,__FILE__,__LINE__)
       end if
    else
       if (mod(time-params%av_start_time,params%time_step) /= 0) then
          write(message,*) 'Unexpected calling of GLINT at time ', time
          call write_log(message,GM_FATAL,__FILE__,__LINE__)
       end if
    end if

    ! Check input fields are correct ----------------------------------------------------------------

    call check_input_fields(params, humid, lwdown, swdown, airpress, zonwind, merwind)

    ! Reset output flag

    if (present(output_flag)) output_flag = .false.
    if (present(ice_tstep))   ice_tstep = .false.

    ! Sort out anomaly coupling

    if (params%anomaly_params%enabled) then
       yearfrac = real(mod(time,days_in_year),rk)/real(days_in_year,rk)
       call anomaly_calc(params%anomaly_params, yearfrac, rawtemp, rawprecip, anomtemp, anomprecip)
       precip => anomprecip
       temp   => anomtemp
    else
       precip => rawprecip
       temp   => rawtemp
    end if

    ! Do averaging and so on...

    call accumulate_averages(params,           &
                             temp,    precip,  &
                             zonwind, merwind, &
                             humid,   lwdown,  &
                             swdown,  airpress)

    ! Increment step counter

    params%av_steps = params%av_steps + 1

    ! ---------------------------------------------------------
    ! If this is a mass balance timestep, prepare global fields, and do a timestep
    ! for each model instance
    ! ---------------------------------------------------------

    if (time - params%av_start_time + params%time_step > params%tstep_mbal) then

       write(message,*) &
            'Incomplete forcing of GLINT mass-balance time-step detected at time ', time
       call write_log(message,GM_FATAL,__FILE__,__LINE__)

    else if (time - params%av_start_time + params%time_step == params%tstep_mbal) then

       ! Set output_flag

       ! At present, outputs are done for each mass-balance timestep, since
       ! that involved least change to the code. However, it might be good
       ! to change the output to occur with user-specified frequency.

       if (present(output_flag)) output_flag = .true.

       ! Allocate output fields

       if (present(orog_out)) then
          allocate(orog_out_temp(size(orog_out,1),size(orog_out,2)))
       else
          allocate(orog_out_temp(params%g_grid_orog%nx, params%g_grid_orog%ny))
       end if
       allocate(albedo_temp(size(orog,1),size(orog,2)))
       allocate(if_temp(size(orog,1),size(orog,2)))
       allocate(vf_temp(size(orog,1),size(orog,2)))
       allocate(sif_temp(size(orog,1),size(orog,2)))
       allocate(svf_temp(size(orog,1),size(orog,2)))
       allocate(sd_temp(size(orog,1),size(orog,2)))
       allocate(wout_temp(size(orog,1),size(orog,2)))
       allocate(win_temp(size(orog,1),size(orog,2)))

       ! Populate output flag derived type

       call populate_output_flags(out_f,                           &
                                  orog_out,       albedo,          &
                                  ice_frac,       veg_frac,        &
                                  snowice_frac,   snowveg_frac,    &
                                  snow_depth,                      &
                                  water_in,       water_out,       &
                                  total_water_in, total_water_out, &
                                  ice_volume)

       ! Zero outputs if present

       if (present(orog_out))        orog_out        = 0.0
       if (present(albedo))          albedo          = 0.0
       if (present(ice_frac))        ice_frac        = 0.0
       if (present(veg_frac))        veg_frac        = 0.0
       if (present(snowice_frac))    snowice_frac    = 0.0
       if (present(snowveg_frac))    snowveg_frac    = 0.0
       if (present(snow_depth))      snow_depth      = 0.0
       if (present(water_out))       water_out       = 0.0
       if (present(water_in))        water_in        = 0.0
       if (present(total_water_in))  total_water_in  = 0.0
       if (present(total_water_out)) total_water_out = 0.0
       if (present(ice_volume))      ice_volume      = 0.0

       ! Calculate averages by dividing by number of steps elapsed
       ! since last model timestep.

       call calculate_averages(params)

       ! Calculate total accumulated precipitation - multiply
       ! by time since last model timestep

       params%g_av_precip = params%g_av_precip*params%tstep_mbal*hours2seconds

       ! Calculate temperature half-range

       params%g_temp_range = (params%g_max_temp-params%g_min_temp)/2.0

       if (GLC_DEBUG .and. main_task) then
          i = params%instances(1)%model%numerics%idiag_global
          j = params%instances(1)%model%numerics%jdiag_global
          write(stdout,*) ' '
          write(stdout,*) 'Take a mass balance timestep, time (hr) =', time
          write(stdout,*) 'av_steps =', real(params%av_steps,rk)
          write(stdout,*) 'tstep_mbal (hr) =', params%tstep_mbal
          write(stdout,*) 'i, j =', i, j
          write(stdout,*) 'call glint_i_tstep'
       end if

       ! Do a timestep for each instance

       do i = 1,params%ninstances
          
          call glint_i_tstep(time,                    params%instances(i),       &
                             params%g_av_temp,        params%g_temp_range,       &
                             params%g_av_precip,      params%g_av_zonwind,       &
                             params%g_av_merwind,     params%g_av_humid,         &
                             params%g_av_lwdown,      params%g_av_swdown,        &
                             params%g_av_airpress,                               &
                             orog,                    orog_out_temp,             &
                             albedo_temp,             if_temp,                   &
                             vf_temp,                 sif_temp,                  &
                             svf_temp,                sd_temp,                   &
                             win_temp,                wout_temp,                 &
                             twin_temp,               twout_temp,                &
                             icevol_temp,             out_f,                     &
                             .true.,                  icets)

          if (GLC_DEBUG .and. main_task) then
             write(stdout,*) 'Finished glc_glint_ice tstep, instance =', i
             write(stdout,*) 'Upscale fields to global grid'
          end if

          ! Add this contribution to the global output
          ! Only the main task has valid values for the global output fields
          if (main_task) then
             
             if (present(orog_out)) orog_out=splice_field(orog_out,orog_out_temp, &
                  params%instances(i)%frac_cov_orog,params%cov_norm_orog)

             if (present(albedo)) albedo=splice_field(albedo,albedo_temp, &
                  params%instances(i)%frac_coverage,params%cov_normalise)

             if (present(ice_frac)) ice_frac=splice_field(ice_frac,if_temp, &
                  params%instances(i)%frac_coverage,params%cov_normalise)

             if (present(veg_frac)) veg_frac=splice_field(veg_frac,vf_temp, &
                  params%instances(i)%frac_coverage,params%cov_normalise)

             if (present(snowice_frac))snowice_frac=splice_field(snowice_frac,sif_temp, &
                  params%instances(i)%frac_coverage,params%cov_normalise)

             if (present(snowveg_frac)) snowveg_frac=splice_field(snowveg_frac, &
                  svf_temp,params%instances(i)%frac_coverage, params%cov_normalise)

             if (present(snow_depth)) snow_depth=splice_field(snow_depth, &
                  sd_temp,params%instances(i)%frac_coverage,params%cov_normalise)

             if (present(water_in)) water_in=splice_field(water_in,win_temp, &
                  params%instances(i)%frac_coverage,params%cov_normalise)

             if (present(water_out)) water_out=splice_field(water_out, &
                  wout_temp, params%instances(i)%frac_coverage,params%cov_normalise)

          end if

          ! Add total water variables to running totals
          ! WJS (1-15-13): These fields are only valid in single-task runs; multi-task
          ! runs should generate an error in glint_i_tstep if you try to compute any of
          ! these. But to be safe, we check here, too

          if (present(total_water_in))  then
             if (tasks > 1) call write_log('total_water_in is only valid when running with a single task', &
                                           GM_FATAL, __FILE__, __LINE__)

             total_water_in  = total_water_in  + twin_temp
          end if

          if (present(total_water_out)) then
             if (tasks > 1) call write_log('total_water_out is only valid when running with a single task', &
                                           GM_FATAL, __FILE__, __LINE__)

             total_water_out = total_water_out + twout_temp
          end if

          if (present(ice_volume)) then
             if (tasks > 1) call write_log('ice_volume is only valid when running with a single task', &
                                           GM_FATAL, __FILE__, __LINE__)

             ice_volume      = ice_volume      + icevol_temp
          end if

          ! Set flag
          if (present(ice_tstep)) then
             ice_tstep = (ice_tstep .or. icets)
          end if

       enddo    ! ninstances

       ! Scale output water fluxes to be in mm/s

       if (present(water_in))  water_in  = water_in/ &
                                           (params%tstep_mbal*hours2seconds)

       if (present(water_out)) water_out = water_out/ &
                                           (params%tstep_mbal*hours2seconds)

       ! ---------------------------------------------------------
       ! Reset averaging fields, flags and counters
       ! ---------------------------------------------------------

       params%g_av_temp    = 0.0
       params%g_av_precip  = 0.0
       params%g_av_zonwind = 0.0
       params%g_av_merwind = 0.0
       params%g_av_humid   = 0.0
       params%g_av_lwdown  = 0.0
       params%g_av_swdown  = 0.0
       params%g_av_airpress = 0.0
       params%g_temp_range = 0.0
       params%g_max_temp   = -1000.0
       params%g_min_temp   = 1000.0

       params%av_steps      = 0
       params%new_av        = .true.
       params%next_av_start = time+params%time_step

       deallocate(albedo_temp,if_temp,vf_temp,sif_temp,svf_temp,sd_temp,wout_temp,win_temp,orog_out_temp)

    endif    ! time - params%av_start_time + params%time_step > params%tstep_mbal

  end subroutine glint

  !===================================================================

  subroutine glint_gcm(params,         time,            &
                       qsmb,           tsfc,            &
                       topo,                            &
                       output_flag,    ice_tstep,       &
                       gfrac,          gtopo,           &
                       grofi,          grofl,           &
                       ghflx)

    ! Main Glint subroutine for GCM coupling.
    !
    ! It does all necessary spatial and temporal averaging, 
    ! and calls the dynamic ice sheet model when required. 
    !
    ! Input fields should be taken as means over the period since the last call.
    ! See the user documentation for more information.
    !
    ! Global output fields are only valid on the main task. Fields that are integrated
    ! over the whole domain are only valid in single-task runs; trying to compute these 
    ! in multi-task runs will generate a fatal error.

    use glimmer_utils
    use glint_interp
    use glint_timestep
    use glimmer_log
    use glimmer_paramets, only: scyr
    use parallel, only: main_task, tasks

    implicit none

    ! Subroutine argument declarations -------------------------------------------------------------

    type(glint_params),              intent(inout) :: params          !*FD parameters for this run
    integer,                         intent(in)    :: time            !*FD Current model time        (hours)

    real(rk),dimension(:,:,:),intent(in)    :: qsmb          ! input surface mass balance of glacier ice (kg/m^2/s)
    real(rk),dimension(:,:,:),intent(in)    :: tsfc          ! input surface ground temperature (deg C)
    real(rk),dimension(:,:,:),intent(in)    :: topo          ! input surface elevation (m)

    !TODO - Do we need both of these?
    logical,                optional,intent(out)   :: output_flag     ! Set true if outputs are set
    logical,                optional,intent(out)   :: ice_tstep       ! Set when an ice dynamic timestep has been done
                                                                      !  and new output is available

    real(rk),dimension(:,:,:),optional,intent(inout) :: gfrac         ! output ice fractional area [0,1]
    real(rk),dimension(:,:,:),optional,intent(inout) :: gtopo         ! output surface elevation (m)
    real(rk),dimension(:,:,:),optional,intent(inout) :: grofi         ! output ice runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(inout) :: grofl         ! output liquid runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(inout) :: ghflx         ! output heat flux (W/m^2, positive down)

    ! Internal variables ----------------------------------------------------------------------------

    integer :: i, j, ig, jg, n

    logical :: icets
    character(250) :: message

    real(rk), dimension(:,:,:), allocatable ::   &
       gfrac_temp    ,&! gfrac for a single instance
       gtopo_temp    ,&! gtopo for a single instance
       grofi_temp    ,&! grofi for a single instance
       grofl_temp    ,&! grofl for a single instance
       ghflx_temp      ! ghflx for a single instance

    if (GLC_DEBUG .and. main_task) then
!       write (stdout,*) 'In subroutine glint_gcm, current time (hr) =', time
!       write (stdout,*) 'av_start_time =', params%av_start_time
!       write (stdout,*) 'next_av_start =', params%next_av_start
!       write (stdout,*) 'new_av =', params%new_av
!       write (stdout,*) 'tstep_mbal =', params%tstep_mbal
    end if

    ! Check we're expecting a call now --------------------------------------------------------------

    if (params%new_av) then
       if (time == params%next_av_start) then
          params%av_start_time = time
          params%new_av = .false.
       else
          write(message,*) 'Unexpected calling of GLINT at time ', time
          call write_log(message,GM_FATAL,__FILE__,__LINE__)
       end if
    else
       if (mod (time - params%av_start_time, params%time_step) /= 0) then
          write(message,*) 'Unexpected calling of GLINT at time ', time
          call write_log(message,GM_FATAL,__FILE__,__LINE__)
       end if
    end if

    ! Check input fields are correct ----------------------------------------------------------------

    ! Reset output flag

    if (present(output_flag)) output_flag = .false.
    if (present(ice_tstep))   ice_tstep = .false.

    ! Accumulate input fields for later averaging

    call accumulate_averages_gcm(params, qsmb, tsfc, topo)

    ! Increment step counter

    params%av_steps = params%av_steps + 1

    ! ---------------------------------------------------------
    ! If this is a mass balance timestep, prepare global fields, and do a timestep
    ! for each model instance
    ! ---------------------------------------------------------

    if (time - params%av_start_time + params%time_step > params%tstep_mbal) then

       write(message,*) &
            'Incomplete forcing of GLINT mass-balance time-step detected at time ', time
       call write_log(message,GM_FATAL,__FILE__,__LINE__)

    else if (time - params%av_start_time + params%time_step == params%tstep_mbal) then

       ! Set output_flag

       ! At present, outputs are done for each mass-balance timestep, since
       ! that involved least change to the code. However, it might be good
       ! to change the output to occur with user-specified frequency.

       if (present(output_flag)) output_flag = .true.

       ! Allocate output fields
       ! Each *_temp field contains the output for one ice sheet instance.
       ! If there are multiple instances, the various *_temp fields are spliced together.

       allocate(gfrac_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
       allocate(gtopo_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
       allocate(grofi_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
       allocate(grofl_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
       allocate(ghflx_temp(params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))

       ! Zero global outputs if present

       !TODO - Change to dp
       if (present(gfrac)) gfrac = 0.0
       if (present(gtopo)) gtopo = 0.0
       if (present(grofi)) grofi = 0.0
       if (present(grofl)) grofl = 0.0
       if (present(ghflx)) ghflx = 0.0

       ! Calculate averages by dividing by number of steps elapsed
       ! since last model timestep.

       call calculate_averages_gcm(params)

       if (GLC_DEBUG .and. main_task) then
          i = params%instances(1)%model%numerics%idiag_global
          j = params%instances(1)%model%numerics%jdiag_global
          write(stdout,*) 'Take a mass balance timestep, time (hr) =', time
          write(stdout,*) 'av_steps =', real(params%av_steps,rk)
          write(stdout,*) 'tstep_mbal (hr) =', params%tstep_mbal
          write(stdout,*) 'i, j =', i, j
          if (params%gcm_smb) then
             do n = 1, params%g_grid%nec
                write (stdout,*) ' '
                write (stdout,*) 'n =', n
                write (stdout,*) 'g_av_topo (m)    =', params%g_av_topo(i,j,n)
                write (stdout,*) 'g_av_tsfc (degC) =', params%g_av_tsfc(i,j,n)
                write (stdout,*) 'g_av_qsmb (m/yr) =', params%g_av_qsmb(i,j,n)*scyr/1000.d0
!                write (stdout,*) 'g_av_qsmb (kg/m2/s) =', params%g_av_qsmb(i,j,n)
             enddo
          endif
          write(stdout,*) 'call glint_i_tstep'
       end if

       ! Calculate total surface mass balance - multiply by time since last model timestep
       ! Note on units: We want g_av_qsmb to have units of m per time step.
       ! Divide by 1000 to convert from mm to m.
       ! Multiply by hours2seconds = 3600 to convert from 1/s to 1/hr.  (tstep_mbal has units of hours)

       params%g_av_qsmb(:,:,:) = params%g_av_qsmb(:,:,:) * params%tstep_mbal * hours2seconds / 1000._rk

       ! Do a timestep for each instance

       do i = 1, params%ninstances
          
          call glint_i_tstep_gcm(time,                  &
                                 params%instances(i),   &
                                 icets,                 &
                                 params%g_av_qsmb,      &
                                 params%g_av_tsfc,      &
                                 params%g_av_topo,      &
                                 gmask  = params%g_grid%mask,    &
                                 gfrac  = gfrac_temp,            &
                                 gtopo  = gtopo_temp,            &
                                 grofi  = grofi_temp,            &
                                 grofl  = grofl_temp,            &
                                 ghflx  = ghflx_temp  )

          if (GLC_DEBUG .and. main_task) then
             write(stdout,*) 'Finished glc_glint_ice tstep, instance =', i
             write(stdout,*) 'Upscale fields to global grid'
          end if

          ! Set flag
          if (present(ice_tstep)) then
             ice_tstep = (ice_tstep .or. icets)
          end if

          ! Upscale the output to elevation classes on the global grid

          call get_i_upscaled_fields_gcm(params%instances(i), params%g_grid%nec, &
                                         params%instances(i)%lgrid%size%pt(1),   &
                                         params%instances(i)%lgrid%size%pt(2),   &
                                         params%g_grid%nx,    params%g_grid%ny,  &
                                         gfrac_temp,          gtopo_temp,        &
                                         grofi_temp,          grofl_temp,        &
                                         ghflx_temp )

          if (GLC_DEBUG .and. main_task) then
             ig = iglint_global  ! in glint_type; make sure values there are appropriate
             jg = jglint_global
             write(stdout,*) ' '
             write(stdout,*) 'After upscaling, instance', i
             do n = 1, params%g_grid%nec
                write(stdout,*) ' '
                write(stdout,*) 'n =', n
                write(stdout,*) 'gfrac_temp(n) =', gfrac_temp(ig,jg,n)
                write(stdout,*) 'gtopo_temp(n) =', gtopo_temp(ig,jg,n)
                write(stdout,*) 'grofi_temp(n) =', grofi_temp(ig,jg,n)
                write(stdout,*) 'grofl_temp(n) =', grofl_temp(ig,jg,n)
                write(stdout,*) 'ghflx_temp(n) =', ghflx_temp(ig,jg,n)
             enddo
          end if

          ! Add this contribution to the global output

          call splice_fields_gcm(gfrac_temp, gtopo_temp,    &
                                 grofi_temp, grofl_temp,    &
                                 ghflx_temp,                &
                                 gfrac,      gtopo,         &
                                 grofi,      grofl,         &
                                 ghflx,                     &
                                 params%g_grid%nec,         &
                                 params%instances(i)%frac_coverage, &
                                 params%cov_normalise)
                         
       enddo    ! ninstances

       ! ---------------------------------------------------------
       ! Reset averaging fields, flags and counters
       ! ---------------------------------------------------------

       !TODO - Change to dp
       params%g_av_qsmb    = 0.0
       params%g_av_tsfc    = 0.0
       params%g_av_topo    = 0.0

       params%av_steps      = 0
       params%new_av        = .true.
       params%next_av_start = time + params%time_step

       deallocate(gfrac_temp, gtopo_temp, grofi_temp, grofl_temp, ghflx_temp)

    endif    ! time - params%av_start_time + params%time_step > params%tstep_mbal

  end subroutine glint_gcm

  !===================================================================

  subroutine end_glint(params,close_logfile)

    !*FD tidy-up operations for Glint
    use glint_initialise
    use glimmer_log
    implicit none

    type(glint_params),intent(inout) :: params          ! parameters for this run
    logical, intent(in), optional    :: close_logfile   ! if true, then close the log file
                                                        ! (GCM may do this elsewhere)                                  
    integer :: i

    ! end individual instances

    do i = 1, params%ninstances
       call glint_i_end(params%instances(i))
    enddo

    if (present(close_logfile)) then
       if (close_logfile) call close_log
    else
       call close_log
    endif

  end subroutine end_glint

  !=====================================================

  integer function glint_coverage_map(params, coverage, cov_orog)

    !*FD Retrieve ice model fractional coverage map. 
    !*FD This function is provided so that glimmer may
    !*FD be restructured without altering the interface.
    !*FD This is currently only valid on the main task.
    !*RV Three return values are possible:
    !*RV \begin{description}
    !*RV \item[0 ---] Successful return
    !*RV \item[1 ---] Coverage map not calculated yet (fail)
    !*RV \item[2 ---] Coverage array is the wrong size (fail)
    !*RV \end{description}

    implicit none

    type(glint_params),intent(in) :: params                  !*FD ice model parameters
    real(rk),dimension(:,:),intent(out) :: coverage          !*FD array to hold coverage map
    real(rk),dimension(:,:),intent(out),optional :: cov_orog !*FD Orography coverage

    if (.not. params%coverage_calculated) then
       glint_coverage_map = 1
       return
    endif

    if ( size(coverage,1) /= params%g_grid%nx .or. &
         size(coverage,2) /= params%g_grid%ny) then
       glint_coverage_map = 2
       return
    endif

    glint_coverage_map = 0
    coverage = params%total_coverage
    if (present(cov_orog)) cov_orog = params%total_cov_orog

  end function glint_coverage_map

  !=====================================================

  !----------------------------------------------------------------------
  ! PRIVATE INTERNAL GLIMMER SUBROUTINES FOLLOW.............
  !----------------------------------------------------------------------

  subroutine glint_allocate_arrays(params)

    !*FD allocates glimmer arrays

    implicit none

    type(glint_params),intent(inout) :: params !*FD ice model parameters

    allocate(params%g_av_precip (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_temp   (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_max_temp  (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_min_temp  (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_temp_range(params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_zonwind(params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_merwind(params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_humid  (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_lwdown (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_swdown (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_airpress(params%g_grid%nx,params%g_grid%ny))

    allocate(params%total_coverage(params%g_grid%nx, params%g_grid%ny))
    allocate(params%cov_normalise (params%g_grid%nx, params%g_grid%ny))
    allocate(params%total_cov_orog(params%g_grid_orog%nx, params%g_grid_orog%ny))
    allocate(params%cov_norm_orog (params%g_grid_orog%nx, params%g_grid_orog%ny))

  end subroutine glint_allocate_arrays

  !========================================================

  subroutine splice_fields_gcm(gfrac_temp, gtopo_temp,    &
                               grofi_temp, grofl_temp,    &
                               ghflx_temp,                &
                               gfrac,      gtopo,         &
                               grofi,      grofl,         &
                               ghflx,                     &
                               nec,                       &
                               frac_coverage,             &
                               cov_normalise)

     use parallel, only: main_task

     ! Add the output for this instance to the global output

     real(rk), dimension(:,:,:), intent(in) :: gfrac_temp  ! output fields for this instance
     real(rk), dimension(:,:,:), intent(in) :: gtopo_temp  ! output fields for this instance
     real(rk), dimension(:,:,:), intent(in) :: grofi_temp  ! output fields for this instance
     real(rk), dimension(:,:,:), intent(in) :: grofl_temp  ! output fields for this instance
     real(rk), dimension(:,:,:), intent(in) :: ghflx_temp  ! output fields for this instance

     real(rk), dimension(:,:,:), intent(inout) :: gfrac    ! spliced global output field
     real(rk), dimension(:,:,:), intent(inout) :: gtopo    ! spliced global output field
     real(rk), dimension(:,:,:), intent(inout) :: grofi    ! spliced global output field
     real(rk), dimension(:,:,:), intent(inout) :: grofl    ! spliced global output field
     real(rk), dimension(:,:,:), intent(inout) :: ghflx    ! spliced global output field

     integer, intent(in) :: nec   ! number of elevation classes

     real(rk), dimension(:,:), intent(in) :: frac_coverage  ! map of fractional coverage of global gridcells
                                                ! by local gridcells
     real(rk), dimension(:,:), intent(in) :: cov_normalise  ! normalisation values
     
     integer :: n

     ! Only the main task has valid values for the global output fields

     if (main_task) then

        do n = 1, nec

           gfrac(:,:,n) = splice_field(gfrac(:,:,n),            &
                                       gfrac_temp(:,:,n),       &
                                       frac_coverage,           &
                                       cov_normalise)

           gtopo(:,:,n) = splice_field(gtopo(:,:,n),            &
                                       gtopo_temp(:,:,n),       &
                                       frac_coverage,           &
                                       cov_normalise)

           grofi(:,:,n) = splice_field(grofi(:,:,n),            &
                                       grofi_temp(:,:,n),       &
                                       frac_coverage,           &
                                       cov_normalise)

           grofl(:,:,n) = splice_field(grofl(:,:,n),            &
                                       grofl_temp(:,:,n),       &
                                       frac_coverage,           &
                                       cov_normalise)

           ghflx(:,:,n) = splice_field(ghflx(:,:,n),            &
                                       ghflx_temp(:,:,n),       &
                                       frac_coverage,           &
                                       cov_normalise)

        enddo   ! nec

     endif  ! main_task

  end subroutine splice_fields_gcm

  !========================================================

  function splice_field(global, local, coverage, normalise)

    !*FD Splices an upscaled field into a global field

    real(rk),dimension(:,:),intent(in) :: global    !*FD Field to receive the splice
    real(rk),dimension(:,:),intent(in) :: local     !*FD The field to be spliced in
    real(rk),dimension(:,:),intent(in) :: coverage  !*FD The coverage fraction
    real(rk),dimension(:,:),intent(in) :: normalise !*FD The normalisation field

    real(rk),dimension(size(global,1),size(global,2)) :: splice_field

    where (coverage == 0.0)
       splice_field = global
    elsewhere
       splice_field = (global*(1.0-coverage/normalise)) + (local*coverage/normalise)
    end where

  end function splice_field

  !========================================================

  !TODO - Move this subroutine to a glint_setup module?

  subroutine glint_readconfig(config, ninstances, fnames, infnames)

    !*FD Determine whether a given config file is a
    !*FD top-level glint config file, and return parameters
    !*FD accordingly.

    use glimmer_config
    use glimmer_log
    implicit none

    ! Arguments -------------------------------------------

    type(ConfigSection),      pointer :: config !*FD structure holding sections of configuration file
    integer,              intent(out) :: ninstances !*FD Number of instances to create
    character(fname_length),dimension(:),pointer :: fnames !*FD list of filenames (output)
    character(fname_length),dimension(:) :: infnames !*FD list of filenames (input)

    ! Internal variables ----------------------------------

    type(ConfigSection), pointer :: section
    character(len=100) :: message
    integer :: i

    if (associated(fnames)) nullify(fnames)

    call GetSection(config,section,'GLINT')
    if (associated(section)) then
       call GetValue(section,'n_instance',ninstances)
       allocate(fnames(ninstances))
       do i=1,ninstances
          call GetSection(section%next,section,'GLINT instance')
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

!!$    call write_log('GLINT global')
!!$    call write_log('------------')
!!$    write(message,*) 'number of instances :',params%ninstances
!!$    call write_log(message)
!!$    call write_log('')

  end subroutine glint_readconfig

  !========================================================

  subroutine calc_bounds(lon, lat, lonb, latb)

    !*FD Calculates the boundaries between
    !*FD global grid-boxes. Note that we assume that the boundaries lie 
    !*FD half-way between the 
    !*FD points, both latitudinally and longitudinally, although 
    !*FD this isn't strictly true for a Gaussian grid.

    use glimmer_map_trans, only: loncorrect

    implicit none

    real(rk),dimension(:),intent(in) :: lon,lat    !*FD locations of global grid-points (degrees)
    real(rk),dimension(:),intent(out) :: lonb,latb !*FD boundaries of grid-boxes (degrees)

    real(rk) :: dlon

    integer :: nxg,nyg,i,j

    nxg=size(lon) ; nyg=size(lat)

    ! Latitudes first - we assume the boundaries of the first and 
    ! last boxes coincide with the poles. Not sure how to
    ! handle it if they don't...

    latb(1)=90.0
    latb(nyg+1)=-90.0

    do j=2,nyg
       latb(j)=lat(j-1)-(lat(j-1)-lat(j))/2.0
    enddo

    ! Longitudes

    if (lon(1)<lon(nxg)) then
       dlon=lon(1)-lon(nxg)+360.0
    else
       dlon=lon(1)-lon(nxg)
    endif
    lonb(1)=lon(nxg)+dlon/2
    lonb(1)=loncorrect(lonb(1),0.0_rk)      

    lonb(nxg+1)=lonb(1)

    do i=2,nxg
       if (lon(i)<lon(i-1)) then
          dlon=lon(i)-lon(i-1)+360.0
       else
          dlon=lon(i)-lon(i-1)
       endif
       lonb(i)=lon(i-1)+dlon/2
       lonb(i)=loncorrect(lonb(i),0.0_rk)      
    enddo

  end subroutine calc_bounds

  !========================================================

  integer function check_mbts(timesteps)

    !*FD Checks to see that all mass-balance time-steps are
    !*FD the same. Flags a fatal error if not, else assigns that
    !*FD value to the output

    use glimmer_log

    implicit none

    integer,dimension(:) :: timesteps !*FD Array of mass-balance timsteps

    integer :: n,i

    n=size(timesteps)
    if (n==0) then
       check_mbts=0
       return
    endif

    check_mbts=timesteps(1)

    do i=2,n
       if (timesteps(i)/=check_mbts) then
          call write_log('All instances must have the same mass-balance and ice timesteps', &
               GM_FATAL,__FILE__,__LINE__)
       endif
    enddo

  end function check_mbts

  !========================================================

  subroutine check_init_args(orog_lats, orog_longs, orog_latb, orog_lonb)

    !*FD Checks which combination arguments have been supplied to
    !*FD define the global grid, and rejects unsuitable combinations

    use glimmer_log

    real(rk),dimension(:),optional,intent(in) :: orog_lats 
    real(rk),dimension(:),optional,intent(in) :: orog_longs 
    real(rk),dimension(:),optional,intent(in) :: orog_latb
    real(rk),dimension(:),optional,intent(in) :: orog_lonb 

    integer :: args
    integer,dimension(5) :: allowed=(/0,3,7,11,15/)

    args=0

    if (present(orog_lats))  args=args+1
    if (present(orog_longs)) args=args+2
    if (present(orog_latb))  args=args+4
    if (present(orog_lonb))  args=args+8

    if (.not.any(args==allowed)) then
       call write_log('Unexpected combination of arguments to initialise_glint', &
            GM_FATAL,__FILE__,__LINE__)
    end if

  end subroutine check_init_args

  !========================================================

  subroutine check_input_fields(params, humid, lwdown, swdown, airpress, zonwind, merwind)

    use glimmer_log

    type(glint_params),              intent(inout) :: params   !*FD parameters for this run
    real(rk),dimension(:,:),optional,intent(in)    :: humid    !*FD Surface humidity (%)
    real(rk),dimension(:,:),optional,intent(in)    :: lwdown   !*FD Downwelling longwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: swdown   !*FD Downwelling shortwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: airpress !*FD surface air pressure (Pa)
    real(rk),dimension(:,:),optional,intent(in)    :: zonwind  !*FD Zonal component of the wind field (m/s)
    real(rk),dimension(:,:),optional,intent(in)    :: merwind  !*FD Meridional component of the wind field (m/s)

    if (params%enmabal) then
       if (.not.(present(humid).and.present(lwdown).and. &
            present(swdown).and.present(airpress).and. &
            present(zonwind).and.present(merwind))) &
            call write_log('Necessary fields not supplied for Energy Balance Mass Balance model',GM_FATAL, &
            __FILE__,__LINE__)
    end if

    if (params%need_winds) then
       if (.not.(present(zonwind).and.present(merwind))) &
          call write_log('Need to supply zonal and meridional wind fields to GLINT',GM_FATAL, &
            __FILE__,__LINE__)
    end if

  end subroutine check_input_fields

  !========================================================

  subroutine accumulate_averages(params, temp, precip, zonwind, merwind, humid, lwdown, swdown, airpress)

    type(glint_params),              intent(inout) :: params   !*FD parameters for this run
    real(rk),dimension(:,:),         intent(in)    :: temp     !*FD Surface temperature field (celsius)
    real(rk),dimension(:,:),         intent(in)    :: precip   !*FD Precipitation rate        (mm/s)
    real(rk),dimension(:,:),optional,intent(in)    :: zonwind  !*FD Zonal component of the wind field (m/s)
    real(rk),dimension(:,:),optional,intent(in)    :: merwind  !*FD Meridional component of the wind field (m/s)
    real(rk),dimension(:,:),optional,intent(in)    :: humid    !*FD Surface humidity (%)
    real(rk),dimension(:,:),optional,intent(in)    :: lwdown   !*FD Downwelling longwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: swdown   !*FD Downwelling shortwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: airpress !*FD surface air pressure (Pa)

    params%g_av_temp    = params%g_av_temp    + temp
    params%g_av_precip  = params%g_av_precip  + precip

    if (params%need_winds) params%g_av_zonwind = params%g_av_zonwind + zonwind
    if (params%need_winds) params%g_av_merwind = params%g_av_merwind + merwind

    if (params%enmabal) then
       params%g_av_humid    = params%g_av_humid    + humid
       params%g_av_lwdown   = params%g_av_lwdown   + lwdown
       params%g_av_swdown   = params%g_av_swdown   + swdown
       params%g_av_airpress = params%g_av_airpress + airpress
    endif

    ! Ranges of temperature

    where (temp > params%g_max_temp) params%g_max_temp = temp
    where (temp < params%g_min_temp) params%g_min_temp = temp

  end subroutine accumulate_averages

  !========================================================

  subroutine accumulate_averages_gcm(params, qsmb, tsfc, topo)

    type(glint_params), intent(inout)   :: params     ! model parameters
    real(rk),dimension(:,:,:),intent(in)  :: qsmb     ! flux of glacier ice (kg/m^2/s)
    real(rk),dimension(:,:,:),intent(in)  :: tsfc     ! surface ground temperature (C)
    real(rk),dimension(:,:,:),intent(in)  :: topo     ! surface elevation (m)

    params%g_av_qsmb(:,:,:) = params%g_av_qsmb(:,:,:) + qsmb(:,:,:)
    params%g_av_tsfc(:,:,:) = params%g_av_tsfc(:,:,:) + tsfc(:,:,:)
    params%g_av_topo(:,:,:) = params%g_av_topo(:,:,:) + topo(:,:,:)

    !TODO - No need to accumulate and average topo?
    !       topo should not change during a mass balance timestep.

  end subroutine accumulate_averages_gcm

  !========================================================

  subroutine calculate_averages(params)

    type(glint_params),              intent(inout) :: params   !*FD parameters for this run

    params%g_av_temp    = params%g_av_temp   /real(params%av_steps)
    params%g_av_precip  = params%g_av_precip /real(params%av_steps)
    if (params%need_winds) params%g_av_zonwind = params%g_av_zonwind/real(params%av_steps)
    if (params%need_winds) params%g_av_merwind = params%g_av_merwind/real(params%av_steps)
    if (params%enmabal) then
       params%g_av_humid    = params%g_av_humid   /real(params%av_steps)
       params%g_av_lwdown   = params%g_av_lwdown  /real(params%av_steps)
       params%g_av_swdown   = params%g_av_swdown  /real(params%av_steps)
       params%g_av_airpress = params%g_av_airpress/real(params%av_steps)
    endif

  end subroutine calculate_averages

  !========================================================

  subroutine calculate_averages_gcm(params)

    type(glint_params),              intent(inout) :: params   !*FD parameters for this run

    !TODO - Do not average topo?
    params%g_av_qsmb(:,:,:) = params%g_av_qsmb(:,:,:) / real(params%av_steps,rk)
    params%g_av_tsfc(:,:,:) = params%g_av_tsfc(:,:,:) / real(params%av_steps,rk)
    params%g_av_topo(:,:,:) = params%g_av_topo(:,:,:) / real(params%av_steps,rk)

  end subroutine calculate_averages_gcm

  !========================================================

  subroutine populate_output_flags(out_f,                     &
                                   orog_out,       albedo,    &
                                   ice_frac,       veg_frac,  &
                                   snowice_frac,   snowveg_frac,  &
                                   snow_depth,                &
                                   water_in,       water_out, &
                                   total_water_in, total_water_out, &
                                   ice_volume)

    type(output_flags),intent(inout) :: out_f

    real(rk),dimension(:,:),optional,intent(inout) :: orog_out        !*FD The fed-back, output orography (m)
    real(rk),dimension(:,:),optional,intent(inout) :: albedo          !*FD surface albedo
    real(rk),dimension(:,:),optional,intent(inout) :: ice_frac        !*FD grid-box ice-fraction
    real(rk),dimension(:,:),optional,intent(inout) :: veg_frac        !*FD grid-box veg-fraction
    real(rk),dimension(:,:),optional,intent(inout) :: snowice_frac    !*FD grid-box snow-covered ice fraction
    real(rk),dimension(:,:),optional,intent(inout) :: snowveg_frac    !*FD grid-box snow-covered veg fraction
    real(rk),dimension(:,:),optional,intent(inout) :: snow_depth      !*FD grid-box mean snow depth (m water equivalent)
    real(rk),dimension(:,:),optional,intent(inout) :: water_in        !*FD Input water flux          (mm)
    real(rk),dimension(:,:),optional,intent(inout) :: water_out       !*FD Output water flux         (mm)
    real(rk),               optional,intent(inout) :: total_water_in  !*FD Area-integrated water flux in (kg)
    real(rk),               optional,intent(inout) :: total_water_out !*FD Area-integrated water flux out (kg)
    real(rk),               optional,intent(inout) :: ice_volume      !*FD Total ice volume (m$^3$)


    out_f%orog         = present(orog_out)
    out_f%albedo       = present(albedo)
    out_f%ice_frac     = present(ice_frac)
    out_f%veg_frac     = present(veg_frac)
    out_f%snowice_frac = present(snowice_frac)
    out_f%snowveg_frac = present(snowveg_frac)
    out_f%snow_depth   = present(snow_depth)
    out_f%water_out    = present(water_out)
    out_f%water_in     = present(water_in)
    out_f%total_win    = present(total_water_in)
    out_f%total_wout   = present(total_water_out)
    out_f%ice_vol      = present(ice_volume)

  end subroutine populate_output_flags

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glint_main

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
