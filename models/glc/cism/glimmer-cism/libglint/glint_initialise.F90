! WJS (1-30-12): The following (turning optimization off) is needed as a workaround for an
! xlf compiler bug, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

#ifdef CPRIBM
@PROCESS ALIAS_SIZE(107374182)
#endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_initialise.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module glint_initialise

  !*FD Initialise GLINT model instance

  use glint_type

  implicit none

  private
  public glint_i_initialise, glint_i_initialise_gcm, glint_i_end, calc_coverage

contains

  subroutine glint_i_initialise(config,           instance,         &
                                grid,             grid_orog,        &
                                mbts,             idts,             &
                                need_winds,       enmabal,          &
                                force_start,      force_dt,         &
                                gcm_restart,      gcm_restart_file, &
                                gcm_config_unit)

    !*FD Initialise a GLINT ice model instance

    use glimmer_paramets, only: GLC_DEBUG
    use glimmer_log
    use glimmer_config
    use glint_global_grid
    use glint_mbal_coupling
    use glint_io          , only: glint_io_createall     , glint_io_writeall
    use glint_mbal_io     , only: glint_mbal_io_createall, glint_mbal_io_writeall
    use glimmer_ncio
    use glide_nc_custom   , only: glide_nc_fillall
    use glide
    use glissade
    use glint_constants
    use glimmer_restart_gcm
    use glide_diagnostics
    use parallel, only: main_task

    implicit none

    ! Arguments
    type(ConfigSection), pointer         :: config           !*FD structure holding sections of configuration file   
    type(glint_instance),  intent(inout) :: instance         !*FD The instance being initialised.
    type(global_grid),     intent(in)    :: grid             !*FD Global grid to use
    type(global_grid),     intent(in)    :: grid_orog        !*FD Global grid to use for orography
    integer,               intent(out)   :: mbts             !*FD mass-balance time-step (hours)
    integer,               intent(out)   :: idts             !*FD ice dynamics time-step (hours)
    logical,               intent(inout) :: need_winds       !*FD Set if this instance needs wind input
    logical,               intent(inout) :: enmabal          !*FD Set if this instance uses the energy balance
                                                             !    mass-bal model
    integer,               intent(in)    :: force_start      !*FD glint forcing start time (hours)
    integer,               intent(in)    :: force_dt         !*FD glint forcing time step (hours)
    logical,     optional, intent(in)    :: gcm_restart      !*FD logical flag to read from a restart file
    character(*),optional, intent(in)    :: gcm_restart_file !*FD restart filename for restart
    integer,     optional, intent(in)    :: gcm_config_unit  !*FD fileunit for reading config files

    ! Internal
    real(sp),dimension(:,:),allocatable :: thk
    integer :: config_fileunit, restart_fileunit
    real(dp) :: timeyr       ! time in years

    config_fileunit = 99
    if (present(gcm_config_unit)) then
       config_fileunit = gcm_config_unit
    endif

    ! initialise model

    call glide_config(instance%model, config, config_fileunit)

    ! if this is a continuation run, then set up to read restart
    ! (currently assumed to be a CESM restart file)

    if (present(gcm_restart)) then

      if (gcm_restart) then

         if (present(gcm_restart_file)) then

            ! read the restart file
            call glimmer_read_restart_gcm(instance%model, gcm_restart_file)
            instance%model%options%is_restart = 1
 
         else

            call write_log('Missing gcm_restart_file when gcm_restart is true',&
                           GM_FATAL,__FILE__,__LINE__)

         endif

      endif
    endif

    if (instance%model%options%whichdycore == DYCORE_GLIDE) then  ! SIA dycore

       ! initialise the model
       call glide_initialise(instance%model)

       !TODO - Remove the oldglide option?

       ! compute the initial diagnostic state
       !WHL - Do not call this if comparing to oldglide (cism1) results
      if (.not. oldglide) then
       call glide_init_state_diagnostic(instance%model)
      endif

    else       ! glam/glissade HO dycore     

       ! initialise the model
       call glissade_initialise(instance%model)

       ! compute the initial diagnostic state
       call glissade_diagnostic_variable_solve(instance%model)

    endif

    instance%ice_tstep = get_tinc(instance%model)*years2hours
    instance%glide_time = instance%model%numerics%tstart
    idts = instance%ice_tstep

    ! read glint configuration

    call glint_i_readconfig(instance, config)    
    call glint_i_printconfig(instance)    

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file (specific to glint - other configs,
    ! e.g. glide, isos, are handled separately by their setup routines).
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/17/13
    ! Note: the corresponding call for glide is placed within *_readconfig, which is probably more appropriate,
    ! but putting this call into glint_i_readconfig creates a circular dependency.  

    call define_glint_restart_variables(instance)
 
    ! create glint variables for the glide output files
    call glint_io_createall(instance%model, data=instance)

    ! create instantaneous glint variables

    call openall_out(instance%model, outfiles=instance%out_first)
    call glint_mbal_io_createall(instance%model, data=instance, outfiles=instance%out_first)

    ! fill dimension variables

    call glide_nc_fillall(instance%model)
    call glide_nc_fillall(instance%model, outfiles=instance%out_first)

    ! Check we've used all the config sections

    call CheckSections(config)

    ! New grid (grid on this task)

    ! WJS (1-11-13): I'm not sure if it's correct to set the origin to (0,0) when running
    ! on multiple tasks, with a decomposed grid. However, as far as I can tell, the
    ! origin of this variable isn't important, so I'm not trying to fix it right now.

    instance%lgrid = coordsystem_new(0.d0, 0.d0, &
                                     get_dew(instance%model), &
                                     get_dns(instance%model), &
                                     get_ewn(instance%model), &
                                     get_nsn(instance%model))

    ! Allocate arrays appropriately

    call glint_i_allocate(instance, grid%nx, grid%ny, grid_orog%nx, grid_orog%ny)

    ! Read data and initialise climate

    call glint_i_readdata(instance)

    ! Create grid spanning full domain and other information needed for downscaling &
    ! upscaling. Note that, currently, these variables only have valid data on the main
    ! task, since all downscaling & upscaling is done there
    
    call setup_lgrid_fulldomain(instance, grid, grid_orog)

    ! initialise the mass-balance accumulation

    call glint_mbc_init(instance%mbal_accum, &
                        instance%lgrid, &
                        config,         &
                        instance%whichacab, &
                        instance%snowd, &
                        instance%siced, &
                        instance%lgrid%size%pt(1), &
                        instance%lgrid%size%pt(2), &
                        real(instance%lgrid%delta%pt(1),rk))

    instance%mbal_tstep = instance%mbal_accum%mbal%tstep
    mbts = instance%mbal_tstep

    instance%next_time = force_start-force_dt+instance%mbal_tstep

    if (GLC_DEBUG .and. main_task) then
       write (6,*) 'Called glint_mbc_init'
       write (6,*) 'mbal tstep =', mbts
       write (6,*) 'next_time =', instance%next_time
       write (6,*) 'start_time =', instance%mbal_accum%start_time
    end if

    ! Mass-balance accumulation length

    if (instance%mbal_accum_time == -1) then
       instance%mbal_accum_time = max(instance%ice_tstep,instance%mbal_tstep)
    end if

    if (instance%mbal_accum_time < instance%mbal_tstep) then
       call write_log('Mass-balance accumulation timescale must be as '//&
            'long as mass-balance time-step',GM_FATAL,__FILE__,__LINE__)
    end if

    if (mod(instance%mbal_accum_time,instance%mbal_tstep) /= 0) then
       call write_log('Mass-balance accumulation timescale must be an '// &
            'integer multiple of the mass-balance time-step',GM_FATAL,__FILE__,__LINE__)
    end if

    if (.not.(mod(instance%mbal_accum_time,instance%ice_tstep)==0 .or. &
         mod(instance%ice_tstep,instance%mbal_accum_time)==0)) then
       call write_log('Mass-balance accumulation timescale and ice dynamics '//&
            'timestep must divide into one another',GM_FATAL,__FILE__,__LINE__)
    end if

    if (instance%ice_tstep_multiply/=1 .and. mod(instance%mbal_accum_time,int(years2hours))/=0.0) then
       call write_log('For ice time-step multiplication, mass-balance accumulation timescale '//&
            'must be an integer number of years',GM_FATAL,__FILE__,__LINE__)
    end if

    ! Initialise some other stuff

    if (instance%mbal_accum_time>instance%ice_tstep) then
       instance%n_icetstep = instance%ice_tstep_multiply*instance%mbal_accum_time/instance%ice_tstep
    else
       instance%n_icetstep = instance%ice_tstep_multiply
    end if

    !This was commented out because it destroys exact restart
    !TODO - Find another way to set thk to snowd?
    ! Copy snow-depth to thickness if no thickness is present

!!    allocate(thk(get_ewn(instance%model),get_nsn(instance%model)))
!!    call glide_get_thk(instance%model,thk)
!!    where (instance%snowd>0.0 .and. thk==0.0)
!!       thk=instance%snowd
!!    elsewhere
!!       thk=thk
!!    endwhere
!!    call glide_set_thk(instance%model,thk)
!!    deallocate(thk)

   ! Write initial ice sheet diagnostics for this instance

    call glide_write_diagnostics(instance%model,                  &
                                 instance%model%numerics%time,    &
                                 tstep_count = instance%model%numerics%timecounter)

    ! Write netCDF output for this instance

    call glide_io_writeall(instance%model, instance%model)
    call glint_io_writeall(instance, instance%model)
    call glint_mbal_io_writeall(instance, instance%model, outfiles=instance%out_first)

    if (instance%whichprecip == 2) need_winds=.true.
    if (instance%whichacab == 3) then
       need_winds = .true.
       enmabal = .true.
    end if

  end subroutine glint_i_initialise

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_initialise_gcm(config,           instance,         &
                                    grid,             &
                                    mbts,             idts,             &
                                    force_start,      force_dt,         &
                                    gcm_restart,      gcm_restart_file, &
                                    gcm_config_unit)

    ! Initialise a GLINT ice model instance for GCM coupling

    use glimmer_paramets, only: GLC_DEBUG
    use glimmer_log
    use glimmer_config
    use glint_global_grid
    use glint_mbal_coupling
    use glint_io          , only: glint_io_createall     , glint_io_writeall
    use glint_mbal_io     , only: glint_mbal_io_createall, glint_mbal_io_writeall
    use glimmer_ncio
    use glide_nc_custom   , only: glide_nc_fillall
    use glide
    use glissade
    use glint_constants
    use glimmer_restart_gcm
    use glide_diagnostics
    use parallel, only: main_task

    implicit none

    ! Arguments
    type(ConfigSection), pointer         :: config           ! structure holding sections of configuration file   
    type(glint_instance),  intent(inout) :: instance         ! The instance being initialised.
    type(global_grid),     intent(in)    :: grid             ! Global grid to use
    integer,               intent(out)   :: mbts             ! mass-balance time-step (hours)
    integer,               intent(out)   :: idts             ! ice dynamics time-step (hours)

    !TODO - Change names of force_start and force_dt?
    integer,               intent(in)    :: force_start      ! glint forcing start time (hours)
    integer,               intent(in)    :: force_dt         ! glint forcing time step (hours)

    logical,     optional, intent(in)    :: gcm_restart      ! logical flag to read from a restart file
    character(*),optional, intent(in)    :: gcm_restart_file ! restart filename for restart
    integer,     optional, intent(in)    :: gcm_config_unit  ! fileunit for reading config files

    ! Internal

    integer :: config_fileunit, restart_fileunit

    config_fileunit = 99
    if (present(gcm_config_unit)) then
       config_fileunit = gcm_config_unit
    endif

    ! initialise model

    call glide_config(instance%model, config, config_fileunit)

    ! if this is a continuation run, then set up to read restart
    ! (currently assumed to be a CESM restart file)

    if (present(gcm_restart)) then

      if (gcm_restart) then

         if (present(gcm_restart_file)) then

            ! read the restart file
            call glimmer_read_restart_gcm(instance%model, gcm_restart_file)
            instance%model%options%is_restart = 1
 
         else

            call write_log('Missing gcm_restart_file when gcm_restart is true',&
                           GM_FATAL,__FILE__,__LINE__)

         endif

      endif
    endif

    if (instance%model%options%whichdycore == DYCORE_GLIDE) then  ! SIA dycore

       ! initialise the model
       call glide_initialise(instance%model)

       !TODO - Remove the oldglide option?

       ! compute the initial diagnostic state
!WHL - Do not call this if comparing to oldglide (cism1) results
      if (.not. oldglide) then
       call glide_init_state_diagnostic(instance%model)
      endif

    else       ! glam/glissade HO dycore     

       ! initialise the model
       call glissade_initialise(instance%model)

       ! compute the initial diagnostic state
       call glissade_diagnostic_variable_solve(instance%model)

    endif

    instance%ice_tstep = get_tinc(instance%model)*years2hours
    idts = instance%ice_tstep

    instance%glide_time = instance%model%numerics%tstart

    ! read glint configuration

    call glint_i_readconfig(instance, config)    
    call glint_i_printconfig(instance)    

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file (specific to glint - other configs,
    ! e.g. glide, isos, are handled separately by their setup routines).
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/17/13
    ! Note: the corresponding call for glide is placed within *_readconfig, which is probably more appropriate,
    ! but putting this call into glint_i_readconfig creates a circular dependency.  

    call define_glint_restart_variables(instance)
 
    ! create glint variables for the glide output files
    call glint_io_createall(instance%model, data=instance)

    ! create instantaneous glint variables
    call openall_out(instance%model, outfiles=instance%out_first)
    call glint_mbal_io_createall(instance%model, data=instance, outfiles=instance%out_first)

    ! fill dimension variables
    call glide_nc_fillall(instance%model)
    call glide_nc_fillall(instance%model, outfiles=instance%out_first)

    ! Check we've used all the config sections

    call CheckSections(config)

    ! New grid (grid on this task)

    ! WJS (1-11-13): I'm not sure if it's correct to set the origin to (0,0) when running
    ! on multiple tasks, with a decomposed grid. However, as far as I can tell, the
    ! origin of this variable isn't important, so I'm not trying to fix it right now.

    instance%lgrid = coordsystem_new(0.d0, 0.d0, &
                                     get_dew(instance%model), &
                                     get_dns(instance%model), &
                                     get_ewn(instance%model), &
                                     get_nsn(instance%model))

    ! Allocate arrays appropriately

    call glint_i_allocate_gcm(instance, grid%nx, grid%ny)

    ! Read data and initialise climate

    call glint_i_readdata(instance)

    ! Create grid spanning full domain and other information needed for downscaling &
    ! upscaling. Note that, currently, these variables only have valid data on the main
    ! task, since all downscaling & upscaling is done there
    
    call setup_lgrid_fulldomain(instance, grid)

    ! initialise the mass-balance accumulation

    call glint_mbc_init_gcm(instance%mbal_accum, &
                            instance%lgrid,      &
                            instance%whichacab)

    !TODO - Do we need two copies of this tstep variable?
    instance%mbal_tstep = instance%mbal_accum%mbal%tstep

    mbts = instance%mbal_tstep

    instance%next_time = force_start - force_dt + instance%mbal_tstep

    if (GLC_DEBUG .and. main_task) then
       write (6,*) 'Called glint_mbc_init'
       write (6,*) 'mbal tstep =', mbts
       write (6,*) 'next_time =', instance%next_time
       write (6,*) 'start_time =', instance%mbal_accum%start_time
    end if

    ! Mass-balance accumulation length

    if (instance%mbal_accum_time == -1) then
       instance%mbal_accum_time = max(instance%ice_tstep,instance%mbal_tstep)
    end if

    if (instance%mbal_accum_time < instance%mbal_tstep) then
       call write_log('Mass-balance accumulation timescale must be as '//&
                      'long as mass-balance time-step',GM_FATAL,__FILE__,__LINE__)
    end if

    if (mod(instance%mbal_accum_time,instance%mbal_tstep) /= 0) then
       call write_log('Mass-balance accumulation timescale must be an '// &
                      'integer multiple of the mass-balance time-step',GM_FATAL,__FILE__,__LINE__)
    end if

    if (.not. (mod(instance%mbal_accum_time, instance%ice_tstep)==0 .or.   &
               mod(instance%ice_tstep, instance%mbal_accum_time)==0)) then
       call write_log('Mass-balance accumulation timescale and ice dynamics '//&
                      'timestep must divide into one another',GM_FATAL,__FILE__,__LINE__)
    end if

    if (instance%ice_tstep_multiply/=1 .and. mod(instance%mbal_accum_time,int(years2hours)) /= 0.0) then
       call write_log('For ice time-step multiplication, mass-balance accumulation timescale '//&
                      'must be an integer number of years',GM_FATAL,__FILE__,__LINE__)
    end if

    ! Initialise some other stuff

    if (instance%mbal_accum_time>instance%ice_tstep) then
       instance%n_icetstep = instance%ice_tstep_multiply*instance%mbal_accum_time/instance%ice_tstep
    else
       instance%n_icetstep = instance%ice_tstep_multiply
    end if

   ! Write initial ice sheet diagnostics for this instance

    call glide_write_diagnostics(instance%model,                  &
                                 instance%model%numerics%time,    &
                                 tstep_count = instance%model%numerics%timecounter)

    ! Write netCDF output for this instance

    call glide_io_writeall(instance%model, instance%model)
    call glint_io_writeall(instance, instance%model)
    call glint_mbal_io_writeall(instance, instance%model, outfiles=instance%out_first)

  end subroutine glint_i_initialise_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_end(instance)

    !*FD Tidy up 

    use glide
    use glimmer_ncio
    implicit none
    type(glint_instance),  intent(inout) :: instance    !*FD The instance being initialised.

    call glide_finalise(instance%model)
    call closeall_out(instance%model,outfiles=instance%out_first)
    instance%out_first => null()

  end subroutine glint_i_end

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_readdata(instance)
    !*FD read data from netCDF file and initialise climate

    use glint_io
    use glide_thck, only: glide_calclsrf
    implicit none

    type(glint_instance),intent(inout)   :: instance    !*FD Instance whose elements are to be allocated.

    ! read data
    call glint_io_readall(instance,instance%model)

    call glide_calclsrf(instance%model%geometry%thck,instance%model%geometry%topg, &
         instance%model%climate%eus,instance%model%geometry%lsrf)
    instance%model%geometry%usrf = instance%model%geometry%thck + instance%model%geometry%lsrf

  end subroutine glint_i_readdata

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine setup_lgrid_fulldomain(instance, grid, grid_orog)

    !*FD Set up the local (icesheet) grid spanning the full domain (i.e., across all tasks).
    !*FD This also sets up auxiliary variables that depend on this full domain lgrid,
    !*FD such as the downscaling and upscaling derived types.
    !*FD This routine is required because we currently do downscaling and upscaling just
    !*FD on the main task, with the appropriate gathers / scatters.
    !*FD Thus, this creates a lgrid spanning the full domain on the main task;
    !*FD other tasks are left with an uninitialized lgrid_fulldomain.
    !*FD Tasks other than the main task also have uninitialized ups, downs and frac_coverage
    !*FD (along with the similar *_orog variables).
    
    use glint_type         , only : glint_instance
    use glint_global_grid  , only : global_grid
    use parallel           , only : main_task, global_ewn, global_nsn, distributed_gather_var
    use glimmer_coordinates, only : coordsystem_new
    use glide_types        , only : get_dew, get_dns

    implicit none

    ! Arguments

    type(glint_instance), intent(inout) :: instance
    type(global_grid)   , intent(in)    :: grid
    type(global_grid)   , intent(in), optional :: grid_orog

    ! Internal variables

    integer, dimension(:,:), allocatable :: out_mask_fulldomain

    ! Beginning of code

    call distributed_gather_var(instance%out_mask, out_mask_fulldomain)

    if (main_task) then

       instance%lgrid_fulldomain = coordsystem_new(0.d0, 0.d0, &
                                                   get_dew(instance%model), &
                                                   get_dns(instance%model), &
                                                   global_ewn, &
                                                   global_nsn)
    
       call new_downscale(instance%downs, instance%model%projection, grid, &
                          instance%lgrid_fulldomain, mpint=(instance%use_mpint==1))

       call new_upscale(instance%ups, grid,  instance%model%projection, &
                        out_mask_fulldomain, instance%lgrid_fulldomain) ! Initialise upscaling parameters

       if (present(grid_orog)) then
          call new_upscale(instance%ups_orog, grid_orog, instance%model%projection, &
                           out_mask_fulldomain, instance%lgrid_fulldomain) ! Initialise upscaling parameters
       endif

       call calc_coverage(instance%lgrid_fulldomain, &
                          instance%ups,   &             
                          grid,           &
                          out_mask_fulldomain, &
                          instance%frac_coverage)

       if (present(grid_orog)) then
          call calc_coverage(instance%lgrid_fulldomain, &               
                             instance%ups_orog,  &             
                             grid_orog,     &
                             out_mask_fulldomain, &
                             instance%frac_cov_orog)
       endif

    end if

  end subroutine setup_lgrid_fulldomain

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_coverage(lgrid_fulldomain, ups,     grid,   &
                           mask_fulldomain,  frac_coverage)

    ! Calculates the fractional coverage of the global grid-boxes by the ice model domain

    use glimmer_map_types
    use glimmer_coordinates
    use glint_global_grid

    ! Arguments

    type(coordsystem_type), intent(in)  :: lgrid_fulldomain  !*FD Local grid, spanning full domain (all tasks)
    type(upscale),          intent(in)  :: ups               !*FD Upscaling used
    type(global_grid),      intent(in)  :: grid              !*FD Global grid used
    integer, dimension(:,:),intent(in)  :: mask_fulldomain   !*FD Mask of points for upscaling, spanning full domain (all tasks)
    real(rk),dimension(:,:),intent(out) :: frac_coverage     !*FD Map of fractional 
                                                             !*FD coverage of global by local grid-boxes.
    ! Internal variables

    integer,dimension(grid%nx,grid%ny) :: tempcount
    integer :: i,j

    ! Beginning of code

    tempcount=0

    do i=1,lgrid_fulldomain%size%pt(1)
       do j=1,lgrid_fulldomain%size%pt(2)
          tempcount(ups%gboxx(i,j),ups%gboxy(i,j))=tempcount(ups%gboxx(i,j),ups%gboxy(i,j))+mask_fulldomain(i,j)
       enddo
    enddo

    do i=1,grid%nx
       do j=1,grid%ny
          if (tempcount(i,j) == 0) then
             frac_coverage(i,j) = 0.0
          else
             frac_coverage(i,j) = (tempcount(i,j)*lgrid_fulldomain%delta%pt(1)*lgrid_fulldomain%delta%pt(2))/ &
                                  (lon_diff(grid%lon_bound(i+1),grid%lon_bound(i))*D2R*EQ_RAD**2*    &
                                  (sin(grid%lat_bound(j)*D2R)-sin(grid%lat_bound(j+1)*D2R)))
          endif
       enddo
    enddo

    ! Fix points that should be 1.0 by checking their surroundings

    ! Interior points first

    do i=2,grid%nx-1
       do j=2,grid%ny-1
          if ((frac_coverage(i,j)   /= 0).and. &
              (frac_coverage(i+1,j) /= 0).and. &
              (frac_coverage(i,j+1) /= 0).and. &
              (frac_coverage(i-1,j) /= 0).and. &
              (frac_coverage(i,j-1) /= 0)) &
                   frac_coverage(i,j) = 1.0
       enddo
    enddo

    ! top and bottom edges

    do i=2,grid%nx/2
       if ((frac_coverage(i,1)   /= 0).and. &
           (frac_coverage(i+1,1) /= 0).and. &
           (frac_coverage(i,2)   /= 0).and. &
           (frac_coverage(i-1,1) /= 0).and. &
           (frac_coverage(i+grid%nx/2,1) /= 0)) &
                frac_coverage(i,1) = 1.0
    enddo

    do i=grid%nx/2+1,grid%nx-1
       if ((frac_coverage(i,1)   /= 0).and. &
           (frac_coverage(i+1,1) /= 0).and. &
           (frac_coverage(i,2)   /= 0).and. &
           (frac_coverage(i-1,1) /= 0).and. &
           (frac_coverage(i-grid%nx/2,1) /= 0)) &
                frac_coverage(i,1) = 1.0
    enddo

    do i=2,grid%nx/2
       if ((frac_coverage(i,grid%ny)   /= 0).and. &
           (frac_coverage(i+1,grid%ny) /= 0).and. &
           (frac_coverage(i+grid%nx/2,grid%ny) /= 0).and. &
           (frac_coverage(i-1,grid%ny) /= 0).and. &
           (frac_coverage(i,grid%ny-1) /= 0)) &
                frac_coverage(i,grid%ny) = 1.0
    enddo

    do i=grid%nx/2+1,grid%nx-1
       if ((frac_coverage(i,grid%ny) /= 0).and. &
           (frac_coverage(i+1,grid%ny) /= 0).and. &
           (frac_coverage(i-grid%nx/2,grid%ny) /= 0).and. &
           (frac_coverage(i-1,grid%ny) /= 0).and. &
           (frac_coverage(i,grid%ny-1) /= 0)) &
                frac_coverage(i,grid%ny) = 1.0
    enddo

    ! left and right edges

    do j=2,grid%ny-1
       if ((frac_coverage(1,j) /= 0).and. &
           (frac_coverage(2,j) /= 0).and. &
           (frac_coverage(1,j+1) /= 0).and. &
           (frac_coverage(grid%nx,j) /= 0).and. &
           (frac_coverage(1,j-1) /= 0)) &
                frac_coverage(1,j) = 1.0
       if ((frac_coverage(grid%nx,j) /= 0).and. &
           (frac_coverage(1,j) /= 0).and. &
           (frac_coverage(grid%nx,j+1) /= 0).and. &
           (frac_coverage(grid%nx-1,j) /= 0).and. &
           (frac_coverage(grid%nx,j-1) /= 0)) &
                frac_coverage(grid%nx,j) = 1.0
    enddo

    ! corners

    if ((frac_coverage(1,1) /= 0).and. &
        (frac_coverage(2,1) /= 0).and. &
        (frac_coverage(1,2) /= 0).and. &
        (frac_coverage(grid%nx,1) /= 0).and. &
        (frac_coverage(grid%nx/2+1,1) /= 0)) &
             frac_coverage(1,1) = 1.0

    if ((frac_coverage(1,grid%ny) /= 0).and. &
        (frac_coverage(2,grid%ny) /= 0).and. &
        (frac_coverage(grid%nx/2+1,grid%ny) /= 0).and. &
        (frac_coverage(grid%nx,grid%ny) /= 0).and. &
        (frac_coverage(1,grid%ny-1) /= 0)) &
             frac_coverage(1,grid%ny) = 1.0

    if ((frac_coverage(grid%nx,1) /= 0).and. &
        (frac_coverage(1,1) /= 0).and. &
        (frac_coverage(grid%nx,2) /= 0).and. &
        (frac_coverage(grid%nx-1,1) /= 0).and. &
        (frac_coverage(grid%nx/2,1) /= 0)) &
             frac_coverage(grid%nx,1) = 1.0

    if ((frac_coverage(grid%nx,grid%ny) /= 0).and. &
        (frac_coverage(1,grid%ny) /= 0).and. &
        (frac_coverage(grid%nx/2,grid%ny) /= 0).and. &
        (frac_coverage(grid%nx-1,grid%ny) /= 0).and. &
        (frac_coverage(grid%nx,grid%ny-1) /= 0)) &
             frac_coverage(grid%nx,grid%ny) = 1.0

    ! Finally fix any rogue points > 1.0

    where (frac_coverage>1.0) frac_coverage=1.0

  end subroutine calc_coverage

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function lon_diff(a,b)

    implicit none

    real(rk),intent(in) :: a,b
    real(rk) :: aa,bb

    aa=a ; bb=b

    do
       if (aa > bb) exit
       aa = aa + 360.0
    enddo

    lon_diff = aa - bb

  end function lon_diff

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine define_glint_restart_variables(instance)

    ! This subroutine analyzes the glint options input by the user in the config file
    ! and determines which variables are necessary for an exact restart.  MJH 1/11/2013

    ! Please comment thoroughly the reasons why a particular variable needs to be a restart variable for a given config.

    use glint_io, only: glint_add_to_restart_variable_list
    use glint_mbal_io, only: glint_mbal_add_to_restart_variable_list
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glint_instance), intent (in) :: instance  !*FD Derived type that includes all glint options

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    ! Variables needed for restart with glint.
    ! TODO I am simply inserting outmask because it was the only variable with hot=1 in glint_vars.def
    !      Not sure it is needed
    call glint_add_to_restart_variable_list('outmask')

    ! Variables needed for restart with glint_mbal
    ! No variables had hot=1 in glint_mbal_vars.def, so I am not adding any restart variables here.
    ! call glint_mbal_add_to_restart_variable_list('')

  end subroutine define_glint_restart_variables

end module glint_initialise
