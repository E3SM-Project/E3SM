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
!   glad_initialise.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
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

module glad_initialise

  !> Initialise GLAD model instance

  use glad_type
  use glimmer_global, only: dp
  implicit none

  private
  public glad_i_initialise_gcm, glad_i_end

contains

  subroutine glad_i_initialise_gcm(config,           instance,         &
                                    force_start,      force_dt,         &
                                    gcm_restart,      gcm_restart_file, &
                                    gcm_config_unit)

    ! Initialise a GLAD ice model instance for GCM coupling

    use glimmer_paramets, only: GLC_DEBUG
    use glimmer_log
    use glimmer_config
    use glimmer_coordinates, only : coordsystem_new
    use glad_mbal_coupling, only : glad_mbc_init
    use glad_io          , only: glad_io_createall     , glad_io_writeall
    use glad_mbal_io     , only: glad_mbal_io_createall, glad_mbal_io_writeall
    use glimmer_ncio
    use glide_nc_custom   , only: glide_nc_fillall
    use glide
    use glissade
    use glad_constants
    use glad_restart_gcm
    use glide_diagnostics
    use parallel, only: main_task

    implicit none

    ! Arguments
    type(ConfigSection), pointer         :: config           ! structure holding sections of configuration file   
    type(glad_instance),  intent(inout) :: instance         ! The instance being initialised.

    integer,               intent(in)    :: force_start      ! glad forcing start time (hours)
    integer,               intent(in)    :: force_dt         ! glad forcing time step (hours)

    logical,     optional, intent(in)    :: gcm_restart      ! logical flag to read from a restart file
    character(*),optional, intent(in)    :: gcm_restart_file ! restart filename for restart
    integer,     optional, intent(in)    :: gcm_config_unit  ! fileunit for reading config files

    ! Internal

    integer :: config_fileunit

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
            call glad_read_restart_gcm(instance%model, gcm_restart_file)
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

       ! compute the initial diagnostic state
       call glide_init_state_diagnostic(instance%model)

    else       ! glam/glissade HO dycore     

       ! initialise the model
       call glissade_initialise(instance%model)

       ! compute the initial diagnostic state
       call glissade_diagnostic_variable_solve(instance%model)

    endif

    instance%ice_tstep = get_tinc(instance%model)*nint(years2hours)

    instance%glide_time = instance%model%numerics%tstart

    ! read glad configuration

    call glad_i_readconfig(instance, config)    
    call glad_i_printconfig(instance)    

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file (specific to glad - other configs,
    ! e.g. glide, isos, are handled separately by their setup routines).
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/17/13
    ! Note: the corresponding call for glide is placed within *_readconfig, which is probably more appropriate,
    ! but putting this call into glad_i_readconfig creates a circular dependency.  

    call define_glad_restart_variables(instance)
 
    ! create glad variables for the glide output files
    call glad_io_createall(instance%model, data=instance)

    ! create instantaneous glad variables
    call openall_out(instance%model, outfiles=instance%out_first)
    call glad_mbal_io_createall(instance%model, data=instance, outfiles=instance%out_first)

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

    call glad_i_allocate_gcm(instance, force_start)

    ! Read data and initialise climate

    call glad_i_readdata(instance)

    ! initialise the mass-balance accumulation

    call glad_mbc_init(instance%mbal_accum, instance%lgrid)

    ! If flag set to force frequent coupling (for testing purposes),
    ! then decrease all coupling timesteps to very short intervals
    if (instance%test_coupling) then
       instance%mbal_accum%tstep = 24
       instance%mbal_accum_time =  24
       instance%ice_tstep =        24
    endif

    instance%mbal_tstep = instance%mbal_accum%tstep

    instance%next_time = force_start - force_dt + instance%mbal_tstep

    if (GLC_DEBUG .and. main_task) then
       write (6,*) 'Called glad_mbc_init'
       write (6,*) 'mbal tstep =', instance%mbal_tstep
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

    if (instance%ice_tstep_multiply/=1 .and. mod(instance%mbal_accum_time,nint(years2hours)) /= 0.d0) then
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
    call glad_io_writeall(instance, instance%model)
    call glad_mbal_io_writeall(instance, instance%model, outfiles=instance%out_first)

  end subroutine glad_i_initialise_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glad_i_end(instance)

    !> Tidy up 

    use glide
    use glimmer_ncio
    implicit none
    type(glad_instance),  intent(inout) :: instance    !> The instance being initialised.

    call glide_finalise(instance%model)
    call closeall_out(instance%model,outfiles=instance%out_first)
    instance%out_first => null()

  end subroutine glad_i_end

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glad_i_readdata(instance)
    !> read data from netCDF file and initialise climate

    use glad_io
    use glide_thck, only: glide_calclsrf
    implicit none

    type(glad_instance),intent(inout)   :: instance    !> Instance whose elements are to be allocated.

    ! read data
    call glad_io_readall(instance,instance%model)

    call glide_calclsrf(instance%model%geometry%thck,instance%model%geometry%topg, &
         instance%model%climate%eus,instance%model%geometry%lsrf)
    instance%model%geometry%usrf = instance%model%geometry%thck + instance%model%geometry%lsrf

  end subroutine glad_i_readdata

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine define_glad_restart_variables(instance)

    ! This subroutine analyzes the glad options input by the user in the config file
    ! and determines which variables are necessary for an exact restart.  MJH 1/11/2013

    ! Please comment thoroughly the reasons why a particular variable needs to be a restart variable for a given config.

    use glad_io, only: glad_add_to_restart_variable_list
    use glad_mbal_io, only: glad_mbal_add_to_restart_variable_list
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glad_instance), intent (in) :: instance  !> Derived type that includes all glad options

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    ! lat and lon need to be on the input file. Since a restart run only reads the
    ! restart file (and not the original input file) we need to write lat and lon back to
    ! the restart file so they will be available for the following run segment.
    
    call glad_add_to_restart_variable_list('lat lon')
    
    ! The variables rofi_tavg, rofl_tavg, and hflx_tavg are time-averaged fluxes on the local grid
    !  from the previous coupling interval. They are included here so that the coupler can be sent
    !  the correct fluxes after restart; otherwise these fluxes would have values of zero.
    !TODO - Add av_count_output so we can restart in the middle of a mass balance timestep?
   
    call glad_add_to_restart_variable_list('rofi_tavg rofl_tavg hflx_tavg')

  end subroutine define_glad_restart_variables


end module glad_initialise
