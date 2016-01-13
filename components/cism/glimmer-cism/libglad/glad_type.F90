!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_type.F90 - part of the Community Ice Sheet Model (CISM)  
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

#define NCO outfile%nc
#define NCI infile%nc

module glad_type

  !> contains type definitions for GLAD

  use glimmer_global, only: dp
  use glide_types
  use glad_input_averages, only : glad_input_averages_type, initialize_glad_input_averages
  use glad_mbal_coupling, only : glad_mbc
  use glad_output_fluxes, only : glad_output_fluxes_type, initialize_glad_output_fluxes
  
  implicit none

  ! Constants that describe the options available

  ! basic Glad options

  integer, parameter :: EVOLVE_ICE_FALSE = 0   ! do not let the ice sheet evolve
                                               ! (hold the ice state fixed at initial condition)
  integer, parameter :: EVOLVE_ICE_TRUE  = 1   ! let the ice sheet evolve

  integer, parameter :: ZERO_GCM_FLUXES_FALSE = 0 ! send true fluxes to the GCM
  integer, parameter :: ZERO_GCM_FLUXES_TRUE  = 1 ! zero out all fluxes sent to the GCM

  !TODO - Add other Glad options here to avoid hardwiring of case numbers?

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glad_instance

     !> Derived type holding information about ice model instance. 
     !> Note that variables used for downscaling & upscaling are only valid on the main task,
     !> since all downscaling and upscaling is done there.

     type(coordsystem_type)           :: lgrid              !> Local grid for interfacing with glide (grid on this task)
                                                            !> (WJS: Note that origin may be incorrect with multiple tasks;
                                                            !>  as far as I can tell, this isn't currently a problem)
     type(glad_input_averages_type)   :: glad_inputs        !> Time-averaged inputs from the climate model
     type(glide_global_type)          :: model              !> The instance and all its arrays.
     character(fname_length)          :: paramfile          !> The name of the configuration file.
     integer                          :: ice_tstep          !> Ice timestep in hours
     integer                          :: mbal_tstep         !> Mass-balance timestep in hours
     integer                          :: mbal_accum_time    !> Accumulation time for mass-balance (hours)
                                                            !> (defaults to ice time-step)
     integer                          :: ice_tstep_multiply=1 !> Ice time multiplier (non-dimensional)
     integer                          :: n_icetstep         !> Number of ice time-steps per mass-balance accumulation
     real(dp)                         :: glide_time         !> Time as seen by glide (years)
     integer                          :: next_time          !> The next time we expect to be called (hours)

     ! Climate inputs, on the local grid -------------------------

     real(dp),dimension(:,:),pointer :: artm => null() !> Annual mean air temperature
     real(dp),dimension(:,:),pointer :: acab => null() !> Annual mass balance (m/y water equiv)

     ! Arrays to accumulate mass-balance quantities --------------

     type(glad_mbc) :: mbal_accum
     
     ! Climate options -------------------------------------------

     integer :: evolve_ice = 1

     !> Whether the ice sheet can evolve:
     !> \begin{description}
     !> \item[0] The ice sheet cannot evolve; hold fixed at initial state
     !> \item[1] The ice sheet can evolve

     logical :: test_coupling = .false.

     integer :: zero_gcm_fluxes = ZERO_GCM_FLUXES_FALSE
     
     !> Whether to zero out the fluxes (e.g., calving flux) sent to the GCM
     !> \begin{description}
     !> \item[0] send true fluxes to the GCM
     !> \item[1] zero out all fluxes sent to the GCM
     !> \end{description}

     ! Latitude & longitude of model grid points
     real(dp), dimension(:,:), pointer :: lat(:,:) => null()
     real(dp), dimension(:,:), pointer :: lon(:,:) => null()

     ! Fields for averaging dycore output
     type(glad_output_fluxes_type) :: glad_output_fluxes
     real(dp), dimension(:,:), pointer :: hflx_tavg => null()  ! conductive heat flux at top surface (W m-2)
     real(dp), dimension(:,:), pointer :: rofi_tavg => null()  ! solid ice runoff (kg m-2 s-1)
     real(dp), dimension(:,:), pointer :: rofl_tavg => null()  ! liquid runoff from basal/interior melting (kg m-2 s-1)
     
     ! Pointers to file input and output

     type(glimmer_nc_output),pointer :: out_first => null() !> first element of linked list defining netCDF outputs
     type(glimmer_nc_input), pointer :: in_first => null()  !> first element of linked list defining netCDF inputs

  end type glad_instance

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

  subroutine glad_i_allocate_gcm(instance, force_start)

    ! Allocate top-level arrays in the model instance, and ice model arrays.

    implicit none

    type(glad_instance),intent(inout) :: instance    !> Instance whose elements are to be allocated.
    integer,            intent(in)    :: force_start !> glad forcing start time (hours)
    
    integer :: ewn,nsn    ! dimensions of local grid

    ewn = get_ewn(instance%model)
    nsn = get_nsn(instance%model)

    ! First deallocate if necessary

    if (associated(instance%artm))          deallocate(instance%artm)
    if (associated(instance%acab))          deallocate(instance%acab)

    if (associated(instance%lat))           deallocate(instance%lat)
    if (associated(instance%lon))           deallocate(instance%lon)
    
    if (associated(instance%rofi_tavg))     deallocate(instance%rofi_tavg)
    if (associated(instance%rofl_tavg))     deallocate(instance%rofl_tavg)
    if (associated(instance%hflx_tavg))     deallocate(instance%hflx_tavg)


    ! Then reallocate and zero...

    allocate(instance%artm(ewn,nsn));          instance%artm = 0.d0
    allocate(instance%acab(ewn,nsn));          instance%acab = 0.d0

    allocate(instance%lat(ewn,nsn));           instance%lat  = 0.d0
    allocate(instance%lon(ewn,nsn));           instance%lon  = 0.d0

    allocate(instance%rofi_tavg(ewn,nsn));    instance%rofi_tavg = 0.d0
    allocate(instance%rofl_tavg(ewn,nsn));    instance%rofl_tavg = 0.d0
    allocate(instance%hflx_tavg(ewn,nsn));    instance%hflx_tavg = 0.d0
    
    call initialize_glad_input_averages(instance%glad_inputs, ewn=ewn, nsn=nsn, &
         next_av_start=force_start)

    call initialize_glad_output_fluxes(instance%glad_output_fluxes, ewn=ewn, nsn=nsn)
    
  end subroutine glad_i_allocate_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Move the next two subroutines to a new module called glad_setup?
  !       This would be analogous to the organization of Glide.

  subroutine glad_i_readconfig(instance,config)

    !> read glad configuration

    use glimmer_config
    use glimmer_log
    use glad_constants, only: years2hours

    implicit none

    ! Arguments

    type(ConfigSection), pointer       :: config      !> structure holding sections of configuration file   
    type(glad_instance),intent(inout) :: instance    !> The instance being initialised.

    ! Internals

    type(ConfigSection), pointer :: section
    real(dp) :: mbal_time_temp ! Accumulation time in years

    mbal_time_temp = -1.d0

    call GetSection(config,section,'GLAD climate')
    if (associated(section)) then
       call GetValue(section,'evolve_ice',instance%evolve_ice)
       call GetValue(section,'test_coupling',instance%test_coupling)       
       call GetValue(section,'mbal_accum_time',mbal_time_temp)
       call GetValue(section,'ice_tstep_multiply',instance%ice_tstep_multiply)
       call GetValue(section,'zero_gcm_fluxes',instance%zero_gcm_fluxes)
    end if

    if (mbal_time_temp > 0.0) then
       instance%mbal_accum_time = mbal_time_temp * years2hours
    else
       instance%mbal_accum_time = -1
    end if

    call glad_nc_readparams(instance,config)

  end subroutine glad_i_readconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine glad_nc_readparams(instance,config)

    !> read netCDF I/O related configuration file
    !> based on glimmer_ncparams

    use glimmer_config
    use glimmer_ncparams, only: handle_output, handle_input, configstring
    implicit none

    type(glad_instance)         :: instance  !> GLAD instance
    type(ConfigSection), pointer :: config !> structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section
    type(glimmer_nc_output), pointer :: output
    type(glimmer_nc_input), pointer :: input

    ! Initialise local pointers 
    output => null()
    input => null()

    ! setup outputs
    call GetSection(config,section,'GLAD output')
    do while(associated(section))
       output => handle_output(section,output,0.d0,configstring)
       if (.not.associated(instance%out_first)) then
          instance%out_first => output
       end if
       call GetSection(section%next,section,'GLAD output')
    end do

    ! setup inputs
    call GetSection(config,section,'GLAD input')
    do while(associated(section))
       input => handle_input(section,input)
       if (.not.associated(instance%in_first)) then
          instance%in_first => input
       end if
       call GetSection(section%next,section,'GLAD input')
    end do
    
    output => null()
    input => null()

  end subroutine glad_nc_readparams

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glad_i_printconfig(instance)

    use glimmer_log
    use glad_constants, only: hours2years
    use parallel, only: tasks

    implicit none

    ! Argument

    type(glad_instance),intent(inout) :: instance    !> The instance to be printed

    ! Internal

    character(len=100) :: message

    call write_log(' ')
    call write_log('GLAD climate')
    call write_log('-------------')
    write(message,*) 'evolve_ice (0=fixed, 1=evolve):  ',instance%evolve_ice
    call write_log(message)
    write(message,*) 'test_coupling:                   ',instance%test_coupling
    call write_log(message)    

    if (instance%evolve_ice == EVOLVE_ICE_FALSE) then
       call write_log('The ice sheet state will not evolve after initialization')
    endif

    if (instance%mbal_accum_time == -1) then
       call write_log('Mass-balance accumulation time will be set to max(ice timestep, mbal timestep)')
    else
       write(message,*) 'Mass-balance accumulation time:',instance%mbal_accum_time * hours2years,' years'
       call write_log(message)
    end if

    write(message,*) 'ice_tstep_multiply:',instance%ice_tstep_multiply
    call write_log(message)

    write(message,*) 'zero_gcm_fluxes: ', instance%zero_gcm_fluxes
    call write_log(message)

  end subroutine glad_i_printconfig

end module glad_type
