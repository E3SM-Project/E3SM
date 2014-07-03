!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_type.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#define NCO outfile%nc
#define NCI infile%nc

module glint_type

  !*FD contains type definitions for GLINT

  use glimmer_global
  use glint_interp
  use glide_types
  use glint_mbal_coupling, only: glint_mbc, mbal_has_snow_model

  implicit none

  ! Constants that describe the options available

  ! basic Glint options

  integer, parameter :: EVOLVE_ICE_FALSE = 0   ! do not let the ice sheet evolve
                                               ! (hold the ice state fixed at initial condition)
  integer, parameter :: EVOLVE_ICE_TRUE  = 1   ! let the ice sheet evolve

  !TODO - Add other Glint options here to avoid hardwiring of case numbers

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!TODO - glint_instance includes information that is not needed if the SMB is received from a GCM.
!       Maybe we should create a new derived type (glint_instance_gcm?) without the extra information.

  type glint_instance

     !*FD Derived type holding information about ice model instance. 
     !*FD Note that variables used for downscaling & upscaling are only valid on the main task,
     !*FD since all downscaling and upscaling is done there.

     type(coordsystem_type)           :: lgrid              !*FD Local grid for interfacing with glide (grid on this task)
                                                            !*FD (WJS: Note that origin may be incorrect with multiple tasks;
                                                            !*FD  as far as I can tell, this isn't currently a problem)
     type(coordsystem_type)           :: lgrid_fulldomain   !*FD Local grid on the full domain (across all tasks),
                                                            !*FD used for downscaling & upscaling
                                                            !*FD (ONLY VALID ON MAIN TASK)
     type(downscale)                  :: downs              !*FD Downscaling parameters.
                                                            !*FD (ONLY VALID ON MAIN TASK)
     type(upscale)                    :: ups                !*FD Upscaling parameters
                                                            !*FD (ONLY VALID ON MAIN TASK)
     type(upscale)                    :: ups_orog           !*FD Upscaling parameters for orography (to cope
                                                            !*FD with need to convert to spectral form).
                                                            !*FD (ONLY VALID ON MAIN TASK)
     type(glide_global_type)          :: model              !*FD The instance and all its arrays.
     character(fname_length)          :: paramfile          !*FD The name of the configuration file.
     integer                          :: ice_tstep          !*FD Ice timestep in hours
     integer                          :: mbal_tstep         !*FD Mass-balance timestep in hours
     integer                          :: mbal_accum_time    !*FD Accumulation time for mass-balance (hours)
                                                            !*FD (defaults to ice time-step)
     integer                          :: ice_tstep_multiply=1 !*FD Ice time multiplier (non-dimensional)
     integer                          :: n_icetstep         !*FD Number of ice time-steps per mass-balance accumulation
     real(dp)                         :: glide_time         !*FD Time as seen by glide (years)
     integer                          :: next_time          !*FD The next time we expect to be called (hours)

     ! Climate inputs from global model --------------------------

     !TODO - Change to dp
     real(sp),dimension(:,:),pointer :: artm        => null() !*FD Annual mean air temperature
     real(sp),dimension(:,:),pointer :: arng        => null() !*FD Annual air temperature half-range
     real(sp),dimension(:,:),pointer :: prcp        => null() !*FD Precipitation (mm or m)
     real(sp),dimension(:,:),pointer :: snowd       => null() !*FD Snow depth (m)
     real(sp),dimension(:,:),pointer :: siced       => null() !*FD Superimposed ice depth (m)
     real(rk),dimension(:,:),pointer :: xwind       => null() !*FD $x$-component of surface winds (m/s)
     real(rk),dimension(:,:),pointer :: ywind       => null() !*FD $y$-component of surface winds (m/s)
     real(rk),dimension(:,:),pointer :: humid       => null() !*FD Surface humidity (%)
     real(rk),dimension(:,:),pointer :: lwdown      => null() !*FD Downwelling longwave (W/m^2)
     real(rk),dimension(:,:),pointer :: swdown      => null() !*FD Downwelling shortwave (W/m^2)
     real(rk),dimension(:,:),pointer :: airpress    => null() !*FD Surface air pressure (Pa)
     real(dp),dimension(:,:),pointer :: global_orog => null() !*FD Global orography (m)
     real(sp),dimension(:,:),pointer :: local_orog  => null() !*FD Local orography (m)

     ! Locally calculated climate/mass-balance fields ------------

     real(sp),dimension(:,:),pointer :: ablt => null() !*FD Annual ablation
     real(sp),dimension(:,:),pointer :: acab => null() !*FD Annual mass-balance

     ! Arrays to accumulate mass-balance quantities --------------

     type(glint_mbc) :: mbal_accum

     ! Fractional coverage information ---------------------------

     real(rk) ,dimension(:,:),pointer :: frac_coverage => null() 
     !*FD Fractional coverage of each global gridbox by the projected grid.
     !*FD (ONLY VALID ON MAIN TASK)

     real(rk) ,dimension(:,:),pointer :: frac_cov_orog => null() 
     !*FD Fractional coverage of each global gridbox by the projected grid (orography).
     !*FD (ONLY VALID ON MAIN TASK)

     ! Output masking --------------------------------------------

     integer, dimension(:,:),pointer :: out_mask => null() 

     !*FD Array indicating whether a point should be considered or ignored 
     !*FD when upscaling data for output. 1 means use, 0 means ignore.

     ! Climate options -------------------------------------------

     integer :: evolve_ice = 1

     !*FD Whether the ice sheet can evolve:
     !*FD \begin{description}
     !*FD \item[0] The ice sheet cannot evolve; hold fixed at initial state
     !*FD \item[1] The ice sheet can evolve

     integer :: whichacab = 1

     !*FD Which mass-balance scheme: 
     !*FD \begin{description}
     !*FD \item[0] Receive surface mass balance from climate model
     !*FD \item[1] PDD mass-balance model
     !*FD \item[2] Accumulation only 
     !*FD \item[3] RAPID energy balance model
     !*FD \end{description}

     integer :: whichprecip = 1

     !*FD Source of precipitation:
     !*FD \begin{description}
     !*FD \item[1] Use large-scale precip as is.
     !*FD \item[2] Use parameterization of \emph{Roe and Lindzen} 
     !*FD \end{description}

     integer :: use_mpint = 0
   
     !*FD Flag to control if mean-preserving interpolation is used

     ! Climate parameters ----------------------------------------------------------

     !TODO - Change to dp
     real(sp) :: ice_albedo   =   0.4 !*FD Ice albedo. (fraction)
     real(sp) :: lapse_rate   =   8.0 !*FD Uniform lapse rate in deg C/km 
     !*FD (N.B. This should be \emph{positive} for temperature falling with height!)
     real(sp) :: data_lapse_rate = 8.0 !*FD Implied lapse rate in large-scale data (used for
     !*FD tuning). Set equal to lapse\_rate if not supplied.

     ! Counter for averaging temperature input --------------------------------------

     integer  :: av_count = 0 !*FD Counter for averaging temperature input

     ! Pointers to file input and output

     type(glimmer_nc_output),pointer :: out_first => null() !*FD first element of linked list defining netCDF outputs
     type(glimmer_nc_input), pointer :: in_first => null()  !*FD first element of linked list defining netCDF inputs

  end type glint_instance

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type output_flags
     !*FD A derived type used internally to communicate the outputs which need
     !*FD to be upscaled, thus avoiding unnecessary calculation

     logical :: orog         !*FD Set if we need to upscale the orography
     logical :: albedo       !*FD Set if we need to upscale the albedo
     logical :: ice_frac     !*FD Set if we need to upscale the ice fraction
     logical :: veg_frac     !*FD Set if we need to upscale the veg fraction
     logical :: snowice_frac !*FD Set if we need to upscale the snow-covered ice fraction
     logical :: snowveg_frac !*FD Set if we need to upscale the snow-covered veg fraction
     logical :: snow_depth   !*FD Set if we need to upscale the snow depth
     logical :: water_in     !*FD Set if we need to upscale the input water flux
     logical :: water_out    !*FD Set if we need to upscale the output water flux
     logical :: total_win    !*FD Set if we need to sum the total water taken up by ice sheet
     logical :: total_wout   !*FD Set if we need to sum the total ablation by the ice sheet
     logical :: ice_vol      !*FD Set if we need to calculate the total ice volume
  end type output_flags


  !TODO - These diagnostic points are grid-specific.  
  !       For GCM coupling, they should be set based on the global GCM grid.

  ! diagnostic points on global grid, useful for debugging

  integer, parameter :: iglint_global = 56     ! SW Greenland point on 64 x 32 glint_example global grid
  integer, parameter :: jglint_global = 4      ! j increases from north to south

contains

  subroutine glint_i_allocate(instance,nxg,nyg,nxgo,nygo)

    !*FD Allocate top-level arrays in the model instance, and ice model arrays.

    implicit none

    type(glint_instance),intent(inout) :: instance  !*FD Instance whose elements are to be allocated.
    integer,             intent(in)    :: nxg       !*FD Longitudinal size of global grid (grid-points).
    integer,             intent(in)    :: nyg       !*FD Latitudinal size of global grid (grid-points).
    integer,             intent(in)    :: nxgo      !*FD Longitudinal size of global orog grid (grid-points).
    integer,             intent(in)    :: nygo      !*FD Latitudinal size of global orog grid (grid-points).

    integer ewn,nsn

    ewn=get_ewn(instance%model)
    nsn=get_nsn(instance%model)

    ! First deallocate if necessary
    ! Downscaled global arrays

    if (associated(instance%artm))          deallocate(instance%artm)
    if (associated(instance%arng))          deallocate(instance%arng)
    if (associated(instance%prcp))          deallocate(instance%prcp)
    if (associated(instance%snowd))         deallocate(instance%snowd)
    if (associated(instance%siced))         deallocate(instance%siced)
    if (associated(instance%xwind))         deallocate(instance%xwind)
    if (associated(instance%ywind))         deallocate(instance%ywind)
    if (associated(instance%humid))         deallocate(instance%humid)
    if (associated(instance%lwdown))        deallocate(instance%lwdown)
    if (associated(instance%swdown))        deallocate(instance%swdown)
    if (associated(instance%airpress))      deallocate(instance%airpress)
    if (associated(instance%global_orog))   deallocate(instance%global_orog) 
    if (associated(instance%local_orog))    deallocate(instance%local_orog)

    ! Local climate arrays

    if (associated(instance%ablt))          deallocate(instance%ablt)
    if (associated(instance%acab))          deallocate(instance%acab)

    ! Fractional coverage

    if (associated(instance%frac_coverage)) deallocate(instance%frac_coverage)
    if (associated(instance%frac_cov_orog)) deallocate(instance%frac_cov_orog)

    ! Output mask

    if (associated(instance%out_mask))      deallocate(instance%out_mask)

    ! Then reallocate and zero...
    ! Global input fields

    allocate(instance%artm(ewn,nsn));        instance%artm        = 0.0
    allocate(instance%arng(ewn,nsn));        instance%arng        = 0.0
    allocate(instance%prcp(ewn,nsn));        instance%prcp        = 0.0
    allocate(instance%snowd(ewn,nsn));       instance%snowd       = 0.0
    allocate(instance%siced(ewn,nsn));       instance%siced       = 0.0
    allocate(instance%xwind(ewn,nsn));       instance%xwind       = 0.0
    allocate(instance%ywind(ewn,nsn));       instance%ywind       = 0.0
    allocate(instance%humid(ewn,nsn));       instance%humid       = 0.0
    allocate(instance%lwdown(ewn,nsn));      instance%lwdown      = 0.0
    allocate(instance%swdown(ewn,nsn));      instance%swdown      = 0.0
    allocate(instance%airpress(ewn,nsn));    instance%airpress    = 0.0
    allocate(instance%global_orog(ewn,nsn)); instance%global_orog = 0.0
    allocate(instance%local_orog(ewn,nsn));  instance%local_orog  = 0.0

    ! Local fields

    allocate(instance%ablt(ewn,nsn)); instance%ablt = 0.0
    allocate(instance%acab(ewn,nsn)); instance%acab = 0.0

    ! Fractional coverage map

    allocate(instance%frac_coverage(nxg,nyg)); instance%frac_coverage = 0.0
    allocate(instance%frac_cov_orog(nxgo,nygo)); instance%frac_cov_orog = 0.0

    ! Output mask

    allocate(instance%out_mask(ewn,nsn)); instance%out_mask = 1

  end subroutine glint_i_allocate

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_allocate_gcm(instance,nxg,nyg)

    ! Allocate top-level arrays in the model instance, and ice model arrays.

    implicit none

    type(glint_instance),intent(inout) :: instance  !*FD Instance whose elements are to be allocated.
    integer,             intent(in)    :: nxg       !*FD Longitudinal size of global grid (grid-points).
    integer,             intent(in)    :: nyg       !*FD Latitudinal size of global grid (grid-points).

    integer ewn,nsn

    ewn = get_ewn(instance%model)
    nsn = get_nsn(instance%model)

    ! First deallocate if necessary

    if (associated(instance%artm))          deallocate(instance%artm)
    if (associated(instance%acab))          deallocate(instance%acab)
    if (associated(instance%frac_coverage)) deallocate(instance%frac_coverage)
    if (associated(instance%out_mask))      deallocate(instance%out_mask)  !TODO - Is this needed?

    ! Then reallocate and zero...

    allocate(instance%artm(ewn,nsn));          instance%artm = 0.0
    allocate(instance%acab(ewn,nsn));          instance%acab = 0.0
    allocate(instance%frac_coverage(nxg,nyg)); instance%frac_coverage = 0.0
    allocate(instance%out_mask(ewn,nsn));      instance%out_mask = 1   !TODO - Is this needed?

  end subroutine glint_i_allocate_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Should the next two subroutines go in a new module called glint_setup?
  !       This would be analogous to the organization of Glide.

  subroutine glint_i_readconfig(instance,config)

    !*FD read glint configuration

    use glimmer_config
    use glimmer_log
    use glint_constants, only: years2hours

    implicit none

    ! Arguments

    type(ConfigSection), pointer       :: config      !*FD structure holding sections of configuration file   
    type(glint_instance),intent(inout) :: instance    !*FD The instance being initialised.

    ! Internals

    type(ConfigSection), pointer :: section
    real(dp) :: mbal_time_temp ! Accumulation time in years

    mbal_time_temp = -1.d0

    call GetSection(config,section,'GLINT climate')
    if (associated(section)) then
       call GetValue(section,'evolve_ice',instance%evolve_ice)
       call GetValue(section,'precip_mode',instance%whichprecip)
       call GetValue(section,'acab_mode',instance%whichacab)
       call GetValue(section,'ice_albedo',instance%ice_albedo)
       call GetValue(section,'lapse_rate',instance%lapse_rate)
       instance%data_lapse_rate=instance%lapse_rate
       call GetValue(section,'data_lapse_rate',instance%data_lapse_rate)
       call GetValue(section,'mbal_accum_time',mbal_time_temp)
       call GetValue(section,'ice_tstep_multiply',instance%ice_tstep_multiply)
       call GetValue(section,'mean_preserving',instance%use_mpint)
    end if

    if (mbal_time_temp > 0.0) then
       instance%mbal_accum_time = mbal_time_temp * years2hours
    else
       instance%mbal_accum_time = -1
    end if

    call glint_nc_readparams(instance,config)

  end subroutine glint_i_readconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine glint_nc_readparams(instance,config)

    !*FD read netCDF I/O related configuration file
    !*FD based on glimmer_ncparams

    use glimmer_config
    use glimmer_ncparams, only: handle_output, handle_input, configstring
    implicit none

    type(glint_instance)         :: instance  !*FD GLINT instance
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section
    type(glimmer_nc_output), pointer :: output
    type(glimmer_nc_input), pointer :: input

    ! Initialise local pointers 
    output => null()
    input => null()

    ! setup outputs
    call GetSection(config,section,'GLINT output')
    do while(associated(section))
       output => handle_output(section,output,0.d0,configstring)
       if (.not.associated(instance%out_first)) then
          instance%out_first => output
       end if
       call GetSection(section%next,section,'GLINT output')
    end do

    ! setup inputs
    call GetSection(config,section,'GLINT input')
    do while(associated(section))
       input => handle_input(section,input)
       if (.not.associated(instance%in_first)) then
          instance%in_first => input
       end if
       call GetSection(section%next,section,'GLINT input')
    end do
    
    output => null()
    input => null()

  end subroutine glint_nc_readparams

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_printconfig(instance)

    use glimmer_log
    use glint_constants, only: hours2years

    implicit none

    ! Argument

    type(glint_instance),intent(inout) :: instance    !*FD The instance to be printed

    ! Internal

    character(len=100) :: message

    !TODO - Make these messages more informative 
    !       (Print a description instead of just the option number.)
    !TODO - Remove hardwired option numbers

    call write_log(' ')
    call write_log('GLINT climate')
    call write_log('-------------')
    write(message,*) 'evolve_ice  ',instance%evolve_ice
    call write_log(message)
    write(message,*) 'precip_mode ',instance%whichprecip
    call write_log(message)
    write(message,*) 'acab_mode   ',instance%whichacab
    call write_log(message)

    if (instance%evolve_ice == EVOLVE_ICE_FALSE) then
       call write_log('The ice sheet state will not evolve after initialization')
    endif

    if (instance%whichacab /= 0) then  ! not getting SMB from GCM
       write(message,*) 'ice_albedo  ',instance%ice_albedo
       call write_log(message)
       write(message,*) 'lapse_rate  ',instance%lapse_rate
       call write_log(message)
       write(message,*) 'data_lapse_rate',instance%data_lapse_rate
       call write_log(message)
    endif

    if (instance%mbal_accum_time == -1) then
       call write_log('Mass-balance accumulation time == ice time-step')
    else
       write(message,*) 'Mass-balance accumulation time:',instance%mbal_accum_time * hours2years,' years'
       call write_log(message)
    end if

    write(message,*) 'ice_tstep_multiply',instance%ice_tstep_multiply
    call write_log(message)

    select case(instance%use_mpint)
    case(1)
       write(message,*) 'Using mean-preserving interpolation'
       call write_log(message)
    case(0)
       write(message,*) 'Using normal interpolation'
       call write_log(message)
    case default
       write(message,*) 'Unrecognised value of instance%use_mpint'
       call write_log(message,GM_FATAL)
    end select

  end subroutine glint_i_printconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Can we remove this function and only use mbal_has_snow_model?

  logical function glint_has_snow_model(instance)

    use glint_mbal, only: mbal_has_snow_model

    type(glint_instance),            intent(in)  :: instance

    glint_has_snow_model = mbal_has_snow_model(instance%mbal_accum%mbal)

  end function glint_has_snow_model

end module glint_type
