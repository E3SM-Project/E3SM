!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_setup.F90 - part of the Community Ice Sheet Model (CISM)  
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

#include "glide_mask.inc"

module glide_setup

  ! general routines for initialisation, etc, called from top-level glimmer subroutines

  use glimmer_global, only: dp

  implicit none

  private
  public :: glide_readconfig, glide_printconfig, glide_scale_params, &
            glide_load_sigma, glide_read_sigma, glide_calc_sigma

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine glide_readconfig(model,config)

    ! read GLIDE configuration file
    ! Note: sigma coordinates are handled by a subsequent call to glide_read_sigma
 
    use glide_types
    use glimmer_config
    implicit none
    type(glide_global_type) :: model        !> model instance
    type(ConfigSection), pointer :: config  !> structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section

    ! read grid size  parameters
    call GetSection(config,section,'grid')
    if (associated(section)) then
       call handle_grid(section, model)
    end if

    ! read time parameters
    call GetSection(config,section,'time')
    if (associated(section)) then
       call handle_time(section, model)
    end if

    ! read options parameters
    call GetSection(config,section,'options')
    if (associated(section)) then
       call handle_options(section, model)
    end if

    !read options for higher-order computation
    call GetSection(config,section,'ho_options')
    if (associated(section)) then
        call handle_ho_options(section, model)
    end if

     !read options for computation using an external dycore -- Doug Ranken 04/20/12
    call GetSection(config,section,'external_dycore_options')
    if (associated(section)) then
        call handle_dycore_options(section, model)
    end if

    ! read parameters
    call GetSection(config,section,'parameters')
    if (associated(section)) then
       call handle_parameters(section, model)
    end if

    ! read GTHF 
    ! NOTE: The [GTHF] section is ignored unless model%options%gthf = GTHF_COMPUTE
    if (model%options%gthf == GTHF_COMPUTE) then
       call GetSection(config,section,'GTHF')
       if (associated(section)) then
          call handle_gthf(section, model)
       end if
    endif

    ! read isostasy
    ! NOTE: The [isostasy] section is ignored unless model%options%isostasy = ISOSTASY_COMPUTE
    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call GetSection(config,section,'isostasy')
       if (associated(section)) then
          call handle_isostasy(section, model)
       end if
    endif

    ! Till options are not currently supported
    ! read till parameters
!!    call GetSection(config,section,'till_options')
!!    if (associated(section)) then
!!       call handle_till_options(section, model)
!!    end if

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file.
    ! (Glint restart variables are handled separately by Glint setup routines.)
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/17/13

    call define_glide_restart_variables(model%options)

  end subroutine glide_readconfig

!-------------------------------------------------------------------------

  subroutine glide_printconfig(model)

    !> print model configuration to log
    use glimmer_log
    use glide_types
    implicit none
    type(glide_global_type)  :: model !> model instance

    call write_log_div
    call print_grid(model)
    call print_time(model)
    call print_options(model)
    call print_parameters(model)
    call print_gthf(model)
    call print_isostasy(model)
!!    call print_till_options(model)  ! disabled for now

  end subroutine glide_printconfig

!-------------------------------------------------------------------------
    
  subroutine glide_scale_params(model)
    !> scale parameters
    use glide_types
    use glimmer_physcon,  only: scyr, gn
    use glimmer_paramets, only: thk0, tim0, len0, vel0, vis0, acc0, tau0

    implicit none

    type(glide_global_type)  :: model !> model instance

    model%numerics%dttem = model%numerics%ntem * model%numerics%tinc   

    ! convert dt and dttem from yr to scaled time units
    model%numerics%dt     = model%numerics%tinc * scyr / tim0   
    model%numerics%dttem  = model%numerics%dttem * scyr / tim0   

    ! allow for subcycling of ice transport
    model%numerics%dt_transport = model%numerics%dt / real(model%numerics%subcyc, dp)

    model%numerics%thklim = model%numerics%thklim  / thk0       
    model%numerics%thklim_temp = model%numerics%thklim_temp  / thk0

    model%numerics%dew = model%numerics%dew / len0
    model%numerics%dns = model%numerics%dns / len0

    !TODO - Scale eus for calving?
    !       Currently the scaling for eus (like relx and topg) is handled automatically in glide_io.F90.
    !       Would need to handle eus scaling separately if reading from config file.

    model%calving%marine_limit = model%calving%marine_limit / thk0
    model%calving%calving_minthck = model%calving%calving_minthck / thk0
    model%calving%calving_timescale = model%calving%calving_timescale * scyr / tim0

    model%numerics%periodic_offset_ew = model%numerics%periodic_offset_ew / thk0
    model%numerics%periodic_offset_ns = model%numerics%periodic_offset_ns / thk0

    model%velowk%trc0   = vel0 * len0 / (thk0**2)
    model%velowk%btrac_const = model%paramets%btrac_const/model%velowk%trc0/scyr
    model%velowk%btrac_max   = model%paramets%btrac_max / model%velowk%trc0/scyr    
    model%velowk%btrac_slope = model%paramets%btrac_slope*acc0/model%velowk%trc0

    ! scale beta parameters
    model%velocity%ho_beta_const = model%velocity%ho_beta_const / (tau0/(vel0*scyr))
    model%velocity%beta_grounded_min = model%velocity%beta_grounded_min / (tau0/(vel0*scyr))

    ! scale basal melting parameters (yr^{-1} -> s^{-1})
    model%temper%bmlt_float_omega = model%temper%bmlt_float_omega / scyr
    model%temper%bmlt_float_rate  = model%temper%bmlt_float_rate / scyr

  end subroutine glide_scale_params

!-------------------------------------------------------------------------

  subroutine glide_read_sigma(model,config)

    ! read sigma levels from configuration file, if present
    ! called immediately after glide_readconfig

    use glide_types
    use glimmer_config
    use glimmer_log
    implicit none

    type(glide_global_type) :: model        !> model instance
    type(ConfigSection), pointer :: config  !> structure holding sections of configuration file
        
    ! local variables
    type(ConfigSection), pointer :: section

    ! read sigma levels
    ! NOTE: The [sigma] section is ignored unless model%options%which_sigma = SIGMA_CONFIG

    if (model%options%which_sigma == SIGMA_CONFIG) then
       call GetSection(config,section,'sigma')
       if (associated(section)) then
          call handle_sigma(section, model)
       else
          model%options%which_sigma = SIGMA_COMPUTE_GLIDE  ! default to standard sigma levels
          call write_log('No [sigma] section present; will compute standard Glide sigma levels')
       end if
    endif

  end subroutine glide_read_sigma

!-------------------------------------------------------------------------

  subroutine glide_load_sigma(model,unit)

    ! Compute sigma coordinates or read them from a file
    ! Note: This subroutine is called from glide_initialise or glissade_initialise.
    !       If sigma levels are provided in the config file, then they are read
    !        in by glide_read_sigma, and model%options%which_sigma is set to
    !        SIGMA_CONFIG, in which case this subroutine does nothing.

    use glide_types
    use glimmer_log
    use glimmer_filenames
    use parallel

    implicit none

    ! Arguments
    type(glide_global_type),intent(inout) :: model !> Ice model to use
    integer,               intent(in)    :: unit   !> Logical file unit to use. 
                                                   !> (Must not already be in use)

    ! Internal variables

    integer :: up,upn
    logical :: there
    real(dp) :: level

    ! Beginning of code

    upn=model%general%upn

    select case(model%options%which_sigma)

    case(SIGMA_COMPUTE_GLIDE)   !  compute standard Glide sigma levels

       do up = 1,upn
          level = real(up-1,kind=dp) / real(upn-1,kind=dp)
          model%numerics%sigma(up) = glide_calc_sigma(level, 2.d0)
       end do

       call write_log('Computing Glide sigma levels')

    case(SIGMA_EXTERNAL)        ! read from external file

       if (main_task) inquire (exist=there, file=process_path(model%funits%sigfile))
       call broadcast(there)
       if (.not.there) then
          call write_log('Sigma levels file: '//trim(process_path(model%funits%sigfile))// &
               ' does not exist',GM_FATAL)
       end if
       call write_log('Reading sigma file: '//process_path(model%funits%sigfile))
       if (main_task) then
          open(unit,file=process_path(model%funits%sigfile))
          read(unit,'(f9.7)',err=10,end=10) (model%numerics%sigma(up), up=1,upn)
          close(unit)
       end if
       call broadcast(model%numerics%sigma)

    case(SIGMA_CONFIG)          ! read from config file

       ! sigma levels have already been read from glide_read_sigma

       call write_log('Getting sigma levels from configuration file')

    case(SIGMA_COMPUTE_EVEN)

       do up = 1,upn
          model%numerics%sigma(up) = real(up-1,kind=dp) / real(upn-1,kind=dp)
       enddo

       call write_log('Computing evenly spaced sigma levels')

    case(SIGMA_COMPUTE_PATTYN)

       do up = 1,upn
          if (up == 1) then
             model%numerics%sigma(up) = 0.d0
          else if (up == upn) then
             model%numerics%sigma(up) = 1.d0
          else
             level = real(up-1,kind=dp) / real(upn-1,kind=dp)
             model%numerics%sigma(up) = glide_calc_sigma_pattyn(level)
          end if
       enddo

       call write_log('Computing Pattyn sigma levels')

    end select


    !NOTE: Glam will always use evenly spaced levels,
    !      overriding other values of which_sigma 
    !      (including sigma levels in config file)

    if (model%options%whichdycore == DYCORE_GLAM) then   ! evenly spaced levels are required

       do up = 1,upn
          model%numerics%sigma(up) = real(up-1,kind=dp) / real(upn-1,kind=dp)
       enddo

       call write_log('Using evenly spaced sigma levels for Glam as required')

    endif

    ! Compute stagsigma (= sigma values at layers midpoints)

    model%numerics%stagsigma(1:upn-1) =   &
            (model%numerics%sigma(1:upn-1) + model%numerics%sigma(2:upn)) / 2.0_dp

    ! Compute stagwbndsigma, adding the boundaries to stagsigma

    model%numerics%stagwbndsigma(1:upn-1) = model%numerics%stagsigma(1:upn-1)
    model%numerics%stagwbndsigma(0) = 0.d0
    model%numerics%stagwbndsigma(upn) = 1.d0        

    call print_sigma(model)

    return

10  call write_log('something wrong with sigma coord file',GM_FATAL)
    
  end subroutine glide_load_sigma

!--------------------------------------------------------------------------------

  function glide_calc_sigma(x,n)

     implicit none
     real(dp) :: glide_calc_sigma, x, n
      
     glide_calc_sigma = (1-(x+1)**(-n)) / (1-2**(-n))

  end function glide_calc_sigma

!--------------------------------------------------------------------------------

  function glide_calc_sigma_pattyn(x)

     ! Implements an alternate set of sigma levels that encourages better
     ! convergence for higher-order velocities

     implicit none
     real(dp) :: glide_calc_sigma_pattyn, x

     glide_calc_sigma_pattyn =   &
         (-2.5641025641d-4)*(41d0*x)**2+3.5256410256d-2*(41d0*x)-8.0047080075d-13

  end function glide_calc_sigma_pattyn

!--------------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private procedures
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! grid sizes

  subroutine handle_grid(section, model)
    use glimmer_config
    use glide_types
    use glimmer_filenames
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'ewn',model%general%ewn)
    call GetValue(section,'nsn',model%general%nsn)
    call GetValue(section,'upn',model%general%upn)
    call GetValue(section,'dew',model%numerics%dew)
    call GetValue(section,'dns',model%numerics%dns)
    call GetValue(section,'sigma_file',model%funits%sigfile)
    call GetValue(section,'global_bc',model%general%global_bc)

    ! We set this flag to one to indicate we've got a sigfile name.
    ! A warning/error is generated if sigma levels are specified in some other way
    ! and mangle the name
    if (trim(model%funits%sigfile) /= '') then
       model%funits%sigfile = filenames_inputname(model%funits%sigfile)
       model%options%which_sigma = SIGMA_EXTERNAL
    end if

  end subroutine handle_grid

!--------------------------------------------------------------------------------

  subroutine print_grid(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=512) :: message


    call write_log('Grid specification')
    call write_log('------------------')
    write(message,*) 'ewn             : ',model%general%ewn
    call write_log(trim(message))
    write(message,*) 'nsn             : ',model%general%nsn
    call write_log(trim(message))
    write(message,*) 'upn             : ',model%general%upn
    call write_log(trim(message))
    write(message,*) 'EW grid spacing : ',model%numerics%dew
    call write_log(trim(message))
    write(message,*) 'NS grid spacing : ',model%numerics%dns
    call write_log(trim(message))
    if (model%general%global_bc==GLOBAL_BC_PERIODIC) then
       write(message,*) 'Periodic global boundary conditions'
       call write_log(trim(message))
    elseif (model%general%global_bc==GLOBAL_BC_OUTFLOW) then
       write(message,*) 'Outflow global boundary conditions; scalars in global halo will be set to zero'
       call write_log(trim(message))
    elseif (model%general%global_bc==GLOBAL_BC_NO_PENETRATION) then
       write(message,*) 'No-penetration global boundary conditions; outflow set to zero at global boundaries'
       call write_log(trim(message))
    endif

    write(message,*) 'sigma file      : ',trim(model%funits%sigfile)
    call write_log(trim(message))
    call write_log('')

  end subroutine print_grid

!--------------------------------------------------------------------------------

  ! time
  subroutine handle_time(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

!TODO - Make the ice dynamic timestep more flexible.
!       To handle timesteps both greater and less than one year, we may want to
!        define ice_dt_option and ice_dt_count in place of the current dt.
!       For instance, ice_dt_option could be either 'nyears' or 'steps_per_year'.
!       For timesteps < 1 year, we would use ice_dt_option = 'steps_per_year'.
!       This would ensure that the ice sheet dynamic timestep divides evenly
!        into the mass balance timestep (= 1 year) when running with Glint.
    call GetValue(section,'tstart',model%numerics%tstart)
    call GetValue(section,'tend',model%numerics%tend)
    call GetValue(section,'dt',model%numerics%tinc)
    call GetValue(section,'subcyc',model%numerics%subcyc)
    call GetValue(section,'ntem',model%numerics%ntem)
    call GetValue(section,'profile',model%numerics%profile_period)

    call GetValue(section,'idiag',model%numerics%idiag)
    call GetValue(section,'jdiag',model%numerics%jdiag)

    !Note: Either dt_diag or ndiag can be specified in the config file.
    !      If dt_diag is specified, it is used to compute ndiag. (Output is written every ndiag timesteps.)
    call GetValue(section,'dt_diag',model%numerics%dt_diag)
    call GetValue(section,'ndiag',model%numerics%ndiag)

  end subroutine handle_time
  
!--------------------------------------------------------------------------------

  subroutine print_time(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log('Time steps')
    call write_log('----------')
    write(message,*) 'start time (yr)     : ',model%numerics%tstart
    call write_log(message)
    write(message,*) 'end time (yr)       : ',model%numerics%tend
    call write_log(message)
    write(message,*) 'time step (yr)      : ',model%numerics%tinc
    call write_log(message)
    write(message,*) 'thermal dt factor   : ',model%numerics%ntem
    call write_log(message)
    if ( (model%numerics%ntem < 1.0d0) .or. & 
       (floor(model%numerics%ntem) /= model%numerics%ntem) ) then
       call write_log('ntem is a multiplier on the basic time step.  It should be a positive integer.  Aborting.',GM_FATAL)
    endif
    write(message,*) 'profile frequency   : ',model%numerics%profile_period
    call write_log(message)

    if (model%numerics%dt_diag > 0.d0) then
       write(message,*) 'diagnostic interval (years):',model%numerics%dt_diag
       call write_log(message)
    elseif (model%numerics%ndiag > 0) then
       write(message,*) 'diagnostic interval (steps):',model%numerics%ndiag
       call write_log(message)
    endif

     !WHL - Written to log in glide_init_diag
!    write(message,*) 'idiag               : ',model%numerics%idiag
!    call write_log(message)
!    write(message,*) 'jdiag               : ',model%numerics%jdiag
!    call write_log(message)

    call write_log('')

  end subroutine print_time

!--------------------------------------------------------------------------------

  ! options
  subroutine handle_options(section, model)

    use glimmer_config
    use glide_types

    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'dycore',model%options%whichdycore)
    call GetValue(section,'evolution',model%options%whichevol)
    call GetValue(section,'temperature',model%options%whichtemp)
    call GetValue(section,'temp_init',model%options%temp_init)
    call GetValue(section,'flow_law',model%options%whichflwa)
    call GetValue(section,'slip_coeff',model%options%whichbtrc)
    call GetValue(section,'basal_water',model%options%whichbwat)
    call GetValue(section,'bmlt_float',model%options%whichbmlt_float)
    call GetValue(section,'basal_mass_balance',model%options%basal_mbal)
    call GetValue(section,'gthf',model%options%gthf)
    call GetValue(section,'isostasy',model%options%isostasy)
    call GetValue(section,'marine_margin',model%options%whichcalving)
    call GetValue(section,'calving_init',model%options%calving_init)
    call GetValue(section,'calving_domain',model%options%calving_domain)
    call GetValue(section,'vertical_integration',model%options%whichwvel)
    call GetValue(section,'topo_is_relaxed',model%options%whichrelaxed)
    call GetValue(section,'periodic_ew',model%options%periodic_ew)
    call GetValue(section,'sigma',model%options%which_sigma)
    call GetValue(section,'ioparams',model%funits%ncfile)

    ! Both terms 'hotstart' and 'restart' are supported in the config file, 
    ! but if they are both supplied for some reason, then restart will be used.
    ! 'restart' is the preferred term moving forward.  
    ! 'hotstart' is retained for backward compatability.
    call GetValue(section,'hotstart',model%options%is_restart)
    call GetValue(section,'restart',model%options%is_restart)

    call GetValue(section,'restart_extend_velo',model%options%restart_extend_velo)

    ! These are not currently supported
    !call GetValue(section, 'use_plume',model%options%use_plume)
    !call GetValue(section,'basal_proc',model%options%which_bproc)

  end subroutine handle_options

!--------------------------------------------------------------------------------

  !Higher order options
  subroutine handle_ho_options(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type) :: model
    
    call GetValue(section, 'which_ho_efvs',      model%options%which_ho_efvs)
    call GetValue(section, 'which_ho_disp',      model%options%which_ho_disp)
    call GetValue(section, 'which_ho_babc',      model%options%which_ho_babc)
    call GetValue(section, 'which_ho_resid',     model%options%which_ho_resid)
    call GetValue(section, 'which_ho_nonlinear', model%options%which_ho_nonlinear)
    call GetValue(section, 'which_ho_sparse',    model%options%which_ho_sparse)
    call GetValue(section, 'which_ho_approx',    model%options%which_ho_approx)
    call GetValue(section, 'which_ho_precond',   model%options%which_ho_precond)
    call GetValue(section, 'which_ho_gradient',  model%options%which_ho_gradient)
    call GetValue(section, 'which_ho_gradient_margin', model%options%which_ho_gradient_margin)
    call GetValue(section, 'which_ho_vertical_remap',  model%options%which_ho_vertical_remap)
    call GetValue(section, 'which_ho_assemble_beta',   model%options%which_ho_assemble_beta)
    call GetValue(section, 'which_ho_assemble_taud',   model%options%which_ho_assemble_taud)
    call GetValue(section, 'which_ho_assemble_bfric',  model%options%which_ho_assemble_bfric)
    call GetValue(section, 'which_ho_ground',    model%options%which_ho_ground)
    call GetValue(section, 'which_ho_flotation_function', model%options%which_ho_flotation_function)
    call GetValue(section, 'which_ho_ice_age',   model%options%which_ho_ice_age)
    call GetValue(section, 'glissade_maxiter',   model%options%glissade_maxiter)

  end subroutine handle_ho_options

!--------------------------------------------------------------------------------

  ! Handles external dycore options -- Doug Ranken 03/26/12
  subroutine handle_dycore_options(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type) :: model
    
    call GetValue(section, 'external_dycore_type', model%options%external_dycore_type)
    call GetValue(section, 'dycore_input_file',  model%options%dycore_input_file)
    if (model%options%external_dycore_type .eq. 1) model%options%whichdycore = 4 ! DYCORE_BISICLES
    if (model%options%external_dycore_type .eq. 2) model%options%whichdycore = 3 ! DYCORE_ALBANYFELIX

    print *,"In handle_dycore_options, external dycore type, input file = ", &
             model%options%external_dycore_type,model%options%dycore_input_file 
    ! print *,"In handle_dycore_options, whichdycore = ",model%options%whichdycore
  end subroutine handle_dycore_options

!--------------------------------------------------------------------------------

  subroutine print_options(model)

    use glide_types
    use glimmer_log

    use parallel

    implicit none

    type(glide_global_type)  :: model
    character(len=500) :: message

    ! basic Glide/Glimmer options

    character(len=*), dimension(0:4), parameter :: dycore = (/ &
         'glide              ', &  ! Glimmer SIA
         'glam               ', &  ! Payne-Price finite difference
         'glissade           ', &  ! prototype finite element
         'albany-felix       ', &  ! External Albany-FELIX finite element
         'bisicles           ' /)  ! External BISICLES-Chombo FVM 

    character(len=*), dimension(0:5), parameter :: evolution = (/ &
         'pseudo-diffusion                      ', &
         'ADI scheme                            ', &
         'iterated diffusion                    ', &
         'incremental remapping                 ', &   
         '1st order upwind                      ', &
         'thickness fixed at initial value      '  /)

    character(len=*), dimension(0:3), parameter :: temperature = (/ &
         'isothermal            ', &
         'prognostic temperature', &
         'constant in time      ', &
         'prognostic enthalpy   ' /)

    character(len=*), dimension(0:2), parameter :: temp_init = (/ &
         'set to 0 C             ', &
         'set to surface air temp', &
         'linear vertical profile' /)

    character(len=*), dimension(0:2), parameter :: flow_law = (/ &
         'const 1e-16 Pa^-n a^-1      ', &
         'Paterson and Budd (T = -5 C)', &
         'Paterson and Budd           ' /)

    !TODO - Rename slip_coeff to something like which_btrc?
    character(len=*), dimension(0:5), parameter :: slip_coeff = (/ &
         'no basal sliding       ', &
         'constant basal traction', &
         'constant where bwat > 0', &
         'constant where T = Tpmp', &
         'linear function of bmlt', &
         'tanh function of bwat  ' /)

    character(len=*), dimension(0:4), parameter :: basal_water = (/ &
         'none                     ', &
         'local water balance      ', &
         'local + steady-state flux', &
         'Constant value (= 10 m)  ', &
         'ocean penetration        ' /)
!!         'From basal proc model    '/) ! not supported

      ! basal proc model is disabled for now.
!!    character(len=*), dimension(0:2), parameter :: which_bproc = (/ &
!!         'Basal proc mod disabled '  , &
!!         'Basal proc, high res.   '   , &
!!         'Basal proc, fast calc.  ' /)
    character(len=*), dimension(0:0), parameter :: which_bproc = (/ &
         'Basal process model disabled ' /)

    character(len=*), dimension(0:1), parameter :: b_mbal = (/ &
         'not in continuity eqn', &
         'in continuity eqn    ' /)

    character(len=*), dimension(0:2), parameter :: which_bmlt_float = (/ &
         'none                     ', &
         'constant                 ', &
         'MISMIP+ melt rate profile' /)

    ! NOTE: Set gthf = 1 in the config file to read the geothermal heat flux from an input file.
    !       Otherwise it will be overwritten, even if the 'bheatflx' field is present.

    character(len=*), dimension(0:2), parameter :: gthf = (/ &
         'uniform geothermal flux         ', &
         'read flux from file, if present ', &
         'compute flux from diffusion eqn ' /)

    ! NOTE: This option has replaced the old do_isos option
    character(len=*), dimension(0:1), parameter :: isostasy = (/ &
         'no isostasy calculation         ', &
         'compute isostasy with model     ' /)

    !TODO - Change 'marine_margin' to 'calving'?  Would have to modify standard config files
    character(len=*), dimension(0:7), parameter :: marine_margin = (/ &
         'do nothing at marine margin     ', &
         'remove all floating ice         ', &
         'remove fraction of floating ice ', &
         'relaxed bedrock threshold       ', &
         'present bedrock threshold       ', &
         'Huybrechts grounding-line scheme', &
         'ice thickness threshold         ', &
         'damage-based calving scheme     ' /) 

    character(len=*), dimension(0:1), parameter :: init_calving = (/ &
         'no calving at initialization    ', &
         'ice calves at initialization    ' /)

    character(len=*), dimension(0:2), parameter :: domain_calving = (/ &
         'calving only at the ocean edge             ',  &
         'calving in all cells where criterion is met',  &
         'calving in cells connected to ocean edge   '/)

    character(len=*), dimension(0:1), parameter :: vertical_integration = (/ &
         'standard     ', &
         'obey upper BC' /)

    ! higher-order options

    character(len=*), dimension(0:2), parameter :: ho_whichefvs = (/ &
         'constant value                 ', & 
         'multiple of flow factor        ', &
         'nonlinear, from eff strain rate' /)

    character(len=*), dimension(-1:1), parameter :: ho_whichdisp = (/ &
         'no dissipation                    ', &
         '0-order SIA                       ', &
         'first-order model (Blatter-Pattyn)' /)

    character(len=*), dimension(0:13), parameter :: ho_whichbabc = (/ &
         'constant beta                                    ', &
         'beta depends on basal temp (melting or frozen)   ', &
         'till yield stress (Picard)                       ', &
         'beta is a function of bwat                       ', &
         'no slip (using large B^2)                        ', &
         'beta passed from CISM                            ', &
         'no slip (Dirichlet implementation)               ', &
         'till yield stress (Newton)                       ', &
         'beta as in ISMIP-HOM test C                      ', &
         'power law using effective pressure               ', &
         'Coulomb friction law w/ effec press              ', &
         'Coulomb friction law w/ effec press, const flwa_b', &
         'min of Coulomb stress and power-law stress (Tsai)', &
         'simple pattern of beta                           ' /)

    character(len=*), dimension(0:1), parameter :: which_ho_nonlinear = (/ &
         'use standard Picard iteration  ', &
         'use JFNK                       '/)

    character(len=*), dimension(0:4), parameter :: ho_whichresid = (/ &
         'max value                   ', &
         'max value ignoring ubas     ', &
         'mean value                  ', &
         'L2 norm of Ax-b = resid     ', &
         'relative L2 norm, |Ax-b|/|b|' /)

    character(len=*), dimension(-1:4), parameter :: ho_whichsparse = (/ &
         'PCG with incomplete Cholesky preconditioner', &
         'BiCG with LU preconditioner                ', &
         'GMRES with LU preconditioner               ', &
         'Native PCG solver, standard                ', &
         'Native PCG solver, Chronopoulos-Gear       ', &
         'Trilinos interface                         '/)

    character(len=*), dimension(-1:4), parameter :: ho_whichapprox = (/ &
         'SIA only (glissade_velo_sia)                     ', &
         'SIA only (glissade_velo_higher)                  ', &
         'SSA only (glissade_velo_higher)                  ', &
         'Blatter-Pattyn HO (glissade_velo_higher)         ', &
         'Depth-integrated L1L2 (glissade_velo_higher)     ', &
         'Depth-integrated viscosity (glissade_velo_higher)' /)

    character(len=*), dimension(0:2), parameter :: ho_whichprecond = (/ &
         'No preconditioner (glissade PCG)        ', &
         'Diagonal preconditioner (glissade PCG)  ', &
         'SIA preconditioner (glissade PCG)       ' /)

    character(len=*), dimension(0:1), parameter :: ho_whichgradient = (/ &
         'centered gradient (glissade dycore)      ', &
         'upstream gradient (glissade dycore)      ' /)

    character(len=*), dimension(0:2), parameter :: ho_whichgradient_margin = (/ &
         'all neighbor cells in gradient (glissade dycore)         ', &
         'ice-covered &/or land cells in gradient (glissade dycore)', &
         'only ice-covered cells in gradient (glissade dycore)     ' /)

    character(len=*), dimension(0:1), parameter :: ho_whichvertical_remap = (/ &
         'first-order accurate  ', &
         'second-order accurate ' /)

    character(len=*), dimension(0:1), parameter :: ho_whichassemble_beta = (/ &
         'standard finite-element assembly (glissade dycore) ', &
         'use local beta at each vertex (glissade dycore)    '  /)

    character(len=*), dimension(0:1), parameter :: ho_whichassemble_taud = (/ &
         'standard finite-element assembly (glissade dycore)       ', &
         'use local driving stress at each vertex (glissade dycore)'  /)

    character(len=*), dimension(0:1), parameter :: ho_whichassemble_bfric = (/ &
         'standard finite-element assembly (glissade dycore)       ', &
         'use local basal friction at each vertex (glissade dycore)'  /)

    character(len=*), dimension(0:2), parameter :: ho_whichground = (/ &
         'f_ground = 0 or 1; no GLP  (glissade dycore)       ', &
         '0 <= f_ground <= 1, based on GLP (glissade dycore) ', &
         'f_ground = 1 for all active cells (glissade dycore)' /)

    character(len=*), dimension(0:2), parameter :: ho_whichflotation_function = (/ &
         'f_pattyn = (-rhow*b)/(rhoi*H)  ', &
         '1/fpattyn = (rhoi*H)/(-rhow*b) ', &
         'linear = -rhow*b - rhoi*H      ' /)

    character(len=*), dimension(0:1), parameter :: ho_whichice_age = (/ &
         'ice age computation off', &
         'ice age computation on ' /)

    call write_log('GLIDE options')
    call write_log('-------------')

    write(message,*) 'I/O parameter file      : ',trim(model%funits%ncfile)
    call write_log(message)

    if (model%options%whichdycore < 0 .or. model%options%whichdycore >= size(dycore)) then
       call write_log('Error, dycore option out of range',GM_FATAL)
    end if
    write(message,*) 'Dycore                  : ',model%options%whichdycore,dycore(model%options%whichdycore)
    call write_log(message)

    ! unsupported dycore options
    if (model%options%whichdycore == DYCORE_GLAM) then
      call write_log('Glam dycore is not currently scientifically supported.  USE AT YOUR OWN RISK.', GM_WARNING)
    endif
    if (model%options%whichdycore == DYCORE_ALBANYFELIX) then
      call write_log('Albany-FELIX dycore is not currently scientifically supported.  USE AT YOUR OWN RISK.', GM_WARNING)
    endif
    if (model%options%whichdycore == DYCORE_BISICLES) then
      call write_log('BISICLES dycore is not currently scientifically supported.  USE AT YOUR OWN RISK.', GM_WARNING)
    endif

    ! Forbidden options associated with the Glide dycore
    if (model%options%whichdycore == DYCORE_GLIDE) then

       if (model%options%whichevol==EVOL_INC_REMAP     .or.  &
           model%options%whichevol==EVOL_UPWIND        .or.  &
           model%options%whichevol==EVOL_NO_THICKNESS) then
          call write_log('Error, Glam/glissade thickness evolution options cannot be used with Glide dycore', GM_FATAL)
       endif

       if (model%options%whichtemp == TEMP_ENTHALPY) then
          call write_log('Error, Enthalpy scheme cannot be used with Glide dycore', GM_FATAL)
       endif

       if (tasks > 1) then
          call write_log('Error, Glide dycore not supported for runs with more than one processor', GM_FATAL)
       end if

       if (model%options%whichevol==EVOL_ADI) then
          call write_log('Warning, exact restarts are not currently possible with ADI evolution', GM_WARNING)
       endif

    else   ! forbidden evolution options with dycores other than Glide

       if (model%options%whichevol==EVOL_PSEUDO_DIFF .or.  &
           model%options%whichevol==EVOL_ADI         .or.  &
           model%options%whichevol==EVOL_DIFFUSION) then
          call write_log('Error, Glide thickness evolution options cannot be used with glam/glissade dycore', GM_FATAL)
       endif

    endif

    ! Forbidden options for running in parallel
    if (tasks > 1 .and. (model%options%which_ho_sparse==HO_SPARSE_BICG  .or.  &
                         model%options%which_ho_sparse==HO_SPARSE_GMRES .or.  &
                         model%options%which_ho_sparse==HO_SPARSE_PCG_INCH) ) then
       call write_log('Error, SLAP solver not supported for more than one processor', GM_FATAL)
    end if

    if (tasks > 1 .and. model%options%which_ho_babc==HO_BABC_ISHOMC) then
       call write_log('Error, ISHOM C basal BCs not supported for more than one processor', GM_FATAL)
    endif

    if (tasks > 1 .and. model%options%whichbwat==BWATER_FLUX) then
       call write_log('Error, flux-based basal water option not supported for more than one processor', GM_FATAL)
    endif

    ! Forbidden options associated with Glam and Glissade dycores
   
    if (model%options%whichdycore == DYCORE_GLISSADE) then 

       if ( (model%options%which_ho_approx == HO_APPROX_SSA  .or.  &
             model%options%which_ho_approx == HO_APPROX_L1L2 .or.  &
             model%options%which_ho_approx == HO_APPROX_DIVA)   &
                                .and.                            &
             (model%options%which_ho_sparse == HO_SPARSE_PCG_STANDARD .or.    &
              model%options%which_ho_sparse == HO_SPARSE_PCG_CHRONGEAR) ) then
          if (model%options%which_ho_precond == HO_PRECOND_SIA) then
             call write_log('Error, cannot use SIA preconditioning for 2D solve', GM_FATAL)
          endif
       endif

       if (model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
          
          if (model%options%which_ho_disp == HO_DISP_FIRSTORDER ) then
             call write_log('Error, cannot use first-order dissipation with local SIA solver', GM_FATAL)
          endif

          if (model%general%global_bc == GLOBAL_BC_NO_PENETRATION) then
             call write_log('Error, cannot use no-penetration BC with local SIA solver', GM_FATAL)
          endif

       endif  ! Glissade local SIA solver

    endif

    if (model%options%whichdycore /= DYCORE_GLISSADE) then 

       if (model%options%which_ho_sparse == HO_SPARSE_PCG_STANDARD .or.   &
           model%options%which_ho_sparse == HO_SPARSE_PCG_CHRONGEAR) then
          call write_log('Error, native PCG solver requires glissade dycore', GM_FATAL)
       endif

       if (model%general%global_bc == GLOBAL_BC_NO_PENETRATION) then
          call write_log('Error, no-penetration BC requires glissade dycore', GM_FATAL)
       endif

    endif

    if (model%options%whichdycore == DYCORE_GLAM) then
       if (model%options%which_ho_approx /= HO_APPROX_BP) then
          call write_log('Error, Glam dycore must use Blatter-Pattyn approximation', GM_FATAL)
       endif
    endif

    ! Config specific to Albany-Felix dycore   
    if (model%options%whichdycore == DYCORE_ALBANYFELIX) then
       call write_log('Warning, Albany-FELIX dycore requires external libraries, and it is still in development!!!', GM_WARNING)
    endif

    !NOTE : Old option 3 (TEMP_REMAP_ADV) has been removed.
    ! If this has been set, then change to option 1 (TEMP_PROGNOSTIC), which applies to any dycore.

    if (model%options%whichtemp < 0 .or. model%options%whichtemp >= size(temperature)) then
       call write_log('Error, temperature option out of range',GM_FATAL)
    end if
    write(message,*) 'temperature calculation : ',model%options%whichtemp,temperature(model%options%whichtemp)
    call write_log(message)

    ! unsupported temperature options
    if (model%options%whichtemp == TEMP_ENTHALPY) then
      call write_log('Enthalpy-based formulation for solving temperature evolution is not currently scientifically supported.  &
           &USE AT YOUR OWN RISK.', GM_WARNING)
    endif

    if (model%options%temp_init < 0 .or. model%options%temp_init >= size(temp_init)) then
       call write_log('Error, temp_init option out of range',GM_FATAL)
    end if
    ! Note: If reading temperature from an input or restart file, the temp_init option is overridden,
    !        in which case it could be confusing here to write the option to the log file.
    !       The method actually used is written to the log file by glide_init_temp. 

    if (model%options%whichflwa < 0 .or. model%options%whichflwa >= size(flow_law)) then
       call write_log('Error, flow_law out of range',GM_FATAL)
    end if
    write(message,*) 'flow law                : ',model%options%whichflwa,flow_law(model%options%whichflwa)
    call write_log(message)

    if (model%options%whichbwat < 0 .or. model%options%whichbwat >= size(basal_water)) then
       call write_log('Error, basal_water out of range',GM_FATAL)
    end if
    write(message,*) 'basal_water             : ',model%options%whichbwat,basal_water(model%options%whichbwat)
    call write_log(message)

    ! unsupported basal_water options
    if (model%options%whichbwat == BWATER_FLUX) then
      call write_log('Steady state routing basal_water option is not currently scientifically supported.  &
           &USE AT YOUR OWN RISK.', GM_WARNING)
    endif
    if (model%options%whichbwat == BWATER_OCEAN_PENETRATION) then
      call write_log('Ocean penetration basal_water option is not currently scientifically supported.  &
           &USE AT YOUR OWN RISK.', GM_WARNING)
    endif

    if (model%options%whichcalving < 0 .or. model%options%whichcalving >= size(marine_margin)) then
       call write_log('Error, marine_margin out of range',GM_FATAL)
    end if
    write(message,*) 'marine_margin           : ', model%options%whichcalving, marine_margin(model%options%whichcalving)
    call write_log(message)

    if (model%options%calving_init < 0 .or. model%options%calving_init >= size(init_calving)) then
       call write_log('Error, calving_init out of range',GM_FATAL)
    end if
    write(message,*) 'calving_init            : ', model%options%calving_init, init_calving(model%options%calving_init)
    call write_log(message)

    if (model%options%calving_domain < 0 .or. model%options%calving_domain >= size(domain_calving)) then
       call write_log('Error, calving_domain out of range',GM_FATAL)
    end if
    write(message,*) 'calving_domain          : ', model%options%calving_domain, domain_calving(model%options%calving_domain)
    call write_log(message)

    ! unsupported calving options

    if (model%options%whichdycore == DYCORE_GLISSADE) then
       if (model%options%whichcalving == CALVING_FLOAT_FRACTION) then
          write(message,*) 'Warning, calving float fraction option deprecated with Glissade_dycore; set calving_timescale instead'
          call write_log(message)
       endif
    else   ! not Glissade
       if (model%options%whichcalving == CALVING_THCK_THRESHOLD) then
          call write_log('Error, calving thickness threshold model is supported for Glissade dycore only', GM_FATAL)
       endif
       if (model%options%whichcalving == CALVING_DAMAGE) then
          call write_log('Error, calving damage model is supported for Glissade dycore only', GM_FATAL)
       endif
       if (model%options%calving_domain /= CALVING_DOMAIN_OCEAN_EDGE) then
          write(message,*) 'Warning, calving domain can be selected for Glissade dycore only; user selection ignored'
          call write_log(message)
       endif
       if (model%calving%calving_timescale > 0.0d0) then
          write(message,*) 'Warning, calving timescale option suppored for Glissade dycore only; user selection ignored'
          call write_log(message)
       endif
    endif

    if (model%options%whichbtrc < 0 .or. model%options%whichbtrc >= size(slip_coeff)) then
       call write_log('Error, slip_coeff out of range',GM_FATAL)
    end if

    !WHL - Currently, not all basal traction options are supported for the Glissade SIA solver
    if (model%options%whichdycore == DYCORE_GLISSADE .and. model%options%which_ho_approx == HO_APPROX_LOCAL_SIA) then
       if (model%options%whichbtrc > BTRC_CONSTANT_TPMP) then
          call write_log('Error, slip_coeff out of range for Glissade dycore',GM_FATAL)
       end if
    endif

    write(message,*) 'slip_coeff              : ', model%options%whichbtrc, slip_coeff(model%options%whichbtrc)
    call write_log(message)

    if (model%options%whichevol < 0 .or. model%options%whichevol >= size(evolution)) then
       call write_log('Error, evolution out of range',GM_FATAL)
    end if

    write(message,*) 'evolution               : ', model%options%whichevol, evolution(model%options%whichevol)
    call write_log(message)

    if (model%options%whichwvel < 0 .or. model%options%whichwvel >= size(vertical_integration)) then
       call write_log('Error, vertical_integration out of range',GM_FATAL)
    end if

    if (model%options%whichwvel /= VERTINT_STANDARD .and. model%options%whichdycore /= DYCORE_GLIDE) then
       call write_log('Error, only standard vertical velocity calculation is supported for higher-order dycores.',GM_FATAL)
    end if

    write(message,*) 'vertical_integration    : ',model%options%whichwvel,vertical_integration(model%options%whichwvel)
    call write_log(message)

    if (model%options%whichbmlt_float < 0 .or. model%options%whichbmlt_float >= size(which_bmlt_float)) then
       call write_log('Error, bmlt_float out of range',GM_FATAL)
    end if

    write(message,*) 'basal melt, floating ice: ',model%options%whichbmlt_float, which_bmlt_float(model%options%whichbmlt_float)
    call write_log(message)

    if (model%options%basal_mbal < 0 .or. model%options%basal_mbal >= size(b_mbal)) then
       call write_log('Error, basal_mass_balance out of range',GM_FATAL)
    end if

    write(message,*) 'basal mass balance      : ',model%options%basal_mbal,b_mbal(model%options%basal_mbal)
    call write_log(message)

    if (model%options%gthf < 0 .or. model%options%gthf >= size(gthf)) then
       print*, 'gthf =', model%options%gthf
       call write_log('Error, geothermal flux option out of range',GM_FATAL)
    end if

    write(message,*) 'geothermal heat flux    : ',model%options%gthf,gthf(model%options%gthf)
    call write_log(message)

    if (model%options%isostasy < 0 .or. model%options%isostasy >= size(isostasy)) then
       print*, 'isostasy =', model%options%isostasy
       call write_log('Error, isostasy option out of range',GM_FATAL)
    end if

    write(message,*) 'isostasy                : ',model%options%isostasy,isostasy(model%options%isostasy)
    call write_log(message)

    if (model%options%whichrelaxed==1) then
       call write_log('First topo time slice has relaxed bedrock topography')
    end if

    if (model%options%periodic_ew) then
       if (model%options%whichevol == EVOL_ADI) then
          call write_log('Periodic boundary conditions not implemented in ADI scheme',GM_FATAL)
       end if
       call write_log('Periodic EW lateral boundary condition')
       call write_log('  Slightly cheated with how temperature is implemented.',GM_WARNING)
    end if

    if (model%options%is_restart == RESTART_TRUE) then
       call write_log('Restarting model from a previous run')
       if (model%options%restart_extend_velo == RESTART_EXTEND_VELO_TRUE) then
          call write_log('Using extended velocity fields for restart')
       endif
    end if

!!     This option is not currently supported
!!    if (model%options%which_bproc < 0 .or. model%options%which_bproc >= size(which_bproc)) then
!!       call write_log('Error, basal_proc out of range',GM_FATAL)
!!    end if
!!    write(message,*) 'basal_proc              : ',model%options%which_bproc,which_bproc(model%options%which_bproc)
!!    call write_log(message)

    !HO options

    if (model%options%whichdycore /= DYCORE_GLIDE) then   ! glam/glissade higher-order

       call write_log(' ')
       call write_log('Higher-order options:')
       call write_log('----------')

       write(message,*) 'ho_whichefvs            : ',model%options%which_ho_efvs,  &
                         ho_whichefvs(model%options%which_ho_efvs)
       call write_log(message)
       if (model%options%which_ho_efvs < 0 .or. model%options%which_ho_efvs >= size(ho_whichefvs)) then
          call write_log('Error, HO effective viscosity input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_whichdisp            : ',model%options%which_ho_disp,  &
                         ho_whichdisp(model%options%which_ho_disp)
       call write_log(message)
       if (model%options%which_ho_disp < -1 .or. model%options%which_ho_disp >= size(ho_whichdisp)-1) then
          call write_log('Error, HO dissipation input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_whichbabc            : ',model%options%which_ho_babc,  &
                         ho_whichbabc(model%options%which_ho_babc)
       call write_log(message)
       if (model%options%which_ho_babc < 0 .or. model%options%which_ho_babc >= size(ho_whichbabc)) then
          call write_log('Error, HO basal BC input out of range', GM_FATAL)
       end if

       ! unsupported ho-babc options
       !TODO - Decide if some of these are now supported?
       if (model%options%which_ho_babc == HO_BABC_YIELD_NEWTON) then
         call write_log('Yield stress higher-order basal boundary condition is not currently scientifically supported.  &
              &USE AT YOUR OWN RISK.', GM_WARNING)
       endif
       if (model%options%which_ho_babc == HO_BABC_POWERLAW) then
         call write_log('Weertman-style power law higher-order basal boundary condition is not currently scientifically &
              &supported.  USE AT YOUR OWN RISK.', GM_WARNING)
       endif
       if (model%options%which_ho_babc == HO_BABC_COULOMB_FRICTION          .or.  &
           model%options%which_ho_babc == HO_BABC_COULOMB_CONST_BASAL_FLWA) then
         call write_log('Coulomb friction law higher-order basal boundary condition is not currently scientifically supported.  &
              &USE AT YOUR OWN RISK.', GM_WARNING)
       endif

       write(message,*) 'which_ho_nonlinear      : ',model%options%which_ho_nonlinear,  &
                         which_ho_nonlinear(model%options%which_ho_nonlinear)
       call write_log(message)
       if (model%options%which_ho_nonlinear < 0 .or. model%options%which_ho_nonlinear >= size(which_ho_nonlinear)) then
          call write_log('Error, HO nonlinear solution input out of range', GM_FATAL)
       end if
       ! unsupported nonlinear options
       if (model%options%which_ho_nonlinear == HO_NONLIN_JFNK) then
         call write_log('JFNK treatment of nonlinearity in momentum balance is not currently scientifically supported.  &
              &USE AT YOUR OWN RISK.', GM_WARNING)
       endif

       write(message,*) 'ho_whichresid           : ',model%options%which_ho_resid,  &
                         ho_whichresid(model%options%which_ho_resid)
       call write_log(message)
       if (model%options%which_ho_resid < 0 .or. model%options%which_ho_resid >= size(ho_whichresid)) then
          call write_log('Error, HO residual input out of range', GM_FATAL)
       end if
       ! unsupported resid options
       if (model%options%which_ho_resid == HO_RESID_MAXU) then
         call write_log('Residual as max. value of normalized velocity vector update is not currently scientifically supported.  &
              &USE AT YOUR OWN RISK.', GM_WARNING)
       endif
       if (model%options%which_ho_resid == HO_RESID_MAXU_NO_UBAS) then
         call write_log('Residual as max. value of normalized velocity vector update with basal velocity omitted is not currently &
              &scientifically supported.  USE AT YOUR OWN RISK.', GM_WARNING)
       endif
       if (model%options%which_ho_resid == HO_RESID_MEANU) then
         call write_log('Residual as mean value of normalized velocity vector update is not currently scientifically supported.  &
              &USE AT YOUR OWN RISK.', GM_WARNING)
       endif

       write(message,*) 'ho_whichsparse          : ',model%options%which_ho_sparse,  &
                         ho_whichsparse(model%options%which_ho_sparse)
       call write_log(message)
       if (model%options%which_ho_sparse < -1 .or. model%options%which_ho_sparse >= size(ho_whichsparse)) then
          call write_log('Error, HO sparse solver input out of range', GM_FATAL)
       end if

       if (model%options%whichdycore == DYCORE_GLISSADE) then

          write(message,*) 'ho_whichapprox          : ',model%options%which_ho_approx,  &
                            ho_whichapprox(model%options%which_ho_approx)
          call write_log(message)
          if (model%options%which_ho_approx < -1 .or. model%options%which_ho_approx >= size(ho_whichapprox)-1) then
             call write_log('Error, Stokes approximation out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichgradient        : ',model%options%which_ho_gradient,  &
                            ho_whichgradient(model%options%which_ho_gradient)
          call write_log(message)
          if (model%options%which_ho_gradient < 0 .or. model%options%which_ho_gradient >= size(ho_whichgradient)) then
             call write_log('Error, gradient option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichgradient_margin : ',model%options%which_ho_gradient_margin,  &
                            ho_whichgradient_margin(model%options%which_ho_gradient_margin)
          call write_log(message)
          if (model%options%which_ho_gradient_margin < 0 .or. &
              model%options%which_ho_gradient_margin >= size(ho_whichgradient_margin)) then
             call write_log('Error, gradient margin option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichvertical_remap  : ',model%options%which_ho_vertical_remap,  &
                            ho_whichvertical_remap(model%options%which_ho_vertical_remap)
          call write_log(message)
          if (model%options%which_ho_vertical_remap < 0 .or. &
              model%options%which_ho_vertical_remap >= size(ho_whichvertical_remap)) then
             call write_log('Error, vertical remap option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichassemble_beta   : ',model%options%which_ho_assemble_beta,  &
                            ho_whichassemble_beta(model%options%which_ho_assemble_beta)
          call write_log(message)
          if (model%options%which_ho_assemble_beta < 0 .or. &
              model%options%which_ho_assemble_beta >= size(ho_whichassemble_beta)) then
             call write_log('Error, beta assembly option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichassemble_taud   : ',model%options%which_ho_assemble_taud,  &
                            ho_whichassemble_taud(model%options%which_ho_assemble_taud)
          call write_log(message)
          if (model%options%which_ho_assemble_taud < 0 .or. &
              model%options%which_ho_assemble_taud >= size(ho_whichassemble_taud)) then
             call write_log('Error, driving-stress assembly option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichassemble_bfric  : ',model%options%which_ho_assemble_bfric,  &
                            ho_whichassemble_bfric(model%options%which_ho_assemble_bfric)
          call write_log(message)
          if (model%options%which_ho_assemble_bfric < 0 .or. &
              model%options%which_ho_assemble_bfric >= size(ho_whichassemble_bfric)) then
             call write_log('Error, basal-friction assembly option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichground          : ',model%options%which_ho_ground,  &
                            ho_whichground(model%options%which_ho_ground)
          call write_log(message)
          if (model%options%which_ho_ground < 0 .or. model%options%which_ho_ground >= size(ho_whichground)) then
             call write_log('Error, ground option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichflotation_function:',model%options%which_ho_flotation_function,  &
                            ho_whichflotation_function(model%options%which_ho_flotation_function)
          call write_log(message)
          if (model%options%which_ho_flotation_function < 0 .or. &
               model%options%which_ho_flotation_function >= size(ho_whichflotation_function)) then
             call write_log('Error, flotation_function option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'ho_whichice_age         : ',model%options%which_ho_ice_age,  &
                            ho_whichice_age(model%options%which_ho_ice_age)
          call write_log(message)
          if (model%options%which_ho_ice_age < 0 .or. model%options%which_ho_ice_age >= size(ho_whichice_age)) then
             call write_log('Error, ice_age option out of range for glissade dycore', GM_FATAL)
          end if

          write(message,*) 'glissade_maxiter        : ',model%options%glissade_maxiter
          call write_log(message)

       end if

       if (model%options%whichdycore == DYCORE_GLISSADE .and.   &
           model%options%which_ho_ground == HO_GROUND_NO_GLP .and. &
           model%options%which_ho_flotation_function == HO_FLOTATION_FUNCTION_PATTYN) then
          write(message,*) 'Warning, Pattyn flotation function with no GLP tends to be unstable; inverse Pattyn is more stable'
          call write_log(message)
       endif

       if (model%options%whichdycore == DYCORE_GLISSADE .and.   &
           (model%options%which_ho_sparse == HO_SPARSE_PCG_STANDARD .or.  &
            model%options%which_ho_sparse == HO_SPARSE_PCG_CHRONGEAR) ) then 
          write(message,*) 'ho_whichprecond         : ',model%options%which_ho_precond,  &
                            ho_whichprecond(model%options%which_ho_precond)
          call write_log(message)
          if (model%options%which_ho_precond < 0 .or. model%options%which_ho_precond >= size(ho_whichprecond)) then
             call write_log('Error, glissade preconditioner out of range', GM_FATAL)
          end if
       end if

    endif   ! whichdycore

  end subroutine print_options

!--------------------------------------------------------------------------------

  ! parameters
  subroutine handle_parameters(section, model)

    use glimmer_config
    use glide_types
    use glimmer_log
    use glimmer_physcon, only: rhoi, rhoo, grav, shci, lhci, trpt

    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model
    real(dp), pointer, dimension(:) :: tempvar => NULL()
    integer :: loglevel

    !NOTE: The following physical constants have default values in glimmer_physcon.F90.
    !      Some test cases (e.g., MISMIP) specify different values. The default values
    !      can therefore be overridden by the user in the config file (except that certain
    !      constants in CESM's shr_const_mod cannot be overridden when CISM is coupled to CESM).
    !      These constants are not part of the model derived type.

#ifndef CCSMCOUPLED
    call GetValue(section,'rhoi', rhoi)
    call GetValue(section,'rhoo', rhoo)
    call GetValue(section,'grav', grav)
    call GetValue(section,'shci', shci)
    call GetValue(section,'lhci', lhci)
    call GetValue(section,'trpt', trpt)
#endif

    loglevel = GM_levels-GM_ERROR
    call GetValue(section,'log_level',loglevel)
    call glimmer_set_msg_level(loglevel)
    call GetValue(section,'ice_limit',        model%numerics%thklim)
    call GetValue(section,'ice_limit_temp',   model%numerics%thklim_temp)
    call GetValue(section,'marine_limit',     model%calving%marine_limit)
    call GetValue(section,'calving_fraction', model%calving%calving_fraction)
    call GetValue(section,'calving_timescale',model%calving%calving_timescale)
    call GetValue(section,'calving_minthck',  model%calving%calving_minthck)
    call GetValue(section,'damage_threshold', model%calving%damage_threshold)
    call GetValue(section,'geothermal',       model%paramets%geot)
    !TODO - Change default_flwa to flwa_constant?  Would have to change config files.
    !       Change flow_factor to flow_enhancement_factor?  Would have to change many SIA config files
    call GetValue(section,'flow_factor',      model%paramets%flow_enhancement_factor)
    call GetValue(section,'default_flwa',     model%paramets%default_flwa)
    call GetValue(section,'efvs_constant',    model%paramets%efvs_constant)
    call GetValue(section,'hydro_time',       model%paramets%hydtim)
    call GetValue(section,'max_slope',        model%paramets%max_slope)

    ! NOTE: bpar is used only for BTRC_TANH_BWAT
    !       btrac_max and btrac_slope are used (with btrac_const) for BTRC_LINEAR_BMLT
    !       btrac_const is used for several options

    call GetValue(section,'basal_tract_const', model%paramets%btrac_const)
    call GetValue(section,'basal_tract_max',   model%paramets%btrac_max)
    call GetValue(section,'basal_tract_slope', model%paramets%btrac_slope)

    !WHL - Changed this so that bpar can be read correctly from config file.
    !      This parameter is now called 'basal_tract_tanh' instead of 'basal_tract'.
    call GetValue(section,'basal_tract_tanh',  tempvar, 5)
    if (associated(tempvar)) then
!!       model%paramets%btrac_const = tempvar(1)  ! old code
       model%paramets%bpar(:) = tempvar(:)
       deallocate(tempvar)
    end if

!!    call GetValue(section,'sliding_constant',  model%climate%slidconst)  ! not currently used

    call GetValue(section,'ho_beta_const', model%velocity%ho_beta_const)
    call GetValue(section,'beta_grounded_min', model%velocity%beta_grounded_min)

    ! Friction law parameters
    call GetValue(section, 'friction_powerlaw_k', model%basal_physics%friction_powerlaw_k)
    call GetValue(section, 'coulomb_c', model%basal_physics%Coulomb_C)
    call GetValue(section, 'coulomb_bump_max_slope', model%basal_physics%Coulomb_Bump_max_slope)
    call GetValue(section, 'coulomb_bump_wavelength', model%basal_physics%Coulomb_bump_wavelength)
    call GetValue(section, 'flwa_basal', model%basal_physics%flwa_basal)
    call GetValue(section, 'powerlaw_c', model%basal_physics%powerlaw_C)
    call GetValue(section, 'powerlaw_m', model%basal_physics%powerlaw_m)

    ! ocean penetration parameterization parameter
    call GetValue(section,'p_ocean_penetration', model%paramets%p_ocean_penetration)

    ! ISMIP-HOM parameters
    call GetValue(section,'periodic_offset_ew',model%numerics%periodic_offset_ew)
    call GetValue(section,'periodic_offset_ns',model%numerics%periodic_offset_ns)

    ! MISMIP+ basal melting parameters
    call GetValue(section,'bmlt_float_rate', model%temper%bmlt_float_rate)
    call GetValue(section,'bmlt_float_omega', model%temper%bmlt_float_omega)
    call GetValue(section,'bmlt_float_h0', model%temper%bmlt_float_h0)
    call GetValue(section,'bmlt_float_z0', model%temper%bmlt_float_z0)

  end subroutine handle_parameters

!--------------------------------------------------------------------------------

  subroutine print_parameters(model)

    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log(' ')
    call write_log('Parameters')
    call write_log('----------')

    write(message,*) 'ice limit for dynamics (m)    : ', model%numerics%thklim
    call write_log(message)

    !Note: The Glissade dycore is known to crash for thklim = 0, but has not
    !      been extensively tested for small values of thklim.
    !      Values smaller than 1 mm may be OK, but no guarantees.
    if (model%options%whichdycore == DYCORE_GLISSADE .and.   &
        model%numerics%thklim < 1.d-3) then   ! 1 mm
       call write_log('ice limit (thklim) is too small for Glissade dycore', GM_FATAL)
    endif

    if (model%options%whichdycore /= DYCORE_GLIDE) then
       write(message,*) 'ice limit for temperature (m) : ', model%numerics%thklim_temp
       call write_log(message)
    endif

    if (model%options%whichcalving == CALVING_FLOAT_FRACTION) then
       write(message,*) 'ice fraction lost due to calving : ', model%calving%calving_fraction
       call write_log(message)
    end if

    if (model%options%whichcalving == CALVING_RELX_THRESHOLD .or.  &
        model%options%whichcalving == CALVING_TOPG_THRESHOLD) then
       write(message,*) 'marine depth limit (m)        : ', model%calving%marine_limit
       call write_log(message)
    endif

    if (model%options%whichcalving == CALVING_THCK_THRESHOLD) then
       write(message,*) 'calving thickness limit (m)   : ', model%calving%calving_minthck
       call write_log(message)
    endif

    if (model%options%whichcalving == CALVING_DAMAGE) then
       write(message,*) 'calving damage threshold      : ', model%calving%damage_threshold
       call write_log(message)
    end if

    if (model%calving%calving_timescale > 0.0d0) then
       write(message,*) 'calving time scale (yr)       : ', model%calving%calving_timescale
       call write_log(message)
    endif

    write(message,*) 'ice density (kg/m^3)          : ', rhoi
    call write_log(message)

    write(message,*) 'ocean density (kg/m^3)        : ', rhoo
    call write_log(message)

    write(message,*) 'gravitational accel (m/s^2)   : ', grav
    call write_log(message)

    write(message,*) 'heat capacity of ice (J/kg/K) : ', shci
    call write_log(message)

    write(message,*) 'latent heat of ice (J/kg)     : ', lhci
    call write_log(message)

    write(message,*) 'triple point of water (K)     : ', trpt
    call write_log(message)

    write(message,*) 'geothermal flux  (W/m^2)      : ', model%paramets%geot
    call write_log(message)

    write(message,*) 'flow enhancement factor       : ', model%paramets%flow_enhancement_factor
    call write_log(message)

    write(message,*) 'basal hydro time constant (yr): ', model%paramets%hydtim
    call write_log(message)

    if (model%options%whichdycore == DYCORE_GLISSADE) then
       write(message,*) 'max surface slope             : ', model%paramets%max_slope
       call write_log(message)
    end if       
 
    if (model%options%whichflwa == FLWA_CONST_FLWA) then
       write(message,*) 'constant flow factor (Pa^-n yr^-1) :', model%paramets%default_flwa
       call write_log(message)
    end if

    if (model%options%which_ho_efvs == HO_EFVS_CONSTANT) then
       write(message,*) 'constant effec viscosity (Pa yr)   :', model%paramets%efvs_constant
       call write_log(message)
    end if

    if (model%options%whichbtrc == BTRC_CONSTANT      .or.  &
        model%options%whichbtrc == BTRC_CONSTANT_BWAT .or.  &
        model%options%whichbtrc == BTRC_LINEAR_BMLT   .or.  &
        model%options%whichbtrc == BTRC_CONSTANT_TPMP) then
       write(message,*) 'basal traction param (m/yr/Pa)      : ', model%paramets%btrac_const
       call write_log(message)
    end if

    if (model%options%whichbtrc == BTRC_TANH_BWAT) then
       write(message,*) 'basal traction tanh factors: ',model%paramets%bpar(1)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(2)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(3)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(4)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(5)
       call write_log(message)
    end if

    if (model%options%whichbtrc == BTRC_LINEAR_BMLT) then
       write(message,*) 'basal traction max            : ',model%paramets%btrac_max
       call write_log(message)
       write(message,*) 'basal traction slope          : ',model%paramets%btrac_slope
       call write_log(message)
    end if

    if (model%options%which_ho_babc == HO_BABC_CONSTANT) then
       write(message,*) 'uniform beta (Pa yr/m)        : ',model%velocity%ho_beta_const
       call write_log(message)
    end if

    if (model%options%which_ho_babc == HO_BABC_ISHOMC) then
       if (model%general%ewn /= model%general%nsn) then
          call write_log('Error, must have ewn = nsn for ISMIP-HOM test C', GM_FATAL)
       endif
    endif

    if (model%options%which_ho_babc == HO_BABC_POWERLAW) then
       write(message,*) 'roughness parameter, k, for power-law friction law : ',model%basal_physics%friction_powerlaw_k
       call write_log(message)
    end if

    if (model%options%which_ho_babc == HO_BABC_COULOMB_FRICTION          .or.  &
        model%options%which_ho_babc == HO_BABC_COULOMB_CONST_BASAL_FLWA) then
       write(message,*) 'C coefficient for Coulomb friction law       : ', model%basal_physics%Coulomb_C
       call write_log(message)
       write(message,*) 'bed bump max. slope for Coulomb friction law : ', model%basal_physics%Coulomb_Bump_max_slope
       call write_log(message)
       write(message,*) 'bed bump wavelength for Coulomb friction law : ', model%basal_physics%Coulomb_bump_wavelength
       call write_log(message)
       if (model%options%which_ho_babc == HO_BABC_COULOMB_CONST_BASAL_FLWA) then
          write(message,*) 'constant basal flwa for Coulomb friction law : ', model%basal_physics%flwa_basal
          call write_log(message)
       endif
    end if

    if (model%options%which_ho_babc == HO_BABC_COULOMB_POWERLAW_TSAI) then
       write(message,*) 'C coefficient for Coulomb friction law       : ', model%basal_physics%Coulomb_C
       call write_log(message)
       write(message,*) 'C coefficient for power law, Pa (m/yr)^(-1/3): ', model%basal_physics%powerlaw_C
       call write_log(message)
       write(message,*) 'm exponent for power law                     : ', model%basal_physics%powerlaw_m
       call write_log(message)
    endif

    if (model%options%whichbwat == BWATER_OCEAN_PENETRATION) then
       write(message,*) 'p_ocean_penetration                : ', model%paramets%p_ocean_penetration
       call write_log(message)
    endif

    if (model%velocity%beta_grounded_min > 0.d0) then
       write(message,*) 'min beta for grounded ice (Pa yr/m): ', model%velocity%beta_grounded_min
       call write_log(message)
    endif
    
    if (model%numerics%idiag < 1 .or. model%numerics%idiag > model%general%ewn     &
                                        .or.                                                     &
        model%numerics%jdiag < 1 .or. model%numerics%jdiag > model%general%nsn) then
        call write_log('Error, global diagnostic point (idiag, jdiag) is out of bounds', GM_FATAL)
    endif

    ! ISMIP-HOM parameters
    if (model%numerics%periodic_offset_ew /= 0.d0) then
       write(message,*) 'periodic offset_ew (m)  :  ',model%numerics%periodic_offset_ew
       call write_log(message)
    endif

    if (model%numerics%periodic_offset_ns /= 0.d0) then
       write(message,*) 'periodic offset_ns (m)  :  ',model%numerics%periodic_offset_ns
       call write_log(message)
    endif

    ! parameters for basal melting of floating ice (including MISMIP+)
    if (model%options%whichbmlt_float == BMLT_FLOAT_CONSTANT) then
       write(message,*) 'bmlt_float_rate (m/yr)   :  ', model%temper%bmlt_float_rate * scyr
       call write_log(message)
    elseif (model%options%whichbmlt_float == BMLT_FLOAT_MISMIP) then
       write(message,*) 'bmlt_float_omega (yr^-1) :  ', model%temper%bmlt_float_omega * scyr
       call write_log(message)
       write(message,*) 'bmlt_float_h0 (m)        :  ', model%temper%bmlt_float_h0
       call write_log(message)
       write(message,*) 'bmlt_float_z0 (m)        :  ', model%temper%bmlt_float_z0
       call write_log(message)
    endif

    call write_log('')

  end subroutine print_parameters

!--------------------------------------------------------------------------------

  ! Sigma levels
  subroutine handle_sigma(section, model)

    use glimmer_config
    use glide_types
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    if (model%options%which_sigma==SIGMA_EXTERNAL) then
       call write_log('Sigma levels specified twice - use only'// &
            ' config file or separate file, not both',GM_FATAL)
    else
       call GetValue(section,'sigma_levels',model%numerics%sigma,model%general%upn)
    end if

  end subroutine handle_sigma

!--------------------------------------------------------------------------------

  subroutine print_sigma(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message,temp
    integer :: i

    call write_log('Sigma levels:')
    call write_log('------------------')
    message=''
    do i=1,model%general%upn
       write(temp,'(f6.3)') model%numerics%sigma(i)
       message=trim(message)//trim(temp)
    enddo
    call write_log(trim(message))
    call write_log('')
    
  end subroutine print_sigma

!--------------------------------------------------------------------------------

  ! geothermal heat flux calculations
  subroutine handle_gthf(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'num_dim',model%lithot%num_dim)
    call GetValue(section,'nlayer',model%lithot%nlayer)
    call GetValue(section,'surft',model%lithot%surft)
    call GetValue(section,'rock_base',model%lithot%rock_base)
    call GetValue(section,'numt',model%lithot%numt)
    call GetValue(section,'rho',model%lithot%rho_r)
    call GetValue(section,'shc',model%lithot%shc_r)
    call GetValue(section,'con',model%lithot%con_r)
  end subroutine handle_gthf

!--------------------------------------------------------------------------------

  subroutine print_gthf(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message
    
    if (model%options%gthf == GTHF_COMPUTE) then
       call write_log('Geothermal heat flux configuration')
       call write_log('----------------------------------')
       if (model%lithot%num_dim==1) then
          call write_log('solve 1D diffusion equation')
       else if (model%lithot%num_dim==3) then          
          call write_log('solve 3D diffusion equation')
       else
          call write_log('Wrong number of dimensions.',GM_FATAL,__FILE__,__LINE__)
       end if
       write(message,*) 'number of layers                     : ',model%lithot%nlayer
       call write_log(message)
       write(message,*) 'initial surface temperature          : ',model%lithot%surft
       call write_log(message)
       write(message,*) 'rock base                            : ',model%lithot%rock_base
       call write_log(message)
       write(message,*) 'density of rock layer                : ',model%lithot%rho_r
       call write_log(message)
       write(message,*) 'specific heat capacity of rock layer : ',model%lithot%shc_r
       call write_log(message)
       write(message,*) 'thermal conductivity of rock layer   : ',model%lithot%con_r
       call write_log(message)
       write(message,*) 'number of time steps for spin-up     : ',model%lithot%numt
       call write_log(message)
       call write_log('')
    end if
  end subroutine print_gthf

!--------------------------------------------------------------------------------

  subroutine handle_isostasy(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'lithosphere',model%isostasy%lithosphere)
    call GetValue(section,'asthenosphere',model%isostasy%asthenosphere)
    call GetValue(section,'relaxed_tau',model%isostasy%relaxed_tau)
    call GetValue(section,'update',model%isostasy%period)

    !NOTE: This value used to be in a separate section ('elastic lithosphere')
    !      Now part of 'isostasy' section
    call GetValue(section,'flexural_rigidity',model%isostasy%rbel%d)

!!    call GetSection(config,section,'elastic lithosphere')
!!    if (associated(section)) then
!!       call GetValue(section,'flexural_rigidity',isos%rbel%d)
!!    end if

  end subroutine handle_isostasy

!--------------------------------------------------------------------------------

  subroutine print_isostasy(model)
    use glide_types
    use glimmer_log
    use parallel, only: tasks
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message
    
    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call write_log('Isostasy')
       call write_log('--------')

       if (model%isostasy%lithosphere==LITHOSPHERE_LOCAL) then
          call write_log('using local lithosphere approximation')
       else if (model%isostasy%lithosphere==LITHOSPHERE_ELASTIC) then
          if (tasks > 1) then
             call write_log('Error, elastic lithosphere not supported for multiple processors',GM_FATAL)
          endif
          call write_log('using elastic lithosphere approximation')
          write(message,*) ' flexural rigidity : ', model%isostasy%rbel%d
          call write_log(message)
          write(message,*) ' update period (yr): ', model%isostasy%period
          call write_log(message)
       else
          call write_log('Error, unknown lithosphere option',GM_FATAL)
       end if

       if (model%isostasy%asthenosphere==ASTHENOSPHERE_FLUID) then
          call write_log('using fluid mantle')
       else if (model%isostasy%asthenosphere==ASTHENOSPHERE_RELAXING) then
          call write_log('using relaxing mantle')
          write(message,*) ' characteristic time constant (yr): ', model%isostasy%relaxed_tau
          call write_log(message)
       else
          call write_log('Error, unknown asthenosphere option',GM_FATAL)
       end if
       call write_log('')
    endif   ! compute isostasy

  end subroutine print_isostasy

!--------------------------------------------------------------------------------

! These options are disabled for now.

!!  subroutine handle_till_options(section,model)
!!    !Till options
!!    use glimmer_config
!!    use glide_types
!!    implicit none
!!    type(ConfigSection), pointer :: section
!!    type(glide_global_type) :: model

!!    if (model%options%which_bproc==1) then
!!        call GetValue(section, 'fric',  model%basalproc%fric)
!!        call GetValue(section, 'etillo',  model%basalproc%etillo)
!!        call GetValue(section, 'No',  model%basalproc%No)
!!        call GetValue(section, 'Comp',  model%basalproc%Comp)
!!        call GetValue(section, 'Cv',  model%basalproc%Cv)
!!        call GetValue(section, 'Kh',  model%basalproc%Kh)
!!    else if (model%options%which_bproc==2) then
!!        call GetValue(section, 'aconst',  model%basalproc%aconst)
!!        call GetValue(section, 'bconst',  model%basalproc%bconst)
!!    end if
!!    if (model%options%which_bproc > 0) then
!!        call GetValue(section, 'Zs',  model%basalproc%Zs)
!!        call GetValue(section, 'tnodes',  model%basalproc%tnodes)
!!        call GetValue(section, 'till_hot', model%basalproc%till_hot)
!!    end if  
!!  end subroutine handle_till_options    

!!  subroutine print_till_options(model)
!!    use glide_types
!!    use glimmer_log
!!    implicit none
!!    type(glide_global_type)  :: model
!!    character(len=100) :: message

!!    if (model%options%which_bproc > 0) then 
!!        call write_log('Till options')
!!        call write_log('----------')
!!        if (model%options%which_bproc==1) then
!!            write(message,*) 'Internal friction           : ',model%basalproc%fric
!!            call write_log(message)
!!            write(message,*) 'Reference void ratio        : ',model%basalproc%etillo
!!            call write_log(message)
!!            write(message,*) 'Reference effective Stress  : ',model%basalproc%No
!!            call write_log(message)
!!            write(message,*) 'Compressibility             : ',model%basalproc%Comp
!!            call write_log(message)
!!            write(message,*) 'Diffusivity                 : ',model%basalproc%Cv
!!            call write_log(message)
!!            write(message,*) 'Hyd. conductivity           : ',model%basalproc%Kh
!!            call write_log(message)
!!        end if
!!        if (model%options%which_bproc==2) then
!!            write(message,*) 'aconst  : ',model%basalproc%aconst
!!            call write_log(message)
!!            write(message,*) 'bconst  : ',model%basalproc%aconst
!!            call write_log(message)
!!        end if
!!        write(message,*) 'Solid till thickness : ',model%basalproc%Zs
!!        call write_log(message)
!!        write(message,*) 'Till nodes number : ',model%basalproc%tnodes
!!        call write_log(message)
!!        write(message,*) 'till_hot  :',model%basalproc%till_hot
!!        call write_log(message)
!!    end if
!!  end subroutine print_till_options

!--------------------------------------------------------------------------------

  subroutine define_glide_restart_variables(options)
    !> This subroutine analyzes the glide/glissade options input by the user in the config file
    !> and determines which variables are necessary for an exact restart.  MJH 1/11/2013

    ! Please comment thoroughly the reasons why a particular variable needs to be a restart variable for a given config.
    ! Note: this subroutine assumes that any restart variables you add you loadable.  Check glide_vars.def to make sure any variables you add have load: 1

    use glide_types
    use glide_io, only: glide_add_to_restart_variable_list

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glide_options), intent (in) :: options  !> Derived type holding all model options

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    !This was the restart list as of 1/11/13 using the old hot=1 system in glide_vars.def:
    !restart_variable_list=' lat  relx  tauf  thk  thkmask  topg  bheatflx  bmlt_ground  bwat  uvel  vvel  wgrd  flwa  temp  litho_temp  age '

    ! Start with a few variables that we always want - prognostic variables and b.c.
    ! topg - needed to reconstruct all other geometry fields
    ! thk - prognostic variable
    ! temp - prognostic variable
    ! Note: the conversion from temp/flwa to tempstag/flwastag (if necessary) happens in glide_io.F90
    ! bheatflx, artm, acab - boundary conditions.  Of course if these fields are 0 they don't need 
    !        to be in the restart file, but without adding a check for that we cannot assume any of them are.
    !        There are some options where artm would not be needed.  Logic could be added to make that distinction.
    !        Note that bheatflx may not be an input variable but can also be assigned as a parameter in the config file!
    call glide_add_to_restart_variable_list('topg thk temp bheatflx artm acab')

    ! add dycore specific restart variables
    select case (options%whichdycore)

      case (DYCORE_GLIDE)
        ! thkmask - TODO is this needed?
        ! wgrd & wvel - temp driver calculates weff = f(wgrd, wvel) so both are needed by temp code.
        !               It looks possible to calculate wvel on a restart from wgrd because wvel does not 
        !               appear to require a time derivative (see subroutine wvelintg).  
        !               wgrd does require time derivatives and therefore should be
        !               calculated at the end of each time step and stored as a restart variable
        !               so that the time derivatives do not need to be restart variables.
        !               For now I am calculating wvel at the same time (end of glide time step) 
        !               and then saving both as restart variables.  This has the advantage of
        !               them being on consistent time levels in the output file.  
        !               (If we waited to calculate wvel in the temp driver, we would not need to
        !                add it as a restart variable, been then in the output wgrd and wvel would
        !                be based on different time levels.)
        ! flwa - in principal this could be reconstructed from temp.  However in the current 
        !        implementation of glide the flwa calculation occurs after temp evolution but 
        !        before thk evolution.  This means flwa is calculated from the current temp and 
        !        the old thk.  The old thk is not available on a restart (just the current thk).
        !        (thk is needed to calculate flwa for 1) a mask for where ice is, 2) correction for pmp.)
        call glide_add_to_restart_variable_list('thkmask wgrd wvel flwa uvel vvel')
    
        ! slip option for SIA
        select case (options%whichbtrc)
          case (0)
            ! no restart variable needed when no-slip is chosen
          case default
            ! when a slip option is chosen, ubas & vbas are needed by the temperature solver
            ! for calculating basal heating prior to the first calculation of velocity.
            ! Rather than recalculate the sliding field on restart, it is easier and 
            ! less error-prone to have them be restart variables.  
            ! This could either be done by making ubas, vbas restart variables or
            ! having them assigned from the bottom level of uvel,vvel on init
            ! Note that btrc and soft are not needed as restart variables because
            ! their current implementation is as a scalar ('basal_tract_const' config parameter).
            ! If they are ever implemented as 2-d fields, then they (probably just one of them)
            ! should become restart variables.
            
            ! Nothing needs to happen because ubas,vbas are assigned from uvel,vel in glide_init_state_diagnostic()
        end select

      case (DYCORE_GLAM, DYCORE_GLISSADE)
        ! beta - b.c. needed for runs with sliding - could add logic to only include in that case
        ! flwa is not needed for glissade.
        ! TODO not sure if thkmask is needed for HO

        call glide_add_to_restart_variable_list('thkmask kinbcmask bfricflx dissip')

        ! uvel,vvel: These are needed for an exact restart because we can only recalculate
        !            them to within the picard/jfnk convergence tolerance.
        ! uvel/vvel_extend - These are identical to uvel and vvel, except that the mesh includes
        !                     points along the north and east boundaries of the domain.
        !                    CISM requires these fields for exact restart if the boundary velocities are nonzero,
        !                     as in MISMIP test problems with periodic BCs.
        !                    To output these fields, the user must set restart_extend_velo = 1 in the config file.
        ! Note: It never hurts to write uvel/vvel_extend in place of uvel/vvel. But for most cases where restart
        !       is required (e.g., whole-ice-sheet simulations), velocities are zero along the boundaries
        !       and uvel/vvel are sufficient.

        if (options%restart_extend_velo == RESTART_EXTEND_VELO_TRUE) then
           call glide_add_to_restart_variable_list('uvel_extend vvel_extend')
        else
           call glide_add_to_restart_variable_list('uvel vvel')
        endif

        ! Glissade approximation options
        select case (options%which_ho_approx)

           case (HO_APPROX_DIVA)
              ! DIVA requires the 2D velocity, basal traction and effective viscosity for exact restart.
              ! Since the 2D velocity and basal traction are located on the staggered grid, these fields
              !  must be written to and read from the extended grid if velocites are nonzero at the boundaries.
              !
              ! Note: The 2D velocity is needed if the DIVA scheme solves for the mean velocity.
              !       If DIVA is configured to solve for the velocity at a specific level (e.g., the surface),
              !       then the 2D velocity could instead be copied from the 3D velocity array.
              ! Note: In addition to uvel/vvel_2D, DIVA requires the full 3D velocity field for exact restart,
              !       because horizontal transport is done before updating the velocity.
              
              if (options%restart_extend_velo == RESTART_EXTEND_VELO_TRUE) then
                 call glide_add_to_restart_variable_list('uvel_2d_extend vvel_2d_extend btractx_extend btracty_extend efvs')
              else
                 call glide_add_to_restart_variable_list('uvel_2d vvel_2d btractx btracty efvs')
              endif
              
           case default
              ! Other approximations (including SSA and L1L2) use the 3D uvel and vvel to initialize the velocity

        end select   ! which_ho_approx
           
      end select ! which_dycore

    ! ==== Other non-dycore specific options ====

    ! basal water option
    select case (options%whichbwat)
      case (BWATER_NONE, BWATER_CONST)
        ! no restart variables needed
      case default
        ! restart needs to know bwat value
        call glide_add_to_restart_variable_list('bwat')
    end select

    ! internal water option (for enthalpy scheme)
    select case (options%whichtemp)
      case (TEMP_ENTHALPY)
        ! restart needs to know internal water fraction
        call glide_add_to_restart_variable_list('waterfrac')
      case default
        ! no restart variables needed
    end select

    ! basal melting option
    select case (options%whichbmlt_float)
       case (BMLT_FLOAT_CONSTANT)
          ! bmlt_float_mask needs to be in restart file if not = 0 everywhere
          ! TODO - Add bmlt_float_mask to the restart file only if not = 0 everywhere?
          call glide_add_to_restart_variable_list('bmlt_float_mask')
       case default
          ! no restart variables needed
    end select

    ! basal sliding option
    select case (options%which_ho_babc)
      case (HO_BABC_POWERLAW, HO_BABC_COULOMB_FRICTION, HO_BABC_COULOMB_CONST_BASAL_FLWA)
        ! These friction laws need effective pressure
         !TODO - Does effecpress need to be a restart variable?
        call glide_add_to_restart_variable_list('effecpress')
        ! C_space_factor needs to be in restart file if not = 1 everywhere
        !TODO - Add C_space_factor to the restart file only if not = 1 everywhere?
        call glide_add_to_restart_variable_list('C_space_factor')
      case(HO_BABC_COULOMB_POWERLAW_TSAI)
        call glide_add_to_restart_variable_list('effecpress')
      case default
        ! Other HO basal boundary conditions may need the external beta field  (although there are a few that don't)
        !Note: If using beta from an external file, then 'beta' here needs to be the fixed, external field,
        !      and not the internal beta field that may have been weighted by the grounded fraction or otherwise adjusted.
        call glide_add_to_restart_variable_list('beta')
    end select

    ! geothermal heat flux option
    select case (options%gthf)
      case(GTHF_COMPUTE)
         ! restart needs to know lithosphere temperature
         call glide_add_to_restart_variable_list('litho_temp')
      case default
         ! no restart variables needed
    end select

    !WHL - added isostasy option
    select case (options%isostasy)
      case(ISOSTASY_COMPUTE)
         ! restart needs to know relaxation depth
         ! TODO MJH: I suspect that relx is only needed when asthenosphere=1 (relaxing mantle), but I'm not sure -
         !      this should be tested when isostasy implementation is finalized/tested.
         call glide_add_to_restart_variable_list('relx')
      case default
         ! no new restart variables needed
    end select

    !WHL - added ice_age option
    !      Note: Ice age is a diagnostic field, not part of the prognostic ice state.
    !      Omitting it from restart will only break the diagnostic.
    select case (options%which_ho_ice_age)
       case(HO_ICE_AGE_COMPUTE)
          call glide_add_to_restart_variable_list('ice_age')
       case default
          ! no restart variables needed
    end select
    !
    ! basal processes module - requires tauf for a restart
!!    if (options%which_bproc /= BAS_PROC_DISABLED ) then
!!        call glide_add_to_restart_variable_list('tauf')
!!    endif

    ! TODO bmlt was set as a restart variable, but I'm not sure when or if it is needed.

    ! TODO age should be a restart variable if it is an input variable.  
    ! Same goes for b.c. (bheatflxm, artm, acab) and any other tracers that get introduced.
    ! These could be included all the time (as I have down above for b.c.), or 
    ! we could add logic to only include them when they were in the input file.
    ! To do this, this subroutine would have to be moved to after where input files are read,
    ! glide_io_readall(), but before the output files are created, glide_io_createall()

    ! TODO lat is only needed for some climate drivers.  It is not needed for cism_driver.
    ! Need to add logic that will add it only when those drivers are used.

  end subroutine define_glide_restart_variables

!--------------------------------------------------------------------------------

end module glide_setup

!--------------------------------------------------------------------------------
