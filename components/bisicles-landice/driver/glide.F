! WJS (1-30-12): The following (turning optimization off) is needed as a workaround for an
! xlf compiler bug, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

!CLEANUP - glide.F90
! Moved higher-order computations to a new module, glissade.F90.
! Simplified glide.F90 to include only SIA computations.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide.F90 - part of the Community Ice Sheet Model (CISM)  
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

!=======================================================================

module glide

  ! Driver for Glide (serial, SIA) dynamical core
  
  use glide_types
  use glide_stop
  use glide_io
  use glide_lithot
  use glide_profile
  use glimmer_config
  use glimmer_global, only: dp

  use glimmer_paramets, only: oldglide

  implicit none

  integer, private, parameter :: dummyunit=99

!WHL - debug
  logical, parameter :: verbose_glide = .false.

contains

!=======================================================================

  subroutine glide_config(model,config,fileunit)

    ! Read glide configuration from file and print it to the log

    use glide_setup
    use isostasy
    use glimmer_ncparams
    use glimmer_config
    use glimmer_map_init
    use glimmer_filenames

    implicit none

    type(glide_global_type), intent(inout) :: model  ! model instance
    type(ConfigSection), pointer  :: config          ! structure holding sections of configuration file
    integer, intent(in), optional :: fileunit        ! fileunit for reading config file 

    type(ConfigSection), pointer :: ncconfig
    integer :: unit

    unit = 99
    if (present(fileunit)) then
       unit = fileunit
    endif

    ! read configuration file
    call glide_readconfig (model,config)
    call glide_printconfig(model)

    ! read sigma levels from config file, if present
    call glide_read_sigma(model,config)

    !WHL - Moved isostasy configuration to glide_setup
!    call isos_readconfig(model%isos,config)
!    call isos_printconfig(model%isos)

    ! read mapping from config file
    ! **** Use of dew and dns here is an ugly fudge that
    ! **** allows the use of old [GLINT projection] config section
    ! **** for backwards compatibility. It will be deleted soon.
    ! **** (You have been warned!)
    ! **** N.B. Here, dew and dns are unscaled - i.e. real distances in m

    call glimmap_readconfig(model%projection,   config,   &
                            model%numerics%dew, model%numerics%dns)

    ! netCDF I/O
    if (trim(model%funits%ncfile) == '') then
       ncconfig => config
    else
       call ConfigRead(process_path(model%funits%ncfile), ncconfig, unit)
    end if

    call glimmer_nc_readparams(model, ncconfig)

  end subroutine glide_config

!=======================================================================

  subroutine glide_initialise(model)

    ! Initialise Glide model instance

    use glide_setup
    use glimmer_ncio
    use glide_velo, only: init_velo
    use glide_thck
    use glide_temp
    use glimmer_log
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glimmer_coordinates, only: coordsystem_new
    use glide_diagnostics, only: glide_init_diag
    use glide_bwater

    use parallel, only: distributed_grid

    type(glide_global_type), intent(inout) :: model     ! model instance

!TODO - Is glimmer_version_char still used?
!       Old Glide does not include this variable.
    character(len=100), external :: glimmer_version_char

    integer, parameter :: my_nhalo = 0   ! no halo layers for Glide dycore

!!!Old Glide has this:
!!!    call write_log(glimmer_version)  

    call write_log(trim(glimmer_version_char()))

    ! initialise scales
    call glimmer_init_scales

    ! scale parameters (some conversions to SI units)
    call glide_scale_params(model)

    ! set up coordinate systems

    ! Note: nhalo = 0 is included in call to distributed_grid to set other halo
    !  variables (lhalo, uhalo, etc.) to 0 instead of default values

!WHL - distributed_grid is not in old glide
      
    call distributed_grid(model%general%ewn, model%general%nsn,  &
                          nhalo_in=my_nhalo)

    model%general%ice_grid = coordsystem_new(0.d0,               0.d0, &
                                             model%numerics%dew, model%numerics%dns, &
                                             model%general%ewn,  model%general%nsn)

    model%general%velo_grid = coordsystem_new(model%numerics%dew/2.d0, model%numerics%dns/2.d0, &
                                              model%numerics%dew,      model%numerics%dns,      &
                                              model%general%ewn-1,     model%general%nsn-1)

    ! allocate arrays
    call glide_allocarr(model)

!TODO - Eliminate the bed softness parameter and set btrc to model%velowo%btrac_const in glide_velo?
    ! initialise bed softness to uniform parameter
    model%velocity%bed_softness = model%velowk%btrac_const

    ! set uniform basal heat flux (positive down)
    !NOTE: This value will be overridden if we read bheatflx from an input file 
    !      (model%options%gthf = 1) or compute it (model%options%gthf = 2)
    model%temper%bheatflx = model%paramets%geot

    ! compute sigma levels or load from external file
    ! (if not already read from config file)
    call glide_load_sigma(model,dummyunit)

    ! open all input files and forcing files
    call openall_in(model)

    ! read first time slice
    call glide_io_readall(model,model)

    ! write projection info to log
    call glimmap_printproj(model%projection)
   
    !WHL - Should have been read from glide_io_readall
    ! read lithot if required
!!    if (model%options%gthf > 0) then
!    if (model%options%gthf == GTHF_COMPUTE) then
!       call glide_lithot_io_readall(model,model)
!    end if

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first
    call init_isostasy(model)

    select case(model%options%whichrelaxed)

    case(RELAXED_TOPO_INPUT)   ! Supplied topography is relaxed
       model%isostasy%relx = model%geometry%topg
    case(RELAXED_TOPO_COMPUTE) ! Supplied topography is in equilibrium
                               !TODO - test case RELAXED_TOPO_COMPUTE
       call isos_relaxed(model)
    end select

    ! open all output files
    call openall_out(model)

    ! create glide variables
    call glide_io_createall(model, model)

!WHL - debug
!    print*, ' '
!    print*, 'Created Glide variables'
!    print*, 'max, min bheatflx (W/m2)=', maxval(model%temper%bheatflx), minval(model%temper%bheatflx)

    ! If a 2D bheatflx field is present in the input file, it will have been written 
    !  to model%temper%bheatflx.  For the case model%options%gthf = 0, we want to use
    !  a uniform heat flux instead.
    ! If no bheatflx field is present in the input file, then we default to the 
    !  prescribed uniform value, model%paramets%geot.

    if (model%options%gthf == GTHF_UNIFORM) then

       ! Check to see if this flux was present in the input file
       ! (by checking whether the flux is nonuniform over the domain)
       if (abs(maxval(model%temper%bheatflx) - minval(model%temper%bheatflx)) > 1.d-6) then  
          call write_log('Setting uniform prescribed geothermal flux')
          call write_log('(Set gthf = 1 to read geothermal flux field from input file)')
       endif

       ! set uniform basal heat flux (positive down)
       model%temper%bheatflx = model%paramets%geot

!WHL - debug
!       print*, ' '
!       print*, 'Use uniform bheatflx'
!       print*, 'max, min bheatflx (W/m2)=', maxval(model%temper%bheatflx), minval(model%temper%bheatflx)

    endif
 
    !TODO - Change subroutine names to glide_init_velo, glide_init_thck

    ! initialise velocity calc
    call init_velo(model)

!WHL - old glide has a call to init_temp, which is similar to glide_init_temp
!      but does not set the temperature or compute flwa until later call to timeevoltemp
!WHL - In old glide I added artm as a hotstart variable

    ! Initialize temperature field - this needs to happen after input file is
    !  read so we can assign artm (which could possibly be read in) if temp has not been input.
    !
    ! Note: If the temperature field has not been read already from an input or restart file, 
    !        then temperature is initialized by this subroutine based on model%options%temp_init.  
    !       If the temperature has been read already, this subroutine will *not* overwrite it.
  
    call glide_init_temp(model)

    ! Initialize basal hydrology model, if enabled
    call bwater_init(model)

    ! initialise thickness evolution calc
    call init_thck(model)

    if (model%options%gthf == GTHF_COMPUTE) then
!!       call glide_lithot_io_createall(model)  !WHL - Variables should have been created by glide_io_createall
       call init_lithot(model)
    end if

!WHL - This call will set the ice column temperature to artm as in old glide,
!       regardless of the value of model%options%temp_init
!      Commented out at least for now.  To reproduce results of old_glide, make sure
!       model%options%temp_init = TEMP_INIT_ARTM.
!!  if (oldglide) then
!!    if (model%options%hotstart.ne.1) then
!!       ! initialise Glen's flow parameter A using an isothermal temperature distribution
!!       call glide_temp_driver(model,0)
!!    endif
!!  endif  ! oldglide
    
!WHL - This option is disabled for now.
    ! *mb* added; initialization of basal proc. module
!!    if (model%options%which_bproc == BAS_PROC_FULLCALC .or. &
!!        model%options%which_bproc == BAS_PROC_FASTCALC) then
!!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!!                              model%numerics%dttem)
!!    end if      

    call glide_set_mask(model%numerics,                                   &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)
 
    ! calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)

    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

    ! initialise thckwk variables; used in timeders subroutine
    model%thckwk%olds(:,:,1) = model%geometry%thck(:,:)
    model%thckwk%olds(:,:,2) = model%geometry%usrf(:,:)

    ! initialise standard glide profiling
    call glide_prof_init(model)

    !TODO - Unclear on how subroutine register_model is used - Is it needed for serial code?
    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error

    call register_model(model)

    ! initialise model diagnostics                                                                                                                 \
                                                                                                                                                  
    call glide_init_diag(model)

!WHL - debug
!    print*, 'After glide_initialise:'
!    print*, 'max, min thck (m)=', maxval(model%geometry%thck)*thk0, minval(model%geometry%thck)*thk0
!    print*, 'max, min usrf (m)=', maxval(model%geometry%usrf)*thk0, minval(model%geometry%usrf)*thk0
!    print*, 'max, min artm =', maxval(model%climate%artm), minval(model%climate%artm)
!    print*, 'max, min temp =', maxval(model%temper%temp), minval(model%temper%temp)
!    print*, 'max, min flwa =', maxval(model%temper%flwa), minval(model%temper%flwa)

!    print*, ' '
!    print*, 'thck:'
!    do j = model%general%nsn, 1, -1
!       write(6,'(30f5.0)') thk0 * model%geometry%thck(:,j)
!    enddo
!    print*, ' '
!    print*, 'temp, k = 2:'
!    do j = model%general%nsn+1, 0, -1
!       write(6,'(32f5.0)') model%temper%temp(2,:,j)
!    enddo
!    print*, 'basal temp:'
!    do j = model%general%nsn+1, 0, -1
!       write(6,'(32f5.0)') model%temper%temp(model%general%upn,:,j)
!    enddo

  end subroutine glide_initialise

!=======================================================================

  subroutine glide_init_state_diagnostic(model)

    ! Calculate diagnostic variables for the initial model state
    ! This provides calculation of output fields at time 0
    ! This is analagous to glissade_diagnostic_variable_solve but is only 
    ! called from init.  The glide tstep routines take care of these calculations
    ! during time stepping.  
    ! Note that none of this is needed on a restart - this code ensures a complete 
    ! set of diagnostic output fields for the initial state.

    use glide_thck
    use glide_velo
    use glide_mask
    use glimmer_paramets, only: tim0
    use glimmer_physcon, only: scyr
    use glide_ground, only: glide_calve_ice
    use glide_bwater, only: calcbwat
    use glide_temp, only: glide_calcbmlt, glide_calcbpmp
    use glide_grid_operators

    type(glide_global_type), intent(inout) :: model     ! model instance

    integer :: i, j

    if (model%options%is_restart == RESTART_TRUE) then
       ! On a restart, just assign the basal velocity from uvel/vvel (which are restart variables)
       ! to ubas/vbas which are used by the temperature solver to calculate basal heating.
       ! During time stepping ubas/vbas are calculated by slipvelo during thickness evolution or below on a cold start.
       model%velocity%ubas = model%velocity%uvel(model%general%upn,:,:)
       model%velocity%vbas = model%velocity%vvel(model%general%upn,:,:)

    else   ! Only make the calculations on a cold start

    ! ------------------------------------------------------------------------ 
    ! ***Part 1: Make geometry consistent with calving law, if desired
    ! ------------------------------------------------------------------------       

       if (model%options%calving_init == CALVING_INIT_ON) then

          !WHL - debug
          print*, 'Calving at initialization, whichcalving =', model%options%whichcalving

          ! ------------------------------------------------------------------------ 
          ! Remove ice which should calve, depending on the value of whichcalving
          ! ------------------------------------------------------------------------ 
       
          !  Note: On a cold start, glide_calve_ice needs the mask to be calculated, but a call to 
          !  glide_set_mask occurs in glide_initialise and does not need to be repeated here.
          
          call glide_calve_ice(model%options%whichcalving,     &
                               model%geometry%thck,            &
                               model%isostasy%relx,            &
                               model%geometry%topg,            &
                               model%geometry%thkmask,         &
                               model%calving%marine_limit,     &
                               model%calving%calving_fraction, &
                               model%climate%eus,              &
                               model%calving%calving_thck)
          
          ! We now need to recalculate the mask because glide_calve_ice may have modified the geometry.
          call glide_set_mask(model%numerics,                                &
                              model%geometry%thck,  model%geometry%topg,     &
                              model%general%ewn,    model%general%nsn,       &
                              model%climate%eus,    model%geometry%thkmask,  &
                              model%geometry%iarea, model%geometry%ivol)
          
       endif  ! calving_init

       ! Compute total areas of grounded and floating ice
       call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                               model%geometry%thkmask,                        &
                               model%geometry%iareaf, model%geometry%iareag)

    ! ------------------------------------------------------------------------ 
    ! ***Part 2: Calculate geometry related fields
    ! ------------------------------------------------------------------------    

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------

       call glide_calclsrf(model%geometry%thck, model%geometry%topg, &
                        model%climate%eus,   model%geometry%lsrf)

       model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives
    !
    ! This call is needed here to make sure stagthck is calculated 
    ! the same way as in thck_lin_evolve/thck_nonlin_evolve
    ! ------------------------------------------------------------------------     

       call glide_prof_start(model,model%glide_prof%geomderv)

       call glide_geometry_derivs(model)   ! stagvarb, geomders as in old Glide

       call glide_prof_stop(model,model%glide_prof%geomderv)

       call glide_prof_start(model,model%glide_prof%ice_mask1)

    !TREY This sets local values of dom, mask, totpts, and empty
    !EIB! call veries between lanl and gc2, this is lanl version
    !magi a hack, someone explain what whichthck=5 does

!WHL - Modified this subroutine so that ice can accumulate in regions with
!      a small positive mass balance.

       call glide_thck_index(model%geometry% thck,      &
                          model%climate%  acab,      &
                          model%geometry% thck_index,  &
                          model%geometry% totpts,    &
                          .true.,                    &
                          model%geometry% empty)

       call glide_prof_stop(model,model%glide_prof%ice_mask1)


    ! ------------------------------------------------------------------------ 
    ! Part 3: Solve velocity
    ! ------------------------------------------------------------------------    

       ! initial value for flwa should already be calculated as part of glide_init_temp()
       ! calculate the part of the vertically averaged velocity field which solely depends on the temperature

        call velo_integrate_flwa(model%velowk,             &
                                 model%geomderv%stagthck,  &
                                 model%temper%flwa)

        ! Calculate diffusivity

       call velo_calc_diffu(model%velowk,              &
                            model%geomderv%stagthck,   &
                            model%geomderv%dusrfdew,   &
                            model%geomderv%dusrfdns,   &
                            model%velocity%diffu)

       ! If necessary, compute staggered variables required for basal traction calculation

       if (model%options%whichbtrc == BTRC_CONSTANT_BWAT) then

          !TODO - I think the next two calls are not needed, given that bwat should be in restart file for this option.

          ! Calculate basal melt rate --------------------------------------------------
          ! Note: For the initial state, we won't have values for ubas/vbas (unless they were 
          !       supplied in the input file) to get an initial guess of sliding heating.
          !       We could iterate on this, but for simplicity that is not done.

          call glide_calcbmlt(model, &
                              model%temper%temp, &
                              model%geometry%thck, &
                              model%geomderv%stagthck, &
                              model%geomderv%dusrfdew, &
                              model%geomderv%dusrfdns, &
                              model%velocity%ubas, &
                              model%velocity%vbas, &
                              model%temper%bmlt_ground, &
                              GLIDE_IS_FLOAT(model%geometry%thkmask))

          ! Note: calcbwat computes stagbwat
          call calcbwat(model, &
                        model%options%whichbwat, &
                        model%temper%bmlt_ground, &
                        model%temper%bwat, &
                        model%temper%bwatflx, &
                        model%geometry%thck, &
                        model%geometry%topg, &
                        model%temper%temp(model%general%upn,:,:), &
                        GLIDE_IS_FLOAT(model%geometry%thkmask), &
                        model%tempwk%wphi)


          ! This call is redundant for now, but is needed if the call to calcbwat is removed
          call stagvarb(model%temper%bwat, &
                        model%temper%stagbwat ,&
                        model%general%ewn, &
                        model%general%nsn)

       elseif (model%options%whichbtrc == BTRC_CONSTANT_TPMP) then

          call stagvarb(model%temper%temp(model%general%upn,1:model%general%ewn,1:model%general%nsn), &
                        model%temper%stagbtemp ,&
                        model%general%  ewn, &
                        model%general%  nsn)
       
          call glide_calcbpmp(model,               &
                              model%geometry%thck, &
                              model%temper%bpmp)
 
          call stagvarb(model%temper%bpmp, &
                        model%temper%stagbpmp ,&
                        model%general%  ewn, &
                        model%general%  nsn)

       endif   ! whichbtrc

       !------------------------------------------------------------------------ 
       ! Calculate basal traction factor
       !------------------------------------------------------------------------ 

       do j = 1, model%general%nsn-1
          do i = 1, model%general%ewn-1
             if (model%geomderv%stagthck(i,j)*thk0 < 1000.d0) then
                model%temper%stagbtemp(i,j) = model%temper%stagbpmp(i,j)
             else
                model%temper%stagbtemp(i,j) = -20.d0
             endif
          enddo
       enddo

       call calc_btrc(model,                    &
                      model%options%whichbtrc,  &
                      model%velocity%btrc)

       call slipvelo(model,                     &
                     0,                         &
                     model%velocity%btrc,       &
                     model%velocity%ubas,       &
                     model%velocity%vbas)
       

       ! Calculate velocity
       call velo_calc_velo(model%velowk,            model%geomderv%stagthck,  &
                           model%geomderv%dusrfdew, model%geomderv%dusrfdns,  &
                           model%temper%flwa,       model%velocity%diffu,     &
                           model%velocity%ubas,     model%velocity%vbas,      &
                           model%velocity%uvel,     model%velocity%vvel,      &
                           model%velocity%uflx,     model%velocity%vflx,      &
                           model%velocity%velnorm)    

    endif  ! if a restart


    ! MJH: I have left these calls outside of the restart if-construct so that there will
    ! always be a velnorm field calculated, which can be helpful for debugging.

    ! ------------------------------------------------------------------------ 
    ! Part 4: Calculate other diagnostic fields that depend on velocity
    ! ------------------------------------------------------------------------    

    ! ------------------------------------------------------------------------
    ! basal shear stress calculation
    ! ------------------------------------------------------------------------

    call calc_basal_shear(model%geomderv%stagthck,                          &
                          model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
                          model%velocity%tau_x,    model%velocity%tau_y)

    ! velocity norm
    model%velocity%velnorm = sqrt(model%velocity%uvel**2 + model%velocity%vvel**2)

!WHL - debug
       if (verbose_glide) then

          print*, ' '
          print*, 'stagthck:'
          do i = 1, model%general%ewn-1
             write(6,'(i7)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i3)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f7.2)',advance='no') model%geomderv%stagthck(i,j)*thk0
             enddo 
             print*, ' '
          enddo

          print*, ' '
          print*, 'diffu (m^2/yr):'
          do i = 1, model%general%ewn-1
             write(6,'(i8)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i3)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f8.0)',advance='no') -model%velocity%diffu(i,j) * vel0*len0*scyr
             enddo
             print*, ' '
          enddo

          print*, ' '
          print*, 'ubas:'
          do i = 1, model%general%ewn-1
             write(6,'(i7)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i4)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f7.2)',advance='no') model%velocity%uvel(model%general%upn,i,j)*(vel0*scyr)
             enddo
             print*, ' '
          enddo

          print*, ' '
          print*, 'vbas:'
          do i = 1, model%general%ewn-1
             write(6,'(i7)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i4)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f7.2)',advance='no') model%velocity%vvel(model%general%upn,i,j)*(vel0*scyr)
             enddo
             print*, ' '
          enddo
          
          print*, ' '
          print*, 'uvel, k = 1:'
          do i = 1, model%general%ewn-1
             write(6,'(i8)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i4)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f8.2)',advance='no') model%velocity%uvel(1,i,j) * (vel0*scyr)
             enddo 
             print*, ' '
          enddo

          print*, ' '
          print*, 'u=vvel, k = 1:'
          do i = 1, model%general%ewn-1
             write(6,'(i8)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i4)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f8.2)',advance='no') model%velocity%vvel(1,i,j) * (vel0*scyr)
             enddo 
             print*, ' '
          enddo

       endif  ! verbose_glide

  end subroutine glide_init_state_diagnostic

!=======================================================================

  subroutine glide_tstep_p1(model,time)

    ! Perform first part of time-step of an ice model instance:
    ! temperature advection, vertical conduction, and internal dissipation.

    use glide_thck
    use glide_velo
    use glide_temp
    use glide_mask
    use glimmer_paramets, only: tim0
    use glimmer_physcon, only: scyr
    use glide_grid_operators
    use glide_bwater

    type(glide_global_type), intent(inout) :: model     ! model instance
    real(dp),  intent(in)   :: time                     ! current time in years

    ! Update internal clock
    model%numerics%time = time  
    model%temper%newtemps = .false.

    model%thckwk%oldtime = model%numerics%time - (model%numerics%dt * tim0/scyr)

    call glide_prof_start(model,model%glide_prof%geomderv)

    ! Update geometric quantities: stagthck, dusrfdew/dns, dthckdew/dns

    call glide_geometry_derivs(model)  ! compute stagthck, dusrfdew/dns, dthckdew/dns

    call glide_prof_stop(model,model%glide_prof%geomderv)

    call glide_prof_start(model,model%glide_prof%ice_mask1)

    !WHL - Modified this subroutine so that ice can accumulate in regions with
    !      a small positive mass balance.

    call glide_thck_index(model%geometry% thck,        &
                          model%climate%  acab,        &
                          model%geometry% thck_index,  &
                          model%geometry% totpts,      &
                          .true.,                      &
                          model%geometry% empty)

    call glide_prof_stop(model,model%glide_prof%ice_mask1)

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 

    if (model%options%gthf == GTHF_COMPUTE) then
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glen's A, if necessary
    ! ------------------------------------------------------------------------ 

    ! Note: These times have units of years.
    !       dttem has scaled units, so multiply by tim0/scyr to convert to years

    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%dttem*tim0/scyr)) then

       call glide_prof_start(model,model%glide_prof%temperature)

       if (oldglide) then   ! compute vertical velocity in glide_tstep_p1 
                          ! In new glide, this is called in glide_tstep_p3
         
          call glide_velo_vertical(model)

       endif   ! oldglide = T

       ! temperature advection, vertical conduction, and internal dissipation

       call glide_temp_driver(model, model%options%whichtemp)

       model%temper%newtemps = .true.

       call glide_prof_stop(model,model%glide_prof%temperature)

       ! Update hydrology, if needed ------------------------------------------------
       call calcbwat(model, &
                     model%options%whichbwat, &
                     model%temper%bmlt_ground, &
                     model%temper%bwat, &
                     model%temper%bwatflx, &
                     model%geometry%thck, &
                     model%geometry%topg, &
                     model%temper%temp(model%general%upn,:,:), &
                     GLIDE_IS_FLOAT(model%geometry%thkmask), &
                     model%tempwk%wphi)

    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate basal traction factor
    ! ------------------------------------------------------------------------ 

    call calc_btrc(model,                    &
                   model%options%whichbtrc,  &
                   model%velocity%btrc)

!WHL - debug
!    print*, ' '
!    print*, 'After glide_tstep_p1:'
!    print*, 'max, min temp =', maxval(model%temper%temp), minval(model%temper%temp)
!    print*, 'max, min flwa =', maxval(model%temper%flwa), minval(model%temper%flwa)

!    print*, ' '
!    print*, 'temp, k = 2:'
!    do j = model%general%nsn+1, 0, -1
!       write(6,'(14f12.7)') model%temper%temp(2,3:16,j)
!    enddo
!    print*, 'basal temp:'
!    do j = model%general%nsn+1, 0, -1
!       write(6,'(14f12.7)') model%temper%temp(model%general%upn,3:16,j)
!    enddo

  end subroutine glide_tstep_p1

!=======================================================================

  subroutine glide_tstep_p2(model)

    ! Perform second part of time-step of an ice model instance:
    ! thickness evolution by one of several methods.

    use glide_thck
    use glide_velo
    use glide_temp
    use glide_mask
    use isostasy
    use glide_ground, only: glide_calve_ice

    type(glide_global_type), intent(inout) :: model    ! model instance

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 

    call glide_prof_start(model,model%glide_prof%ice_evo)

    select case(model%options%whichevol)

    case(EVOL_PSEUDO_DIFF) ! Use precalculated uflx, vflx -----------------------------------

       call thck_lin_evolve(model,model%temper%newtemps)

    case(EVOL_ADI) ! Use explicit leap frog method with uflx,vflx -------------------

       call stagleapthck(model,model%temper%newtemps)

    case(EVOL_DIFFUSION) ! Use non-linear calculation that incorporates velocity calc -----

       call thck_nonlin_evolve(model,model%temper%newtemps)

    end select

    call glide_prof_stop(model,model%glide_prof%ice_evo)

    ! ------------------------------------------------------------------------ 
    ! get new mask
    ! Note: A call to glide_set_mask is needed before glide_calve_ice.
    ! ------------------------------------------------------------------------ 

    call glide_prof_start(model,model%glide_prof%ice_mask2)

    !TODO - Calculate area and vol separately from glide_set_mask?

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

    call glide_prof_stop(model,model%glide_prof%ice_mask2)

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichcalving
    ! ------------------------------------------------------------------------ 

    !TODO - Some arguments for glide_calve_ice may not be needed.
    !       Old glide includes only arguments through model%climate%calving.

    call glide_calve_ice(model%options%whichcalving,     &
                         model%geometry%thck,            &
                         model%isostasy%relx,            &
                         model%geometry%topg,            &
                         model%geometry%thkmask,         &
                         model%calving%marine_limit,     &
                         model%calving%calving_fraction, &
                         model%climate%eus,              &
                         model%calving%calving_thck)

    ! Recalculate the mask following calving
    ! Note - This call to glide_set_mask is not in old Glide, but should have been.

 if (.not. oldglide) then   ! recalculate the thickness mask after calving
    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)
 endif   ! oldglide = F

 if (.not. oldglide) then   ! calculate area of floating and grounded ice
    call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                            model%geometry%thkmask,                        &
                            model%geometry%iareaf, model%geometry%iareag)
 endif   ! oldglide = F

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------

    call glide_prof_start(model,model%glide_prof%isos_water)

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       if (model%numerics%time >= model%isostasy%next_calc) then
          model%isostasy%next_calc = model%isostasy%next_calc + model%isostasy%period
          call isos_icewaterload(model)
          model%isostasy%new_load = .true.
       end if
    end if

    call glide_prof_stop(model,model%glide_prof%isos_water)
    
    ! ------------------------------------------------------------------------
    ! basal shear stress calculation
    ! ------------------------------------------------------------------------

! Old glide just passes 'model'

    call calc_basal_shear(model%geomderv%stagthck,                          &
                          model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
                          model%velocity%tau_x,    model%velocity%tau_y)

! not in old glide, but this is a useful diagnostic

    ! velocity norm
    model%velocity%velnorm = sqrt(model%velocity%uvel**2 + model%velocity%vvel**2)

!WHL - debug
!    print*, ' '
!    print*, 'After tstep_p2:'
!    print*, 'max, min thck (m)=', maxval(model%geometry%thck)*thk0, minval(model%geometry%thck)*thk0
!    print*, 'max, min usrf (m)=', maxval(model%geometry%usrf)*thk0, minval(model%geometry%usrf)*thk0
!    print*, 'max uvel, vvel =', maxval(model%velocity%uvel), maxval(model%velocity%vvel)

!    print*, ' '
!    print*, 'thck:'
!    do j = model%general%nsn, 1, -1
!       write(6,'(14f12.7)') thk0 * model%geometry%thck(3:16,j)
!    enddo
!    print*, 'sfc uvel:'
!    do j = model%general%nsn-1, 1, -1
!       write(6,'(14f12.7)') model%velocity%uvel(1,3:16,j)
!    enddo
!    print*, 'sfc vvel:'
!    do j = model%general%nsn-1, 1, -1
!       write(6,'(14f12.7)') model%velocity%vvel(1,3:16,j)
!    enddo

  end subroutine glide_tstep_p2

!=======================================================================

  subroutine glide_tstep_p3(model)

    ! Perform third part of time-step of an ice model instance:
    ! calculate isostatic adjustment and upper and lower ice surface

    use isostasy
    use glide_setup
    use glide_velo, only: glide_velo_vertical
    use glide_thck, only: glide_calclsrf
    implicit none

    type(glide_global_type), intent(inout) :: model     ! model instance
    
    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 

    call glide_prof_start(model,model%glide_prof%isos)

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call isos_compute(model)
    end if

    call glide_prof_stop(model,model%glide_prof%isos)

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------

    call glide_calclsrf(model%geometry%thck, model%geometry%topg, &
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

    !TODO - Move timecounter to a driver routine?
    !CESM Glimmer code has this after the netCDF write.

    ! increment time counter
    model%numerics%timecounter = model%numerics%timecounter + 1

    !TODO - Combine these timeders and vert velo calls into a subroutine?

    ! For exact restart, compute wgrd here and write it to the restart file.
    ! (This is easier than writing thckwk quantities to the restart file.)

 if (.not. oldglide) then  ! compute vertical velocity in glide_tstep_p3

    ! compute vertical velocity
         
    call t_startf('vertical_velo')

    call glide_velo_vertical(model)

    call t_stopf('vertical_velo')

 endif  ! oldglide = F

 !WHL - Moved netCDF output to simple_glide
 !!       call glide_io_writeall(model,model)

  end subroutine glide_tstep_p3

!=======================================================================

end module glide
