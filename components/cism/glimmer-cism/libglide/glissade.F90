!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

! WJS (1-30-12): The following (turning optimization off) is needed as a workaround for an
! xlf compiler bug, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

!CLEANUP - glissade.F90
!
! NOTE: MJH Lines that start with !### are ones I have identified to be deleted.
!
! This is a new module, originally copied from glide.F90 (William Lipscomb, June 2012)
! Removed SIA-specific code, leaving only the HO code with remapping transport
! Whenever possible, parallel_halo updates should go in this module rather
!  than at lower levels.
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glissade.f90 - part of the Glimmer-CISM ice model        + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade

  ! Driver for Glissade (parallel, higher-order) dynamical core

  use glimmer_global, only: dp
  use glimmer_log
  use glide_types
  use glide_io
  use glide_lithot
  use glimmer_config

  implicit none

  integer :: &
     ntracer = 1   ! number of tracers to transport (just temperature for now)
                   ! BDM change ntracer from parameter so it's not a constant
                   !TODO - Declare elsewhere?

  integer, private, parameter :: dummyunit=99

contains

!=======================================================================

! Note: There is no glissade_config subroutine; glide_config works for all dycores.

!=======================================================================

  subroutine glissade_initialise(model)

    ! initialise Glissade model instance

!TODO - Are all of these needed?
    use parallel
    use glide_stop, only: register_model
    use glide_setup
    use glimmer_ncio
    use glide_velo, only: init_velo  !TODO - Remove this
    use glissade_temp, only: glissade_init_temp
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glide_ground
    use glide_thck, only : glide_calclsrf
    use glam_strs2, only : glam_velo_init
    use glimmer_coordinates, only: coordsystem_new

    use glissade_velo_higher, only: glissade_velo_higher_init

!!    use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    !TODO - build glimmer_vers file or put this character elsewhere?
    character(len=100), external :: glimmer_version_char

!WHL - debug
!    logical, parameter :: test_parallel = .false.   ! if true, call test_parallel subroutine
    integer :: i, j, nx, ny

!WHL - for artificial adjustment to ismip-hom surface elevation
    logical, parameter :: ismip_hom_adjust_usrf = .false.
    real(dp) :: usrf_ref

    call write_log(trim(glimmer_version_char()))

    ! initialise scales
    call glimmer_init_scales

    ! scale parameters (some conversions to SI units)
    call glide_scale_params(model)

    ! set up coordinate systems
    ! time to change to the parallel values of ewn and nsn

    call distributed_grid(model%general%ewn,model%general%nsn)

    model%general%ice_grid = coordsystem_new(0.d0,               0.d0,               &
                                             model%numerics%dew, model%numerics%dns, &
                                             model%general%ewn,  model%general%nsn)

    !TODO - Change 2. to 2.d0
    model%general%velo_grid = coordsystem_new(model%numerics%dew/2., model%numerics%dns/2., &
                                              model%numerics%dew,    model%numerics%dns,    &
                                              model%general%ewn-1,   model%general%nsn-1)


    ! allocate arrays
    call glide_allocarr(model)

    !TODO - May be able to eliminate the bed softness parameter 
    !       and set btrc to model%velowo%btrac_const in glide_velo

    ! initialise bed softness to uniform parameter
    model%velocity%bed_softness = model%velowk%btrac_const

    ! set uniform basal heat flux (positive down)
    model%temper%bheatflx = model%paramets%geot

    ! compute sigma levels or load from external file
    ! (if not already read from config file)
    call glide_load_sigma(model,dummyunit)

    ! open all input files
    call openall_in(model)

    ! and read first time slice
    call glide_io_readall(model,model)

    ! Write projection info to log
    call glimmap_printproj(model%projection)

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first
    call init_isostasy(model)

    select case(model%options%whichrelaxed)
    case(RELAXED_TOPO_INPUT)   ! Supplied topography is relaxed
       model%isostasy%relx = model%geometry%topg
    case(RELAXED_TOPO_COMPUTE) ! Supplied topography is in equilibrium
                               !TODO - Test this case
       call not_parallel(__FILE__,__LINE__)
       call isos_relaxed(model)
    end select

    ! open all output files
    call openall_out(model)

    ! create glide variables
    call glide_io_createall(model)

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

    endif

    ! initialise glissade components

    !TODO - Most of what's done in init_velo is needed for SIA only
    !       Can remove this call provided velowk is not used elsewhere (e.g., to call wvelintg)
    call init_velo(model)

    !BDM - Just call glissade_init_temp, which will call glissade_init_enthalpy
    !      if necessary

    call glissade_init_temp(model) 

    if (model%options%gthf == GTHF_COMPUTE) then
       call not_parallel(__FILE__,__LINE__)
       call init_lithot(model)
    end if

    if (model%options%whichdycore == DYCORE_GLAM ) then  ! glam finite-difference

       call glam_velo_init(model%general%ewn,    model%general%nsn,  &
                           model%general%upn,                        &
                           model%numerics%dew,   model%numerics%dns, &
                           model%numerics%sigma)

    elseif (model%options%whichdycore == DYCORE_GLISSADE ) then  ! glissade finite-element

       call glissade_velo_higher_init

    endif

!WHL - This option is disabled for now.
    ! *mb* added; initialization of basal proc. module
!!    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
!!        model%options%which_bmod == BAS_PROC_FASTCALC) then        
!!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!!                              model%numerics%ntem)
!!    end if      

    ! calculate mask
       call glide_set_mask(model%numerics,                                &
                           model%geometry%thck,  model%geometry%topg,     &
                           model%general%ewn,    model%general%nsn,       &
                           model%climate%eus,    model%geometry%thkmask,  &
                           model%geometry%iarea, model%geometry%ivol)

    !TODO- Not sure why this needs to be called here.
    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

    ! and calculate lower and upper ice surface

    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

!WHL - The following is a hack to ensure that usrf is uniform (to double precision) for 
!      a given value of i for the ISMIP-HOM tests.  We take one value in each column as 
!      the benchmark and set all other values in that column to the same value.
!      Then we correct the thickness to ensure that lsrf + thck = usrf to double precision.
!      
!      A better way would be to read double-precision data from the input files.

    if (ismip_hom_adjust_usrf) then
       do i = 1, model%general%ewn
          usrf_ref = model%geometry%usrf(i,nhalo+1)
          do j = 1, model%general%nsn
             model%geometry%usrf(i,:) = usrf_ref  ! same usrf for all j in this column
          enddo
       enddo
       model%geometry%thck(:,:) = model%geometry%usrf(:,:) - model%geometry%lsrf(:,:) 
    endif

    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)

     !WHL - debug
!    if (test_parallel) then
!       call glissade_test_parallel (model)
!       call parallel_finalise
!    endif
     
  end subroutine glissade_initialise
  
!=======================================================================

  subroutine glissade_tstep(model, time, no_write)

    ! Perform time-step of an ice model instance with glissade dycore
    !TODO - Reorganize to put isostasy and calving at start of step?

    use parallel
!!    use glimmer_horiz_bcs, only: horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns, &
!!                                 horiz_bcs_unstag_scalar

    use glimmer_paramets, only: tim0, len0, vel0, thk0
    use glimmer_physcon, only: scyr
    use glissade_temp, only: glissade_temp_driver
    use glide_mask, only: glide_set_mask, calc_iareaf_iareag
    use glide_ground, only: glide_marinlim
    use glide_grid_operators
    use isostasy
    use glissade_enthalpy
    use glissade_transport, only: glissade_transport_driver

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    real(dp),  intent(in)   :: time         !*FD Current time in years

    !TODO - Remove this argument; it is not used 
    logical, optional, intent(in) :: no_write

    logical nw

    ! --- Local Variables ---

    integer :: sc  ! subcycling index

    ! temporary thck array in SI units (m)
    !TODO - SCALING - Remove if and when scaling is removed
    real(dp), dimension(model%general%ewn,model%general%nsn) :: thck_unscaled

    ! temporary variables needed to reset geometry for the EVOL_NO_THICKNESS option
    real(dp), dimension(model%general%ewn,model%general%nsn) :: thck_old
    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) :: stagthck_old

    ! temporary bmlt array
    real(dp), dimension(model%general%ewn,model%general%nsn) :: &
       bmlt_continuity  ! = bmlt if basal mass balance is included in continuity equation
                        ! else = 0
    
    !WHL - debug
    integer :: j

    ! ========================

    ! Update internal clock
    model%numerics%time = time  
    model%numerics%timecounter = model%numerics%timecounter + 1
    model%temper%newtemps = .false.

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    !TODO Not sure if this is in the right place.  G1=f(G0,T0) and T1=g(G0,T0)  
    !     If we update G1 now, then we will be doing T1=g(G1,T0).
    if (model%options%gthf == GTHF_COMPUTE) then
       call not_parallel(__FILE__,__LINE__)
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glen's A, if necessary
    ! Vertical diffusion and strain heating only; no advection
    ! ------------------------------------------------------------------------ 

    ! Note: These times have units of years.

    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then

!HALO TODO - Modify glissade_temp_driver to compute over locally owned cells only?

      call t_startf('glissade_temp_driver')
      call glissade_temp_driver(model, model%options%whichtemp)
      call t_stopf('glissade_temp_driver')

       model%temper%newtemps = .true.

    end if

    ! ------------------------------------------------------------------------ 
    ! Halo updates
    ! ------------------------------------------------------------------------ 

    call parallel_halo(model%temper%bwat)    !HALO: not sure if this is needed

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
    ! MJH: This now uses velocity from the previous time step, which is appropriate for a Forward Euler time-stepping scheme
    ! WHL: We used to have EVOL_NO_THICKNESS = -1 as a Glide option, used to hold the ice surface elevation fixed during CESM runs.  
    !      This option has been replaced by a Glint option, evolve_ice.
    !      We now have EVOL_NO_THICKESS = 5 as a glam/glissade option.  It is used to hold the ice surface elevation fixed
    !       while allowing temperature to evolve, which can be useful for model spinup.  This option might need more testing.

    select case(model%options%whichevol)

       case(EVOL_INC_REMAP, EVOL_NO_THICKNESS) 

       ! Use incremental remapping scheme for advecting ice thickness ---
       ! (and temperature too, if whichtemp = TEMP_PROGNOSTIC)
       ! MJH: I put the no thickness evolution option here so that it is still possible 
       ! (but not required) to use IR to advect temperature when thickness evolution is turned off.

       ! TODO  MJH If we really want to support no evolution, then we may want to implement it so that IR does not occur 
       !       at all - right now a run can fail because of a CFL violation in IR even if evolution is turned off.  Do we want
       !       to support temperature evolution without thickness evolution?  If so, then the current implementation may be preferred approach.

       if (model%options%whichevol == EVOL_NO_THICKNESS) then
          ! store old thickness
          thck_old = model%geometry%thck
          stagthck_old = model%geomderv%stagthck
       endif

      call t_startf('inc_remap_driver')

       if (main_task) then
          print *, ' '
          print *, 'Compute dH/dt'
       endif

         !WHL - Testing a new subroutine that updates all the key scalars (thck, temp, etc.) at once
         !TODO - Do we need updates of lsrf, usrf, or topg?

!!       call parallel_halo_scalars(model%geometry%thck,   &
!!                                  model%temper%temp)


!       Velocity update may be needed if velo was not updated in halo at the end of the previous diagnostic solve
!        (just to be on the safe side).

        ! Halo updates for velocities, thickness and tracers
      call t_startf('new_remap_halo_upds')

       call staggered_parallel_halo(model%velocity%uvel)
!       call horiz_bcs_stag_vector_ew(model%velocity%uvel)

       call staggered_parallel_halo(model%velocity%vvel)
!       call horiz_bcs_stag_vector_ns(model%velocity%vvel)

       call parallel_halo(model%geometry%thck)
!       call horiz_bcs_unstag_scalar(model%geometry%thck)

       if (model%options%whichtemp == TEMP_PROGNOSTIC) then
          call parallel_halo(model%temper%temp)
!          call horiz_bcs_unstag_scalar(model%temper%temp)
       endif

       if (model%options%whichtemp == TEMP_ENTHALPY) then
          call parallel_halo(model%temper%temp)
          call parallel_halo(model%temper%waterfrac)
!          call horiz_bcs_unstag_scalar(model%temper%temp)
!          call horiz_bcs_unstag_scalar(model%temper%waterfrac)
       endif

      call t_stopf('new_remap_halo_upds')

      call t_startf('glissade_transport_driver')

       !TODO  It would be less confusing to just store the subcycling dt in a local/module variable - 
       !       really only needs to be calculated once on init

       model%numerics%dt = model%numerics%dt / model%numerics%subcyc

       if (model%options%basal_mbal == BASAL_MBAL_CONTINUITY) then    ! include bmlt in continuity equation
          bmlt_continuity(:,:) = model%temper%bmlt(:,:) * thk0/tim0   ! convert to m/s
       else                                                           ! do not include bmlt in continuity equation
          bmlt_continuity(:,:) = 0.d0
       endif

       ! Call the transport driver.
       ! Note: This subroutine assumes SI units:
       !       * dt (s)
       !       * dew, dns, thck (m)
       !       * uvel, vvel, acab, blmt (m/s)
       !       Since thck has intent(inout), we create and pass a temporary array with units of m.

       !TODO: Test this driver for 1st order upwind transport.

       do sc = 1 , model%numerics%subcyc
          if (model%numerics%subcyc > 1) write(*,*) 'Subcycling transport: Cycle ',sc

          if (model%options%whichtemp == TEMP_PROGNOSTIC) then  ! Use IR to transport thickness, temperature
                                                                ! (and other tracers, if present)
                                                                ! Note: We are passing arrays in SI units.

             ! temporary thickness array in SI units (m)                               
             thck_unscaled(:,:) = model%geometry%thck(:,:) * thk0

             call glissade_transport_driver(model%numerics%dt * tim0,                             &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      &
                                            nhalo,                     ntracer,                   &
                                            model%velocity%uvel(:,:,:) * vel0,                    &
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            thck_unscaled(:,:),                                   &
                                            dble(model%climate%acab(:,:)) * thk0/tim0,            &
                                            bmlt_continuity(:,:),                                 &
                                            model%temper%temp(:,:,:) )

             ! convert thck back to scaled units
             model%geometry%thck(:,:) = thck_unscaled(:,:) / thk0

          elseif (model%options%whichtemp == TEMP_ENTHALPY) then  ! Use IR to transport thickness, temperature,
                                                                ! and waterfrac.  Also set ntracer = 3
                                                                ! Note: We are passing arrays in SI units.

             ! temporary thickness array in SI units (m)                               
             thck_unscaled(:,:) = model%geometry%thck(:,:) * thk0

             ! BDM Set ntracer = 3 since temp and waterfrac need to be passed (and maybe ice age).
             !     If no ice age is present, a dummy array will be passed for tracer = 2
             ntracer = 3

             call glissade_transport_driver(model%numerics%dt * tim0,                             &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      &
                                            nhalo,                     ntracer,                 &
                                            model%velocity%uvel(:,:,:) * vel0,                    &
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            thck_unscaled(:,:),                                   &
                                            dble(model%climate%acab(:,:)) * thk0/tim0,            &
                                            bmlt_continuity(:,:),                                 &
                                            model%temper%temp(:,:,:),                             & 
                                            model%temper%waterfrac(:,:,:)  )

             ! convert thck back to scaled units
             model%geometry%thck(:,:) = thck_unscaled(:,:) / thk0
         
          !TODO - Will we continue to support this option?  May not be needed.

          else  ! Use IR to transport thickness only
                ! Note: In glissade_transport_driver, the ice thickness is transported layer by layer,
                !       which is inefficient if no tracers are being transported.  (It would be more
                !       efficient to transport thickness in one layer only, using a vertically
                !       averaged velocity.)  But this option probably will not be used in practice;
                !       it is left in the code just to ensure backward compatibility with an
                !       older remapping scheme for transporting thickness only.

             ! temporary thickness array in SI units (m)                               
             thck_unscaled(:,:) = model%geometry%thck(:,:) * thk0

             call glissade_transport_driver(model%numerics%dt * tim0,                             &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      &
                                            nhalo,                     ntracer,                   &
                                            model%velocity%uvel(:,:,:) * vel0,                    & 
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            thck_unscaled(:,:),                                   &
                                            dble(model%climate%acab(:,:)) * thk0/tim0,            & 
                                            bmlt_continuity(:,:) )

             ! convert thck back to scaled units
             model%geometry%thck(:,:) = thck_unscaled(:,:) / thk0

          endif  ! whichtemp

!WHL - debug
!    print*, ' '
!    print*, 'After glissade_transport_driver:'
!    print*, 'max, min thck (m)=', maxval(model%geometry%thck)*thk0, minval(model%geometry%thck)*thk0
!    print*, 'max, min acab (m/yr) =', maxval(model%climate%acab)*scale_acab, minval(model%climate%acab)*scale_acab
!    print*, 'max, min artm =', maxval(model%climate%artm), minval(model%climate%artm)
!    print*, 'thklim =', model%numerics%thklim * thk0
!    print*, 'max, min temp =', maxval(model%temper%temp), minval(model%temper%temp)
!    print*, ' '
!    print*, 'thck:'
!    do j = model%general%nsn, 1, -1
!       write(6,'(23f7.2)') model%geometry%thck(7:29,j)*thk0
!    enddo

          ! Update halos of modified fields

         call t_startf('after_remap_haloupds')

         call parallel_halo(model%geometry%thck)
!         call horiz_bcs_unstag_scalar(model%geometry%thck)

         call parallel_halo(model%temper%temp)
!         call horiz_bcs_unstag_scalar(model%temper%temp)

          ! Halo updates of other tracers, if present, would need to go here
         if (model%options%whichtemp == TEMP_ENTHALPY) then
            call parallel_halo(model%temper%waterfrac)
!            call horiz_bcs_unstag_scalar(model%temper%waterfrac)
         endif

         call t_stopf('after_remap_haloupds')

       enddo     ! subcycling

!TODO: Don't divide and multiply this variable
       model%numerics%dt = model%numerics%dt * model%numerics%subcyc

      call t_stopf('glissade_transport_driver')

      call t_stopf('inc_remap_driver')

       if (model%options%whichevol == EVOL_NO_THICKNESS) then
          ! restore old thickness
          model%geometry%thck = thck_old
          model%geomderv%stagthck = stagthck_old
       endif

    end select

!HALO TODO - What halo updates needed here?
!       We could put the various geometry halo updates here, after thickness and temperature evolution.
!       This would include thck and tracer at a minimum.
!       Also include topg if basal topography is evolving.
!       
!       TODO: Make sure a call to calc_flwa (based on post-remap temperature) 
!             is not needed for glide_marinlim.  Should be based on post-remap temperature.

    call parallel_halo(model%geometry%topg)
!    call horiz_bcs_unstag_scalar(model%geometry%topg)

    ! --- Calculate updated mask because marinlim calculation needs a mask.

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask)


    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 

!HALO - Look at marinlim more carefully and see which fields need halo updates before it is called.
!       It appears that marinlim only needs the halo of thkmask for case 5 (which was removed).  
!       If that case is removed, a thkmask halo update does not need to occur here.

    ! TODO: glide_set_mask includes a halo update of model%geometry%thkmask; move it here?
    call parallel_halo(model%geometry%thkmask) 
!    call horiz_bcs_unstag_scalar(model%geometry%thkmask)

    call parallel_halo(model%isostasy%relx)
!    call horiz_bcs_unstag_scalar(model%isostasy%relx)

    call glide_marinlim(model%options%whichmarn, &
         model%geometry%thck,  &
         model%isostasy%relx,      &
         model%geometry%topg,  &
         model%geometry%thkmask,    &
         model%numerics%mlimit,     &
         model%numerics%calving_fraction, &
         model%climate%eus,         &
         model%climate%calving,  &
         model%ground, &
         model%numerics%dew,    &
         model%numerics%dns, &
         model%general%nsn, &
         model%general%ewn)
         ! model%geometry%usrf) not used in routine

!HALO TODO: Test the effect of these updates with a nonzero calving field

    ! halo updates

    call parallel_halo(model%geometry%thck)    ! Updated halo values of thck are needed below in calc_lsrf
!    call horiz_bcs_unstag_scalar(model%geometry%thck)   

!TODO - Remove this call to glide_set_mask?
!      This subroutine is called at the beginning of glissade_velo_driver,
!       so a call here is not needed for the velo diagnostic solve.
!      The question is whether it is needed for the isostasy.
!      And the isostasy may itself be in the wrong place.

    ! --- marinlim adjusts thickness for calved ice.  Therefore the mask needs to be recalculated.
    ! --- This time we want to calculate the optional arguments iarea and ivol because thickness 
    ! --- will not change further during this time step.

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

!HALO TODO - glide_set_mask includes a halo update of model%geometry%thkmask at end of call
!       That update should be moved here if needed later (but may not be needed).

    ! call parallel_halo(model%geometry%thkmask) in previous glide_set_mask call
    ! call horiz_bcs_unstag_scalar(model%geometry%thkmask)

    ! --- Calculate area of ice that is floating and grounded.
    !TODO This subroutine does not use iarea - remove from the call/subroutine.
    !TODO May want to only calculate iarea, iareaf, iareag in glide_write_diag() and remove those calculations here.  

    call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                            model%geometry%iarea,  model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

    !TODO - Need a global sum here? (currently done inside calc_iareaf_iareag)

!TODO These isostasy calls may be in the wrong place.  
! Consider for a forward Euler time step:
! With a relaxing mantle model, topg is a prognostic (time-evolving) variable (I think):
! topg1 = f(topg0, thk0, ...) 
! However, for a fluid mantle where the adjustment is instantaneous, topg is a diagnostic variable 
!(comparable to calculating floatation height of ice in the ocean):
! topg1 = f(thk1)
! In either case, the topg update should be separate from the thickness evolution (because thk1 = f(thk0, vel0=g(topg0,...)).
! However, if the isostasy calculation needs topg0, the icewaterload call should be made BEFORE thck is updated.  
! If the isostasy calculation needs topg1, the icewaterload call should be made AFTER thck is updated.  
! Also, we need think carefully about when marinlim, usrf, lsrf, derivatives should be calculated relative to the topg update via isostasy.

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------

    if (model%options%isostasy == ISOSTASY_COMPUTE) then

       call not_parallel(__FILE__, __LINE__)

       if (model%numerics%time >= model%isostasy%next_calc) then
          model%isostasy%next_calc = model%isostasy%next_calc + model%isostasy%period
          call isos_icewaterload(model)
          model%isostasy%new_load = .true.
       end if
    end if
   
      ! calculate isostatic adjustment and upper and lower ice surface

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 

    !TODO - Test the local isostasy schemes in the parallel model.
    !       The elastic lithosphere scheme is not expected to work in parallel.

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call not_parallel(__FILE__, __LINE__)
       call isos_compute(model)
    end if

    ! ------------------------------------------------------------------------
    ! Calculate diagnostic variables, including velocity
    ! ------------------------------------------------------------------------

    call glissade_diagnostic_variable_solve(model)

    !TODO - Any halo updates needed at the end?  

  end subroutine glissade_tstep

!=======================================================================
!MJH added this diagnostic solve subroutine so it can be called from init.  

  subroutine glissade_diagnostic_variable_solve(model) 

     ! Solve diagnostic (not time-dependent) variables.  This is needed at the end of each time step once the 
     !  prognostic variables (thickness, tracers) have been updated.  
     ! It is also needed to fill out the initial state from the fields that have been read in.

    use parallel

    use glimmer_paramets, only: tim0, len0, vel0
    use glimmer_physcon, only: scyr
    use glide_thck, only: glide_calclsrf
    use glissade_temp, only: glissade_calcflwa
    use glam_velo, only: glam_velo_driver
    use glissade_velo, only: glissade_velo_driver
    use glide_stress, only : glide_calcstrsstr
!!    use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar, horiz_bcs_stag_scalar, horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns

    use glam_grid_operators, only: glam_geometry_derivs

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Local variables

!WHL - debug
    integer :: i, j

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 1. First part of diagnostic solve: 
    !    Now that advection is done, update geometry- and temperature-related 
    !    diagnostic fields that are needed for the velocity solve.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    !WHL - Moved this calculation here from glissade_temp, since flwa is a diagnostic variable.

    ! Calculate Glen's A --------------------------------------------------------
    !
    ! Note: because flwa is not a restart variable in glissade, no check is included 
    !       here for whether to calculate it on initial time (as is done in glide).
    ! 
    ! Note: We are passing in only vertical elements (1:upn-1) of the temp array,
    !       so that it has the same vertical dimensions as flwa.

    ! BDM - adding a call for whichtemp = TEMP_ENTHALPY that includes waterfrac as input
    ! TODO - May be OK to use a single call (with waterfrac) for either TEMP option.

    if (model%options%whichtemp == TEMP_ENTHALPY) then

       call glissade_calcflwa(model%numerics%stagsigma,    &
                              model%numerics%thklim,       &
                              model%temper%flwa,           &
                              model%temper%temp(1:model%general%upn-1,:,:),  &
                              model%geometry%thck,         &
                              model%paramets%flow_factor,  &
                              model%paramets%default_flwa, &
                              model%options%whichflwa,      &
                              model%temper%waterfrac(:,:,:))
    else

       call glissade_calcflwa(model%numerics%stagsigma,    &
                              model%numerics%thklim,       &
                              model%temper%flwa,           &
                              model%temper%temp(1:model%general%upn-1,:,:),  &
                              model%geometry%thck,         &
                              model%paramets%flow_factor,  &
                              model%paramets%default_flwa, &
                              model%options%whichflwa)
    endif

    ! Halo update for flwa
    call parallel_halo(model%temper%flwa)

    ! ------------------------------------------------------------------------
    ! Halo updates for ice topography and thickness
    !
    !WHL - Note the optional argument periodic_offset_ew for topg.
    !      This is for ismip-hom experiments. A positive EW offset means that 
    !       the topography in west halo cells will be raised, and the topography 
    !       in east halo cells will be lowered.  This ensures that the topography
    !       and upper surface elevation are continuous between halo cells
    !       and locally owned cells at the edge of the global domain.
    !      In other cases (anything but ismip-hom), periodic_offset_ew = periodic_offset_ns = 0, 
    !       and this argument will have no effect.
    ! ------------------------------------------------------------------------

    call parallel_halo(model%geometry%thck)
    call parallel_halo(model%geometry%topg, periodic_offset_ew = model%numerics%periodic_offset_ew)

    ! ------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    ! ------------------------------------------------------------------------

    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       & 
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))


    ! ------------------------------------------------------------------------
    ! Update some geometry derivatives
    ! ------------------------------------------------------------------------
    !TODO - The fields computed by glam_geometry_derivs are not required by
    !       the glissade solver (which computes them internally).  
    !       However, some of the fields (stagthck, dusrfdew and dusrfdns) 
    !       are needed during the next timestep by glissade_temp
    !       if we're doing shallow-ice dissipation.  If dissipation is always
    !       higher-order, then I don't think this call is needed here.
    !       (The glam_velo driver includes its own call to glam_geometry_derivs.) 

    call glam_geometry_derivs(model)

    !WHL - Moved glam-specific geometry calculations to glam_velo_driver in glam_velo.F90.

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 2. Second part of diagnostic solve: 
    !    Now that geometry- and temperature-related diagnostic fields are updated, 
    !    solve velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! Do not solve velocity for initial time on a restart because that breaks an exact restart.

    if ( (model%options%is_restart == RESTART_TRUE) .and. &
         (model % numerics % time == model % numerics % tstart) ) then
  
       call write_log('Using uvel, vvel from restart file at initial time')

    else

       ! If this is not a restart or we are not at the initial time, then proceed normally.

       if ( (model % numerics % time == model % numerics % tstart) .and. &
         ( (maxval(abs(model%velocity%uvel))/=0.0d0) .or. & 
           (maxval(abs(model%velocity%vvel))/=0.0d0) ) ) then
          ! If velocity was input and this is NOT a restart, then use the input field as the first guess at the initial time.
          ! This happens automatically, but let the user know.
          ! Using this value versus not will only change the answer within the tolerance of the nonlinear solve.  
          ! If a user already has a good guess from a previous run, they may wish to start things off with it to speed the initial solution.
          call write_log('Using uvel, vvel from input file as initial guess at initial time.')
          call write_log('If this is not desired, please remove those fields from the input file.')
       endif

       if (main_task) then
          print *, ' '
          print *, 'Compute higher-order ice velocities, time =', model%numerics%time
       endif


       !WHL - Broke up into separate calls to glam_velo_driver and glissade_velo_driver
       !      Previously had a single call to glissade_velo_driver

       if (model%options%whichdycore == DYCORE_GLAM) then    ! glam finite-difference dycore

         call t_startf('glam_velo_driver')
          call glam_velo_driver(model)
         call t_stopf('glam_velo_driver')

       else

         call t_startf('glissade_velo_driver')
          call glissade_velo_driver(model)
         call t_stopf('glissade_velo_driver')

       endif  ! whichdycore
 
    endif     ! is_restart

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 3. Third part of diagnostic solve: 
    ! Now that velocity is solved, calculate any diagnostic fields that are
    ! a function of velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! compute the velocity norm (for diagnostic output)

    model%velocity%velnorm = sqrt(model%velocity%uvel**2 + model%velocity%vvel**2)

    ! WHL - Copy uvel and vvel to arrays uvel_icegrid and vvel_icegrid.
    !       These arrays have horizontal dimensions (nx,ny) instead of (nx-1,ny-1).
    !       Thus they are better suited for I/O if we have periodic BC,
    !        where the velocity field we are solving for has global dimensions (nx,ny).
    !       Since uvel and vvel are not defined for i = nx or for j = ny, the
    !        uvel_icegrid and vvel_icegrid arrays will have values of zero at these points.
    !       But these are halo points, so when we write netCDF I/O it shouldn't matter;
    !        we should have the correct values at physical points.
    
    model%velocity%uvel_icegrid(:,:,:) = 0.d0
    model%velocity%vvel_icegrid(:,:,:) = 0.d0

    do j = 1, model%general%nsn-1
       do i = 1, model%general%ewn-1
          model%velocity%uvel_icegrid(:,i,j) = model%velocity%uvel(:,i,j)
          model%velocity%vvel_icegrid(:,i,j) = model%velocity%vvel(:,i,j)             
       enddo
    enddo
        
!TODO - Don't think we need to update ubas, vbas, or velnorm,
!       because these can be derived directly from the 3D uvel and vvel arrays

    call staggered_parallel_halo(model%velocity%velnorm)
    !       call horiz_bcs_stag_scalar(model%velocity%velnorm)
    call staggered_parallel_halo(model%velocity%ubas)
    !       call horiz_bcs_stag_vector_ew(model%velocity%ubas)
    call staggered_parallel_halo(model%velocity%vbas)
    !       call horiz_bcs_stag_vector_ns(model%velocity%vbas)

    call parallel_halo(model%stress%efvs)
    !       call horiz_bcs_unstag_scalar(model%stress%efvs)

    !Tau is calculated in glide_stress and initialized in glide_types.

    call glide_calcstrsstr( model )       !*sfp* added for populating stress tensor w/ HO fields

    ! Includes halo updates of 
    ! model%stress%tau%xx, model%stress%tau%yy, model%stress%tau%xy,
    ! model%stress%tau%scalar, model%stress%tau%xz, model%stress%tau%yz

    !HALO TODO - If the stress%tau halo updates are needed, they should go here (in glissade.F90)
    !            But I think they are not needed.

    ! --- A calculation of wvel could go here if we want to calculate it.
    ! --- For now, such a calculation is not needed.

  end subroutine glissade_diagnostic_variable_solve

!=======================================================================

  ! Not currently called; may not be worth having a separate subroutine for this

  subroutine parallel_halo_scalars(thck,     temp,   &
                                   lsrf,     usrf,   &
                                   topg,     tracers)

    ! Do parallel halo updates for the main scalar state variables

    use parallel

    real(dp), intent(inout), dimension(:,:) :: thck   
    real(dp), intent(inout), dimension(:,:,:) :: temp

    real(dp), intent(inout), dimension(:,:), optional ::  &
       lsrf,       & ! lower ice surface
       usrf,       & ! upper ice surface
       topg          ! basal topography

    real(dp), intent(inout), dimension(:,:,:,:), optional :: tracers

    integer :: nt   ! tracer index

    call parallel_halo(thck)
!    call horiz_bcs_unstag_scalar(thck)

    call parallel_halo(temp)
!!   call horiz_bcs_unstag_scalar(temp)

    ! optional updates for geometry variables

    if (present(lsrf)) then
       call parallel_halo(lsrf)
!       call horiz_bcs_unstag_scalar(lsrf)
    endif

    if (present(usrf)) then
       if (present(lsrf)) then   ! compute usrf = lsrf + thck; no halo call needed
          usrf(:,:) = lsrf(:,:) + thck(:,:)
       else
          call parallel_halo(usrf)
!          call horiz_bcs_unstag_scalar(usrf)
       endif
    endif

    if (present(topg)) then
       call parallel_halo(topg)
!       call horiz_bcs_unstag_scalar(topg)
    endif

    ! optional update for 3D tracers (e.g., ice age)

    if (present(tracers)) then

       do nt = 1, ntracer
          call parallel_halo(tracers(:,:,:,nt))
!          call horiz_bcs_unstag_scalar(tracers(:,:,:,nt))
       enddo

    endif

  end subroutine parallel_halo_scalars

!=======================================================================

    subroutine glissade_test_parallel(model)

    use parallel
    use glissade_transport, only: glissade_transport_driver
    use glimmer_paramets, only: len0
    use glimmer_physcon, only: pi

    ! various tests of parallel model

    type(glide_global_type), intent(inout) :: model      ! model instance

    integer, dimension (:,:), allocatable ::  pgID    ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable ::  pgIDr    ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDr3    ! unique global ID for parallel runs  

    integer, dimension (:,:), allocatable ::  pgIDstagi    ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable ::  pgIDstagr    ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDstagr3    ! unique global ID for parallel runs  

    real(dp), dimension(:,:,:), allocatable :: uvel, vvel   ! uniform velocity field

    logical, dimension(:,:), allocatable :: logvar
 
    integer :: i, j, k, n
    integer :: nx, ny, nz

    integer, parameter :: rdiag = 0    ! rank for diagnostic prints 

    real(dp), parameter :: dt = 1.0       ! time step in yr
    integer, parameter  :: ntstep = 10     ! run for this number of timesteps

!    real(dp), parameter :: umag = 100.    ! uniform speed (m/yr)
    real(dp), parameter :: umag = 1000.   ! uniform speed (m/yr)

    !WHL - Tested all three of these angles (eastward, northward, and northeastward)
!    real(dp), parameter :: theta = 0.d0     ! eastward
    real(dp), parameter :: theta = pi/4.d0   ! northeastward
!    real(dp), parameter :: theta = pi/2.d0  ! northward


    real(dp) :: global_row, global_col, global_ID

    print*, ' '
    print*, 'In test_parallel, this_rank =', this_rank

    nx = model%general%ewn
    ny = model%general%nsn
    nz = model%general%upn

    allocate(logvar(nx,ny))
    allocate(pgID(nx,ny))
    allocate(pgIDr(nx,ny))
    allocate(pgIDr3(nz,nx,ny))
    allocate(pgIDstagi(nx-1,ny-1))
    allocate(pgIDstagr(nx-1,ny-1))
    allocate(pgIDstagr3(nz,nx-1,ny-1))
    allocate(uvel(nz,nx-1,ny-1), vvel(nz,nx-1,ny-1))

    if (main_task) then
       print*, ' '
       print*, 'nx, ny, nz =', nx, ny, nz
       print*, 'uhalo, lhalo =', uhalo, lhalo
       print*, 'global_ewn, global_nsn =', global_ewn, global_nsn
       print*, ' '
    endif

    print*, 'this_rank, global_row/col offset =', this_rank, global_row_offset, global_col_offset

    ! logical 2D field

    logvar(:,:) = .false.

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       logvar(i,j) = .true.
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Logical field, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,*) logvar(1:34,j)
       enddo
    endif

    call parallel_halo(logvar)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,*) logvar(:,j)
       enddo
    endif

    ! integer 2D field

    ! Compute parallel global ID for each grid cell

    pgID(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgID(i,j) = parallel_globalID_scalar(i,j,nz)    ! function in parallel_mpi.F90
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (integer), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34i5)') pgID(:,j)
       enddo
    endif

    call parallel_halo(pgID)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34i5)') pgID(:,j)
       enddo
    endif

    ! real 2D
    
    pgIDr(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDr(i,j) = real(parallel_globalID_scalar(i,j,nz), dp)
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (real 2D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr(:,j)
       enddo
    endif

    call parallel_halo(pgIDr)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr(:,j)
       enddo
    endif

    ! real 3D

    pgIDr3(:,:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          pgIDr3(k,i,j) = real(parallel_globalID_scalar(i,j,nz),dp) + real(k,dp)    ! function in parallel_mpi.F90
       enddo
    enddo
    enddo

    k = 1

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (real 3D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr3(k,:,j)
       enddo
    endif

    call parallel_halo(pgIDr3)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr3(k,:,j)
       enddo
    endif

    ! Repeat for staggered variables

    ! First for an integer 2D field

    pgIDstagi(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDstagi(i,j) = parallel_globalID_scalar(i,j,nz)    ! function in parallel_mpi.F90
    enddo
    enddo

    ! Print
    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Staggered parallel global ID (integer), this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33i5)') pgIDstagi(:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagi)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33i5)') pgIDstagi(:,j)
       enddo
    endif

    ! Then for a real 2D field

    pgIDstagr(:,:) = 0.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDstagr(i,j) = real(parallel_globalID_scalar(i,j,nz),dp)    ! function in parallel_mpi.F90
    enddo
    enddo

    ! Print
    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Staggered parallel global ID (real 2D), this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr(:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagr)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr(:,j)
       enddo
    endif

    ! Then for a real 3D field

    pgIDstagr3(:,:,:) = 0.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          pgIDstagr3(k,i,j) = real(parallel_globalID_scalar(i,j,nz),dp) + real(k,dp)    ! function in parallel_mpi.F90
       enddo
    enddo
    enddo

    k = 1

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Staggered parallel global ID (real 3D), k, this_rank =', k, this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr3(k,:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagr3)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr3(k,:,j)
       enddo
    endif

    ! Run remapping routine

    uvel(:,:,:) = 0.d0
    vvel(:,:,:) = 0.d0

    ! Set velocity in locally owned cells

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          uvel(k,i,j) = umag * cos(theta)
          vvel(k,i,j) = umag * sin(theta)
       enddo
    enddo
    enddo

       if (this_rank == rdiag) then
          write(6,*) ' '
          write(6,*) 'Before halo update:'
          write(6,*) 'uvel, this_rank =', this_rank
          do j = ny-1, 1, -1
             write(6,'(33f6.1)') uvel(1,:,j)
          enddo
          write(6,*) ' '
          write(6,*) 'vvel, this_rank =', this_rank
          do j = ny-1, 1, -1
             write(6,'(33f6.1)') vvel(1,:,j)
          enddo
       endif

    ! staggered halo update

    call staggered_parallel_halo(uvel)
    call staggered_parallel_halo(vvel)

       if (this_rank == rdiag) then
          write(6,*) ' '
          write(6,*) 'After halo update:'
          write(6,*) 'uvel, this_rank =', this_rank
          do j = ny-1, 1, -1
             write(6,'(33f6.1)') uvel(1,:,j)
          enddo
          write(6,*) ' '
          write(6,*) 'vvel, this_rank =', this_rank
          do j = ny-1, 1, -1
             write(6,'(33f6.1)') vvel(1,:,j)
          enddo
       endif

    do n = 1, ntstep

       call glissade_transport_driver(dt,                                                   &
                                      model%numerics%dew * len0, model%numerics%dns * len0, &
                                      model%general%ewn,         model%general%nsn,         &
                                      model%general%upn-1,       model%numerics%sigma,      &
                                      nhalo,                     ntracer,                   &
                                      uvel(:,:,:),               vvel(:,:,:),               &
                                      model%geometry%thck(:,:),                             &
                                      dble(model%climate%acab(:,:)),                        &
                                      model%temper%bmlt(:,:),                               &
                                      model%temper%temp(:,:,:) )

       if (this_rank == rdiag) then
          write(6,*) ' '
          write(6,*) 'New thck, n =', n
          do j = ny, 1, -1
             write(6,'(19e10.2)') model%geometry%thck(1:19,j)
          enddo
          write(6,*) ' '
          write(6,*) 'New layer 1 temp, n =', n
          do j = ny, 1, -1
             write(6,'(19f10.2)') model%temper%temp(1,1:19,j)
          enddo
       endif

    enddo  ! ntstep

    if (main_task) print*, 'Done in parallel diagnostic test'

    deallocate(logvar)
    deallocate(pgID)
    deallocate(pgIDr)
    deallocate(pgIDr3)
    deallocate(pgIDstagi)
    deallocate(pgIDstagr)
    deallocate(pgIDstagr3)
    deallocate(uvel)
    deallocate(vvel)

    end subroutine glissade_test_parallel

!=======================================================================


end module glissade
