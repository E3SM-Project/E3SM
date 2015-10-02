#ifdef CPRIBM
@PROCESS ALIAS_SIZE(107374182)
#endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_timestep.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#include "glide_mask.inc"

module glint_timestep
  !*FD timestep of a GLINT instance

  use glint_type
  use glint_constants
  use glimmer_global, only: dp
  implicit none

  private
  public glint_i_tstep, glint_i_tstep_gcm

contains

  subroutine glint_i_tstep(time,            instance,       &
                           g_orog_out,      g_albedo,       &
                           g_ice_frac,      g_veg_frac,     &
                           g_snowice_frac,  g_snowveg_frac, &
                           g_snow_depth,                    &
                           g_water_in,      g_water_out,    &
                           t_win,           t_wout,         &
                           ice_vol,         out_f,          &
                           ice_tstep)

    !*FD Performs time-step of an ice model instance. 
    !*FD Note that input quantities here are accumulated/average totals since the last call.
    !*FD Global output arrays are only valid on the main task.
    !
    
    use glimmer_paramets
    use glimmer_physcon, only: rhow,rhoi
    use glimmer_log
    use glimmer_coordinates, only: coordsystem_allocate
    use glide
    use glissade
    use glide_io
    use glint_upscale, only: glint_upscaling
    use glint_mbal_coupling
    use glint_io
    use glint_mbal_io
    use glint_routing
    use glide_diagnostics
    use parallel, only: tasks, main_task

    implicit none

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    integer,                intent(in)   :: time         !*FD Current time in hours
    type(glint_instance), intent(inout)  :: instance     !*FD Model instance
    real(dp),dimension(:,:),intent(out)  :: g_orog_out   !*FD Output orography (m)
    real(dp),dimension(:,:),intent(out)  :: g_albedo     !*FD Output surface albedo 
    real(dp),dimension(:,:),intent(out)  :: g_ice_frac   !*FD Output ice fraction
    real(dp),dimension(:,:),intent(out)  :: g_veg_frac   !*FD Output veg fraction
    real(dp),dimension(:,:),intent(out)  :: g_snowice_frac !*FD Output snow-ice fraction
    real(dp),dimension(:,:),intent(out)  :: g_snowveg_frac !*FD Output snow-veg fraction
    real(dp),dimension(:,:),intent(out)  :: g_snow_depth !*FD Output snow depth (m)
    real(dp),dimension(:,:),intent(out)  :: g_water_in   !*FD Input water flux (m)
    real(dp),dimension(:,:),intent(out)  :: g_water_out  !*FD Output water flux (m)
    real(dp),               intent(out)  :: t_win        !*FD Total water input (kg)
    real(dp),               intent(out)  :: t_wout       !*FD Total water output (kg)
    real(dp),               intent(out)  :: ice_vol      !*FD Output ice volume (m$^3$)
    type(output_flags),     intent(in)   :: out_f        !*FD Flags to tell us whether to do output   
    logical,                intent(out)  :: ice_tstep    !*FD Set if we have done an ice time step

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    real(dp),dimension(:,:),pointer :: upscale_temp => null() ! temporary array for upscaling
    real(dp),dimension(:,:),pointer :: routing_temp => null() ! temporary array for flow routing
    real(dp),dimension(:,:),pointer :: accum_temp   => null() ! temporary array for accumulation
    real(dp),dimension(:,:),pointer :: ablat_temp   => null() ! temporary array for ablation
    integer, dimension(:,:),pointer :: fudge_mask   => null() ! temporary array for fudging
    real(dp),dimension(:,:),pointer :: thck_temp    => null() ! temporary array for volume calcs
    real(dp),dimension(:,:),pointer :: calve_temp   => null() ! temporary array for calving flux
    real(dp) :: start_volume,end_volume,flux_fudge            ! note: only valid for single-task runs
    integer :: i, j, k, ii, jj, nx, ny, il, jl, ig, jg

    logical :: gcm_smb   ! true if getting sfc mass balance from a GCM

    !WHL - Moved glint_downscaling call up to glint_main

    ice_tstep = .false.

    ! Assume we always need this, as it's too complicated to work out when we do and don't

    call coordsystem_allocate(instance%lgrid, thck_temp)
    call coordsystem_allocate(instance%lgrid, calve_temp)

    ! ------------------------------------------------------------------------  
    ! Sort out some local orography and remove bathymetry. This relies on the 
    ! point 1,1 being underwater. However, it's a better method than just 
    ! setting all points < 0.0 to zero
    ! ------------------------------------------------------------------------  

    call glide_get_usurf(instance%model, instance%local_orog)
    call glint_remove_bath(instance%local_orog,1,1)

    ! ------------------------------------------------------------------------  
    ! Adjust the surface temperatures using the lapse-rate, by reducing to
    ! sea-level and then back up to high-res orography
    ! ------------------------------------------------------------------------  

    call glint_lapserate(instance%artm, real(instance%global_orog,dp), real(-instance%data_lapse_rate,dp))
    call glint_lapserate(instance%artm, real(instance%local_orog, dp), real( instance%lapse_rate,dp))

    ! Process the precipitation field if necessary ---------------------------
    ! and convert from mm/s to m/s

    call glint_calc_precip(instance)

    ! Get ice thickness ----------------------------------------

    call glide_get_thk(instance%model,thck_temp)

    ! Accumulate mass balance fields -----------------------------------------

    call glint_accumulate_mbal(instance%mbal_accum, time, instance%artm, instance%arng, instance%prcp, &
                               instance%snowd, instance%siced, instance%xwind, instance%ywind, &
                               instance%local_orog, real(thck_temp,dp), instance%humid,    &
                               instance%swdown, instance%lwdown, instance%airpress)

    ! Initialise water budget quantities to zero. These will be over-ridden if
    ! there's an ice-model time-step

    t_win = 0.d0       ; t_wout = 0.d0
    g_water_out = 0.d0 ; g_water_in = 0.d0

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'In glint_i_tstep, time =', time
       write(stdout,*) 'next_time =', instance%next_time
       write(stdout,*) 'Check for ice dynamics timestep'
       write(stdout,*) 'start_time =', instance%mbal_accum%start_time
       write(stdout,*) 'mbal_step =', instance%mbal_tstep
       write(stdout,*) 'mbal_accum_time =', instance%mbal_accum_time
       write(stdout,*) 'time-start_time+mbal_tstep =', time - instance%mbal_accum%start_time + instance%mbal_tstep
       write(stdout,*) 'ice_tstep =', instance%ice_tstep
       write(stdout,*) 'n_icetstep =', instance%n_icetstep
    end if

    ! ------------------------------------------------------------------------  
    ! ICE TIMESTEP begins HERE ***********************************************
    ! ------------------------------------------------------------------------  

    if (time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time) then

       if (instance%mbal_accum_time < instance%ice_tstep) then 
          instance%next_time = instance%next_time + instance%ice_tstep - instance%mbal_tstep

          !WHL - debug
          print*, 'Adjust next_time:', instance%next_time

       end if

       ice_tstep = .true.

       ! Prepare arrays for water budgeting

       if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then
          call coordsystem_allocate(instance%lgrid, accum_temp)
          call coordsystem_allocate(instance%lgrid, ablat_temp)
          accum_temp = 0.d0
          ablat_temp = 0.d0
       end if

       ! Calculate the initial ice volume (scaled and converted to water equivalent)
       ! start_volume is only valid for single-task runs (this is checked in the place
       ! where it is used)

       call glide_get_thk(instance%model, thck_temp)
       thck_temp = thck_temp * rhoi/rhow
       start_volume = sum(thck_temp)

       ! ---------------------------------------------------------------------
       ! Timestepping for the dynamic ice sheet model
       ! ---------------------------------------------------------------------

       do i = 1, instance%n_icetstep

          if (GLC_DEBUG .and. main_task) then
             write (stdout,*) 'Ice sheet timestep, iteration =', i
          end if

          ! Calculate the initial ice volume (scaled and converted to water equivalent)
          call glide_get_thk(instance%model,thck_temp)
          thck_temp = thck_temp * rhoi/rhow

          ! Get latest upper-surface elevation (needed for masking)
          call glide_get_usurf(instance%model, instance%local_orog)

          call glint_remove_bath(instance%local_orog,1,1)

          ! Get the average mass-balance, as m water/year 
          call glint_average_mbal(instance%mbal_accum,    &
                                  instance%artm,  instance%prcp,   &
                                  instance%ablt,  instance%acab,   &
                                  instance%snowd, instance%siced,  &
                                  instance%mbal_accum_time)

          ! Mask out non-accumulation in ice-free areas

          where(thck_temp <= 0.d0 .and. instance%acab < 0.d0)
             instance%acab = 0.d0
             instance%ablt = instance%prcp
          end where

          ! Set acab to zero for ocean cells (bed below sea level, no ice present)

          where (GLIDE_IS_OCEAN(instance%model%geometry%thkmask))
             instance%acab = 0.d0
          endwhere

          ! Put climate inputs in the appropriate places, with conversion ----------

          ! Note on units: 
          ! For this subroutine, input acab is in m/yr; this value is divided 
          !  by scale_acab = scyr*thk0/tim0 and copied to data%climate%acab.
          ! Input artm is in deg C; this value is copied to data%climate%artm (no unit conversion).

          call glide_set_acab(instance%model, instance%acab * rhow/rhoi)
          call glide_set_artm(instance%model, instance%artm)

          ! This will work only for single-processor runs
          if (GLC_DEBUG .and. tasks==1) then
             il = instance%model%numerics%idiag
             jl = instance%model%numerics%jdiag
             write (stdout,*) ' '
             write (stdout,*) 'After glide_set_acab, glide_set_artm: i, j =', il, jl
             write (stdout,*) 'acab (m/y), artm (C) =', instance%acab(il,jl)*rhow/rhoi, instance%artm(il,jl)
          end if

          ! Adjust glint acab and ablt for output
 
          where (instance%acab < -thck_temp .and. thck_temp > 0.d0)
             instance%acab = -thck_temp
             instance%ablt =  thck_temp
          end where

          instance%glide_time = instance%glide_time + instance%model%numerics%tinc

          ! call the dynamic ice sheet model (provided the ice is allowed to evolve)

          if (instance%evolve_ice == EVOLVE_ICE_TRUE) then

             if (instance%model%options%whichdycore == DYCORE_GLIDE) then

                call glide_tstep_p1(instance%model,instance%glide_time)

                call glide_tstep_p2(instance%model)

                call glide_tstep_p3(instance%model)

             else   ! glam/glissade dycore

                call glissade_tstep(instance%model,instance%glide_time)

             endif

          endif   ! evolve_ice

          ! Add the calved ice to the ablation field

          call glide_get_calving(instance%model, calve_temp)
          calve_temp = calve_temp * rhoi/rhow

          instance%ablt = instance%ablt + calve_temp/instance%model%numerics%tinc
          instance%acab = instance%acab - calve_temp/instance%model%numerics%tinc

          ! Accumulate for water-budgeting
          if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then
             accum_temp = accum_temp + instance%prcp*instance%model%numerics%tinc
             ablat_temp = ablat_temp + instance%ablt*instance%model%numerics%tinc
          endif

          ! write ice sheet diagnostics at specified interval

          call glide_write_diagnostics(instance%model,                  &
                                       instance%model%numerics%time,    &
                                       tstep_count = instance%model%numerics%timecounter)

          ! write netCDf output

          call glide_io_writeall(instance%model,instance%model)
          call glint_io_writeall(instance,instance%model)

       end do   ! n_icestep

       ! Calculate flux fudge factor --------------------------------------------

       if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then

          ! WJS (1-15-13): I am pretty sure (but not positive) that the stuff in this
          ! conditional will only work right with a single task
          if (tasks > 1) then
             call write_log('The sums in the computation of a flux fudge factor only work with a single task', &
                            GM_FATAL, __FILE__, __LINE__)
          end if
          
          call coordsystem_allocate(instance%lgrid,fudge_mask)

          call glide_get_thk(instance%model,thck_temp)
          end_volume = sum(thck_temp)

          where (thck_temp > 0.d0)
             fudge_mask = 1
          elsewhere
             fudge_mask = 0
          endwhere

          flux_fudge = (start_volume + sum(accum_temp) - sum(ablat_temp) - end_volume) / sum(fudge_mask)

          ! Apply fudge_factor

          where(thck_temp > 0.d0)
             ablat_temp = ablat_temp + flux_fudge
          endwhere
          
          deallocate(fudge_mask)
          fudge_mask => null()

       endif

       ! Upscale water flux fields ----------------------------------------------
       ! First water input (i.e. mass balance + ablation)

       if (out_f%water_in) then
       
          call coordsystem_allocate(instance%lgrid, upscale_temp)

          where (thck_temp > 0.d0)
             upscale_temp = accum_temp
          elsewhere
             upscale_temp = 0.d0
          endwhere

          call local_to_global_avg(instance%ups,   &
                                   upscale_temp,   &
                                   g_water_in,     &
                                   instance%out_mask)
          deallocate(upscale_temp)
          upscale_temp => null()
       endif

       ! Now water output (i.e. ablation) - and do routing

       if (out_f%water_out) then    
       
          ! WJS (1-15-13): The flow_router routine (called bolew) currently seems to
          ! assume that it's working on the full (non-decomposed) domain. I'm not sure
          ! what the best way is to fix this, so for now we only allow this code to be
          ! executed if tasks==1.
          if (tasks > 1) then
             call write_log('water_out computation assumes a single task', &
                            GM_FATAL, __FILE__, __LINE__)
          end if

          call coordsystem_allocate(instance%lgrid, upscale_temp)
          call coordsystem_allocate(instance%lgrid, routing_temp)

          where (thck_temp > 0.d0)
             upscale_temp = ablat_temp
          elsewhere
             upscale_temp = 0.d0
          endwhere

          call glide_get_usurf(instance%model, instance%local_orog)
          call flow_router(instance%local_orog, &
                           upscale_temp,        &
                           routing_temp,        &
                           instance%out_mask,   &
                           real(instance%lgrid%delta%pt(1),dp), &
                           real(instance%lgrid%delta%pt(2),dp))

          call local_to_global_avg(instance%ups,   &
                                   routing_temp,   &
                                   g_water_out,    &
                                   instance%out_mask)

          deallocate(upscale_temp,routing_temp)
          upscale_temp => null()
          routing_temp => null()

       endif

       ! Sum water fluxes and convert if necessary ------------------------------

       if (out_f%total_win) then
          if (tasks > 1) call write_log('t_win sum assumes a single task', &
                                        GM_FATAL, __FILE__, __LINE__)

          t_win  = sum(accum_temp) * instance%lgrid%delta%pt(1)* &
                                     instance%lgrid%delta%pt(2)
       endif

       if (out_f%total_wout) then
          if (tasks > 1) call write_log('t_wout sum assumes a single task', &
                                        GM_FATAL, __FILE__, __LINE__)

          t_wout = sum(ablat_temp) * instance%lgrid%delta%pt(1)* &
                                     instance%lgrid%delta%pt(2)
       endif

    end if  ! time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time

    ! Output instantaneous values

    call glint_mbal_io_writeall(instance, instance%model,       &
                                outfiles = instance%out_first,  &
                                time = time*hours2years)

    ! ------------------------------------------------------------------------ 
    ! Upscaling of output
    ! ------------------------------------------------------------------------ 

    ! We now upscale all fields at once...

    call glint_upscaling(instance, g_orog_out, g_albedo, g_ice_frac, g_veg_frac, &
                         g_snowice_frac, g_snowveg_frac, g_snow_depth)

    ! Calculate ice volume ---------------------------------------------------

    if (out_f%ice_vol) then
       if (tasks > 1) call write_log('ice_vol sum assumes a single task', &
                                     GM_FATAL, __FILE__, __LINE__)

       call glide_get_thk(instance%model, thck_temp)
       ice_vol = sum(thck_temp) * instance%lgrid%delta%pt(1)* &
                                  instance%lgrid%delta%pt(2)
    endif

    ! Tidy up ----------------------------------------------------------------

    if (associated(accum_temp)) then 
       deallocate(accum_temp)
       accum_temp => null()
    end if

    if (associated(ablat_temp)) then
       deallocate(ablat_temp)
       ablat_temp => null()
    end if

    if (associated(calve_temp)) then
       deallocate(calve_temp)
       calve_temp => null()
    end if

    if (associated(thck_temp)) then
       deallocate(thck_temp)
       thck_temp => null()
    endif

  end subroutine glint_i_tstep

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_tstep_gcm(time,            instance,       &
                               ice_tstep)

    ! Performs time-step of an ice model instance. 
    ! Input quantities here are accumulated/average totals since the last call.
    ! Global output arrays are only valid on the main task.
    !
    use glimmer_paramets
    use glimmer_physcon, only: rhow, rhoi
    use glimmer_log
    use glimmer_coordinates, only: coordsystem_allocate
    use glide
    use glissade
    use glide_io
    use glint_downscale, only: glint_accumulate_input_gcm, glint_average_input_gcm
    use glint_upscale, only: glint_accumulate_output_gcm
    use glint_io
    use glint_mbal_io
    use glide_diagnostics
    use parallel, only: tasks, main_task

    implicit none

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    integer,                intent(in)   :: time         ! Current time in hours
    type(glint_instance), intent(inout)  :: instance     ! Model instance
    logical,                intent(out)  :: ice_tstep    ! Set if we have done an ice time step

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    !TODO - Are these needed?
    real(dp),dimension(:,:),pointer :: thck_temp    => null() ! temporary array for volume calcs

    integer :: i, j, k, nx, ny, il, jl, ig, jg

    !WHL - Moved glint_downscaling call up to glint_main

    ice_tstep = .false.

    call coordsystem_allocate(instance%lgrid, thck_temp)

    ! ------------------------------------------------------------------------
    ! Sort out some local orography and remove bathymetry. This relies on the 
    ! point 1,1 being underwater. However, it's a better method than just 
    ! setting all points < 0.0 to zero
    ! ------------------------------------------------------------------------  

!TODO: Determine if glint_remove_bath is needed in a CESM run. If so, fix it to work with
!      multiple tasks. 

!!    call glide_get_usurf(instance%model, instance%local_orog)
!!    call glint_remove_bath(instance%local_orog,1,1)

    ! Get ice thickness ----------------------------------------

    call glide_get_thk(instance%model,thck_temp)

    ! Accumulate Glide input fields, acab and artm
    ! Note: At this point, instance%acab has units of m
    !       Upon averaging (in glint_mbal_gcm), units are converted to m/yr

    call glint_accumulate_input_gcm(instance%mbal_accum,   time,        &
                                    instance%acab,         instance%artm)


    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'In glint_i_tstep_gcm, time =', time
       write(stdout,*) 'next_time =', instance%next_time
       write(stdout,*) 'Check for ice dynamics timestep'
       write(stdout,*) 'time =', time
       write(stdout,*) 'start_time =', instance%mbal_accum%start_time
       write(stdout,*) 'mbal_step =', instance%mbal_tstep
       write(stdout,*) 'mbal_accum_time =', instance%mbal_accum_time
       write(stdout,*) 'time-start_time+mbal_tstep =', time - instance%mbal_accum%start_time + instance%mbal_tstep
       write(stdout,*) 'ice_tstep =', instance%ice_tstep
       write(stdout,*) 'n_icetstep =', instance%n_icetstep
    end if

    ! ------------------------------------------------------------------------  
    ! ICE TIMESTEP begins HERE ***********************************************
    ! ------------------------------------------------------------------------  

    if (time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time) then

       if (instance%mbal_accum_time < instance%ice_tstep) then 
          instance%next_time = instance%next_time + instance%ice_tstep - instance%mbal_tstep
       end if

       ice_tstep = .true.

       ! ---------------------------------------------------------------------
       ! Timestepping for ice sheet model
       ! ---------------------------------------------------------------------

       do i = 1, instance%n_icetstep

          if (GLC_DEBUG .and. main_task) then
             write (stdout,*) 'Ice sheet timestep, iteration =', i
          end if

          ! Get average values of acab and artm during mbal_accum_time
          ! instance%acab has units of m/yr w.e. after averaging

          call glint_average_input_gcm(instance%mbal_accum, instance%mbal_accum_time,  &
                                       instance%acab,       instance%artm)
                                  
          ! Calculate the initial ice volume (scaled and converted to water equivalent)
          call glide_get_thk(instance%model,thck_temp)
          thck_temp = thck_temp * rhoi/rhow

          !TODO: Determine if glint_remove_bath is needed in a CESM run. If so, fix it to work with
          !      multiple tasks.  (And decide whether the call is needed both here and above) 
          ! Get latest upper-surface elevation (needed for masking)
!!          call glide_get_usurf(instance%model, instance%local_orog)
!!          call glint_remove_bath(instance%local_orog,1,1)

          ! Mask out non-accumulation in ice-free areas

          where(thck_temp <= 0.d0 .and. instance%acab < 0.d0)
             instance%acab = 0.d0
          end where

          ! Set acab to zero for ocean cells (bed below sea level, no ice present)

          where (GLIDE_IS_OCEAN(instance%model%geometry%thkmask))
             instance%acab = 0.d0
          endwhere

          ! Put climate inputs in the appropriate places, with conversion ----------

          ! Note on units:
          ! For subroutine glide_set_acab, input acab is in m/yr ice; this value is multiplied 
          !  by tim0/(scyr*thk0) and copied to data%climate%acab.
          ! Input artm is in deg C; this value is copied to data%climate%artm (no unit conversion).

          !TODO - It is confusing to have units of m/yr w.e. for instance%acab, compared to units
          !       of m/yr ice for Glide. It would be better to have the same units consistently.
          !       E.g., switch to w.e. in Glide

          call glide_set_acab(instance%model, instance%acab * rhow/rhoi)
          call glide_set_artm(instance%model, instance%artm)

          ! This will work only for single-processor runs
          if (GLC_DEBUG .and. tasks==1) then
             il = instance%model%numerics%idiag
             jl = instance%model%numerics%jdiag
             write (stdout,*) ' '
             write (stdout,*) 'After glide_set_acab, glide_set_artm: i, j =', il, jl
             write (stdout,*) 'acab (m/y), artm (C) =', instance%acab(il,jl)*rhow/rhoi, instance%artm(il,jl)
          end if

          ! Adjust glint acab and ablt for output
 
          where (instance%acab < -thck_temp .and. thck_temp > 0.d0)
             instance%acab = -thck_temp
          end where

          instance%glide_time = instance%glide_time + instance%model%numerics%tinc

          ! call the dynamic ice sheet model (provided the ice is allowed to evolve)

          if (instance%evolve_ice == EVOLVE_ICE_TRUE) then

             if (instance%model%options%whichdycore == DYCORE_GLIDE) then
 
                call glide_tstep_p1(instance%model, instance%glide_time)

                call glide_tstep_p2(instance%model)
 
                call glide_tstep_p3(instance%model)

             else   ! glam/glissade dycore

                call glissade_tstep(instance%model, instance%glide_time)

             endif

          endif  ! evolve_ice

          ! write ice sheet diagnostics at specified interval (model%numerics%dt_diag)

          call glide_write_diagnostics(instance%model,                  &
                                       instance%model%numerics%time,    &
                                       tstep_count = instance%model%numerics%timecounter)

          ! write netCDF output

          call glide_io_writeall(instance%model,instance%model)
          call glint_io_writeall(instance,instance%model)

          ! Accumulate Glide output fields to be sent to GCM

          call glint_accumulate_output_gcm(instance%model,            &
                                           instance%av_count_output,  &
                                           instance%new_tavg_output,  &
                                           instance%rofi_tavg,        &
                                           instance%rofl_tavg,        &
                                           instance%hflx_tavg )

       end do   ! instance%n_icetstep

    end if   ! time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time

    ! Output instantaneous values

    call glint_mbal_io_writeall(instance, instance%model,       &
                                outfiles = instance%out_first,  &
                                time = time*hours2years)

    ! Deallocate

    if (associated(thck_temp)) then
       deallocate(thck_temp)
       thck_temp => null()
    endif

  end subroutine glint_i_tstep_gcm

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Rewrite to support multiple tasks?

  subroutine glint_remove_bath(orog,x,y)

    ! Sets ocean areas to zero height, working recursively from
    ! a known ocean point.

    use glimmer_log
    use parallel, only : tasks

    real(dp),dimension(:,:),intent(inout) :: orog !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y  !*FD Location of starting point (index)

    integer :: nx,ny

    ! Currently, this routine is called assuming point 1,1 is ocean... this won't be true
    ! when running on multiple processors, with a distributed grid
    ! This can't be made a fatal error, because this is currently called even if we have
    ! more than one task... the hope is just that the returned data aren't needed in CESM.
    if (tasks > 1) then
       call write_log('Use of glint_remove_bath currently assumes the use of only one task', &
                      GM_WARNING, __FILE__, __LINE__)
    end if

    nx=size(orog,1) ; ny=size(orog,2)

    if (orog(x,y) < 0.d0) orog(x,y) = 0.d0
    call glint_find_bath(orog,x,y,nx,ny)

  end subroutine glint_remove_bath

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  recursive subroutine glint_find_bath(orog,x,y,nx,ny)

    !*FD Recursive subroutine called by {\tt glimmer\_remove\_bath}.

    real(dp),dimension(:,:),intent(inout) :: orog  !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y   !*FD Starting point
    integer,                intent(in)    :: nx,ny !*FD Size of array {\tt orography}

    integer,dimension(4) :: xi = (/ -1,1,0,0 /)
    integer,dimension(4) :: yi = (/ 0,0,-1,1 /)
    integer :: ns = 4
    integer :: i

    do i=1,ns
       if (x+xi(i) <= nx .and. x+xi(i) > 0 .and. &
           y+yi(i) <= ny .and. y+yi(i) > 0) then
          if (orog(x+xi(i),y+yi(i)) < 0.d0) then
             orog(x+xi(i),y+yi(i)) = 0.d0
             call glint_find_bath(orog,x+xi(i),y+yi(i),nx,ny)
          endif
       endif
    enddo

  end subroutine glint_find_bath

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!WHL - This used to be called glint_lapserate_dp

  subroutine glint_lapserate(temp,topo,lr)

    !*FD Corrects the temperature field
    !*FD for height, using a constant lapse rate.
    !*FD
    !*FD This is the double-precision version, aliased as \texttt{glimmer\_lapserate}.

    implicit none

    real(dp),dimension(:,:), intent(inout) :: temp !*FD temperature at sea-level in $^{\circ}$C
                                                   !*FD used for input and output
    real(dp),dimension(:,:), intent(in)    :: topo !*FD topography field (m above msl)
    real(dp),                intent(in)    :: lr   !*FD Lapse rate ($^{\circ}\mathrm{C\,km}^{-1}$).
                                                   !*FD
                                                   !*FD NB: the lapse rate is positive for 
                                                   !*FD falling temp with height\ldots

    temp = temp-(lr*topo/1000.d0)                  ! The lapse rate calculation.

  end subroutine glint_lapserate

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!WHL - Removed subroutine glint_lapserate_sp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_calc_precip(instance)

    use glint_precip_param
    use glimmer_log

    !*FD Process precip if necessary

    type(glint_instance) :: instance

    select case (instance%whichprecip)

    case(PRECIP_STANDARD)
       ! Do nothing to the precip field

    case(PRECIP_RL)
       ! Use the Roe/Lindzen parameterisation
       call glint_precip(instance%prcp, &
                         instance%xwind, &
                         instance%ywind, &
                         instance%artm, &
                         instance%local_orog, &
                         real(instance%lgrid%delta%pt(1),dp), &
                         real(instance%lgrid%delta%pt(2),dp), &
                         fixed_a=.true.)

    case default

       call write_log('Invalid value of whichprecip',GM_FATAL,__FILE__,__LINE__)

    end select

    ! Convert from mm/s to m/s - very important!

    instance%prcp = instance%prcp * 1.d-3

  end subroutine glint_calc_precip

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glint_timestep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
