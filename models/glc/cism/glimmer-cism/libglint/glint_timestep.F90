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

  implicit none

  private
  public glint_i_tstep, glint_i_tstep_gcm, &
         get_i_upscaled_fields, get_i_upscaled_fields_gcm

  interface glint_lapserate
     module procedure glint_lapserate_dp, glint_lapserate_sp
  end interface

contains

  subroutine glint_i_tstep(time,            instance,       &
                           g_temp,          g_temp_range,   &
                           g_precip,        g_zonwind,      &
                           g_merwind,       g_humid,        &
                           g_lwdown,        g_swdown,       &
                           g_airpress,      g_orog,         &
                           g_orog_out,      g_albedo,       &
                           g_ice_frac,      g_veg_frac,     &
                           g_snowice_frac,  g_snowveg_frac, &
                           g_snow_depth,                    &
                           g_water_in,      g_water_out,    &
                           t_win,           t_wout,         &
                           ice_vol,         out_f,          &
                           orogflag,        ice_tstep)

    !*FD Performs time-step of an ice model instance. 
    !*FD Note that input quantities here are accumulated/average totals since the
    !*FD last call.
    !*FD Global output arrays are only valid on the main task.
    !
    
    use glimmer_paramets
    use glimmer_physcon, only: rhow,rhoi
    use glimmer_log
    use glide
    use glissade
    use glide_io
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
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)

    real(rk),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_humid      !*FD Global surface humidity (%)
    real(rk),dimension(:,:),intent(in)   :: g_lwdown     !*FD Global downwelling longwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_swdown     !*FD Global downwelling shortwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_airpress   !*FD Global surface air pressure (Pa)

    real(rk),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_orog_out   !*FD Output orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_albedo     !*FD Output surface albedo 
    real(rk),dimension(:,:),intent(out)  :: g_ice_frac   !*FD Output ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_veg_frac   !*FD Output veg fraction
    real(rk),dimension(:,:),intent(out)  :: g_snowice_frac !*FD Output snow-ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_snowveg_frac !*FD Output snow-veg fraction
    real(rk),dimension(:,:),intent(out)  :: g_snow_depth !*FD Output snow depth (m)
    real(rk),dimension(:,:),intent(out)  :: g_water_in   !*FD Input water flux (m)
    real(rk),dimension(:,:),intent(out)  :: g_water_out  !*FD Output water flux (m)
    real(rk),               intent(out)  :: t_win        !*FD Total water input (kg)
    real(rk),               intent(out)  :: t_wout       !*FD Total water output (kg)
    real(rk),               intent(out)  :: ice_vol      !*FD Output ice volume (m$^3$)
    type(output_flags),     intent(in)   :: out_f        !*FD Flags to tell us whether to do output   
    logical,                intent(in)   :: orogflag     !*FD Set if we have new global orog
    logical,                intent(out)  :: ice_tstep    !*FD Set if we have done an ice time step

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    real(rk),dimension(:,:),pointer :: upscale_temp => null() ! temporary array for upscaling
    real(rk),dimension(:,:),pointer :: routing_temp => null() ! temporary array for flow routing
    real(rk),dimension(:,:),pointer :: accum_temp   => null() ! temporary array for accumulation
    real(rk),dimension(:,:),pointer :: ablat_temp   => null() ! temporary array for ablation
    integer, dimension(:,:),pointer :: fudge_mask   => null() ! temporary array for fudging
    real(sp),dimension(:,:),pointer :: thck_temp    => null() ! temporary array for volume calcs
    real(sp),dimension(:,:),pointer :: calve_temp   => null() ! temporary array for calving flux
    real(rk) :: start_volume,end_volume,flux_fudge            ! note: only valid for single-task runs
    integer :: i, j, k, ii, jj, nx, ny, il, jl, ig, jg

    logical :: gcm_smb   ! true if getting sfc mass balance from a GCM

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'In glint_i_tstep, time =', time
       write(stdout,*) 'next_time =', instance%next_time
    end if

    ! Check whether we're doing anything this time.

    if (time /= instance%next_time) then
       return
    else
       instance%next_time = instance%next_time + instance%mbal_tstep
    end if

    ! Assume we always need this, as it's too complicated to work out when we do and don't

    call coordsystem_allocate(instance%lgrid, thck_temp)
    call coordsystem_allocate(instance%lgrid, calve_temp)
    ice_tstep = .false.

    ! Downscale input fields from global to local grid
    ! This subroutine computes instance%acab and instance%artm, the key inputs to GLIDE.

    call glint_downscaling(instance,                  &
                           g_temp,     g_temp_range,  &
                           g_precip,   g_orog,        &
                           g_zonwind,  g_merwind,     &
                           g_humid,    g_lwdown,      &
                           g_swdown,   g_airpress,    &
                           orogflag)

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

    call glint_lapserate(instance%artm, real(instance%global_orog,rk), real(-instance%data_lapse_rate,rk))
    call glint_lapserate(instance%artm, real(instance%local_orog,rk),  real(instance%lapse_rate,rk))

    ! Process the precipitation field if necessary ---------------------------
    ! and convert from mm/s to m/s

    call glint_calc_precip(instance)

    ! Get ice thickness ----------------------------------------

    call glide_get_thk(instance%model,thck_temp)

    ! Do accumulation --------------------------------------------------------

    call glint_accumulate(instance%mbal_accum, time, instance%artm, instance%arng, instance%prcp, &
                          instance%snowd, instance%siced, instance%xwind, instance%ywind, &
                          instance%local_orog, real(thck_temp,rk), instance%humid,    &
                          instance%swdown, instance%lwdown, instance%airpress)

    ! Initialise water budget quantities to zero. These will be over-ridden if
    ! there's an ice-model time-step

    t_win=0.0       ; t_wout=0.0
    g_water_out=0.0 ; g_water_in=0.0

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'Check for ice dynamics timestep'
       write(stdout,*) 'time =', time
       write(stdout,*) 'start_time =', instance%mbal_accum%start_time
       write(stdout,*) 'mbal_step =', instance%mbal_tstep
       write(stdout,*) 'mbal_accum_time =', instance%mbal_accum_time
    end if

    ! ------------------------------------------------------------------------  
    ! ICE TIMESTEP begins HERE ***********************************************
    ! ------------------------------------------------------------------------  

    if (time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time) then

       if (instance%mbal_accum_time < instance%ice_tstep) then 
          instance%next_time = instance%next_time + instance%ice_tstep - instance%mbal_tstep
       end if

       ice_tstep = .true.

       ! Prepare arrays for water budgeting

       if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then
          call coordsystem_allocate(instance%lgrid, accum_temp)
          call coordsystem_allocate(instance%lgrid, ablat_temp)
          accum_temp = 0.0
          ablat_temp = 0.0
       end if

       ! Calculate the initial ice volume (scaled and converted to water equivalent)
       ! start_volume is only valid for single-task runs (this is checked in the place
       ! where it is used)

       call glide_get_thk(instance%model, thck_temp)
       thck_temp = thck_temp*real(rhoi/rhow)
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
          thck_temp = thck_temp*real(rhoi/rhow)

          ! Get latest upper-surface elevation (needed for masking)
          call glide_get_usurf(instance%model, instance%local_orog)

          call glint_remove_bath(instance%local_orog,1,1)

          ! Get the mass-balance, as m water/year 
          call glint_get_mbal(instance%mbal_accum, instance%artm, instance%prcp, instance%ablt, &
                              instance%acab, instance%snowd, instance%siced, instance%mbal_accum_time)

          ! Mask out non-accumulation in ice-free areas

          where(thck_temp <= 0.0 .and. instance%acab < 0.0)
             instance%acab = 0.0
             instance%ablt = instance%prcp
          end where

          ! Set acab to zero for ocean cells (bed below sea level, no ice present)

          where (GLIDE_IS_OCEAN(instance%model%geometry%thkmask))
             instance%acab = 0.0
          endwhere

          ! Put climate inputs in the appropriate places, with conversion ----------

          ! Note on units: 
          ! For this subroutine, input acab is in m/yr; this value is divided 
          !  by scale_acab = scyr*thk0/tim0 and copied to data%climate%acab.
          ! Input artm is in deg C; this value is copied to data%climate%artm (no unit conversion).

          !TODO - Change to dp
          call glide_set_acab(instance%model, instance%acab*real(rhow/rhoi))
          call glide_set_artm(instance%model, instance%artm)

          ! This will work only for single-processor runs
          if (GLC_DEBUG .and. tasks==1) then
             il = instance%model%numerics%idiag_global
             jl = instance%model%numerics%jdiag_global
             write (stdout,*) ' '
             write (stdout,*) 'After glide_set_acab, glide_set_artm: i, j =', il, jl
             write (stdout,*) 'acab (m/y), artm (C) =', instance%acab(il,jl)*rhow/rhoi, instance%artm(il,jl)
          end if

          ! Adjust glint acab and ablt for output
 
          where (instance%acab < -thck_temp .and. thck_temp > 0.0)
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
          calve_temp = calve_temp * real(rhoi/rhow)

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

          where (thck_temp > 0.0)
             fudge_mask = 1
          elsewhere
             fudge_mask = 0
          endwhere

          flux_fudge = (start_volume + sum(accum_temp) - sum(ablat_temp) - end_volume) / sum(fudge_mask)

          ! Apply fudge_factor

          where(thck_temp > 0.0)
             ablat_temp = ablat_temp + flux_fudge
          endwhere
          
          deallocate(fudge_mask)
          fudge_mask => null()

       endif

       ! Upscale water flux fields ----------------------------------------------
       ! First water input (i.e. mass balance + ablation)

       if (out_f%water_in) then
          call coordsystem_allocate(instance%lgrid, upscale_temp)

          where (thck_temp > 0.0)
             upscale_temp = accum_temp
          elsewhere
             upscale_temp = 0.0
          endwhere

          call mean_to_global(instance%ups,   &
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

          where (thck_temp > 0.0)
             upscale_temp = ablat_temp
          elsewhere
             upscale_temp = 0.0
          endwhere

          call glide_get_usurf(instance%model, instance%local_orog)
          call flow_router(instance%local_orog, &
                           upscale_temp, &
                           routing_temp, &
                           instance%out_mask, &
                           real(instance%lgrid%delta%pt(1),rk), &
                           real(instance%lgrid%delta%pt(2),rk))

          call mean_to_global(instance%ups,   &
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
    !TODO - This subroutine is called here and also from subroutine glint.
    !       Are both calls needed?

    call get_i_upscaled_fields(instance, g_orog_out, g_albedo, g_ice_frac, g_veg_frac, &
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
                               ice_tstep,                       &
                               qsmb_g,          tsfc_g,         &
                               topo_g,          gmask,          &
                               gfrac,           gtopo,          &
                               grofi,           grofl,          &
                               ghflx )

    ! Performs time-step of an ice model instance. 
    ! Input quantities here are accumulated/average totals since the last call.
    ! Global output arrays are only valid on the main task.
    !
    use glimmer_paramets
    use glimmer_physcon, only: rhow, rhoi
    use glimmer_log
    use glide
    use glissade
    use glide_io
    use glint_mbal_coupling
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

    real(rk),dimension(:,:,:),intent(in)  :: qsmb_g    ! Depth of new ice (m)
    real(rk),dimension(:,:,:),intent(in)  :: tsfc_g    ! Surface temperature (C)
    real(rk),dimension(:,:,:),intent(in)  :: topo_g    ! Surface elevation (m)

    integer, dimension(:,:),  optional,intent(in)  :: gmask     ! = 1 where global data are valid, else = 0

    real(rk),dimension(:,:,:),optional,intent(out) :: gfrac     ! ice fractional area [0,1]
    real(rk),dimension(:,:,:),optional,intent(out) :: gtopo     ! surface elevation (m)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofi     ! ice runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofl     ! liquid runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: ghflx     ! heat flux (W/m^2, positive down)

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    !TODO - Are these needed?
    real(rk),dimension(:,:),pointer :: upscale_temp => null() ! temporary array for upscaling
    real(sp),dimension(:,:),pointer :: thck_temp    => null() ! temporary array for volume calcs
    real(sp),dimension(:,:),pointer :: calve_temp   => null() ! temporary array for calving flux

    integer :: i, j, k, nx, ny, il, jl, ig, jg

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'In glint_i_tstep, time =', time
       write(stdout,*) 'next_time =', instance%next_time
    end if

    ! Zero outputs

    !TODO - Change to dp
    if (present(gfrac)) gfrac(:,:,:) = 0._rk
    if (present(gtopo)) gtopo(:,:,:) = 0._rk
    if (present(grofi)) grofi(:,:,:) = 0._rk
    if (present(grofl)) grofl(:,:,:) = 0._rk
    if (present(ghflx)) ghflx(:,:,:) = 0._rk

    ! Check whether we're doing anything this time.

    if (time /= instance%next_time) then
       return
    else
       instance%next_time = instance%next_time + instance%mbal_tstep
    end if

    !TODO - Are these needed?
    call coordsystem_allocate(instance%lgrid, thck_temp)
    call coordsystem_allocate(instance%lgrid, calve_temp)

    ice_tstep = .false.

    ! Downscale input fields from global to local grid
    ! This subroutine computes instance%acab and instance%artm, the key inputs to GLIDE.

       call glint_downscaling_gcm (instance,              &
                                   qsmb_g,      tsfc_g,   &
                                   topo_g,      gmask)

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

    ! Accumulate acab and artm

    call glint_accumulate_gcm(instance%mbal_accum,   time,        &
                              instance%acab,         instance%artm)

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'Check for ice dynamics timestep'
       write(stdout,*) 'time =', time
       write(stdout,*) 'start_time =', instance%mbal_accum%start_time
       write(stdout,*) 'mbal_step =', instance%mbal_tstep
       write(stdout,*) 'mbal_accum_time =', instance%mbal_accum_time
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

          ! Calculate the initial ice volume (scaled and converted to water equivalent)
          call glide_get_thk(instance%model,thck_temp)
          thck_temp = thck_temp * real(rhoi/rhow)

          !TODO: Determine if glint_remove_bath is needed in a CESM run. If so, fix it to work with
          !      multiple tasks.  (And decide whether the call is needed both here and above) 
          ! Get latest upper-surface elevation (needed for masking)
!!          call glide_get_usurf(instance%model, instance%local_orog)
!!          call glint_remove_bath(instance%local_orog,1,1)

          call glint_get_mbal_gcm(instance%mbal_accum, instance%mbal_accum_time,  &
                                  instance%acab,       instance%artm)
                                  
          ! Mask out non-accumulation in ice-free areas

          !TODO - Change to dp
          where(thck_temp <= 0.0 .and. instance%acab < 0.0)
             instance%acab = 0.0
          end where

          ! Set acab to zero for ocean cells (bed below sea level, no ice present)

          where (GLIDE_IS_OCEAN(instance%model%geometry%thkmask))
             instance%acab = 0.0
          endwhere

          ! Put climate inputs in the appropriate places, with conversion ----------

          ! Note on units: 
          ! For this subroutine, input acab is in m/yr; this value is multiplied 
          !  by tim0/(scyr*thk0) and copied to data%climate%acab.
          ! Input artm is in deg C; this value is copied to data%climate%artm (no unit conversion).

          !TODO - Just rhow/rhoi without 'real'?
          call glide_set_acab(instance%model, instance%acab*real(rhow/rhoi))
          call glide_set_artm(instance%model, instance%artm)

          ! This will work only for single-processor runs
          if (GLC_DEBUG .and. tasks==1) then
             il = instance%model%numerics%idiag_global
             jl = instance%model%numerics%jdiag_global
             write (stdout,*) ' '
             write (stdout,*) 'After glide_set_acab, glide_set_artm: i, j =', il, jl
             write (stdout,*) 'acab (m/y), artm (C) =', instance%acab(il,jl)*rhow/rhoi, instance%artm(il,jl)
          end if

          ! Adjust glint acab and ablt for output
 
          where (instance%acab < -thck_temp .and. thck_temp > 0.0)
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

          ! Add the calved ice to the ablation field

          !TODO - Use this to compute the solid ice runoff,grofi?
          !       Also add basal melting (bmlt) to the liquid runoff, grofl.

          call glide_get_calving(instance%model, calve_temp)
          calve_temp = calve_temp * real(rhoi/rhow)

          ! write ice sheet diagnostics at specified interval (model%numerics%dt_diag)

          call glide_write_diagnostics(instance%model,                  &
                                       instance%model%numerics%time,    &
                                       tstep_count = instance%model%numerics%timecounter)

          ! write netCDF output

          call glide_io_writeall(instance%model,instance%model)
          call glint_io_writeall(instance,instance%model)

       end do   ! instance%n_icetstep

    end if   ! time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time

    ! Output instantaneous values

    call glint_mbal_io_writeall(instance, instance%model,       &
                                outfiles = instance%out_first,  &
                                time = time*hours2years)

    ! Deallocate

    if (associated(calve_temp)) then
       deallocate(calve_temp)
       calve_temp => null()
    end if

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

    real(sp),dimension(:,:),intent(inout) :: orog !*FD Orography --- used for input and output
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

    if (orog(x,y) < 0.0) orog(x,y)=0.0
    call glint_find_bath(orog,x,y,nx,ny)

  end subroutine glint_remove_bath

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  recursive subroutine glint_find_bath(orog,x,y,nx,ny)

    !*FD Recursive subroutine called by {\tt glimmer\_remove\_bath}.

    real(sp),dimension(:,:),intent(inout) :: orog  !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y   !*FD Starting point
    integer,                intent(in)    :: nx,ny !*FD Size of array {\tt orography}

    integer,dimension(4) :: xi=(/ -1,1,0,0 /)
    integer,dimension(4) :: yi=(/ 0,0,-1,1 /)
    integer :: ns=4,i

    do i=1,ns
       if (x+xi(i) <= nx.and.x+xi(i) > 0.and. &
            y+yi(i) <= ny.and.y+yi(i) > 0) then
          if (orog(x+xi(i),y+yi(i)) < 0.0) then
             orog(x+xi(i),y+yi(i))=0.0
             call glint_find_bath(orog,x+xi(i),y+yi(i),nx,ny)
          endif
       endif
    enddo

  end subroutine glint_find_bath

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_lapserate_dp(temp,topo,lr)

    !*FD Corrects the temperature field
    !*FD for height, using a constant lapse rate.
    !*FD
    !*FD This the double-precision version, aliased as \texttt{glimmer\_lapserate}.

    implicit none

    real(dp),dimension(:,:), intent(inout) :: temp !*FD temperature at sea-level in $^{\circ}$C
                                                   !*FD used for input and output
    real(rk),dimension(:,:), intent(in)    :: topo !*FD topography field (m above msl)
    real(rk),                intent(in)    :: lr   !*FD Lapse rate ($^{\circ}\mathrm{C\,km}^{-1}$).
                                                   !*FD
                                                   !*FD NB: the lapse rate is positive for 
                                                   !*FD falling temp with height\ldots

    temp=temp-(lr*topo/1000.0)                     ! The lapse rate calculation.

  end subroutine glint_lapserate_dp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Remove when we switch to dp

  subroutine glint_lapserate_sp(temp,topo,lr)

    !*FD Corrects the temperature field for height, using a constant lapse rate.
    !*FD
    !*FD This is the single-precision version, aliased as \texttt{glimmer\_lapserate}.

    implicit none

    real(sp),dimension(:,:),intent(inout) :: temp  !*FD temperature at sea-level in $^{\circ}$C
                                                   !*FD used for input and output
    real(rk),dimension(:,:), intent(in)    :: topo !*FD topography field (m above msl)
    real(rk),                intent(in)    :: lr   !*FD Lapse rate ($^{\circ}\mathrm{C\,km}^{-1}$).
                                                   !*FD
                                                   !*FD NB: the lapse rate is positive for 
                                                   !*FD falling temp with height\ldots

    temp=temp-(lr*topo/1000.0)                     ! The lapse rate calculation.

  end subroutine glint_lapserate_sp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_calc_precip(instance)

    use glint_precip_param
    use glimmer_log

    !*FD Process precip if necessary

    type(glint_instance) :: instance

    select case (instance%whichprecip)

    case(1)
       ! Do nothing to the precip field

    case(2)
       ! Use the Roe/Lindzen parameterisation
       call glint_precip(instance%prcp, &
                         instance%xwind, &
                         instance%ywind, &
                         instance%artm, &
                         instance%local_orog, &
                         real(instance%lgrid%delta%pt(1),rk), &
                         real(instance%lgrid%delta%pt(2),rk), &
                         fixed_a=.true.)

    case default

       call write_log('Invalid value of whichprecip',GM_FATAL,__FILE__,__LINE__)

    end select

    ! Convert from mm/s to m/s - very important!

    instance%prcp = instance%prcp*0.001

  end subroutine glint_calc_precip

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_downscaling(instance,                  &
                               g_temp,     g_temp_range,  &
                               g_precip,   g_orog,        &
                               g_zonwind,  g_merwind,     &
                               g_humid,    g_lwdown,      &
                               g_swdown,   g_airpress,    &
                               orogflag)

    use glint_interp

    !*FD Downscale global fields to the local ice sheet grid

    type(glint_instance) :: instance
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)
    real(rk),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(rk),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_humid      !*FD Global surface humidity (%)
    real(rk),dimension(:,:),intent(in)   :: g_lwdown     !*FD Global downwelling longwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_swdown     !*FD Global downwelling shortwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_airpress   !*FD Global surface air pressure (Pa)
    logical,                intent(in)   :: orogflag

    call interp_to_local(instance%lgrid_fulldomain,g_temp,      instance%downs,localsp=instance%artm)
    call interp_to_local(instance%lgrid_fulldomain,g_temp_range,instance%downs,localsp=instance%arng,z_constrain=.true.)
    call interp_to_local(instance%lgrid_fulldomain,g_precip,    instance%downs,localsp=instance%prcp,z_constrain=.true.)

    if (instance%whichacab==3) then
       call interp_to_local(instance%lgrid_fulldomain,g_humid,   instance%downs,localrk=instance%humid,z_constrain=.true.)
       call interp_to_local(instance%lgrid_fulldomain,g_lwdown,  instance%downs,localrk=instance%lwdown)
       call interp_to_local(instance%lgrid_fulldomain,g_swdown,  instance%downs,localrk=instance%swdown)
       call interp_to_local(instance%lgrid_fulldomain,g_airpress,instance%downs,localrk=instance%airpress,z_constrain=.true.)
    end if

    if (orogflag) call interp_to_local(instance%lgrid_fulldomain,g_orog,instance%downs,localdp=instance%global_orog,z_constrain=.true.)

    if (instance%whichprecip==2 .or. instance%whichacab==3) &
         call interp_wind_to_local(instance%lgrid_fulldomain,g_zonwind,g_merwind,instance%downs,instance%xwind,instance%ywind)

  end subroutine glint_downscaling

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_downscaling_gcm (instance,            &
                                    qsmb_g,     tsfc_g,  &
                                    topo_g,     gmask)
 
    use glimmer_paramets, only: thk0, GLC_DEBUG

    use glint_type
    use glint_interp, only: interp_to_local
    use parallel, only: tasks

    ! Downscale fields from the global grid (with multiple elevation classes)
    ! to the local ice sheet grid.
    ! 
    ! This routine is used for downscaling when the surface mass balance is
    ! computed in the GCM land surface model.

    type(glint_instance), intent(inout) :: instance
    real(dp),dimension(:,:,:),intent(in) :: qsmb_g       ! Surface mass balance (m)
    real(dp),dimension(:,:,:),intent(in) :: tsfc_g       ! Surface temperature (C)
    real(dp),dimension(:,:,:),intent(in) :: topo_g       ! Surface elevation (m)
    integer ,dimension(:,:),  intent(in),optional :: gmask ! = 1 where global data are valid
                                                           ! = 0 elsewhere

    real(dp), parameter :: maskval = 0.0_dp    ! value written to masked out gridcells

    integer ::       &
       nec,          &      ! number of elevation classes
       nxl, nyl             ! local grid dimensions

    integer :: i, j, n, ig, jg
 
    real(dp), dimension(:,:,:), allocatable ::   &
       qsmb_l,    &! interpolation of global mass balance to local grid
       tsfc_l,    &! interpolation of global sfc temperature to local grid
       topo_l      ! interpolation of global topography in each elev class to local grid

    real(dp) :: fact, usrf

    real(dp), parameter :: lapse = 0.0065_dp   ! atm lapse rate, deg/m
                                               ! used only for extrapolating temperature outside
                                               !  the range provided by the climate model
    nec = size(qsmb_g,3)
    nxl = instance%lgrid%size%pt(1)
    nyl = instance%lgrid%size%pt(2)

    allocate(qsmb_l(nxl,nyl,nec))
    allocate(tsfc_l(nxl,nyl,nec))
    allocate(topo_l(nxl,nyl,nec))

    !   Downscale global fields for each elevation class to local grid.

    if (present(gmask)) then   ! set local field = maskval where the global field is masked out
                               ! (i.e., where instance%downs%lmask = 0).
       do n = 1, nec
          call interp_to_local(instance%lgrid_fulldomain, qsmb_g(:,:,n), instance%downs, localdp=qsmb_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
          call interp_to_local(instance%lgrid_fulldomain, tsfc_g(:,:,n), instance%downs, localdp=tsfc_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
          call interp_to_local(instance%lgrid_fulldomain, topo_g(:,:,n), instance%downs, localdp=topo_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
       enddo

    else    ! global field values are assumed to be valid everywhere

       do n = 1, nec
          call interp_to_local(instance%lgrid_fulldomain, qsmb_g(:,:,n), instance%downs, localdp=qsmb_l(:,:,n))
          call interp_to_local(instance%lgrid_fulldomain, tsfc_g(:,:,n), instance%downs, localdp=tsfc_l(:,:,n))
          call interp_to_local(instance%lgrid_fulldomain, topo_g(:,:,n), instance%downs, localdp=topo_l(:,:,n))
       enddo

    endif

    ! The following output only works correctly if running with a single task
    if (GLC_DEBUG .and. tasks==1) then
       ig = iglint_global   ! in glint_type; make sure values are appropriate
       jg = jglint_global
       write (stdout,*) ' ' 
       write (stdout,*) 'Interpolate fields to local grid'
       write (stdout,*) 'Global cell =', ig, jg
       do n = 1, nec
          write(stdout,*) n, topo_g(ig,jg, n)
       enddo

       do j = 1, nyl
       do i = 1, nxl
           if ( (instance%downs%xloc(i,j,1) == ig .and. instance%downs%yloc(i,j,1) == jg) .or.  &
                (instance%downs%xloc(i,j,2) == ig .and. instance%downs%yloc(i,j,2) == jg) .or.  &
                (instance%downs%xloc(i,j,3) == ig .and. instance%downs%yloc(i,j,3) == jg) .or.  &
                (instance%downs%xloc(i,j,4) == ig .and. instance%downs%yloc(i,j,4) == jg) ) then
               write(stdout,*) i, j, thk0 * instance%model%geometry%usrf(i,j)
           endif
       enddo
       enddo
    
       i = instance%model%numerics%idiag_global
       j = instance%model%numerics%jdiag_global
       write (stdout,*) ' ' 
       write (stdout,*) 'Interpolated to local cells: i, j =', i, j
       do n = 1, nec
          write (stdout,*) ' '
          write (stdout,*) 'n =', n
          write (stdout,*) 'qsmb_l =', qsmb_l(i,j,n)
          write (stdout,*) 'tsfc_l =', tsfc_l(i,j,n)
          write (stdout,*) 'topo_l =', topo_l(i,j,n)
       enddo

    end if ! GLC_DEBUG

!   Interpolate tsfc and qsmb to local topography using values in the neighboring 
!    elevation classes.
!   If the local topography is outside the bounds of the global elevations classes,
!    extrapolate the temperature using the prescribed lapse rate.

    do j = 1, nyl
    do i = 1, nxl

       usrf = instance%model%geometry%usrf(i,j) * thk0   ! actual sfc elevation (m)

       if (usrf <= topo_l(i,j,1)) then
          instance%acab(i,j) = qsmb_l(i,j,1)
          instance%artm(i,j) = tsfc_l(i,j,1) + lapse*(topo_l(i,j,1)-usrf)
       elseif (usrf > topo_l(i,j,nec)) then
          instance%acab(i,j) = qsmb_l(i,j,nec)
          instance%artm(i,j) = tsfc_l(i,j,nec) - lapse*(usrf-topo_l(i,j,nec))
       else
          do n = 2, nec
             if (usrf > topo_l(i,j,n-1) .and. usrf <= topo_l(i,j,n)) then
                fact = (topo_l(i,j,n) - usrf) / (topo_l(i,j,n) - topo_l(i,j,n-1)) 
                instance%acab(i,j) = fact*qsmb_l(i,j,n-1) + (1._dp-fact)*qsmb_l(i,j,n)
                instance%artm(i,j) = fact*tsfc_l(i,j,n-1) + (1._dp-fact)*tsfc_l(i,j,n)
                exit
             endif
          enddo
       endif   ! usrf

       ! The following output only works correctly if running with a single task
       if (GLC_DEBUG .and. tasks==1) then
          if (i==instance%model%numerics%idiag_global .and. j==instance%model%numerics%jdiag_global) then
             n = 4  
             write (stdout,*) ' '
             write (stdout,*) 'Interpolated values, i, j, n =', i, j, n
             write (stdout,*) 'usrf =', usrf
             write (stdout,*) 'acab =', instance%acab(i,j)
             write (stdout,*) 'artm =', instance%artm(i,j)
             write (stdout,*) 'topo(n-1) =', topo_l(i,j,n-1)
             write (stdout,*) 'topo(n) =', topo_l(i,j,n)
             write (stdout,*) 'qsmb(n-1) =', qsmb_l(i,j,n-1)
             write (stdout,*) 'qsmb(n) =', qsmb_l(i,j,n)
             write (stdout,*) 'tsfc(n-1) =', tsfc_l(i,j,n-1)
             write (stdout,*) 'tsfc(n) =', tsfc_l(i,j,n)
             write (stdout,*) 'fact = ', (topo_l(i,j,n) - usrf) / (topo_l(i,j,n) - topo_l(i,j,n-1)) 
          endif
       end if

    enddo  ! i
    enddo  ! j

    deallocate(qsmb_l, tsfc_l, topo_l)

  end subroutine glint_downscaling_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Is there a better place to put this subroutine?

  subroutine get_i_upscaled_fields(instance,                    &
                                   orog,         albedo,        &
                                   ice_frac,     veg_frac,      &
                                   snowice_frac, snowveg_frac,  &
                                   snow_depth)

    !*FD Upscales and returns certain fields
    !*FD Output fields are only valid on the main task
    !*FD 
    !*FD \begin{itemize}
    !*FD \item \texttt{orog} --- the orographic elevation (m)
    !*FD \item \texttt{albedo} --- the albedo of ice/snow (this is only a notional value --- need to do
    !*FD some work here)
    !*FD \item \texttt{ice\_frac} --- The fraction covered by ice
    !*FD \item \texttt{veg\_frac} --- The fraction of exposed vegetation
    !*FD \item \texttt{snowice\_frac} --- The fraction of snow-covered ice
    !*FD \item \texttt{snowveg\_frac} --- The fraction of snow-covered vegetation
    !*FD \item \texttt{snow_depth} --- The mean snow-depth over those parts covered in snow (m w.e.)
    !*FD \end{itemize}

    use glimmer_paramets

    ! Arguments ----------------------------------------------------------------------------------------

    type(glint_instance),   intent(in)  :: instance      !*FD the model instance

    real(rk),dimension(:,:),intent(out) :: orog          !*FD the orographic elevation (m)
    real(rk),dimension(:,:),intent(out) :: albedo        !*FD the albedo of ice/snow
    real(rk),dimension(:,:),intent(out) :: ice_frac      !*FD The fraction covered by ice
    real(rk),dimension(:,:),intent(out) :: veg_frac      !*FD The fraction of exposed vegetation
    real(rk),dimension(:,:),intent(out) :: snowice_frac  !*FD The fraction of snow-covered ice
    real(rk),dimension(:,:),intent(out) :: snowveg_frac  !*FD The fraction of snow-covered vegetation
    real(rk),dimension(:,:),intent(out) :: snow_depth    !*FD The mean snow-depth over those 
    !*FD parts covered in snow (m w.e.)

    ! Internal variables -------------------------------------------------------------------------------

    real(rk),dimension(:,:),pointer :: temp => null()

    ! --------------------------------------------------------------------------------------------------
    ! Orography

    call mean_to_global(instance%ups_orog, &
                        instance%model%geometry%usrf, &
                        orog,    &
                        instance%out_mask)
    orog=thk0*orog

    call coordsystem_allocate(instance%lgrid,temp)

    !TODO - Change to dp
    ! Ice-no-snow fraction
    where (instance%mbal_accum%snowd==0.0.and.instance%model%geometry%thck>0.0)
       temp=1.0
    elsewhere
       temp=0.0
    endwhere
    call mean_to_global(instance%ups, &
                        temp, &
                        ice_frac,    &
                        instance%out_mask)

    ! Ice-with-snow fraction
    where (instance%mbal_accum%snowd>0.0.and.instance%model%geometry%thck>0.0)
       temp=1.0
    elsewhere
       temp=0.0
    endwhere
    call mean_to_global(instance%ups, &
                        temp, &
                        snowice_frac,    &
                        instance%out_mask)

    ! Veg-with-snow fraction (if ice <10m thick)
    where (instance%mbal_accum%snowd>0.0.and.instance%model%geometry%thck<=(10.0/thk0))
       temp=1.0
    elsewhere
       temp=0.0
    endwhere
    call mean_to_global(instance%ups, &
                        temp, &
                        snowveg_frac,    &
                        instance%out_mask)

    ! Remainder is veg only
    veg_frac=1.0-ice_frac-snowice_frac-snowveg_frac

    ! Snow depth

    call mean_to_global(instance%ups, &
                        instance%mbal_accum%snowd, &
                        snow_depth,    &
                        instance%out_mask)

    ! Albedo

    where ((ice_frac+snowice_frac)>0.0)
       albedo=instance%ice_albedo
    elsewhere
       albedo=0.0
    endwhere

    deallocate(temp)
    temp => null()

  end subroutine get_i_upscaled_fields

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Is there a better place to put this subroutine?

  subroutine get_i_upscaled_fields_gcm(instance,    nec,      &
                                       nxl,         nyl,      &
                                       nxg,         nyg,      &
                                       gfrac,       gtopo,    &
                                       grofi,       grofl,    &
                                       ghflx)

    ! Upscale fields from the local grid to the global grid (with multiple elevation classes).
    ! Output fields are only valid on the main task.
    ! The upscaled fields are passed to the GCM land surface model, which has the option
    !  of updating the fractional area and surface elevation of glaciated gridcells.

    use glimmer_paramets, only: thk0, GLC_DEBUG
    use glimmer_log
    use parallel, only: tasks, main_task

    ! Arguments ----------------------------------------------------------------------------
 
    type(glint_instance),     intent(in)  :: instance      ! the model instance
    integer,                  intent(in)  :: nec           ! number of elevation classes
    integer,                  intent(in)  :: nxl,nyl       ! local grid dimensions 
    integer,                  intent(in)  :: nxg,nyg       ! global grid dimensions 

    !TODO - Should these be inout?
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gfrac   ! ice-covered fraction [0,1]
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gtopo   ! surface elevation (m)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: grofi   ! ice runoff (calving) flux (kg/m^2/s)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: grofl   ! liquid runoff (basal melt) flux (kg/m^2/s)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: ghflx   ! heat flux (m)
 
    ! Internal variables ----------------------------------------------------------------------
 
    real(dp),dimension(nxl,nyl) :: local_field
    real(dp),dimension(nxl,nyl) :: local_topo   ! local surface elevation (m)
    real(dp),dimension(nxl,nyl) :: local_thck   ! local ice thickness (m)

    !TODO - Put this parameter elsewhere?  Make it equal to thkmin?
    real(dp), parameter :: min_thck = 1.0_dp    ! min thickness (m) for setting gfrac = 1

    integer :: i, j            ! indices
 
    integer :: il, jl, ig, jg
    character(len=100) :: message

    !TODO - Pass in topomax as an argument instead of hardwiring it here
    real(dp), dimension(0:nec) :: topomax   ! upper elevation limit of each class

    ! Given the value of nec, specify the upper and lower elevation boundaries of each class.
    ! Note: These must be consistent with the values in the GCM.  Better to pass as an argument.
    if (nec == 1) then
       topomax = (/ 0._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 3) then
       topomax = (/ 0._dp,  1000._dp,  2000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 5) then
       topomax = (/ 0._dp,   500._dp,  1000._dp,  1500._dp,  2000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 10) then
       topomax = (/ 0._dp,   200._dp,   400._dp,   700._dp,  1000._dp,  1300._dp,  &
                            1600._dp,  2000._dp,  2500._dp,  3000._dp, 10000._dp /)
    elseif (nec == 36) then
       topomax = (/ 0._dp,   200._dp,   400._dp,   600._dp,   800._dp,  &
                 1000._dp,  1200._dp,  1400._dp,  1600._dp,  1800._dp,  &
                 2000._dp,  2200._dp,  2400._dp,  2600._dp,  2800._dp,  &
                 3000._dp,  3200._dp,  3400._dp,  3600._dp,  3800._dp,  &
                 4000._dp,  4200._dp,  4400._dp,  4600._dp,  4800._dp,  &
                 5000._dp,  5200._dp,  5400._dp,  5600._dp,  5800._dp,  &
                 6000._dp,  6200._dp,  6400._dp,  6600._dp,  6800._dp,  &
                 7000._dp, 10000._dp /)
    else
       if (GLC_DEBUG .and. main_task) then
          write(message,'(a6,i3)') 'nec =', nec
          call write_log(trim(message), GM_DIAGNOSTIC)
       end if
       call write_log('ERROR: Current supported values of nec (no. of elevation classes) are 1, 3, 5, 10, or 36', &
                       GM_FATAL,__FILE__,__LINE__)
    endif

    local_topo(:,:) = thk0 * instance%model%geometry%usrf(:,:)
    local_thck(:,:) = thk0 * instance%model%geometry%thck(:,:)

    ! The following output only works correctly if running with a single task
    if (GLC_DEBUG .and. tasks==1) then
       ig = iglint_global    ! defined in glint_type
       jg = jglint_global
       il = instance%model%numerics%idiag_global
       jl = instance%model%numerics%jdiag_global
       write(stdout,*) 'In get_i_upscaled_fields_gcm'
       write(stdout,*) 'il, jl =', il, jl
       write(stdout,*) 'ig, jg =', ig, jg
       write(stdout,*) 'nxl, nyl =', nxl,nyl
       write(stdout,*) 'nxg, nyg =', nxg,nyg
       write(stdout,*) 'topo =', local_topo(il,jl) 
       write(stdout,*) 'thck =', local_thck(il,jl) 
       write(stdout,*) 'local out_mask =', instance%out_mask(il,jl)
    end if

    ! temporary field: = 1 where ice thickness exceeds threshold, else = 0
    !TODO - Use > 0 as threshold?

    do j = 1, nyl
    do i = 1, nxl
       if (local_thck(i,j) > min_thck) then
          local_field(i,j) = 1.d0
       else
          local_field(i,j) = 0.d0
       endif
    enddo
    enddo

    ! ice fraction
    !TODO - gfrac should be fraction of total grid cell with ice in each elevation class
    !       Currently is ice-covered fraction of cells in a given elevation class

    call mean_to_global_mec(instance%ups,                       &
                            nxl,                nyl,            &
                            nxg,                nyg,            &
                            nec,                topomax,        &
                            local_field,        gfrac,          &
                            local_topo,         instance%out_mask)

    ! surface elevation

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_topo,          gtopo,     &
                            local_topo,          instance%out_mask)

    !TODO - For upscaling, need to copy the appropriate Glide fields into the local_field array

    ! ice runoff

    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         grofi,     &
                            local_topo,          instance%out_mask)

    ! liquid runoff

    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         grofl,     &
                            local_topo,          instance%out_mask)

    ! heat flux

    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         ghflx,     &
                            local_topo,          instance%out_mask)
    
    if (GLC_DEBUG .and. main_task) then

!       write(stdout,*) ' '
!       write(stdout,*) 'global ifrac:'
!       do n = 1, nec
!          write(stdout,*) n, gfrac(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global gtopo:'
!       do n = 1, nec
!          write(stdout,*) n, gtopo(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global grofi:'
!       do n = 1, nec
!          write(stdout,*) n, grofi(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global grofl:'
!       do n = 1, nec
!          write(stdout,*) n, grofl(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global ghflx:'
!       do n = 1, nec
!          write(stdout,*) n, ghflx(ig, jg, n)
!       enddo

    end if

  end subroutine get_i_upscaled_fields_gcm

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glint_timestep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


