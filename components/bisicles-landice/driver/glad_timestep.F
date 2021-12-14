#ifdef CPRIBM
@PROCESS ALIAS_SIZE(107374182)
#endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_timestep.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glad_timestep
  !> timestep of a GLAD instance

  use glad_type
  use glad_constants
  use glimmer_global, only: dp
  implicit none

  private
  public glad_i_tstep_gcm

contains


  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glad_i_tstep_gcm(time,            instance,       &
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
    use glad_mbal_coupling, only : glad_accumulate_input_gcm, glad_average_input_gcm
    use glad_io
    use glad_mbal_io
    use glide_diagnostics
    use parallel, only: tasks, main_task, this_rank
    use glad_output_fluxes, only : accumulate_output_fluxes, reset_output_fluxes
    
    implicit none

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    integer,                intent(in)   :: time         ! Current time in hours
    type(glad_instance), intent(inout)  :: instance     ! Model instance
    logical,                intent(out)  :: ice_tstep    ! Set if we have done an ice time step

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    real(dp),dimension(:,:),pointer :: thck_temp => null() ! temporary array for volume calcs

    integer :: i, il, jl

    if (GLC_DEBUG .and. main_task) then
       print*, 'In glad_i_tstep_gcm'
    endif

    ice_tstep = .false.

    call coordsystem_allocate(instance%lgrid, thck_temp)

    ! ------------------------------------------------------------------------
    ! Sort out some local orography and remove bathymetry. This relies on the 
    ! point 1,1 being underwater. However, it's a better method than just 
    ! setting all points < 0.0 to zero
    ! ------------------------------------------------------------------------  

    !Note: Call to glad_remove_bath is commented out for now.  Not sure if it is needed in GCM runs.
!!    call glide_get_usurf(instance%model, instance%local_orog)
!!    call glad_remove_bath(instance%local_orog,1,1)

    ! Get ice thickness ----------------------------------------

    call glide_get_thk(instance%model,thck_temp)

    ! Accumulate Glide input fields, acab and artm
    ! Note: At this point, instance%acab has units of m
    !       Upon averaging (in glad_average_input_gcm), units are converted to m/yr

    call glad_accumulate_input_gcm(instance%mbal_accum,   time,        &
                                    instance%acab,         instance%artm)


    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) ' '
       write(stdout,*) 'In glad_i_tstep_gcm, time =', time
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

       call reset_output_fluxes(instance%glad_output_fluxes)

       ! ---------------------------------------------------------------------
       ! Timestepping for ice sheet model
       ! ---------------------------------------------------------------------

       do i = 1, instance%n_icetstep

          if (GLC_DEBUG .and. main_task) then
             write (stdout,*) 'Ice sheet timestep, iteration =', i
          end if

          ! Get average values of acab and artm during mbal_accum_time
          ! instance%acab has units of m/yr w.e. after averaging

          call glad_average_input_gcm(instance%mbal_accum, instance%mbal_accum_time,  &
                                       instance%acab,       instance%artm)
                                  
          ! Calculate the initial ice volume (scaled and converted to water equivalent)
          call glide_get_thk(instance%model,thck_temp)
          thck_temp = thck_temp * rhoi/rhow

          !Note: Call to glad_remove_bath is commented out for now.  Not sure if it is needed in GCM runs.
          ! Get latest upper-surface elevation (needed for masking)
!!          call glide_get_usurf(instance%model, instance%local_orog)
!!          call glad_remove_bath(instance%local_orog,1,1)

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

          !TODO - It is confusing to have units of m/yr w.e. for instance%acab, compared to units m/yr ice for Glide. 
          !       Change to use the same units consistently?  E.g., switch to w.e. in Glide

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

          ! Adjust glad acab for output
 
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
          call glad_io_writeall(instance,instance%model)

          ! Accumulate Glide output fields to be sent to GCM

          call accumulate_output_fluxes(instance%glad_output_fluxes, instance%model)

       end do   ! instance%n_icetstep

    end if   ! time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time

!WHL - debug
    print*, 'output instantaneous values'

    ! Output instantaneous values

    call glad_mbal_io_writeall(instance, instance%model,       &
                                outfiles = instance%out_first,  &
                                time = time*hours2years)

    ! Deallocate

    if (associated(thck_temp)) then
       deallocate(thck_temp)
       thck_temp => null()
    endif

    if (GLC_DEBUG .and. main_task) then
       write(stdout,*) 'Done in glad_i_tstep_gcm'
    endif

  end subroutine glad_i_tstep_gcm

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Rewrite glad_remove_bath to support multiple tasks?
  !       Calls to this subroutine are currently commented out.

  subroutine glad_remove_bath(orog,x,y)

    ! Sets ocean areas to zero height, working recursively from
    ! a known ocean point.

    use glimmer_log
    use parallel, only : tasks

    real(dp),dimension(:,:),intent(inout) :: orog !> Orography --- used for input and output
    integer,                intent(in)    :: x,y  !> Location of starting point (index)

    integer :: nx,ny

    ! Currently, this routine is called assuming point 1,1 is ocean... this won't be true
    ! when running on multiple processors, with a distributed grid
    ! This can't be made a fatal error, because this is currently called even if we have
    ! more than one task... the hope is just that the returned data aren't needed in CESM.
    if (tasks > 1) then
       call write_log('Use of glad_remove_bath currently assumes the use of only one task', &
                      GM_WARNING, __FILE__, __LINE__)
    end if

    nx=size(orog,1) ; ny=size(orog,2)

    if (orog(x,y) < 0.d0) orog(x,y) = 0.d0
    call glad_find_bath(orog,x,y,nx,ny)

  end subroutine glad_remove_bath

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  recursive subroutine glad_find_bath(orog,x,y,nx,ny)

    !> Recursive subroutine called by {\tt glimmer\_remove\_bath}.

    real(dp),dimension(:,:),intent(inout) :: orog  !> Orography --- used for input and output
    integer,                intent(in)    :: x,y   !> Starting point
    integer,                intent(in)    :: nx,ny !> Size of array {\tt orography}

    integer,dimension(4) :: xi = (/ -1,1,0,0 /)
    integer,dimension(4) :: yi = (/ 0,0,-1,1 /)
    integer :: ns = 4
    integer :: i

    do i=1,ns
       if (x+xi(i) <= nx .and. x+xi(i) > 0 .and. &
           y+yi(i) <= ny .and. y+yi(i) > 0) then
          if (orog(x+xi(i),y+yi(i)) < 0.d0) then
             orog(x+xi(i),y+yi(i)) = 0.d0
             call glad_find_bath(orog,x+xi(i),y+yi(i),nx,ny)
          endif
       endif
    enddo

  end subroutine glad_find_bath

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module glad_timestep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
