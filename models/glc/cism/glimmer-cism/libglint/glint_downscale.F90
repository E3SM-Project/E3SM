!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_downscale.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

  module glint_downscale

  ! This module contains subroutines for downscaling fields from the global to the local grid.
  ! Much of the actual work is done at a lower level, in glint_interp.F90.

  use glint_type
  use glint_constants
  use glimmer_global, only: dp
  implicit none

  private
  public glint_downscaling, glint_downscaling_gcm,  &
         glint_init_input_gcm, glint_accumulate_input_gcm, glint_average_input_gcm

  !Note: The three subroutines glint_*_input_gcm are based on old Glint subroutines
  !      in glint_mbal_coupling.F90.

contains

  subroutine glint_downscaling(instance,                  &
                               g_temp,     g_temp_range,  &
                               g_precip,   g_orog,        &
                               g_zonwind,  g_merwind,     &
                               g_humid,    g_lwdown,      &
                               g_swdown,   g_airpress,    &
                               orogflag)

    use glint_interp

    !*FD Downscale global input fields to the local ice sheet grid

    type(glint_instance) :: instance
    real(dp),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(dp),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(dp),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)
    real(dp),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(dp),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(dp),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(dp),dimension(:,:),intent(in)   :: g_humid      !*FD Global surface humidity (%)
    real(dp),dimension(:,:),intent(in)   :: g_lwdown     !*FD Global downwelling longwave (W/m^2)
    real(dp),dimension(:,:),intent(in)   :: g_swdown     !*FD Global downwelling shortwave (W/m^2)
    real(dp),dimension(:,:),intent(in)   :: g_airpress   !*FD Global surface air pressure (Pa)
    logical,                intent(in)   :: orogflag

    call interp_to_local(instance%lgrid_fulldomain,g_temp,      instance%downs,localdp=instance%artm)
    call interp_to_local(instance%lgrid_fulldomain,g_temp_range,instance%downs,localdp=instance%arng,z_constrain=.true.)
    call interp_to_local(instance%lgrid_fulldomain,g_precip,    instance%downs,localdp=instance%prcp,z_constrain=.true.)

    if (instance%whichacab==MASS_BALANCE_EBM) then
       call interp_to_local(instance%lgrid_fulldomain,g_humid,   instance%downs,localdp=instance%humid,z_constrain=.true.)
       call interp_to_local(instance%lgrid_fulldomain,g_lwdown,  instance%downs,localdp=instance%lwdown)
       call interp_to_local(instance%lgrid_fulldomain,g_swdown,  instance%downs,localdp=instance%swdown)
       call interp_to_local(instance%lgrid_fulldomain,g_airpress,instance%downs,localdp=instance%airpress,z_constrain=.true.)
    end if

    if (orogflag) call interp_to_local(instance%lgrid_fulldomain,g_orog,instance%downs,localdp=instance%global_orog,z_constrain=.true.)

    if (instance%whichprecip==PRECIP_RL .or. instance%whichacab==MASS_BALANCE_EBM) &
         call interp_wind_to_local(instance%lgrid_fulldomain,g_zonwind,g_merwind,instance%downs,instance%xwind,instance%ywind)

  end subroutine glint_downscaling

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_downscaling_gcm (instance,&
                                    qsmb_g,     tsfc_g,  &
                                    topo_g,     gmask)
 
    use glimmer_paramets, only: thk0, GLC_DEBUG
    use glint_constants, only: lapse
    use glint_type
    use glint_interp, only: interp_to_local, copy_to_local
    use parallel, only: tasks

    ! Downscale global input fields from the global grid (with multiple elevation classes)
    ! to the local ice sheet grid.
    ! 
    ! This routine is used for downscaling when the surface mass balance is
    ! computed in the GCM land surface model.
    
    type(glint_instance), intent(inout) :: instance
    real(dp),dimension(:,:,0:),intent(in) :: qsmb_g       ! Surface mass balance (m)
    real(dp),dimension(:,:,0:),intent(in) :: tsfc_g       ! Surface temperature (C)
    real(dp),dimension(:,:,0:),intent(in) :: topo_g       ! Surface elevation (m)
    integer ,dimension(:,:),   intent(in),optional :: gmask   ! = 1 where global data are valid
                                                             ! = 0 elsewhere
    real(dp), parameter :: maskval = 0.0_dp    ! value written to masked out gridcells

    integer :: nxl, nyl, nec             ! local grid dimensions

    integer :: i, j, n
 
    real(dp), dimension(:,:,:), allocatable ::   &
       qsmb_l,    &! interpolation of global mass balance to local grid
       tsfc_l,    &! interpolation of global sfc temperature to local grid
       topo_l      ! interpolation of global topography in each elev class to local grid

    real(dp) :: fact, usrf, thck
                                               
    nxl = instance%lgrid%size%pt(1)
    nyl = instance%lgrid%size%pt(2)
    nec = ubound(qsmb_g,3)

    allocate(qsmb_l(nxl,nyl,0:nec))
    allocate(tsfc_l(nxl,nyl,0:nec))
    allocate(topo_l(nxl,nyl,0:nec))

    !   Downscale global fields for each elevation class to local grid (horizontal interpolation).

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

    ! For elevation class 0 (bare land), simply set the values to the values of the
    ! global parent cell. No vertical/horizontal interpolation is used, since these
    ! elevation-dependent values are not constrained to a discrete elevation band. Also
    ! note that we do not consider gmask here.

    call copy_to_local(instance%lgrid_fulldomain, qsmb_g(:,:,0), instance%downs, qsmb_l(:,:,0))
    call copy_to_local(instance%lgrid_fulldomain, tsfc_g(:,:,0), instance%downs, tsfc_l(:,:,0))

    ! topo_l(:,:,0) isn't used right now, but compute it anyway for consistency
    call copy_to_local(instance%lgrid_fulldomain, topo_g(:,:,0), instance%downs, topo_l(:,:,0))

!   Interpolate tsfc and qsmb to local topography using values in the neighboring 
!    elevation classes (vertical interpolation).

!   If the local topography is outside the bounds of the global elevations classes,
!    extrapolate the temperature using the prescribed lapse rate.

    do j = 1, nyl
       do i = 1, nxl

          usrf = instance%model%geometry%usrf(i,j) * thk0
          thck = instance%model%geometry%thck(i,j) * thk0

          if (thck <= min_thck) then !if ice-free...
             if (usrf > 0.) then !and on land (not ocean)...
                ! As noted above, no vertical interpolation is done for ice-free land
                instance%acab(i,j) = qsmb_l(i,j,0)
                instance%artm(i,j) = tsfc_l(i,j,0)
                if (instance%acab(i,j) < 0.d0) then
                   write (stdout,*)'ERROR: SMB is negative over bare land point'
                   write (stdout,*)'instance%acab(i,j) = ',instance%acab(i,j)
                   write (stdout,*)'instance%artm(i,j) = ',instance%artm(i,j)          
                   write (stdout,*)'qsmb_l(i,j,0) = ',qsmb_l(i,j,0)
                   write (stdout,*)'usrf=',usrf
                   write (stdout,*)'thck=',thck
                   write (stdout,*)'Stopping in glint_downscale.'
                   stop      
                endif
             else  ! usrf <= 0  -- assumed to be ocean
                ! CISM assumes any point with usrf <= 0 is ocean, and thus can't form ice
                ! (actually, this isn't exactly the cutoff used elsewhere in CISM - we may
                ! want to change this conditional to use GLIDE_IS_OCEAN). So it makes no
                ! sense to pass acab and artm there. However, this could lead to a loss of
                ! conservation, e.g., if CLM thinks a grid cell is bare ground with some
                ! positive SMB (glacial inception) yet CISM says it's ocean (so ignores
                ! the SMB). We eventually want to handle this by keeping CLM consistent
                ! with CISM in terms of its breakdown into land vs "ocean" (e.g., wetland
                ! in CLM). In that case, if CISM says a point is ocean, then it would
                ! tell CLM that that point is ocean, and so CLM wouldn't try to generate
                ! SMB there.

                instance%acab(i,j) = 0.d0
                instance%artm(i,j) = 0.d0
             endif
          else !if ice-covered...
             if (usrf <= topo_l(i,j,1)) then
                instance%acab(i,j) = qsmb_l(i,j,1)
                instance%artm(i,j) = tsfc_l(i,j,1) + lapse*(topo_l(i,j,1)-usrf)
             elseif (usrf > topo_l(i,j,nec)) then
                instance%acab(i,j) = qsmb_l(i,j,nec)
                instance%artm(i,j) = tsfc_l(i,j,nec) - lapse*(usrf-topo_l(i,j,nec))
             else
                do n = 2,nec
                   if (usrf > topo_l(i,j,n-1) .and. usrf <= topo_l(i,j,n)) then
                      fact = (topo_l(i,j,n) - usrf) / (topo_l(i,j,n) - topo_l(i,j,n-1)) 
                      instance%acab(i,j) = fact*qsmb_l(i,j,n-1) + (1._dp-fact)*qsmb_l(i,j,n)
                      instance%artm(i,j) = fact*tsfc_l(i,j,n-1) + (1._dp-fact)*tsfc_l(i,j,n)
                      exit
                   endif
                enddo
             endif ! usrf, inner
          endif ! thck  
       enddo ! i
    enddo ! j

    deallocate(qsmb_l, tsfc_l, topo_l)
    
  end subroutine glint_downscaling_gcm

  ! +++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_init_input_gcm(params,  &
                                  lgrid,   &
                                  whichacab)

    ! Simplified version of glint_mbc_init, used when coupling
    ! to a GCM that provides the surface mass balance and temperature

    use glimmer_coordinates
    use glint_constants, only: years2hours

    type(glint_mbc)        :: params      ! mass balance parameters
    type(coordsystem_type) :: lgrid       ! local grid
    integer, intent(in)    :: whichacab   ! mass balance method
                                          ! = 0 for GCM coupling

    ! Deallocate if necessary

    if (associated(params%acab_save))  deallocate(params%acab_save)
    if (associated(params%artm_save))  deallocate(params%artm_save)
    if (associated(params%acab))       deallocate(params%acab)
    if (associated(params%artm))       deallocate(params%artm)

    ! Allocate arrays and zero

    call coordsystem_allocate(lgrid,params%acab_save);  params%acab_save = 0.d0
    call coordsystem_allocate(lgrid,params%artm_save);  params%artm_save = 0.d0
    call coordsystem_allocate(lgrid,params%acab);       params%acab = 0.d0
    call coordsystem_allocate(lgrid,params%artm);       params%artm = 0.d0

    ! Set the mbal method and tstep

    params%mbal%which = whichacab     ! = 0 for GCM coupling
    params%mbal%tstep = years2hours   ! no. of hours in 1 year

  end subroutine glint_init_input_gcm

  ! +++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_accumulate_input_gcm(params, time, acab, artm)

    type(glint_mbc)  :: params
    integer :: time

    real(dp),dimension(:,:),intent(in) :: acab   ! Surface mass balance (m)
    real(dp),dimension(:,:),intent(in) :: artm   ! Mean air temperature (degC)

    ! Things to do the first time

    if (params%new_accum) then

       params%new_accum = .false.
       params%av_count  = 0

       ! Initialise

       params%acab_save = 0.d0
       params%artm_save = 0.d0
       params%start_time = time

    end if

    params%av_count = params%av_count + 1

    ! Accumulate

    params%acab_save = params%acab_save + acab
    params%artm_save = params%artm_save + artm

    ! Copy instantaneous fields

    params%acab = acab
    params%artm = artm

  end subroutine glint_accumulate_input_gcm

  ! +++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_average_input_gcm(params, dt, acab, artm)

    use glint_constants, only: hours2years

    type(glint_mbc)  :: params
    integer,                intent(in)    :: dt     !*FD mbal accumulation time (hours)
    real(dp),dimension(:,:),intent(out)   :: artm   !*FD Mean air temperature (degC)
    real(dp),dimension(:,:),intent(out)   :: acab   !*FD Mass-balance (m/yr)

    if (.not. params%new_accum) then
       params%artm_save = params%artm_save / real(params%av_count,dp)
    end if
    artm  = params%artm_save

    ! Note: acab_save has units of m, but acab has units of m/yr
    acab  = params%acab_save / real(dt*hours2years,dp)

    params%new_accum = .true.

  end subroutine glint_average_input_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module glint_downscale

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
