!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_upscale.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

  module glint_upscale

  ! This module contains subroutines for upscaling fields from the local to the global grid.
  ! Much of the actual work is done at a lower level, in glint_interp.F90.

  use glint_type
  use glint_constants
  use glimmer_global, only: dp
  implicit none

  private
  public glint_upscaling, glint_upscaling_gcm,  &
         glint_accumulate_output_gcm

contains

  subroutine glint_upscaling(instance,                    &
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
    use glimmer_coordinates, only: coordsystem_allocate

    ! Arguments ----------------------------------------------------------------------------------------

    type(glint_instance),   intent(in)  :: instance      !*FD the model instance

    real(dp),dimension(:,:),intent(out) :: orog          !*FD the orographic elevation (m)
    real(dp),dimension(:,:),intent(out) :: albedo        !*FD the albedo of ice/snow
    real(dp),dimension(:,:),intent(out) :: ice_frac      !*FD The fraction covered by ice
    real(dp),dimension(:,:),intent(out) :: veg_frac      !*FD The fraction of exposed vegetation
    real(dp),dimension(:,:),intent(out) :: snowice_frac  !*FD The fraction of snow-covered ice
    real(dp),dimension(:,:),intent(out) :: snowveg_frac  !*FD The fraction of snow-covered vegetation
    real(dp),dimension(:,:),intent(out) :: snow_depth    !*FD The mean snow-depth over those 
    !*FD parts covered in snow (m w.e.)

    ! Internal variables -------------------------------------------------------------------------------

    real(dp),dimension(:,:),pointer :: temp => null()

    ! --------------------------------------------------------------------------------------------------
    ! Orography

    call local_to_global_avg(instance%ups_orog, &
                             instance%model%geometry%usrf, &
                             orog,    &
                             instance%out_mask)
    orog=thk0*orog

    call coordsystem_allocate(instance%lgrid,temp)

    ! Ice-no-snow fraction
    where (instance%mbal_accum%snowd == 0.d0 .and. instance%model%geometry%thck > 0.d0)
       temp = 1.d0
    elsewhere
       temp = 0.d0
    endwhere

    call local_to_global_avg(instance%ups, &
                             temp, &
                             ice_frac,    &
                             instance%out_mask)

    ! Ice-with-snow fraction
    where (instance%mbal_accum%snowd > 0.d0 .and. instance%model%geometry%thck > 0.d0)
       temp = 1.d0
    elsewhere
       temp = 0.d0
    endwhere
    call local_to_global_avg(instance%ups, &
                             temp, &
                             snowice_frac,    &
                             instance%out_mask)

    ! Veg-with-snow fraction (if ice <10m thick)
    where (instance%mbal_accum%snowd > 0.d0 .and. instance%model%geometry%thck <= (10.d0/thk0))
       temp = 1.d0
    elsewhere
       temp = 0.d0
    endwhere
    call local_to_global_avg(instance%ups, &
                             temp, &
                             snowveg_frac,    &
                             instance%out_mask)

    ! Remainder is veg only
    veg_frac = 1.d0 - ice_frac - snowice_frac - snowveg_frac

    ! Snow depth

    call local_to_global_avg(instance%ups, &
                             instance%mbal_accum%snowd, &
                             snow_depth,    &
                             instance%out_mask)

    ! Albedo

    where ((ice_frac+snowice_frac) > 0.d0)
       albedo = instance%ice_albedo
    elsewhere
       albedo = 0.d0
    endwhere

    deallocate(temp)
    temp => null()

  end subroutine glint_upscaling

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_upscaling_gcm(instance,    nec,      &
                                 nxl,         nyl,      &
                                 nxg,         nyg,      &
                                 box_areas,             &
                                 gfrac,       gtopo,    &
                                 grofi,       grofl,    &
                                 ghflx,                 &
                                 init_call)

    ! Upscale fields from the local grid to the global grid (with multiple elevation classes).
    ! Output fields are only valid on the main task.
    ! The upscaled fields are passed to the GCM land surface model, which has the option
    !  of updating the fractional area and surface elevation of glaciated gridcells.
    ! If instance%zero_gcm_fluxes is true, then the upscaled versions of grofi, grofl and
    ! ghflx are zeroed out

    use glimmer_paramets, only: thk0, GLC_DEBUG
    use glimmer_log
    use parallel, only: tasks, main_task

!WHL - debug
    use glimmer_paramets, only: tim0

    ! Arguments ----------------------------------------------------------------------------
 
    type(glint_instance), intent(inout) :: instance      ! the model instance
    integer,              intent(in)    :: nec           ! number of elevation classes
    integer,              intent(in)    :: nxl,nyl       ! local grid dimensions
    integer,              intent(in)    :: nxg,nyg       ! local grid dimensions    
    real(dp),dimension(nxg,nyg),      intent(in)  :: box_areas ! global grid cell areas (m^2)
    real(dp),dimension(nxg,nyg,0:nec),intent(out) :: gfrac   ! ice/land-covered fraction [0,1]
    real(dp),dimension(nxg,nyg,0:nec),intent(out) :: gtopo   ! surface elevation (m)
    real(dp),dimension(nxg,nyg,0:nec),intent(out) :: ghflx   ! heat flux (m)
    real(dp),dimension(nxg,nyg),      intent(out) :: grofi   ! ice runoff (calving) flux (kg/m^2/s)
    real(dp),dimension(nxg,nyg),      intent(out) :: grofl   ! liquid runoff (basal melt) flux (kg/m^2/s)
 
    logical, intent(in), optional :: init_call   ! true if called during initialization

    ! Internal variables ----------------------------------------------------------------------

    integer :: i, j, n      ! indices
    integer :: il, jl, ig, jg

    character(len=100) :: message
    real(dp) :: dew, dns    ! gridcell dimensions
    real(dp) :: usrf, thck, topg  ! surface elevation, ice thickness, bed elevation (m)

    real(dp), dimension(nxl,nyl) ::  &
       area_l,             &! local gridcell area
       area_rofi_l,        &! area*rofi on local grid
       area_rofl_l          ! area*rofl on local grid       

    real(dp), dimension(nxg,nyg) ::  &
       area_g               ! global gridcell area (including ocean)

    real(dp), dimension(nxl,nyl,0:nec) ::  &
       area_frac_l,                        &! area*frac per elevation class on local grid
       area_topo_l,                        &! area*topo per elevation class on local grid
       area_hflx_l                          ! area*hflx per elevation class on local grid
       
    integer, dimension(nxl,nyl,0:nec) ::   &   
       area_mask_l                          !binary mask, defined at all elevation classes,
                                            !that defines whether some ice elevation is present.
                                            !For the 0-indexed bed information, this is 1 
                                            !for all land points (ice or ice-free).
    
    !TODO - Pass in topomax as an argument instead of hardwiring it here
    real(dp), dimension(0:nec) :: topomax   ! upper elevation limit of each class

    logical :: first_call   ! if calling the first time, then do not average the accumulated fluxes
                            ! use values from restart file if available

    first_call = .false.
    if (present(init_call)) then
       if (init_call) first_call = .true.
    endif 

    dew = get_dew(instance%model)
    dns = get_dns(instance%model)

    ! Given the value of nec, specify the upper and lower elevation boundaries of each class.
    ! TODO: These must be consistent with the values in the GCM.  Better to pass as an argument.

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
     
    gfrac(:,:,0:nec) = 0.d0
    gtopo(:,:,0:nec) = 0.d0
    ghflx(:,:,0:nec) = 0.d0
    grofi(:,:)   = 0.d0
    grofl(:,:)   = 0.d0
    area_l(:,:) = 0.d0
    area_g(:,:) = 0.d0
    area_frac_l(:,:,0:nec) = 0.d0
    area_topo_l(:,:,0:nec) = 0.d0
    area_hflx_l(:,:,0:nec) = 0.d0
    area_mask_l(:,:,0:nec)   = 0
    area_rofi_l(:,:) = 0.d0
    area_rofl_l(:,:) = 0.d0
    
    ! Compute time-average fluxes (unless called during initialization)
    if (first_call) then
       ! do nothing; use values from restart file if restarting
    else
       if (instance%av_count_output > 0) then
          instance%rofi_tavg(:,:) = instance%rofi_tavg(:,:) / real(instance%av_count_output,dp)
          instance%rofl_tavg(:,:) = instance%rofl_tavg(:,:) / real(instance%av_count_output,dp)
          instance%hflx_tavg(:,:) = instance%hflx_tavg(:,:) / real(instance%av_count_output,dp)
       else
          instance%rofi_tavg(:,:) = 0.d0
          instance%rofl_tavg(:,:) = 0.d0
          instance%hflx_tavg(:,:) = 0.d0
       endif
    endif

    ! Reset the logical variable for averaging output

    instance%new_tavg_output = .true.

    !To account for the potential for variable grid sizes, multiply all values by the area of each grid
    !cell, upscale, then divide by the upscaled area afterwards.  Note that this is redundant for the current, constant-
    !sized grid.  Also in this loop: bin ice sheet grid cell values by elevation.

    do j = 1, nyl
      do i = 1, nxl

         usrf = thk0 * instance%model%geometry%usrf(i,j)
         thck = thk0 * instance%model%geometry%thck(i,j)
         topg = thk0 * instance%model%geometry%topg(i,j)
         
         if (usrf > 0.d0) then!if not at sea level (assume a land point)...
            if (thck <= min_thck) then !and is not ice covered...
               area_frac_l(i,j,0) = dew*dns !accumulate bare land area fraction in 0-indexed cell.
            else !and is ice covered...
               do n = 1, nec       
                  if (usrf >= topomax(n-1) .and. usrf < topomax(n)) then    !local cell is in elev class n
                    area_frac_l(i,j,n) = dew*dns                            !accumulate ice area fraction
                    area_topo_l(i,j,n) = dew*dns * usrf                     !accumulate topography
                    ! Setting hflx to 0 for now to avoid giving the impression that it's
                    ! being used, since currently CLM doesn't handle it
!                    area_hflx_l(i,j,n) = dew*dns * instance%hflx_tavg(i,j) !accumulate heat flux
                    area_hflx_l(i,j,n) = 0.d0
                    area_mask_l(i,j,n) = 1                                  !accumulate ice area
                    exit
                  endif
               enddo   ! nec       
            endif  
            !for upscaled bed topographies and heat fluxes, include values under ice.
            ! The rationale for topography is that it gets hard to analyze / explain
            ! results if bare land topography changes whenever ice expands or retreats - 
            ! so we're using a formulation that results in bare land topography being
            ! constant in time (as long as topg is constant - e.g., neglecting isostasy).
            area_topo_l(i,j,0) = dew*dns * topg
            ! Setting hflx to 0 for now to avoid giving the impression that it's
            ! being used, since currently CLM doesn't handle it
!           area_hflx_l(i,j,0) = dew*dns * instance%hflx_tavg(i,j)
            area_hflx_l(i,j,0) = 0.d0
            area_mask_l(i,j,0) = 1
            area_l(i,j) = dew*dns            
         else
            !grid cell is presumed ocean point.    
            ! Ocean points are unglaciated, so we treat this as a 0-area point. We
            ! eventually want to handle this by keeping CLM consistent with CISm in terms
            ! of its breakdown into land vs "ocean" (e.g., wetland in CLM). In that case,
            ! if CISM says a point is ocean, then it would tell CLM that that point is
            ! ocean, and so CLM wouldn't try to generate SMB there.

            area_l(i,j) = 0.d0
         endif

         ! Runoff fluxes (note, these fluxes can be nonzero for cells with no ice in current timestep)
         area_rofi_l(i,j) = dew*dns * instance%rofi_tavg(i,j)
         area_rofl_l(i,j) = dew*dns * instance%rofl_tavg(i,j)
         
      enddo
    enddo

    ! Map the area-weighted local values to the global grid.  Note that there are three ways to do this.
    ! Ensure that your choice makes physical sense. 
    ! 1) Sum up all values of 'children' ice sheet grid points, for each 'parent' climate grid cell.  
    !    Example usage: get total volume of ice associated with each parent climate grid cell, 
    ! 2) Average all values of 'children' ice sheet grid points, for each 'parent' climate grid cell.
    !    Example: get average surface elevation of children ice sheet points, for each parent climate grid cell.
    ! 3) Minimum of all values of 'children' ice sheet grid points, for each 'parent' climate grid cell.

    ! Total area of non-ocean ice cells within global grid cell
    call local_to_global_sum(instance%ups, &
                             area_l,       &
                             area_g)

    !Total solid ice flux                            
    call local_to_global_sum(instance%ups,       &
                             area_rofi_l(:,:),   &
                             grofi(:,:))

    !Total basal runoff
    call local_to_global_sum(instance%ups,       &
                             area_rofl_l(:,:),   &
                             grofl(:,:))
                     
    ! Loop over elevation classes to generate global-grid-based gfrac, gtopo, ghflx fields.
    ! Note the values in the zero index refer to bare-land values.
    do n = 0, nec

       call local_to_global_sum(instance%ups,         &
                                area_frac_l(:,:,n),   &
                                gfrac(:,:,n))

       if (n==0) then !for bare land topography, use minimum elevation of child grid cell, as the value for the parent grid cell.
       
         call local_to_global_min(instance%ups,         &                                
                                area_topo_l(:,:,n),   &
                                gtopo(:,:,n),   &
                                area_mask_l(:,:,n))
                                
       else

         call local_to_global_avg(instance%ups,         &
                                area_topo_l(:,:,n),   &
                                gtopo(:,:,n),   &
                                area_mask_l(:,:,n))
       endif

       call local_to_global_avg(instance%ups,         &
                                area_hflx_l(:,:,n),   &
                                ghflx(:,:,n),   &
                                area_mask_l(:,:,n))

       do j = 1, nyg
          do i = 1, nxg
             if (area_g(i,j) > 0.d0) then
                gfrac(i,j,n) = gfrac(i,j,n) / area_g(i,j) !fraction of elevation class, relative to total land area 
             else
               gfrac(i,j,n) = 0.d0
             endif
             if (n==0) then !non-ice-class values: use model-derived ice-free topography/heat flux, regardless of whether ice exists.
                gtopo(i,j,n) = gtopo(i,j,n) / (dew*dns)
                gtopo(i,j,n) = max(0.d0,gtopo(i,j,n)) !Keep elevations to equal or greater than sea level.
                ghflx(i,j,n) = ghflx(i,j,n) / (dew*dns)
             else  !ice class values
               if (gfrac(i,j,n) > 0.d0) then !ice ice area exists, use model-derived topography/heat flux values
                  gtopo(i,j,n) = gtopo(i,j,n) / (dew*dns) 
                  ghflx(i,j,n) = ghflx(i,j,n) / (dew*dns)              
               else !use prescribed, idealized topography/heat flux values
                  gtopo(i,j,n) = mean_elevation_virtual(n, nec, topomax)
                  ghflx(i,j,n) = 0.d0
               endif
             endif
          enddo
       enddo

    enddo

    do j = 1, nyg
       do i = 1, nxg
          ! Find mean ice runoff from calving 
          ! Note: Here we divide by box_areas (the area of the global grid cell), which in general is not equal
          ! to area_g (the sum over area of the local grid cells associated with the global cell).
          ! We do this to ensure conservation of ice mass (=flux*area) when multiplying later by the
          ! area of the global grid cell.
          if (box_areas(i,j) > 0.d0) then
             grofi(i,j) = grofi(i,j) / box_areas(i,j)     
             grofl(i,j) = grofl(i,j) / box_areas(i,j)
          endif
       enddo
    enddo

    if (instance%zero_gcm_fluxes == ZERO_GCM_FLUXES_TRUE) then
       ghflx(:,:,0:nec) = 0.d0
       grofi(:,:)   = 0.d0
       grofl(:,:)   = 0.d0
    end if
    
  end subroutine glint_upscaling_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_accumulate_output_gcm(model,            &
                                         av_count_output,  &
                                         new_tavg_output,  &
                                         rofi_tavg,        &
                                         rofl_tavg,        &
                                         hflx_tavg)

    ! Given the calving, basal melting, and conductive heat flux fields from the dycore,
    ! accumulate contributions to the rofi, rofl, and hflx fields to be sent to the coupler.

    use glimmer_paramets, only: thk0, tim0

    use glimmer_scales, only: scale_acab  ! for testing

!WHL - debug - Set to inout if specifying the model fields for testing
    type(glide_global_type), intent(in)  :: model

    integer,  intent(inout) :: av_count_output     ! step counter 
    logical,  intent(inout) :: new_tavg_output     ! if true, start new averaging
    real(dp), dimension(:,:), intent(inout) :: rofi_tavg    ! solid ice runoff (kg m-2 s-1)
    real(dp), dimension(:,:), intent(inout) :: rofl_tavg    ! liquid runoff from basal/interior melting (kg m-2 s-1)
    real(dp), dimension(:,:), intent(inout) :: hflx_tavg    ! conductive heat flux at top surface (W m-2)

    ! things to do the first time

    if (new_tavg_output) then

       new_tavg_output = .false.
       av_count_output  = 0

       ! Initialise
       rofi_tavg(:,:) = 0.d0
       rofl_tavg(:,:) = 0.d0
       hflx_tavg(:,:) = 0.d0

    end if

    av_count_output = av_count_output + 1

    !--------------------------------------------------------------------
    ! Accumulate solid runoff (calving)
    !--------------------------------------------------------------------
                       
    ! Note on units: model%climate%calving has dimensionless ice thickness units
    !                Multiply by thk0 to convert to meters of ice
    !                Multiply by rhoi to convert to kg/m^2 water equiv.
    !                Divide by (dt*tim0) to convert to kg/m^2/s

    ! Convert to kg/m^2/s
    rofi_tavg(:,:) = rofi_tavg(:,:)  &
                   + model%climate%calving(:,:) * thk0 * rhoi / (model%numerics%dt * tim0)

    !--------------------------------------------------------------------
    ! Accumulate liquid runoff (basal melting)
    !--------------------------------------------------------------------
    !TODO - Add internal melting for enthalpy case
                       
    ! Note on units: model%temper%bmlt has dimensionless units of ice thickness per unit time
    !                Multiply by thk0/tim0 to convert to meters ice per second
    !                Multiply by rhoi to convert to kg/m^2/s water equiv.

    ! Convert to kg/m^2/s
    rofl_tavg(:,:) = rofl_tavg(:,:)  &
                   + model%temper%bmlt(:,:) * thk0/tim0 * rhoi

    !--------------------------------------------------------------------
    ! Accumulate basal heat flux
    !--------------------------------------------------------------------

    ! Note on units: model%temper%ucondflx has units of W/m^2, positive down
    !                Flip the sign so that hflx_tavg is positive up.

    hflx_tavg(:,:) = hflx_tavg(:,:) &
                   - model%temper%ucondflx(:,:)

  end subroutine glint_accumulate_output_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function mean_elevation_virtual(ec, nec, topomax)
    
    ! For a "virtual" elevation class (that is, an elevation class that has 0 area),
    ! return the "mean" elevation of the given elevation class.

    use glimmer_log

    ! Arguments ----------------------------------------------------------------------------
    
    integer  , intent(in) :: ec          ! elevation class
    integer  , intent(in) :: nec         ! number of elevation classes
    real(dp) , intent(in) :: topomax(0:) ! upper elevation limit of each class

    ! TODO: replace this with a call to shr_assert, if/when glimmer-cism pulls in csm_share
    if (ubound(topomax, 1) /= nec) then
       call write_log('ERROR: upper bound of topomax does not match nec', &
            GM_FATAL, __FILE__, __LINE__)
    end if

    if (ec < nec) then
       mean_elevation_virtual = (topomax(ec-1) + topomax(ec))/2.0_dp
    else if (ec == nec) then
       ! In the top elevation class; in this case, assignment of a "mean" elevation is
       ! somewhat arbitrary
       
       if (nec > 1) then
          mean_elevation_virtual = 2.0_dp * topomax(ec-1) - topomax(ec-2)
       else
          ! entirely arbitrary
          mean_elevation_virtual = 1000._dp
       end if
    else
       call write_log('ERROR: class out of bounds', &
            GM_FATAL, __FILE__, __LINE__)
    end if

  end function mean_elevation_virtual

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module glint_upscale

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
