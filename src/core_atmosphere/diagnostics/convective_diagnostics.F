! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module convective_diagnostics

    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type, MPAS_LOG_ERR, MPAS_LOG_CRIT
    use mpas_kind_types, only : RKIND
    use mpas_log, only : mpas_log_write

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag

    type (MPAS_clock_type), pointer :: clock

    public :: convective_diagnostics_setup, &
              convective_diagnostics_update, &
              convective_diagnostics_compute, &
              convective_diagnostics_reset

    private

    !
    ! For any fields where we need the min/max/mean, it is helpful to determine
    ! in the setup routine whether these fields will actually appear in any
    ! output streams; if so, we need to handle the field in the update/reset
    ! routines, but if not, we can save some computation
    !
    logical :: is_needed_updraft_helicity
    logical :: is_needed_w_max
    logical :: is_needed_lml_wsp_max
  

    contains


    !-----------------------------------------------------------------------
    !  routine convective_diagnostics_setup
    !
    !> \brief Set-up the convective diagnostics module
    !> \author Michael Duda
    !> \date   14 October 2016
    !> \details
    !>  To avoid later work in dereferencing pointers to various pools,
    !>  this routine saves pool pointers for use by
    !>  the convective_diagnostics_compute routine.
    !
    !-----------------------------------------------------------------------
    subroutine convective_diagnostics_setup(all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type, MPAS_STREAM_OUTPUT, MPAS_STREAM_INPUT, &
                                       MPAS_STREAM_INPUT_OUTPUT
        use mpas_pool_routines, only : mpas_pool_get_subpool
        use mpas_atm_diagnostics_utils, only : mpas_stream_inclusion_count

        implicit none

        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)

        clock => simulation_clock

        is_needed_updraft_helicity = .false.
        is_needed_w_max = .false.
        is_needed_lml_wsp_max = .false.

        if (mpas_stream_inclusion_count('updraft_helicity_max', direction=MPAS_STREAM_OUTPUT) > 0) then
            is_needed_updraft_helicity = .true.
        end if
        if (mpas_stream_inclusion_count('w_velocity_max', direction=MPAS_STREAM_OUTPUT) > 0) then
            is_needed_w_max = .true.
        end if
        if (mpas_stream_inclusion_count('wind_speed_level1_max', direction=MPAS_STREAM_OUTPUT) > 0) then
            is_needed_lml_wsp_max = .true.
        end if
   
    end subroutine convective_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine convective_diagnostics_update
    !
    !> \brief Updates the maximum pointwise values for convective diagnostics
    !> \author Michael Duda
    !> \date   18 October 2016
    !> \details
    !>  Updates the maximum pointwise values for updraft helicity, w velocity,
    !>  and lowest-model-level wind speed.
    !
    !-----------------------------------------------------------------------
    subroutine convective_diagnostics_update()

        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

        implicit none

!  - subroutine to compute max values of convective diagnostics over some time period
!  - called at the end of every timestep
!  - values must be zeroed out after output, and this determines the period for the max.
!
!  WCS, March 2015, for the SPC spring experiment forecasts

        integer :: iCell, k, nVertLevelsP1, i, j
        integer, pointer :: nCells, nVertLevels, nCellsSolve, nVertexDegree, nVertices
        integer, pointer :: moist_start, moist_end
        integer, dimension(:,:), pointer :: cellsOnVertex
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: cqu

        real (kind=RKIND), allocatable, dimension(:,:) :: updraft_helicity
        real (kind=RKIND), allocatable, dimension(:) :: z_agl
        real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex, vorticity, w, zgrid
        real (kind=RKIND), dimension(:), pointer :: areaCell, updraft_helicity_max, w_velocity_max
        real (kind=RKIND), dimension(:), pointer :: wind_speed_level1_max
        real (kind=RKIND), dimension(:,:), pointer :: uzonal, umeridional
        real (kind=RKIND) :: uph

        if (is_needed_updraft_helicity .or. is_needed_w_max .or. is_needed_lml_wsp_max) then

            call mpas_pool_get_dimension(mesh, 'nCells', nCells)
            call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
            call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
            call mpas_pool_get_dimension(mesh, 'vertexDegree', nVertexDegree)
            call mpas_pool_get_array(mesh, 'areaCell', areaCell)
            call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
            call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
            call mpas_pool_get_array(mesh, 'zgrid', zgrid)

            call mpas_pool_get_array(diag, 'vorticity', vorticity)
            call mpas_pool_get_array(diag, 'updraft_helicity_max', updraft_helicity_max)
            call mpas_pool_get_array(diag, 'w_velocity_max', w_velocity_max)
            call mpas_pool_get_array(diag, 'wind_speed_level1_max', wind_speed_level1_max)
            call mpas_pool_get_array(diag, 'uReconstructZonal', uzonal)
            call mpas_pool_get_array(diag, 'uReconstructMeridional', umeridional)

            call mpas_pool_get_array(state, 'w', w, 1)

            nVertLevelsP1 = nVertLevels + 1

            allocate(z_agl(nVertLevelsP1))
!
!  updraft helicity, first compute vorticity at cell center, then mulitply by w
!  *** where does cellOnVertex point if cell is outside block and block halo?  -> cell nCell+1
!
            if (is_needed_w_max) then
                do iCell=1, nCellsSolve
                  w_velocity_max(iCell) = max( w_velocity_max(iCell), maxval(w(1:nVertLevels,iCell)))
                end do
            end if

            if (is_needed_lml_wsp_max) then
                do iCell=1, nCellsSolve
                  wind_speed_level1_max(iCell) = max(wind_speed_level1_max(iCell), sqrt(uzonal(1,iCell)**2+umeridional(1,iCell)**2))
                end do
            end if

            if (is_needed_updraft_helicity) then
                allocate(updraft_helicity(nVertLevels,nCells+1))
                updraft_helicity(:,:) = 0.
                do i=1,nVertices
                  do j=1,nvertexDegree
                    iCell = cellsOnVertex(j,i)
                    updraft_helicity(1:nVertLevels,iCell) = updraft_helicity(1:nVertLevels,iCell) + kiteAreasOnVertex(j,i)*vorticity(1:nVertLevels,i)
                  end do
                end do
                do iCell=1,nCellsSolve
                  do k=1,nVertLevels
                    updraft_helicity(k,iCell) =  max(0.,0.5*(w(k,iCell)+w(k+1,iCell)))  &                  
                                                           * max(0.,updraft_helicity(k,iCell)/areaCell(iCell))                     
                  end do
                end do
!
!  compute diagnostics
!
                do iCell=1,nCellsSolve
    
                  ! compute above ground level (AGL) heights
                  z_agl(1:nVertLevelsP1) = zgrid(1:nVertLevelsP1,iCell) - zgrid(1,iCell)
                  uph = integral_zstaggered(updraft_helicity(1:nVertLevels,iCell),z_agl,2000.0_RKIND,5000.0_RKIND,nVertLevels,nVertLevelsP1)
                  updraft_helicity_max(iCell) = max( updraft_helicity_max(iCell),uph)
        
                end do
                deallocate(updraft_helicity)
            end if
    
            deallocate(z_agl)
        end if
   
    end subroutine convective_diagnostics_update


    !-----------------------------------------------------------------------
    !  routine convective_diagnostics_compute
    !
    !> \brief Computes convective diagnostic
    !> \author Michael Duda
    !> \date   14 October 2016
    !> \details
    !>  This routine computes several diagnostics used in Spring Experiment
    !>  runs and originally added by WCS in March 2015.
    !>  The following fields are computed by this routine:
    !>    cape
    !>    cin
    !>    lcl
    !>    lfc
    !>    srh_0_1km
    !>    srh_0_3km
    !>    uzonal_surface
    !>    uzonal_1km
    !>    uzonal_6km
    !>    umeridional_surface
    !>    umeridional_1km
    !>    umeridional_6km
    !>    temperature_surface
    !>    dewpoint_surface
    !
    !-----------------------------------------------------------------------
    subroutine convective_diagnostics_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
        use mpas_constants, only : rvord
        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

        implicit none

        integer :: iCell, k
        integer, pointer :: nCells, nCellsSolve, nVertLevels, nVertLevelsP1

        ! Fields that are computed in this routine
        real (kind=RKIND), dimension(:), pointer :: cape
        real (kind=RKIND), dimension(:), pointer :: cin
        real (kind=RKIND), dimension(:), pointer :: lcl
        real (kind=RKIND), dimension(:), pointer :: lfc
        real (kind=RKIND), dimension(:), pointer :: srh_0_1km
        real (kind=RKIND), dimension(:), pointer :: srh_0_3km
        real (kind=RKIND), dimension(:), pointer :: uzonal_surface
        real (kind=RKIND), dimension(:), pointer :: uzonal_1km
        real (kind=RKIND), dimension(:), pointer :: uzonal_6km
        real (kind=RKIND), dimension(:), pointer :: umeridional_surface
        real (kind=RKIND), dimension(:), pointer :: umeridional_1km
        real (kind=RKIND), dimension(:), pointer :: umeridional_6km
        real (kind=RKIND), dimension(:), pointer :: temperature_surface
        real (kind=RKIND), dimension(:), pointer :: dewpoint_surface

        ! Other fields used in the computation of convective diagnostics 
        ! defined above
        real (kind=RKIND), dimension(:,:), pointer :: height
        real (kind=RKIND), dimension(:,:), pointer :: uzonal
        real (kind=RKIND), dimension(:,:), pointer :: umeridional
        real (kind=RKIND), dimension(:,:), pointer :: relhum
        real (kind=RKIND), dimension(:,:), pointer :: exner
        real (kind=RKIND), dimension(:,:), pointer :: theta_m
        real (kind=RKIND), dimension(:,:), pointer :: pressure_p
        real (kind=RKIND), dimension(:,:), pointer :: pressure_base
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        integer, pointer :: index_qv

        real (kind=RKIND), parameter :: dev_motion = 7.5, z_bunker_bot = 0., z_bunker_top = 6000.
        real (kind=RKIND)            :: u_storm, v_storm, u_srh_bot, v_srh_bot, u_srh_top, v_srh_top
        real (kind=RKIND)            :: u_mean, v_mean, u_shear, v_shear, shear_magnitude
        real (kind=RKIND)            :: b_term, cape_out, cin_out
        real (kind=RKIND), dimension(:), allocatable :: dudz, dvdz, zp
        real (kind=RKIND), dimension(:), allocatable :: zrel, srh
        real (kind=RKIND), dimension(:), allocatable :: p_in, t_in, td_in

        real (kind=RKIND), dimension(:,:), allocatable :: temperature, dewpoint

        real (kind=RKIND) :: evp

        logical :: need_cape, need_cin, need_lcl, need_lfc, need_srh_01km, need_srh_03km, need_uzonal_sfc, need_uzonal_1km, &
                   need_uzonal_6km, need_umerid_sfc, need_umerid_1km, need_umerid_6km, need_tsfc, need_tdsfc
        logical :: need_any_diags

        need_any_diags = .false.

        need_cape = MPAS_field_will_be_written('cape')
        need_any_diags = need_any_diags .or. need_cape
        need_cin = MPAS_field_will_be_written('cin')
        need_any_diags = need_any_diags .or. need_cin
        need_lcl = MPAS_field_will_be_written('lcl')
        need_any_diags = need_any_diags .or. need_lcl
        need_lfc = MPAS_field_will_be_written('lfc')
        need_any_diags = need_any_diags .or. need_lfc
        need_srh_01km = MPAS_field_will_be_written('srh_0_1km')
        need_any_diags = need_any_diags .or. need_srh_01km
        need_srh_03km = MPAS_field_will_be_written('srh_0_3km')
        need_any_diags = need_any_diags .or. need_srh_03km
        need_uzonal_sfc = MPAS_field_will_be_written('uzonal_surface')
        need_any_diags = need_any_diags .or. need_uzonal_sfc
        need_uzonal_1km = MPAS_field_will_be_written('uzonal_1km')
        need_any_diags = need_any_diags .or. need_uzonal_1km
        need_uzonal_6km = MPAS_field_will_be_written('uzonal_6km')
        need_any_diags = need_any_diags .or. need_uzonal_6km
        need_umerid_sfc = MPAS_field_will_be_written('umeridional_surface')
        need_any_diags = need_any_diags .or. need_umerid_sfc
        need_umerid_1km = MPAS_field_will_be_written('umeridional_1km')
        need_any_diags = need_any_diags .or. need_umerid_1km
        need_umerid_6km = MPAS_field_will_be_written('umeridional_6km')
        need_any_diags = need_any_diags .or. need_umerid_6km
        need_tsfc = MPAS_field_will_be_written('temperature_surface')
        need_any_diags = need_any_diags .or. need_tsfc
        need_tdsfc = MPAS_field_will_be_written('dewpoint_surface')
        need_any_diags = need_any_diags .or. need_tdsfc

        if (need_any_diags) then

            call mpas_pool_get_dimension(mesh, 'nCells', nCells)
            call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
            call mpas_pool_get_dimension(mesh, 'nVertLevelsP1', nVertLevelsP1)

            call mpas_pool_get_array(diag, 'cape', cape)
            call mpas_pool_get_array(diag, 'cin',  cin)
            call mpas_pool_get_array(diag, 'lcl',  lcl)
            call mpas_pool_get_array(diag, 'lfc',  lfc)
            call mpas_pool_get_array(diag, 'srh_0_1km', srh_0_1km)
            call mpas_pool_get_array(diag, 'srh_0_3km', srh_0_3km)
            call mpas_pool_get_array(diag, 'uzonal_surface', uzonal_surface)
            call mpas_pool_get_array(diag, 'uzonal_1km', uzonal_1km)
            call mpas_pool_get_array(diag, 'uzonal_6km', uzonal_6km)
            call mpas_pool_get_array(diag, 'umeridional_surface', umeridional_surface)
            call mpas_pool_get_array(diag, 'umeridional_1km', umeridional_1km)
            call mpas_pool_get_array(diag, 'umeridional_6km', umeridional_6km)
            call mpas_pool_get_array(diag, 'temperature_surface', temperature_surface)
            call mpas_pool_get_array(diag, 'dewpoint_surface', dewpoint_surface)

            call mpas_pool_get_array(mesh, 'zgrid', height)
            call mpas_pool_get_array(diag, 'uReconstructMeridional', umeridional)
            call mpas_pool_get_array(diag, 'uReconstructZonal', uzonal)
            call mpas_pool_get_array(state, 'theta_m', theta_m, 1)
            call mpas_pool_get_array(state, 'scalars', scalars, 1)
            call mpas_pool_get_array(diag, 'exner', exner)
            call mpas_pool_get_array(diag, 'relhum', relhum)
            call mpas_pool_get_array(diag, 'pressure_base', pressure_base)
            call mpas_pool_get_array(diag, 'pressure_p', pressure_p)

            call mpas_pool_get_dimension(state, 'index_qv', index_qv)

            allocate(temperature(nVertLevels,nCells))
            allocate(dewpoint(nVertLevels,nCells))

            allocate(dudz(nVertLevels))
            allocate(dvdz(nVertLevels))
            allocate(zp(nVertLevels))
            allocate(zrel(nVertLevels+1))
            allocate(srh(nVertLevels+1))

            allocate(p_in(nVertLevels))
            allocate(t_in(nVertLevels))
            allocate(td_in(nVertLevels))

            do iCell = 1,nCells
            do k = 1,nVertLevels
                temperature(k,iCell) = (theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))*exner(k,iCell)

                ! Vapor pressure
                evp = 0.01_RKIND * (pressure_base(k,iCell) + pressure_p(k,iCell)) &
                                 * scalars(index_qv,k,iCell) / (scalars(index_qv,k,iCell) + 0.622_RKIND)
                evp = max(evp, 1.0e-8_RKIND)

                ! Dewpoint temperature following Bolton (1980)
                dewpoint(k,iCell) = (243.5_RKIND * log(evp/6.112_RKIND)) / (17.67_RKIND - log(evp/6.112_RKIND))
                dewpoint(k,iCell) = dewpoint(k,iCell) + 273.15
            enddo
            enddo


            ! first the shear values.  We will use lowest model level velocity for surface velocity
            do iCell=1,nCellsSolve

              zp(1:nVertLevels) = 0.5*(height(1:nVertLevels,iCell)+height(2:nVertlevels+1,iCell)) - height(1,iCell)
              zrel(1:nVertLevels+1) = height(1:nVertLevels+1,iCell) - height(1,iCell)

              uzonal_surface(iCell) = uzonal(1,iCell)
              umeridional_surface(iCell) = umeridional(1,iCell)
              temperature_surface(iCell) = temperature(1,iCell)
              dewpoint_surface(iCell) = dewpoint(1,iCell)
              if (need_uzonal_1km) then
                  uzonal_1km(iCell)      = column_height_value(uzonal(1:nVertLevels,iCell),      zp, 1000.0_RKIND, nVertLevels)
              end if
              if (need_umerid_1km) then
                  umeridional_1km(iCell) = column_height_value(umeridional(1:nVertLevels,iCell), zp, 1000.0_RKIND, nVertLevels)
              end if
              if (need_uzonal_6km) then
                  uzonal_6km(iCell)      = column_height_value(uzonal(1:nVertLevels,iCell),      zp, 6000.0_RKIND, nVertLevels)
              end if
              if (need_umerid_6km) then
                  umeridional_6km(iCell) = column_height_value(umeridional(1:nVertLevels,iCell), zp, 6000.0_RKIND, nVertLevels)
              end if

           !  storm-relative helicity
           !  first, calculate storm motion, using Bunkers formula for right-moving storms

               if (need_srh_01km .or. need_srh_03km) then
                   u_srh_bot = uzonal(1,iCell)
                   v_srh_bot = umeridional(1,iCell)
                   if(z_bunker_bot .gt. zp(1)) then
                     u_srh_bot = column_height_value( uzonal(1:nVertLevels,iCell), zp, z_bunker_bot, nVertLevels)
                     v_srh_bot = column_height_value( umeridional(1:nVertLevels,iCell), zp, z_bunker_bot, nVertLevels)
                   end if
                   u_srh_top = column_height_value( uzonal(1:nVertLevels,iCell), zp, z_bunker_top, nVertLevels)
                   v_srh_top = column_height_value( umeridional(1:nVertLevels,iCell), zp, z_bunker_top, nVertLevels)
                   u_shear = u_srh_top - u_srh_bot
                   v_shear = v_srh_top - v_srh_bot
                   u_mean = integral_zstaggered(uzonal(1:nVertLevels,iCell),zrel,z_bunker_bot,z_bunker_top,nVertLevels,nVertLevelsP1)/(z_bunker_top-z_bunker_bot)
                   v_mean = integral_zstaggered(umeridional(1:nVertLevels,iCell),zrel,z_bunker_bot,z_bunker_top,nVertLevels,nVertLevelsP1)/(z_bunker_top-z_bunker_bot)
                   shear_magnitude = max(0.0001,sqrt(u_shear**2 + v_shear**2))
                   u_storm = u_mean + dev_motion * v_shear/shear_magnitude
                   v_storm = v_mean - dev_motion * u_shear/shear_magnitude

                   !  calculate horizontal vorticity

                   do k=2, nVertLevels-1
                       dudz(k) = (uzonal(k,iCell)     -uzonal(k-1,iCell)     )/(0.5*(height(k+1,iCell)-height(k-1,iCell)))      
                       dvdz(k) = (umeridional(k,iCell)-umeridional(k-1,iCell))/(0.5*(height(k+1,iCell)-height(k-1,iCell)))
                   end do
                   dudz(1) = dudz(2)
                   dvdz(1) = dvdz(2)
                   dudz(nVertLevels) = dudz(nVertLevels-1)
                   dvdz(nVertLevels) = dvdz(nVertLevels-1)

                   do k=2,nVertLevels
                     srh(k) = - (0.5*(uzonal(k,iCell)      + uzonal(k-1,iCell)      )-u_storm)*dvdz(k) &
                              + (0.5*(umeridional(k,iCell) + umeridional(k-1,iCell) )-v_storm)*dudz(k)
                   end do
                   srh(1) = - (uzonal(1,iCell)     - u_storm)*dvdz(1) &
                            + (umeridional(1,iCell)- v_storm)*dudz(1)
                   srh(nVertLevelsP1) = srh(nVertLevels)

                   do k=1, nVertLevels+1
                     srh(k) = max(0.,srh(k))  ! discounting negative SRH
                   end do

                   if (need_srh_01km) then
                       srh_0_1km(iCell) = integral_zpoint(srh, zrel, 0.0_RKIND, 1000.0_RKIND, nVertLevelsP1)
                   end if
                   if (need_srh_03km) then
                       srh_0_3km(iCell) = integral_zpoint(srh, zrel, 0.0_RKIND, 3000.0_RKIND, nVertLevelsP1)
                   end if
               end if

            end do

            !  calculate cape and cin
            if (need_cape .or. need_cin) then
                do iCell=1, nCellsSolve
                    p_in(1:nVertLevels) = (pressure_p(1:nVertLevels,iCell) + pressure_base(1:nVertLevels,iCell)) / 100.0_RKIND
                    t_in(1:nVertLevels) = temperature(1:nVertLevels,iCell) - 273.15_RKIND
                    td_in(1:nVertLevels) = dewpoint(1:nVertLevels,iCell) - 273.15_RKIND

                 !   do k=1,nVertLevels
                 !     relhum(k,iCell) = max(1.e-08,min(1.,relhum(k,iCell)))
                 !     td_in(k) = 243.04*(log(relhum(k,iCell))+((17.625*t_in(k))/(243.04+t_in(k)))) & 
                 !                      /(17.625-log(relhum(k,iCell))-((17.625*t_in(k))/(243.04+t_in(k)))) 
                 !   end do

                    call getcape( nVertLevels, p_in, t_in, td_in, cape_out, cin_out )

                    cape(iCell) = cape_out
                    cin(iCell) = cin_out

                end do
            end if

            deallocate(temperature)
            deallocate(dewpoint)

            deallocate(dudz)
            deallocate(dvdz)
            deallocate(zp)
            deallocate(zrel)
            deallocate(srh)
            deallocate(p_in)
            deallocate(t_in)
            deallocate(td_in)
        end if

    end subroutine convective_diagnostics_compute


    !-----------------------------------------------------------------------
    !  routine convective_diagnostics_reset
    !
    !> \brief Reset maximum values for convective diagnostics
    !> \author 
    !> \date   
    !> \details
    !>  Resets the maximum pointwise values for updraft helicity, w velocity,
    !>  and lowest-model-level wind speed.
    !
    !-----------------------------------------------------------------------
    subroutine convective_diagnostics_reset()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
        use mpas_pool_routines, only : mpas_pool_get_array

        implicit none

        real (kind=RKIND), dimension(:), pointer :: updraft_helicity_max
        real (kind=RKIND), dimension(:), pointer :: w_velocity_max
        real (kind=RKIND), dimension(:), pointer :: wind_speed_level1_max

        if (MPAS_field_will_be_written('updraft_helicity_max')) then
            call mpas_pool_get_array(diag, 'updraft_helicity_max', updraft_helicity_max)
            updraft_helicity_max(:) = 0.
        end if
        if (MPAS_field_will_be_written('w_velocity_max')) then
            call mpas_pool_get_array(diag, 'w_velocity_max', w_velocity_max)
            w_velocity_max(:) = 0.
        end if
        if (MPAS_field_will_be_written('wind_speed_level1_max')) then
            call mpas_pool_get_array(diag, 'wind_speed_level1_max', wind_speed_level1_max)
            wind_speed_level1_max(:) = 0.
        end if
   
    end subroutine convective_diagnostics_reset


      real (kind=RKIND) function column_height_value( column_values, z, z_interp, n )
      implicit none
      integer n
      real (kind=RKIND) :: column_values(n), z(n), z_interp, wz, wzp1
      integer :: kz, k
!  we assume height increases monotonically with n
      kz = 1
      do k=1,n
        if(z(k) <= z_interp) kz = k
      end do
      kz = min(kz,n-1)
      
      wz = (z(kz+1)-z_interp)/(z(kz+1)-z(kz))
      wzp1 = 1. - wz
      column_height_value = wz*column_values(kz) + wzp1*column_values(kz+1)

      end function column_height_value

!---------------------------

      real (kind=RKIND) function integral_zstaggered( column_values, z, zbot, ztop, n, np1 )
      implicit none
      integer n, np1
      real (kind=RKIND) :: column_values(n), z(np1), zbot, ztop
      real (kind=RKIND) :: zb, zt

      integer :: k

!  integral from z_bot to z_top, assume cell-average values (first-order integration)
!  z increases monotonically

      integral_zstaggered = 0.
      do k=1,n
        zb = max(z(k), zbot)
        zt = min(z(k+1), ztop)
        integral_zstaggered = integral_zstaggered + column_values(k)*max(0.,(zt-zb))
      end do
      end function integral_zstaggered

!---------------------------------

      real (kind=RKIND) function integral_zpoint( column_values, z, zbot, ztop, n )
      implicit none
      integer n
      real(kind=RKIND) :: column_values(n), z(n), zbot, ztop
      real(kind=RKIND) :: zb, zt, dz, zr_midpoint, midpoint_value

      integer :: k

!  integral from z_bot to z_top, assume point values (second-order integration)
!  z increases monotonically

      integral_zpoint = 0.
      do k=1,n-1
        zb = max(z(k), zbot)
        zt = min(z(k+1), ztop)
        dz = max(0.,zt-zb)
        zr_midpoint = (0.5*(zt+zb) - z(k))/(z(k+1)-z(k))
        midpoint_value = column_values(k) + (column_values(k+1)-column_values(k))*zr_midpoint
        integral_zpoint = integral_zpoint + dz*midpoint_value
      end do
      end function integral_zpoint

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    subroutine getcape( nk , p_in , t_in , td_in, cape , cin )
    implicit none

    integer, intent(in) :: nk
    real (kind=RKIND), dimension(nk), intent(in) :: p_in,t_in,td_in
    real (kind=RKIND), intent(out) :: cape,cin

!-----------------------------------------------------------------------
!
!  getcape - a fortran90 subroutine to calculate Convective Available
!            Potential Energy (CAPE) from a sounding.
!
!  Version 1.02                           Last modified:  10 October 2008
!
!  Author:  George H. Bryan
!           Mesoscale and Microscale Meteorology Division
!           National Center for Atmospheric Research
!           Boulder, Colorado, USA
!           gbryan@ucar.edu
!
!  Disclaimer:  This code is made available WITHOUT WARRANTY.
!
!  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
!               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
!
!-----------------------------------------------------------------------
!
!  Input:     nk - number of levels in the sounding (integer)
!
!           p_in - one-dimensional array of pressure (mb) (real)
!
!           t_in - one-dimensional array of temperature (C) (real)
!
!          td_in - one-dimensional array of dewpoint temperature (C) (real)
!
!  Output:  cape - Convective Available Potential Energy (J/kg) (real)
!
!            cin - Convective Inhibition (J/kg) (real)
!
!-----------------------------------------------------------------------
!  User options:

    real (kind=RKIND), parameter :: pinc = 100.0   ! Pressure increment (Pa)
                                      ! (smaller number yields more accurate
                                      !  results,larger number makes code 
                                      !  go faster)

    integer, parameter :: source = 2    ! Source parcel:
                                        ! 1 = surface
                                        ! 2 = most unstable (max theta-e)
                                        ! 3 = mixed-layer (specify ml_depth)

    real (kind=RKIND), parameter :: ml_depth =  200.0  ! depth (m) of mixed layer 
                                          ! for source=3

    integer, parameter :: adiabat = 1   ! Formulation of moist adiabat:
                                        ! 1 = pseudoadiabatic, liquid only
                                        ! 2 = reversible, liquid only
                                        ! 3 = pseudoadiabatic, with ice
                                        ! 4 = reversible, with ice

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!            No need to modify anything below here:
!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    logical :: doit,ice,cloud,not_converged
    integer :: k,kmax,n,nloop,i,orec
    real (kind=RKIND), dimension(nk) :: p,t,td,pi,q,th,thv,z,pt,pb,pc,pn,ptv

    real (kind=RKIND) :: the,maxthe,parea,narea,lfc
    real (kind=RKIND) :: th1,p1,t1,qv1,ql1,qi1,b1,pi1,thv1,qt,dp,dz,ps,frac
    real (kind=RKIND) :: th2,p2,t2,qv2,ql2,qi2,b2,pi2,thv2
    real (kind=RKIND) :: thlast,fliq,fice,tbar,qvbar,qlbar,qibar,lhv,lhs,lhf,rm,cpm
    real*8 :: avgth,avgqv
!    real (kind=RKIND) :: getqvs,getqvi,getthe

!-----------------------------------------------------------------------

    real (kind=RKIND), parameter :: g     = 9.81
    real (kind=RKIND), parameter :: p00   = 100000.0
    real (kind=RKIND), parameter :: cp    = 1005.7
    real (kind=RKIND), parameter :: rd    = 287.04
    real (kind=RKIND), parameter :: rv    = 461.5
    real (kind=RKIND), parameter :: xlv   = 2501000.0
    real (kind=RKIND), parameter :: xls   = 2836017.0
    real (kind=RKIND), parameter :: t0    = 273.15
    real (kind=RKIND), parameter :: cpv   = 1875.0
    real (kind=RKIND), parameter :: cpl   = 4190.0
    real (kind=RKIND), parameter :: cpi   = 2118.636
    real (kind=RKIND), parameter :: lv1   = xlv+(cpl-cpv)*t0
    real (kind=RKIND), parameter :: lv2   = cpl-cpv
    real (kind=RKIND), parameter :: ls1   = xls+(cpi-cpv)*t0
    real (kind=RKIND), parameter :: ls2   = cpi-cpv

    real (kind=RKIND), parameter :: rp00  = 1.0/p00
    real (kind=RKIND), parameter :: eps   = rd/rv
    real (kind=RKIND), parameter :: reps  = rv/rd
    real (kind=RKIND), parameter :: rddcp = rd/cp
    real (kind=RKIND), parameter :: cpdrd = cp/rd
    real (kind=RKIND), parameter :: cpdg  = cp/g

    real (kind=RKIND), parameter :: converge = 0.0002

    integer, parameter :: debug_level =   0

!-----------------------------------------------------------------------

!---- convert p,t,td to mks units; get pi,q,th,thv ----!

    do k=1,nk
        p(k) = 100.0*p_in(k)
        t(k) = 273.15+t_in(k)
       td(k) = 273.15+td_in(k)
       pi(k) = (p(k)*rp00)**rddcp
        q(k) = getqvs(p(k),td(k))
       th(k) = t(k)/pi(k)
      thv(k) = th(k)*(1.0+reps*q(k))/(1.0+q(k))
    enddo

!---- get height using the hydrostatic equation ----!

    z(1) = 0.0
    do k=2,nk
      dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))
      z(k) = z(k-1) + dz
    enddo

!---- find source parcel ----!

  IF(source.eq.1)THEN
    ! use surface parcel
    kmax = 1

  ELSEIF(source.eq.2)THEN
    ! use most unstable parcel (max theta-e)

    IF(p(1).lt.50000.0)THEN
      ! first report is above 500 mb ... just use the first level reported
      kmax = 1
      maxthe = getthe(p(1),t(1),td(1),q(1))
    ELSE
      ! find max thetae below 500 mb
      maxthe = 0.0
      do k=1,nk
        if(p(k).ge.50000.0)then
          the = getthe(p(k),t(k),td(k),q(k))
          if( the.gt.maxthe )then
            maxthe = the
            kmax = k
          endif
        endif
      enddo
    ENDIF
    if(debug_level.ge.100) call mpas_log_write('  kmax,maxthe = $i $r', intArgs=(/kmax/), realArgs=(/maxthe/))

  ELSEIF(source.eq.3)THEN
    ! use mixed layer

    IF( (z(2)-z(1)).gt.ml_depth )THEN
      ! the second level is above the mixed-layer depth:  just use the
      ! lowest level

      avgth = th(1)
      avgqv = q(1)
      kmax = 1

    ELSEIF( z(nk).lt.ml_depth )THEN
      ! the top-most level is within the mixed layer:  just use the
      ! upper-most level

      avgth = th(nk)
      avgqv = q(nk)
      kmax = nk

    ELSE
      ! calculate the mixed-layer properties:

      avgth = 0.0
      avgqv = 0.0
      k = 2
      if(debug_level.ge.100) call mpas_log_write('  ml_depth = $r', realArgs=(/ml_depth/))
      if(debug_level.ge.100) call mpas_log_write('  k,z,th,q: $i $r $r $r', intArgs=(/1/), realArgs=(/z(1),th(1),q(1)/))

      do while( (z(k).le.ml_depth) .and. (k.le.nk) )

        if(debug_level.ge.100) call mpas_log_write('$i $r $r $r', intArgs=(/k/), realArgs=(/z(k),th(k),q(k)/))

        avgth = avgth + 0.5*(z(k)-z(k-1))*(th(k)+th(k-1))
        avgqv = avgqv + 0.5*(z(k)-z(k-1))*(q(k)+q(k-1))

        k = k + 1

      enddo

      th2 = th(k-1)+(th(k)-th(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))
      qv2 =  q(k-1)+( q(k)- q(k-1))*(ml_depth-z(k-1))/(z(k)-z(k-1))

      if(debug_level.ge.100) call mpas_log_write('999 $r $r $r', realArgs=(/ml_depth,th2,qv2/))

      avgth = avgth + 0.5*(ml_depth-z(k-1))*(th2+th(k-1))
      avgqv = avgqv + 0.5*(ml_depth-z(k-1))*(qv2+q(k-1))

      if(debug_level.ge.100) call mpas_log_write('$i $r $r $r', intArgs=(/k/), realArgs=(/z(k),th(k),q(k)/))

      avgth = avgth/ml_depth
      avgqv = avgqv/ml_depth

      kmax = 1

    ENDIF

    if(debug_level.ge.100) call mpas_log_write('$r $r', realArgs=real((/avgth,avgqv/),kind=RKIND))

  ELSE

!    call mpas_log_write('')
!    call mpas_log_write('  Unknown value for source')
!    call mpas_log_write('')
!    call mpas_log_write('  source = $i', intArgs=(/source/))
!    call mpas_log_write('', messageType=MPAS_LOG_CRIT)
    call mpas_log_write('getcape: unknown value for source', messageType=MPAS_LOG_ERR)
    return

  ENDIF

!---- define parcel properties at initial location ----!
    narea = 0.0

  if( (source.eq.1).or.(source.eq.2) )then
    k    = kmax
    th2  = th(kmax)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = t(kmax)
    thv2 = thv(kmax)
    qv2  = q(kmax)
    b2   = 0.0
  elseif( source.eq.3 )then
    k    = kmax
    th2  = avgth
    qv2  = avgqv
    thv2 = th2*(1.0+reps*qv2)/(1.0+qv2)
    pi2  = pi(kmax)
    p2   = p(kmax)
    t2   = th2*pi2
    b2   = g*( thv2-thv(kmax) )/thv(kmax)
  endif

    ql2 = 0.0
    qi2 = 0.0
    qt  = qv2

    cape = 0.0
    cin  = 0.0
    lfc  = 0.0

    doit = .true.
    cloud = .false.
    if(adiabat.eq.1.or.adiabat.eq.2)then
      ice = .false.
    else
      ice = .true.
    endif

      the = getthe(p2,t2,t2,qv2)
      if(debug_level.ge.100) call mpas_log_write('  the = $r', realArgs=(/the/))

!---- begin ascent of parcel ----!

      if(debug_level.ge.100)then
        call mpas_log_write('  Start loop:')
        call mpas_log_write('  p2,th2,qv2 = $r $r $r', realArgs=(/p2,th2,qv2/))
      endif

    do while( doit .and. (k.lt.nk) )

        k = k+1
       b1 =  b2

       dp = p(k-1)-p(k)

      if( dp.lt.pinc )then
        nloop = 1
      else
        nloop = 1 + int( dp/pinc )
        dp = dp/float(nloop)
      endif

      do n=1,nloop

         p1 =  p2
         t1 =  t2
        pi1 = pi2
        th1 = th2
        qv1 = qv2
        ql1 = ql2
        qi1 = qi2
        thv1 = thv2

        p2 = p2 - dp
        pi2 = (p2*rp00)**rddcp

        thlast = th1
        i = 0
        not_converged = .true.

        do while( not_converged )
          i = i + 1
          t2 = thlast*pi2
          if(ice)then
            fliq = max(min((t2-233.15)/(273.15-233.15),1.0),0.0)
            fice = 1.0-fliq
          else
            fliq = 1.0
            fice = 0.0
          endif
          qv2 = min( qt , fliq*getqvs(p2,t2) + fice*getqvi(p2,t2) )
          qi2 = max( fice*(qt-qv2) , 0.0 )
          ql2 = max( qt-qv2-qi2 , 0.0 )

          tbar  = 0.5*(t1+t2)
          qvbar = 0.5*(qv1+qv2)
          qlbar = 0.5*(ql1+ql2)
          qibar = 0.5*(qi1+qi2)

          lhv = lv1-lv2*tbar
          lhs = ls1-ls2*tbar
          lhf = lhs-lhv

          rm=rd+rv*qvbar
          cpm=cp+cpv*qvbar+cpl*qlbar+cpi*qibar
          th2=th1*exp(  lhv*(ql2-ql1)/(cpm*tbar)     &
                       +lhs*(qi2-qi1)/(cpm*tbar)     &
                       +(rm/cpm-rd/cp)*log(p2/p1) )

          if(i .gt. 90 .and. debug_level .gt. 0) call mpas_log_write('$i $r $r $r', intArgs=(/i/), realArgs=(/th2,thlast,th2-thlast/))
          if(i .gt. 100)then
!            call mpas_log_write('')
!            call mpas_log_write('  Error:  lack of convergence')
!            call mpas_log_write('')
!            call mpas_log_write('  ... stopping iteration ')
!            call mpas_log_write('')
!            stop 1001
            if (debug_level .gt. 0) then
               call mpas_log_write('getcape: lack of convergence', messageType=MPAS_LOG_ERR)
            end if
            return
          endif
          if( abs(th2-thlast).gt.converge )then
            thlast=thlast+0.3*(th2-thlast)
          else
            not_converged = .false.
          endif
        enddo

        ! Latest pressure increment is complete.  Calculate some
        ! important stuff:

        if( ql2.ge.1.0e-10 ) cloud = .true.

        IF(adiabat.eq.1.or.adiabat.eq.3)THEN
          ! pseudoadiabat
          qt  = qv2
          ql2 = 0.0
          qi2 = 0.0
        ELSEIF(adiabat.le.0.or.adiabat.ge.5)THEN
!          call mpas_log_write('')
!          call mpas_log_write('  Undefined adiabat')
!          call mpas_log_write('')
!          stop 10000
          call mpas_log_write('getcape: Undefined adiabat', messageType=MPAS_LOG_ERR)
          return
        ENDIF

      enddo

      thv2 = th2*(1.0+reps*qv2)/(1.0+qv2+ql2+qi2)
        b2 = g*( thv2-thv(k) )/thv(k)
        dz = -cpdg*0.5*(thv(k)+thv(k-1))*(pi(k)-pi(k-1))

      the = getthe(p2,t2,t2,qv2)

      ! Get contributions to CAPE and CIN:

      if( (b2.ge.0.0) .and. (b1.lt.0.0) )then
        ! first trip into positive area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b2/(b2-b1)
        parea =  0.5*b2*dz*frac
        narea = narea-0.5*b1*dz*(1.0-frac)
        if(debug_level.ge.200)then
          call mpas_log_write('      b1,b2 = $r $r', realArgs=(/b1,b2/))
          call mpas_log_write('      p1,ps,p2 = $r $r $r', realArgs=(/p(k-1),ps,p(k)/))
          call mpas_log_write('      frac = $r', realArgs=(/frac/))
          call mpas_log_write('      parea = $r', realArgs=(/parea/))
          call mpas_log_write('      narea = $r', realArgs=(/narea/))
        endif
        cin  = cin  + narea
        narea = 0.0
      elseif( (b2.lt.0.0) .and. (b1.gt.0.0) )then
        ! first trip into neg area
        ps = p(k-1)+(p(k)-p(k-1))*(0.0-b1)/(b2-b1)
        frac = b1/(b1-b2)
        parea =  0.5*b1*dz*frac
        narea = -0.5*b2*dz*(1.0-frac)
        if(debug_level.ge.200)then
          call mpas_log_write('      b1,b2 = $r $r', realArgs=(/b1,b2/))
          call mpas_log_write('      p1,ps,p2 = $r $r $r', realArgs=(/p(k-1),ps,p(k)/))
          call mpas_log_write('      frac = $r', realArgs=(/frac/))
          call mpas_log_write('      parea = $r', realArgs=(/parea/))
          call mpas_log_write('      narea = $r', realArgs=(/narea/))
        endif
      elseif( b2.lt.0.0 )then
        ! still collecting negative buoyancy
        parea =  0.0
        narea = narea-0.5*dz*(b1+b2)
      else
        ! still collecting positive buoyancy
        parea =  0.5*dz*(b1+b2)
        narea =  0.0
      endif

      cape = cape + max(0.0,parea)

      if(debug_level.ge.200)then
        call mpas_log_write('$r $r $r $r $r $l', realArgs=(/p2,b1,b2,cape,cin/), logicArgs=(/cloud/))
      endif

      if( (p(k).le.10000.0).and.(b2.lt.0.0) )then
        ! stop if b < 0 and p < 100 mb
        doit = .false.
      endif

    enddo

!---- All done ----!

    return
    end subroutine getcape

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real (kind=RKIND) function getqvs(p,t)
    implicit none

    real (kind=RKIND) :: p,t,es

    real (kind=RKIND), parameter :: eps = 287.04/461.5

    es = 611.2*exp(17.67*(t-273.15)/(t-29.65))
    getqvs = eps*es/(p-es)

    return
    end function getqvs

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real (kind=RKIND) function getqvi(p,t)
    implicit none

    real (kind=RKIND) :: p,t,es

    real (kind=RKIND), parameter :: eps = 287.04/461.5

    es = 611.2*exp(21.8745584*(t-273.15)/(t-7.66))
    getqvi = eps*es/(p-es)

    return
    end function getqvi

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

    real (kind=RKIND) function getthe(p,t,td,q)
    implicit none

    real (kind=RKIND) :: p,t,td,q
    real (kind=RKIND) :: tlcl

    if( (td-t).ge.-0.1 )then
      tlcl = t
    else
      tlcl = 56.0 + ( (td-56.0)**(-1) + 0.00125*log(t/td) )**(-1)
    endif

    getthe=t*( (100000.0/p)**(0.2854*(1.0-0.28*q)) )   &
            *exp( ((3376.0/tlcl)-2.54)*q*(1.0+0.81*q) )

    return
    end function getthe

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------

end module convective_diagnostics
