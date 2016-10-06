! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_atm_boundaries

    use mpas_derived_types, only : mpas_pool_type, mpas_clock_type, block_type, mpas_time_type, mpas_timeInterval_type, MPAS_NOW, &
                                   MPAS_STREAM_LATEST_BEFORE, MPAS_STREAM_EARLIEST_STRICTLY_AFTER, &
                                   MPAS_streamManager_type
    use mpas_pool_routines, only : mpas_pool_get_array, mpas_pool_get_dimension, mpas_pool_get_subpool, mpas_pool_shift_time_levels
    use mpas_kind_types, only : RKIND, StrKIND
    use mpas_timekeeping, only : mpas_get_clock_time, mpas_get_timeInterval, mpas_set_time, operator(-)
    use mpas_stream_manager, only : mpas_stream_mgr_read

    ! Important note: At present, the code in mpas_atm_setup_bdy_masks for
    ! deriving the nearestRelaxationCell field assumes that nSpecZone == 2
    integer, parameter :: nSpecZone = 2
    integer, parameter :: nRelaxZone = 5
    integer, parameter :: nBdyZone = nSpecZone + nRelaxZone

    public :: mpas_atm_update_bdy_tend, &
              mpas_atm_get_bdy_tend, &
              mpas_atm_get_bdy_state, &
              mpas_atm_setup_bdy_masks

    public :: nBdyZone, nSpecZone, nRelaxZone

    private

    type (MPAS_Time_Type) :: LBC_intv_end


    contains


    !***********************************************************************
    !
    !  routine mpas_atm_update_bdy_tend
    !
    !> \brief   Reads new boundary data and updates the LBC tendencies
    !> \author  Michael Duda
    !> \date    27 September 2016
    !> \details 
    !>  This routine reads from the 'lbc_in' stream all variables in the 'lbc'
    !>  pool. When called with firstCall=.true., the latest time before the
    !>  present is read into time level 2 of the lbc pool; otherwise, the
    !>  contents of time level 2 are shifted to time level 1, the earliest
    !>  time strictly later than the present is read into time level 2, and
    !>  the tendencies for all fields in the lbc pool are computed and stored
    !>  in time level 1.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_update_bdy_tend(clock, streamManager, block, firstCall)

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (MPAS_streamManager_type), intent(inout) :: streamManager
        type (block_type), intent(inout) :: block
        logical, intent(in) :: firstCall

        type (mpas_pool_type), pointer :: mesh
        type (mpas_pool_type), pointer :: state
        type (mpas_pool_type), pointer :: lbc
        real (kind=RKIND) :: dt

        real (kind=RKIND), dimension(:,:), pointer :: u
        real (kind=RKIND), dimension(:,:), pointer :: ru
        real (kind=RKIND), dimension(:,:), pointer :: rho_edge
        real (kind=RKIND), dimension(:,:), pointer :: w
        real (kind=RKIND), dimension(:,:), pointer :: theta_m
        real (kind=RKIND), dimension(:,:), pointer :: rtheta_m
        real (kind=RKIND), dimension(:,:), pointer :: rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_u
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_ru
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_edge
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_w
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_theta_m
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rtheta_m
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_zz
        real (kind=RKIND), dimension(:,:,:), pointer :: lbc_tend_scalars

        integer :: ierr
        integer :: dd_intv, s_intv, sn_intv, sd_intv
        type (MPAS_Time_Type) :: currTime
        type (MPAS_TimeInterval_Type) :: lbc_interval
        character(len=StrKIND) :: read_time


        call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
        call mpas_pool_get_subpool(block % structs, 'state', state)
        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        if (firstCall) then
            call MPAS_stream_mgr_read(streamManager, streamID='lbc_in', timeLevel=2, whence=MPAS_STREAM_LATEST_BEFORE, &
                                      actualWhen=read_time, ierr=ierr)
        else
            call mpas_pool_shift_time_levels(lbc)
            call MPAS_stream_mgr_read(streamManager, streamID='lbc_in', timeLevel=2, whence=MPAS_STREAM_EARLIEST_STRICTLY_AFTER, &
                                      actualWhen=read_time, ierr=ierr)
        end if
        call mpas_set_time(currTime, dateTimeString=trim(read_time))

        !
        ! Compute any derived fields from those that were read from the lbc_in stream
        !
        call mpas_pool_get_array(lbc, 'lbc_u', u, 2)
        call mpas_pool_get_array(lbc, 'lbc_ru', ru, 2)
        call mpas_pool_get_array(lbc, 'lbc_rho_edge', rho_edge, 2)
        call mpas_pool_get_array(lbc, 'lbc_theta_m', theta_m, 2)
        call mpas_pool_get_array(lbc, 'lbc_rtheta_m', rtheta_m, 2)
        call mpas_pool_get_array(lbc, 'lbc_rho_zz', rho_zz, 2)
        ru(:,:) = u(:,:) * rho_edge(:,:)
        rtheta_m(:,:) = theta_m(:,:) * rho_zz(:,:)

        if (.not. firstCall) then
            lbc_interval = currTime - LBC_intv_end
            call mpas_get_timeInterval(interval=lbc_interval, DD=dd_intv, S=s_intv, S_n=sn_intv, S_d=sd_intv, ierr=ierr)
            dt = 86400.0_RKIND * real(dd_intv, kind=RKIND) + real(s_intv, kind=RKIND) &
                 + (real(sn_intv, kind=RKIND) / real(sd_intv, kind=RKIND))

            call mpas_pool_get_array(lbc, 'lbc_u', u, 2)
            call mpas_pool_get_array(lbc, 'lbc_ru', ru, 2)
            call mpas_pool_get_array(lbc, 'lbc_rho_edge', rho_edge, 2)
            call mpas_pool_get_array(lbc, 'lbc_w', w, 2)
            call mpas_pool_get_array(lbc, 'lbc_theta_m', theta_m, 2)
            call mpas_pool_get_array(lbc, 'lbc_rtheta_m', rtheta_m, 2)
            call mpas_pool_get_array(lbc, 'lbc_rho_zz', rho_zz, 2)
            call mpas_pool_get_array(lbc, 'lbc_scalars', scalars, 2)

            call mpas_pool_get_array(lbc, 'lbc_u', lbc_tend_u, 1)
            call mpas_pool_get_array(lbc, 'lbc_ru', lbc_tend_ru, 1)
            call mpas_pool_get_array(lbc, 'lbc_rho_edge', lbc_tend_rho_edge, 1)
            call mpas_pool_get_array(lbc, 'lbc_w', lbc_tend_w, 1)
            call mpas_pool_get_array(lbc, 'lbc_theta_m', lbc_tend_theta_m, 1)
            call mpas_pool_get_array(lbc, 'lbc_rtheta_m', lbc_tend_rtheta_m, 1)
            call mpas_pool_get_array(lbc, 'lbc_rho_zz', lbc_tend_rho_zz, 1)
            call mpas_pool_get_array(lbc, 'lbc_scalars', lbc_tend_scalars, 1)


            dt = 1.0_RKIND / dt
            lbc_tend_u(:,:) = (u(:,:) - lbc_tend_u(:,:)) * dt
            lbc_tend_ru(:,:) = (ru(:,:) - lbc_tend_ru(:,:)) * dt
            lbc_tend_rho_edge(:,:) = (rho_edge(:,:) - lbc_tend_rho_edge(:,:)) * dt
            lbc_tend_w(:,:) = (w(:,:) - lbc_tend_w(:,:)) * dt
            lbc_tend_theta_m(:,:) = (theta_m(:,:) - lbc_tend_theta_m(:,:)) * dt
            lbc_tend_rtheta_m(:,:) = (rtheta_m(:,:) - lbc_tend_rtheta_m(:,:)) * dt
            lbc_tend_rho_zz(:,:) = (rho_zz(:,:) - lbc_tend_rho_zz(:,:)) * dt
            lbc_tend_scalars(:,:,:) = (scalars(:,:,:) - lbc_tend_scalars(:,:,:)) * dt
        end if

        LBC_intv_end = currTime

    end subroutine mpas_atm_update_bdy_tend


    !***********************************************************************
    !
    !  routine mpas_atm_get_bdy_tend
    !
    !> \brief   Returns LBC tendencies a specified delta-t in the future
    !> \author  Michael Duda
    !> \date    28 September 2016
    !> \details 
    !>  This function returns an array providing the tendency for the requested
    !>  progostic variable delta_t in the future from the current time known
    !>  by the simulation clock (which is typically the time at the start of
    !>  the current timestep).
    !>
    !>  The vertDim and horizDim should match the nominal block dimensions of
    !>  the field to be returned by the call; for example, a call to retrieve
    !>  the tendency for the 'u' field would set vertDim=nVertLevels and 
    !>  horizDim=nEdges. This function internally adds 1 to the horizontal
    !>  dimension to account for the "garbage" element.
    !>
    !>  The field is identified by the 'field' argument, and this argument is
    !>  prefixed with 'lbc_' before attempting to retrieve the field from
    !>  the 'lbc' pool. For scalars, the field argument should give the name 
    !>  of the constituent, e.g., 'qv'.
    !>
    !>  Example calls to this function:
    !>  
    !>   tend_u(:,:) = mpas_atm_get_bdy_tend(clock, domain % blocklist, nVertLevels, nEdges, 'u', 0.0_RKIND)
    !>   tend_w(:,:) = mpas_atm_get_bdy_tend(clock, domain % blocklist, nVertLevels+1, nCells, 'w', 0.0_RKIND)
    !>   tend_rho_zz(:,:) = mpas_atm_get_bdy_tend(clock, domain % blocklist, nVertLevels, nCells, 'rho_zz', 0.0_RKIND)
    !>   tend_theta_m(:,:) = mpas_atm_get_bdy_tend(clock, domain % blocklist, nVertLevels, nCells, 'theta_m', 0.0_RKIND)
    !>   tend_scalars(1,:,:) = mpas_atm_get_bdy_tend(clock, domain % blocklist, nVertLevels, nCells, 'qv', 0.0_RKIND)
    !
    !-----------------------------------------------------------------------
    function mpas_atm_get_bdy_tend(clock, block, vertDim, horizDim, field, delta_t) result(return_tend)

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (block_type), intent(inout) :: block
        integer, intent(in) :: vertDim, horizDim
        character(len=*), intent(in) :: field
        real (kind=RKIND), intent(in) :: delta_t

        real (kind=RKIND), dimension(vertDim,horizDim+1) :: return_tend

        type (mpas_pool_type), pointer :: lbc
        integer, pointer :: idx
        real (kind=RKIND), dimension(:,:), pointer :: tend
        real (kind=RKIND), dimension(:,:,:), pointer :: tend_scalars
        integer :: ierr


        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        nullify(tend)
        call mpas_pool_get_array(lbc, 'lbc_'//trim(field), tend, 1)

        if (associated(tend)) then
            return_tend(:,:) = tend(:,:)
        else
            call mpas_pool_get_array(lbc, 'lbc_scalars', tend_scalars, 1)
            call mpas_pool_get_dimension(lbc, 'index_'//trim(field), idx)

            return_tend(:,:) = tend_scalars(idx,:,:)
        end if

    end function mpas_atm_get_bdy_tend


    !***********************************************************************
    !
    !  routine mpas_atm_get_bdy_state
    !
    !> \brief   Returns LBC state at a specified delta-t in the future
    !> \author  Michael Duda
    !> \date    28 September 2016
    !> \details 
    !>  This function returns an array providing the state for the requested
    !>  progostic variable delta_t in the future from the current time known
    !>  by the simulation clock (which is typically the time at the start of
    !>  the current timestep).
    !>
    !>  The vertDim and horizDim should match the nominal block dimensions of
    !>  the field to be returned by the call; for example, a call to retrieve
    !>  the state of the 'u' field would set vertDim=nVertLevels and 
    !>  horizDim=nEdges. This function internally adds 1 to the horizontal
    !>  dimension to account for the "garbage" element.
    !>
    !>  The field is identified by the 'field' argument, and this argument is
    !>  prefixed with 'lbc_' before attempting to retrieve the field from
    !>  the 'lbc' pool. For scalars, the field argument should give the name 
    !>  of the constituent, e.g., 'qv'.
    !>
    !>  Example calls to this function:
    !>  
    !>   u(:,:) = mpas_atm_get_bdy_state(clock, domain % blocklist, nVertLevels, nEdges, 'u', 0.0_RKIND)
    !>   w(:,:) = mpas_atm_get_bdy_state(clock, domain % blocklist, nVertLevels+1, nCells, 'w', 0.0_RKIND)
    !>   rho_zz(:,:) = mpas_atm_get_bdy_state(clock, domain % blocklist, nVertLevels, nCells, 'rho_zz', 0.0_RKIND)
    !>   theta_m(:,:) = mpas_atm_get_bdy_state(clock, domain % blocklist, nVertLevels, nCells, 'theta_m', 0.0_RKIND)
    !>   scalars(1,:,:) = mpas_atm_get_bdy_state(clock, domain % blocklist, nVertLevels, nCells, 'qv', 0.0_RKIND)
    !
    !-----------------------------------------------------------------------
    function mpas_atm_get_bdy_state(clock, block, vertDim, horizDim, field, delta_t) result(return_state)

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (block_type), intent(inout) :: block
        integer, intent(in) :: vertDim, horizDim
        character(len=*), intent(in) :: field
        real (kind=RKIND), intent(in) :: delta_t

        real (kind=RKIND), dimension(vertDim,horizDim+1) :: return_state

        type (mpas_pool_type), pointer :: lbc
        integer, pointer :: idx
        real (kind=RKIND), dimension(:,:), pointer :: tend
        real (kind=RKIND), dimension(:,:), pointer :: state
        real (kind=RKIND), dimension(:,:,:), pointer :: tend_scalars
        real (kind=RKIND), dimension(:,:,:), pointer :: state_scalars
        type (MPAS_Time_Type) :: currTime
        type (MPAS_TimeInterval_Type) :: lbc_interval
        integer :: dd_intv, s_intv, sn_intv, sd_intv
        real (kind=RKIND) :: dt
        integer :: ierr


        currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)

        lbc_interval = LBC_intv_end - currTime

        call mpas_get_timeInterval(interval=lbc_interval, DD=dd_intv, S=s_intv, S_n=sn_intv, S_d=sd_intv, ierr=ierr)
        dt = 86400.0_RKIND * real(dd_intv, kind=RKIND) + real(s_intv, kind=RKIND) &
             + (real(sn_intv, kind=RKIND) / real(sd_intv, kind=RKIND))

        dt = dt - delta_t

        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        nullify(tend)
        nullify(state)
        call mpas_pool_get_array(lbc, 'lbc_'//trim(field), tend, 1)
        call mpas_pool_get_array(lbc, 'lbc_'//trim(field), state, 2)

        if (associated(tend) .and. associated(state)) then
            return_state(:,:) = state(:,:) - dt * tend(:,:)
        else
            call mpas_pool_get_array(lbc, 'lbc_scalars', tend_scalars, 1)
            call mpas_pool_get_array(lbc, 'lbc_scalars', state_scalars, 2)
            call mpas_pool_get_dimension(lbc, 'index_'//trim(field), idx)

            return_state(:,:) = state_scalars(idx,:,:) - dt * tend_scalars(idx,:,:)
        end if

    end function mpas_atm_get_bdy_state


    !***********************************************************************
    !
    !  routine mpas_atm_setup_bdy_masks
    !
    !> \brief   Prepares various fields for boundaries of limited-area
    !> \author  Michael Duda
    !> \date    28 September 2016
    !> \details 
    !>  This routine prepares (1) the mask field needed to distinguish cells in 
    !>  the specified zone from those in the relaxation zone, and (2) a field
    !>  of indices identifying the closest relaxation cell to each cell in 
    !>  the specified zone..
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_setup_bdy_masks(mesh, configs)

        implicit none

        type (mpas_pool_type), intent(inout) :: mesh
        type (mpas_pool_type), intent(in) :: configs

        integer :: iCell, i, j, ii, jj
        real (kind=RKIND) :: d, dmin

        integer, pointer :: nCells
        integer, dimension(:), pointer :: bdyMaskCell, bdyMaskEdge, bdyMaskVertex
        integer, dimension(:), pointer :: nearestRelaxationCell
        integer, dimension(:), pointer :: nEdgesOnCell
        integer, dimension(:,:), pointer :: cellsOnCell
        real (kind=RKIND), dimension(:), pointer :: specZoneMaskCell, specZoneMaskEdge, specZoneMaskVertex
        real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell

        call mpas_pool_get_dimension(mesh, 'nCells', nCells)

        call mpas_pool_get_array(mesh, 'bdyMaskCell', bdyMaskCell)
        call mpas_pool_get_array(mesh, 'bdyMaskEdge', bdyMaskEdge)
        call mpas_pool_get_array(mesh, 'bdyMaskVertex', bdyMaskVertex)
        call mpas_pool_get_array(mesh, 'specZoneMaskCell', specZoneMaskCell)
        call mpas_pool_get_array(mesh, 'specZoneMaskEdge', specZoneMaskEdge)
        call mpas_pool_get_array(mesh, 'specZoneMaskVertex', specZoneMaskVertex)
        call mpas_pool_get_array(mesh, 'nearestRelaxationCell', nearestRelaxationCell)
        call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
        call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
        call mpas_pool_get_array(mesh, 'xCell', xCell)
        call mpas_pool_get_array(mesh, 'yCell', yCell)
        call mpas_pool_get_array(mesh, 'zCell', zCell)

        !
        ! Setup mask identifying cells/edges/vertices in the specified zone
        !
        specZoneMaskCell(:) = 0.0_RKIND
        where (bdyMaskCell(:) > nRelaxZone) specZoneMaskCell(:) = 1.0_RKIND

        specZoneMaskEdge(:) = 0.0_RKIND
        where (bdyMaskEdge(:) > nRelaxZone) specZoneMaskEdge(:) = 1.0_RKIND

        specZoneMaskVertex(:) = 0.0_RKIND
        where (bdyMaskVertex(:) > nRelaxZone) specZoneMaskVertex(:) = 1.0_RKIND


        nearestRelaxationCell(:) = nCells+1

        !
        ! For nearest relaxation cell to inner specified zone, just search
        ! all cellsOnCell with bdyMaskCell == nRelaxZone
        !
        do iCell=1,nCells
            if (bdyMaskCell(iCell) == (nRelaxZone+1)) then
                dmin = 1.0e36
                do j=1,nEdgesOnCell(iCell)
                    i = cellsOnCell(j,iCell)
                    if (bdyMaskCell(i) == nRelaxZone) then
                        d = (xCell(i) - xCell(iCell))**2 + (yCell(i) - yCell(iCell))**2 + (zCell(i) - zCell(iCell))**2
                        if (d < dmin) then
                            dmin = d
                            nearestRelaxationCell(iCell) = i
                        end if
                    end if
                end do
            end if
        end do

        !
        ! For nearest relaxation cell to outer specified zone, search
        ! all cellsOnCell of cellsOnCell
        !
        do iCell=1,nCells
            if (bdyMaskCell(iCell) == (nRelaxZone+2)) then
                dmin = 1.0e36
                do j=1,nEdgesOnCell(iCell)
                    i = cellsOnCell(j,iCell)
                    if (bdyMaskCell(i) == (nRelaxZone+1)) then

                        do jj=1,nEdgesOnCell(i)
                            ii = cellsOnCell(jj,i)
                            if (bdyMaskCell(ii) == nRelaxZone) then

                                d = (xCell(ii) - xCell(iCell))**2 + (yCell(ii) - yCell(iCell))**2 + (zCell(ii) - zCell(iCell))**2
                                if (d < dmin) then
                                    dmin = d
                                    nearestRelaxationCell(iCell) = ii
                                end if

                            end if
                        end do

                    end if
                end do
            end if
        end do

    end subroutine mpas_atm_setup_bdy_masks

end module mpas_atm_boundaries
