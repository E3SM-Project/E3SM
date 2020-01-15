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

    ! Important note: At present, the code in mpas_atm_setup_bdy_masks for
    ! deriving the nearestRelaxationCell field assumes that nSpecZone == 2
    integer, parameter :: nSpecZone = 2
    integer, parameter :: nRelaxZone = 5
    integer, parameter :: nBdyZone = nSpecZone + nRelaxZone

    public :: mpas_atm_update_bdy_tend, &
              mpas_atm_get_bdy_tend, &
              mpas_atm_get_bdy_state, &
              mpas_atm_setup_bdy_masks, &
              mpas_atm_bdy_checks

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
    subroutine mpas_atm_update_bdy_tend(clock, streamManager, block, firstCall, ierr)

        use mpas_constants, only : rvord
        use mpas_stream_manager, only : mpas_stream_mgr_read
        use mpas_log, only : mpas_log_write
        use mpas_derived_types, only : MPAS_STREAM_MGR_NOERR, MPAS_LOG_ERR
        use mpas_timekeeping, only : mpas_get_time

        implicit none

        type (mpas_clock_type), intent(in) :: clock
        type (MPAS_streamManager_type), intent(inout) :: streamManager
        type (block_type), intent(inout) :: block
        logical, intent(in) :: firstCall
        integer, intent(out) :: ierr

        character(len=StrKIND) :: lbc_intv_start_string
        character(len=StrKIND) :: lbc_intv_end_string

        type (mpas_pool_type), pointer :: mesh
        type (mpas_pool_type), pointer :: state
        type (mpas_pool_type), pointer :: lbc
        real (kind=RKIND) :: dt

        integer, pointer :: nCells
        integer, pointer :: nEdges
        integer, pointer :: index_qv

        real (kind=RKIND), dimension(:,:), pointer :: u
        real (kind=RKIND), dimension(:,:), pointer :: ru
        real (kind=RKIND), dimension(:,:), pointer :: rho_edge
        real (kind=RKIND), dimension(:,:), pointer :: w
        real (kind=RKIND), dimension(:,:), pointer :: theta
        real (kind=RKIND), dimension(:,:), pointer :: rtheta_m
        real (kind=RKIND), dimension(:,:), pointer :: rho_zz
        real (kind=RKIND), dimension(:,:), pointer :: rho
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_u
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_ru
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_edge
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_w
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_theta
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rtheta_m
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho_zz
        real (kind=RKIND), dimension(:,:), pointer :: lbc_tend_rho
        real (kind=RKIND), dimension(:,:,:), pointer :: lbc_tend_scalars

        integer, dimension(:,:), pointer :: cellsOnEdge
        real (kind=RKIND), dimension(:,:), pointer :: zz

        integer :: dd_intv, s_intv, sn_intv, sd_intv
        type (MPAS_Time_Type) :: currTime
        type (MPAS_TimeInterval_Type) :: lbc_interval
        character(len=StrKIND) :: read_time
        integer :: iEdge
        integer :: cell1, cell2


        ierr = 0

        call mpas_pool_get_subpool(block % structs, 'mesh', mesh)
        call mpas_pool_get_subpool(block % structs, 'state', state)
        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        if (firstCall) then
            call MPAS_stream_mgr_read(streamManager, streamID='lbc_in', timeLevel=2, whence=MPAS_STREAM_LATEST_BEFORE, &
                                      actualWhen=read_time, ierr=ierr)
            if (ierr /= MPAS_STREAM_MGR_NOERR) then
               call mpas_log_write('Could not read from ''lbc_in'' stream on or before the current date '// &
                                   'to update lateral boundary tendencies', messageType=MPAS_LOG_ERR)
               ierr = 1
            end if
        else
            call mpas_pool_shift_time_levels(lbc)
            call MPAS_stream_mgr_read(streamManager, streamID='lbc_in', timeLevel=2, whence=MPAS_STREAM_EARLIEST_STRICTLY_AFTER, &
                                      actualWhen=read_time, ierr=ierr)
            if (ierr /= MPAS_STREAM_MGR_NOERR) then
               call mpas_log_write('Could not read from ''lbc_in'' stream after the current date '// &
                                   'to update lateral boundary tendencies', messageType=MPAS_LOG_ERR)
               ierr = 1
            end if
        end if
        if (ierr /= 0) then
           return
        end if

        call mpas_set_time(currTime, dateTimeString=trim(read_time))

        !
        ! Compute any derived fields from those that were read from the lbc_in stream
        !
        call mpas_pool_get_array(lbc, 'lbc_u', u, 2)
        call mpas_pool_get_array(lbc, 'lbc_ru', ru, 2)
        call mpas_pool_get_array(lbc, 'lbc_rho_edge', rho_edge, 2)
        call mpas_pool_get_array(lbc, 'lbc_theta', theta, 2)
        call mpas_pool_get_array(lbc, 'lbc_rtheta_m', rtheta_m, 2)
        call mpas_pool_get_array(lbc, 'lbc_rho_zz', rho_zz, 2)
        call mpas_pool_get_array(lbc, 'lbc_rho', rho, 2)
        call mpas_pool_get_array(lbc, 'lbc_scalars', scalars, 2)

        call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
        call mpas_pool_get_dimension(lbc, 'index_qv', index_qv)
        call mpas_pool_get_array(mesh, 'zz', zz)

        ! Compute lbc_rho_zz
        zz(:,nCells+1) = 1.0_RKIND          ! Avoid potential division by zero in the following line
        rho_zz(:,:) = rho(:,:) / zz(:,:)

        ! Average lbc_rho_zz to edges
        do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 > 0 .and. cell2 > 0) then
                rho_edge(:,iEdge) = 0.5_RKIND * (rho_zz(:,cell1) + rho_zz(:,cell2))
            end if
        end do

        ru(:,:) = u(:,:) * rho_edge(:,:)
        rtheta_m(:,:) = theta(:,:) * rho_zz(:,:) * (1.0_RKIND + rvord * scalars(index_qv,:,:))

        if (.not. firstCall) then
            lbc_interval = currTime - LBC_intv_end
            call mpas_get_timeInterval(interval=lbc_interval, DD=dd_intv, S=s_intv, S_n=sn_intv, S_d=sd_intv, ierr=ierr)
            dt = 86400.0_RKIND * real(dd_intv, kind=RKIND) + real(s_intv, kind=RKIND) &
                 + (real(sn_intv, kind=RKIND) / real(sd_intv, kind=RKIND))

            call mpas_pool_get_array(lbc, 'lbc_u', u, 2)
            call mpas_pool_get_array(lbc, 'lbc_ru', ru, 2)
            call mpas_pool_get_array(lbc, 'lbc_rho_edge', rho_edge, 2)
            call mpas_pool_get_array(lbc, 'lbc_w', w, 2)
            call mpas_pool_get_array(lbc, 'lbc_theta', theta, 2)
            call mpas_pool_get_array(lbc, 'lbc_rtheta_m', rtheta_m, 2)
            call mpas_pool_get_array(lbc, 'lbc_rho_zz', rho_zz, 2)
            call mpas_pool_get_array(lbc, 'lbc_rho', rho, 2)
            call mpas_pool_get_array(lbc, 'lbc_scalars', scalars, 2)

            call mpas_pool_get_array(lbc, 'lbc_u', lbc_tend_u, 1)
            call mpas_pool_get_array(lbc, 'lbc_ru', lbc_tend_ru, 1)
            call mpas_pool_get_array(lbc, 'lbc_rho_edge', lbc_tend_rho_edge, 1)
            call mpas_pool_get_array(lbc, 'lbc_w', lbc_tend_w, 1)
            call mpas_pool_get_array(lbc, 'lbc_theta', lbc_tend_theta, 1)
            call mpas_pool_get_array(lbc, 'lbc_rtheta_m', lbc_tend_rtheta_m, 1)
            call mpas_pool_get_array(lbc, 'lbc_rho_zz', lbc_tend_rho_zz, 1)
            call mpas_pool_get_array(lbc, 'lbc_rho', lbc_tend_rho, 1)
            call mpas_pool_get_array(lbc, 'lbc_scalars', lbc_tend_scalars, 1)


            dt = 1.0_RKIND / dt
            lbc_tend_u(:,:) = (u(:,:) - lbc_tend_u(:,:)) * dt
            lbc_tend_ru(:,:) = (ru(:,:) - lbc_tend_ru(:,:)) * dt
            lbc_tend_rho_edge(:,:) = (rho_edge(:,:) - lbc_tend_rho_edge(:,:)) * dt
            lbc_tend_w(:,:) = (w(:,:) - lbc_tend_w(:,:)) * dt
            lbc_tend_theta(:,:) = (theta(:,:) - lbc_tend_theta(:,:)) * dt
            lbc_tend_rtheta_m(:,:) = (rtheta_m(:,:) - lbc_tend_rtheta_m(:,:)) * dt
            lbc_tend_rho_zz(:,:) = (rho_zz(:,:) - lbc_tend_rho_zz(:,:)) * dt
            lbc_tend_rho(:,:) = (rho(:,:) - lbc_tend_rho(:,:)) * dt
            lbc_tend_scalars(:,:,:) = (scalars(:,:,:) - lbc_tend_scalars(:,:,:)) * dt

            !
            ! Logging the lbc start and end times appears to be backwards, but
            ! until the end of this function, LBC_intv_end == the last interval
            ! time and currTime == the next interval time.
            !
            call mpas_get_time(LBC_intv_end, dateTimeString=lbc_intv_start_string)
            call mpas_get_time(currTime, dateTimeString=lbc_intv_end_string)
            call mpas_log_write('----------------------------------------------------------------------')
            call mpas_log_write('Updated lateral boundary conditions. LBCs are now valid')
            call mpas_log_write('from '//trim(lbc_intv_start_string)//' to '//trim(lbc_intv_end_string))
            call mpas_log_write('----------------------------------------------------------------------')

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
    !>   tend_theta(:,:) = mpas_atm_get_bdy_tend(clock, domain % blocklist, nVertLevels, nCells, 'theta', 0.0_RKIND)
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
    !>   theta(:,:) = mpas_atm_get_bdy_state(clock, domain % blocklist, nVertLevels, nCells, 'theta', 0.0_RKIND)
    !>   scalars(1,:,:) = mpas_atm_get_bdy_state(clock, domain % blocklist, nVertLevels, nCells, 'qv', 0.0_RKIND)
    !
    !-----------------------------------------------------------------------
    function mpas_atm_get_bdy_state(clock, block, vertDim, horizDim, field, delta_t) result(return_state)

        use mpas_pool_routines, only : mpas_pool_get_error_level, mpas_pool_set_error_level
        use mpas_derived_types, only : MPAS_POOL_SILENT

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
        integer :: err_level
        integer :: ierr


        currTime = mpas_get_clock_time(clock, MPAS_NOW, ierr)

        lbc_interval = LBC_intv_end - currTime

        call mpas_get_timeInterval(interval=lbc_interval, DD=dd_intv, S=s_intv, S_n=sn_intv, S_d=sd_intv, ierr=ierr)
        dt = 86400.0_RKIND * real(dd_intv, kind=RKIND) + real(s_intv, kind=RKIND) &
             + (real(sn_intv, kind=RKIND) / real(sd_intv, kind=RKIND))

        dt = dt - delta_t

        call mpas_pool_get_subpool(block % structs, 'lbc', lbc)

        !
        ! The first two calls to mpas_pool_get_array, below, may cause harmless warning
        ! messages, which we can silence by setting the pool error level
        !
        err_level = mpas_pool_get_error_level()
        call mpas_pool_set_error_level(MPAS_POOL_SILENT)

        nullify(tend)
        nullify(state)

        call mpas_pool_get_array(lbc, 'lbc_'//trim(field), tend, 1)
        call mpas_pool_get_array(lbc, 'lbc_'//trim(field), state, 2)

        ! Reset the pool error level to its previous value
        call mpas_pool_set_error_level(err_level)

        !
        ! If we have both a tendency and state for this boundary field from the tend and state
        ! pools, then the requested field was not a scalar constituent; otherwise, we need to
        ! query the field as a scalar constituent
        !
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
        ! NB: The specZoneMask{Cell,Edge,Vertex} fields receive a default value of 0.0 from the Registry,
        !     so no need to initialize them here
        !
        where (bdyMaskCell(:) > nRelaxZone) specZoneMaskCell(:) = 1.0_RKIND
        where (bdyMaskEdge(:) > nRelaxZone) specZoneMaskEdge(:) = 1.0_RKIND
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


    !***********************************************************************
    !
    !  routine mpas_atm_bdy_checks
    !
    !> \brief   Checks compatibility of limited-area settings
    !> \author  Michael Duda
    !> \date    12 May 2019
    !> \details 
    !>  This routine checks that settings related to limited-area simulations
    !>  are compatible. Specifically, the following are checked by this routine:
    !>
    !>  1) If config_apply_lbcs = true, the bdyMaskCell field must have non-zero elements
    !>  2) If config_apply_lbcs = false, the bdyMaskCell field must not have non-zero elements
    !>  3) If config_apply_lbcs = true, the lbc_in stream must have a valid input interval
    !>
    !>  If any of the above are not true, this routine prints an error message and
    !>  returns a non-zero value in ierr; otherwise, a value of 0 is returned.
    !
    !-----------------------------------------------------------------------
    subroutine mpas_atm_bdy_checks(dminfo, blockList, streamManager, ierr)

        use mpas_log, only : mpas_log_write
        use mpas_kind_types, only : StrKIND
        use mpas_derived_types, only : dm_info, block_type, mpas_pool_type, MPAS_LOG_ERR, MPAS_STREAM_PROPERTY_RECORD_INTV, &
                                       MPAS_STREAM_MGR_NOERR, MPAS_STREAM_INPUT
        use mpas_stream_manager, only : mpas_stream_mgr_get_property
        use mpas_pool_routines, only : mpas_pool_get_config, mpas_pool_get_dimension, mpas_pool_get_subpool, mpas_pool_get_array
        use mpas_dmpar, only : mpas_dmpar_max_int

        implicit none

        type (dm_info), pointer :: dminfo
        type (block_type), pointer :: blockList
        type (MPAS_streamManager_type), pointer :: streamManager
        integer, intent(out) :: ierr

        character(len=StrKIND) :: input_interval
        logical, pointer :: config_apply_lbcs => null()
        integer, pointer :: nCellsSolve => null()
        type (mpas_pool_type), pointer :: meshPool => null()
        type (block_type), pointer :: block => null()
        integer, dimension(:), pointer :: bdyMaskCell => null()
        integer :: maxvar2d_local, maxvar2d_global


        call mpas_pool_get_config(blocklist % configs, 'config_apply_lbcs', config_apply_lbcs)

        call mpas_log_write('')
        call mpas_log_write('Checking consistency of limited-area settings...')
        call mpas_log_write(' - config_apply_lbcs = $l', logicArgs=(/config_apply_lbcs/))

        !
        ! Check whether any elements of bdyMaskCell have non-zero values
        !
        maxvar2d_local = -huge(maxvar2d_local)
        block => blockList
        do while (associated(block))
            call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
            call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
            call mpas_pool_get_array(meshPool, 'bdyMaskCell', bdyMaskCell)

            maxvar2d_local = max(maxvar2d_local, maxval(bdyMaskCell(1:nCellsSolve)))

            block => block % next
        end do

        call mpas_dmpar_max_int(dminfo, maxvar2d_local, maxvar2d_global)
        call mpas_log_write(' - Maximum value in bdyMaskCell = $i', intArgs=(/maxvar2d_global/))

        !
        ! If there are boundary cells, config_apply_lbcs must be set to true
        !
        if (.not. config_apply_lbcs .and. maxvar2d_global > 0) then
            call mpas_log_write('Boundary cells found in the bdyMaskCell field, but config_apply_lbcs = false.', &
                                messageType=MPAS_LOG_ERR)
            call mpas_log_write('Please ensure that config_apply_lbcs = true for limited-area simulations.', &
                                messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        end if

        !
        ! If there are no boundary cells, config_apply_lbcs must be set to false
        !
        if (config_apply_lbcs .and. maxvar2d_global == 0) then
            call mpas_log_write('config_apply_lbcs = true, but no boundary cells found in the bdyMaskCell field.', &
                                messageType=MPAS_LOG_ERR)
            call mpas_log_write('Please ensure that config_apply_lbcs = false for global simulations.', &
                                messageType=MPAS_LOG_ERR)
            ierr = 1
            return
        end if

        !
        ! If config_apply_lbcs = true, check that the 'lbc_in' stream has a valid input interval
        !
        if (config_apply_lbcs) then
            call mpas_stream_mgr_get_property(streamManager, 'lbc_in', MPAS_STREAM_PROPERTY_RECORD_INTV, &
                                              input_interval, MPAS_STREAM_INPUT, ierr)
            if (ierr /= MPAS_STREAM_MGR_NOERR) then
                call mpas_log_write('Unable to retrieve input interval for the ''lbc_in'' stream.', messageType=MPAS_LOG_ERR)
                ierr = 1
                return
            end if

            call mpas_log_write(' - Input interval for ''lbc_in'' stream = '''//trim(input_interval)//'''')
            if (trim(input_interval) == 'none') then
                call mpas_log_write('Input interval for the ''lbc_in'' stream must be a valid interval '// &
                                    'when config_apply_lbcs = true.', messageType=MPAS_LOG_ERR)
                ierr = 1
                return
            end if
        end if

        call mpas_log_write(' ----- done checking limited-area settings -----')
        call mpas_log_write('')
        ierr = 0

    end subroutine mpas_atm_bdy_checks

end module mpas_atm_boundaries
