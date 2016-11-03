! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module soundings

    use mpas_kind_types, only : RKIND, StrKIND
    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type

    public :: soundings_setup, &
              soundings_compute, &
              soundings_cleanup

    private

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag

    integer :: nSoundings = 0
    logical, allocatable, dimension(:) :: stationOwned
    real (kind=RKIND), allocatable, dimension(:) :: stationLats
    real (kind=RKIND), allocatable, dimension(:) :: stationLons
    integer, allocatable, dimension(:) :: stationCells
    character (len=StrKIND), allocatable, dimension(:) :: stationNames

    type (MPAS_clock_type), pointer :: simulationClock

    real (kind=RKIND) :: pi_const = 2.0_RKIND * asin(1.0_RKIND)


    contains


    !-----------------------------------------------------------------------
    !  routine soundings_setup
    !
    !> \brief Reads sounding locations and sets sounding alarm
    !> \author Michael Duda
    !> \date   20 April 2016
    !> \details
    !>  This routine checks on the existence of a 'sounding_locations.txt' file,
    !>  and, if present, reads sounding locations from this file and determines
    !>  which grid cell contains each of the locations.
    !
    !-----------------------------------------------------------------------
    subroutine soundings_setup(configs, all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_config
        use mpas_io_units, only : stderrUnit, mpas_new_unit, mpas_release_unit
        use mpas_timekeeping, only : MPAS_timeInterval_type, MPAS_time_type, MPAS_set_timeInterval, &
                                     MPAS_get_clock_time, MPAS_add_clock_alarm, MPAS_NOW

        implicit none

        type (MPAS_pool_type), pointer :: configs
        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        character(len=StrKIND), pointer :: soundingInterval

        integer :: i, ierr
        integer :: sndUnit
        real (kind=RKIND) :: station_lat, station_lon
        character (len=StrKIND) :: station_name
        logical :: exists
        integer :: nearestCell
        integer, pointer :: nCells, nCellsSolve, maxEdges
        integer, dimension(:), pointer :: nEdgesOnCell
        integer, dimension(:,:), pointer :: cellsOnCell
        real (kind=RKIND), dimension(:), pointer :: latCell, lonCell
        type (MPAS_timeInterval_type) :: intv
        type (MPAS_time_type) :: now


        simulationClock => simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)

        call mpas_pool_get_config(configs, 'config_sounding_interval', soundingInterval)

        if (trim(soundingInterval) == 'none') then
            return
        end if

        inquire(file='sounding_locations.txt', exist=exists)

        if (.not. exists) then
            write(stderrUnit,*) 'Sounding location file ''sounding_locations.txt'' not found.'
            return
        end if
      
        call MPAS_set_timeInterval(intv, timeString=trim(soundingInterval))
        now = MPAS_get_clock_time(simulationClock, MPAS_NOW)
        call MPAS_add_clock_alarm(simulationClock, 'soundingAlarm', now, alarmTimeInterval=intv) 

        call mpas_new_unit(sndUnit)

        open(sndUnit, file='sounding_locations.txt', form='formatted', status='old', iostat=ierr)

        if (ierr /= 0) then
            write(stderrUnit,*) 'Error opening sounding location file ''sounding_locations.txt'''
            call mpas_release_unit(sndUnit)
            return
        end if

        do
            read(sndUnit, fmt=*, iostat=ierr) station_lat, station_lon, station_name
            if (ierr == 0) then
                nSoundings = nSoundings + 1
            else
                exit
            end if
        end do

        if (nSoundings == 0) then
            close(sndUnit)
            call mpas_release_unit(sndUnit)
            return
        end if

        write(stderrUnit,*) 'Read ', nSoundings, ' sounding locations from ''sounding_locations.txt'''

        rewind(sndUnit)

        allocate(stationOwned(nSoundings))
        allocate(stationLats(nSoundings))
        allocate(stationLons(nSoundings))
        allocate(stationCells(nSoundings))
        allocate(stationNames(nSoundings))

        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_dimension(mesh, 'maxEdges', maxEdges)
        call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
        call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
        call mpas_pool_get_array(mesh, 'latCell', latCell)
        call mpas_pool_get_array(mesh, 'lonCell', lonCell)

        nearestCell = nCells

        do i=1,nSoundings
            read(sndUnit, fmt=*) stationLats(i), stationLons(i), stationNames(i)
            nearestCell = nearest_cell((stationLats(i) * pi_const / 180.0_RKIND), (stationLons(i) * pi_const / 180.0_RKIND), &
                                       nearestCell, nCells, maxEdges, nEdgesOnCell, cellsOnCell, latCell, lonCell)
            if (nearestCell <= nCellsSolve) then
                stationOwned(i) = .true.
                stationCells(i) = nearestCell
            else
                stationOwned(i) = .false.
            end if
        end do

        close(sndUnit)

        call mpas_release_unit(sndUnit)

    end subroutine soundings_setup


    !-----------------------------------------------------------------------
    !  routine soundings_compute
    !
    !> \brief If soundings alarm is ringing, writes soundings text files
    !> \author Michael Duda
    !> \date   20 April 2016
    !> \details
    !>  If this routine is called when the 'soundingAlarm' alarm is ringing,
    !>  each calling task will write the sounding locations within its blocks
    !>  to text files on disk.
    !
    !-----------------------------------------------------------------------
    subroutine soundings_compute()

        use mpas_io_units, only : stderrUnit
        use mpas_derived_types, only : MPAS_pool_type
        use mpas_pool_routines, only : MPAS_pool_get_dimension, MPAS_pool_get_array
        use mpas_timekeeping, only : MPAS_time_type, MPAS_is_alarm_ringing, MPAS_reset_clock_alarm, MPAS_get_clock_time, &
                                     MPAS_get_time, MPAS_NOW
        use mpas_constants, only : rvord

        implicit none

        integer :: iStn, k
        integer, pointer :: nVertLevels, index_qv
        real (kind=RKIND), dimension(:,:), pointer :: pressure_base, pressure_p, uReconstructZonal, uReconstructMeridional, zgrid, &
                                                      theta_m, exner
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
        real (kind=RKIND) :: tmpc, tdpc, dir, spd, rh, log_rh, qvs, pres
        type (MPAS_time_type) :: now
        character(len=StrKIND) :: nowString
        integer :: yyyy, mm, dd, h, m, s
        character(len=StrKIND) :: fname
        character(len=10) :: stid


        if (MPAS_is_alarm_ringing(simulationClock, 'soundingAlarm')) then

            now = MPAS_get_clock_time(simulationClock, MPAS_NOW)
            call mpas_get_time(now, YYYY=yyyy, MM=mm, DD=dd, H=h, M=m, S=s, dateTimeString=nowString)

            call MPAS_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
            call MPAS_pool_get_dimension(state, 'index_qv', index_qv)
            call MPAS_pool_get_array(mesh, 'zgrid', zgrid)
            call MPAS_pool_get_array(state, 'scalars', scalars, 1)
            call MPAS_pool_get_array(state, 'theta_m', theta_m, 1)
            call MPAS_pool_get_array(diag, 'pressure_base', pressure_base)
            call MPAS_pool_get_array(diag, 'pressure_p', pressure_p)
            call MPAS_pool_get_array(diag, 'exner', exner)
            call MPAS_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
            call MPAS_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)

!            write(stderrUnit,*) '--- Writing soundings at '//trim(nowString)//'---'

            do iStn=1,nSoundings
                if (stationOwned(iStn)) then
!                    write(stderrUnit,*) 'Writing sounding for station '//trim(stationNames(iStn))

                    write(fname,'(a,i4.4,i2.2,i2.2,i2.2,i2.2,a)') trim(stationNames(iStn))//'.', yyyy, mm, dd, h, m, '.snd'
                    open(97,file=trim(fname),form='formatted',status='replace')

                    write(stid,'(a)') trim(stationNames(iStn))

                    write(97,'(a)') '  SNPARM = PRES;HGHT;TMPC;DWPC;DRCT;SPED;'
                    write(97,'(a)') ''
                    write(97,'(a,i2.2,i2.2,i2.2,a,i2.2,i2.2)') ' STID = '//stid//'  STNM = 99999       TIME = ', mod(yyyy,100), mm, dd,'/', h, m
                    write(97,'(a,f6.2,a,f7.2,a)') ' SLAT = ', stationLats(iStn), '      SLON = ', stationLons(iStn), '     SELV =    -999'
                    write(97,'(a)') ''
                    write(97,'(a)') '      PRES      HGHT     TMPC     DWPC     DRCT     SPED'

                    do k=1,nVertLevels
                        tmpc = theta_m(k,stationCells(iStn)) / (1.0_RKIND + rvord * scalars(index_qv,k,stationCells(iStn))) * exner(k,stationCells(iStn))
                        pres = pressure_base(k,stationCells(iStn)) + pressure_p(k,stationCells(iStn))
!                        if (tmpc >= 273.15_RKIND) then
                            qvs = rslf(pres, tmpc)
!                        else
!                            qvs = rsif(pres, tmpc)
!                        end if
                        rh = max(1.0e-8,min(1.0,scalars(index_qv,k,stationCells(iStn))/qvs))
                        log_rh = log(rh)
                        tmpc = tmpc - 273.15_RKIND
                        pres = pres * 0.01
                        tdpc = 243.04*(log_rh+(17.625*tmpc/(243.04+tmpc))) / (17.625-log_rh-((17.625*tmpc)/(243.04+tmpc))) 
                        spd = sqrt(uReconstructZonal(k,stationCells(iStn))**2 + uReconstructMeridional(k,stationCells(iStn))**2)
                        if (spd == 0.0) then
                            dir = 0.0
                        else
                            dir = acos(-uReconstructMeridional(k,stationCells(iStn)) / spd)
                            if (uReconstructZonal(k,stationCells(iStn)) > 0.0) then
                                dir = 2.0 * pi_const - dir
                            end if
                            dir = dir * 180.0_RKIND / pi_const
                        end if
                        write(97,'(f10.2,f10.2,f9.2,f9.2,f9.2,f9.2)') &
                                    pres, &
                                    0.5 * (zgrid(k,stationCells(iStn)) + zgrid(k+1,stationCells(iStn))), &               ! Avg to layer midpoint
                                    tmpc, &
                                    tdpc, &
                                    dir, &
                                    spd
                    end do

                    close(97)
                end if
            end do

            call MPAS_reset_clock_alarm(simulationClock, 'soundingAlarm')
        else
!            write(stderrUnit,*) '--- Not yet time to write soundings ---'
        end if
   
    end subroutine soundings_compute


    !-----------------------------------------------------------------------
    !  routine soundings_cleanup
    !
    !> \brief Deallocates memory used to handle soundings
    !> \author Michael Duda
    !> \date   20 April 2016
    !> \details
    !>  This routine should be called last among public routines in this module
    !>  to release any allocated memory used in the storage of sounding
    !>  locations.
    !
    !-----------------------------------------------------------------------
    subroutine soundings_cleanup()

        implicit none

        if (nSoundings > 0) then
            deallocate(stationOwned)
            deallocate(stationLats)
            deallocate(stationLons)
            deallocate(stationCells)
            deallocate(stationNames)
            nSoundings = 0
        end if
   
    end subroutine soundings_cleanup


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finds the MPAS grid cell nearest to (target_lat, target_lon)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer function nearest_cell(target_lat, target_lon, start_cell, nCells, maxEdges, &
                                  nEdgesOnCell, cellsOnCell, latCell, lonCell)
    
        implicit none
    
        real (kind=RKIND), intent(in) :: target_lat, target_lon
        integer, intent(in) :: start_cell
        integer, intent(in) :: nCells, maxEdges
        integer, dimension(nCells), intent(in) :: nEdgesOnCell
        integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
        real (kind=RKIND), dimension(nCells), intent(in) :: latCell, lonCell
    
        integer :: i
        integer :: iCell
        integer :: current_cell
        real (kind=RKIND) :: current_distance, d
        real (kind=RKIND) :: nearest_distance
    
        nearest_cell = start_cell
        current_cell = -1
    
        do while (nearest_cell /= current_cell)
            current_cell = nearest_cell
            current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                               target_lon, 1.0_RKIND)
            nearest_cell = current_cell
            nearest_distance = current_distance
            do i = 1, nEdgesOnCell(current_cell)
                iCell = cellsOnCell(i,current_cell)
                if (iCell <= nCells) then
                    d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0_RKIND)
                    if (d < nearest_distance) then
                        nearest_cell = iCell
                        nearest_distance = d
                    end if
                end if
            end do
        end do
    
    end function nearest_cell
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) 
    !    on a sphere with given radius.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real (kind=RKIND) function sphere_distance(lat1, lon1, lat2, lon2, radius)
    
        implicit none
    
        real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius
        real (kind=RKIND) :: arg1
    
        arg1 = sqrt( sin(0.5*(lat2-lat1))**2 + cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
        sphere_distance = 2.0 * radius * asin(arg1)
    
    end function sphere_distance


!=============================================================================================
!NOTE: functions rslf and rsif are taken from module_mp_thompson temporarily for computing
!      the diagnostic relative humidity. These two functions will be removed from this module
!      when the Thompson cloud microphysics scheme will be restored to MPAS-Dev.
!      Laura D. Fowler (laura@ucar.edu) / 2013-07-11.

!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE LIQUID SATURATION VAPOR MIXING RATIO AS
! A FUNCTION OF TEMPERATURE AND PRESSURE
!
      REAL FUNCTION RSLF(P,T)

      IMPLICIT NONE
      REAL, INTENT(IN):: P, T
      REAL:: ESL,X
      REAL, PARAMETER:: C0= .611583699E03
      REAL, PARAMETER:: C1= .444606896E02
      REAL, PARAMETER:: C2= .143177157E01
      REAL, PARAMETER:: C3= .264224321E-1
      REAL, PARAMETER:: C4= .299291081E-3
      REAL, PARAMETER:: C5= .203154182E-5
      REAL, PARAMETER:: C6= .702620698E-8
      REAL, PARAMETER:: C7= .379534310E-11
      REAL, PARAMETER:: C8=-.321582393E-13

      X=MAX(-80.,T-273.16)

!      ESL=612.2*EXP(17.67*X/(T-29.65))
      ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
      RSLF=.622*ESL/(P-ESL)

!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!    ESL = EXP(54.842763 - 6763.22 / T - 4.210 * ALOG(T) + 0.000367 * T
!        + TANH(0.0415 * (T - 218.8)) * (53.878 - 1331.22
!        / T - 9.44523 * ALOG(T) + 0.014025 * T))

      END FUNCTION RSLF
!+---+-----------------------------------------------------------------+
! THIS FUNCTION CALCULATES THE ICE SATURATION VAPOR MIXING RATIO AS A
! FUNCTION OF TEMPERATURE AND PRESSURE
!
      REAL FUNCTION RSIF(P,T)

      IMPLICIT NONE
      REAL, INTENT(IN):: P, T
      REAL:: ESI,X
      REAL, PARAMETER:: C0= .609868993E03
      REAL, PARAMETER:: C1= .499320233E02
      REAL, PARAMETER:: C2= .184672631E01
      REAL, PARAMETER:: C3= .402737184E-1
      REAL, PARAMETER:: C4= .565392987E-3
      REAL, PARAMETER:: C5= .521693933E-5
      REAL, PARAMETER:: C6= .307839583E-7
      REAL, PARAMETER:: C7= .105785160E-9
      REAL, PARAMETER:: C8= .161444444E-12

      X=MAX(-80.,T-273.16)
      ESI=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
      RSIF=.622*ESI/(P-ESI)

!    ALTERNATIVE
!  ; Source: Murphy and Koop, Review of the vapour pressure of ice and
!             supercooled water for atmospheric applications, Q. J. R.
!             Meteorol. Soc (2005), 131, pp. 1539-1565.
!     ESI = EXP(9.550426 - 5723.265/T + 3.53068*ALOG(T) - 0.00728332*T)

      END FUNCTION RSIF

end module soundings
