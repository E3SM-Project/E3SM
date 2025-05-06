! Copyright (c) 2013,  Pacific Northwest National Laboratory, (PNNL)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_scaled_sfwf
!> \brief MPAS ocean scale data freshwater fluxes
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
!>  Computes a precipitation factor to multiply the fresh water flux
!>  due to precipitation uniformly to ensure a balance of fresh water
!>  at the ocean surface.
!>  Designed following precip_adjustment subroutine located at:
!>  https://github.com/ESCOMP/POP2-CESM/blob/master/source/forcing_sfwf.F90
!
!-----------------------------------------------------------------------

module ocn_scaled_sfwf

    use mpas_kind_types
    use mpas_derived_types
    use mpas_global_sum_mod
    use mpas_timekeeping
    use mpas_timer

    use ocn_config
    use ocn_mesh

    use shr_const_mod

    use ocn_constants, only: & 
        pi, &      
        rho_fw, & 
        rho_sw, &
        ocn_ref_salinity

    implicit none
    private
    save

    !--------------------------------------------------------------------
    !
    ! Public parameters
    !
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    !
    ! Public member functions
    !
    !--------------------------------------------------------------------

    public :: ocn_scaled_sfwf_init, &
              ocn_update_scaling_factor

    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------

    logical :: scaledSFWFOn
    character (len=*), parameter :: alarmID = 'scaledSFWFUpdateAlarm'

<<<<<<< HEAD
    real (kind=RKIND) :: earthArea ! surface area of earth (m2)

    real(kind=RKIND) :: rhow, areaCellGlobal, volumeCellGlobal

=======
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)
!***********************************************************************

contains


!***********************************************************************
!
<<<<<<< HEAD
!  routine ocn_scaled_sfwf_init
=======
!  routine ocn_init_scaled_sfwf
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)
!
!> \brief   Initialize scaling of data freshwater fluxes
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
<<<<<<< HEAD
!>  Set alarms needed to compute yearly and running means of water fluxes 
!
!-----------------------------------------------------------------------

    subroutine ocn_scaled_sfwf_init(domain)!{{{
=======
!>  Set alarms needed to compute daily and running means of water fluxes 
!
!-----------------------------------------------------------------------

    subroutine ocn_init_scaled_sfwf(domain)!{{{
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)

        type (domain_type), intent(inout) :: domain

        ! Alarm variables
        type (mpas_Time_Type) :: alarmTime
        type (mpas_TimeInterval_type) :: alarmTimeStep

        integer :: err_tmp

<<<<<<< HEAD
        character(len=19) :: timeString

        if (config_use_precip_scaling) then
=======
        if (config_sfwf_balance) then
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)
            scaledSFWFOn = .true.
        else
            scaledSFWFOn = .false.
            return
        endif

<<<<<<< HEAD
        if (trim(config_precip_scaling_mode) /= 'time-dependent' .and.  &
            trim(config_precip_scaling_mode) /= 'constant') then
            call mpas_log_write('scaledSFWFOn = .true. requires config_precip_scaling_mode = time-dependent or constant', &
                                mpas_LOG_CRIT)
        end if

        if (config_precip_scaling_history_years < 1) then
            call mpas_log_write('scaledSFWFOn = .true. requires config_precip_scaling_history_years >= 1', &
                                mpas_LOG_CRIT)
        endif
        
        ! Setup time interval string based on config_precip_scaling_history_years (treated as "years")
        write(timeString, '( I4.4,"-00-00_00:00:00" )') config_precip_scaling_history_years

        ! Setup Alarm for updating of scaling factors for precipitation flux
        alarmTime = mpas_get_clock_time(domain % clock, &
                                        mpas_START_TIME, &
                                        ierr=err_tmp)
        
        call mpas_set_timeInterval(alarmTimeStep, &
                                   timeString=timeString,& 
                                   ierr=err_tmp)

        call mpas_add_clock_alarm(domain % clock, alarmID, alarmTime + alarmTimeStep, &
                                  alarmTimeInterval=alarmTimeStep, ierr=err_tmp)
        ! Density of sea water
        rhow = rho_sw

        ! Compute Earth's surface area (used for normalizing fluxes)
        earthArea = 4.0_RKIND * pi * sphereRadius * sphereRadius


    end subroutine ocn_scaled_sfwf_init!}}}
=======
        if (config_sfwf_balance_type < 0) then
            call mpas_log_write('scaledSFWFOn = .true. requires config_sfwf_balance_type >= 0', &
                                mpas_LOG_CRIT)
        endif

        ! Setup Alarm for updating of scaled SFWF
        alarmTime = mpas_get_clock_time(domain % clock, &
                                        mpas_START_TIME, &
                                        ierr=err_tmp)
        call mpas_set_timeInterval(alarmTimeStep, &
                                   timeString='0000-00-01_00:00:00', &
                                   ierr=err_tmp)
        call mpas_add_clock_alarm(domain % clock, alarmID, alarmTime + alarmTimeStep, &
                                  alarmTimeInterval=alarmTimeStep, ierr=err_tmp)

    end subroutine ocn_init_scaled_sfwf!}}}
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)


!***********************************************************************
!
!  routine ocn_update_scaling_factor
!
!> \brief  Update scaling factors for data freshwater fluxes
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
<<<<<<< HEAD
!>  Accumulate yearly mean of surface freshwater fluxes. 
!>  If we are at the end of a year, used it to derive a factor to scale 
=======
!>  Accumulate daily mean of surface freshwater fluxes. 
!>  If we are at the end of a day, used it to derive a factor to scale 
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)
!>  down/up precipitation fluxes and maintain freshwater balance.  
!
!-----------------------------------------------------------------------

<<<<<<< HEAD
    subroutine ocn_update_scaling_factor(domain,timeLevel)!{{{

       type (domain_type), intent(inout) :: domain
       
       integer, intent(in), optional :: timeLevel

       integer :: timeLevelLocal
       
       integer :: err_tmp

       if (.not. scaledSFWFOn) then
           return
       end if

       if (trim(config_precip_scaling_mode) == 'constant') then
           return
       end if 

       if ( present(timeLevel) ) then
          timeLevelLocal = timeLevel
       else
          timeLevelLocal = 1
       end if

       ! since the current clock time is for the end of the accumulation
       ! intereval, first, accumulate the yearly mean
       call accumulate_mean_fw_fluxes(domain,timeLevelLocal)

       ! then, compute update the history and running mean if we are at the
       ! end of a year
       if(mpas_is_alarm_ringing(domain % clock, alarmID, ierr=err_tmp)) then
#ifdef mpas_DEBUG
            call mpas_log_write('       Computing Scaling Factor For Water Fluxes')
#endif
            call compute_scaling_factor(domain,timeLevelLocal)
            call mpas_reset_clock_alarm(domain % clock, alarmID, ierr=err_tmp)
       endif

    end subroutine ocn_update_scaling_factor!}}}

!***********************************************************************
!
!  routine compute_mass_from_ssh
!
!> \brief  Derive ssh inferred freshwater mass 
!> \author Shixuan Zhang
!> \date   May 2025
!> \details
!>  calculate ssh equivalent total mass over ocean grid point 
!
!-----------------------------------------------------------------------
    subroutine compute_mass_from_ssh(domain,timeLevel,sshMass)!{{{

        type (domain_type), intent(inout) :: domain

        integer, intent(in) :: timeLevel

        real(kind=RKIND), intent(out) :: &
                sshMass

        type (block_type), pointer :: & 
                block_ptr

        type (mpas_pool_type), pointer ::  &
                statePool, &
                meshPool

        integer, pointer :: &
                nCellsSolve, &
                nVertLevels

        real(kind=RKIND), dimension(:), pointer :: &
                areaCell, & 
                ssh

        real (kind=RKIND), dimension(:), allocatable :: &
                localArrayForsshSum

        integer :: iCell, k

        character(len=160) :: m
        
        real (kind=RKIND) :: A, c

        block_ptr => domain % blocklist

        do while (associated(block_ptr))

           call mpas_pool_get_dimension(block_ptr % dimensions, "nCellsSolve", nCellsSolve)
           call mpas_pool_get_dimension(block_ptr % dimensions, "nVertLevels", nVertLevels)

           call mpas_pool_get_subpool(block_ptr % structs, "mesh", meshPool)
           call mpas_pool_get_array(meshPool, "areaCell", areaCell)

           call mpas_pool_get_subpool(block_ptr % structs, "state", statePool)
           call mpas_pool_get_array(statePool, "ssh", ssh, timeLevel) 

           allocate(localArrayForsshSum(nCellsSolve))
           
           localArrayForsshSum(:) = 0.0_RKIND
           
           !$omp parallel
           !$omp do schedule(runtime)
           do iCell=1,nCellsSolve
               localArrayForsshSum(iCell) = ssh(iCell) * areaCell(iCell) * rhow
           enddo
           !$omp end do
           !$omp end parallel

            block_ptr => block_ptr % next
        end do ! block_ptr

        sshMass = mpas_global_sum(localArrayForsshSum, domain % dminfo % comm)

        !-------------------------------------------------------------
        ! Output to log file
        !-------------------------------------------------------------
        if (config_precip_scaling_write_to_logfile) then
            A = earthArea
            c = 1._RKIND / A 
            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write('WATER CONSERVATION CHECK')
            call mpas_log_write(' ')
            call mpas_log_write('Total Mass Equivalent to SSH')
            call mpas_log_write('MPAS-Ocean name          Mass(kg)          Mass(kg/m2)')
            write(m,"('sshMass           ',es16.8,'          ',es16.8)") sshMass,sshMass * c
            call mpas_log_write(m)
            call mpas_log_write(' ')
            call mpas_log_write('----------------------------------------------------------')
        end if

    end subroutine compute_mass_from_ssh!}}}
=======
    subroutine ocn_update_scaling_factor(domain)!{{{

        type (domain_type), intent(inout) :: domain

        integer :: err_tmp

        if (.not. scaledSFWFOn) then
            return
        end if

        ! since the current clock time is for the end of the accumulation
        ! intereval, first, accumulate the daily mean
        call accumulate_mean_fw_fluxes(domain)

       ! then, compute update the history and running mean if we are at the
       ! end of a day
       if(mpas_is_alarm_ringing(domain % clock, alarmID, ierr=err_tmp)) then
#ifdef mpas_DEBUG
            call mpas_log_write('       Computing Scaling Factor For SFWF')
#endif
            call update_scaling_factor(domain)
            call mpas_reset_clock_alarm(domain % clock, alarmID, ierr=err_tmp)
        endif

    end subroutine ocn_update_scaling_factor!}}}

>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)

!***********************************************************************
!
!  routine accumulate_mean_fw_fluxes
!
!> \brief  Accumulate surface water fluxes 
!> \author Shixuan Zhang
!> \date   May 2025
!> \details
!>  Accumulate current surface freshwater and precipitation fluxes 
<<<<<<< HEAD
!>  into the yearly mean value
!
!-----------------------------------------------------------------------

    subroutine accumulate_mean_fw_fluxes(domain,timeLevel)!{{{

        type (domain_type), intent(inout) :: domain

        integer, intent(in) :: timeLevel

=======
!>  into the daily mean value
!
!-----------------------------------------------------------------------

    subroutine accumulate_mean_fw_fluxes(domain)!{{{

        type (domain_type), intent(inout) :: domain

>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)
        type (block_type), pointer :: block_ptr

        type (mpas_pool_type), pointer ::  & 
                forcingPool, & 
<<<<<<< HEAD
                meshPool

        real(kind=RKIND), dimension(:), pointer :: &
                areaCell, &
                rainFlux, &
                snowFlux

        real (kind=RKIND), pointer :: &
                avgprecFlux,  &
                SSHinitial

        integer, pointer :: & 
                nCellsSolve, & 
                nAccumulated

        real (kind=RKIND) :: & 
                precFlux

        real (kind=RKIND), dimension(:), allocatable :: & 
             localArrayForPrecSum

        integer :: iCell, k
=======
                statePool, &
                tracersPool, &
                meshPool

        type (MPAS_timeInterval_type) :: &
                timeStepESMF

        real(kind=RKIND), dimension(:), pointer :: &
                areaCell, &
                rainFlux, &
                snowFlux, &
                sshOld, &
                sshNew

        real(kind=RKIND), dimension(:,:,:), pointer :: &
                activeTracersNew, &
                activeTracersOld
  
        integer, dimension(:), pointer :: & 
                minLevelCell, & 
                maxLevelCell

        real (kind=RKIND), pointer :: &
                avgprecFlux,  &
                avgfwFlux, &
                avgsalFlux

        integer, pointer :: & 
                nCellsSolve, & 
                nVertLevels, &
                indexSalinity, &
                nAccumulated

        real (kind=RKIND) :: & 
                precFlux, & 
                fwFlux, &
                salFlux

        real(kind=RKIND) :: dt

        real (kind=RKIND), dimension(:), allocatable :: & 
             localArrayForPrecSum, &
             localArrayForSSHSum, &
             localArrayForSalinitySum

        integer :: &
             iCell, &
             k

        call mpas_get_timeInterval(timeStepESMF, dt=dt)
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)

        block_ptr => domain % blocklist

        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        ! independent of blocks
        call mpas_pool_get_array(forcingPool, 'avgprecFlux', avgprecFlux)
<<<<<<< HEAD
=======
        call mpas_pool_get_array(forcingPool, 'avgfwFlux', avgfwFlux)
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)

        call mpas_pool_get_dimension(forcingPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedfwFlux', nAccumulated)

<<<<<<< HEAD
        if ( nAccumulated == 0 ) then 
           call mpas_pool_get_array(forcingPool, "SSHinitial", SSHinitial)
           call compute_mass_from_ssh(domain,timeLevel,SSHinitial)
        end if

        do while (associated(block_ptr))

           call mpas_pool_get_dimension(block_ptr % dimensions, "nCellsSolve", nCellsSolve)

           call mpas_pool_get_subpool(block_ptr % structs, "mesh", meshPool)
           call mpas_pool_get_array(meshPool, "areaCell", areaCell)

           call mpas_pool_get_subpool(block_ptr % structs, "forcing", forcingPool)
           call mpas_pool_get_array(forcingPool, 'rainFlux', rainFlux)
           call mpas_pool_get_array(forcingPool, 'snowFlux', snowFlux)
           call mpas_pool_get_array(forcingPool, 'avgprecFlux', avgprecFlux)
           call mpas_pool_get_array(forcingPool, 'nAccumulatedfwFlux', nAccumulated)

           allocate(localArrayForPrecSum(nCellsSolve))

           localArrayForPrecSum(:) = 0.0_RKIND
=======
        do while (associated(block_ptr))

           call mpas_pool_get_dimension(block_ptr % dimensions, "nCellsSolve", nCellsSolve)
           call mpas_pool_get_subpool(block_ptr % structs, "mesh", meshPool)
           call mpas_pool_get_subpool(block_ptr % structs, "forcing", forcingPool)
           call mpas_pool_get_subpool(block_ptr % structs, "state", statePool)

           call mpas_pool_get_array(meshPool, "areaCell", areaCell)
           call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
           call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

           call mpas_pool_get_array(statePool, "layerThickness", layerThickness, 1)
           call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)
           call mpas_pool_get_array(statePool, 'ssh', sshOld, 1)

           call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
           call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)
           call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersNew, 2) 
           call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersOld, 1) 

           call mpas_pool_get_array(forcingPool, 'rainFlux',               rainFlux)
           call mpas_pool_get_array(forcingPool, 'snowFlux',               snowFlux)

           call mpas_pool_get_array(forcingPool, 'avgprecFlux', avgprecFlux)
           call mpas_pool_get_array(forcingPool, 'avgfwFlux', avgfwFlux)

           call mpas_pool_get_array(forcingPool, 'nAccumulatedfwFlux', nAccumulated)

           allocate(localArrayForPrecSum(nCellsSolve))
           allocate(localArrayForSSHSum(nCellsSolve))
           localArrayForPrecSum(:) = 0.0_RKIND
           localArrayForSSHSum(:) = 0.0_RKIND
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)

           !precipitation flux and ssh 
           !$omp parallel
           !$omp do schedule(runtime)
           do iCell=1,nCellsSolve
<<<<<<< HEAD
               localArrayForPrecSum(iCell) = (rainFlux(iCell) + snowFlux(iCell)) * areaCell(iCell) ! unit: kg/s
=======
               localArrayForPrecSum(iCell) = (rainFlux(iCell)+snowFlux(iCell)) * areaCell(iCell) ! unit: kg/s
               localArrayForSSHSum(iCell) = (sshNew(iCell)-sshOld(iCell))* areaCell(iCell) * rho_fw / dt ! unit: kg/s
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)
           enddo
           !$omp end do
           !$omp end parallel

<<<<<<< HEAD
           block_ptr => block_ptr % next
        end do ! block_ptr

        precFlux = mpas_global_sum(localArrayForPrecSum, domain % dminfo % comm) 

        avgprecFlux = (avgprecFlux * nAccumulated + precFlux) / ( nAccumulated + 1 )

        deallocate (localArrayForPrecSum)
=======
           allocate(localArraySalinitySum(nCellsSolve))
           localArraySalinitySum(:) = 0.0_RKIND
           
           !salinity flux
           !convert to 
           !$omp parallel
           !$omp do schedule(runtime) private(k)
           do iCell = 1,nCells
             do k = minLevelCell(iCell), maxLevelCell(iCell)
                localArraySalinitySum(iCell) = localArraySalinitySum(iCell) - & 
                     areaCell(iCell) * layerThickness(k,iCell) * & ! cell volume
                     (activeTracersNew(indexSalinity,k,iCell) - activeTracersOld(indexSalinity,k,iCell)) * & ! * ! cell salinity
                     rho_fw * 1.e3_RKIND / ocn_ref_salinity / dt ! F = - rho * h * 1/S * dS/dt 
             end do
           end do
           !$omp end do
           !$omp end parallel

            block_ptr => block_ptr % next
        end do ! block_ptr

        precFlux = mpas_global_sum(localArrayForPrecSum, domain % dminfo % comm)
        fwFlux = mpas_global_sum(localArrayForSSHSum, domain % dminfo % comm)
        salFlux = mpas_global_sum(localArraySalinitySum, domain % dminfo % comm)

        deallocate (localArrayForPrecSum)
        deallocate (localArrayForSSHSum)
        deallocate (localArraySalinitySum)

        avgprecFlux = (avgprecFlux * nAccumulated + precFlux) / ( nAccumulated + 1 )
        avgfwFlux = (avgfwFlux * nAccumulated + fwFlux) / ( nAccumulated + 1 )
        avgsalFlux = (avgsalFlux * nAccumulated + salFlux) / ( nAccumulated + 1 )
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)

        nAccumulated = nAccumulated + 1

    end subroutine accumulate_mean_fw_fluxes!}}}

<<<<<<< HEAD
!***********************************************************************
!
!  routine compute_scaling_factor
=======

!***********************************************************************
!
!  routine update_scaling_factor
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)
!
!> \brief  Update scaling factors for data freshwater fluxes 
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
!>  Update the net mass fluxes, and use it to scale the freshwater 
!>  fluxes based on the water flux balance relationship
!
!-----------------------------------------------------------------------

<<<<<<< HEAD
    subroutine compute_scaling_factor(domain,timeLevel)!{{{

        type (domain_type), intent(inout) :: domain

        integer, intent(in) :: timeLevel

        type (block_type), pointer :: block_ptr

        type (mpas_pool_type), pointer :: forcingPool
        
        type (mpas_timeInterval_type) :: &
                timeStepESMF

        real (kind=RKIND), pointer :: & 
                avgprecFlux, & 
                avgfwFlux, &
                SSHinitial, &
                SSHfinal, & 
=======
    subroutine update_scaling_factor(domain)!{{{

        type (domain_type), intent(inout) :: domain

        type (block_type), pointer :: block_ptr
        type (mpas_pool_type), pointer :: forcingPool

        real (kind=RKIND), pointer :: & 
                avgprecFlux, & 
                avgfwFlux, & 
                avgsalFlux, & 
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)
                SFWFScalingFactor

        real (kind=RKIND), dimension(:), pointer :: & 
                totalAccumulatedprecFluxHistory, &
<<<<<<< HEAD
                totalAccumulatedfwFluxHistory
=======
                totalAccumulatedfwFluxHistory, &
                totalAccumulatedsalFluxHistory
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)

        integer, pointer :: nValidHistory, nAccumulated

        real (kind=RKIND), dimension(:), allocatable :: &
                tmpprecHistory, &
<<<<<<< HEAD
                tmpfwHistory

        real (kind=RKIND) :: previousprecFlux, previousfwFlux
        real (kind=RKIND) :: precFluxMean, fwFluxMean

        integer :: nHistory, err_tmp        

        real(kind=RKIND) :: dt, dtAvg, v, A, c

        character(len=160) :: m
  
        call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err_tmp)
        call mpas_get_timeInterval(timeStepESMF, dt=dt)

        ! independent of blocks
=======
                tmpfwHistory, &
                tmpsalHistory 

        real (kind=RKIND) :: timeInterval
        real (kind=RKIND) :: previousprecFlux, previousfwFlux, previoussalFlux
        real (kind=RKIND) :: precFluxMean, fwFluxMean, salFluxMean
        integer :: nHistory

>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)
        block_ptr => domain % blocklist
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedfwFlux', nAccumulated)
        call mpas_pool_get_array(forcingPool, 'nValidTotalfwFluxHistory', nValidHistory)
        call mpas_pool_get_array(forcingPool, 'avgprecFlux', avgprecFlux)
        call mpas_pool_get_array(forcingPool, 'avgfwFlux', avgfwFlux)
<<<<<<< HEAD
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedprecFluxHistory', totalAccumulatedprecFluxHistory)
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedfwFluxHistory', totalAccumulatedfwFluxHistory)
        call mpas_pool_get_array(forcingPool, 'SFWFScalingFactor', SFWFScalingFactor)

        call mpas_pool_get_array(forcingPool, "SSHinitial", SSHinitial)
        call mpas_pool_get_array(forcingPool, "SSHfinal", SSHfinal)
        call compute_mass_from_ssh(domain,timeLevel,SSHfinal)

        nHistory = config_precip_scaling_history_years
        dtAvg = dt * nAccumulated

        avgfwFlux = (SSHfinal - SSHinitial) / dtAvg ! unit: kg / s 

        previousprecFlux = totalAccumulatedprecFluxHistory(nValidHistory) 
        previousfwFlux = totalAccumulatedfwFluxHistory(nValidHistory) 
        
        if (nValidHistory == 0) then
            ! Followed the module mpas_ocn_scaled_dismf.F
            ! We keep zero as the first entry in totalAccumulated fluxes 
            nValidHistory = 1
        end if

        if (nHistory > 1) then 
            if (nValidHistory == nHistory) then
                ! We need to shift the history, since it's full
                allocate(tmpprecHistory(nHistory))
                allocate(tmpfwHistory(nHistory))
                tmpprecHistory(:) = totalAccumulatedprecFluxHistory(:)
                tmpfwHistory(:) = totalAccumulatedfwFluxHistory(:)
                totalAccumulatedprecFluxHistory(1:nHistory - 1) = tmpprecHistory(2:nHistory)
                totalAccumulatedfwFluxHistory(1:nHistory - 1) = tmpfwHistory(2:nHistory)
            else
                ! the history isn't full yet, so we just add to the end
                nValidHistory = nValidHistory + 1
            end if
            ! add the new total, the previous total plus the new yearly average
            totalAccumulatedprecFluxHistory(nValidHistory) = avgprecFlux
            totalAccumulatedfwFluxHistory(nValidHistory) = avgfwFlux 
        else 
            totalAccumulatedprecFluxHistory(nValidHistory) = avgprecFlux 
            totalAccumulatedfwFluxHistory(nValidHistory) = avgfwFlux 
        end if


        if (trim(config_precip_scaling_mode) == 'constant') then 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! constant factors setup by user via namelist!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          SFWFScalingFactor = config_precip_scaling_constant_factor

        else if (trim(config_precip_scaling_mode) == 'time-dependent') then
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! computed as the ratio of ssh-inferred      !
          ! mass flux to precipitation flux.
          ! Average values are used when nHistory > 1  !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          precFluxMean = sum(totalAccumulatedprecFluxHistory(1:nHistory)) / nHistory
          fwFluxMean = sum(totalAccumulatedfwFluxHistory(1:nHistory)) / nHistory

          SFWFScalingFactor = SFWFScalingFactor - fwFluxMean / precFluxMean

        end if 

        ! reset yearly averaging of the precipitation and freshwater fluxes 
        nAccumulated = 0
        avgprecFlux = 0.0_RKIND 

        !-------------------------------------------------------------
        ! Output to log file
        !-------------------------------------------------------------
        if (config_precip_scaling_write_to_logfile) then

            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write('WATER CONSERVATION CHECK')
            A = earthArea
            c = 1.0_RKIND / A
            call mpas_log_write('WATER FLUXES')
            call mpas_log_write('MPAS-Ocean name          kg/s (F)                coupler name     short name     kg/m^2/s (F/A)')
            write(m,"('precFluxMean            ',es16.8,'       x2o_Faxa_(rain/snow)     prect     ',f16.8)") & 
                    precFluxMean,precFluxMean*c; call mpas_log_write(m)
            write(m,"('fwFluxMean              ',es16.8,'       ssh                      ssh       ',f16.8)") &
                    fwFluxMean,fwFluxMean*c; call mpas_log_write(m)

            call mpas_log_write('SCALING FACTOR ON PRECIPITATION FLUXES')

            call mpas_log_write(' ')
            write(m, '(A)') 'SFWFScalingMethod       '//trim(config_precip_scaling_mode); call mpas_log_write(m) 

            call mpas_log_write(' ')
            write(m,"('SFWFScalingFactor       ',es16.8,'       precip_fact              fact      ',f16.8)") &
                    SFWFScalingFactor,1.0_RKIND - fwFluxMean / precFluxMean; call mpas_log_write(m)
            call mpas_log_write(' ')
            call mpas_log_write('----------------------------------------------------------')

        end if 

    end subroutine compute_scaling_factor!}}}
=======
        call mpas_pool_get_array(forcingPool, 'avgsalFlux', avgsalFlux)
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedprecFluxHistory', totalAccumulatedprecFluxHistory)
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedfwFluxHistory', totalAccumulatedfwFluxHistory)
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedsalFluxHistory', totalAccumulatedsalFluxHistory)
        call mpas_pool_get_array(forcingPool, 'SFWFScalingFactor', SFWFScalingFactor)

        nHistory = config_sfwf_history_days 

        previousprecFlux = totalAccumulatedprecFluxHistory(nValidHistory)
        previousfwFlux = totalAccumulatedfwFluxHistory(nValidHistory)
        previoussalFlux = totalAccumulatedsalFluxHistory(nValidHistory)

        if (nValidHistory == 0) then
            ! here followed the module mpas_ocn_scaled_dismf.F
            nValidHistory = 1
        end if

        if (nValidHistory == nHistory) then
            ! we need to shift the history, since it's full
            allocate(tmpprecHistory(nHistory))
            allocate(tmpfwHistory(nHistory))
            allocate(tmpsalHistory(nHistory))
            tmpprecHistory(:) = totalAccumulatedprecFluxHistory(:)
            tmpfwHistory(:) = totalAccumulatedfwFluxHistory(:)
            tmpsalHistory(:) = totalAccumulatedsalFluxHistory(:)
            totalAccumulatedprecFluxHistory(1:nHistory - 1) = tmpprecHistory(2:nHistory)
            totalAccumulatedfwFluxHistory(1:nHistory - 1) = tmpfwHistory(2:nHistory)
            totalAccumulatedsalFluxHistory(1:nHistory - 1) = tmpsalHistory(2:nHistory)
        else
            ! the history isn't full yet, so we just add to the end
            nValidHistory = nValidHistory + 1
        end if

        ! add the new total, the previous total plus the new daily average
        totalAccumulatedprecFluxHistory(nValidHistory) = previousprecFlux + SHR_CONST_CDAY * avgprecFlux
        totalAccumulatedfwFluxHistory(nValidHistory) = previousfwFlux + SHR_CONST_CDAY * avgfwFlux
        totalAccumulatedsalFluxHistory(nValidHistory) = previoussalFlux + SHR_CONST_CDAY * avgsalFlux

        timeInterval = SHR_CONST_CDAY * (nValidHistory - 1)
        ! the running mean is the difference between the newest and oldest
        ! totals divided by the time between them
        precFluxMean = (totalAccumulatedprecFluxHistory(nValidHistory) - totalAccumulatedprecFluxHistory(1)) &
                       / timeInterval
        fwFluxMean = (totalAccumulatedfwFluxHistory(nValidHistory) - totalAccumulatedfwFluxHistory(1)) &
                      / timeInterval
        salFluxMean = (totalAccumulatedsalFluxHistory(nValidHistory) - totalAccumulatedsalFluxHistory(1)) &
                      / timeInterval

        SFWFScalingFactor = SFWFScalingFactor - (fwFluxMean + salFluxMean) / precFluxMean

        ! reset daily averaging of the precipitation and freshwater fluxes 
        nAccumulated = 0
        avgprecFlux = 0.0_RKIND 
        avgfwFlux = 0.0_RKIND
        avgsalFlux = 0.0_RKIND

    end subroutine update_scaling_factor!}}}
>>>>>>> f16cc1705a (Code change to dynamically compute a scaling factor for the fresh water flux)

end module ocn_scaled_sfwf
