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
!>    https://github.com/ESCOMP/POP2-CESM/blob/master/source/forcing_sfwf.F90
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

    use shr_kind_mod, only: & 
        SHR_KIND_R8

    use shr_const_mod

    use ocn_constants, only: & 
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

    public :: ocn_init_scaled_sfwf, &
              ocn_update_scaling_factor

    !--------------------------------------------------------------------
    !
    ! Private module variables
    !
    !--------------------------------------------------------------------

    logical :: scaledSFWFOn
    character (len=*), parameter :: alarmID = 'scaledSFWFUpdateAlarm'

    real (kind=RKIND) :: earthAreaE3SM ! surface area of earth in coupler

    real(kind=RKIND) :: rhow, areaCellGlobal, volumeCellGlobal

!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_init_scaled_sfwf
!
!> \brief   Initialize scaling of data freshwater fluxes
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
!>  Set alarms needed to compute yearly and running means of water fluxes 
!
!-----------------------------------------------------------------------

    subroutine ocn_init_scaled_sfwf(domain)!{{{

        type (domain_type), intent(inout) :: domain

        ! Alarm variables
        type (mpas_Time_Type) :: alarmTime
        type (mpas_TimeInterval_type) :: alarmTimeStep

        integer :: err_tmp

        ! taking PI from SHR_CONST_PI in share/util/shr_const_mod.F90 to match coupler
        real (kind=RKIND), parameter :: piE3SM = SHR_CONST_PI
        ! taking earth radius from SHR_CONST_REARTH in share/util/shr_const_mod.F90 to match coupler
        real (kind=RKIND), parameter :: earthRadiusE3SM = SHR_CONST_REARTH ! radius of earth, m

        character(len=19) :: timeString

        if (config_sfwf_balance .and. config_sfwf_balance_type > 0) then
            scaledSFWFOn = .true.
        else
            scaledSFWFOn = .false.
            return
        endif

        if (config_sfwf_history_years < 1) then
            call mpas_log_write('scaledSFWFOn = .true. requires config_sfwf_history_years >= 1', &
                                mpas_LOG_CRIT)
        endif
        
        ! Setup time interval string based on config_sfwf_history_years (treated as "years")
        write(timeString, '( I4.4,"-00-00_00:00:00" )') config_sfwf_history_years

        ! Setup Alarm for updating of scaling factors for precipitation flux
        alarmTime = mpas_get_clock_time(domain % clock, &
                                        mpas_START_TIME, &
                                        ierr=err_tmp)
        
        call mpas_set_timeInterval(alarmTimeStep, &
                                   timeString=timeString,& 
                                   ierr=err_tmp)

        call mpas_add_clock_alarm(domain % clock, alarmID, alarmTime + alarmTimeStep, &
                                  alarmTimeInterval=alarmTimeStep, ierr=err_tmp)
                          
        ! Compute Earth's surface area (used for normalizing fluxes)
        earthAreaE3SM = 4.0_RKIND * piE3SM * earthRadiusE3SM*earthRadiusE3SM

        ! Density of sea water
        rhow = rho_sw

    end subroutine ocn_init_scaled_sfwf!}}}


!***********************************************************************
!
!  routine ocn_update_scaling_factor
!
!> \brief  Update scaling factors for data freshwater fluxes
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
!>  Accumulate yearly mean of surface freshwater fluxes. 
!>  If we are at the end of a year, used it to derive a factor to scale 
!>  down/up precipitation fluxes and maintain freshwater balance.  
!
!-----------------------------------------------------------------------

    subroutine ocn_update_scaling_factor(domain,timeLevel)!{{{

       type (domain_type), intent(inout) :: domain
       
       integer, intent(in), optional :: timeLevel

       integer :: timeLevelLocal
       
       integer :: err_tmp

       if (.not. scaledSFWFOn) then
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
            call update_scaling_factor(domain,timeLevelLocal)
            call mpas_reset_clock_alarm(domain % clock, alarmID, ierr=err_tmp)
       endif

    end subroutine ocn_update_scaling_factor!}}}

!***********************************************************************
!
!  routine compute_ssh_mass
!
!> \brief  Derive ssh inferred freshwater mass 
!> \author Shixuan Zhang
!> \date   May 2025
!> \details
!>  calculate ssh equivalent total mass over ocean grid point 
!
!-----------------------------------------------------------------------
    subroutine compute_ssh_mass(domain,timeLevel,sshMass)!{{{

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
        if (config_sfwf_balance_write_to_logfile) then
            A = earthAreaE3SM
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

    end subroutine compute_ssh_mass!}}}

!***********************************************************************
!
!  routine accumulate_mean_fw_fluxes
!
!> \brief  Accumulate surface water fluxes 
!> \author Shixuan Zhang
!> \date   May 2025
!> \details
!>  Accumulate current surface freshwater and precipitation fluxes 
!>  into the yearly mean value
!
!-----------------------------------------------------------------------

    subroutine accumulate_mean_fw_fluxes(domain,timeLevel)!{{{

        type (domain_type), intent(inout) :: domain

        integer, intent(in) :: timeLevel

        type (block_type), pointer :: block_ptr

        type (mpas_pool_type), pointer ::  & 
                forcingPool, & 
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

        block_ptr => domain % blocklist

        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        ! independent of blocks
        call mpas_pool_get_array(forcingPool, 'avgprecFlux', avgprecFlux)

        call mpas_pool_get_dimension(forcingPool, 'nCellsSolve', nCellsSolve)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedfwFlux', nAccumulated)

        if ( nAccumulated == 0 ) then 
           call mpas_pool_get_array(forcingPool, "SSHinitial", SSHinitial)
           call compute_ssh_mass(domain,timeLevel,SSHinitial)
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

           !precipitation flux and ssh 
           !$omp parallel
           !$omp do schedule(runtime)
           do iCell=1,nCellsSolve
               localArrayForPrecSum(iCell) = (rainFlux(iCell) + snowFlux(iCell)) * areaCell(iCell) ! unit: kg/s
           enddo
           !$omp end do
           !$omp end parallel

           block_ptr => block_ptr % next
        end do ! block_ptr

        precFlux = mpas_global_sum(localArrayForPrecSum, domain % dminfo % comm) 

        avgprecFlux = (avgprecFlux * nAccumulated + precFlux) / ( nAccumulated + 1 )

        deallocate (localArrayForPrecSum)

        nAccumulated = nAccumulated + 1

    end subroutine accumulate_mean_fw_fluxes!}}}

!***********************************************************************
!
!  routine update_scaling_factor
!
!> \brief  Update scaling factors for data freshwater fluxes 
!> \author Shixuan Zhang 
!> \date   May 2025
!> \details
!>  Update the net mass fluxes, and use it to scale the freshwater 
!>  fluxes based on the water flux balance relationship
!
!-----------------------------------------------------------------------

    subroutine update_scaling_factor(domain,timeLevel)!{{{

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
                SFWFScalingFactor

        real (kind=RKIND), dimension(:), pointer :: & 
                totalAccumulatedprecFluxHistory, &
                totalAccumulatedfwFluxHistory

        integer, pointer :: nValidHistory, nAccumulated

        real (kind=RKIND), dimension(:), allocatable :: &
                tmpprecHistory, &
                tmpfwHistory

        real (kind=RKIND) :: previousprecFlux, previousfwFlux
        real (kind=RKIND) :: precFluxMean, fwFluxMean

        integer :: nHistory, err_tmp        

        real(kind=RKIND) :: dt, dtAvg, v, A, c

        character(len=160) :: m
  
        call mpas_set_timeInterval(timeStepESMF, timeString=config_dt, ierr=err_tmp)
        call mpas_get_timeInterval(timeStepESMF, dt=dt)

        ! independent of blocks
        block_ptr => domain % blocklist
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_array(forcingPool, 'nAccumulatedfwFlux', nAccumulated)
        call mpas_pool_get_array(forcingPool, 'nValidTotalfwFluxHistory', nValidHistory)
        call mpas_pool_get_array(forcingPool, 'avgprecFlux', avgprecFlux)
        call mpas_pool_get_array(forcingPool, 'avgfwFlux', avgfwFlux)
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedprecFluxHistory', totalAccumulatedprecFluxHistory)
        call mpas_pool_get_array(forcingPool, 'totalAccumulatedfwFluxHistory', totalAccumulatedfwFluxHistory)
        call mpas_pool_get_array(forcingPool, 'SFWFScalingFactor', SFWFScalingFactor)

        call mpas_pool_get_array(forcingPool, "SSHinitial", SSHinitial)
        call mpas_pool_get_array(forcingPool, "SSHfinal", SSHfinal)
        call compute_ssh_mass(domain,timeLevel,SSHfinal)

        nHistory = config_sfwf_history_years
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

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Scaling factors are computed as the ratio of ssh-inferred mass flux to precipitation flux.
        ! Average values are used when nHistory > 1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        precFluxMean = sum(totalAccumulatedprecFluxHistory(1:nHistory)) / nHistory
        fwFluxMean = sum(totalAccumulatedfwFluxHistory(1:nHistory)) / nHistory 

        if (config_sfwf_use_prev_factor) then 
           SFWFScalingFactor = SFWFScalingFactor - fwFluxMean / precFluxMean
        else
           SFWFScalingFactor = 1.0_RKIND - fwFluxMean / precFluxMean 
        end if 

        ! reset yearly averaging of the precipitation and freshwater fluxes 
        nAccumulated = 0
        avgprecFlux = 0.0_RKIND 

        !-------------------------------------------------------------
        ! Output to log file
        !-------------------------------------------------------------
        if (config_sfwf_balance_write_to_logfile) then
            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write('WATER CONSERVATION CHECK')
            A = earthAreaE3SM
            c = 1.0_RKIND / A
            call mpas_log_write('WATER FLUXES')
            call mpas_log_write('MPAS-Ocean name          kg/s (F)                coupler name     short name     kg/m^2/s (F/A)')
            write(m,"('precFluxMean            ',es16.8,'       x2o_Faxa_(rain/snow)     prect     ',f16.8)") & 
                    precFluxMean,precFluxMean*c; call mpas_log_write(m)
            write(m,"('fwFluxMean              ',es16.8,'       ssh                      ssh       ',f16.8)") &
                    fwFluxMean,fwFluxMean*c; call mpas_log_write(m)

            call mpas_log_write('SCALING FACTOR ON PRECIPITATION FLUXES')
            call mpas_log_write(' ')
            write(m,"('SFWFScalingFactor       ',es16.8,'       precip_fact              fact      ',f16.8)") &
                    SFWFScalingFactor,1.0_RKIND - fwFluxMean / precFluxMean; call mpas_log_write(m)
            call mpas_log_write(' ')
            call mpas_log_write('----------------------------------------------------------')
        end if 

    end subroutine update_scaling_factor!}}}

end module ocn_scaled_sfwf
