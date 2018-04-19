! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_conservation_check
!
!> \brief MPAS sea ice analysis mode member: conservation_check
!> \author Adrian K. Turner
!> \date   9th September 2015
!> \details
!>  MPAS sea ice analysis mode member: conservation_check
!>
!-----------------------------------------------------------------------

module seaice_conservation_check

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager
   use mpas_log, only: mpas_log_write

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

   public :: seaice_bootstrap_conservation_check, &
             seaice_init_conservation_check, &
             seaice_precompute_conservation_check, &
             seaice_compute_conservation_check, &
             seaice_restart_conservation_check, &
             seaice_finalize_conservation_check

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine seaice_bootstrap_conservation_check
!
!> \brief   Bootstrap MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    10th November 2015
!> \details
!>  This routine conducts all boostraps required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_bootstrap_conservation_check(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_bootstrap_conservation_check!}}}

!***********************************************************************
!
!  routine seaice_init_conservation_check
!
!> \brief   Initialize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_init_conservation_check(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_init_conservation_check!}}}

!***********************************************************************
!
!  routine seaice_precompute_conservation_check
!
!> \brief   Precompute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  This routine conducts all pre-computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_precompute_conservation_check(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type(MPAS_pool_type), pointer :: &
           conservationCheckAMPool, &
           conservationCheckEnergyAMPool, &
           conservationCheckMassAMPool, &
           conservationCheckSaltAMPool

      integer, pointer :: &
           performConservationPrecompute

      real(kind=RKIND), pointer :: &
           initialEnergy, &
           initialMass, &
           initialSalt

      err = 0

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
      call MPAS_pool_get_array(conservationCheckAMPool, "performConservationPrecompute", performConservationPrecompute)

      if (performConservationPrecompute == 1) then

         ! zero the accumulated fluxes
         call reset_accumulated_variables(domain)

         ! initial total energy
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "initialEnergy", initialEnergy)

         call compute_total_energy(domain, initialEnergy)

         ! initial total mass
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckMassAM", conservationCheckMassAMPool)
         call MPAS_pool_get_array(conservationCheckMassAMPool, "initialMass", initialMass)

         call compute_total_mass(domain, initialMass)

         ! initial total salt
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "initialSalt", initialSalt)

         call compute_total_salt(domain, initialSalt)

         performConservationPrecompute = 0

      endif

   end subroutine seaice_precompute_conservation_check!}}}

!***********************************************************************
!
!  routine seaice_compute_conservation_check
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_compute_conservation_check(domain, instance, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      logical, pointer :: &
           config_AM_conservationCheck_write_to_logfile

      integer :: &
           ierr

      type(MPAS_pool_type), pointer :: &
           conservationCheckAMPool

      integer, pointer :: &
           performConservationPrecompute

      type(MPAS_Time_type) :: &
           currentTime

      character(len=strKIND) :: &
           timeStr

      err = 0

      call MPAS_pool_get_config(domain % blocklist % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                               config_AM_conservationCheck_write_to_logfile)

      if (config_AM_conservationCheck_write_to_logfile .and. &
          MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then
         call mpas_log_write('==========================================================')
         currentTime = MPAS_get_clock_time(domain % clock, MPAS_NOW, ierr=ierr)
         call MPAS_get_time(currentTime, dateTimeString=timeStr, ierr=ierr)
         call mpas_log_write(' Conservation checks: '//trim(timeStr))
      endif

      ! energy conservation check
      call energy_conservation(domain, err)

      ! mass conservation check
      call mass_conservation(domain, err)

      ! salt conservation check
      call salt_conservation(domain, err)

      if (config_AM_conservationCheck_write_to_logfile .and. &
         MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then
         call mpas_log_write('==========================================================')

         ! set precompute to happen next timestep
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
         call MPAS_pool_get_array(conservationCheckAMPool, "performConservationPrecompute", performConservationPrecompute)
         performConservationPrecompute = 1

      endif

   end subroutine seaice_compute_conservation_check!}}}

!***********************************************************************
!
!  routine energy_conservation
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine energy_conservation(domain, err)

      use ice_constants_colpkg, only: &
           Lfresh, &
           Lvap

      type(domain_type), intent(inout) :: &
           domain

      integer, intent(out) :: &
           err !< Output: error flag

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           conservationCheckEnergyAMPool

      real(kind=RKIND), pointer :: &
           initialEnergy, &
           finalEnergy, &
           energyChange, &
           netEnergyFlux, &
           absoluteEnergyError, &
           relativeEnergyError

      real(kind=RKIND), pointer :: &
           accumulatedSurfaceHeatFlux, &
           accumulatedOceanHeatFlux, &
           accumulatedFreezingPotential, &
           accumulatedSnowfallHeat, &
           accumulatedLatentHeat

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           tracersAggregatePool, &
           icestatePool, &
           shortwavePool, &
           oceanFluxesPool, &
           atmosFluxesPool, &
           atmosCouplingPool, &
           diagnosticsPool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           iceAreaCell, &
           iceAreaCellInitial, &
           absorbedShortwaveFlux, &
           oceanShortwaveFlux, &
           sensibleHeatFlux, &
           longwaveUp, &
           longwaveDown, &
           surfaceHeatFlux, &
           latentHeatFlux, &
           evaporativeWaterFlux, &
           oceanHeatFluxArea, &
           freezingMeltingPotentialInitial, &
           snowfallRate

      real(kind=RKIND), pointer :: &
           dt

      logical, pointer :: &
           config_calc_surface_temperature, &
           config_AM_conservationCheck_write_to_logfile

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell, &
           ierr

      integer, parameter :: &
           nSums = 5

      character(len=17) :: &
           formatString

      err = 0

      call MPAS_pool_get_config(domain % blocklist % configs, "config_dt", dt)
      call MPAS_pool_get_config(domain % blocklist % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                               config_AM_conservationCheck_write_to_logfile)

      !-------------------------------------------------------------
      ! Net heat flux to ice
      !-------------------------------------------------------------

      allocate(sumArray(nSums))
      allocate(sumArrayOut(nSums))

      sumArray = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_config(block % configs, "config_calc_surface_temperature", config_calc_surface_temperature)

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
         call MPAS_pool_get_subpool(block % structs, "shortwave", shortwavePool)
         call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxesPool)
         call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxesPool)
         call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)
         call MPAS_pool_get_subpool(block % structs, "diagnostics", diagnosticsPool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
         call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
         call MPAS_pool_get_array(shortwavePool, "absorbedShortwaveFlux", absorbedShortwaveFlux)
         call MPAS_pool_get_array(oceanFluxesPool, "oceanShortwaveFlux", oceanShortwaveFlux)
         call MPAS_pool_get_array(oceanFluxesPool, "oceanHeatFluxArea", oceanHeatFluxArea)
         call MPAS_pool_get_array(atmosFluxesPool, "sensibleHeatFlux", sensibleHeatFlux)
         call MPAS_pool_get_array(atmosFluxesPool, "longwaveUp", longwaveUp)
         call MPAS_pool_get_array(atmosFluxesPool, "surfaceHeatFlux", surfaceHeatFlux)
         call MPAS_pool_get_array(atmosFluxesPool, "latentHeatFlux", latentHeatFlux)
         call MPAS_pool_get_array(atmosFluxesPool, "evaporativeWaterFlux", evaporativeWaterFlux)
         call MPAS_pool_get_array(atmosCouplingPool, "longwaveDown", longwaveDown)
         call MPAS_pool_get_array(atmosCouplingPool, "snowfallRate", snowfallRate)
         call MPAS_pool_get_array(diagnosticsPool, "freezingMeltingPotentialInitial", freezingMeltingPotentialInitial)

         ! surface heat flux
         if (config_calc_surface_temperature) then

            do iCell = 1, nCellsSolve

               sumArray(1) = sumArray(1) + &
                    (absorbedShortwaveFlux(iCell) - oceanShortwaveFlux(iCell) + &
                     sensibleHeatFlux(iCell) + longwaveUp(iCell)) * iceAreaCell(iCell) * areaCell(iCell) + &
                    longwaveDown(iCell) * iceAreaCellInitial(iCell) * areaCell(iCell)

            enddo ! iCell

         else

            do iCell = 1, nCellsSolve

               sumArray(1) = sumArray(1) + &
                    (surfaceHeatFlux(iCell) - latentHeatFlux(iCell)) * iceAreaCell(iCell) * areaCell(iCell)

            enddo ! iCell

         endif

         do iCell = 1, nCellsSolve

            ! ocean heat flux
            sumArray(2) = sumArray(2) + oceanHeatFluxArea(iCell) * areaCell(iCell)

            ! freezing potential
            sumArray(3) = sumArray(3) + max(0.0_RKIND, freezingMeltingPotentialInitial(iCell)) * areaCell(iCell)

            ! snowfall heat input
            sumArray(4) = sumArray(4) - snowfallRate(iCell) * iceAreaCellInitial(iCell) * areaCell(iCell) * Lfresh

            ! latent heat
            sumArray(5) = sumArray(5) + evaporativeWaterFlux(iCell) * iceAreaCell(iCell) * areaCell(iCell) * Lvap

         enddo ! iCell

         block => block % next
      enddo

      ! perform the sums over processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

      ! accumulate fluxes
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)

      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSurfaceHeatFlux", accumulatedSurfaceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedOceanHeatFlux", accumulatedOceanHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedFreezingPotential", accumulatedFreezingPotential)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSnowfallHeat", accumulatedSnowfallHeat)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLatentHeat", accumulatedLatentHeat)

      accumulatedSurfaceHeatFlux   = accumulatedSurfaceHeatFlux   + sumArrayOut(1)
      accumulatedOceanHeatFlux     = accumulatedOceanHeatFlux     + sumArrayOut(2)
      accumulatedFreezingPotential = accumulatedFreezingPotential + sumArrayOut(3)
      accumulatedSnowfallHeat      = accumulatedSnowfallHeat      + sumArrayOut(4)
      accumulatedLatentHeat        = accumulatedLatentHeat        + sumArrayOut(5)

      ! cleanup
      deallocate(sumArray)
      deallocate(sumArrayOut)

      !-------------------------------------------------------------
      ! Energy conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! get initial energy
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "initialEnergy", initialEnergy)

         ! get final energy
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "finalEnergy", finalEnergy)
         call compute_total_energy(domain, finalEnergy)

         ! compute the energy change
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyChange", energyChange)
         energyChange = finalEnergy - initialEnergy

         ! calculate the final net energy flux to the ice
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "netEnergyFlux", netEnergyFlux)

         netEnergyFlux = &
               accumulatedSurfaceHeatFlux &
             + accumulatedSnowfallHeat &
             + accumulatedLatentHeat &
             - accumulatedOceanHeatFlux &
             - accumulatedFreezingPotential

         ! compute the final energy error
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "absoluteEnergyError", absoluteEnergyError)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "relativeEnergyError", relativeEnergyError)

         absoluteEnergyError = netEnergyFlux * dt - energyChange
         relativeEnergyError = absoluteEnergyError / (finalEnergy - 1.0_RKIND) ! why the minus 1????

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            formatString = "(a32,2x,1pe24.17)"

            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write(' Energy conservation check')
            call mpas_log_write(' ')
            call mpas_log_write(' Initial energy           (J) = $r', realArgs=(/initialEnergy/))
            call mpas_log_write(' Final energy             (J) = $r', realArgs=(/finalEnergy/))
            call mpas_log_write(' Energy change            (J) = $r', realArgs=(/energyChange/))
            call mpas_log_write(' ')
            call mpas_log_write(' Surface heat flux        (W) = $r', realArgs=(/accumulatedSurfaceHeatFlux/))
            call mpas_log_write(' Ocean heat flux          (W) = $r', realArgs=(/accumulatedOceanHeatFlux/))
            call mpas_log_write(' Freezing heat flux       (W) = $r', realArgs=(/accumulatedFreezingPotential/))
            call mpas_log_write(' Snowfall heat flux       (W) = $r', realArgs=(/accumulatedSnowfallHeat/))
            call mpas_log_write(' Latent heat flux         (W) = $r', realArgs=(/accumulatedLatentHeat/))
            call mpas_log_write(' Net energy flux          (W) = $r', realArgs=(/netEnergyFlux/))
            call mpas_log_write(' Net energy flux          (J) = $r', realArgs=(/netEnergyFlux * dt/))
            call mpas_log_write(' ')
            call mpas_log_write(' Absolute energy error    (J) = $r', realArgs=(/absoluteEnergyError/))
            call mpas_log_write(' Absolute energy error    (W) = $r', realArgs=(/absoluteEnergyError / dt/))
            call mpas_log_write(' Relative energy error        = $r', realArgs=(/relativeEnergyError/))

         endif

      endif

    end subroutine energy_conservation

!***********************************************************************
!
!  routine mass_conservation
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    11th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine mass_conservation(domain, err)

      use ice_constants_colpkg, only: &
           rhoi

      type (domain_type), intent(inout) :: &
           domain

      integer, intent(out) :: &
           err !< Output: error flag

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           conservationCheckMassAMPool

      real(kind=RKIND), pointer :: &
           initialMass, &
           finalMass, &
           massChange, &
           netMassFlux, &
           absoluteMassError, &
           relativeMassError

      real(kind=RKIND), pointer :: &
           accumulatedRainfallRate, &
           accumulatedSnowfallRate, &
           accumulatedEvaporation, &
           accumulatedFreshWater, &
           accumulatedFrazilWater

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           tracersAggregatePool, &
           icestatePool, &
           atmosCouplingPool, &
           atmosFluxesPool, &
           oceanFluxesPool, &
           meltGrowthRatesPool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           iceAreaCell, &
           iceAreaCellInitial, &
           rainfallRate, &
           snowfallRate, &
           evaporativeWaterFlux, &
           oceanFreshWaterFluxArea, &
           frazilFormation

      real(kind=RKIND), pointer :: &
           dt

      logical, pointer :: &
           config_update_ocean_fluxes, &
           config_AM_conservationCheck_write_to_logfile

      character(len=strKIND), pointer :: &
           config_thermodynamics_type

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell, &
           ierr

      integer, parameter :: &
           nSums = 5

      character(len=17) :: &
           formatString

      err = 0

      call MPAS_pool_get_config(domain % blocklist % configs, "config_dt", dt)
      call MPAS_pool_get_config(domain % blocklist % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                               config_AM_conservationCheck_write_to_logfile)

      !-------------------------------------------------------------
      ! Net mass flux to ice
      !-------------------------------------------------------------

      allocate(sumArray(nSums))
      allocate(sumArrayOut(nSums))

      sumArray = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_config(block % configs, "config_update_ocean_fluxes", config_update_ocean_fluxes)
         call MPAS_pool_get_config(block % configs, "config_thermodynamics_type", config_thermodynamics_type)

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
         call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)
         call MPAS_pool_get_subpool(block % structs, "atmos_fluxes", atmosFluxesPool)
         call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxesPool)
         call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", meltGrowthRatesPool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
         call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
         call MPAS_pool_get_array(atmosCouplingPool, "rainfallRate", rainfallRate)
         call MPAS_pool_get_array(atmosCouplingPool, "snowfallRate", snowfallRate)
         call MPAS_pool_get_array(atmosFluxesPool, "evaporativeWaterFlux", evaporativeWaterFlux)
         call MPAS_pool_get_array(oceanFluxesPool, "oceanFreshWaterFluxArea", oceanFreshWaterFluxArea)
         call MPAS_pool_get_array(meltGrowthRatesPool, "frazilFormation", frazilFormation)

         do iCell = 1, nCellsSolve

            ! rainfall
            sumArray(1) = sumArray(1) + &
                 rainfallRate(iCell) * iceAreaCellInitial(iCell) * areaCell(iCell)

            ! snowfall
            sumArray(2) = sumArray(2) + &
                 snowfallRate(iCell) * iceAreaCellInitial(iCell) * areaCell(iCell)

            ! evaporation
            sumArray(3) = sumArray(3) + &
                 evaporativeWaterFlux(iCell) * iceAreaCell(iCell) * areaCell(iCell)

            ! fresh water flux to ocean
            sumArray(4) = sumArray(4) + &
                 oceanFreshWaterFluxArea(iCell) * areaCell(iCell)

         enddo ! iCell

         if (config_update_ocean_fluxes .and. trim(config_thermodynamics_type) == "mushy") then

            do iCell = 1, nCellsSolve

               ! frazil ice
               sumArray(5) = sumArray(5) + &
                    (frazilFormation(iCell) * areaCell(iCell) * rhoi) / dt

            enddo ! iCell

         endif

         block => block % next
      enddo

      ! perform the sums over processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

      ! accumulate fluxes
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckMassAM", conservationCheckMassAMPool)

      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedRainfallRate", accumulatedRainfallRate)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedSnowfallRate", accumulatedSnowfallRate)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedEvaporation", accumulatedEvaporation)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedFreshWater", accumulatedFreshWater)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedFrazilWater", accumulatedFrazilWater)

      accumulatedRainfallRate = accumulatedRainfallRate + sumArrayOut(1)
      accumulatedSnowfallRate = accumulatedSnowfallRate + sumArrayOut(2)
      accumulatedEvaporation  = accumulatedEvaporation  + sumArrayOut(3)
      accumulatedFreshWater   = accumulatedFreshWater   + sumArrayOut(4)
      accumulatedFrazilWater  = accumulatedFrazilWater  + sumArrayOut(5)

      ! cleanup
      deallocate(sumArray)
      deallocate(sumArrayOut)

      !-------------------------------------------------------------
      ! Mass conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! get initial mass
         call MPAS_pool_get_array(conservationCheckMassAMPool, "initialMass", initialMass)

         ! get final mass
         call MPAS_pool_get_array(conservationCheckMassAMPool, "finalMass", finalMass)
         call compute_total_mass(domain, finalMass)

         ! compute the energy change
         call MPAS_pool_get_array(conservationCheckMassAMPool, "massChange", massChange)
         massChange = finalMass - initialMass

         ! calculate the final net energy flux to the ice
         call MPAS_pool_get_array(conservationCheckMassAMPool, "netMassFlux", netMassFlux)

         netMassFlux = &
                accumulatedRainfallRate &
              + accumulatedSnowfallRate &
              + accumulatedEvaporation &
              - accumulatedFreshWater &
              + accumulatedFrazilWater

         ! compute the final energy error
         call MPAS_pool_get_array(conservationCheckMassAMPool, "absoluteMassError", absoluteMassError)
         call MPAS_pool_get_array(conservationCheckMassAMPool, "relativeMassError", relativeMassError)

         absoluteMassError = netMassFlux * dt - massChange
         relativeMassError = absoluteMassError / (finalmass + 1.0_RKIND) ! why the plus 1????

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            formatString = "(a32,2x,1pe24.17)"

            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write(' Mass conservation check')
            call mpas_log_write(' ')
            call mpas_log_write(' Initial mass            (kg) = $r', realArgs=(/initialMass/))
            call mpas_log_write(' Final mass              (kg) = $r', realArgs=(/finalMass/))
            call mpas_log_write(' Mass change             (kg) = $r', realArgs=(/massChange/))
            call mpas_log_write(' ')
            call mpas_log_write(' Rainfall mass flux    (kg/s) = $r', realArgs=(/accumulatedRainfallRate/))
            call mpas_log_write(' Snowfall mass flux    (kg/s) = $r', realArgs=(/accumulatedSnowfallRate/))
            call mpas_log_write(' Evaporative mass flux (kg/s) = $r', realArgs=(/accumulatedEvaporation/))
            call mpas_log_write(' Fresh water mass flux (kg/s) = $r', realArgs=(/accumulatedFreshWater/))
            call mpas_log_write(' Frazil water flux     (kg/s) = $r', realArgs=(/accumulatedFrazilWater/))
            call mpas_log_write(' Net mass flux         (kg/s) = $r', realArgs=(/netMassFlux/))
            call mpas_log_write(' Net mass flux           (kg) = $r', realArgs=(/netMassFlux * dt/))
            call mpas_log_write(' ')
            call mpas_log_write(' Absolute mass error     (kg) = $r', realArgs=(/absoluteMassError/))
            call mpas_log_write(' Absolute mass error   (kg/s) = $r', realArgs=(/absoluteMassError / dt/))
            call mpas_log_write(' Relative mass error          = $r', realArgs=(/relativeMassError/))

         endif

      endif

    end subroutine mass_conservation

!***********************************************************************
!
!  routine salt_conservation
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    11th September 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine salt_conservation(domain, err)

      use ice_constants_colpkg, only: &
           rhoi, &
           ice_ref_salinity

      type (domain_type), intent(inout) :: &
           domain

      integer, intent(out) :: &
           err !< Output: error flag

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           conservationCheckSaltAMPool

      real(kind=RKIND), pointer :: &
           initialSalt, &
           finalSalt, &
           saltChange, &
           netSaltFlux, &
           absoluteSaltError, &
           relativeSaltError

      real(kind=RKIND), pointer :: &
           accumulatedOceanSaltFlux, &
           accumulatedFrazilSaltFlux

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           oceanFluxesPool, &
           meltGrowthRatesPool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           oceanSaltFluxArea, &
           frazilFormation

      real(kind=RKIND), pointer :: &
           dt

      logical, pointer :: &
           config_update_ocean_fluxes, &
           config_AM_conservationCheck_write_to_logfile

      character(len=strKIND), pointer :: &
           config_thermodynamics_type

      integer, pointer :: &
           nCellsSolve

      integer :: &
           iCell, &
           ierr

      integer, parameter :: &
           nSums = 2

      character(len=17) :: &
           formatString

      err = 0

      call MPAS_pool_get_config(domain % blocklist % configs, "config_dt", dt)
      call MPAS_pool_get_config(domain % blocklist % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                               config_AM_conservationCheck_write_to_logfile)

      !-------------------------------------------------------------
      ! Net salt flux to ice
      !-------------------------------------------------------------

      allocate(sumArray(nSums))
      allocate(sumArrayOut(nSums))

      sumArray = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_config(block % configs, "config_update_ocean_fluxes", config_update_ocean_fluxes)
         call MPAS_pool_get_config(block % configs, "config_thermodynamics_type", config_thermodynamics_type)

         call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "ocean_fluxes", oceanFluxesPool)
         call MPAS_pool_get_subpool(block % structs, "melt_growth_rates", meltGrowthRatesPool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(oceanFluxesPool, "oceanSaltFluxArea", oceanSaltFluxArea)
         call MPAS_pool_get_array(meltGrowthRatesPool, "frazilFormation", frazilFormation)

         do iCell = 1, nCellsSolve

            ! salt flux to ocean
            sumArray(1) = sumArray(1) + &
                 oceanSaltFluxArea(iCell) * areaCell(iCell)

         enddo ! iCell

         if (config_update_ocean_fluxes .and. trim(config_thermodynamics_type) == "mushy") then

            do iCell = 1, nCellsSolve

               ! frazil ice
               sumArray(2) = sumArray(2) + &
                    (frazilFormation(iCell) * areaCell(iCell) * rhoi * ice_ref_salinity * 0.001_RKIND) / dt

            enddo ! iCell

         endif

         block => block % next
      enddo

      ! perform the sums over processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

      ! accumulate fluxes
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)

      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedOceanSaltFlux", accumulatedOceanSaltFlux)
      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedFrazilSaltFlux", accumulatedFrazilSaltFlux)

      accumulatedOceanSaltFlux  = accumulatedOceanSaltFlux  + sumArrayOut(1)
      accumulatedFrazilSaltFlux = accumulatedFrazilSaltFlux + sumArrayOut(2)

      ! cleanup
      deallocate(sumArray)
      deallocate(sumArrayOut)

      !-------------------------------------------------------------
      ! Salt conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! get initial salt content
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "initialSalt", initialSalt)

         ! get final salt content
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "finalSalt", finalSalt)
         call compute_total_salt(domain, finalSalt)

         ! compute the salt content change
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "saltChange", saltChange)
         saltChange = finalSalt - initialSalt

         ! calculate the final net salt flux to the ice
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "netSaltFlux", netSaltFlux)

         netSaltFlux = &
              - accumulatedOceanSaltFlux &
              + accumulatedFrazilSaltFlux

         ! compute the final salt error
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "absoluteSaltError", absoluteSaltError)
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "relativeSaltError", relativeSaltError)

         absoluteSaltError = netSaltFlux * dt - saltChange
         relativeSaltError = absoluteSaltError / (finalSalt - 1.0_RKIND) ! why the minus 1????

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            formatString = "(a32,2x,1pe24.17)"

            call mpas_log_write('----------------------------------------------------------')
            call mpas_log_write(' Salt conservation check')
            call mpas_log_write(' ')
            call mpas_log_write(' Initial salt            (kg) = $r', realArgs=(/initialSalt/))
            call mpas_log_write(' Final salt              (kg) = $r', realArgs=(/finalSalt/))
            call mpas_log_write(' Salt change             (kg) = $r', realArgs=(/saltChange/))
            call mpas_log_write(' ')
            call mpas_log_write(' Ocean salt flux       (kg/s) = $r', realArgs=(/accumulatedOceanSaltFlux/))
            call mpas_log_write(' Frazil salt flux      (kg/s) = $r', realArgs=(/accumulatedFrazilSaltFlux/))
            call mpas_log_write(' Net salt flux         (kg/s) = $r', realArgs=(/netSaltFlux/))
            call mpas_log_write(' Net salt flux           (kg) = $r', realArgs=(/netSaltFlux * dt/))
            call mpas_log_write(' ')
            call mpas_log_write(' Absolute mass error     (kg) = $r', realArgs=(/absoluteSaltError/))
            call mpas_log_write(' Absolute mass error   (kg/s) = $r', realArgs=(/absoluteSaltError / dt/))
            call mpas_log_write(' Relative salt error          = $r', realArgs=(/relativeSaltError/))

         endif

      endif

    end subroutine salt_conservation

!***********************************************************************
!
!  routine compute_total_energy
!
!> \brief   Compute total energy of sea-ice system
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  Calculate the total energy of the sea-ice system
!
!-----------------------------------------------------------------------

    subroutine compute_total_energy(domain, totalEnergy)

      type (domain_type), intent(inout) :: &
           domain

      real(kind=RKIND), intent(out) :: &
           totalEnergy

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           tracersPool

      integer, pointer :: &
           nCellsSolve, &
           nCategories, &
           nIceLayers, &
           nSnowLayers

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           iceEnthalpy, &
           snowEnthalpy, &
           iceVolumeCategory, &
           snowVolumeCategory

      real(kind=RKIND) :: &
           nIceLayersInverse, &
           nSnowLayersInverse

      integer :: &
           iCell, &
           iCategory, &
           iIceLayer, &
           iSnowLayer

      real(kind=RKIND) :: &
           energy

      energy = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)

         call MPAS_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, 'nCategories', nCategories)
         call MPAS_pool_get_dimension(block % dimensions, 'nIceLayers', nIceLayers)
         call MPAS_pool_get_dimension(block % dimensions, 'nSnowLayers', nSnowLayers)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(tracersPool, "iceEnthalpy", iceEnthalpy, 1)
         call MPAS_pool_get_array(tracersPool, "snowEnthalpy", snowEnthalpy, 1)
         call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)
         call MPAS_pool_get_array(tracersPool, "snowVolumeCategory", snowVolumeCategory, 1)

         nIceLayersInverse  = 1.0_RKIND / real(nIceLayers,RKIND)
         nSnowLayersInverse = 1.0_RKIND / real(nSnowLayers,RKIND)

         do iCell = 1, nCellsSolve
            do iCategory = 1, nCategories

               do iIceLayer = 1, nIceLayers

                  energy = energy + &
                       iceEnthalpy(iIceLayer,iCategory,iCell) * &
                       iceVolumeCategory(1,iCategory,iCell) * &
                       nIceLayersInverse * &
                       areaCell(iCell)

               enddo ! iIceLayer

               do iSnowLayer = 1, nSnowLayers

                  energy = energy + &
                       snowEnthalpy(iSnowLayer,iCategory,iCell) * &
                       snowVolumeCategory(1,iCategory,iCell) * &
                       nSnowLayersInverse * &
                       areaCell(iCell)

               enddo ! iIceLayer

            enddo ! iCategory
         enddo ! iCell

         block => block % next
      enddo

      ! sum across processors
      call MPAS_dmpar_sum_real(domain % dminfo, energy, totalEnergy)

    end subroutine compute_total_energy

!***********************************************************************
!
!  routine compute_total_mass
!
!> \brief   Compute total mass of sea-ice system
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  Calculate the total mass of the sea-ice system
!
!-----------------------------------------------------------------------

    subroutine compute_total_mass(domain, totalMass)

      use ice_constants_colpkg, only: &
           rhoi, &
           rhos

      type (domain_type), intent(inout) :: &
           domain

      real(kind=RKIND), intent(out) :: &
           totalMass

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           tracersPool, &
           tracersAggregatePool

      integer, pointer :: &
           nCellsSolve, &
           nCategories

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           iceVolumeCell, &
           snowVolumeCell

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           iceAreaCategory, &
           iceVolumeCategory, &
           snowVolumeCategory, &
           levelIceArea, &
           levelIceVolume

      logical, pointer :: &
           config_use_topo_meltponds

      integer :: &
           iCell, &
           iCategory

      real(kind=RKIND) :: &
           mass

      mass = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_config(block % configs, "config_use_topo_meltponds", config_use_topo_meltponds)

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)
         call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

         call MPAS_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call MPAS_pool_get_dimension(block % dimensions, 'nCategories', nCategories)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)
         call MPAS_pool_get_array(tracersAggregatePool, "snowVolumeCell", snowVolumeCell)
         call MPAS_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
         call MPAS_pool_get_array(tracersPool, "levelIceArea", levelIceArea, 1)
         call MPAS_pool_get_array(tracersPool, "levelIceVolume", levelIceVolume, 1)

         do iCell = 1, nCellsSolve

            ! ice and snow mass
            mass = mass + &
                 (iceVolumeCell(iCell)  * rhoi + &
                  snowVolumeCell(iCell) * rhos) * areaCell(iCell)

         enddo ! iCell

         if (config_use_topo_meltponds) then

            do iCell = 1, nCellsSolve

               do iCategory = 1, nCategories

                  ! pond mass
                  mass = mass + &
                       iceAreaCategory(1,icategory,iCell) * levelIceArea(1,icategory,iCell) * &
                       levelIceVolume(1,icategory,iCell)  * areaCell(iCell)

               enddo ! iCategory

            enddo ! iCell

         endif

         block => block % next
      enddo

      ! sum across processors
      call MPAS_dmpar_sum_real(domain % dminfo, mass, totalMass)

    end subroutine compute_total_mass

!***********************************************************************
!
!  routine compute_total_salt
!
!> \brief   Compute total salt of sea-ice system
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  Calculate the total salt of the sea-ice system
!
!-----------------------------------------------------------------------

    subroutine compute_total_salt(domain, totalSalt)

      use ice_constants_colpkg, only: &
           ice_ref_salinity, &
           rhoi, &
           rhos

      type (domain_type), intent(inout) :: &
           domain

      real(kind=RKIND), intent(out) :: &
           totalSalt

      type(block_type), pointer :: &
           block

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           tracersAggregatePool

      integer, pointer :: &
           nCellsSolve

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           iceVolumeCell

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           iceVolumeCategory

      integer :: &
           iCell, &
           iCategory

      real(kind=RKIND) :: &
           salt

      salt = 0.0_RKIND

      block => domain % blocklist
      do while (associated(block))

         call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

         call MPAS_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)

         do iCell = 1, nCellsSolve

            ! ice and snow mass
            salt = salt + &
                 iceVolumeCell(iCell) * areaCell(iCell)

         enddo ! iCell

         block => block % next
      enddo

      salt = salt * rhoi * ice_ref_salinity * 0.001_RKIND

      ! sum across processors
      call MPAS_dmpar_sum_real(domain % dminfo, salt, totalSalt)

    end subroutine compute_total_salt

!***********************************************************************
!
!  routine reset_accumulated_variables
!
!> \brief   Reset the accumulated fluxes
!> \author  Adrian K. Turner
!> \date    7th September 2015
!> \details This routine resets accumulated fluxes after the
!> conservation calculation has been performed
!
!-----------------------------------------------------------------------

    subroutine reset_accumulated_variables(domain)

      type(domain_type), intent(inout) :: &
           domain

      type(MPAS_pool_type), pointer :: &
           conservationCheckEnergyAMPool, &
           conservationCheckMassAMPool, &
           conservationCheckSaltAMPool

      real(kind=RKIND), pointer :: &
           accumulatedSurfaceHeatFlux, &
           accumulatedOceanHeatFlux, &
           accumulatedFreezingPotential, &
           accumulatedSnowfallHeat, &
           accumulatedLatentHeat, &
           accumulatedRainfallRate, &
           accumulatedSnowfallRate, &
           accumulatedEvaporation, &
           accumulatedFreshWater, &
           accumulatedFrazilWater, &
           accumulatedOceanSaltFlux, &
           accumulatedFrazilSaltFlux

      ! heat
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)

      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSurfaceHeatFlux", accumulatedSurfaceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedOceanHeatFlux", accumulatedOceanHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedFreezingPotential", accumulatedFreezingPotential)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedSnowfallHeat", accumulatedSnowfallHeat)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "accumulatedLatentHeat", accumulatedLatentHeat)

      accumulatedSurfaceHeatFlux   = 0.0_RKIND
      accumulatedOceanHeatFlux     = 0.0_RKIND
      accumulatedFreezingPotential = 0.0_RKIND
      accumulatedSnowfallHeat      = 0.0_RKIND
      accumulatedLatentHeat        = 0.0_RKIND

      ! mass
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckMassAM", conservationCheckMassAMPool)

      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedRainfallRate", accumulatedRainfallRate)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedSnowfallRate", accumulatedSnowfallRate)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedEvaporation", accumulatedEvaporation)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedFreshWater", accumulatedFreshWater)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "accumulatedFrazilWater", accumulatedFrazilWater)

      accumulatedRainfallRate = 0.0_RKIND
      accumulatedSnowfallRate = 0.0_RKIND
      accumulatedEvaporation  = 0.0_RKIND
      accumulatedFreshWater   = 0.0_RKIND
      accumulatedFrazilWater  = 0.0_RKIND

      ! salt
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)

      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedOceanSaltFlux", accumulatedOceanSaltFlux)
      call MPAS_pool_get_array(conservationCheckSaltAMPool, "accumulatedFrazilSaltFlux", accumulatedFrazilSaltFlux)

      accumulatedOceanSaltFlux  = 0.0_RKIND
      accumulatedFrazilSaltFlux = 0.0_RKIND

    end subroutine reset_accumulated_variables

!***********************************************************************
!
!  routine seaice_restart_conservation_check
!
!> \brief   Save restart for MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_restart_conservation_check(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_restart_conservation_check!}}}

!***********************************************************************
!
!  routine seaice_finalize_conservation_check
!
!> \brief   Finalize MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    9th September 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Seaice analysis member.
!
!-----------------------------------------------------------------------

   subroutine seaice_finalize_conservation_check(domain, instance, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      character(len=*), intent(in) :: instance

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine seaice_finalize_conservation_check!}}}

!-----------------------------------------------------------------------

end module seaice_conservation_check

! vim: foldmethod=marker
