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

     use seaice_constants, only: &
          pii

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

      type(block_type), pointer :: &
           blockPtr

      type(MPAS_pool_type), pointer :: &
           conservationCheckAMPool, &
           conservationCheckAreaAMPool, &
           meshPool

      real(kind=RKIND), dimension(:), pointer :: &
           latCell, &
           areaCell

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      real(kind=RKIND), pointer :: &
           earthRadius, &
           earthArea

      real(kind=RKIND), dimension(:), pointer :: &
           domainArea

      integer, pointer :: &
           nCellsSolve, &
           nHemispheres

      integer :: &
           iCell, &
           iHemisphere

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      integer, parameter :: &
           nVars = 1

      integer :: &
           nSums

      err = 0

      call MPAS_pool_get_config(domain % configs, "config_earth_radius", earthRadius)

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)
      nSums = nHemispheres * nVars

      allocate(sumArray(nSums))
      allocate(sumArrayOut(nSums))

      sumArray = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_dimension(blockPtr % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAreaAM", conservationCheckAreaAMPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)

         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)
         call MPAS_pool_get_array(conservationCheckAreaAMPool, "earthArea", earthArea)
         call MPAS_pool_get_array(conservationCheckAreaAMPool, "domainArea", domainArea)
         call MPAS_pool_get_array(meshPool, "latCell", latCell)
         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

         earthArea = 4.0_RKIND * pii * earthRadius**2

         do iCell = 1, nCellsSolve

            ! global
            cellInHemisphere(1,iCell) = 1

            ! northern hemisphere
            if (latCell(iCell) >= 0.0) then
               cellInHemisphere(2,iCell) = 1
            else
               cellInHemisphere(2,iCell) = 0
            endif

            ! southern hemisphere
            if (latCell(iCell) < 0.0) then
               cellInHemisphere(3,iCell) = 1
            else
               cellInHemisphere(3,iCell) = 0
            endif

            ! domain area
            do iHemisphere = 1, nHemispheres
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  sumArray(iHemisphere) = sumArray(iHemisphere) + areaCell(iCell)
               endif
            enddo ! iHemisphere

         enddo ! iCell

         blockPtr => blockPtr % next
      enddo

      ! perform the sums over processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAreaAM", conservationCheckAreaAMPool)

      call MPAS_pool_get_array(conservationCheckAreaAMPool, "domainArea", domainArea)

      do iHemisphere = 1, nHemispheres
         domainArea(iHemisphere) = sumArrayOut(iHemisphere)
      enddo ! iHemisphere

      ! cleanup
      deallocate(sumArray)
      deallocate(sumArrayOut)

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
           conservationCheckSaltAMPool, &
           conservationCheckCarbonAMPool

      integer, pointer :: &
           performConservationPrecompute

      logical, pointer :: &
           config_use_column_biogeochemistry

      real(kind=RKIND), dimension(:), pointer :: &
           initialEnergy, &
           initialMass, &
           initialSalt, &
           initialCarbon

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

         ! initial total carbon
         call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)

         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckCarbonAM", conservationCheckCarbonAMPool)
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "initialCarbon", initialCarbon)

         if (config_use_column_biogeochemistry) then
            call compute_total_carbon(domain, initialCarbon)
         else
            initialCarbon = 0.0_RKIND
         endif

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
           nAccumulate, &
           performConservationPrecompute

      logical, pointer :: &
           config_use_column_biogeochemistry

      type(MPAS_Time_type) :: &
           currentTime

      character(len=strKIND) :: &
           timeStr

      err = 0

      ! number of accumulations
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
      call MPAS_pool_get_array(conservationCheckAMPool, "nAccumulate", nAccumulate)
      nAccumulate = nAccumulate + 1

      call MPAS_pool_get_config(domain % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                   config_AM_conservationCheck_write_to_logfile)

      if (config_AM_conservationCheck_write_to_logfile .and. &
          MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then
         call mpas_log_write('===================================================================================')
         currentTime = MPAS_get_clock_time(domain % clock, MPAS_NOW, ierr=ierr)
         call MPAS_get_time(currentTime, dateTimeString=timeStr, ierr=ierr)
         call mpas_log_write(' Conservation checks: '//trim(timeStr))
      endif

      ! area analysis
      call area_analysis(domain, err)

      ! energy conservation check
      call energy_conservation(domain, err)

      ! mass conservation check
      call mass_conservation(domain, err)

      ! salt conservation check
      call salt_conservation(domain, err)

      ! initial total carbon
      call MPAS_pool_get_config(domain % configs, "config_use_column_biogeochemistry", config_use_column_biogeochemistry)
      if (config_use_column_biogeochemistry) call carbon_conservation(domain, err)

      if (config_AM_conservationCheck_write_to_logfile .and. &
         MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then
         call mpas_log_write('===================================================================================')
      endif

      ! set precompute to happen next timestep after output
      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
         call MPAS_pool_get_array(conservationCheckAMPool, "performConservationPrecompute", performConservationPrecompute)
         performConservationPrecompute = 1
      endif

   end subroutine seaice_compute_conservation_check!}}}

!***********************************************************************
!
!  routine area_analysis
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Adrian K. Turner
!> \date    1st April 2021
!> \details
!>  Analyses areas used in the model.
!
!-----------------------------------------------------------------------

   subroutine area_analysis(domain, err)

      type(domain_type), intent(inout) :: &
           domain

      integer, intent(out) :: &
           err !< Output: error flag

      type(block_type), pointer :: &
           blockPtr

      type(MPAS_pool_type), pointer :: &
           conservationCheckAMPool, &
           conservationCheckAreaAMPool, &
           meshPool, &
           tracersAggregatePool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           latCell, &
           iceAreaCell

      real(kind=RKIND), pointer :: &
           earthRadius, &
           earthArea

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      real(kind=RKIND), dimension(:), pointer :: &
           domainArea, &
           accumulatedSeaIceArea

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      real(kind=RKIND) :: &
           seaiceAreaCell

      logical, pointer :: &
           config_AM_conservationCheck_write_to_logfile

      integer, pointer :: &
           nCellsSolve,   &
           nHemispheres,  &
           nAccumulate

      integer :: &
           iCell, &
           iHemisphere, &
           ierr

      integer, parameter :: &
           nVars = 1

      integer :: &
           nSums

      call MPAS_pool_get_config(domain % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                   config_AM_conservationCheck_write_to_logfile)

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)
      nSums = nHemispheres * nVars

      allocate(sumArray(nSums))
      allocate(sumArrayOut(nSums))

      sumArray = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_dimension(blockPtr % dimensions, "nCellsSolve", nCellsSolve)
         call MPAS_pool_get_dimension(blockPtr % dimensions, "nHemispheres", nHemispheres)

         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(meshPool, "latCell", latCell)
         call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)

         do iCell = 1, nCellsSolve

            seaiceAreaCell = areaCell(iCell) * iceAreaCell(iCell)

            do iHemisphere = 1, nHemispheres
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  sumArray(iHemisphere) = sumArray(iHemisphere) + seaiceAreaCell
               endif
            enddo ! iHemisphere

         enddo ! iCell

         blockPtr => blockPtr % next
      enddo

      ! perform the sums over processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

      ! accumulate fluxes
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAreaAM", conservationCheckAreaAMPool)

      call MPAS_pool_get_array(conservationCheckAreaAMPool, "accumulatedSeaIceArea", accumulatedSeaIceArea)

      do iHemisphere = 1, nHemispheres
         accumulatedSeaIceArea(iHemisphere) = accumulatedSeaIceArea(iHemisphere) + sumArrayOut(iHemisphere)
      enddo ! iHemisphere

      ! cleanup
      deallocate(sumArray)
      deallocate(sumArrayOut)

      !-------------------------------------------------------------
      ! Area analysis
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         call MPAS_pool_get_array(conservationCheckAMPool, "nAccumulate", nAccumulate)

         accumulatedSeaIceArea(:) = accumulatedSeaIceArea(:)/nAccumulate

         call MPAS_pool_get_config(domain % configs, "config_earth_radius", earthRadius)

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         call MPAS_pool_get_array(conservationCheckAreaAMPool, "earthArea", earthArea)
         call MPAS_pool_get_array(conservationCheckAreaAMPool, "domainArea", domainArea)


         if (config_AM_conservationCheck_write_to_logfile) then

            call mpas_log_write('-----------------------------------------------------------------------------------')
            call mpas_log_write(' Area analysis                                 Global             NH             SH')
            call mpas_log_write(' ')
            call mpas_log_write(' Earth radius                   (m) = '//trim(hemisphere_format((/earthRadius/))))
            call mpas_log_write(' Earth area                    (m2) = '//trim(hemisphere_format((/earthArea/))))
            call mpas_log_write(' Domain area                   (m2) = '//trim(hemisphere_format(domainArea)))
            call mpas_log_write(' Sea-ice area                  (m2) = '//trim(hemisphere_format(accumulatedSeaIceArea)))

         endif

      endif

    end subroutine area_analysis

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
           blockPtr

      type(MPAS_pool_type), pointer :: &
           conservationCheckEnergyAMPool, &
           conservationCheckAMPool, &
           conservationCheckAreaAMPool

      real(kind=RKIND), dimension(:), pointer :: &
           initialEnergy, &
           finalEnergy, &
           energyChange, &
           energyChangeFlux, &
           netEnergyFlux, &
           absoluteEnergyError, &
           relativeEnergyError

      real(kind=RKIND), dimension(:), pointer :: &
           energyConsSurfaceHeatFlux, &
           energyConsAbsorbedShortwaveFlux, &
           energyConsOceanShortwaveFlux, &
           energyConsSensibleHeatFlux, &
           energyConsLongwaveUp, &
           energyConsLongwaveDown, &
           energyConsOceanHeatFlux, &
           energyConsFreezingPotential, &
           energyConsSnowfallHeat, &
           energyConsLatentHeat

      real(kind=RKIND) :: &
           surfaceHeatCell, &
           surfaceHeatFluxCell, &
           absorbedShortwaveFluxCell, &
           oceanShortwaveFluxCell, &
           sensibleHeatFluxCell, &
           longwaveUpCell, &
           longwaveDownCell, &
           oceanHeatCell, &
           freezingPotentialCell, &
           snowfallHeatCell, &
           latentHeatCell

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

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      real(kind=RKIND), pointer :: &
           dt, &
           earthArea

      logical, pointer :: &
           config_calc_surface_temperature, &
           config_AM_conservationCheck_write_to_logfile

      integer, pointer :: &
           nCellsSolve, &
           nHemispheres, &
           nAccumulate

      integer :: &
           iCell, &
           iHemisphere, &
           iSumPrev, &
           ierr

      integer :: &
           nVars

      integer :: &
           nSums

      real(kind=RKIND) :: &
           fluxScale

      logical, parameter :: &
           defaultLatentHeatFormulation = .true.

      err = 0

      call MPAS_pool_get_config(domain % blocklist % configs, "config_calc_surface_temperature", config_calc_surface_temperature)
      if (config_calc_surface_temperature) then
         nVars = 9
      else
         nVars = 5
      endif

      call MPAS_pool_get_config(domain % configs, "config_dt", dt)
      call MPAS_pool_get_config(domain % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                   config_AM_conservationCheck_write_to_logfile)

      !-------------------------------------------------------------
      ! Net heat flux to ice
      !-------------------------------------------------------------

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)
      nSums = nHemispheres * nVars

      allocate(sumArray(nSums))
      allocate(sumArrayOut(nSums))

      sumArray = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_dimension(blockPtr % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(blockPtr % structs, "icestate", icestatePool)
         call MPAS_pool_get_subpool(blockPtr % structs, "shortwave", shortwavePool)
         call MPAS_pool_get_subpool(blockPtr % structs, "ocean_fluxes", oceanFluxesPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "atmos_fluxes", atmosFluxesPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "atmos_coupling", atmosCouplingPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "diagnostics", diagnosticsPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)

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
         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)

         do iCell = 1, nCellsSolve

            oceanHeatCell         = -oceanHeatFluxArea(iCell) * areaCell(iCell)
            freezingPotentialCell = -max(0.0_RKIND, freezingMeltingPotentialInitial(iCell)) * areaCell(iCell)
            snowfallHeatCell      = -snowfallRate(iCell) * iceAreaCellInitial(iCell) * areaCell(iCell) * Lfresh

            do iHemisphere = 1, nHemispheres
               iSumPrev = (iHemisphere-1) * nVars
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  sumArray(iSumPrev+1) = sumArray(iSumPrev+1) + oceanHeatCell
                  sumArray(iSumPrev+2) = sumArray(iSumPrev+2) + freezingPotentialCell
                  sumArray(iSumPrev+3) = sumArray(iSumPrev+3) + snowfallHeatCell
               endif
            enddo ! iHemisphere

         enddo ! iCell

         if (defaultLatentHeatFormulation) then

            do iCell = 1, nCellsSolve

               latentHeatCell = latentHeatFlux(iCell) * iceAreaCell(iCell) * areaCell(iCell)

               do iHemisphere = 1, nHemispheres
                  iSumPrev = (iHemisphere-1) * nVars
                  if (cellInHemisphere(iHemisphere,iCell) == 1) then
                     sumArray(iSumPrev+4) = sumArray(iSumPrev+4) + latentHeatCell
                  endif
               enddo ! iHemisphere

            enddo ! iCell

         else

            do iCell = 1, nCellsSolve

               latentHeatCell = evaporativeWaterFlux(iCell) * iceAreaCell(iCell) * areaCell(iCell) * Lvap

               do iHemisphere = 1, nHemispheres
                  iSumPrev = (iHemisphere-1) * nVars
                  if (cellInHemisphere(iHemisphere,iCell) == 1) then
                     sumArray(iSumPrev+4) = sumArray(iSumPrev+4) + latentHeatCell
                  endif
               enddo ! iHemisphere

            enddo ! iCell

         endif

         ! surface heat flux
         if (config_calc_surface_temperature) then

            do iCell = 1, nCellsSolve

               absorbedShortwaveFluxCell = absorbedShortwaveFlux(iCell) * iceAreaCell(iCell) * areaCell(iCell)
               oceanShortwaveFluxCell    = -oceanShortwaveFlux(iCell)   * iceAreaCell(iCell) * areaCell(iCell)
               sensibleHeatFluxCell      = sensibleHeatFlux(iCell)      * iceAreaCell(iCell) * areaCell(iCell)
               longwaveUpCell            = longwaveUp(iCell)            * iceAreaCell(iCell) * areaCell(iCell)
               longwaveDownCell          = longwaveDown(iCell)          * iceAreaCellInitial(iCell) * areaCell(iCell)

               do iHemisphere = 1, nHemispheres
                  iSumPrev = (iHemisphere-1) * nVars
                  if (cellInHemisphere(iHemisphere,iCell) == 1) then
                     sumArray(iSumPrev+5) = sumArray(iSumPrev+5) + absorbedShortwaveFluxCell
                     sumArray(iSumPrev+6) = sumArray(iSumPrev+6) + oceanShortwaveFluxCell
                     sumArray(iSumPrev+7) = sumArray(iSumPrev+7) + sensibleHeatFluxCell
                     sumArray(iSumPrev+8) = sumArray(iSumPrev+8) + longwaveUpCell
                     sumArray(iSumPrev+9) = sumArray(iSumPrev+9) + longwaveDownCell
                  endif
               enddo ! iHemisphere

            enddo ! iCell

         else

            do iCell = 1, nCellsSolve

               surfaceHeatFluxCell = (surfaceHeatFlux(iCell) - latentHeatFlux(iCell)) * iceAreaCell(iCell) * areaCell(iCell)

               do iHemisphere = 1, nHemispheres
                  iSumPrev = (iHemisphere-1) * nVars
                  if (cellInHemisphere(iHemisphere,iCell) == 1) then
                     sumArray(iSumPrev+5) = sumArray(iSumPrev+5) + surfaceHeatFluxCell
                  endif
               enddo ! iHemisphere

            enddo ! iCell

         endif

         blockPtr => blockPtr % next
      enddo

      ! perform the sums over processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

      ! accumulate fluxes
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)

      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsOceanHeatFlux",     energyConsOceanHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsFreezingPotential", energyConsFreezingPotential)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsSnowfallHeat",      energyConsSnowfallHeat)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsLatentHeat",        energyConsLatentHeat)
      do iHemisphere = 1, nHemispheres
         iSumPrev = (iHemisphere-1) * nVars
         energyConsOceanHeatFlux    (iHemisphere) = energyConsOceanHeatFlux    (iHemisphere) + sumArrayOut(iSumPrev+1)
         energyConsFreezingPotential(iHemisphere) = energyConsFreezingPotential(iHemisphere) + sumArrayOut(iSumPrev+2)
         energyConsSnowfallHeat     (iHemisphere) = energyConsSnowfallHeat     (iHemisphere) + sumArrayOut(iSumPrev+3)
         energyConsLatentHeat       (iHemisphere) = energyConsLatentHeat       (iHemisphere) + sumArrayOut(iSumPrev+4)
      enddo ! iHemisphere

      if (config_calc_surface_temperature) then
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsAbsorbedShortwaveFlux", energyConsAbsorbedShortwaveFlux)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsOceanShortwaveFlux",    energyConsOceanShortwaveFlux)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsSensibleHeatFlux",      energyConsSensibleHeatFlux)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsLongwaveUp",            energyConsLongwaveUp)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsLongwaveDown",          energyConsLongwaveDown)
         do iHemisphere = 1, nHemispheres
            iSumPrev = (iHemisphere-1) * nVars
            energyConsAbsorbedShortwaveFlux(iHemisphere) = energyConsAbsorbedShortwaveFlux(iHemisphere) + sumArrayOut(iSumPrev+5)
            energyConsOceanShortwaveFlux   (iHemisphere) = energyConsOceanShortwaveFlux   (iHemisphere) + sumArrayOut(iSumPrev+6)
            energyConsSensibleHeatFlux     (iHemisphere) = energyConsSensibleHeatFlux     (iHemisphere) + sumArrayOut(iSumPrev+7)
            energyConsLongwaveUp           (iHemisphere) = energyConsLongwaveUp           (iHemisphere) + sumArrayOut(iSumPrev+8)
            energyConsLongwaveDown         (iHemisphere) = energyConsLongwaveDown         (iHemisphere) + sumArrayOut(iSumPrev+9)
         enddo ! iHemisphere
      else
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsSurfaceHeatFlux", energyConsSurfaceHeatFlux)
         do iHemisphere = 1, nHemispheres
            iSumPrev = (iHemisphere-1) * nVars
            energyConsSurfaceHeatFlux(iHemisphere) = energyConsSurfaceHeatFlux(iHemisphere) + sumArrayOut(iSumPrev+5)
         enddo ! iHemisphere
      endif

      ! cleanup
      deallocate(sumArray)
      deallocate(sumArrayOut)

      !-------------------------------------------------------------
      ! Energy conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! convert fluxes to fluxes per m2 of whole earth and per second
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
         call MPAS_pool_get_array(conservationCheckAMPool, "nAccumulate", nAccumulate)
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAreaAM", conservationCheckAreaAMPool)
         call MPAS_pool_get_array(conservationCheckAreaAMPool, "earthArea", earthArea)
         fluxScale = 1.0_RKIND / (earthArea * real(nAccumulate))

         ! get initial energy
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "initialEnergy", initialEnergy)

         ! get final energy
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "finalEnergy", finalEnergy)
         call compute_total_energy(domain, finalEnergy)

         ! compute the energy change
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyChange", energyChange)
         energyChange(:) = finalEnergy(:) - initialEnergy(:)
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyChangeFlux", energyChangeFlux)
         energyChangeFlux(:) = energyChange(:) * (fluxScale / dt)

         ! surface heat flux
         if (config_calc_surface_temperature) then
            call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsSurfaceHeatFlux", energyConsSurfaceHeatFlux)
            energyConsSurfaceHeatFlux(:) = &
                 energyConsAbsorbedShortwaveFlux(:) + &
                 energyConsOceanShortwaveFlux(:) + &
                 energyConsSensibleHeatFlux(:) + &
                 energyConsLongwaveUp(:) + &
                 energyConsLongwaveDown(:)
         endif

         ! calculate the final net energy flux to the ice
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "netEnergyFlux", netEnergyFlux)
         netEnergyFlux(:) = &
               energyConsSurfaceHeatFlux(:) &
             + energyConsSnowfallHeat(:) &
             + energyConsLatentHeat(:) &
             + energyConsOceanHeatFlux(:) &
             + energyConsFreezingPotential(:)

         ! absolute error
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "absoluteEnergyError", absoluteEnergyError)
         absoluteEnergyError(:) = netEnergyFlux(:) * dt - energyChange(:)

         ! rescale fluxes
         energyConsSurfaceHeatFlux(:)       = energyConsSurfaceHeatFlux(:)       * fluxScale
         energyConsAbsorbedShortwaveFlux(:) = energyConsAbsorbedShortwaveFlux(:) * fluxScale
         energyConsOceanShortwaveFlux(:)    = energyConsOceanShortwaveFlux(:)    * fluxScale
         energyConsSensibleHeatFlux(:)      = energyConsSensibleHeatFlux(:)      * fluxScale
         energyConsLongwaveUp(:)            = energyConsLongwaveUp(:)            * fluxScale
         energyConsLongwaveDown(:)          = energyConsLongwaveDown(:)          * fluxScale
         energyConsSnowfallHeat(:)          = energyConsSnowfallHeat(:)          * fluxScale
         energyConsLatentHeat(:)            = energyConsLatentHeat(:)            * fluxScale
         energyConsOceanHeatFlux(:)         = energyConsOceanHeatFlux(:)         * fluxScale
         energyConsFreezingPotential(:)     = energyConsFreezingPotential(:)     * fluxScale
         netEnergyFlux(:)                   = netEnergyFlux(:)                   * fluxScale

         ! relative flux error
         call MPAS_pool_get_array(conservationCheckEnergyAMPool, "relativeEnergyError", relativeEnergyError)
         do iHemisphere = 1, nHemispheres
            if (abs(finalEnergy(iHemisphere)) > 0.0) then
               relativeEnergyError(iHemisphere) = absoluteEnergyError(iHemisphere) / finalEnergy(iHemisphere)
            else
               relativeEnergyError(iHemisphere) = 0.0_RKIND
            endif
         enddo ! iHemisphere

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            call mpas_log_write('-----------------------------------------------------------------------------------')
            call mpas_log_write(' Energy conservation check')
            call mpas_log_write(' ')
            call mpas_log_write(' Initial energy                 (J) = '//trim(hemisphere_format(initialEnergy)))
            call mpas_log_write(' Final energy                   (J) = '//trim(hemisphere_format(finalEnergy)))
            call mpas_log_write(' Energy change                  (J) = '//trim(hemisphere_format(energyChange)))
            call mpas_log_write(' Energy change flux          (W/m2) = '//trim(hemisphere_format(energyChangeFlux)))
            call mpas_log_write(' ')
            call mpas_log_write(' Surface heat flux           (W/m2) = '//trim(hemisphere_format(energyConsSurfaceHeatFlux)))
            call mpas_log_write('    Absorbed shortwave flux  (W/m2) = '//trim(hemisphere_format(energyConsAbsorbedShortwaveFlux)))
            call mpas_log_write('    Ocean Shortwave flux     (W/m2) = '//trim(hemisphere_format(energyConsOceanShortwaveFlux)))
            call mpas_log_write('    Sensible heat flux       (W/m2) = '//trim(hemisphere_format(energyConsSensibleHeatFlux)))
            call mpas_log_write('    Longwave up              (W/m2) = '//trim(hemisphere_format(energyConsLongwaveUp)))
            call mpas_log_write('    Longwave down            (W/m2) = '//trim(hemisphere_format(energyConsLongwaveDown)))
            call mpas_log_write(' Ocean heat flux             (W/m2) = '//trim(hemisphere_format(energyConsOceanHeatFlux)))
            call mpas_log_write(' Freezing heat flux          (W/m2) = '//trim(hemisphere_format(energyConsFreezingPotential)))
            call mpas_log_write(' Snowfall heat flux          (W/m2) = '//trim(hemisphere_format(energyConsSnowfallHeat)))
            call mpas_log_write(' Latent heat flux            (W/m2) = '//trim(hemisphere_format(energyConsLatentHeat)))
            call mpas_log_write(' ')
            call mpas_log_write(' Net energy change              (J) = '//trim(hemisphere_format((netEnergyFlux * dt) / fluxScale)))
            call mpas_log_write(' Net energy flux             (W/m2) = '//trim(hemisphere_format(netEnergyFlux)))
            call mpas_log_write(' ')
            call mpas_log_write(' Absolute energy error          (J) = '//trim(hemisphere_format(absoluteEnergyError)))
            call mpas_log_write(' Absolute energy error       (W/m2) = '//trim(hemisphere_format((absoluteEnergyError * fluxScale) / dt)))
            call mpas_log_write(' Relative energy error              = '//trim(hemisphere_format(relativeEnergyError)))

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
           blockPtr

      type(MPAS_pool_type), pointer :: &
           conservationCheckMassAMPool, &
           conservationCheckAMPool, &
           conservationCheckAreaAMPool

      real(kind=RKIND), dimension(:), pointer :: &
           initialMass, &
           finalMass, &
           massChange, &
           massChangeFlux, &
           netMassFlux, &
           absoluteMassError, &
           relativeMassError

      real(kind=RKIND), dimension(:), pointer :: &
           massConsRainfallRate, &
           massConsSnowfallRate, &
           massConsEvaporation, &
           massConsFreshWater, &
           massConsFrazilWater

      real(kind=RKIND) :: &
           rainfallCell, &
           snowfallCell, &
           evaporativeWaterCell, &
           oceanFreshWaterCell, &
           frazilCell

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

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      real(kind=RKIND), pointer :: &
           dt, &
           earthArea

      logical, pointer :: &
           config_update_ocean_fluxes, &
           config_AM_conservationCheck_write_to_logfile

      character(len=strKIND), pointer :: &
           config_thermodynamics_type

      integer, pointer :: &
           nCellsSolve, &
           nHemispheres, &
           nAccumulate

      integer :: &
           iCell, &
           iHemisphere, &
           iSumPrev, &
           ierr

      integer, parameter :: &
           nVars = 5

      integer :: &
           nSums

      real(kind=RKIND) :: &
           fluxScale

      err = 0

      call MPAS_pool_get_config(domain % configs, "config_dt", dt)
      call MPAS_pool_get_config(domain % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                   config_AM_conservationCheck_write_to_logfile)

      !-------------------------------------------------------------
      ! Net mass flux to ice
      !-------------------------------------------------------------

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)
      nSums = nHemispheres * nVars

      allocate(sumArray(nSums))
      allocate(sumArrayOut(nSums))

      sumArray = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_config(blockPtr % configs, "config_update_ocean_fluxes", config_update_ocean_fluxes)
         call MPAS_pool_get_config(blockPtr % configs, "config_thermodynamics_type", config_thermodynamics_type)

         call MPAS_pool_get_dimension(blockPtr % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(blockPtr % structs, "icestate", icestatePool)
         call MPAS_pool_get_subpool(blockPtr % structs, "atmos_coupling", atmosCouplingPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "atmos_fluxes", atmosFluxesPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "ocean_fluxes", oceanFluxesPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "melt_growth_rates", meltGrowthRatesPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
         call MPAS_pool_get_array(icestatePool, "iceAreaCellInitial", iceAreaCellInitial)
         call MPAS_pool_get_array(atmosCouplingPool, "rainfallRate", rainfallRate)
         call MPAS_pool_get_array(atmosCouplingPool, "snowfallRate", snowfallRate)
         call MPAS_pool_get_array(atmosFluxesPool, "evaporativeWaterFlux", evaporativeWaterFlux)
         call MPAS_pool_get_array(oceanFluxesPool, "oceanFreshWaterFluxArea", oceanFreshWaterFluxArea)
         call MPAS_pool_get_array(meltGrowthRatesPool, "frazilFormation", frazilFormation)
         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)

         do iCell = 1, nCellsSolve

            rainfallCell         = rainfallRate(iCell) * iceAreaCellInitial(iCell) * areaCell(iCell)
            snowfallCell         = snowfallRate(iCell) * iceAreaCellInitial(iCell) * areaCell(iCell)
            evaporativeWaterCell = evaporativeWaterFlux(iCell) * iceAreaCell(iCell) * areaCell(iCell)
            oceanFreshWaterCell  = -oceanFreshWaterFluxArea(iCell) * areaCell(iCell)

            do iHemisphere = 1, nHemispheres
               iSumPrev = (iHemisphere-1) * nVars
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  sumArray(iSumPrev+1) = sumArray(iSumPrev+1) + rainfallCell
                  sumArray(iSumPrev+2) = sumArray(iSumPrev+2) + snowfallCell
                  sumArray(iSumPrev+3) = sumArray(iSumPrev+3) + evaporativeWaterCell
                  sumArray(iSumPrev+4) = sumArray(iSumPrev+4) + oceanFreshWaterCell
               endif
            enddo ! iHemisphere

         enddo ! iCell

         ! frazil ice
         if (config_update_ocean_fluxes .and. trim(config_thermodynamics_type) == "mushy") then

            do iCell = 1, nCellsSolve

               frazilCell = (frazilFormation(iCell) * areaCell(iCell) * rhoi) / dt
               do iHemisphere = 1, nHemispheres
                  iSumPrev = (iHemisphere-1) * nVars
                  if (cellInHemisphere(iHemisphere,iCell) == 1) then
                     sumArray(iSumPrev+5) = sumArray(iSumPrev+5) + frazilCell
                  endif
               enddo ! iHemisphere

            enddo ! iCell

         endif

         blockPtr => blockPtr % next
      enddo

      ! perform the sums over processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

      ! accumulate fluxes
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckMassAM", conservationCheckMassAMPool)

      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsRainfallRate", massConsRainfallRate)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsSnowfallRate", massConsSnowfallRate)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsEvaporation", massConsEvaporation)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsFreshWater", massConsFreshWater)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsFrazilWater", massConsFrazilWater)

      do iHemisphere = 1, nHemispheres
         iSumPrev = (iHemisphere-1) * nVars
         massConsRainfallRate(iHemisphere) = massConsRainfallRate(iHemisphere) + sumArrayOut(iSumPrev+1)
         massConsSnowfallRate(iHemisphere) = massConsSnowfallRate(iHemisphere) + sumArrayOut(iSumPrev+2)
         massConsEvaporation (iHemisphere) = massConsEvaporation (iHemisphere) + sumArrayOut(iSumPrev+3)
         massConsFreshWater  (iHemisphere) = massConsFreshWater  (iHemisphere) + sumArrayOut(iSumPrev+4)
         massConsFrazilWater (iHemisphere) = massConsFrazilWater (iHemisphere) + sumArrayOut(iSumPrev+5)
      enddo ! iHemisphere

      ! cleanup
      deallocate(sumArray)
      deallocate(sumArrayOut)

      !-------------------------------------------------------------
      ! Mass conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! convert to fluxes per m2 of whole earth and per second
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
         call MPAS_pool_get_array(conservationCheckAMPool, "nAccumulate", nAccumulate)
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAreaAM", conservationCheckAreaAMPool)
         call MPAS_pool_get_array(conservationCheckAreaAMPool, "earthArea", earthArea)
         fluxScale = 1.0_RKIND / (earthArea * real(nAccumulate))

         ! get initial mass
         call MPAS_pool_get_array(conservationCheckMassAMPool, "initialMass", initialMass)

         ! get final mass
         call MPAS_pool_get_array(conservationCheckMassAMPool, "finalMass", finalMass)
         call compute_total_mass(domain, finalMass)

         ! compute the energy change
         call MPAS_pool_get_array(conservationCheckMassAMPool, "massChange", massChange)
         massChange(:) = finalMass(:) - initialMass(:)
         call MPAS_pool_get_array(conservationCheckMassAMPool, "massChangeFlux", massChangeFlux)
         massChangeFlux(:) = massChange(:) * (fluxScale / dt)

         ! calculate the final net energy flux to the ice
         call MPAS_pool_get_array(conservationCheckMassAMPool, "netMassFlux", netMassFlux)
         netMassFlux(:) = &
                massConsRainfallRate(:) &
              + massConsSnowfallRate(:) &
              + massConsEvaporation(:) &
              + massConsFreshWater(:) &
              + massConsFrazilWater(:)

         ! absolute error
         call MPAS_pool_get_array(conservationCheckMassAMPool, "absoluteMassError", absoluteMassError)
         absoluteMassError(:) = netMassFlux(:) * dt - massChange(:)

         ! rescale fluxes
         massConsRainfallRate(:) = massConsRainfallRate(:) * fluxScale
         massConsSnowfallRate(:) = massConsSnowfallRate(:) * fluxScale
         massConsEvaporation(:)  = massConsEvaporation(:)  * fluxScale
         massConsFreshWater(:)   = massConsFreshWater(:)   * fluxScale
         massConsFrazilWater(:)  = massConsFrazilWater(:)  * fluxScale
         netMassFlux(:)          = netMassFlux(:)          * fluxScale

         ! compute the final energy error
         call MPAS_pool_get_array(conservationCheckMassAMPool, "relativeMassError", relativeMassError)
         do iHemisphere = 1, nHemispheres
            if (abs(finalMass(iHemisphere)) > 0.0) then
               relativeMassError(iHemisphere) = absoluteMassError(iHemisphere) / finalMass(iHemisphere)
            else
               relativeMassError(iHemisphere) = 0.0_RKIND
            endif
         enddo ! iHemisphere

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            call mpas_log_write('-----------------------------------------------------------------------------------')
            call mpas_log_write(' Mass conservation check')
            call mpas_log_write(' ')
            call mpas_log_write(' Initial mass                  (kg) = '//trim(hemisphere_format(initialMass)))
            call mpas_log_write(' Final mass                    (kg) = '//trim(hemisphere_format(finalMass)))
            call mpas_log_write(' Mass change                   (kg) = '//trim(hemisphere_format(massChange)))
            call mpas_log_write(' Mass change flux          (kg/m2s) = '//trim(hemisphere_format(massChangeFlux)))
            call mpas_log_write(' ')
            call mpas_log_write(' Rainfall mass flux        (kg/m2s) = '//trim(hemisphere_format(massConsRainfallRate)))
            call mpas_log_write(' Snowfall mass flux        (kg/m2s) = '//trim(hemisphere_format(massConsSnowfallRate)))
            call mpas_log_write(' Evaporative mass flux     (kg/m2s) = '//trim(hemisphere_format(massConsEvaporation)))
            call mpas_log_write(' Fresh water mass flux     (kg/m2s) = '//trim(hemisphere_format(massConsFreshWater)))
            call mpas_log_write(' Frazil water flux         (kg/m2s) = '//trim(hemisphere_format(massConsFrazilWater)))
            call mpas_log_write(' ')
            call mpas_log_write(' Net mass change               (kg) = '//trim(hemisphere_format((netMassFlux * dt) / fluxScale)))
            call mpas_log_write(' Net mass flux             (kg/m2s) = '//trim(hemisphere_format(netMassFlux)))
            call mpas_log_write(' ')
            call mpas_log_write(' Absolute mass error           (kg) = '//trim(hemisphere_format(absoluteMassError)))
            call mpas_log_write(' Absolute mass error       (kg/m2s) = '//trim(hemisphere_format((absoluteMassError * fluxScale) / dt)))
            call mpas_log_write(' Relative mass error                = '//trim(hemisphere_format(relativeMassError)))

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
           blockPtr

      type(MPAS_pool_type), pointer :: &
           conservationCheckSaltAMPool, &
           conservationCheckAMPool, &
           conservationCheckAreaAMPool

      real(kind=RKIND), dimension(:), pointer :: &
           initialSalt, &
           finalSalt, &
           saltChange, &
           saltChangeFlux, &
           netSaltFlux, &
           absoluteSaltError, &
           relativeSaltError

      real(kind=RKIND), dimension(:), pointer :: &
           saltConsOceanSaltFlux, &
           saltConsFrazilSaltFlux

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

      real(kind=RKIND) :: &
           frazilCell, &
           oceanSaltCell

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      real(kind=RKIND), pointer :: &
           dt, &
           earthArea

      logical, pointer :: &
           config_update_ocean_fluxes, &
           config_AM_conservationCheck_write_to_logfile

      character(len=strKIND), pointer :: &
           config_thermodynamics_type

      integer, pointer :: &
           nCellsSolve, &
           nHemispheres, &
           nAccumulate

      integer :: &
           iCell, &
           iHemisphere, &
           iSumPrev, &
           ierr

      integer, parameter :: &
           nVars = 2

      integer :: &
           nSums

      real(kind=RKIND) :: &
           fluxScale

      err = 0

      call MPAS_pool_get_config(domain % configs, "config_dt", dt)
      call MPAS_pool_get_config(domain % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                   config_AM_conservationCheck_write_to_logfile)

      !-------------------------------------------------------------
      ! Net salt flux to ice
      !-------------------------------------------------------------

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)
      nSums = nHemispheres * nVars

      allocate(sumArray(nSums))
      allocate(sumArrayOut(nSums))

      sumArray = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_config(blockPtr % configs, "config_update_ocean_fluxes", config_update_ocean_fluxes)
         call MPAS_pool_get_config(blockPtr % configs, "config_thermodynamics_type", config_thermodynamics_type)

         call MPAS_pool_get_dimension(blockPtr % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "ocean_fluxes", oceanFluxesPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "melt_growth_rates", meltGrowthRatesPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(oceanFluxesPool, "oceanSaltFluxArea", oceanSaltFluxArea)
         call MPAS_pool_get_array(meltGrowthRatesPool, "frazilFormation", frazilFormation)
         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)

         do iCell = 1, nCellsSolve

            oceanSaltCell = -oceanSaltFluxArea(iCell) * areaCell(iCell)

            do iHemisphere = 1, nHemispheres
               iSumPrev = (iHemisphere-1) * nVars
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  sumArray(iSumPrev+1) = sumArray(iSumPrev+1) + oceanSaltCell
               endif
            enddo ! iHemisphere

         enddo ! iCell

         ! frazil
         if (config_update_ocean_fluxes .and. trim(config_thermodynamics_type) == "mushy") then

            do iCell = 1, nCellsSolve

               frazilCell = (frazilFormation(iCell) * areaCell(iCell) * rhoi * ice_ref_salinity * 0.001_RKIND) / dt
               do iHemisphere = 1, nHemispheres
                  iSumPrev = (iHemisphere-1) * nVars
                  if (cellInHemisphere(iHemisphere,iCell) == 1) then
                     sumArray(iSumPrev+2) = sumArray(iSumPrev+2) + frazilCell
                  endif
               enddo ! iHemisphere


            enddo ! iCell

         endif

         blockPtr => blockPtr % next
      enddo

      ! perform the sums over processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

      ! accumulate fluxes
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)

      call MPAS_pool_get_array(conservationCheckSaltAMPool, "saltConsOceanSaltFlux", saltConsOceanSaltFlux)
      call MPAS_pool_get_array(conservationCheckSaltAMPool, "saltConsFrazilSaltFlux", saltConsFrazilSaltFlux)

      do iHemisphere = 1, nHemispheres
         iSumPrev = (iHemisphere-1) * nVars
         saltConsOceanSaltFlux (iHemisphere) = saltConsOceanSaltFlux (iHemisphere) + sumArrayOut(iSumPrev+1)
         saltConsFrazilSaltFlux(iHemisphere) = saltConsFrazilSaltFlux(iHemisphere) + sumArrayOut(iSumPrev+2)
      enddo ! iHemisphere

      ! cleanup
      deallocate(sumArray)
      deallocate(sumArrayOut)

      !-------------------------------------------------------------
      ! Salt conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! convert fluxes to fluxes per m2 of whole earth and per second
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
         call MPAS_pool_get_array(conservationCheckAMPool, "nAccumulate", nAccumulate)
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAreaAM", conservationCheckAreaAMPool)
         call MPAS_pool_get_array(conservationCheckAreaAMPool, "earthArea", earthArea)
         fluxScale = 1.0_RKIND / (earthArea * real(nAccumulate))

         ! get initial salt content
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "initialSalt", initialSalt)

         ! get final salt content
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "finalSalt", finalSalt)
         call compute_total_salt(domain, finalSalt)

         ! compute the salt content change
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "saltChange", saltChange)
         saltChange(:) = finalSalt(:) - initialSalt(:)
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "saltChangeFlux", saltChangeFlux)
         saltChangeFlux(:) = saltChange(:) * (fluxScale / dt)

         ! calculate the final net salt flux to the ice
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "netSaltFlux", netSaltFlux)
         netSaltFlux(:) = &
              + saltConsOceanSaltFlux(:) &
              + saltConsFrazilSaltFlux(:)

         ! absolute error
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "absoluteSaltError", absoluteSaltError)
         absoluteSaltError(:) = netSaltFlux(:) * dt - saltChange(:)

         ! rescale fluxes
         saltConsOceanSaltFlux(:)  = saltConsOceanSaltFlux(:)  * fluxScale
         saltConsFrazilSaltFlux(:) = saltConsFrazilSaltFlux(:) * fluxScale
         netSaltFlux(:)               = netSaltFlux(:)               * fluxScale

         ! compute the final salt error
         call MPAS_pool_get_array(conservationCheckSaltAMPool, "relativeSaltError", relativeSaltError)
         do iHemisphere = 1, nHemispheres
            if (abs(finalSalt(iHemisphere)) > 0.0) then
               relativeSaltError(iHemisphere) = absoluteSaltError(iHemisphere) / finalSalt(iHemisphere)
            else
               relativeSaltError(iHemisphere) = 0.0_RKIND
            endif
         enddo ! iHemisphere

         !-------------------------------------------------------------
         ! Output to log file
         !-------------------------------------------------------------

         if (config_AM_conservationCheck_write_to_logfile) then

            call mpas_log_write('-----------------------------------------------------------------------------------')
            call mpas_log_write(' Salt conservation check')
            call mpas_log_write(' ')
            call mpas_log_write(' Initial salt                  (kg) = '//trim(hemisphere_format(initialSalt)))
            call mpas_log_write(' Final salt                    (kg) = '//trim(hemisphere_format(finalSalt)))
            call mpas_log_write(' Salt change                   (kg) = '//trim(hemisphere_format(saltChange)))
            call mpas_log_write(' Salt change flux          (kg/m2s) = '//trim(hemisphere_format(saltChangeFlux)))
            call mpas_log_write(' ')
            call mpas_log_write(' Ocean salt flux           (kg/m2s) = '//trim(hemisphere_format(saltConsOceanSaltFlux)))
            call mpas_log_write(' Frazil salt flux          (kg/m2s) = '//trim(hemisphere_format(saltConsFrazilSaltFlux)))
            call mpas_log_write(' ')
            call mpas_log_write(' Net salt change               (kg) = '//trim(hemisphere_format((netSaltFlux * dt) / fluxScale)))
            call mpas_log_write(' Net salt flux             (kg/m2s) = '//trim(hemisphere_format(netSaltFlux)))
            call mpas_log_write(' ')
            call mpas_log_write(' Absolute mass error           (kg) = '//trim(hemisphere_format(absoluteSaltError)))
            call mpas_log_write(' Absolute mass error       (kg/m2s) = '//trim(hemisphere_format((absoluteSaltError * fluxScale) / dt)))
            call mpas_log_write(' Relative salt error                = '//trim(hemisphere_format(relativeSaltError)))

         endif

      endif

    end subroutine salt_conservation

!***********************************************************************
!
!  routine carbon_conservation
!
!> \brief   Compute MPAS-Seaice analysis member
!> \author  Nicole Jeffery
!> \date    27 May 2020
!> \details
!>  This routine conducts all computations to verify carbon conservation
!>  in seaice.
!
!-----------------------------------------------------------------------

   subroutine carbon_conservation(domain, err)

      use seaice_constants, only: &
           gramsCarbonPerMolCarbon

      type (domain_type), intent(inout) :: &
           domain

      integer, intent(out) :: &
           err !< Output: error flag

      integer, pointer :: &
           maxAlgaeType, &
           maxDOCType, &
           maxDICType, &
           maxDONType, &
           maxIronType

      type(block_type), pointer :: &
           blockPtr

      type(MPAS_pool_type), pointer :: &
           conservationCheckCarbonAMPool, &
           conservationCheckAMPool, &
           conservationCheckAreaAMPool

      real(kind=RKIND), dimension(:), pointer :: &
           initialCarbon, &
           finalCarbon, &
           carbonChange, &
           carbonChangeFlux, &
           netCarbonFlux, &
           absoluteCarbonError, &
           relativeCarbonError, &
           iceAreaCell, &
           accumAbsoluteCarbonError, &
           accumRelativeCarbonError

      real(kind=RKIND), dimension(:), pointer :: &
           carbonConsOceanCarbonFlux, &
           carbonConsOceanCarbonFluxCheck, &
           carbonConsDiatomFlux, &
           carbonConsSmallAlgaeFlux, &
           carbonConsDOC1Flux, &
           carbonConsDOC2Flux, &
           carbonConsDICFlux, &
           carbonConsTOCFlux, &
           carbonConsDONCarbonFlux, &
           carbonConsHumicsFlux, &
           carbonConsSemiLabileDOCFlux

      real(kind=RKIND), pointer :: &
           config_ratio_C_to_N_diatoms, &
           config_ratio_C_to_N_small_plankton, &
           config_ratio_C_to_N_phaeocystis, &
           config_ratio_C_to_N_proteins

      real(kind=RKIND) :: &
           totalOceanCarbonCell, &
           totalOceanCarbonCellCheck, &
           totalOceanDiatomCell, &
           totalOceanSmallAlgaeCell, &
           totalOceanDOC1Cell, &
           totalOceanDOC2Cell, &
           totalOceanDONCarbonCell, &
           totalOceanDICCell, &
           totalOceanHumicsCell, &
           totalOceanSemiLabileDOCCell

      real(kind=RKIND), dimension(:), allocatable :: &
           sumArray, &
           sumArrayOut

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           tracersAggregatePool, &
           biogeochemistryPool

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           totalOceanCarbonFlux, &
           oceanHumicsFlux

      real(kind=RKIND), dimension(:,:), pointer :: &
           oceanBioFluxes, &
           oceanAlgaeFlux, &
           oceanDOCFlux, &
           oceanDICFlux, &
           oceanDONFlux

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      real(kind=RKIND), pointer :: &
           dt, &
           earthArea, &
           relativeErrorBounds

      logical, pointer :: &
           config_AM_conservationCheck_write_to_logfile, &
           config_AM_conservationCheck_carbon_failure_abort

      integer, pointer :: &
           nCellsSolve, &
           nHemispheres, &
           nAccumulate

      type (MPAS_Time_Type) :: &
           currTime, &
           startTime

      integer :: &
           iCell, &
           iHemisphere, &
           iSumPrev, &
           ierr

      integer :: &
           nVars, &
           nZBGC, &
           nSums

      real(kind=RKIND) :: &
           fluxScale, &
           relativeErrorStepBound

      character(len=strKIND) :: &
           timeStr, &
           timeStrStart

      character(len=16) :: &
           valStr

      err = 0
      relativeErrorStepBound = 1.0E-7 ! mass conservation relative error bound 1e-5

      call MPAS_pool_get_config(domain % configs, "config_dt", dt)
      call MPAS_pool_get_config(domain % configs, "config_AM_conservationCheck_write_to_logfile", &
                                                   config_AM_conservationCheck_write_to_logfile)
      call MPAS_pool_get_config(domain % configs, "config_AM_conservationCheck_carbon_failure_abort", &
                                                   config_AM_conservationCheck_carbon_failure_abort)

      !-------------------------------------------------------------
      ! Net carbon flux to ice
      !-------------------------------------------------------------

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)

      nVars = 10

      nSums = nHemispheres * nVars

      allocate(sumArray(nSums))
      allocate(sumArrayOut(nSums))

      sumArray = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_dimension(blockPtr % dimensions, "nCellsSolve", nCellsSolve)

         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "biogeochemistry", biogeochemistryPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)

         call MPAS_pool_get_config(blockPtr % configs, "config_ratio_C_to_N_diatoms", config_ratio_C_to_N_diatoms)
         call MPAS_pool_get_config(blockPtr % configs, "config_ratio_C_to_N_small_plankton", config_ratio_C_to_N_small_plankton)
         call MPAS_pool_get_config(blockPtr % configs, "config_ratio_C_to_N_phaeocystis", config_ratio_C_to_N_phaeocystis)
         call MPAS_pool_get_config(blockPtr % configs, "config_ratio_C_to_N_proteins", config_ratio_C_to_N_proteins)
         call MPAS_pool_get_dimension(blockPtr % dimensions, "maxAlgaeType", maxAlgaeType)
         call MPAS_pool_get_dimension(blockPtr % dimensions, "maxDOCType", maxDOCType)
         call MPAS_pool_get_dimension(blockPtr % dimensions, "maxDICType", maxDICType)
         call MPAS_pool_get_dimension(blockPtr % dimensions, "maxDONType", maxDONType)
         call MPAS_pool_get_dimension(blockPtr % dimensions, "maxIronType", maxIronType)


         call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(biogeochemistryPool, "totalOceanCarbonFlux", totalOceanCarbonFlux)
         call MPAS_pool_get_array(biogeochemistryPool, "oceanBioFluxes", oceanBioFluxes)
         call MPAS_pool_get_array(biogeochemistryPool, "oceanAlgaeFlux", oceanAlgaeFlux)
         call MPAS_pool_get_array(biogeochemistryPool, "oceanDOCFlux", oceanDOCFlux)
         call MPAS_pool_get_array(biogeochemistryPool, "oceanDICFlux", oceanDICFlux)
         call MPAS_pool_get_array(biogeochemistryPool, "oceanDONFlux", oceanDONFlux)
         call MPAS_pool_get_array(biogeochemistryPool, "oceanHumicsFlux", oceanHumicsFlux)
         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)

         do iCell = 1, nCellsSolve

            ! compute total carbon flux to ocean from individual fluxes

            nZBGC = 1
            totalOceanDiatomCell = oceanAlgaeFlux(nZBGC,iCell) * areaCell(iCell) * &
               config_ratio_C_to_N_diatoms * gramsCarbonPerMolCarbon * 1.0e-6_RKIND * &
               iceAreaCell(iCell)
            nZBGC = nZBGC+1
            totalOceanSmallAlgaeCell = oceanAlgaeFlux(nZBGC,iCell) * areaCell(iCell) * &
               config_ratio_C_to_N_small_plankton * gramsCarbonPerMolCarbon * 1.0e-6_RKIND * &
               iceAreaCell(iCell)
            nZBGC = 1
            totalOceanDOC1Cell = oceanDOCFlux(nZBGC,iCell) * areaCell(iCell) * &
               gramsCarbonPerMolCarbon * 1.0e-6_RKIND * &
               iceAreaCell(iCell)
            nZBGC = nZBGC+1
            totalOceanDOC2Cell = oceanDOCFlux(nZBGC,iCell) * areaCell(iCell) * &
               gramsCarbonPerMolCarbon * 1.0e-6_RKIND * &
               iceAreaCell(iCell)

            nZBGC = 1
            totalOceanDICCell = oceanDICFlux(nZBGC,iCell) * areaCell(iCell) * &
               gramsCarbonPerMolCarbon * 1.0e-6_RKIND * &
               iceAreaCell(iCell)

            nZBGC = 1
            totalOceanDONCarbonCell = oceanDONFlux(nZBGC,iCell) * areaCell(iCell) * &
               config_ratio_C_to_N_proteins * gramsCarbonPerMolCarbon * 1.0e-6_RKIND * &
               iceAreaCell(iCell)

            totalOceanSemiLabileDOCCell = totalOceanDONCarbonCell + totalOceanDOC1Cell + &
               totalOceanDOC2Cell

            totalOceanHumicsCell = oceanHumicsFlux(iCell) * areaCell(iCell) * &
               gramsCarbonPerMolCarbon * 1.0e-6_RKIND * &
               iceAreaCell(iCell)

            totalOceanCarbonCellCheck = totalOceanHumicsCell + totalOceanDONCarbonCell + &
               totalOceanDOC2Cell + totalOceanDOC1Cell + totalOceanDiatomCell + totalOceanSmallAlgaeCell + &
               totalOceanDICCell

            ! carbon flux to ocean computed in code
            totalOceanCarbonCell = totalOceanCarbonFlux(iCell) * areaCell(iCell) * &
               gramsCarbonPerMolCarbon * 1.0e-6_RKIND

            do iHemisphere = 1, nHemispheres
               iSumPrev = (iHemisphere-1) * nVars
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  sumArray(iSumPrev+1) = sumArray(iSumPrev+1) + totalOceanCarbonCell
                  sumArray(iSumPrev+2) = sumArray(iSumPrev+2) + totalOceanCarbonCellCheck
                  sumArray(iSumPrev+3) = sumArray(iSumPrev+3) + totalOceanDiatomCell
                  sumArray(iSumPrev+4) = sumArray(iSumPrev+4) + totalOceanSmallAlgaeCell
                  sumArray(iSumPrev+5) = sumArray(iSumPrev+5) + totalOceanDOC1Cell
                  sumArray(iSumPrev+6) = sumArray(iSumPrev+6) + totalOceanDOC2Cell
                  sumArray(iSumPrev+7) = sumArray(iSumPrev+7) + totalOceanDICCell
                  sumArray(iSumPrev+8) = sumArray(iSumPrev+8) + totalOceanDONCarbonCell
                  sumArray(iSumPrev+9) = sumArray(iSumPrev+9) + totalOceanHumicsCell
                  sumArray(iSumPrev+10)= sumArray(iSumPrev+10) + totalOceanSemiLabileDOCCell
               endif
            enddo ! iHemisphere

         enddo ! iCell

         blockPtr => blockPtr % next
      enddo

      ! perform the sums over processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nSums, sumArray, sumArrayOut)

      ! accumulate fluxes
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckCarbonAM", conservationCheckCarbonAMPool)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsOceanCarbonFlux", carbonConsOceanCarbonFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsOceanCarbonFluxCheck", carbonConsOceanCarbonFluxCheck)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDiatomFlux", carbonConsDiatomFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsSmallAlgaeFlux", carbonConsSmallAlgaeFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDOC1Flux", carbonConsDOC1Flux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDOC2Flux", carbonConsDOC2Flux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDONCarbonFlux", carbonConsDONCarbonFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDICFlux", carbonConsDICFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsTOCFlux", carbonConsTOCFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsHumicsFlux", carbonConsHumicsFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsSemiLabileDOCFlux", carbonConsSemiLabileDOCFlux)

      do iHemisphere = 1, nHemispheres
         iSumPrev = (iHemisphere-1) * nVars
         carbonConsOceanCarbonFlux(iHemisphere) = carbonConsOceanCarbonFlux(iHemisphere) + sumArrayOut(iSumPrev+1)
         carbonConsOceanCarbonFluxCheck(iHemisphere) = carbonConsOceanCarbonFluxCheck(iHemisphere) + sumArrayOut(iSumPrev+2)
         carbonConsDiatomFlux(iHemisphere) = carbonConsDiatomFlux(iHemisphere) + sumArrayOut(iSumPrev+3)
         carbonConsSmallAlgaeFlux(iHemisphere) = carbonConsSmallAlgaeFlux(iHemisphere) + sumArrayOut(iSumPrev+4)
         carbonConsDOC1Flux(iHemisphere) = carbonConsDOC1Flux(iHemisphere) + sumArrayOut(iSumPrev+5)
         carbonConsDOC2Flux(iHemisphere) = carbonConsDOC2Flux(iHemisphere) + sumArrayOut(iSumPrev+6)
         carbonConsDICFlux(iHemisphere) = carbonConsDICFlux(iHemisphere) + sumArrayOut(iSumPrev+7)
         carbonConsDONCarbonFlux(iHemisphere) = carbonConsDONCarbonFlux(iHemisphere) + sumArrayOut(iSumPrev+8)
         carbonConsHumicsFlux(iHemisphere) = carbonConsHumicsFlux(iHemisphere) + sumArrayOut(iSumPrev+9)
         carbonConsSemiLabileDOCFlux(iHemisphere) = carbonConsSemiLabileDOCFlux(iHemisphere) + sumArrayOut(iSumPrev+10)
      enddo ! iHemisphere

      ! cleanup
      deallocate(sumArray)
      deallocate(sumArrayOut)

      !-------------------------------------------------------------
      ! Carbon conservation error
      !-------------------------------------------------------------

      if (MPAS_stream_mgr_ringing_alarms(domain % streamManager, "conservationCheckOutput", ierr=ierr)) then

         ! convert fluxes to fluxes per m2 of whole earth and per second
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
         call MPAS_pool_get_array(conservationCheckAMPool, "nAccumulate", nAccumulate)
         call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAreaAM", conservationCheckAreaAMPool)
         call MPAS_pool_get_array(conservationCheckAreaAMPool, "earthArea", earthArea)
         fluxScale = 1.0_RKIND / (earthArea * real(nAccumulate))

         ! get initial carbon content
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "initialCarbon", initialCarbon)

         ! get final carbon content
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "finalCarbon", finalCarbon)
         call compute_total_carbon(domain, finalCarbon)

         ! compute the carbon content change
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonChange", carbonChange)
         carbonChange(:) = finalCarbon(:) - initialCarbon(:)
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonChangeFlux", carbonChangeFlux)
         carbonChangeFlux(:) = carbonChange(:) * (fluxScale / dt)

         ! calculate the final net carbon flux to the ice
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "netCarbonFlux", netCarbonFlux)
         netCarbonFlux(:) = &
              - carbonConsOceanCarbonFlux(:)

         ! absolute error
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "absoluteCarbonError", absoluteCarbonError)
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "accumAbsoluteCarbonError", accumAbsoluteCarbonError)
         absoluteCarbonError(:) = netCarbonFlux(:) * dt - carbonChange(:)

         ! rescale fluxes
         carbonConsOceanCarbonFlux(:) = carbonConsOceanCarbonFlux(:) * fluxScale
         carbonConsDiatomFlux(:) = carbonConsDiatomFlux(:) * fluxScale
         carbonConsSmallAlgaeFlux(:) = carbonConsSmallAlgaeFlux(:) * fluxScale
         carbonConsDOC1Flux(:) = carbonConsDOC1Flux(:) * fluxScale
         carbonConsDOC2Flux(:) = carbonConsDOC2Flux(:) * fluxScale
         carbonConsDONCarbonFlux(:) = carbonConsDONCarbonFlux(:) * fluxScale
         carbonConsDICFlux(:) = carbonConsDICFlux(:) * fluxScale
         carbonConsHumicsFlux(:) = carbonConsHumicsFlux(:) * fluxScale
         carbonConsSemiLabileDOCFlux(:) = carbonConsSemiLabileDOCFlux(:) * fluxScale
         carbonConsOceanCarbonFluxCheck(:) = carbonConsOceanCarbonFluxCheck(:) * fluxScale
         netCarbonFlux(:)             = netCarbonFlux(:)             * fluxScale
         carbonConsTOCFlux(:) = carbonConsDiatomFlux(:) + carbonConsSmallAlgaeFlux(:) + &
            carbonConsDOC1Flux(:) + carbonConsDOC2Flux(:) + carbonConsDONCarbonFlux(:) + &
            carbonConsHumicsFlux(:)

         ! compute the final carbon error
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "relativeCarbonError", relativeCarbonError)
         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "accumRelativeCarbonError", accumRelativeCarbonError)
         do iHemisphere = 1, nHemispheres
            if (abs(finalCarbon(iHemisphere)) > 0.0) then
               relativeCarbonError(iHemisphere) = absoluteCarbonError(iHemisphere) / finalCarbon(iHemisphere)
            else
               relativeCarbonError(iHemisphere) = 0.0_RKIND
            endif
         enddo ! iHemisphere

         call MPAS_pool_get_array(conservationCheckCarbonAMPool, "relativeErrorBounds", relativeErrorBounds)
         relativeErrorBounds = relativeErrorBounds + relativeErrorStepBound

         startTime = mpas_get_clock_time(domain % clock, MPAS_START_TIME, ierr=ierr)
         currTime = mpas_get_clock_time(domain % clock, MPAS_NOW, ierr=ierr)

         call MPAS_get_time(currTime, dateTimeString=timeStr, ierr=ierr)
         call MPAS_get_time(startTime, dateTimeString=timeStrStart, ierr=ierr)

         if (trim(timeStr) /= trim(timeStrStart)) then
            accumRelativeCarbonError(:) = accumRelativeCarbonError(:) + relativeCarbonError(:)
            accumAbsoluteCarbonError(:) = accumAbsoluteCarbonError(:) + absoluteCarbonError(:)
            !-------------------------------------------------------------
            ! Output to log file
            !-------------------------------------------------------------

            if (config_AM_conservationCheck_write_to_logfile) then

               call mpas_log_write('-----------------------------------------------------------------------------')
               call mpas_log_write(' Carbon conservation check')
               call mpas_log_write(' ')
               call mpas_log_write(' Initial carbon              (kg) = '//trim(hemisphere_format(initialCarbon)))
               call mpas_log_write(' Final carbon                (kg) = '//trim(hemisphere_format(finalCarbon)))
               call mpas_log_write(' Carbon change               (kg) = '//trim(hemisphere_format(carbonChange)))
               call mpas_log_write(' Carbon change flux      (kg/m2s) = '//trim(hemisphere_format(carbonChangeFlux)))
               call mpas_log_write(' ')
               call mpas_log_write(' Ocean carbon flux       (kg/m2s) = '//trim(hemisphere_format(carbonConsOceanCarbonFlux)))
               call mpas_log_write(' ')
               call mpas_log_write(' Individual ocean carbon fluxes     :')
               call mpas_log_write(' Ocean diatom flux       (kg/m2s) = '//trim(hemisphere_format(carbonConsDiatomFlux)))
               call mpas_log_write(' Ocean small algae flux  (kg/m2s) = '//trim(hemisphere_format(carbonConsSmallAlgaeFlux)))
               call mpas_log_write(' Ocean polysaccharid flux(kg/m2s) = '//trim(hemisphere_format(carbonConsDOC1Flux)))
               call mpas_log_write(' Ocean lipid flux        (kg/m2s) = '//trim(hemisphere_format(carbonConsDOC2Flux)))
               call mpas_log_write(' Ocean DIC flux          (kg/m2s) = '//trim(hemisphere_format(carbonConsDICFlux)))
               call mpas_log_write(' Ocean DON flux          (kg/m2s) = '//trim(hemisphere_format(carbonConsDONCarbonFlux)))
               call mpas_log_write(' Ocean humics flux       (kg/m2s) = '//trim(hemisphere_format(carbonConsHumicsFlux)))
               call mpas_log_write(' Sum of ocean fluxes     (kg/m2s) = '//trim(hemisphere_format(carbonConsOceanCarbonFluxCheck)))
               call mpas_log_write(' ')

               call mpas_log_write(' Total semilabile DOCflux(kg/m2s) = '//trim(hemisphere_format(carbonConsSemiLabileDOCFlux)))
               call mpas_log_write(' Total inorganic C flux  (kg/m2s) = '//trim(hemisphere_format(carbonConsDICFlux)))
               call mpas_log_write(' Total organic C flux    (kg/m2s) = '//trim(hemisphere_format(carbonConsTOCFlux)))
               call mpas_log_write(' ')
               call mpas_log_write(' Net carbon change           (kg) = '//trim(hemisphere_format((netCarbonFlux * dt) / fluxScale)))
               call mpas_log_write(' Net carbon flux         (kg/m2s) = '//trim(hemisphere_format(netCarbonFlux)))
               call mpas_log_write(' ')
               call mpas_log_write(' Absolute carbon error       (kg) = '//trim(hemisphere_format(absoluteCarbonError)))
               call mpas_log_write(' Absolute carbon error/s (kg/m2s) = '//trim(hemisphere_format((absoluteCarbonError * fluxScale) / dt)))
               call mpas_log_write(' Relative carbon error            = '//trim(hemisphere_format(relativeCarbonError)))

               call mpas_log_write(' ')
               call mpas_log_write(' Total run accumulated error ')
               call mpas_log_write(' Accum. absolute carbon error (kg)= '//trim(hemisphere_format(accumAbsoluteCarbonError)))
               call mpas_log_write(' Accum. relative carbon error     = '//trim(hemisphere_format(accumRelativeCarbonError)))
               write(valStr,fmt=' (es15.6)') relativeErrorBounds

               call mpas_log_write(' Accum. relative error bounds     = '//trim(valStr))
            endif
            if ((MAXVAL(abs(accumRelativeCarbonError(:))) > relativeErrorBounds) .or. &
                (MAXVAL(abs(relativeCarbonError(:))) > relativeErrorStepBound) ) then
               if (config_AM_conservationCheck_carbon_failure_abort) then
                  call mpas_log_write('maximum accumulated relative $r carbon error exceeds bounds: $r or relative carbon error $r exceeds step bounds: $r', MPAS_LOG_CRIT, &
                     realArgs=(/MAXVAL(abs(accumRelativeCarbonError(:))), relativeErrorBounds, MAXVAL(abs(relativeCarbonError(:))), relativeErrorStepBound /))
               else
                  call mpas_log_write('WARNING: accumulated relative carbon error exceeds bounds ')
               endif
            endif
         endif
      endif

    end subroutine carbon_conservation

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

      real(kind=RKIND), dimension(:), intent(out) :: &
           totalEnergy

      type(block_type), pointer :: &
           blockPtr

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           tracersPool, &
           conservationCheckAMPool

      integer, pointer :: &
           nCellsSolve, &
           nCategories, &
           nIceLayers, &
           nSnowLayers, &
           nHemispheres

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           iceEnthalpy, &
           snowEnthalpy, &
           iceVolumeCategory, &
           snowVolumeCategory

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      real(kind=RKIND) :: &
           nIceLayersInverse, &
           nSnowLayersInverse

      integer :: &
           iCell, &
           iCategory, &
           iIceLayer, &
           iSnowLayer, &
           iHemisphere

      real(kind=RKIND) :: &
           energyCell

      real(kind=RKIND), dimension(:), allocatable :: &
           energy

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)

      allocate(energy(nHemispheres))
      energy(:) = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "tracers", tracersPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)

         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nCellsSolve', nCellsSolve)
         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nCategories', nCategories)
         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nIceLayers', nIceLayers)
         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nSnowLayers', nSnowLayers)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(tracersPool, "iceEnthalpy", iceEnthalpy, 1)
         call MPAS_pool_get_array(tracersPool, "snowEnthalpy", snowEnthalpy, 1)
         call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)
         call MPAS_pool_get_array(tracersPool, "snowVolumeCategory", snowVolumeCategory, 1)
         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)

         nIceLayersInverse  = 1.0_RKIND / real(nIceLayers,RKIND)
         nSnowLayersInverse = 1.0_RKIND / real(nSnowLayers,RKIND)

         do iCell = 1, nCellsSolve

            energyCell = 0.0;

            do iCategory = 1, nCategories

               do iIceLayer = 1, nIceLayers

                  energyCell = energyCell + &
                       iceEnthalpy(iIceLayer,iCategory,iCell) * &
                       iceVolumeCategory(1,iCategory,iCell) * &
                       nIceLayersInverse * &
                       areaCell(iCell)

               enddo ! iIceLayer

               do iSnowLayer = 1, nSnowLayers

                  energyCell = energyCell + &
                       snowEnthalpy(iSnowLayer,iCategory,iCell) * &
                       snowVolumeCategory(1,iCategory,iCell) * &
                       nSnowLayersInverse * &
                       areaCell(iCell)

               enddo ! iIceLayer

            enddo ! iCategory

            do iHemisphere = 1, nHemispheres
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  energy(iHemisphere) = energy(iHemisphere) + energyCell
               endif
            enddo ! iHemisphere

         enddo ! iCell

         blockPtr => blockPtr % next
      enddo

      ! sum across processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nHemispheres, energy, totalEnergy)

      ! clean up
      deallocate(energy)

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

      real(kind=RKIND), dimension(:), intent(out) :: &
           totalMass

      type(block_type), pointer :: &
           blockPtr

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           tracersPool, &
           tracersAggregatePool, &
           conservationCheckAMPool

      integer, pointer :: &
           nCellsSolve, &
           nCategories, &
           nHemispheres

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

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      logical, pointer :: &
           config_use_topo_meltponds

      integer :: &
           iCell, &
           iCategory, &
           iHemisphere

      real(kind=RKIND) :: &
           massCell

      real(kind=RKIND), dimension(:), allocatable :: &
           mass

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)

      allocate(mass(nHemispheres))
      mass(:) = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_config(blockPtr % configs, "config_use_topo_meltponds", config_use_topo_meltponds)

         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "tracers", tracersPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)

         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nCellsSolve', nCellsSolve)
         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nCategories', nCategories)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)
         call MPAS_pool_get_array(tracersAggregatePool, "snowVolumeCell", snowVolumeCell)
         call MPAS_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
         call MPAS_pool_get_array(tracersPool, "levelIceArea", levelIceArea, 1)
         call MPAS_pool_get_array(tracersPool, "levelIceVolume", levelIceVolume, 1)
         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)

         do iCell = 1, nCellsSolve

            ! ice and snow mass
            massCell = &
                 (iceVolumeCell(iCell)  * rhoi + &
                  snowVolumeCell(iCell) * rhos) * areaCell(iCell)

            do iHemisphere = 1, nHemispheres
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  mass(iHemisphere) = mass(iHemisphere) + massCell
               endif
            enddo ! iHemisphere

         enddo ! iCell

         if (config_use_topo_meltponds) then

            do iCell = 1, nCellsSolve

               massCell = 0.0;

               do iCategory = 1, nCategories

                  ! pond mass
                  massCell = massCell + &
                       iceAreaCategory(1,icategory,iCell) * levelIceArea(1,icategory,iCell) * &
                       levelIceVolume(1,icategory,iCell)  * areaCell(iCell)

               enddo ! iCategory

               do iHemisphere = 1, nHemispheres
                  if (cellInHemisphere(iHemisphere,iCell) == 1) then
                     mass(iHemisphere) = mass(iHemisphere) + massCell
                  endif
               enddo ! iHemisphere

            enddo ! iCell

         endif

         blockPtr => blockPtr % next
      enddo

      ! sum across processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nHemispheres, mass, totalMass)

      ! clean up
      deallocate(mass)

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

      real(kind=RKIND), dimension(:), intent(out) :: &
           totalSalt

      type(block_type), pointer :: &
           blockPtr

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           tracersAggregatePool, &
           conservationCheckAMPool

      integer, pointer :: &
           nCellsSolve, &
           nHemispheres

      real(kind=RKIND), dimension(:), pointer :: &
           areaCell, &
           iceVolumeCell

      real(kind=RKIND), dimension(:,:,:), pointer :: &
           iceVolumeCategory

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      integer :: &
           iCell, &
           iCategory, &
           iHemisphere

      real(kind=RKIND) :: &
           saltCell

      real(kind=RKIND), dimension(:), allocatable :: &
           salt

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)

      allocate(salt(nHemispheres))
      salt(:) = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "tracers_aggregate", tracersAggregatePool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)

         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nCellsSolve', nCellsSolve)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)
         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)

         do iCell = 1, nCellsSolve

            ! ice and snow mass
            saltCell = &
                 iceVolumeCell(iCell) * areaCell(iCell)

            do iHemisphere = 1, nHemispheres
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  salt(iHemisphere) = salt(iHemisphere) + saltCell
               endif
            enddo ! iHemisphere

         enddo ! iCell

         blockPtr => blockPtr % next
      enddo

      salt(:) = salt(:) * rhoi * ice_ref_salinity * 0.001_RKIND

      ! sum across processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nHemispheres, salt, totalSalt)

      ! clean up
      deallocate(salt)

    end subroutine compute_total_salt

!***********************************************************************
!
!  routine compute_total_carbon
!
!> \brief   Compute total carbon of sea ice system
!> \author  Nicole Jeffery
!> \date    27 May 2020
!> \details
!>  Calculate the total carbon of the sea ice system
!
!-----------------------------------------------------------------------

    subroutine compute_total_carbon(domain, totalCarbon)

      use seaice_constants, only: &
           gramsCarbonPerMolCarbon

      type (domain_type), intent(inout) :: &
           domain

      real(kind=RKIND), dimension(:), intent(out) :: &
           totalCarbon

      type(block_type), pointer :: &
           blockPtr

      type(MPAS_pool_type), pointer :: &
           meshPool, &
           biogeochemistryPool, &
           conservationCheckAMPool

      integer, pointer :: &
           nCellsSolve, &
           nHemispheres

      real(kind=RKIND), dimension(:), pointer :: &
           totalCarbonContentCell, &
           areaCell

      integer, dimension(:,:), pointer :: &
           cellInHemisphere

      integer :: &
           iCell, &
           iHemisphere

      real(kind=RKIND) :: &
           carbonCell

      real(kind=RKIND), dimension(:), allocatable :: &
           carbon

      call MPAS_pool_get_dimension(domain % blocklist % dimensions, "nHemispheres", nHemispheres)

      allocate(carbon(nHemispheres))
      carbon(:) = 0.0_RKIND

      blockPtr => domain % blocklist
      do while (associated(blockPtr))

         call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "biogeochemistry", biogeochemistryPool)
         call MPAS_pool_get_subpool(blockPtr % structs, "conservationCheckAM", conservationCheckAMPool)

         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nCellsSolve', nCellsSolve)
         call MPAS_pool_get_dimension(blockPtr % dimensions, 'nHemispheres', nHemispheres)

         call MPAS_pool_get_array(meshPool, "areaCell", areaCell)
         call MPAS_pool_get_array(biogeochemistryPool, "totalCarbonContentCell", totalCarbonContentCell)
         call MPAS_pool_get_array(conservationCheckAMPool, "cellInHemisphere", cellInHemisphere)

         call compute_carbon_cell(blockPtr,totalCarbonContentCell)

         do iCell = 1, nCellsSolve

            ! ice carbon mass (mmols)
            carbonCell = totalCarbonContentCell(iCell) * areaCell(iCell) * &
                         gramsCarbonPerMolCarbon * 1.0e-6_RKIND

            do iHemisphere = 1, nHemispheres
               if (cellInHemisphere(iHemisphere,iCell) == 1) then
                  carbon(iHemisphere) = carbon(iHemisphere) + carbonCell
               endif
            enddo ! iHemisphere

         enddo ! iCell

         blockPtr => blockPtr % next
      enddo

      ! sum across processors
      call MPAS_dmpar_sum_real_array(domain % dminfo, nHemispheres, carbon, totalCarbon)

      ! clean up
      deallocate(carbon)

    end subroutine compute_total_carbon

!***********************************************************************
!
!  compute_carbon_cell
!
!> \brief
!> \author Nicole Jeffery, LANL
!> \date 26 May 2020
!> \details Calculate the total carbon concentration in the sea ice cell
!> by summing the appropriate biogeochemical tracers in units of mmol C
!>
!>      Total carbon = algal nitrogen groups * (C to N ratios) + dissolved carbon groups
!>                   + dissolved inorganic carbon + dissolved organic nitrogen * (C to N ratio)
!>                   + humic material
!
!-----------------------------------------------------------------------

  subroutine compute_carbon_cell(blockPtr,totalCarbonContentCell)

    use seaice_constants, only: &
         skeletalLayerThickness

    real(kind=RKIND), dimension(:), intent(out) :: &
         totalCarbonContentCell

    type(block_type), intent(in) :: &
         blockPtr

    logical, pointer :: &
         config_use_skeletal_biochemistry, &
         config_use_vertical_biochemistry, &
         config_use_vertical_tracers, &
         config_use_carbon, &
         config_use_DON, &
         config_use_humics

    integer, pointer :: &
         nBioLayersP1, &
         nBioLayers, &
         nAlgae, &
         nDOC, &
         nDIC, &
         nDON

    type(MPAS_pool_type), pointer :: &
         mesh, &
         biogeochemistry, &
         tracers_aggregate

    real(kind=RKIND), dimension(:), pointer :: &
         brineFractionCell, &
         iceVolumeCell

    real(kind=RKIND), dimension(:,:), pointer :: &
         skeletalAlgaeConcCell, &
         skeletalDOCConcCell, &
         skeletalDICConcCell, &
         skeletalDONConcCell, &
         skeletalHumicsConcCell, &
         verticalAlgaeConcCell, &
         verticalDOCConcCell, &
         verticalDICConcCell, &
         verticalDONConcCell, &
         verticalHumicsConcCell

    real(kind=RKIND), pointer :: &
         config_ratio_C_to_N_diatoms, &
         config_ratio_C_to_N_small_plankton, &
         config_ratio_C_to_N_phaeocystis, &
         config_ratio_C_to_N_proteins

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    real(kind=RKIND), dimension(:), allocatable :: &
         ratio_C_to_N, &
         verticalGridSpace

    real(kind=RKIND) :: &
         totalOceanCarbonCell

    integer :: &
         iBioTracers, &
         iBioCount, &
         iLayers, &
         iCell

    call MPAS_pool_get_config(blockPtr % configs, "config_use_skeletal_biochemistry", config_use_skeletal_biochemistry)
    call MPAS_pool_get_config(blockPtr % configs, "config_use_vertical_biochemistry", config_use_vertical_biochemistry)
    call MPAS_pool_get_config(blockPtr % configs, "config_use_vertical_tracers", config_use_vertical_tracers)
    call MPAS_pool_get_config(blockPtr % configs, "config_use_carbon", config_use_carbon)
    call MPAS_pool_get_config(blockPtr % configs, "config_use_DON", config_use_DON)
    call MPAS_pool_get_config(blockPtr % configs, "config_use_humics",config_use_humics)
    call MPAS_pool_get_config(blockPtr % configs, "config_ratio_C_to_N_diatoms", config_ratio_C_to_N_diatoms)
    call MPAS_pool_get_config(blockPtr % configs, "config_ratio_C_to_N_small_plankton", config_ratio_C_to_N_small_plankton)
    call MPAS_pool_get_config(blockPtr % configs, "config_ratio_C_to_N_phaeocystis", config_ratio_C_to_N_phaeocystis)
    call MPAS_pool_get_config(blockPtr % configs, "config_ratio_C_to_N_proteins", config_ratio_C_to_N_proteins)

    call MPAS_pool_get_dimension(blockPtr % dimensions, "nBioLayers", nBioLayers)
    call MPAS_pool_get_dimension(blockPtr % dimensions, "nBioLayersP1", nBioLayersP1)
    call MPAS_pool_get_dimension(blockPtr % dimensions, "nAlgae", nAlgae)
    call MPAS_pool_get_dimension(blockPtr % dimensions, "nDOC", nDOC)
    call MPAS_pool_get_dimension(blockPtr % dimensions, "nDIC", nDIC)
    call MPAS_pool_get_dimension(blockPtr % dimensions, "nDON", nDON)

    call MPAS_pool_get_subpool(blockPtr % structs, "tracers_aggregate", tracers_aggregate)
    call MPAS_pool_get_subpool(blockPtr % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(blockPtr % structs, "biogeochemistry", biogeochemistry)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)

    call MPAS_pool_get_array(tracers_aggregate, "skeletalAlgaeConcCell", skeletalAlgaeConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDOCConcCell", skeletalDOCConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDICConcCell", skeletalDICConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalDONConcCell", skeletalDONConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "skeletalHumicsConcCell", skeletalHumicsConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalAlgaeConcCell", verticalAlgaeConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDOCConcCell", verticalDOCConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDICConcCell", verticalDICConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalDONConcCell", verticalDONConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "verticalHumicsConcCell", verticalHumicsConcCell)
    call MPAS_pool_get_array(tracers_aggregate, "brineFractionCell", brineFractionCell)
    call MPAS_pool_get_array(tracers_aggregate, "iceVolumeCell", iceVolumeCell)

    allocate(ratio_C_to_N(3))
    allocate(verticalGridSpace(nBioLayersP1))

    ratio_C_to_N(1) = config_ratio_C_to_N_diatoms
    ratio_C_to_N(2) = config_ratio_C_to_N_small_plankton
    ratio_C_to_N(3) = config_ratio_C_to_N_phaeocystis

    verticalGridSpace(:) = 1.0_RKIND/real(nBioLayers,kind=RKIND)
    verticalGridSpace(1) = verticalGridSpace(1)/2.0_RKIND
    verticalGridSpace(nBioLayersP1) = verticalGridSpace(1)

    totalCarbonContentCell(:) = 0.0_RKIND

    if (config_use_skeletal_biochemistry) then
       do iCell = 1, nCellsSolve

          ! algal nitrogen
          do iBioTracers = 1, nAlgae
             totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) +  skeletalAlgaeConcCell(iBioTracers,iCell)* &
                   skeletalLayerThickness * ratio_C_to_N(iBioTracers)
          enddo

          if (config_use_carbon) then

             ! DOC
             do iBioTracers = 1, nDOC
                totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) +  skeletalDOCConcCell(iBioTracers,iCell)* &
                   skeletalLayerThickness
             enddo

             ! DIC
             do iBioTracers = 1, nDIC
                totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) +  skeletalDICConcCell(iBioTracers,iCell)* &
                   skeletalLayerThickness
             enddo
          endif

          ! DON
          if (config_use_DON) then
             do iBioTracers = 1, nDON
                totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) +  skeletalDONConcCell(iBioTracers,iCell)* &
                   config_ratio_C_to_N_proteins * skeletalLayerThickness
             enddo
          endif

          ! humic material
          if (config_use_humics) &
             totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) +  skeletalHumicsConcCell(1,iCell)* &
                skeletalLayerThickness

       enddo ! iCell

    elseif (config_use_vertical_tracers) then

       do iCell = 1, nCellsSolve

          if (config_use_vertical_biochemistry) then
             iBioCount = 0

             ! algal nitrogen
             do iBioTracers = 1, nAlgae

                do iLayers = 1,nBioLayersP1
                   iBiocount = iBiocount + 1
                   totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) + &
                      verticalAlgaeConcCell(iBioCount,iCell) * ratio_C_to_N(iBioTracers) * &
                      verticalGridSpace(iLayers) * iceVolumeCell(iCell) * brineFractionCell(iCell)
                enddo
                iBiocount = iBioCount + 2
             enddo
          endif

          if (config_use_carbon) then
             iBioCount = 0

             ! DOC
             do iBioTracers = 1, nDOC

                do iLayers = 1,nBioLayersP1
                   iBioCount = iBioCount + 1
                   totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) + &
                      verticalDOCConcCell(iBioCount,iCell) * verticalGridSpace(iLayers) * &
                      iceVolumeCell(iCell) * brineFractionCell(iCell)
                enddo
                iBiocount = iBioCount + 2
             enddo
             iBioCount = 0

             ! DIC
             do iBioTracers = 1, nDIC

                do iLayers = 1,nBioLayersP1
                   iBioCount = iBioCount + 1
                   totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) + &
                      verticalDICConcCell(iBioCount,iCell) * verticalGridSpace(iLayers) * &
                      iceVolumeCell(iCell) * brineFractionCell(iCell)
                enddo
                iBiocount = iBioCount + 2
             enddo
          endif

          if (config_use_DON) then
             iBioCount = 0

             ! dissolve organic nitrogen
             do iBioTracers = 1, nDON

                do iLayers = 1,nBioLayersP1
                   iBiocount = iBiocount + 1
                   totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) + &
                      verticalDONConcCell(iBioCount,iCell) * config_ratio_C_to_N_proteins * &
                      verticalGridSpace(iLayers) * iceVolumeCell(iCell) * brineFractionCell(iCell)
                enddo
                iBiocount = iBioCount + 2
             enddo
          endif

          ! humic material
          if (config_use_humics) then
             do iLayers = 1, nBioLayersP1
                totalCarbonContentCell(iCell) = totalCarbonContentCell(iCell) + &
                   verticalHumicsConcCell(iLayers,iCell) * verticalGridSpace(iLayers) * &
                   iceVolumeCell(iCell) * brineFractionCell(iCell)
             enddo
           endif


        enddo ! iCell
     endif

     deallocate(ratio_C_to_N)
     deallocate(verticalGridSpace)

  end subroutine compute_carbon_cell

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
           conservationCheckAMPool, &
           conservationCheckAreaAMPool, &
           conservationCheckEnergyAMPool, &
           conservationCheckMassAMPool, &
           conservationCheckSaltAMPool, &
           conservationCheckCarbonAMPool

      integer, pointer :: &
           nAccumulate

      real(kind=RKIND), dimension(:), pointer :: &
           accumulatedSeaIceArea

      real(kind=RKIND), dimension(:), pointer :: &
           energyConsSurfaceHeatFlux, &
           energyConsAbsorbedShortwaveFlux, &
           energyConsOceanShortwaveFlux, &
           energyConsSensibleHeatFlux, &
           energyConsLongwaveUp, &
           energyConsLongwaveDown, &
           energyConsOceanHeatFlux, &
           energyConsFreezingPotential, &
           energyConsSnowfallHeat, &
           energyConsLatentHeat

      real(kind=RKIND), dimension(:), pointer :: &
           massConsRainfallRate, &
           massConsSnowfallRate, &
           massConsEvaporation, &
           massConsFreshWater, &
           massConsFrazilWater

      real(kind=RKIND), dimension(:), pointer :: &
           saltConsOceanSaltFlux, &
           saltConsFrazilSaltFlux

      real(kind=RKIND), dimension(:), pointer :: &
           carbonConsOceanCarbonFlux, &
           carbonConsOceanCarbonFluxCheck, &
           carbonConsDiatomFlux, &
           carbonConsSmallAlgaeFlux, &
           carbonConsDOC1Flux, &
           carbonConsDOC2Flux, &
           carbonConsDICFlux, &
           carbonConsDONCarbonFlux, &
           carbonConsHumicsFlux, &
           carbonConsSemiLabileDOCFlux

      ! number of accumulations
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAM", conservationCheckAMPool)
      call MPAS_pool_get_array(conservationCheckAMPool, "nAccumulate", nAccumulate)
      nAccumulate = 0

      ! area
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckAreaAM", conservationCheckAreaAMPool)
      call MPAS_pool_get_array(conservationCheckAreaAMPool, "accumulatedSeaIceArea", accumulatedSeaIceArea)
      accumulatedSeaIceArea(:) = 0.0_RKIND

      ! heat
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckEnergyAM", conservationCheckEnergyAMPool)

      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsSurfaceHeatFlux",       energyConsSurfaceHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsAbsorbedShortwaveFlux", energyConsAbsorbedShortwaveFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsOceanShortwaveFlux",    energyConsOceanShortwaveFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsSensibleHeatFlux",      energyConsSensibleHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsLongwaveUp",            energyConsLongwaveUp)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsLongwaveDown",          energyConsLongwaveDown)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsOceanHeatFlux",         energyConsOceanHeatFlux)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsFreezingPotential",     energyConsFreezingPotential)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsSnowfallHeat",          energyConsSnowfallHeat)
      call MPAS_pool_get_array(conservationCheckEnergyAMPool, "energyConsLatentHeat",            energyConsLatentHeat)

      energyConsSurfaceHeatFlux(:)       = 0.0_RKIND
      energyConsAbsorbedShortwaveFlux(:) = 0.0_RKIND
      energyConsOceanShortwaveFlux(:)    = 0.0_RKIND
      energyConsSensibleHeatFlux(:)      = 0.0_RKIND
      energyConsLongwaveUp(:)            = 0.0_RKIND
      energyConsLongwaveDown(:)          = 0.0_RKIND
      energyConsOceanHeatFlux(:)         = 0.0_RKIND
      energyConsFreezingPotential(:)     = 0.0_RKIND
      energyConsSnowfallHeat(:)          = 0.0_RKIND
      energyConsLatentHeat(:)            = 0.0_RKIND

      ! mass
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckMassAM", conservationCheckMassAMPool)

      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsRainfallRate", massConsRainfallRate)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsSnowfallRate", massConsSnowfallRate)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsEvaporation",  massConsEvaporation)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsFreshWater",   massConsFreshWater)
      call MPAS_pool_get_array(conservationCheckMassAMPool, "massConsFrazilWater",  massConsFrazilWater)

      massConsRainfallRate(:) = 0.0_RKIND
      massConsSnowfallRate(:) = 0.0_RKIND
      massConsEvaporation(:)  = 0.0_RKIND
      massConsFreshWater(:)   = 0.0_RKIND
      massConsFrazilWater(:)  = 0.0_RKIND

      ! salt
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckSaltAM", conservationCheckSaltAMPool)

      call MPAS_pool_get_array(conservationCheckSaltAMPool, "saltConsOceanSaltFlux",  saltConsOceanSaltFlux)
      call MPAS_pool_get_array(conservationCheckSaltAMPool, "saltConsFrazilSaltFlux", saltConsFrazilSaltFlux)

      saltConsOceanSaltFlux(:)  = 0.0_RKIND
      saltConsFrazilSaltFlux(:) = 0.0_RKIND

      ! carbon
      call MPAS_pool_get_subpool(domain % blocklist % structs, "conservationCheckCarbonAM", conservationCheckCarbonAMPool)

      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsOceanCarbonFlux", carbonConsOceanCarbonFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsOceanCarbonFluxCheck", carbonConsOceanCarbonFluxCheck)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDiatomFlux", carbonConsDiatomFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsSmallAlgaeFlux", carbonConsSmallAlgaeFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDOC1Flux", carbonConsDOC1Flux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDOC2Flux", carbonConsDOC2Flux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDICFlux", carbonConsDICFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsDONCarbonFlux", carbonConsDONCarbonFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsHumicsFlux", carbonConsHumicsFlux)
      call MPAS_pool_get_array(conservationCheckCarbonAMPool, "carbonConsSemiLabileDOCFlux", carbonConsSemiLabileDOCFlux)

      carbonConsOceanCarbonFlux(:)  = 0.0_RKIND
      carbonConsOceanCarbonFluxCheck(:)  = 0.0_RKIND
      carbonConsDiatomFlux(:)  = 0.0_RKIND
      carbonConsSmallAlgaeFlux(:)  = 0.0_RKIND
      carbonConsDOC1Flux(:)  = 0.0_RKIND
      carbonConsDOC2Flux(:)  = 0.0_RKIND
      carbonConsDICFlux(:)  = 0.0_RKIND
      carbonConsDONCarbonFlux(:)  = 0.0_RKIND
      carbonConsHumicsFlux(:)  = 0.0_RKIND
      carbonConsSemiLabileDOCFlux(:)  = 0.0_RKIND

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

!***********************************************************************
!
!  routine hemisphere_format
!
!> \brief   Format output string
!> \author  Adrian K. Turner
!> \date    13th April 2021
!> \details
!>  Formats output string for hemispheric values
!
!-----------------------------------------------------------------------

   function hemisphere_format(arrayIn) result(outStr)

     real(kind=RKIND), dimension(:), intent(in) :: &
          arrayIn

     character(len=strKIND) :: &
          outStr

     character(len=16) :: &
          valStr

     integer :: &
          i

     outStr = ""

     do i = 1, size(arrayIn)

        write(valStr,fmt=' (es15.6)') arrayIn(i)
        outStr = trim(outStr)//valStr

     enddo ! i

   end function hemisphere_format

!-----------------------------------------------------------------------

end module seaice_conservation_check

! vim: foldmethod=marker
