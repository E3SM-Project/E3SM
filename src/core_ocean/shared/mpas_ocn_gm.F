! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module ocn_gm

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_threading

   use ocn_constants
   use ocn_config

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

   public :: ocn_GM_compute_Bolus_velocity, &
             ocn_GM_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   private :: tridiagonal_solve

   ! Config options
   real(kind=RKIND), parameter :: epsGM = 1.0e-12_RKIND

   integer :: visbeckTopIndex, visbeckBottomIndex
   logical :: local_config_GM_compute_visbeck
   logical :: local_config_GM_lat_variable_c2
   logical :: local_config_GM_kappa_lat_depth_variable
   real(kind=RKIND), parameter :: local_config_GM_min_phase_speed = 0.1_RKIND
   real(kind=RKIND) :: slopeTaperFactor, sfcTaperFactor, rediGMinitValue
!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_GM_compute_Bolus_velocity
!
!> \brief   Computes GM Bolus velocity
!> \details
!>  This routine is the main driver for the Gent-McWilliams (GM) parameterization.
!>  It computes horizontal and vertical density gradients, the slope
!>  of isopycnal surfaces, and solves a boundary value problem in each column
!>  for the stream function, which is used to compute the Bolus velocity.
!
!-----------------------------------------------------------------------

   subroutine ocn_GM_compute_Bolus_velocity(statePool, diagnosticsPool, &
                                            meshPool, scratchPool, timeLevelIn)
      !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type(mpas_pool_type), intent(in) :: &
         statePool, &! pool containing state variables
         meshPool              ! pool containing mesh quantities

      integer, intent(in), optional :: &
         timeLevelIn          ! time level for state variables

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(mpas_pool_type), intent(inout) :: &
         diagnosticsPool, &! pool containing some diagnostics
         scratchPool           ! pool containing some scratch space

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:, :), pointer :: density, displacedDensity, zMid, normalGMBolusVelocity, &
                                                    layerThicknessEdge, gradDensityEdge, &
                                                    k33, gmStreamFuncTopOfEdge, BruntVaisalaFreqTop, gmStreamFuncTopOfCell, &
                                                    layerThickness, inSituThermalExpansionCoeff, &
                                                    inSituSalineContractionCoeff, RediKappaScaling, gmKappaScaling, RediKappaSfcTaper, &
                                                    zTop, RiTopOfCell

      real(kind=RKIND), dimension(:), pointer   :: c_visbeck, gmBolusKappa, cGMphaseSpeed, &
                ssh, eddyLength, eddyTime, betaEdge, fEdge, gmResolutionTaper, RediKappa
      real(kind=RKIND), dimension(:, :, :), pointer :: slopeTriadUp, slopeTriadDown
      real(kind=RKIND), dimension(:), pointer   :: areaCell, dcEdge, dvEdge, tridiagA, tridiagB, tridiagC, rightHandSide
      integer, dimension(:), pointer   :: maxLevelEdgeTop, maxLevelCell, nEdgesOnCell, indMLD
      integer, dimension(:, :), pointer :: cellsOnEdge, edgesOnCell
      integer                          :: i, k, iEdge, cell1, cell2, iCell, N, iter, iCellSelf, maxLocation
      real(kind=RKIND)                 :: h1, h2, areaEdge, c, BruntVaisalaFreqTopEdge, rtmp
      real(kind=RKIND)                 :: kappaGMEdge, sumN2, countN2, maxN, kappaSum, ltSum
      real(kind=RKIND)                 :: sumRi, RiTopOfEdge, zEdge, zMLD, sshEdge, sfcTaper
      real(kind=RKIND) :: gradDensityTopOfEdge, dcEdgeInv, drhoDx, drhoDT, drhoDS, dTdx, dSdx, BVFcent
      real(kind=RKIND) :: slopeTaperUp, slopeTaperDown, sfcTaperUp, sfcTaperDown, invAreaCell
      real(kind=RKIND), dimension(:), allocatable :: dzTop, dTdzTop, dSdzTop, k33Norm
      ! Dimensions
      integer :: nsmooth, nCells, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      type(field2DReal), pointer :: gradDensityEdgeField

      type(mpas_pool_type), pointer :: tracersPool
      real(kind=RKIND), dimension(:, :, :), pointer :: activeTracers
      integer, pointer :: indexTemperature, indexSalinity
      integer :: timeLevel

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_timer_start('gm bolus velocity')

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperature)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinity)

      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'inSituThermalExpansionCoeff', inSituThermalExpansionCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'inSituSalineContractionCoeff', inSituSalineContractionCoeff)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'indMLD', indMLD)

      call mpas_pool_get_array(diagnosticsPool, 'cGMphaseSpeed', cGMphaseSpeed)
      call mpas_pool_get_array(diagnosticsPool, 'RediKappaScaling', RediKappaScaling)
     call mpas_pool_get_array(diagnosticsPool, 'gmKappaScaling', gmKappaScaling)
      call mpas_pool_get_array(diagnosticsPool, 'RediKappaSfcTaper', RediKappaSfcTaper)
      call mpas_pool_get_array(diagnosticsPool, 'normalGMBolusVelocity', normalGMBolusVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'slopeTriadUp', slopeTriadUp)
      call mpas_pool_get_array(diagnosticsPool, 'slopeTriadDown', slopeTriadDown)
      call mpas_pool_get_array(diagnosticsPool, 'k33', k33)
      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop', BruntVaisalaFreqTop)
      call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfEdge', gmStreamFuncTopOfEdge)
      call mpas_pool_get_array(diagnosticsPool, 'gmStreamFuncTopOfCell', gmStreamFuncTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'gmBolusKappa', gmBolusKappa)
      call mpas_pool_get_array(diagnosticsPool, 'gmResolutionTaper', gmResolutionTaper)
      call mpas_pool_get_array(diagnosticsPool, 'RediKappa', RediKappa)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

      call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_field(scratchPool, 'gradDensityEdge', gradDensityEdgeField)
      call mpas_allocate_scratch_field(gradDensityEdgeField, .True., .false.)

      gradDensityEdge => gradDensityEdgeField%array

      allocate (rightHandSide(nVertLevels))
      allocate (tridiagA(nVertLevels))
      allocate (tridiagB(nVertLevels))
      allocate (tridiagC(nVertLevels))

      nCells = nCellsArray(size(nCellsArray))
      nEdges = nEdgesArray(size(nEdgesArray))

      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            gradDensityEdge(k, iEdge) = 0.0_RKIND
            normalGMBolusVelocity(k, iEdge) = 0.0_RKIND
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            slopeTriadUp(k, :, iEdge) = 0.0_RKIND
            slopeTriadDown(k, :, iEdge) = 0.0_RKIND
         end do
      end do
      !$omp end do

      !$omp do schedule(runtime) private(k)
      do iCell = 1, nCells + 1
         do k = 1, nVertLevels
            k33(k, iCell) = 0.0_RKIND
         end do
      end do
      !$omp end do
      allocate (dzTop(nVertLevels + 1))
      allocate (dTdzTop(nVertLevels + 1))
      allocate (dSdzTop(nVertLevels + 1))
      allocate (k33Norm(nVertLevels + 1))

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         RediKappaScaling(:, iCell) = 1.0_RKIND
         gmKappaScaling(:, iCell) = 1.0_RKIND
         RediKappaSfcTaper(:, iCell) = 1.0_RKIND
      end do
      !$omp end do

      if ( config_Redi_use_N2_based_taper .or. local_config_GM_kappa_lat_depth_variable ) then

         !$omp do schedule(runtime) private(maxLocation, cell1, cell2, k, BruntVaisalaFreqTopEdge, maxN)
         do iCell = 1, nCells
            k = min(maxLevelCell(iCell) - 1, max(1, indMLD(iCell)))
            k = 1
            maxN = max(BruntVaisalaFreqTop(k, iCell), 0.0_RKIND)
            do while (BruntVaisalaFreqTop(k + 1, iCell) > maxN .and. k < maxLevelCell(iCell) - 1)
               k = k + 1
               maxN = max(BruntVaisalaFreqTop(k, iCell), 0.0_RKIND)
            end do

            maxLocation = k
            do k = maxLocation, maxLevelCell(iCell)
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTop(k, iCell), 0.0_RKIND)

               RediKappaScaling(k, iCell) = min(max(config_GM_min_stratification_ratio, &
                                               BruntVaisalaFreqTopEdge/(maxN + 1.0E-10_RKIND)), &
                                           1.0_RKIND)
            end do
         end do
         !$omp end do
      end if

      if (local_config_GM_kappa_lat_depth_variable) then
        !$omp do schedule(runtime)
        do iCell=1,nCells
          gmKappaScaling(:, iCell) = RediKappaScaling(:,iCell)
        enddo
        !$omp end do
      end if

      if (config_Redi_use_surface_taper) then
         !$omp do schedule (runtime) private(zMLD, k)
         do iCell = 1, nCells
            k = max(1, indMLD(iCell))
            zMLD = ssh(iCell) - zMid(k, iCell)

            do k = 1, indMLD(iCell)
               RediKappaSfcTaper(k, iCell) = abs((ssh(iCell) - zMid(k, iCell))/(zMLD))
            end do
         end do
         !$omp end do
      end if

      nCells = nCellsArray(3)
      !$omp do schedule(runtime) private(i, k, iEdge, cell1, cell2, iCellSelf, &
      !$omp& invAreaCell, dcEdgeInv, areaEdge, drhoDT, drhoDS, dTdx, dSdx, drhoDx, &
      !$omp& sfcTaper, sfcTaperUp, sfcTaperDown, slopeTaperDown, slopeTaperUp)
      do iCell = 1, nCells
         invAreaCell = 1.0_RKIND/areaCell(iCell)
         k33(1:maxLevelCell(iCell) + 1, iCell) = 0.0_RKIND
         k33Norm(1:maxLevelCell(iCell) + 1) = epsGM
         ! prep dz, dTdz and dSdz for this column
         do k = 2, maxLevelCell(iCell)
            dzTop(k) = 0.5_RKIND*(layerThickness(k - 1, iCell) + layerThickness(k, iCell))
            dTdzTop(k) = (activeTracers(indexTemperature, k - 1, iCell) &
                          - activeTracers(indexTemperature, k, iCell)) &
                         /dzTop(k)
            dSdzTop(k) = (activeTracers(indexSalinity, k - 1, iCell) &
                          - activeTracers(indexSalinity, k, iCell)) &
                         /dzTop(k)
         end do
         dzTop(1) = -1e-15_RKIND
         dTdzTop(1) = -1e-15_RKIND
         dSdzTop(1) = -1e-15_RKIND
         dzTop(maxLevelCell(iCell) + 1) = -1e-15_RKIND
         dTdzTop(maxLevelCell(iCell) + 1) = -1e-15_RKIND
         dSdzTop(maxLevelCell(iCell) + 1) = -1e-15_RKIND

         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            if (cell1 == iCell) then
               iCellSelf = 1
            else  ! cell2 == iCell
               iCellSelf = 2
            end if
            dcEdgeInv = 1.0_RKIND/dcEdge(iEdge)
            areaEdge = dcEdge(iEdge)*dvEdge(iEdge)

            do k = 1, maxLevelEdgeTop(iEdge)
               drhoDT = -inSituThermalExpansionCoeff(k, iCell)
               drhoDS = inSituSalineContractionCoeff(k, iCell)
               dTdx = (activeTracers(indexTemperature, k, cell2) &
                       - activeTracers(indexTemperature, k, cell1)) &
                      *dcEdgeInv
               dSdx = (activeTracers(indexSalinity, k, cell2) &
                       - activeTracers(indexSalinity, k, cell1)) &
                      *dcEdgeInv
               drhoDx = drhoDT*dTdx + drhoDS*dSdx

               ! Always compute *Up on the top cell and *Down on the bottom
               ! cell, even though they are never used. This avoids an if
               ! statement or separate computation for top and bottom.
               slopeTriadUp(k, iCellSelf, iEdge) = &
                  -drhoDx/ &
                  (drhoDT*dTdzTop(k) &
                   + drhoDS*dSdzTop(k) + 1E-15_RKIND)
               slopeTriadDown(k, iCellSelf, iEdge) = &
                  -drhoDx/ &
                  (drhoDT*dTdzTop(k + 1) &
                   + drhoDS*dSdzTop(k + 1) + 1E-15_RKIND)

               ! set taper of slope ('F' function from Danabasoglu and McWilliams 95)
               if (abs(slopeTriadDown(k, iCellSelf, iEdge)) > 0.6_RKIND*config_redi_maximum_slope) then
                  slopeTaperDown = 0.0_RKIND
               else if (abs(slopeTriadDown(k, iCellSelf, iEdge)) < 0.2_RKIND*config_redi_maximum_slope) then
                  slopeTaperDown = 1.0_RKIND
               else
                  slopeTaperDown = 0.5_RKIND*(1.0_RKIND - ((2.5_RKIND*abs(slopeTriadDown(k, iCellSelf, iEdge)))/ &
                                   config_Redi_maximum_slope - 1.0_RKIND)*(4.0_RKIND - abs(10.0_RKIND* &
                                   abs(slopeTriadDown(k, iCellSelf, iEdge))/config_Redi_maximum_slope - 4.0_RKIND)))
               end if
               if (abs(slopeTriadUp(k, iCellSelf, iEdge)) > 0.6_RKIND*config_redi_maximum_slope) then
                  slopeTaperUp = 0.0_RKIND
               else if (abs(slopeTriadUp(k, iCellSelf, iEdge)) < 0.2_RKIND*config_redi_maximum_slope) then
                  slopeTaperUp = 1.0_RKIND
               else
                  slopeTaperUp = 0.5_RKIND*(1.0_RKIND - ((2.5_RKIND*abs(slopeTriadUp(k, iCellSelf, iEdge)))/ &
                                 config_Redi_maximum_slope - 1.0_RKIND)*(4.0_RKIND - abs(10.0_RKIND* &
                                 abs(slopeTriadUp(k, iCellSelf, iEdge))/config_Redi_maximum_slope - 4.0_RKIND)))
               end if

               slopeTaperUp = 1.0_RKIND + slopeTaperFactor*(slopeTaperUp - 1.0_RKIND)
               slopeTaperDown = 1.0_RKIND + slopeTaperFactor*(slopeTaperDown - 1.0_RKIND)

               sfcTaper = min(RediKappaSfcTaper(k, cell1), RediKappaSfcTaper(k, cell2))
               sfcTaperUp = 1.0_RKIND + sfcTaperFactor*(sfcTaper - 1.0_RKIND)
               sfcTaperDown = 1.0_RKIND + sfcTaperFactor*(sfcTaper - 1.0_RKIND)

               ! Griffies 1998 eqn 34
               if (k > 1) then
                  k33(k, iCell) = k33(k, iCell) + slopeTaperUp*sfcTaperUp* &
                                  areaEdge*dzTop(k)*slopeTriadUp(k, iCellSelf, iEdge)**2
                  k33Norm(k) = k33Norm(k) + areaEdge*dzTop(k)
               end if

               !Taper Redi by tapering the slopes
               k33(k + 1, iCell) = k33(k + 1, iCell) + slopeTaperDown*sfcTaperDown* &
                                   areaEdge*dzTop(k + 1)*slopeTriadDown(k, iCellSelf, iEdge)**2
               k33Norm(k + 1) = k33Norm(k + 1) + areaEdge*dzTop(k + 1)

               slopeTriadUp(k, iCellSelf, iEdge) = &
                  slopeTaperUp*sfcTaperUp*slopeTriadUp(k, iCellSelf, iEdge)
               slopeTriadDown(k, iCellSelf, iEdge) = &
                  slopeTaperDown*sfcTaperDown*slopeTriadDown(k, iCellSelf, iEdge)

            end do ! maxLevelEdgeTop(iEdge)
         end do ! nEdgesOnCell(iCell)

         ! Normalize k33
         do k = 2, maxLevelCell(iCell)
            k33(k, iCell) = k33(k, iCell)/k33Norm(k)*RediKappaSfcTaper(k, iCell)* &
                            RediKappaScaling(k, iCell)
         end do
         k33(1, iCell) = 0.0_RKIND
         k33(maxLevelCell(iCell) + 1, iCell) = 0.0_RKIND
      end do ! iCell
      !$omp end do

      deallocate (dzTop)
      deallocate (dTdzTop)
      deallocate (dSdzTop)
      deallocate (k33Norm)

      ! allow disabling of K33 for testing
      if (config_disable_redi_k33) then
         nCells = nCellsArray(size(nCellsArray))
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            k33(:, iCell) = 0.0_RKIND
         end do
         !$omp end do
      end if

      !--------------------------------------------------------------------
      !
      ! Compute stream function and Bolus velocity for Bolus part of GM
      !
      !--------------------------------------------------------------------

      if (config_use_GM) then
         nEdges = nEdgesArray(3)
         !$omp do schedule(runtime) private(k, cell1, cell2, dcEdgeInv, &
         !$omp& drhoDT, drhoDS, dTdx, dSdx, drhoDx)
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)
            dcEdgeInv = 1.0_RKIND/dcEdge(iEdge)

            do k = 1, maxLevelEdgeTop(iEdge)
               drhoDT = -0.5_RKIND*(inSituThermalExpansionCoeff(k, cell1) + &
                                    inSituThermalExpansionCoeff(k, cell2))
               drhoDS = 0.5_RKIND*(inSituSalineContractionCoeff(k, cell1) + &
                                   inSituSalineContractionCoeff(k, cell2))
               dTdx = (activeTracers(indexTemperature, k, cell2) &
                       - activeTracers(indexTemperature, k, cell1)) &
                      *dcEdgeInv
               dSdx = (activeTracers(indexSalinity, k, cell2) &
                       - activeTracers(indexSalinity, k, cell1)) &
                      *dcEdgeInv
               drhoDx = drhoDT*dTdx + drhoDS*dSdx

               gradDensityEdge(k, iEdge) = drhoDx*rho_sw
            end do
         end do
         !$omp end do

         nEdges = nEdgesArray(3)

         if (local_config_GM_lat_variable_c2) then
            !$omp do schedule(runtime) private(k, cell1, cell2, sumN2, ltSum, countN2, BruntVaisalaFreqTopEdge)
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)
               sumN2 = 0.0
               ltSum = 0.0
               countN2 = 0

               do k = 2, maxLevelEdgeTop(iEdge)

                  BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k, cell1) + BruntVaisalaFreqTop(k, cell2))
                  BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)

                  sumN2 = sumN2 + BruntVaisalaFreqTopEdge*layerThicknessEdge(k, iEdge)
                  ltSum = ltSum + layerThicknessEdge(k, iEdge)
                  countN2 = countN2 + 1

               end do

               if (countN2 > 0) cGMphaseSpeed(iEdge) = max(local_config_GM_min_phase_speed, sqrt(sumN2/ltSum)* &
                            ltSum/(config_GM_baroclinic_mode*3.141592_RKIND))

            end do
            !$omp end do

         else
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               cGMphaseSpeed(iEdge) = config_GM_constant_gravWaveSpeed
            end do
            !$omp end do
         end if

         if (local_config_GM_compute_visbeck) then
           call mpas_pool_get_array(diagnosticsPool, 'betaEdge', betaEdge)
           call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
           call mpas_pool_get_array(diagnosticsPool, 'eddyLength', eddyLength)
           call mpas_pool_get_array(diagnosticsPool, 'eddyTime', eddyTime)
           call mpas_pool_get_array(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)
           call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
           call mpas_pool_get_array(diagnosticsPool, 'c_visbeck', c_visbeck)
           !$omp do schedule(runtime) private(k, cell1, cell2, sumN2, ltSum, sumRi, countN2, &
           !$omp& zEdge, RiTopOfEdge, BruntVaisalaFreqTopEdge)
           do iEdge = 1, nEdges
              k=1
              sumN2 = 0.0_RKIND
              ltSum = 0.0_RKIND
              sumRi = 0.0_RKIND
              countN2 = 0
              cell1 = cellsOnEdge(1, iEdge)
              cell2 = cellsOnEdge(2, iEdge)
              zEdge = 0.5_RKIND*(zTop(k,cell1) + zTop(k,cell2))
              do while(zEdge > -config_GM_visbeck_min_depth .and. k < maxLevelEdgeTop(iEdge))
                k = k+1
                zEdge = 0.5_RKIND*(zTop(k,cell1) + zTop(k,cell2))
              enddo

              k = max(1,k - 1)
              zEdge = 0.5_RKIND*(zTop(k,cell1) + zTop(k,cell2))
              do while(zEdge > -config_GM_visbeck_max_depth .and. k < maxLevelEdgeTop(iEdge))
                BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k, cell1) + BruntVaisalaFreqTop(k, cell2))
                BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
                sumN2 = sumN2 + layerThicknessEdge(k,iEdge)*BruntVaisalaFreqTopEdge
                RiTopOfEdge = 0.5_RKIND*(RiTopOfCell(k,cell1) + RiTopOfCell(k,cell2))
                RiTopOfEdge = max(RiTopOfEdge, 0.0_RKIND)
                sumRi = sumRi + layerThicknessEdge(k,iEdge)*RiTopOfEdge
                ltsum = ltsum + layerThicknessEdge(k,iEdge)
                countN2 = countN2 + 1

                k = k + 1
                zEdge = 0.5_RKIND*(zTop(k,cell1) + zTop(k,cell2))
             end do

             if (countN2 > 0) then
                c_visbeck(iEdge) = sqrt(sumN2*ltsum)
                sumRi = sumRi / (ltsum + 1.0E-11_RKIND)

                eddyLength(iEdge) = min(c_visbeck(iEdge)/(1.0E-15_RKIND + abs(fEdge(iEdge))), &
                                sqrt(c_visbeck(iEdge) / (2.0_RKIND*betaEdge(iEdge))))
                eddyTime(iEdge) = 1.0_RKIND / (max(abs(fEdge(iEdge)), sqrt(2.0_RKIND*c_visbeck(iEdge)*betaEdge(iEdge))) /  &
                                      (1.0E-11_RKIND + sqrt(sumRi)))
                gmBolusKappa(iEdge) = config_GM_visbeck_alpha * eddyLength(iEdge)**2.0_RKIND / (1.0E-15 + eddyTime(iEdge))
                gmBolusKappa(iEdge) = max(config_GM_visbeck_min_kappa, min(gmBolusKappa(iEdge), config_GM_visbeck_max_kappa))
             else !for really shallow columns use min bolus kappa
               gmBolusKappa(iEdge) = config_GM_visbeck_min_kappa
             end if
           end do
           !omp end do
         end if

         if(config_Redi_set_RediKappa_to_GMKappa) then
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               RediKappa(iEdge) = gmResolutionTaper(iEdge)*gmBolusKappa(iEdge)
            end do
            !$omp end do
         endif

         nEdges = nEdgesArray(3)
         !$omp do schedule(runtime) private(cell1, cell2, k,  gradDensityTopOfEdge, &
         !$omp& kappaGMEdge, BruntVaisalaFreqTopEdge, N)
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)

            gmStreamFuncTopOfEdge(:, iEdge) = 0.0_RKIND

            ! Construct the tridiagonal matrix
            if (maxLevelEdgeTop(iEdge) .GE. 3) then
               ! First row
               k = 2
               gradDensityTopOfEdge = (layerThicknessEdge(k - 1, iEdge)*gradDensityEdge(k - 1, iEdge) + &
                                       layerThicknessEdge(k, iEdge)*gradDensityEdge(k, iEdge))/(layerThicknessEdge(k - 1, iEdge) + &
                                                                                                layerThicknessEdge(k, iEdge))
               kappaGMEdge = 0.5_RKIND*(gmKappaScaling(k, cell1) + gmKappaScaling(k, cell2))
               BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k, cell1) + BruntVaisalaFreqTop(k, cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
               tridiagB(k - 1) = -2.0_RKIND*cGMphaseSpeed(iEdge)**2/(layerThicknessEdge(k - 1, iEdge) &
                                                                     *layerThicknessEdge(k, iEdge)) - BruntVaisalaFreqTopEdge
               tridiagC(k - 1) = 2.0_RKIND*cGMphaseSpeed(iEdge)**2/layerThicknessEdge(k, iEdge) &
                                 /(layerThicknessEdge(k - 1, iEdge) + layerThicknessEdge(k, iEdge))
               rightHandSide(k - 1) = gmBolusKappa(iEdge)*kappaGMEdge*gravity/rho_sw &
                                      *gradDensityTopOfEdge

               ! Second to next to the last rows
               do k = 3, maxLevelEdgeTop(iEdge) - 1
                  gradDensityTopOfEdge = (layerThicknessEdge(k - 1, iEdge)*gradDensityEdge(k - 1, iEdge) + &
                                       layerThicknessEdge(k, iEdge)*gradDensityEdge(k, iEdge)) & 
                                       /(layerThicknessEdge(k - 1, iEdge) + &
                                       layerThicknessEdge(k, iEdge))
                  kappaGMEdge = 0.5_RKIND*(gmKappaScaling(k, cell1) + gmKappaScaling(k, cell2))
                  BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k, cell1) + BruntVaisalaFreqTop(k, cell2))
                  BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
                  tridiagA(k - 2) = 2.0_RKIND*cGMphaseSpeed(iEdge)**2/layerThicknessEdge(k - 1, iEdge) &
                                    /(layerThicknessEdge(k - 1, iEdge) + layerThicknessEdge(k, iEdge))
                  tridiagB(k - 1) = -2.0_RKIND*cGMphaseSpeed(iEdge)**2/(layerThicknessEdge(k - 1, iEdge) &
                                                                        *layerThicknessEdge(k, iEdge)) - BruntVaisalaFreqTopEdge
                  tridiagC(k - 1) = 2.0_RKIND*cGMphaseSpeed(iEdge)**2/layerThicknessEdge(k, iEdge) &
                                    /(layerThicknessEdge(k - 1, iEdge) + layerThicknessEdge(k, iEdge))
                  rightHandSide(k - 1) = gmBolusKappa(iEdge)*kappaGMEdge*gravity/rho_sw &
                                         *gradDensityTopOfEdge
               end do

               ! Last row
               k = maxLevelEdgeTop(iEdge)
               gradDensityTopOfEdge = (layerThicknessEdge(k - 1, iEdge)*gradDensityEdge(k - 1, iEdge) + &
                                       layerThicknessEdge(k, iEdge)*gradDensityEdge(k, iEdge))/(layerThicknessEdge(k - 1, iEdge) + &
                                                                                                layerThicknessEdge(k, iEdge))
               kappaGMEdge = 0.5_RKIND*(gmKappaScaling(k, cell1) + gmKappaScaling(k, cell2))
               BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k, cell1) + BruntVaisalaFreqTop(k, cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
               tridiagA(k - 2) = 2.0_RKIND*cGMphaseSpeed(iEdge)**2/layerThicknessEdge(k - 1, iEdge) &
                                 /(layerThicknessEdge(k - 1, iEdge) + layerThicknessEdge(k, iEdge))
               tridiagB(k - 1) = -2.0_RKIND*cGMphaseSpeed(iEdge)**2/(layerThicknessEdge(k - 1, iEdge) &
                                                                     *layerThicknessEdge(k, iEdge)) - BruntVaisalaFreqTopEdge
               rightHandSide(k - 1) = gmBolusKappa(iEdge)*kappaGMEdge*gravity/rho_sw &
                                      *gradDensityTopOfEdge

               ! Total number of rows
               N = maxLevelEdgeTop(iEdge) - 1

               ! Call the tridiagonal solver
               call tridiagonal_solve(tridiagA, tridiagB, tridiagC, rightHandSide, &
                                      gmStreamFuncTopOfEdge(2:maxLevelEdgeTop(iEdge), iEdge), N)
            end if
         end do
         !$omp end do

         nEdges = nEdgesArray(3)
         ! Compute normalGMBolusVelocity from the stream function and apply resolution taper of GM here
         !$omp do schedule(runtime) private(k)
         do iEdge = 1, nEdges
            do k = 1, maxLevelEdgeTop(iEdge)
               normalGMBolusVelocity(k, iEdge) = gmResolutionTaper(iEdge) * (gmStreamFuncTopOfEdge(k, iEdge) &
                      - gmStreamFuncTopOfEdge(k + 1, iEdge)) /layerThicknessEdge(k, iEdge)
            end do
         end do
         !$omp end do

         nCells = nCellsArray(1)

         ! Interpolate gmStreamFuncTopOfEdge to cell centers for visualization
         !$omp do schedule(runtime) private(i, iEdge, areaEdge, k, rtmp)
         do iCell = 1, nCells
            gmStreamFuncTopOfCell(:, iCell) = 0.0_RKIND
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)

               areaEdge = 0.25_RKIND*dcEdge(iEdge)*dvEdge(iEdge)

               do k = 1, maxLevelEdgeTop(iEdge)
                  rtmp = 0.5_RKIND*(gmStreamFuncTopOfEdge(k, iEdge) + gmStreamFuncTopOfEdge(k + 1, iEdge))*areaEdge
                  gmStreamFuncTopOfCell(k, iCell) = gmStreamFuncTopOfCell(k, iCell) + gmResolutionTaper(iEdge) * rtmp
               end do
            end do
         end do
         !$omp end do

         !$omp do schedule(runtime)
         do iCell = 1, nCells
            gmStreamFuncTopOfCell(:, iCell) = gmStreamFuncTopOfCell(:, iCell)/areaCell(iCell)
         end do
         !$omp end do

         deallocate (rightHandSide)
         deallocate (tridiagA)
         deallocate (tridiagB)
         deallocate (tridiagC)

      end if !end config_use_GM

      ! Deallocate scratch variables
      call mpas_deallocate_scratch_field(gradDensityEdgeField, .true.)
      call mpas_timer_stop('gm bolus velocity')

   end subroutine ocn_GM_compute_Bolus_velocity!}}}

!***********************************************************************
!
!  routine tridiagonal_solve
!
!> \brief   Solve the matrix equation Ax=r for x, where A is tridiagonal.
!> \details
!>  Solve the matrix equation Ax=r for x, where A is tridiagonal.
!>  A is an nxn matrix, with:
!>  a sub-diagonal, filled from 1:n-1 (a(1) appears on row 2)
!>  b diagonal, filled from 1:n
!>  c sup-diagonal, filled from 1:n-1  (c(1) apears on row 1)
!
!-----------------------------------------------------------------------
! mrp note:  This subroutine also appears in vmix and should really be put in the framework.
   subroutine tridiagonal_solve(a, b, c, r, x, n) !{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: n
      real(KIND=RKIND), dimension(n), intent(in) :: a, b, c, r

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real(KIND=RKIND), dimension(n), intent(out) :: x

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real(KIND=RKIND), dimension(n) :: bTemp, rTemp
      real(KIND=RKIND) :: m
      integer i

      ! Use work variables for b and r
      bTemp(1) = b(1)
      rTemp(1) = r(1)

      ! First pass: set the coefficients
      do i = 2, n
         m = a(i - 1)/bTemp(i - 1)
         bTemp(i) = b(i) - m*c(i - 1)
         rTemp(i) = r(i) - m*rTemp(i - 1)
      end do

      x(n) = rTemp(n)/bTemp(n)
      ! Second pass: back-substition
      do i = n - 1, 1, -1
         x(i) = (rTemp(i) - c(i)*x(i + 1))/bTemp(i)
      end do

   end subroutine tridiagonal_solve !}}}

!***********************************************************************
!
!  routine ocn_GM_init
!
!> \brief   Initializes ocean momentum horizontal pressure gradient
!> \details
!>  This routine initializes parameters required for the computation of the
!>  horizontal pressure gradient.
!
!-----------------------------------------------------------------------

   subroutine ocn_GM_init(domain, err)!{{{

      type(domain_type), intent(in) :: domain

      integer, intent(out) :: err !< Output: error flag

      real(kind=RKIND) :: sphereRadius

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: meshPool
      type(mpas_pool_type), pointer :: diagnosticsPool
      real(kind=RKIND), dimension(:), pointer :: dcEdge 

      integer :: iEdge
      integer, pointer :: nEdges
      integer, pointer :: nVertLevels
      real(kind=RKIND), pointer :: sphere_radius
      real(kind=RKIND), dimension(:), pointer   :: latEdge, betaEdge, fEdge, &
            gmResolutionTaper, gmBolusKappa
      integer, dimension(:, :), pointer :: cellsOnEdge, edgesOncell

      err = 0

      if ((.not. config_use_GM) .and. (.not. config_use_Redi)) return

      block => domain%blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block%structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(diagnosticsPool, 'gmBolusKappa', gmBolusKappa)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
         call mpas_pool_get_array(diagnosticsPool, 'gmResolutionTaper', gmResolutionTaper)
         if (config_Redi_use_slope_taper) then
            slopeTaperFactor = 1.0_RKIND
         else
            slopeTaperFactor = 0.0_RKIND
         end if

         if (config_Redi_use_surface_taper) then
            if (config_AM_mixedLayerDepths_enable .and. config_AM_mixedLayerDepths_Dthreshold) then
               sfcTaperFactor = 1.0_RKIND
            else
               call mpas_log_write('Redi Surface tapering requires MLD AM enabled with dThresh option selected.', &
                                   MPAS_LOG_CRIT)
               err = 1
               call mpas_dmpar_finalize(domain%dminfo)
            end if
         else
            sfcTaperFactor = 0.0_RKIND
         end if

         RediGMinitValue = 1.0_RKIND
         if (config_GM_closure == 'constant') then
            local_config_GM_lat_variable_c2 = .false.
            local_config_GM_kappa_lat_depth_variable = .false.
            local_config_GM_compute_visbeck = .false.
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               gmBolusKappa(iEdge) = config_GM_kappa
            end do
            !$omp end do
            RediGMinitValue = 1.0_RKIND
         else if (config_GM_closure == 'N2_dependent') then
            local_config_GM_lat_variable_c2 = .true.
            local_config_GM_kappa_lat_depth_variable = .true.
            local_config_GM_compute_visbeck = .false.
            RediGMinitValue = 0.0_RKIND
            ! for N2 dependence, we still assign Kappa as a constant.
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               gmBolusKappa(iEdge) = config_GM_kappa
            end do
            !$omp end do
         else if (config_GM_closure == 'visbeck') then
            local_config_GM_lat_variable_c2 = .true.
            local_config_GM_kappa_lat_depth_variable = .false.
            local_config_GM_compute_visbeck = .true.

            call mpas_pool_get_array(diagnosticsPool, 'betaEdge', betaEdge)
            call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
            call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
            !$omp do schedule(runtime)
            do iEdge=1,nEdges
              betaEdge(iEdge) = 2.0_RKIND*omega*cos(latEdge(iEdge)) / sphere_radius
              gmBolusKappa(iEdge) = config_GM_visbeck_max_kappa
            end do
            !$omp end do
         else
            call mpas_log_write('Invalid choice of config_GM_closure.', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain%dminfo)
         end if

         ! Add resolution taper
         if (config_eddying_resolution_taper == 'none') then
            !$omp do schedule(runtime)
            do iEdge = 1,nEdges
               gmResolutionTaper(iEdge) = 1.0_RKIND
            end do
            !$omp end do
            ! Nothing to do, as we just keep the same assignment as above.
         else if (config_eddying_resolution_taper == 'ramp') then
            !$omp do schedule(runtime) 
            do iEdge = 1, nEdges
               if (dcEdge(iEdge) <= config_eddying_resolution_ramp_min) then
                  gmResolutionTaper(iEdge) = 0.0_RKIND
               else if (dcEdge(iEdge) >= config_eddying_resolution_ramp_max) then
                  gmResolutionTaper(iEdge) = 1.0_RKIND 
               else
                  gmResolutionTaper(iEdge) = 1.0_RKIND &
                                        /(config_eddying_resolution_ramp_max - config_eddying_resolution_ramp_min) &
                                        *(dcEdge(iEdge) - config_eddying_resolution_ramp_min)
               end if
            end do
            !$omp end do
         else
            call mpas_log_write('Invalid choice of config_eddying_resolution_taper.', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain%dminfo)
         end if

         block => block%next
      end do
   end subroutine ocn_GM_init!}}}

!***********************************************************************

end module ocn_gm

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
