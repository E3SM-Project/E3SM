! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
module ocn_gm

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_threading

   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_submesoscale_eddies

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

   public :: ocn_GM_compute_Bolus_velocity,    &
             ocn_GM_init,                      &
             ocn_GM_add_to_transport_vel,      &
             ocn_GM_eddyEnergy_tendency,       &
             ocn_GM_compute_kappa_geometric

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   private :: tridiagonal_solve

   ! Config options
   real(kind=RKIND), parameter :: epsGM = 1.0e-12_RKIND

   ! The following logical variables are used to configure the three
   ! available GM closures (constant, N2_dependent, Visbeck, EdenGreatbatch)
   logical :: local_config_GM_compute_Visbeck
   logical :: local_config_GM_lat_variable_c2
   logical :: local_config_GM_kappa_lat_depth_variable
   logical :: local_config_GM_compute_EdenGreatbatch
   logical :: local_config_GM_diagnostic_geometric
   real(kind=RKIND) :: slopeTaperFactor, sfcTaperFactor, rediGMinitValue
   real(kind=RKIND) :: local_config_GM_constant_bclModeSpeed
   real(kind=RKIND) :: gm_minBclModeSpeed_constant, gm_minBclModeSpeed_compute_on

!***********************************************************************

contains

   function safe_slope(num, den) result(slope)

      real(kind=RKIND), intent(in) :: num, den
      real(kind=RKIND) :: slope

      ! compute isopycnal slopes safely wrt. finite precision
      ! small slope assumption, so don't want |s| > 1

      if (abs(den) > abs(num)) then
         slope = num / den
      else
         slope = sign(1.0_RKIND, den * num)
      end if

   end function

!***********************************************************************
!
!  routine ocn_GM_compute_Bolus_velocity
!
!> \brief   Computes GM Bolus velocity
!> \details
!>  This routine is the main driver for the Gent-McWilliams (GM) parameterization.
!>  It computes GM via the boundary value problem proposed by Ferrari et al (2010)
!>  Ocean Modeling.  It is written as
!>
!>  (c^2 d^2/dz^2 - N^2) psi = grav * kappa_GM / rho_o grad_z rho
!>
!>  Here c is the phase speed of low baroclinic modes, N^2 is the Brunt Vaisala Frequency
!>  the gradient of density is taken at a fixed depth and kappa_GM is the specified
!>  diffusivity.  Herein we have three options to specify the boundary value problem
!>  these options use local flags to specify appropriate values for chosen configuration
!>  The configurations are specified by config_GM_closure, the options are
!>    1. Spatially constant kappa_GM, which also sets the c above to a constant value
!>    2. N2 dependent, this uses three dimensional N2 to modify kappa GM as
!>        kappa_GM = (N/N_ref)^2*kappa_GM(constant
!>       where N_ref is taken as the maximum stratification in the column
!>    3. Visbeck parameterization (Visbeck et al. 1997, JPO)
!>  For options 2 and 3 the phase speed (c) is computed based on column integrated N2, i.e.,
!>   c = (N * H)/(m*pi), where H is the column depth, m is the chosen baroclinic mode
!
!-----------------------------------------------------------------------

   subroutine ocn_GM_compute_Bolus_velocity(statePool, &
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
         scratchPool           ! pool containing some scratch space

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:, :), pointer :: &
                                                    normalVelocity, &
                                                    layerThickness

      real(kind=RKIND), dimension(:), pointer :: totalParamEddyEnergy

      real(kind=RKIND), dimension(:), pointer   :: &
                ssh, fEdge
      real(kind=RKIND), dimension(:), pointer   :: areaCell, dcEdge, dvEdge, tridiagA, tridiagB, tridiagC, rightHandSide
      integer, dimension(:), pointer   :: minLevelEdgeBot, minLevelEdgeTop, maxLevelEdgeTop, minLevelCell, maxLevelCell, nEdgesOnCell
      integer, dimension(:, :), pointer :: cellsOnEdge, edgesOnCell
      integer                          :: i, k, iEdge, cell1, cell2, iCell, N, iCellSelf, maxLocation
      integer                          :: kCell1,kCell2, kEdge
      real(kind=RKIND)                 :: h1, h2, areaEdge, BruntVaisalaFreqTopEdge, rtmp
      real(kind=RKIND)                 :: sum_hN, countN2, BVFref, ltSum
      real(kind=RKIND)                 :: sumRi, RiTopOfEdge, zEdge, zMLD, sfcTaper
      real(kind=RKIND) :: dcEdgeInv, drhoDx, drhoDT, drhoDS, dTdx, dSdx
      real(kind=RKIND) :: slopeTaperUp, slopeTaperDown, sfcTaperUp, sfcTaperDown, invAreaCell
      real(kind=RKIND) :: lt1, lt2, c_min
      real(kind=RKIND) :: sigma, Lr, Length, L_rhines, shearEdgeInv
      real(kind=RKIND) :: eddyLength, eddyTime
      real(kind=RKIND) :: slopeMag, eddyEnergy, kGeometric, drhoDz, hgrad_rho, M2N, N2_ref
      real(kind=RKIND) :: eddyEnergyEdge, M2NInteg, eke,epe
      real(kind=RKIND) :: bvfML, gradBuoy, bvfAv, hML, hAv, indMLDedge
      real(kind=RKIND) :: dzTop1, dTdzTop1, dSdzTop1
      real(kind=RKIND) :: dzTop2, dTdzTop2, dSdzTop2
      real(kind=RKIND), dimension(:), allocatable :: dzTop, dTdzTop, dSdzTop, k33Norm
      real(kind=RKIND), dimension(:), allocatable :: RossbyRadiusTaper
      real(kind=RKIND), dimension(:), allocatable :: structFunction
      real(kind=RKIND) :: c_Visbeck   ! baroclinic wave speed from Visbeck parameterization
      real(kind=RKIND) :: c_mode1, resolution, dc_Lr
      ! Dimensions
      integer :: nCells, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      type(mpas_pool_type), pointer :: tracersPool
      real(kind=RKIND), dimension(:, :, :), pointer :: activeTracers
      integer, pointer :: indexTemperaturePtr, indexSalinityPtr
      integer :: indexTemperature, indexSalinity
      integer :: timeLevel

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_timer_start('gm bolus velocity')

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexTemperaturePtr)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexSalinityPtr)
      indexTemperature = indexTemperaturePtr
      indexSalinity    = indexSalinityPtr

      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeTop', minLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
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

      allocate (rightHandSide(nVertLevels))
      allocate (tridiagA(nVertLevels))
      allocate (tridiagB(nVertLevels))
      allocate (tridiagC(nVertLevels))

      nCells = nCellsArray(size(nCellsArray))
      nEdges = nEdgesArray(size(nEdgesArray))

      !$omp parallel
      !$omp do schedule(runtime) private(k)
      do iEdge = 1, nEdges
         do k = 1, nVertLevels
            normalGMBolusVelocity(k, iEdge) = 0.0_RKIND
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

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         RediKappaSfcTaper(:, iCell) = 1.0_RKIND
      end do
      !$omp end do
      !$omp end parallel

      ! The following code computes scaling for Redi mixing terms and the slope triads
      ! It is only needed when Redi mixing is enabled.
      if ( config_use_Redi ) then

         allocate(dzTop(nVertLevels + 1))
         allocate(dTdzTop(nVertLevels + 1))
         allocate(dSdzTop(nVertLevels + 1))
         allocate(k33Norm(nVertLevels + 1))

         if (config_Redi_use_surface_taper) then
            !$omp parallel
            !$omp do schedule (runtime) private(zMLD, k)
            do iCell = 1, nCells
               k = max(minLevelCell(iCell), indMLD(iCell))
               zMLD = ssh(iCell) - zMid(k, iCell)

               do k = minLevelCell(iCell), indMLD(iCell)
                  RediKappaSfcTaper(k, iCell) = abs((ssh(iCell) - zMid(k, iCell))/(zMLD))
               end do
            end do
            !$omp end do
            !$omp end parallel
         end if

         nCells = nCellsArray(3)
         !$omp parallel
         !$omp do schedule(runtime)  &
         !$omp private(invAreaCell, k33Norm, k, dzTop, dTdzTop, dSdzTop, i, iEdge, cell1, cell2, &
         !$omp         iCellSelf, dcEdgeInv, areaEdge, drhoDT, drhoDS, dTdx, dSdx, drhoDx, &
         !$omp         slopeTaperUp, slopeTaperDown, sfcTaperUp, sfcTaperDown, sfcTaper)
         do iCell = 1, nCells
            invAreaCell = 1.0_RKIND/areaCell(iCell)
            k33(1:maxLevelCell(iCell) + 1, iCell) = 0.0_RKIND
            k33Norm(1:maxLevelCell(iCell) + 1) = epsGM
            ! prep dz, dTdz and dSdz for this column
            do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
               dzTop(k) = 0.5_RKIND*(layerThickness(k - 1, iCell) + layerThickness(k, iCell))
               dTdzTop(k) = (activeTracers(indexTemperature, k - 1, iCell) &
                             - activeTracers(indexTemperature, k, iCell)) &
                            /dzTop(k)
               dSdzTop(k) = (activeTracers(indexSalinity, k - 1, iCell) &
                             - activeTracers(indexSalinity, k, iCell)) &
                            /dzTop(k)
            end do
            dzTop(1:minLevelCell(iCell)) = -1e-15_RKIND
            dTdzTop(1:minLevelCell(iCell)) = -1e-15_RKIND
            dSdzTop(1:minLevelCell(iCell)) = -1e-15_RKIND
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

               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  drhoDT = -thermExpCoeff(k, iCell)
                  drhoDS = salineContractCoeff(k, iCell)
                  dTdx = (activeTracers(indexTemperature, k, cell2) &
                          - activeTracers(indexTemperature, k, cell1)) &
                         *dcEdgeInv
                  dSdx = (activeTracers(indexSalinity, k, cell2) &
                          - activeTracers(indexSalinity, k, cell1)) &
                         *dcEdgeInv
                  drhoDx = drhoDT*dTdx + drhoDS*dSdx


                  sfcTaperUp = drhoDT*dTdzTop(k) + drhoDS*dSdzTop(k)
                  sfcTaperDown = drhoDT*dTdzTop(k+1) + drhoDS*dSdzTop(k+1)

                  ! Always compute *Up on the top cell and *Down on the bottom
                  ! cell, even though they are never used. This avoids an if
                  ! statement or separate computation for top and bottom.
                  slopeTriadUp(k, iCellSelf, iEdge) = &
                     safe_slope( -drhoDx, sfcTaperUp )
                  slopeTriadDown(k, iCellSelf, iEdge) = &
                      safe_slope( -drhoDx, sfcTaperDown )

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
                  if (k < min( indMLD(cell1), indMLD(cell2))) then
                    sfcTaper = 0.0_RKIND
                  else
                    sfcTaper = 1.0_RKIND
                  end if

                  sfcTaperUp = 1.0_RKIND + sfcTaperFactor*(sfcTaper - 1.0_RKIND)
                  sfcTaperDown = 1.0_RKIND + sfcTaperFactor*(sfcTaper - 1.0_RKIND)

!                  ! Griffies 1998 eqn 34
                  if (k > minLevelCell(iCell)) then
                     k33(k, iCell) = k33(k, iCell) +  &
                                     slopeTaperUp*sfcTaperUp*areaEdge*dzTop(k)*slopeTriadUp(k, iCellSelf, iEdge)**2
                     k33Norm(k) = k33Norm(k) + areaEdge*dzTop(k)
                  end if

!                  !Taper Redi by tapering the slopes
                  k33(k + 1, iCell) = k33(k + 1, iCell) +  &
                                      slopeTaperDown*sfcTaperDown*areaEdge*dzTop(k + 1)*slopeTriadDown(k, iCellSelf, iEdge)**2
                  k33Norm(k + 1) = k33Norm(k + 1) + areaEdge*dzTop(k + 1)

                  limiterUp(k,iCellSelf,iEdge) = slopeTaperUp*sfcTaperUp
                  limiterDown(k,iCellSelf,iEdge) = slopeTaperDown*sfcTaperDown

                  slopeTriadUp(k, iCellSelf, iEdge) = &
                     slopeTaperUp*sfcTaperUp*slopeTriadUp(k, iCellSelf, iEdge)
                  slopeTriadDown(k, iCellSelf, iEdge) = &
                     slopeTaperDown*sfcTaperDown*slopeTriadDown(k, iCellSelf, iEdge)

                  !print*, slopeTriadUp(k, iCellSelf, iEdge), slopeTriadDown(k, iCellSelf, iEdge)
                  !print*, k,iEdge,iCellSelf,slopeTriadUp(20,1, iEdge), slopeTriadDown(20,1, iEdge),slopeTriadUp(k, iCellSelf, iEdge), slopeTriadDown(k, iCellSelf, iEdge),slopeTriadDown(k, 1, iEdge)

               end do ! maxLevelEdgeTop(iEdge)
            end do ! nEdgesOnCell(iCell)


!            ! Normalize k33
            do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
               k33(k, iCell) = k33(k, iCell)/k33Norm(k)
            end do
            k33(1:minLevelCell(iCell), iCell) = 0.0_RKIND
            k33(maxLevelCell(iCell) + 1, iCell) = 0.0_RKIND
         end do ! iCell
         !$omp end do
         !$omp end parallel

         deallocate (dzTop)
         deallocate (dTdzTop)
         deallocate (dSdzTop)
         deallocate (k33Norm)

      endif  !config_use_Redi



      ! allow disabling of K33 for testing
      if (config_disable_redi_k33) then
         nCells = nCellsArray(size(nCellsArray))
         !$omp parallel
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            k33(:, iCell) = 0.0_RKIND
         end do
         !$omp end do
         !$omp end parallel
      end if

      !--------------------------------------------------------------------
      !
      ! Compute stream function and Bolus velocity for Bolus part of GM
      !
      !--------------------------------------------------------------------

      if (config_use_GM) then

         nEdges = nEdgesArray(3)

         !--------------------------------------------------------------------------------------------
         ! config_GM_closure = 'geometric'
         !--------------------------------------------------------------------------------------------
         if (local_config_GM_diagnostic_geometric) then
            call ocn_GM_compute_kappa_geometric(statePool, meshPool, timeLevel)
         end if
         !--------------------------------------------------------------------------------------------


         ! For config_GM_closure = 'N2_dependent' use a scaling to taper gmBolusKappa
         ! based on stratification relative to the maximum in the column
         if (local_config_GM_kappa_lat_depth_variable) then
            !$omp parallel
            !$omp do schedule(runtime) private(k, BruntVaisalaFreqTopEdge, BVFref, cell1, cell2)
            do iEdge=1,nEdges
               k=minLevelEdgeBot(iEdge)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k,cell1) + &
                                                    BruntVaisalaFreqTop(k,cell2))
               BVFref = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                 BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k,cell1) + &
                                                      BruntVaisalaFreqTop(k,cell2))
                 BVFref = max(BVFref,max(BruntVaisalaFreqTopEdge, 0.0_RKIND))
               enddo

               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                  BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k,cell1) + &
                                                       BruntVaisalaFreqTop(k,cell2))
                  BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)

                  gmKappaScaling(k, iEdge) = min(max(config_GM_spatially_variable_min_kappa /         &
                                                     config_GM_spatially_variable_max_kappa,          &
                                                     BruntVaisalaFreqTopEdge/(BVFref + 1.0E-10_RKIND)), &
                                                     1.0_RKIND)
                  !if ( k == 10 ) then
                  !   print*, gmKappaScaling(k,iEdge)
                  !endif
                  !if ( k ==   maxLevelEdgeTop(iEdge)-1) then
                  !   print*, gmKappaScaling(k,iEdge)*config_GM_spatially_variable_max_kappa
                  !endif
                  !if ( k == maxLevelEdgeTop(iEdge)-1 ) then
                  !   if ( gmKappaScaling(k, iEdge) * config_GM_spatially_variable_max_kappa > 51.0 ) then
                  !      print*, gmKappaScaling(k, iEdge) * config_GM_spatially_variable_max_kappa
                  !   endif
                  !endif
               end do
           enddo
           !$omp end do
           !$omp end parallel
         end if

         ! For config_GM_closure = 'N2_dependent' or 'Visbeck' or 'EdenGreatbatch' compute a
         ! spatially variable baroclinic phase speed, the mode can be specified by
         ! config_GM_spatially_variable_baroclinic_mode
         if (local_config_GM_lat_variable_c2) then
            !$omp parallel
            !$omp do schedule(runtime) private(k, cell1, cell2, sum_hN, lt1, lt2, ltSum, c_min)
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)
               sum_hN = 0.0_RKIND
               lt1 = 0.0_RKIND
               lt2 = 0.0_RKIND
               ltSum = epsGM
               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)-1

                  lt1 = 0.5_RKIND*(layerThickness(k,cell1) + layerThickness(k-1,cell1))
                  lt2 = 0.5_RKIND*(layerThickness(k,cell2) + layerThickness(k-1,cell2))
                  sum_hN = sum_hN + 0.5_RKIND*(lt1*sqrt(max(BruntVaisalaFreqTop(k,cell1),0.0_RKIND)) + &
                                lt2*sqrt(max(BruntVaisalaFreqTop(k,cell2),0.0_RKIND)))

                  ltSum = ltSum + 0.5*(lt1+lt2)
               end do

               ! Compute the minimum allowed speed of the first baroclinic mode.
               ! If config_GM_minBclModeSpeed_method='constant' use local_config_GM_constant_bclModeSpeed.
               ! If config_GM_minBclModeSpeed_method='computed' use Brunt-Vaisala frequency on this edge.
               ! See initialization of these variables in this modules init routine.
               ! The value 0.5*(lt1+lt2) is meant to be an estimate of the largest thickness
               ! in the column, and c_min is constructed such that we satisfy
               ! the condition c/N < dz.  The condition is an estimate, but
               ! tries to insure that the streamfunction is well resolved even
               ! for the thickest layers.
               c_min = gm_minBclModeSpeed_constant + gm_minBclModeSpeed_compute_on*max(0.01_RKIND,sum_hN/ltSum*(0.5*(lt1+lt2)))

               ! Compute the speed of the first baroclinic mode from the Brunt-Vaisala frequency.
               cGMphaseSpeed(iEdge) = max(c_min,                          &
                                          sum_hN/(config_GM_spatially_variable_baroclinic_mode*pi))

            end do
            !$omp end do
            !$omp end parallel

         else !constant phase speed for config_GM_closure = 'constant'
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               cGMphaseSpeed(iEdge) = local_config_GM_constant_bclModeSpeed
            end do
            !$omp end do
            !$omp end parallel
         end if

         ! When config_GM_closure = 'Visbeck' actually modify the value of gmBolusKappa based on
         ! Visbeck et al (1997), JPO as recast by Cessi (2008), JPO
         if (local_config_GM_compute_Visbeck) then
           call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
           !$omp parallel
           !$omp do schedule(runtime) private(k, cell1, cell2, sum_hN, ltSum, sumRi, countN2, &
           !$omp             zEdge, RiTopOfEdge, BruntVaisalaFreqTopEdge, lt1, lt2, c_Visbeck, &
           !$omp             eddyLength, eddyTime)
           do iEdge = 1, nEdges
              k=minLevelEdgeBot(iEdge)
              sum_hN = 0.0_RKIND
              ltSum = 0.0_RKIND
              sumRi = 0.0_RKIND
              countN2 = 0
              cell1 = cellsOnEdge(1, iEdge)
              cell2 = cellsOnEdge(2, iEdge)

              k = minLevelEdgeBot(iEdge)+1
              zEdge = -layerThickEdgeMean(k-1,iEdge)
              do while(zEdge > -config_GM_Visbeck_max_depth .and. k < maxLevelEdgeTop(iEdge))
                lt1 = 0.5_RKIND*(layerThickness(k,cell1) + layerThickness(k-1,cell1))
                lt2 = 0.5_RKIND*(layerThickness(k,cell2) + layerThickness(k-1,cell2))
                sum_hN = sum_hN + 0.5_RKIND*(lt1*max(BruntVaisalaFreqTop(k,cell1),0.0_RKIND) + &
                                lt2*max(BruntVaisalaFreqTop(k,cell2),0.0_RKIND))
                sumRi = sumRi + 0.5_RKIND*(lt1*max(RiTopOfCell(k,cell1),0.0_RKIND) + &
                                    lt2*max(RiTopOfCell(k,cell2),0.0_RKIND))
                ltsum = ltsum + 0.5_RKIND*(lt1+lt2)
                countN2 = countN2 + 1

                k = k + 1
                zEdge = zEdge - layerThickEdgeMean(k-1,iEdge)
              end do

              if (countN2 > 0) then
                c_Visbeck = sqrt(sum_hN*ltsum)
                sumRi = sumRi / (ltsum + 1.0E-11_RKIND)

                eddyLength = max(dcEdge(iEdge), min(c_Visbeck/(1.0E-15_RKIND + abs(fEdge(iEdge))), &
                                sqrt(c_Visbeck/ (2.0_RKIND*betaEdge(iEdge)))))
                eddyTime = 1.0_RKIND / (1.0E-11_RKIND + max(abs(fEdge(iEdge)), sqrt(2.0_RKIND*c_Visbeck*betaEdge(iEdge))) /  &
                                      (1.0E-11_RKIND + sqrt(sumRi)))
                gmBolusKappa(iEdge) = config_GM_Visbeck_alpha * eddyLength**2.0_RKIND / (1.0E-15 + eddyTime)
                gmBolusKappa(iEdge) = max(config_GM_spatially_variable_min_kappa, min(gmBolusKappa(iEdge),  &
                                          config_GM_spatially_variable_max_kappa))
             else !for really shallow columns use min bolus kappa
                gmBolusKappa(iEdge) = config_GM_spatially_variable_min_kappa
             end if
           end do
           !$omp end do
           !$omp end parallel
         end if

         ! When config_GM_closure = 'EdenGreatbatch' actually modify the value of gmBolusKappa based on
         ! Eden and Greatbach (2008), Ocean Model., but assuming quasi-geostrophy to simplify the EKE budget
         if (local_config_GM_compute_EdenGreatbatch) then
           call mpas_pool_get_array(meshPool, 'fEdge', fEdge)

           !$omp parallel
           !$omp do schedule(runtime) private(cell1, cell2, sigma, RiTopOfEdge, &
           !$omp             L_rhines, Lr, Length, BruntVaisalaFreqTopEdge, shearEdgeInv)
           do iEdge=1,nEdges
              cell1 = cellsOnEdge(1, iEdge)
              cell2 = cellsOnEdge(2, iEdge)
              Lr = min(cGMphaseSpeed(iEdge) / (1.0E-15_RKIND + abs(fEdge(iEdge))), &
                       sqrt(cGMphaseSpeed(iEdge) / (2.0_RKIND*betaEdge(iEdge))))

              do k=minLevelEdgeBot(iEdge)+1,maxLevelEdgeTop(iEdge)
                 BruntVaisalaFreqTopEdge = 0.5_RKIND*(max(BruntVaisalaFreqTop(k,cell1),0.0_RKIND) + &
                                          max(BruntVaisalaFreqTop(k,cell2),0.0_RKIND))
                 shearEdgeInv = (0.5_RKIND*(layerThickEdgeMean(k-1,iEdge) + layerThickEdgeMean(k,iEdge)))**2.0 &
                    / (  (normalVelocity(k-1,iEdge) - normalVelocity(k,iEdge))**2.0 &
                       + (tangentialVelocity(k-1,iEdge) - tangentialVelocity(k,iEdge))**2.0 &
                       + 1.0e-18_RKIND )
                 RiTopOfEdge = BruntVaisalaFreqTopEdge * shearEdgeInv
                 sigma = max(abs(fEdge(iEdge)),sqrt(2.0_RKIND*betaEdge(iEdge)* &
                             cGMphaseSpeed(iEdge))) / sqrt(RiTopOfEdge + config_GM_EG_riMin)
                 L_rhines = sigma / (1.0E-18_RKIND + betaEdge(iEdge))
                 Length = min(config_GM_EG_Rossby_factor*Lr,config_GM_EG_Rhines_factor*L_rhines)
                 !For compatibility with other schemes, we normalize by the 2-D kappa field, which will
                 !be removed below in the computation of the stream function
                 gmKappaScaling(k,iEdge) = max(config_GM_spatially_variable_min_kappa, min(config_GM_EG_kappa_factor* &
                         Length**2.0*sigma, config_GM_spatially_variable_max_kappa)) / gmBolusKappa(iEdge)
              enddo
           enddo
           !$omp end do
           !$omp end parallel
         endif

         nEdges = nEdgesArray(3)
         if (config_GM_horizontal_taper == 'RossbyRadius'.or. &
             config_Redi_horizontal_taper == 'RossbyRadius') then
            allocate (RossbyRadiusTaper(nEdges))
            call mpas_pool_get_array(meshPool, 'fEdge', fEdge)

            !$omp parallel
            !$omp do schedule(runtime) private(k, cell1, cell2, sum_hN, lt1, lt2, c_mode1, resolution, dc_Lr)
            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)
               sum_hN = 0.0

               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)-1

                  lt1 = 0.5_RKIND*(layerThickness(k,cell1) + layerThickness(k-1,cell1))
                  lt2 = 0.5_RKIND*(layerThickness(k,cell2) + layerThickness(k-1,cell2))
                  sum_hN = sum_hN + 0.5_RKIND*(lt1*sqrt(max(BruntVaisalaFreqTop(k,cell1),0.0_RKIND)) + &
                                  lt2*sqrt(max(BruntVaisalaFreqTop(k,cell2),0.0_RKIND)))
               end do

               ! See Chelton (1998)
               ! https://journals.ametsoc.org/view/journals/phoc/28/3/1520-0485_1998_028_0433_gvotfb_2.0.co_2.xml
               ! eqn 2.2a, also Appendix A,
               ! for the first-mode (m=1) baroclinic gravity wave speed in m/s
               ! c_m = 1/pi integral(N dz )

               c_mode1 = max(1.0E-5_RKIND, sum_hN/pi)

               ! See Hallberg (2013) https://doi.org/10.1016/j.ocemod.2013.08.007
               ! just after eqn 4 the defines the deformation
               ! radius as:
               ! L = sqrt( c_g^2/ ( f^2 + 2*beta*c_g) )
               !   = c_g * sqrt(1/( f^2 + 2*beta*c_g) )  (note c_g is guaranteed non-negative)
               ! This comes from a combination of Chelton (1998) eqns 2.3a and
               ! 2.3b, with a smooth transition in between.
               RossbyRadius(iEdge) = c_mode1*sqrt( 1.0_RKIND / &
                  (fEdge(iEdge)*fEdge(iEdge) + 2.0_RKIND*betaEdge(iEdge)*c_mode1 + 1.0E-15_RKIND))

               dc_Lr = dcEdge(iEdge)/RossbyRadius(iEdge)

               ! Ramp the parameterization between the given min and max values
               ! of the ratio of the cell width to the Rossby radius
               if (dc_Lr <= config_GMRedi_Rossby_ramp_min) then
                  RossbyRadiusTaper(iEdge) = 0.0_RKIND
               else if (dc_Lr >= config_GMRedi_Rossby_ramp_max) then
                  RossbyRadiusTaper(iEdge) = 1.0_RKIND
               else
                  ! a linear ramp between the min and max of dc_Lr
                  RossbyRadiusTaper(iEdge) = (dc_Lr - config_GMRedi_Rossby_ramp_min) / &
                     (config_GMRedi_Rossby_ramp_max - config_GMRedi_Rossby_ramp_min)
               end if
            end do
            !$omp end do
            !$omp end parallel

            if (config_GM_horizontal_taper == 'RossbyRadius') then
               !$omp parallel
               !$omp do schedule(runtime)
               do iEdge = 1,nEdges
                  gmHorizontalTaper(iEdge) = RossbyRadiusTaper(iEdge)
               end do
               !$omp end do
               !$omp end parallel
            end if

            if (config_Redi_horizontal_taper == 'RossbyRadius') then
               !$omp parallel
               !$omp do schedule(runtime)
               do iEdge = 1,nEdges
                  RediHorizontalTaper(iEdge) = RossbyRadiusTaper(iEdge)
               end do
               !$omp end do
               !$omp end parallel
            end if

            deallocate (RossbyRadiusTaper)
         end if

         if (config_Redi_closure == 'equalGM') then
            ! Set Redi closure to be the same as the GM closure, i.e. if
            ! config_GM_closure is EdenGreatbatch or Visbeck, you will see those
            ! same closure fields in RediKappa. The horizontal tapering still
            ! comes from the config_Redi_horizontal_taper flag, not GM.
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge=1, nEdges
               RediKappa(iEdge) = gmHorizontalTaper(iEdge) * gmBolusKappa(iEdge)
            end do
            !$omp end do
            !$omp end parallel
         else if (config_Redi_closure == 'constant'.and. &
                  config_Redi_horizontal_taper == 'RossbyRadius') then
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge=1, nEdges
               RediKappa(iEdge) = RediHorizontalTaper(iEdge) * config_Redi_constant_kappa
            end do
            !$omp end do
            !$omp end parallel
         end if

         !!$omp parallel
         !!$omp do schedule(runtime) private(cell1, cell2, k, &
         !!$omp             BruntVaisalaFreqTopEdge, N) &
         !!$omp firstprivate(tridiagB, tridiagC, rightHandSide, tridiagA)
         do iEdge = 1, nEdges
            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)

            gmStreamFuncTopOfEdge(:, iEdge) = 0.0_RKIND

            ! Construct the tridiagonal matrix if there are at least 3 active cells
            if ((maxLevelEdgeTop(iEdge) - minLevelEdgeBot(iEdge)) .GE. 2) then
               ! First row
               k = minLevelEdgeBot(iEdge)+1

               BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k, cell1) + BruntVaisalaFreqTop(k, cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
               tridiagB(k - 1) = -2.0_RKIND*cGMphaseSpeed(iEdge)**2/(layerThickEdgeMean(k - 1, iEdge) &
                                                                     *layerThickEdgeMean(k, iEdge)) - BruntVaisalaFreqTopEdge
               tridiagC(k - 1) = 2.0_RKIND*cGMphaseSpeed(iEdge)**2/layerThickEdgeMean(k, iEdge) &
                                 /(layerThickEdgeMean(k - 1, iEdge) + layerThickEdgeMean(k, iEdge))
               rightHandSide(k - 1) = gmBolusKappa(iEdge)*gmKappaScaling(k,iEdge)*gravity/rho_sw* &
                                            gradDensityEddy(k,iEdge)

               ! Second to next to the last rows
               do k = minLevelEdgeBot(iEdge) + 2, maxLevelEdgeTop(iEdge) - 1
                  BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k, cell1) + BruntVaisalaFreqTop(k, cell2))
                  BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
                  tridiagA(k - 2) = 2.0_RKIND*cGMphaseSpeed(iEdge)**2/layerThickEdgeMean(k - 1, iEdge) &
                                    /(layerThickEdgeMean(k - 1, iEdge) + layerThickEdgeMean(k, iEdge))
                  tridiagB(k - 1) = -2.0_RKIND*cGMphaseSpeed(iEdge)**2/(layerThickEdgeMean(k - 1, iEdge) &
                                                                        *layerThickEdgeMean(k, iEdge)) - BruntVaisalaFreqTopEdge
                  tridiagC(k - 1) = 2.0_RKIND*cGMphaseSpeed(iEdge)**2/layerThickEdgeMean(k, iEdge) &
                                    /(layerThickEdgeMean(k - 1, iEdge) + layerThickEdgeMean(k, iEdge))
                  rightHandSide(k - 1) = gmBolusKappa(iEdge)*gmKappaScaling(k,iEdge)*gravity/rho_sw* &
                                             gradDensityEddy(k,iEdge)
               end do

               ! Last row
               k = maxLevelEdgeTop(iEdge)
               BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k, cell1) + BruntVaisalaFreqTop(k, cell2))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
               tridiagA(k - 2) = 2.0_RKIND*cGMphaseSpeed(iEdge)**2/layerThickEdgeMean(k - 1, iEdge) &
                                 /(layerThickEdgeMean(k - 1, iEdge) + layerThickEdgeMean(k, iEdge))
               tridiagB(k - 1) = -2.0_RKIND*cGMphaseSpeed(iEdge)**2/(layerThickEdgeMean(k - 1, iEdge) &
                                                                     *layerThickEdgeMean(k, iEdge)) - BruntVaisalaFreqTopEdge
               rightHandSide(k - 1) = gmBolusKappa(iEdge)*gmKappaScaling(k,iEdge)*gravity/rho_sw* &
                                            gradDensityEddy(k,iEdge)
               ! Total number of rows
               N = maxLevelEdgeTop(iEdge) - minLevelEdgeBot(iEdge)

               ! Call the tridiagonal solver
               call tridiagonal_solve(tridiagA(minLevelEdgeBot(iEdge):maxLevelEdgeTop(iEdge)), &
                                      tridiagB(minLevelEdgeBot(iEdge):maxLevelEdgeTop(iEdge)), &
                                      tridiagC(minLevelEdgeBot(iEdge):maxLevelEdgeTop(iEdge)), &
                                      rightHandSide(minLevelEdgeBot(iEdge):maxLevelEdgeTop(iEdge)), &
                                      gmStreamFuncTopOfEdge(minLevelEdgeBot(iEdge)+1:maxLevelEdgeTop(iEdge), iEdge), N)
            end if
         end do
         !!$omp end do
         !!$omp end parallel

         ! Compute normalGMBolusVelocity from the stream function and apply resolution taper of GM here
         !$omp parallel
         !$omp do schedule(runtime) private(k)
         do iEdge = 1, nEdges
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               normalGMBolusVelocity(k, iEdge) = gmHorizontalTaper(iEdge) * (gmStreamFuncTopOfEdge(k, iEdge) &
                      - gmStreamFuncTopOfEdge(k + 1, iEdge)) /layerThickEdgeMean(k, iEdge)
            end do
         end do
         !$omp end do
         !$omp end parallel

         deallocate (rightHandSide)
         deallocate (tridiagA)
         deallocate (tridiagB)
         deallocate (tridiagC)

      end if !end config_use_GM

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
!  routine ocn_GM_add_to_transport_velocity
!
!> \brief   Adds GM to transportVel
!> \details
!>  Adds the GM bolus velocity to the normal transport velocity if enabled
!
!-----------------------------------------------------------------------

    subroutine ocn_GM_add_to_transport_vel(normalTransportVelocity, &
                                           nEdges, nVertLevels)!{{{

       real (kind=RKIND), dimension(:,:), intent(inout) :: &
            normalTransportVelocity

       integer, intent(in) :: nEdges, nVertLevels

       integer :: iEdge, k

       if (.not. config_use_gm) return

       !$omp parallel
       !$omp do schedule(runtime) &
       !$omp private(k)
       do iEdge = 1, nEdges
          do k = 1, nVertLevels
             normalTransportVelocity(k,iEdge) = normalTransportVelocity(k,iEdge) + &
                                    normalGMBolusVelocity(k,iEdge)
          end do
       end do
       !$omp end do
       !$omp end parallel

    end subroutine ocn_GM_add_to_transport_vel!}}}

!***********************************************************************
!
!  routine ocn_GM_compute_kappa_geometric
!
!> \brief   Compute tendency of the eddy energy budget equation for
!           the GEOMETRIC parameterization
!> \details
!>  Compute tendency of the eddy energy budget equation for the
!   GEOMETRIC parameterization if 'geometric' GM is enabled.
!
!-----------------------------------------------------------------------

    subroutine ocn_GM_compute_kappa_geometric(statePool, meshPool, &
                                              timeLevelIn)!{{{

       !-----------------------------------------------------------------
       !
       ! input variables
       !
       !-----------------------------------------------------------------

       type(mpas_pool_type), intent(in) :: &
          statePool, &! pool containing state variables
          meshPool    ! pool containing mesh quantities

       integer, intent(in), optional :: &
          timeLevelIn          ! time level for state variables

       !-----------------------------------------------------------------
       !
       ! local variables
       !
       !-----------------------------------------------------------------

       real(kind=RKIND), dimension(:, :), pointer :: layerThickness

       real(kind=RKIND), dimension(:), pointer :: totalParamEddyEnergy

       integer, dimension(:), pointer    :: minLevelEdgeBot, minLevelEdgeTop, &
                                            maxLevelEdgeTop, minLevelCell,    &
                                            maxLevelCell
       integer, dimension(:, :), pointer :: cellsOnEdge
       integer                           :: i, k, iEdge, cell1, cell2, iCell
       integer                           :: kCell1,kCell2, kEdge
       real(kind=RKIND)                  :: BruntVaisalaFreqTopEdge, rtmp
       real(kind=RKIND) :: slopeMag, eddyEnergy, kGeometric, M2N, N2_ref
       real(kind=RKIND) :: eddyEnergyEdge, M2NInteg, eke,epe, BVFref
       real(kind=RKIND) :: gradBuoy, hML, hAv
       real(kind=RKIND), dimension(:), allocatable :: structFunction
       ! Dimensions
       integer :: nCells, nEdges, kMLD
       integer, pointer :: nVertLevels
       integer, dimension(:), pointer :: nCellsArray, nEdgesArray

       integer :: timeLevel

       if (.not. config_use_gm .and. config_GM_closure /= 'geometric') return

       if (present(timeLevelIn)) then
          timeLevel = timeLevelIn
       else
          timeLevel = 1
       end if

       call mpas_pool_get_array(statePool, 'totalParamEddyEnergy', totalParamEddyEnergy, timeLevel)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

       call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
       call mpas_pool_get_array(meshPool, 'minLevelEdgeTop', minLevelEdgeTop)
       call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
       call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

       call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

       !--------------------------------------------------------------------------------------------
       !--------------------------------------------------------------------------------------------
       ! The GEOMETRIC parameterization
       !
       ! K_gm = alpha * int(E)dz * Gamma(z)
       !        / int(Gamma(z)*M^2/N)dz
       !                                   (Eq. (4) in Mak et al., 2018)
       !                                     doi:10.1175/JPO-D-18-0017.1
       !
       ! where,
       !   E = Prognostic eddy energy (Marshall 2012, Mak 2018)
       !   S = slope = drhoDX / drhoDZ
       !   h = layer thickness
       !   M^2 = |grad(b)| ; b=buoyancy
       !   N^2 = Brunt-Vaisala frequency
       !   Gamma(z) = N^2 / N^2_ref (Ferreira et a., 2005)
       !--------------------------------------------------------------------------------------------
       !--------------------------------------------------------------------------------------------

       nEdges = nEdgesArray(3)



       ! Temp arrays to commpute slope (codes from the submesoscale param)

       allocate(structFunction(nVertLevels + 1))

       structFunction(:) = 1.0_RKIND


       do iEdge=1,nEdges
          cell1 = cellsOnEdge(1, iEdge)
          cell2 = cellsOnEdge(2, iEdge)

          ! Build the structure function
          kMLD = min( indMLD(cell1), indMLD(cell2) )
          BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(kMLD,cell1) + &
                                               BruntVaisalaFreqTop(kMLD,cell2))
          ! N_ref at the mixed layer
          BVFref = max(BruntVaisalaFreqTopEdge, 1.0E-14_RKIND)

          do k = minLevelEdgeBot(iEdge),maxLevelEdgeTop(iEdge)
             BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k,cell1) + &
                                                  BruntVaisalaFreqTop(k,cell2))
             BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 1.0E-14_RKIND)

             if (k <= kMLD) then ! Above and at mixed layer
                structFunction(k) = 1.0_RKIND
             else
                structFunction(k) = max(0.0_RKIND,BruntVaisalaFreqTopEdge) &
                                  / BVFref

                !!!!!!!!!!!!!!!!!!!!!!
                !structFunction(k) = max(1.0_RKIND, min(structFunction(k),1.0_RKIND)) ! For vertically constant GEOMETRIC
                !!!!!!!!!!!!!!!!!!!!!!
                structFunction(k) = max(0.01, min(structFunction(k),1.0_RKIND))
                !!!!!!!!!!!!!!!!!!!!!!
             endif
          end do

          k = minLevelEdgeBot(iEdge)
          hML = 0.5_RKIND*layerThickEdgeMean(k,iEdge)

          ! N^2
          BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k,cell1) + &
                                               BruntVaisalaFreqTop(k,cell2))
          BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 1.0E-14_RKIND)
          BVFref = BruntVaisalaFreqTopEdge

          M2NInteg = 0.0_RKIND

          !-------------------------------------------------------------------------
          ! Surface level
          !structFunction(k) = 1.0_RKIND

          ! M^2
          !gradBuoy = abs(gradDensityEddy(k,iEdge) * gravity / rho_sw)

          ! M^2 / N
          !M2N = gradBuoy / sqrt(BruntVaisalaFreqTopEdge)

          ! M^2 / N = slope (=M^2/N^2) * sqrt(BVF) (=N)
          slopeMag =  0.25_RKIND                              &
                   * (  abs(slopeTriadUp(  k,1,iEdge))  &
                      + abs(slopeTriadUp(  k,2,iEdge))  &
                      + abs(slopeTriadDown(k,1,iEdge))  &
                      + abs(slopeTriadDown(k,2,iEdge)) )
          M2N = slopeMag * sqrt(BruntVaisalaFreqTopEdge)

          M2NInteg = structFunction(k) * M2N * hML

          if (minLevelEdgeTop(iEdge) .ge. 1) then
             k = minLevelEdgeTop(iEdge)+1
             hML = 0.5_RKIND*layerThickEdgeMean(k,iEdge)
             BruntVaisalaFreqTopEdge = 0.5_RKIND*(max(BruntVaisalaFreqTop(k,cell1),1.0E-14_RKIND) + &
                                                  max(BruntVaisalaFreqTop(k,cell2),1.0E-14_RKIND))
             BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 1.0e-14_RKIND)
             BVFref = BruntVaisalaFreqTopEdge

             !structFunction(1) = BVFref / BVFref = 1.0
             !structFunction(k) = 1.0_RKIND
             !gradBuoy = gradDensityEddy(minLevelEdgeTop(iEdge)+1,iEdge) * gravity / rho_sw

             ! M^2 / N
             !M2N = gradBuoy / sqrt(BruntVaisalaFreqTopEdge)


             ! M^2 / N = slope (=M^2/N^2) * sqrt(BVF) (=N)
             slopeMag =  0.25_RKIND                              &
                      * (  abs(slopeTriadUp(  k,1,iEdge))  &
                         + abs(slopeTriadUp(  k,2,iEdge))  &
                         + abs(slopeTriadDown(k,1,iEdge))  &
                         + abs(slopeTriadDown(k,2,iEdge)) )
             M2N = slopeMag * sqrt(BruntVaisalaFreqTopEdge)

             M2NInteg = structFunction(k) * M2N * hML
          end if

          !-------------------------------------------------------------------------

          do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
             ! N^2
             BruntVaisalaFreqTopEdge = 0.5_RKIND*(BruntVaisalaFreqTop(k,cell1) + &
                                                  BruntVaisalaFreqTop(k,cell2))
             BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 1.0e-14_RKIND)

             ! Structure function (Gamma(z))
             !structFunction(k) = BruntVaisalaFreqTopEdge / BVFref

             ! Depth integration
             hML = layerThickEdgeMean(k,iEdge)

             ! M^2 = |gradBuoy|
!            gradBuoy = abs(gradDensityEddy(k,iEdge) * gravity / rho_sw)

             ! M^2 / N
!            M2N = gradBuoy / sqrt(BruntVaisalaFreqTopEdge)

             ! M^2 / N = slope^2 (=M^2/N^2) * sqrt(BVF) (=N)
             slopeMag =  0.25_RKIND                              &
                      * (  abs(slopeTriadUp(  k,1,iEdge))  &
                         + abs(slopeTriadUp(  k,2,iEdge))  &
                         + abs(slopeTriadDown(k,1,iEdge))  &
                         + abs(slopeTriadDown(k,2,iEdge)) )
             M2N = slopeMag * sqrt(BruntVaisalaFreqTopEdge)

             M2NInteg = M2NInteg + structFunction(k) * M2N * hML

          end do ! k

          do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
             eddyEnergyEdge = 0.5_RKIND*(totalParamEddyEnergy(cell1) + &
                                         totalParamEddyEnergy(cell2))

             ! K_gm = alpha * (int(E)dz / int(Gamma(z)*M^2/N)dz) * Gamma(z)
             kGeometric = config_GM_geometric_alpha * (eddyEnergyEdge / M2NInteg) &
                                                    * structFunction(k)

             ! Cap kGeometric ( GM_min_kappa < kGeometric < GM_max_kappa )
             if ( kGeometric < config_GM_spatially_variable_min_kappa ) then
                  kGeometric = config_GM_spatially_variable_min_kappa
             else if ( kGeometric > config_GM_spatially_variable_max_kappa ) then
                       kGeometric = config_GM_spatially_variable_max_kappa
             end if

             gmKappaScaling(k,iEdge) = kGeometric / config_GM_spatially_variable_max_kappa

             !gmKappaScaling(k,iEdge) = config_GM_spatially_variable_min_kappa / config_GM_spatially_variable_max_kappa

             !if ( k == 3 ) then
             !  print*,k,kMLD,kGeometric,eddyEnergyEdge, M2NInteg,structFunction(k)
             !endif

          end do ! k

             !-------------------------------------------------------------------------
       end do

       deallocate (structFunction)

    end subroutine ocn_GM_compute_kappa_geometric

!***********************************************************************
!
!  routine ocn_GM_eddyEnergy_tendency
!
!> \brief   Compute tendency of the eddy energy budget equation for
!           the GEOMETRIC parameterization
!> \details
!>  Compute tendency of the eddy energy budget equation for the
!   GEOMETRIC parameterization if 'geometric' GM is enabled.
!
!> \References
!> Marshall et al. (2012), Mak et al. (2018), Mak et al. (2022)
!-----------------------------------------------------------------------

    subroutine ocn_GM_eddyEnergy_tendency(statePool, meshPool, tendPool, &
                                          dt, timeLevelIn)!{{{

       !-----------------------------------------------------------------
       !
       ! input variables
       !
       !-----------------------------------------------------------------

       type(mpas_pool_type), intent(in) :: &
          statePool, &! pool containing state variables
          meshPool,  &! pool containing mesh quantities
          tendPool    ! pool containing tendencies

       integer, intent(in), optional :: &
          timeLevelIn          ! time level for state variables

       real (kind=RKIND), intent(in) :: &
          dt              !< [in] time step (sec) to move forward

       !-----------------------------------------------------------------
       !
       ! input/output variables
       !
       !-----------------------------------------------------------------

       !-----------------------------------------------------------------
       !
       ! local variables
       !
       !-----------------------------------------------------------------

       integer :: nCells, nEdges, nVertices
       integer :: timeLevel
       integer :: i, j, k, iCell, iEdge, iVertex, cell1, cell2
       integer :: kCell1,kCell2, kEdge,eoe
       integer, pointer :: nVertLevels, vertexDegree
       integer, pointer :: nCellsAll, nEdgesAll, nVerticesAll,nCellsOwned
       integer, dimension(:), pointer   :: minLevelEdgeBot, minLevelEdgeTop, &
                                           maxLevelEdgeTop, minLevelCell,    &
                                           maxLevelCell, nEdgesOnCell, nEdgesOnEdge, &
                                           nCellsArray, nEdgesArray, nVerticesArray
       integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnCell
       integer, dimension(:,:), pointer :: edgesOnEdge, edgesOnVertex, edgesOnCell, kiteIndexOnCell
       integer, dimension(:,:), pointer :: edgeSignOnCell, edgeSignOnVertex
       real(kind=RKIND), dimension(:,:), pointer :: layerThickness
       real(kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, kiteAreasOnVertex
       integer, dimension(:,:), pointer :: edgeMask, cellMask
       real(kind=RKIND), dimension(:), pointer :: normalBarotropicVelocity
       real(kind=RKIND), dimension(:), pointer :: totalParamEddyEnergy, tendTotalParamEddyEnergy
       real(kind=RKIND), dimension(:), pointer :: fEdge, latEdge
       real(kind=RKIND), dimension(:), pointer :: areaCell, dcEdge, dvEdge, areaTriangle, angleEdge
       real(kind=RKIND), dimension(:), pointer :: meshScalingDel2
       real(kind=RKIND), dimension(:), allocatable :: BVFSmoothed,tangentialFlux, advVertex
       real(kind=RKIND), dimension(:), allocatable :: longRossbyPhaseSpeed
       real(kind=RKIND) :: areaEdge, BVFref, BruntVaisalaFreqTop1, rtmp
       real(kind=RKIND) :: BVFinteg, c1, cp
       real(kind=RKIND) :: dcEdgeInv, eddyDiff2, gradEddyEnergy, hML, wgt
       real(kind=RKIND) :: invAreaCell, invAreaTriangle, r_tmp_div
       real(kind=RKIND) :: M4N2, gmKappaCell, flux, wgtNorm, gradBuoyCell
       real(kind=RKIND) :: dissipationCoef, dissipationRate
       real(kind=RKIND) :: slopeCell, sourceTerm, advectionTerm

       if (.not. config_use_gm .and. config_GM_closure /= 'geometric') return

       call mpas_timer_start('gm compute geometric tendency')

       if (present(timeLevelIn)) then
          timeLevel = timeLevelIn
       else
          timeLevel = 1
       end if

       call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                            normalBarotropicVelocity, timeLevel)
       call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
       call mpas_pool_get_array(statePool, 'totalParamEddyEnergy', totalParamEddyEnergy, timeLevel)

       call mpas_pool_get_array(tendPool, 'totalParamEddyEnergy', tendTotalParamEddyEnergy)

       call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
       call mpas_pool_get_array(meshPool, 'minLevelEdgeTop', minLevelEdgeTop)
       call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
       call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
       call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
       call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
       call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
       call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(meshPool, 'nEdgesOnEdge', nEdgesOnEdge)
       call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
       call mpas_pool_get_array(meshPool, 'edgesOnEdge',edgesOnEdge)
       call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
       call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)
       call mpas_pool_get_array(meshPool, 'edgeSignOnVertex', edgeSignOnVertex)
       call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)
       call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
       call mpas_pool_get_array(meshPool, 'cellMask', cellMask)
       call mpas_pool_get_array(meshPool, 'weightsOnEdge',weightsOnEdge)
       call mpas_pool_get_array(meshPool, 'areaTriangle',areaTriangle)
       call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex',kiteAreasOnVertex)
       call mpas_pool_get_array(meshPool, 'kiteIndexOnCell',kiteIndexOnCell)
       call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)

       call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
       call mpas_pool_get_dimension(meshPool, 'nEdgesArray', nEdgesArray)
       call mpas_pool_get_dimension(meshPool, 'nVerticesArray', nVerticesArray)
       call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
       call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
       call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
       call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
!      call mpas_pool_get_dimension(meshPool, 'nEdgesAll',nEdgesAll)
!      call mpas_pool_get_dimension(meshPool, 'nVerticesAll',nVerticesAll)
!      call mpas_pool_get_dimension(meshPool, 'nCellsOwned',nCellsOwned)

       nEdges = nEdgesArray(3)
       nVertices = nVerticesArray(2)
       nCells = nCellsArray(2)

       allocate(tangentialFlux(nEdges))
       allocate(advVertex(nVertices))


       !-------------------------------------------------------------------------
       ! Initialization
       !-------------------------------------------------------------------------
       do iCell = 1,nCellsArray(3)
          tendTotalParamEddyEnergy(i) = 0.0_RKIND
       end do

       !-------------------------------------------------------------------------
       ! Advection term
       !-------------------------------------------------------------------------

       ! 2. The intrinsic long Rossby phase speed varies with latitude
       allocate(longRossbyPhaseSpeed(nEdges))

       do i = 1,nEdges
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          ! 2.1. Compute the eigenvalue associated with the first baroclinic mode
          !      (Eq. 6.11.8 of Gill (1982), Eq. 2.2 of Chelton et al. (1998))

          ! Surface level
          hML = 0.5_RKIND*layerThickEdgeMean(minLevelEdgeBot(iEdge),iEdge)
          BruntVaisalaFreqTop1 =  sqrt(0.5_RKIND*(max(BruntVaisalaFreqTop(minLevelCell(cell1)+1,cell1),1.0E-14_RKIND) &
                               +  max(BruntVaisalaFreqTop(minLevelCell(cell2)+1,cell2),1.0E-14_RKIND)))
          BVFinteg = max(BruntVaisalaFreqTop1, 1.0E-14_RKIND) * hML

          do k = minLevelEdgeBot(iEdge)+1,maxLevelEdgeBot(iEdge)
             hML = 0.5_RKIND*(layerThickEdgeMean(k,iEdge) + layerThickEdgeMean(k-1,iEdge))
             BruntVaisalaFreqTop1 =  sqrt(0.5_RKIND*(max(BruntVaisalaFreqTop(k,cell1),1.0E-14_RKIND) &
                                  +  max(BruntVaisalaFreqTop(k,cell2),1.0E-14_RKIND)))
             BVFinteg = BVFinteg + max(BruntVaisalaFreqTop1, 1.0E-14_RKIND) * hML
          end do

          c1 = BVFinteg / pi

          ! 2.2. Compute the long-phase speed that the total eddy energy is to be advected
          !      (e.g., Eq. 12.3.13 of Gill (1982))
          cp = c1**2.0_RKIND * cos(latEdge(iEdge)) &
             / (2.0_RKIND*omega*sphere_radius*sin(latEdge(iEdge))**2.0_RKIND)

          ! 2.3. Cap the long-phase speed bounded by above by the magnitude tropical
          !      planetary wave phase speed (e.g., Eq. 12.3.14 of Gill (1982)), as it diverges
          !      at the equator
          longRossbyPhaseSpeed(iEdge) = min(abs(cp),abs(c1/3.0_RKIND))

          ! 2.4. Zonal vector to normal vector (zonal vector only)
          longRossbyPhaseSpeed(iEdge) = longRossbyPhaseSpeed(iEdge) * cos(angleEdge(iEdge))
       end do ! iEdge


       !-----------------------------------
       ! 1. Barotropic advection
       !    - To obtain smoothed solution for more stable time integration of
       !      the total parameterized eddy energy, compute flux at vertices
       !      and remap the flux at vertices to cells

       !-----------------------------------
       ! 1.1. Computing tangential component of the advection term
       do iEdge = 1,nEdges
          tangentialFlux(iEdge) = 0.0_RKIND
          do i = 1, nEdgesOnEdge(iEdge)
             eoe = edgesOnEdge(i,iEdge)
             if (minLevelEdgeBot(eoe).eq.0) cycle
             cell1 = cellsOnEdge(1,eoe)
             cell2 = cellsOnEdge(2,eoe)

             tangentialFlux(iEdge) = &
             tangentialFlux(iEdge) + weightsOnEdge(i,iEdge)* &
                                     ! [Normal BTR Vel] - [Long Rossby C]
                                     (  normalBarotropicVelocity(eoe)    &
                                      - longRossbyPhaseSpeed(iEdge)  ) * &
                                     edgeMask(minLevelEdgeBot(iEdge),iEdge) !* &
          end do
       end do

       ! 1.2. Compute advection term at vertices (divergence at vertices)
       do iVertex = 1,nVertices
          invAreaTriangle = 1.0_RKIND / areaTriangle(iVertex)
          advVertex(iVertex) = 0.0_RKIND
          do i = 1, vertexDegree
             iEdge = edgesOnVertex(i,iVertex)
             if (minLevelEdgeBot(iEdge).eq.0) cycle
             cell1 = cellsOnEdge(1,iEdge)
             cell2 = cellsOnEdge(2,iEdge)
             rtmp = 0.5_RKIND*(totalParamEddyEnergy(cell1) +  &
                               totalParamEddyEnergy(cell2)) * &
                    edgeMask(minLevelEdgeBot(iEdge),iEdge)
             r_tmp_div = dcEdge(iEdge) * tangentialFlux(iEdge) * invAreaTriangle * rtmp
             advVertex(iVertex) = advVertex(iVertex) &
                                - edgeSignOnVertex(i,iVertex) * r_tmp_div
          end do
       end do

       ! 1.3. Remapping from verticies to cells
       do iCell = 1, nCells
          invAreaCell = 1.0_RKIND/areaCell(iCell)
          do i = 1,nEdgesOnCell(iCell)
             j = kiteIndexOnCell(i,iCell)
             iVertex = verticesOnCell(i, iCell)
             tendTotalParamEddyEnergy(iCell) = &
             tendTotalParamEddyEnergy(iCell) + kiteAreasOnVertex(j,iVertex)* &
                                     advVertex(iVertex)*invAreaCell
          end do
       end do

       !-----------------------------------

       do iCell = 1, nCells
          invAreaCell = 1.0_RKIND/areaCell(iCell)
          sourceTerm = 0.0_RKIND

          !-------------------------------------------------------------------------
          ! Source term
          !-------------------------------------------------------------------------
          !NOTE, gradDensityEddy and BFV are defined at the top cell interface
          !and not valid for minLevelCell, so properties for the top layer
          !rely on the values at (minLevelCell+1).

          ! Surface level
          k = minLevelCell(iCell)
          hML = 0.5_RKIND*layerThickness(k,iCell)
          BruntVaisalaFreqTop1 = BruntVaisalaFreqTop(k,iCell)
          BruntVaisalaFreqTop1 = max(BruntVaisalaFreqTop1, 1.0E-14_RKIND)

          ! M^2
          wgtNorm = 1.0E-14_RKIND
          gradBuoyCell = 0.0_RKIND
          gmKappaCell = 0.0_RKIND
          slopeCell = 0.0_RKIND

          ! Interp edge to cell
          do i = 1, nEdgesOnCell(iCell)
             iEdge = edgesOnCell(i, iCell)
             if (maxLevelEdgeTop(iEdge).eq.0) cycle
             kEdge = minLevelEdgeBot(iEdge)
             wgt = 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge)
             wgtNorm = wgtNorm + wgt


             ! slope (=M^2/N^2)
             slopeCell =  slopeCell + wgt * 0.25_RKIND            &
                       * (  abs(slopeTriadUp(  k,1,iEdge))  &
                          + abs(slopeTriadUp(  k,2,iEdge))  &
                          + abs(slopeTriadDown(k,1,iEdge))  &
                          + abs(slopeTriadDown(k,2,iEdge)) )

             gradBuoyCell = gradBuoyCell + &
                            wgt * abs(gradDensityEddy(kEdge,iEdge) * gravity / rho_sw)
             ! GM Kappa
             gmKappaCell  = gmKappaCell + &
                            wgt * gmKappaScaling(kEdge,iEdge) * gmBolusKappa(iEdge)
          end do ! i

          gradBuoyCell = gradBuoyCell / wgtNorm
          gmKappaCell  =  gmKappaCell / wgtNorm
          slopeCell  =  slopeCell / wgtNorm

          ! M^4 / N^2 = slope^2 (=M^4/N^4) * BVF (=N^2)
          M4N2 = slopeCell**2.0_RKIND * BruntVaisalaFreqTop1
          !M4N2 = gradBuoyCell**2.0_RKIND / BruntVaisalaFreqTop1
          M4N2 = min(M4N2,1E-10_RKIND)

          sourceTerm = (gmKappaCell * M4N2) * hML
          !print*, k, gradBuoyCell**2.0_RKIND,BruntVaisalaFreqTop1,gmKappaCell,M4N2,(gmKappaCell * M4N2) * hML, sourceTerm

          do k = minLevelCell(iCell)+1, maxLevelCell(iCell)
             hML = layerThickness(k,iCell)

             ! N^2
             BruntVaisalaFreqTop1 = BruntVaisalaFreqTop(k,iCell)
             BruntVaisalaFreqTop1 = max(BruntVaisalaFreqTop1, 1.0E-14_RKIND)

             ! M^2
             wgtNorm      = 1.0E-14_RKIND
             gradBuoyCell = 0.0_RKIND
             gmKappaCell  = 0.0_RKIND
             slopeCell    = 0.0_RKIND
             ! Interp edge to cell
             do i = 1, nEdgesOnCell(iCell)
                iEdge = edgesOnCell(i, iCell)
                if (maxLevelEdgeTop(iEdge).eq.0) cycle
                kEdge = k
                wgt = 0.25_RKIND * dcEdge(iEdge) * dvEdge(iEdge)
                wgtNorm = wgtNorm + wgt

                ! slope (=M^2/N^2)
                slopeCell =  slopeCell + wgt * 0.25_RKIND            &
                          * (  abs(slopeTriadUp(  kEdge,1,iEdge))  &
                             + abs(slopeTriadUp(  kEdge,2,iEdge))  &
                             + abs(slopeTriadDown(kEdge,1,iEdge))  &
                             + abs(slopeTriadDown(kEdge,2,iEdge)) )

                gradBuoyCell = gradBuoyCell + &
                               wgt * abs(gradDensityEddy(kEdge,iEdge) * gravity / rho_sw)
                ! GM Kappa
                gmKappaCell  = gmKappaCell + &
                               wgt * gmKappaScaling(kEdge,iEdge) * gmBolusKappa(iEdge)

             !print*, slopeTriadUp(  kEdge,1,iEdge),  slopeTriadUp(  kEdge,2,iEdge), slopeTriadDown(kEdge,1,iEdge), slopeTriadDown(kEdge,2,iEdge)
             end do ! i
             gradBuoyCell = gradBuoyCell / wgtNorm
             gmKappaCell  =  gmKappaCell / wgtNorm
             slopeCell  =  slopeCell / wgtNorm


             ! M^4 / N^2 = slope^2 (=M^4/N^4) * BVF (=N^2)
             M4N2 = slopeCell**2.0_RKIND * BruntVaisalaFreqTop1
             !M4N2 = gradBuoyCell**2.0_RKIND / BruntVaisalaFreqTop1
             M4N2 = min(M4N2,1E-10_RKIND)

             sourceTerm = sourceTerm + (gmKappaCell * M4N2) * hML

          end do ! k

          tendTotalParamEddyEnergy(iCell) = tendTotalParamEddyEnergy(iCell) + sourceTerm

          !-------------------------------------------------------------------------
          ! Dissipation term
          !-------------------------------------------------------------------------

          ! Linear dissipation coefficient
          ! TODO: Need to implement via config
          !dissipationCoef = 100.0_RKIND !v01
          !dissipationCoef = 365.0_RKIND !v02
          dissipationCoef = 90.0_RKIND !v03
          !dissipationCoef = 160.0_RKIND !v04
          dissipationRate = 1.0_RKIND / (dissipationCoef*86400.0_RKIND)
          tendTotalParamEddyEnergy(iCell) = tendTotalParamEddyEnergy(iCell) - &
                                  dissipationRate * (totalParamEddyEnergy(iCell) - 0.1_RKIND)


          !-------------------------------------------------------------------------
          ! Diffusion term
          !-------------------------------------------------------------------------
          ! TODO: Need to implement via config
          !eddyDiff2 = 500.0_RKIND !v01
          !eddyDiff2 = 200.0_RKIND !v02
          eddyDiff2 = 600.0_RKIND !v03
          !eddyDiff2 = 300.0_RKIND !v04


          flux = 0.0_RKIND
          do i = 1, nEdgesOnCell(iCell)
             iEdge = edgesOnCell(i, iCell)
             !if (maxLevelEdgeTop(iEdge).eq.0) cycle
             cell1 = cellsOnEdge(1,iEdge)
             cell2 = cellsOnEdge(2,iEdge)

             rtmp = meshScalingDel2(iEdge) * eddyDiff2 * dvEdge(iEdge) / dcEdge(iEdge)

             ! \kappa_2 \nabla \phi on edge
             gradEddyEnergy = totalParamEddyEnergy(cell2) - totalParamEddyEnergy(cell1)

             ! div(h \kappa_2 \nabla \phi) at cell center
             flux = gradEddyEnergy * rtmp

             tendTotalParamEddyEnergy(iCell) = tendTotalParamEddyEnergy(iCell) - &
                                     edgeSignOnCell(i, iCell) * flux * invAreaCell
          end do

       end do ! iCell

       deallocate(tangentialFlux)
       deallocate(advVertex)
       deallocate(longRossbyPhaseSpeed)

       call mpas_timer_stop('gm compute geometric tendency')

    end subroutine ocn_GM_eddyEnergy_tendency!}}}

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

      real(kind=RKIND) :: sphereRadius, coef

      type(block_type), pointer :: block
      type(mpas_pool_type), pointer :: meshPool, statePool
      real(kind=RKIND), dimension(:), pointer :: dcEdge

      integer :: iCell, iEdge, k
      integer, pointer :: nCells, nEdges
      integer, pointer :: nVertLevels
      real(kind=RKIND), pointer :: sphere_radius
      real(kind=RKIND), dimension(:), pointer   :: latEdge, fEdge
      integer, dimension(:, :), pointer :: cellsOnEdge
      !--------------------
      integer :: cell1, cell2
      real(kind=RKIND) :: eddyEnergyH, BruntVaisalaFreqTopEdge, gradBuoy, hML
      integer, dimension(:), pointer  :: minLevelCell, maxLevelCell
      real(kind=RKIND), dimension(:,:), pointer   :: layerThickness
      real(kind=RKIND), dimension(:), pointer   :: totalParamEddyEnergyCur
      real(kind=RKIND), dimension(:), pointer   :: totalParamEddyEnergyNew
      !--------------------

      err = 0

      if ((.not. config_use_GM) .and. (.not. config_use_Redi)) return

      block => domain%blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block%structs, 'state', statePool)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'fEdge', fEdge)

         if (config_Redi_use_slope_taper) then
            slopeTaperFactor = 1.0_RKIND
         else
            slopeTaperFactor = 0.0_RKIND
         end if

         if (config_Redi_use_surface_taper) then
            sfcTaperFactor = 1.0_RKIND
         else
            sfcTaperFactor = 0.0_RKIND
         end if

         local_config_GM_constant_bclModeSpeed = config_GM_constant_bclModeSpeed

         if (config_GM_minBclModeSpeed_method=='constant') then
            gm_minBclModeSpeed_compute_on = 0.0_RKIND
            gm_minBclModeSpeed_constant = config_GM_constant_bclModeSpeed
         else if (config_GM_minBclModeSpeed_method=='computed') then
            gm_minBclModeSpeed_compute_on = 1.0_RKIND
            gm_minBclModeSpeed_constant = 0.0_RKIND
         else
            call mpas_log_write('config_GM_minBclModeSpeed_method must be constant or computed', &
                                 MPAS_LOG_CRIT)
         end if

         RediGMinitValue = 1.0_RKIND
         !compute beta
         call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
         call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
         !$omp parallel
         !$omp do schedule(runtime)
         do iEdge=1,nEdges
            betaEdge(iEdge) = 2.0_RKIND*omega*cos(latEdge(iEdge)) / sphere_radius
         end do
         !$omp end do
         !$omp end parallel

         if (config_GM_closure == 'constant') then
            local_config_GM_lat_variable_c2 = .false.
            local_config_GM_kappa_lat_depth_variable = .false.
            local_config_GM_compute_Visbeck = .false.
            local_config_GM_compute_EdenGreatbatch = .false.
            local_config_GM_diagnostic_geometric = .false.
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               gmBolusKappa(iEdge) = config_GM_constant_kappa
            end do
            !$omp end do
            !$omp end parallel
            RediGMinitValue = 1.0_RKIND
         else if (config_GM_closure == 'N2_dependent') then
            local_config_GM_lat_variable_c2 = .true.
            local_config_GM_kappa_lat_depth_variable = .true.
            local_config_GM_compute_Visbeck = .false.
            local_config_GM_compute_EdenGreatbatch = .false.
            local_config_GM_diagnostic_geometric = .false.

            RediGMinitValue = 0.0_RKIND
            ! for N2 dependence, we still assign Kappa as a constant to be multiplied by N2 scaling.
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               gmBolusKappa(iEdge) = config_GM_spatially_variable_max_kappa
            end do
            !$omp end do
            !$omp end parallel
         else if (config_GM_closure == 'Visbeck'.or. &
                  config_GM_closure == 'visbeck') then
            local_config_GM_lat_variable_c2 = .true.
            local_config_GM_kappa_lat_depth_variable = .false.
            local_config_GM_compute_Visbeck = .true.
            local_config_GM_compute_EdenGreatbatch = .false.
            local_config_GM_diagnostic_geometric = .false.

            call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
            call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge=1,nEdges
              gmBolusKappa(iEdge) = config_GM_spatially_variable_max_kappa
            end do
            !$omp end do
            !$omp end parallel
         else if (config_GM_closure == 'EdenGreatbatch'.or. &
                  config_GM_closure == 'edenGreatbatch') then
            local_config_GM_lat_variable_c2 = .true.
            local_config_GM_kappa_lat_depth_variable = .false.
            local_config_GM_compute_Visbeck = .false.
            local_config_GM_compute_EdenGreatbatch = .true.
            local_config_GM_diagnostic_geometric = .false.
            call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
            call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge=1,nEdges
              gmBolusKappa(iEdge) = config_GM_spatially_variable_max_kappa
            end do
            !$omp end do
            !$omp end parallel
         else if (config_GM_closure == 'geometric') then
            local_config_GM_lat_variable_c2 = .false.
            local_config_GM_kappa_lat_depth_variable = .false.
            local_config_GM_compute_Visbeck = .false.
            local_config_GM_compute_EdenGreatbatch = .false.
            local_config_GM_diagnostic_geometric = .true.
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               gmBolusKappa(iEdge) = config_GM_spatially_variable_max_kappa
            end do
            !$omp end do
            !$omp end parallel

            ! Initialization for the geometric option --------------------------------
            ! Initialize vertically-integrated eddy energy

            call mpas_pool_get_array(statePool, 'totalParamEddyEnergy', totalParamEddyEnergyCur, 1)
            call mpas_pool_get_array(statePool, 'totalParamEddyEnergy', totalParamEddyEnergyNew, 2)
            call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
            call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            !TODO: Add 'restart' option (if (.not. restart))
            do iCell = 1,nCells
               totalParamEddyEnergyCur(iCell) = 1.0E-1_RKIND
               totalParamEddyEnergyNew(iCell) = 1.0E-1_RKIND
            end do

         else
            call mpas_log_write('Invalid choice of config_GM_closure.', MPAS_LOG_CRIT)
            err = 1
            call mpas_dmpar_finalize(domain%dminfo)
         end if

         ! Initialize horizontal taper
         if (config_GM_horizontal_taper == 'none' .or. &
             config_GM_horizontal_taper == 'RossbyRadius') then
            ! For 'RossbyRadius', the taper is recomputed at every time step.
            ! but just set to one here
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1,nEdges
               gmHorizontalTaper(iEdge) = 1.0_RKIND
            end do
            !$omp end do
            !$omp end parallel
         else if (config_GM_horizontal_taper == 'ramp') then
            coef = 1.0_RKIND &
                    /(config_GM_horizontal_ramp_max - config_GM_horizontal_ramp_min)
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge = 1, nEdges
               if (dcEdge(iEdge) <= config_GM_horizontal_ramp_min) then
                  gmHorizontalTaper(iEdge) = 0.0_RKIND
               else if (dcEdge(iEdge) >= config_GM_horizontal_ramp_max) then
                  gmHorizontalTaper(iEdge) = 1.0_RKIND
               else
                  gmHorizontalTaper(iEdge) = coef &
                     *(dcEdge(iEdge) - config_GM_horizontal_ramp_min)
               end if
            end do
            !$omp end do
            !$omp end parallel
         else
            call mpas_log_write('Invalid choice of config_GM_horizontal_taper.', MPAS_LOG_CRIT)
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
