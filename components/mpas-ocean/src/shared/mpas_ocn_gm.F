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

   use ocn_redi
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

   public :: ocn_GM_compute_Bolus_velocity, &
             ocn_GM_init,                   &
             ocn_GM_add_to_transport_vel

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
   logical :: local_config_GM_kappa_lat_depth_variable_orig
   real(kind=RKIND) :: rediGMinitValue
   real(kind=RKIND) :: local_config_GM_constant_bclModeSpeed
   real(kind=RKIND) :: gm_minBclModeSpeed_constant, gm_minBclModeSpeed_compute_on

   contains
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
                                            meshPool, timeLevelIn)
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

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND), dimension(:, :), pointer :: &
                                                    normalVelocity, &
                                                    layerThickness

      real(kind=RKIND), dimension(:), pointer   :: &
                ssh, fEdge
      real(kind=RKIND), dimension(:), pointer   :: areaCell, dcEdge, dvEdge, tridiagA, tridiagB, tridiagC, rightHandSide
      integer, dimension(:), pointer   :: minLevelEdgeBot, maxLevelEdgeTop, minLevelCell, maxLevelCell, nEdgesOnCell
      integer, dimension(:, :), pointer :: cellsOnEdge, edgesOnCell
      integer                          :: i, k, iEdge, cell1, cell2, iCell, N, iCellSelf, maxLocation
      real(kind=RKIND)                 :: h1, h2, areaEdge, BruntVaisalaFreqTopEdge, rtmp
      real(kind=RKIND)                 :: sum_hN, countN2, maxN, ltSum
      real(kind=RKIND)                 :: sumRi, RiTopOfEdge, zEdge
      real(kind=RKIND) :: dcEdgeInv, invAreaCell
      real(kind=RKIND) :: lt1, lt2, c_min
      real(kind=RKIND) :: sigma, Lr, Length, L_rhines, shearEdgeInv
      real (kind=RKIND) :: eddyLength, eddyTime
      real(kind=RKIND), dimension(:), allocatable :: RossbyRadiusTaper
      real(kind=RKIND) :: c_Visbeck   ! baroclinic wave speed from Visbeck parameterization
      real(kind=RKIND) :: c_mode1, dc_Lr
      ! Dimensions
      integer :: kmin, kmax, nCells, nEdges
      integer, pointer :: nVertLevels
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray

      type(mpas_pool_type), pointer :: tracersPool
      real(kind=RKIND), dimension(:, :, :), pointer :: activeTracers
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

      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
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
      !$omp end parallel

      !--------------------------------------------------------------------
      !
      ! Compute stream function and Bolus velocity for Bolus part of GM
      !
      !--------------------------------------------------------------------

      nEdges = nEdgesArray(3)
      if (config_GM_horizontal_taper == 'RossbyRadius'.or. &
         config_Redi_horizontal_taper == 'RossbyRadius') then
         allocate (RossbyRadiusTaper(nEdges))
         call mpas_pool_get_array(meshPool, 'fEdge', fEdge)

         nCells = nCellsArray(4)

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

      if (config_use_redi) then
         call ocn_redi_compute_isopycnal_slopes(statePool,timeLevel)
      end if

      if (config_use_GM) then

         ! For config_GM_closure = 'N2_dependent' use a scaling to taper gmBolusKappa
         ! based on stratification relative to the maximum in the column
         if (local_config_GM_kappa_lat_depth_variable) then
            !$omp parallel
            !$omp do schedule(runtime) private(i, k, BruntVaisalaFreqTopEdge, maxN, cell1, cell2)
            do iEdge=1,nEdges
               k=minLevelEdgeBot(iEdge)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               i = max(1,min(int(indexBoundaryLayerDepth(cell1)),int(indexBoundaryLayerDepth(cell2))))
               BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTop(i,cell1), &
                                                    BruntVaisalaFreqTop(i,cell2))
               maxN = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)
               !Scan from the boundary layer depth only
               do k = i+1, maxLevelEdgeTop(iEdge)
                 !According to Danabasoglu and Marshall (2007) take max BVF
                 BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTop(k,cell1), &
                                                      BruntVaisalaFreqTop(k,cell2))
                 maxN = max(maxN,max(BruntVaisalaFreqTopEdge, 0.0_RKIND))
               enddo

               do k = minLevelEdgeBot(iEdge), i
                  gmKappaScaling(k,iEdge) = 1.0_RKIND
               end do

               do k = i+1, maxLevelEdgeTop(iEdge)
                  BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTop(k,cell1), &
                                                       BruntVaisalaFreqTop(k,cell2))
                  BruntVaisalaFreqTopEdge = max(BruntVaisalaFreqTopEdge, 0.0_RKIND)

                  gmKappaScaling(k, iEdge) = min(max(config_GM_spatially_variable_min_kappa /         &
                                                     config_GM_spatially_variable_max_kappa,          &
                                                     BruntVaisalaFreqTopEdge/(maxN + 1.0E-10_RKIND)), &
                                                     1.0_RKIND)
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
              do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)-1
                lt1 = 0.5_RKIND*(layerThickness(k,cell1) + layerThickness(k-1,cell1))
                lt2 = 0.5_RKIND*(layerThickness(k,cell2) + layerThickness(k-1,cell2))
                sum_hN = sum_hN + 0.5_RKIND*(lt1*max(BruntVaisalaFreqTop(k,cell1),0.0_RKIND) + &
                                lt2*max(BruntVaisalaFreqTop(k,cell2),0.0_RKIND))
                sumRi = sumRi + 0.5_RKIND*(lt1*max(RiTopOfCell(k,cell1),0.0_RKIND) + &
                                    lt2*max(RiTopOfCell(k,cell2),0.0_RKIND))
                ltsum = ltsum + 0.5_RKIND*(lt1+lt2)
                countN2 = countN2 + 1
              end do

              if (countN2 > 0) then
                c_Visbeck = cGMphaseSpeed(iEdge)
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

      if (config_use_redi) then
         call ocn_redi_compute_isopycnal_slopes(statePool, timeLevel)
      end if

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
      type(mpas_pool_type), pointer :: meshPool
      real(kind=RKIND), dimension(:), pointer :: dcEdge

      integer :: iEdge
      integer, pointer :: nCells, nEdges
      integer, pointer :: nVertLevels
      real(kind=RKIND), pointer :: sphere_radius
      real(kind=RKIND), dimension(:), pointer   :: latEdge, fEdge
      integer, dimension(:, :), pointer :: cellsOnEdge

      err = 0

      if ((.not. config_use_GM) .and. (.not. config_use_Redi)) return

      block => domain%blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'fEdge', fEdge)

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

            call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
            call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
            !$omp parallel
            !$omp do schedule(runtime)
            do iEdge=1,nEdges
              gmBolusKappa(iEdge) = config_GM_spatially_variable_max_kappa
            end do
            !$omp end do
            !$omp end parallel
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
