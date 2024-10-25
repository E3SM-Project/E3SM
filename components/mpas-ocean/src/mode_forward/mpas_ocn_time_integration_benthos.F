! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_benthos
!
!> \brief MPAS ocean  time integration for benthos spinup
!> \author Nicole Jeffery
!> \date   March 2022
!> \details
!>  This module contains the routine for benthos integration only
!
!-----------------------------------------------------------------------

module ocn_time_integration_benthos

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer
   use mpas_threading
   use mpas_timekeeping
   use mpas_log

   use ocn_config
   use ocn_mesh
   use ocn_tendency_benthos
   use ocn_diagnostics_variables
   use ocn_diagnostics
   use ocn_gm

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_time_average_coupled

   use ocn_effective_density_in_land_ice
   use ocn_transport_tests
   use ocn_tracer_benthos

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

   public :: ocn_time_integrator_benthos, &
             ocn_time_integration_benthos_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   character (len=*), parameter :: &
      subcycleGroupName = 'subcycleFields', &! name for cached halo exch
      finalBtrGroupName = 'finalBtrFields'   ! name for cached halo exch

   logical :: &
      unsplit           ! flag for unsplit_explicit time option

   integer, dimension(:), allocatable :: &
      numClinicIterations  ! number of baroclinic iterations for each
                           ! outer timestep iteration

   integer :: &
      neededHalos,     &! number of halo levels needed
      haloDecrement,   &! number of halos to invalidate each cycle
      numTSIterations, &! number of outer timestep iterations
      nBtrSubcycles     ! number of barotropic subcycles

   real (kind=RKIND) :: &
      useVelocityCorrection,&! mask for velocity correction
      splitFact          ! mask for terms in split or unsplit cases

!***********************************************************************

   contains

!***********************************************************************
!
!  ocn_time_integrator_benthos
!
!> \brief MPAS ocean benthos time integration scheme
!> \author Nicole Jeffery
!> \date   March 2022
!> \details
!>  This routine integrates the benthos sub-model forward one master time
!>  step (dt)
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_benthos(domain, dt)!{{{

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         dt              !< [in] time step (sec) to move forward

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain  !< [inout] model state to advance forward

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: &
         block ! structure with subdomain data

      type (mpas_pool_type), pointer :: &
         statePool,         &! structure holding state variables
         tracersPool,       &! structure holding tracers
         meshPool,          &! structure holding mesh variables
         verticalMeshPool,  &! structure holding vertical mesh variables
         tendPool,          &! structure holding tendencies
         tracersTendPool,   &! structure holding tracer tendencies
         forcingPool,       &! structure holding forcing variables
         scratchPool,       &! structure holding temporary variables
         swForcingPool,     &! structure holding short-wave forcing vars
         diagnosticsPool,   &! structure holding tracer diagnostics
         seafloorTendPool,  &! structure holding benthos tendencies
         seafloorTracersPool, &! structure holding benthos tracers
         tracersIdealAgeFieldsPool, & ! structure holding idealAge-related fields
         benthosAuxiliaryPool ! structure holding index to benthic mesh

      logical :: &
         activeTracersOnly   ! only compute tendencies for active tracers

      logical, pointer :: &
         config_use_tracerGroup,    & ! flag for using each tracer group
         config_use_benthosTracers    ! flag for using benthos tracers

      integer :: &
         iCell, iEdge, k, &! loop iterators for cell, edge, vert loops
         iCellBenthic,    & ! loop iterator for benthic cells
         nCells, nEdges,  & ! number of cells or edges (incl halos)
         i,j,             &! generic loop iterators
         cell1, cell2,    &! neighbor cell addresses
         eoe,             &! index for edge on edge
         err,             &! local error flag
         splitExplicitStep, &! loop index for outer timestep loop
         oldBtrSubcycleTime,&! time index for old barotropic value
         newBtrSubcycleTime,&! time index for new barotropic value
         uPerpTime,       &! time index to use for uPerp calculation
         BtrCorIter,      &! loop index for barotropic coriolis cycle
         stage1_tend_time,&! time index for stage 1 tendencies
         edgeHaloComputeCounter, &! halo counters to reduce
         cellHaloComputeCounter   ! halo updates during cycling

      real (kind=RKIND) :: &
         normalThicknessFluxSum, &! sum of thickness flux in column
         thicknessSum,     &! sum of thicknesses in column
         flux,             &! temp for computing flux for barotropic
         sshEdge,          &! sea surface height at edge
         CoriolisTerm,     &! temp for computing coriolis term (fuperp)
         normalVelocityCorrection, &! velocity correction
         temp,             &! temp for holding vars at new time
         temp_h,           &! temporary for phi
         temp_mask,        &! temporary mask (edgeMask)
         lat,              &! cell latitude
         sshCell1,         &! sea sfc height in neighboring cells
         sshCell2

      real (kind=RKIND), dimension(:), allocatable:: &
         uTemp

      real (kind=RKIND), dimension(:), allocatable :: &
         btrvel_temp

      ! State Array Pointers
      real (kind=RKIND), dimension(:), pointer :: &
         sshSubcycleCur,              &! subcycl ssh    at current time
         sshSubcycleNew,              &! subcycl ssh    at new     time
         sshSubcycleCurWithTides,     &! subcycl ssh    at current time
         sshSubcycleNewWithTides,     &! subcycl ssh    at new     time
         normalBarotropicVelocitySubcycleCur,&! barotropic vel subcyc
         normalBarotropicVelocitySubcycleNew,&! at current, new times
         sshCur,                      &! sea sfc height at current time
         sshNew,                      &! sea sfc height at new     time
         normalBarotropicVelocityCur, &! barotropic vel at current time
         normalBarotropicVelocityNew   ! barotropic vel at new     time

      real (kind=RKIND), dimension(:,:), pointer :: &
         normalBaroclinicVelocityCur, &! baroclinic vel at current time
         normalBaroclinicVelocityNew, &! baroclinic vel at new     time
         normalVelocityCur,           &! full velocity  at current time
         normalVelocityNew,           &! full velocity  at new     time
         layerThicknessCur,           &! layer thick    at current time
         layerThicknessNew,           &! layer thick    at new     time
         highFreqThicknessCur,        &! high frq thick at current time
         highFreqThicknessNew,        &! high frq thick at new     time
         lowFreqDivergenceCur,        &! low frq div    at current time
         lowFreqDivergenceNew,        &! low frq div    at new     time
         tracerGroupIdealAgeMask       ! mask for resetting surface ideal age

      type (field1DReal), pointer :: &
         effectiveDensityField         ! field pointer for halo update

      ! State/Tracer arrays, info
      real (kind=RKIND), dimension(:,:,:), pointer ::      &!
         tracersGroupCur,      &! old, new tracer arrays
         tracersGroupNew,      &
         benthosTracersCur,    &! old benthosTracer arrays
         benthosTracersNew      ! new benthosTracer arrays

      integer, pointer :: &
         startIndex, endIndex, &! start, end index for tracer groups
         indexSalinity,        &! tracer index for salinity
         nBenthicVertLevels     ! number of benthic vertical levels

      type (mpas_pool_iterator_type) :: &
         groupItr               ! iterator for tracer groups

      character (len=StrKIND) :: &
         modifiedGroupName,    &! constructed tracer group names
         configName             ! constructed config names for tracers

      real (kind=RKIND), dimension(:), pointer :: &
         benthosInterfaceLayerThickness ! benthos interface-grid layer thickness

      ! Tendency Array Pointers

      real (kind=RKIND), dimension(:), pointer :: &
         sshTend                ! sea-surface height tendency

      real (kind=RKIND), dimension(:,:), pointer :: &
         normalVelocityTend,    &! normal velocity tendency
         highFreqThicknessTend, &! thickness tendency for high-freq iter
         lowFreqDivergenceTend, &! thickness tendency for low freq
         layerThicknessTend      ! main thickness tendency

      real (kind=RKIND), dimension(:,:,:), pointer :: &
         tracersGroupTend,      &! tendencies for tracer groups
         activeTracersTend,     &! active tracer tendencies
         benthosTracersTend      ! benthos tracer tendencies

      ! Forcing pool
      real (kind=RKIND), dimension(:), pointer :: tidalPotentialEta

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Start main time step timer
      call mpas_timer_start("se timestep")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Prep variables before first iteration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("se prep")

      !*** Retrieve model state, pools
      !*** No longer support sub-blocks so this retrieves the only block

      block => domain%blocklist
      call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block%structs, 'verticalMesh', &
                                                 verticalMeshPool)
      call mpas_pool_get_subpool(block%structs, 'state', &
                                                 statePool)
      call mpas_pool_get_subpool(block%structs, 'forcing', &
                                                 forcingPool)
      call mpas_pool_get_subpool(block%structs, 'shortwave', &
                                                 swForcingPool)
      call mpas_pool_get_subpool(block%structs, 'diagnostics', &
                                                 diagnosticsPool)
      call mpas_pool_get_subpool(block%structs, 'tend', &
                                                 tendPool)
      call mpas_pool_get_subpool(block%structs, 'scratch', &
                                                 scratchPool)
      call mpas_pool_get_subpool(statePool, 'tracers', &
                                             tracersPool)
      call mpas_pool_get_subpool(tendPool,  'tracersTend', &
                                             tracersTendPool)

      call mpas_pool_get_config(domain % configs, 'config_use_benthosTracers', config_use_benthosTracers)
      if (config_use_benthosTracers) then    ! need to add benthosStorage
         call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
         call mpas_pool_get_subpool(forcingPool, 'benthosAuxiliary', benthosAuxiliaryPool)
         call mpas_pool_get_dimension(meshPool, 'nBenthicVertLevels', nBenthicVertLevels)
         call mpas_pool_get_subpool(statePool, 'seafloorTracers', seafloorTracersPool)
         call mpas_pool_get_subpool(tendPool, 'seafloorTend', seafloorTendPool)
         call mpas_pool_get_array(seafloorTracersPool, 'benthosTracers', benthosTracersCur,1)
         call mpas_pool_get_array(seafloorTracersPool, 'benthosTracers', benthosTracersNew,2)

         if (associated(benthosTracersCur) .and. associated(benthosTracersNew)) then
            !$omp parallel
            !$omp do schedule(runtime) private(k)
            do iCellBenthic = 1, nBenthicCells ! couple tracers to thickness
               do k = 1, nBenthicVertLevels+1
                  benthosTracersNew(:, k, iCellBenthic) = benthosTracersCur(:, k, iCellBenthic)
               end do
            end do
            !$omp end do
            !$omp end parallel
         end if
      end if

      call mpas_timer_stop("se prep")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large outer timestep iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do splitExplicitStep = 1, 1 !numTSIterations

         call mpas_timer_start('se loop')

         stage1_tend_time = min(splitExplicitStep,2)

         oldBtrSubcycleTime = 1
         newBtrSubcycleTime = 2
         activeTracersOnly = .false.

         call mpas_timer_start('se tracer tend')
         call ocn_tend_benthos(tendPool, statePool, forcingPool, &
                              meshPool, swForcingPool, diagnosticsPool, &
                              dt, activeTracersOnly, 2)

         call mpas_timer_stop('se tracer tend')

         ! update halo for tracer tendencies
         call mpas_timer_start("se halo tracers")

         if (config_use_benthosTracers) then
            call mpas_pool_get_subpool(tendPool, 'seafloorTend', seafloorTendPool)
            call mpas_dmpar_field_halo_exch(domain, 'benthosTracersTend')
         end if

         call mpas_timer_stop("se halo tracers")

         call mpas_timer_start('se loop fini')

         if (config_use_benthosTracers) then
            call mpas_pool_get_subpool(forcingPool, 'benthosAuxiliary', benthosAuxiliaryPool)
            call mpas_pool_get_subpool(statePool, 'seafloorTracers', seafloorTracersPool)
            call mpas_pool_get_subpool(tendPool, 'seafloorTend', seafloorTendPool)
            call mpas_pool_get_dimension(meshPool, 'nBenthicVertLevels', nBenthicVertLevels)
            call mpas_pool_get_array(seafloorTracersPool, 'benthosTracers', benthosTracersCur, 1)
            call mpas_pool_get_array(seafloorTracersPool, 'benthosTracers', benthosTracersNew, 2)
            call mpas_pool_get_array(seafloorTendPool, 'benthosTracersTend', benthosTracersTend)
            call mpas_pool_get_array(seafloorTendPool, 'benthosInterfaceLayerThickness', benthosInterfaceLayerThickness)
            if (associated(benthosTracersNew) .and. associated(benthosTracersCur) .and. associated(benthosTracersTend)) then
               !$omp parallel
               !$omp do schedule(runtime) private(k)
               do iCellBenthic = 1, nBenthicCells
                  do k = 1, nBenthicVertLevels+1
                     benthosTracersNew(:,k,iCellBenthic) = &
                        (benthosTracersCur(:,k,iCellBenthic) * &
                        benthosInterfaceLayerThickness(k) + dt* &
                        benthosTracersTend(:,k,iCellBenthic) )/ &
                        benthosInterfaceLayerThickness(k)
                  end do
               end do
               !$omp end do
               !$omp end parallel
            end if
         end if

         call mpas_timer_stop('se loop fini')
         call mpas_timer_stop('se loop')

      end do  ! outer timestep loop (splitExplicitStep)
      call mpas_timer_start('se fini')
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_stop('se fini')
      call mpas_timer_stop("se timestep")

   end subroutine ocn_time_integrator_benthos!}}}

!***********************************************************************
!
!  routine ocn_time_integration_benthos_init
!
!> \brief   Initialize benthos time stepping within ocean
!> \author  Nicole Jeffery
!> \date    March 2022
!> \details
!>  This routine initializes variables required for
!>  integrating the benthos forward in time.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integration_benthos_init(domain)!{{{

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain      !< [inout] data structure containing most variables

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: &
         statePool,      &! structure containing model state
         meshPool,       &! structure containing mesh fields
         tracersPool      ! structure containing tracer fields

      integer ::         &
         iCell, iEdge, k,&! loop indices for cell, edge and vertical
         kmax,           &! index of deepest active edge
         ierr,           &! local error flag
         cell1, cell2     ! neighbor cell indices across edge

      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) ::       &
         normalThicknessFluxSum, &! vertical sum of thick flux
         layerThicknessSum,      &! vertical sum of layer thickness
         layerThicknessEdge1      ! layer thickness on edge

      real (kind=RKIND), dimension(:), pointer :: &
         refBottomDepth,         &! reference bottom depth
         normalBarotropicVelocity ! normal barotropic velocity

      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThickness,          &! layer thickness cell center
         normalBaroclinicVelocity,&! normal baroclinic velocity
         normalVelocity            ! normal velocity (total)

      type (mpas_time_type) :: &
         nowTime          ! current model time

      type (mpas_timeInterval_type) :: &! various time step calculations
         fullTimeStep,       &! full model time step
         barotropicTimeStep, &! barotropic time step
         remainder,          &! remaining time after interval division
         zeroInterval         ! zero timestep for comparing remainder

      integer, pointer :: nVertLevels, nCells, nEdges

      integer (kind=I8KIND) :: nBtrSubcyclesI8

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      !*** Determine the number of barotropic subcycles based on the
      !*** ratio of time steps

      nowTime = mpas_get_clock_time(domain%clock, MPAS_NOW, ierr)

      call mpas_set_timeInterval( zeroInterval, S=0 )
      call mpas_set_timeInterval( fullTimeStep , timeString=config_dt )
      call mpas_set_timeInterval( barotropicTimeStep, &
                                             timeString=config_btr_dt )

      !*** transfer to I8 for division step
      nBtrSubcyclesI8 = nBtrSubcycles
      call mpas_interval_division( nowTime, fullTimeStep, &
                       barotropicTimeStep, nBtrSubcyclesI8, remainder )
      nBtrSubcycles = nBtrSubcyclesI8

      if ( remainder > zeroInterval ) then
         nBtrSubcycles = nBtrSubcycles + 1
      end if

      !*** Set mask for using velocity correction
      if (config_vel_correction) then
         useVelocityCorrection = 1.0_RKIND
      else
         useVelocityCorrection = 0.0_RKIND
      endif

      !*** Compute some halo needs to reduce the number
      !*** of halo updates during subcycling
      neededHalos = 1 + config_n_btr_cor_iter
      haloDecrement = mod(config_num_halos, neededHalos)

      !*** Determine the time integration type and set associated masks

      select case (trim(config_time_integrator))

      case ('unsplit_explicit')
         unsplit   = .true.
         splitFact = 0.0_RKIND

      case ('split_explicit')
         unsplit   = .false.
         splitFact = 1.0_RKIND
         call mpas_log_write( '*******************************************************************************')
         call mpas_log_write( 'The code is configured to spin-up the benthos using: $i barotropic subcycles', &
            intArgs=(/ nBtrSubcycles /) )
         call mpas_log_write( '*******************************************************************************')

      case default
         call mpas_log_write('Incorrect choice config_time_integrator',&
                             MPAS_LOG_CRIT)

      end select

      !*** set number of baroclinic iterations on each outer
      !*** time step iteration (number can be different on the
      !*** first and last time step iteration)

      numTSIterations     = config_n_ts_iter
      allocate(numClinicIterations(numTSIterations))

      numClinicIterations    = config_n_bcl_iter_mid ! most iterations
      numClinicIterations(1) = config_n_bcl_iter_beg ! first iteration
      numClinicIterations(numTSIterations)=config_n_bcl_iter_end !last

      !*** Initialize some variables if not a restart

      if (.not. config_do_restart) then

         !*** Initialize z-level mesh variables from h, read in from
         !*** input file.

         block => domain % blocklist
         call mpas_pool_get_subpool(block%structs, 'state', statePool)
         call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
         call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)

         call mpas_pool_get_array(statePool, 'layerThickness', &
                                              layerThickness, 1)
         call mpas_pool_get_array(statePool, 'normalVelocity', &
                                              normalVelocity, 1)
         call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                              normalBarotropicVelocity, 1)
         call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                              normalBaroclinicVelocity, 1)

         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
         call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
         call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)

         !*** Compute barotropic velocity at first timestep
         !*** This is only done upon start-up.
         if (unsplit) then

            do iEdge = 1, nEdgesAll
               normalBarotropicVelocity(  iEdge) = 0.0_RKIND
               normalBaroclinicVelocity(:,iEdge) = &
                         normalVelocity(:,iEdge)
            end do

         else ! split explicit

            if (config_filter_btr_mode) then
               do iCell = 1, nCells
                  layerThickness(minLevelCell(iCell),iCell) = refBottomDepth(minLevelCell(iCell))
               enddo
            endif

            do iEdge = 1, nEdges
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               kmax  = maxLevelEdgeTop(iEdge)

               ! normalBarotropicVelocity = sum(h*u)/sum(h) on each edge
               ! ocn_diagnostic_solve has not yet been called, so
               ! compute hEdge just for this edge.

               ! thicknessSum is initialized outside the loop because on
               ! land boundaries maxLevelEdgeTop=0, but we want to
               ! initialize thicknessSum with a nonzero value to avoid
               ! a NaN.
               layerThicknessEdge1 = 0.5_RKIND* &
                                     (layerThickness(minLevelCell(cell1),cell1) + &
                                      layerThickness(minLevelCell(cell2),cell2) )
               normalThicknessFluxSum = layerThicknessEdge1* &
                                        normalVelocity(minLevelEdgeBot(iEdge),iEdge)
               layerThicknessSum = layerThicknessEdge1

               do k=minLevelEdgeBot(iEdge)+1, kmax
                  layerThicknessEdge1 = 0.5_RKIND* &
                                       (layerThickness(k,cell1) + &
                                        layerThickness(k,cell2))

                  normalThicknessFluxSum = normalThicknessFluxSum + &
                                           layerThicknessEdge1* &
                                           normalVelocity(k,iEdge)
                  layerThicknessSum = layerThicknessSum + &
                                      layerThicknessEdge1

               enddo
               normalBarotropicVelocity(iEdge) = &
                     normalThicknessFluxSum/layerThicknessSum

               ! normalBaroclinicVelocity = normalVelocity -
               !     normalBarotropicVelocity
               do k = minLevelEdgeBot(iEdge), kmax
                  normalBaroclinicVelocity(k,iEdge) = &
                            normalVelocity(k,iEdge) - &
                    normalBarotropicVelocity(iEdge)
               enddo

               ! normalBaroclinicVelocity=0,
               ! normalVelocity=0 on land cells
               do k = kmax+1, nVertLevels
                  normalBaroclinicVelocity(k,iEdge) = 0.0_RKIND
                  normalVelocity(k,iEdge) = 0.0_RKIND
               enddo
            enddo ! edge loop

            if (config_filter_btr_mode) then
               ! filter normalBarotropicVelocity out of initial condition

               normalVelocity(:,:) = normalBaroclinicVelocity(:,:)
               normalBarotropicVelocity(:) = 0.0_RKIND

            endif

         endif ! split explicit

      end if ! not restart

      ! Create reusable halo exchange for subcycling
      call mpas_dmpar_exch_group_create(domain, subcycleGroupName)
      call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, &
                                           'sshSubcycle')
      call mpas_dmpar_exch_group_add_field(domain, subcycleGroupName, &
                                  'normalBarotropicVelocitySubcycle')
      call mpas_dmpar_exch_group_build_reusable_buffers(domain, &
                                                     subcycleGroupName)

      !-----------------------------------------------------------------

   end subroutine ocn_time_integration_benthos_init!}}}

!***********************************************************************

end module ocn_time_integration_benthos

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
