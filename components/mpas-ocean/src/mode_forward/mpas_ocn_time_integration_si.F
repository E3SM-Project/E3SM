! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_si
!
!> \brief MPAS ocean split-implicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011 (split explicit base code)
!
!> \author Hyun-Gyu Kang (Oak Ridge National Laboratory)
!> \date   September 2019 (split-implicit code)
!> \details
!>  This module contains the routines for the split-implicit
!>  time integration scheme based on the split-explicit code.
!>  Only stage 2 (barotropic mode) is changed from the explicit
!>  subcycling scheme to the split-implicit scheme.
!
!-----------------------------------------------------------------------

module ocn_time_integration_si

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
   use mpas_global_sum_mod

   use ocn_config
   use ocn_mesh
   use ocn_tendency
   use ocn_diagnostics_variables
   use ocn_diagnostics
   use ocn_gm
   use ocn_submesoscale_eddies

   use ocn_equation_of_state
   use ocn_vmix
   use ocn_vertical_advection
   use ocn_vertical_remap
   use ocn_time_average_coupled

   use ocn_effective_density_in_land_ice
   use ocn_surface_land_ice_fluxes
   use ocn_transport_tests

#ifdef MPAS_OPENACC
#ifdef USE_MAGMA
   use magma
   use iso_c_binding
#endif

#ifdef USE_CUBLAS
   use cublas
#endif
#endif  

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

   public :: ocn_time_integrator_si,                  &
             ocn_time_integration_si_init,            &
             ocn_time_integrator_si_preconditioner,   &
             ocn_time_integrator_si_variable_destroy

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   character (len=*), parameter :: &
      iterGroupName     = 'iterFields',   &! group name for halos
      finalBtrGroupName = 'finalBtrFields' ! group name for halos

   integer, dimension(:), allocatable :: &
      numClinicIterations  ! number of baroclinic iterations for each
                           ! outer timestep iteration

   integer :: &
      numTSIterations ! number of outer timestep iterations

   real (kind=RKIND) :: &
      useVelocityCorrection, & ! mask for velocity correction
      self_attraction_and_loading_beta

   ! Global variables for the split-implicit time stepper ---------------
   real (kind=RKIND), allocatable,dimension(:,:) :: &
      prec_ivmat    ! an inversed preconditioning matrix
   real (kind=RKIND), dimension(2) :: &
      SIcst_allreduce_local2, &! array for local  summations
      SIcst_allreduce_global2  ! array for global summations
   real (kind=RKIND), dimension(3) :: &
      SIcst_allreduce_local3, &! array for local  summations
      SIcst_allreduce_global3  ! array for global summations
   real (kind=RKIND), dimension(9) :: &
      SIcst_allreduce_local9, &! array for local  summations
      SIcst_allreduce_global9  ! array for global summations
   real (kind=RKIND), dimension(:,:), allocatable :: &
      globalReprodSum2fld1, &  ! array for global reproducible sum
      globalReprodSum2fld2, &  ! array for global reproducible sum
      globalReprodSum3fld1, &  ! array for global reproducible sum
      globalReprodSum3fld2, &  ! array for global reproducible sum
      globalReprodSum9fld1, &  ! array for global reproducible sum
      globalReprodSum9fld2     ! array for global reproducible sum
   real (kind=RKIND) :: &
      R1_alpha1s_g_dts, &! pre-computed coefficients
      R1_alpha1s_g_dt,  &! pre-computed coefficients
      R1_alpha1_g,      &! pre-computed coefficients
      alpha1,           &! implicitness parameter (=alpha)
      alpha2,           &! implicitness parameter (=1-alpha)
      tolerance_outer,  &! tolerance for the outer solver iterations
      tolerance_inner,  &! tolerance for the main solver iterations
      total_num_cells,  &! total number of surface cells
      mean_num_cells,   &! mean #cells for each core
      dt_si,            &! dt for the SI btr mode solver
      ! SI constants
      SIcst_alpha0,SIcst_alpha1,SIcst_beta0,SIcst_beta1,           &
      SIcst_omega0,SIcst_gamma0,SIcst_gamma1,SIcst_rho0,SIcst_rho1,&
      SIcst_r0rh0,SIcst_w0rh0,SIcst_r0r0 ,SIcst_r00r0,SIcst_r00w0, &
      SIcst_r00s0,SIcst_r00z0,SIcst_q0y0 ,SIcst_y0y0 ,SIcst_r00q0, &
      SIcst_r00y0,SIcst_r00t0,SIcst_r00v0,SIcst_q0q0,resid
   integer :: &
      nPrecVec,         &! vector lenghts for preconditioner
      si_ismf,          &! SI solver linearization for ISMF cases
      ncpus,            &! Total number of cores
      nSiLargeIter       ! #iteration for the barotropic system
   character(len=5) :: si_algorithm ! Name of the iterative solver algorithm

   type (field1DReal), pointer :: &
      SIvec_zh1_field,SIvec_wh0_field,SIvec_rh1_field   
                         ! field pointer for halo update

   real(kind=RKIND) :: pgf_sal_on

#ifdef USE_GPU_AWARE_MPI
   ! SI - halo exch related variables ----------------------------------

   integer :: nSendCommList,nRecvCommList
   integer :: nMaxSendList,nMaxRecvList

   integer,dimension(:),pointer :: haloLayers

   type si_exchList
      integer,allocatable,dimension(:) :: destList,srcList
   end type si_exchList

   type si_halo
      type(si_exchList),allocatable,dimension(:) :: exch
      integer,allocatable,dimension(:) :: endPointID,nList
   end type si_halo

   type si_local_mpi
      type(si_halo),allocatable,dimension(:) :: halo
      real(kind=RKIND),allocatable,dimension(:) :: buffer
      integer,allocatable,dimension(:) :: nExch,bufferOffset
      integer :: procID,nList,reqID
   end type si_local_mpi
   
   type(si_local_mpi),allocatable,dimension(:),target :: exchSend,exchRecv

   real(kind=RKIND),pointer,dimension(:) :: sbuffer,rbuffer
   real(kind=RKIND),allocatable,dimension(:) :: s1buffer,r1buffer
#endif

#ifdef USE_MAGMA
   integer(kind=8) :: dprec_ivmat,queue,dxvec,dyvec
#endif   

!***********************************************************************

   contains

!***********************************************************************
!
!  ocn_time_integration_si
!
!> \brief MPAS ocean split-implicit time integration scheme
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   September 2011 (split explicit base code)
!> \author Hyun-Gyu Kang (Oak Ridge National Laboratory)
!> \date   JAN 2019 (split-implicit code)
!> \details
!>  This routine integrates a master time step (dt) using a
!>  split-implicit time integrator.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_si(domain, dt)!{{{

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
         tracersIdealAgeFieldsPool     ! structure holding idealAge-related fields

      logical :: &
         activeTracersOnly   ! only compute tendencies for active tracers

      logical, pointer :: &
         config_use_tracerGroup  ! flag for using each tracer group

      integer :: &
         iCell, iEdge, k, &! loop iterators for cell, edge, vert loops
         nCells, nEdges,  &! number of cells or edges (incl halos)
         i,j,             &! generic loop iterators
         cell1, cell2,    &! neighbor cell addresses
         eoe,             &! index for edge on edge
         err,             &! local error flag
         splitImplicitStep, &! loop index for outer timestep loop
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
         thicknessSumCur,  &
         thicknessSumMid,  &
         thicknessSumLag,  &
         flux,             &! temp for computing flux for barotropic
         sshEdge,          &! sea surface height at edge
         sshEdgeCur,       &! 
         sshEdgeMid,       &! 
         sshEdgeLag,       &! 
         sshDiffCur,       &! sea surface height difference
         sshDiffLag,       & 
         CoriolisTerm,     &! temp for computing coriolis term (fuperp)
         normalVelocityCorrection, &! velocity correction
         temp,             &! temp for holding vars at new time
         temp_h,           &! temporary for phi
         temp_mask,        &! temporary mask (edgeMask)
         lat,              &! cell latitude
         sshCell1,         &! sea sfc height in neighboring cells
         sshCell2,         &
         sshCurArea,       &! temp for iterative solver
         sshLagArea,       &
         sshTendb1,        &
         sshTendb2,        &
         sshTendAx,        &
         fluxb1,           &
         fluxb2,           &
         fluxAx

      real (kind=RKIND), dimension(:,:), allocatable:: &
         uTemp

      real (kind=RKIND), dimension(:), allocatable :: &
         btrvel_temp, &
         bottomDepthEdge

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

      type (dm_info) :: dminfo

      ! State/Tracer arrays, info
      real (kind=RKIND), dimension(:,:,:), pointer ::      &!
         tracersGroupCur,      &! old, new tracer arrays
         tracersGroupNew

      integer :: &
         startIndex, endIndex, &! start, end index for tracer groups
         indexSalinityPtr       ! tracer index for salinity

      integer :: indexSalinity
      integer :: iTracer, numTracers, numActiveTracers

      type (mpas_pool_iterator_type) :: &
         groupItr               ! iterator for tracer groups

      character (len=StrKIND) :: &
         modifiedGroupName,    &! constructed tracer group names
         configName             ! constructed config names for tracers

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
         activeTracersTend       ! active tracer tendencies

      ! Forcing pool
      real (kind=RKIND), dimension(:), pointer :: tidalPotentialEta

      real (kind=RKIND), dimension(:), pointer :: &
        seaIcePressure, atmosphericPressure

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure, landIcePressure, landIceDraft

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracersNew

      ! Remap variables
      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThicknessLagNew

      integer ::     &
         siLargeIter  ! Index of the large btr system iteration

      ! These are public variables used from the diagnostics module
      ! indexSurfaceVelocityZonal, indexSurfaceVelocityMeridional
      ! indexSSHGradientZonal, indexSSHGradientMeridional
      ! barotropicForcing, barotropicThicknessFlux
      ! layerThickEdgeFlux, layerThickEdgeMean,
      ! normalTransportVelocity, normalGMBolusVelocity
      ! vertAleTransportTop, normalMLEvelocity
      ! velocityX, velocityY, velocityZ
      ! velocityZonal, velocityMeridional
      ! gradSSH, gradSSHX, gradSSHY, gradSSHZ
      ! gradSSHZonal, gradSSHMeridional
      ! surfaceVelocity, SSHGradient
      ! temperatureShortWaveTendency
      ! activeTracerHorizontalAdvectionTendency
      ! activeTracerVerticalAdvectionTendency
      ! activeTracerSurfaceFluxTendency
      ! activeTracerNonLocalTendency
      ! activeTracerHorMixTendency
      ! activeTracerHorizontalAdvectionEdgeFlux
      ! vertVelocityTop

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      dminfo = domain % dminfo

      call mpas_timer_start("si timestep")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !  Prep variables before first iteration
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call mpas_timer_start("si prep")

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
      call mpas_pool_get_subpool(block%structs, 'tend', &
                                                 tendPool)
      call mpas_pool_get_subpool(block%structs, 'scratch', &
                                                 scratchPool)
      call mpas_pool_get_subpool(statePool, 'tracers', &
                                             tracersPool)
      call mpas_pool_get_subpool(tendPool,  'tracersTend', &
                                             tracersTendPool)

      !*** Retrieve state variables at two time levels

      call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                           normalBaroclinicVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                           normalBarotropicVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', &
                                           normalBarotropicVelocitySubcycleCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocityCur, 1)

      call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                           normalBaroclinicVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                           normalBarotropicVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'normalBarotropicVelocitySubcycle', &
                                           normalBarotropicVelocitySubcycleNew, 2)
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocityNew, 2)

      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1)
      call mpas_pool_get_array(statePool, 'ssh', sshNew, 2)

      call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleCur, 1)
      call mpas_pool_get_array(statePool, 'sshSubcycle', sshSubcycleNew, 2)

      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThicknessNew, 2)

      call mpas_pool_get_array(statePool, 'highFreqThickness', &
                                           highFreqThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'highFreqThickness', &
                                           highFreqThicknessNew, 2)

      call mpas_pool_get_array(statePool, 'lowFreqDivergence', &
                                           lowFreqDivergenceCur, 1)
      call mpas_pool_get_array(statePool, 'lowFreqDivergence', &
                                           lowFreqDivergenceNew, 2)

      call mpas_pool_get_dimension_scalar(tracersPool, 'index_salinity', &
                                                        indexSalinity)

      !*** Retrieve tendency variables

      call mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                          highFreqThicknessTend)
      call mpas_pool_get_array(tendPool, 'normalVelocity', &
                                          normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'ssh', &
                                          sshTend)
      call mpas_pool_get_array(tendPool, 'layerThickness', &
                                          layerThicknessTend)
      call mpas_pool_get_array(tendPool, 'normalVelocity', &
                                          normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'highFreqThickness', &
                                          highFreqThicknessTend)
      call mpas_pool_get_array(tendPool, 'lowFreqDivergence', &
                                          lowFreqDivergenceTend)

      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracersNew, 2)

      allocate(bottomDepthEdge(nEdgesAll+1))

      if (config_transport_tests_flow_id > 0) then
        ! This is a transport test. Write advection velocity from prescribed
        ! flow field.
        call ocn_transport_test_velocity(meshPool, daysSinceStartOfSim, &
          0.5*dt, normalVelocityCur)
      endif

#ifdef MPAS_OPENACC
      !$acc enter data copyin(normalVelocityCur, normalBarotropicVelocityCur, &
      !$acc   sshCur, layerThicknessCur,sshSubcycleCur,sshSubcycleNew)
      !$acc enter data create(normalBaroClinicVelocityCur, normalVelocityNew, &
      !$acc   normalBaroclinicVelocityNew, sshNew, layerThicknessNew)

      if (associated(highFreqThicknessNew)) then
         !$acc enter data copyin(highFreqThicknessCur)
         !$acc enter data create(highFreqThicknessNew)
      endif
      if (associated(lowFreqDivergenceNew)) then
         !$acc enter data copyin(lowFreqDivergenceCur)
         !$acc enter data create(lowFreqDivergenceNew)
      endif

#endif

      ! Initialize * variables that are used to compute baroclinic
      ! tendencies below.

      ! The baroclinic velocity needs be recomputed at the beginning
      ! of a timestep because the implicit vertical mixing is
      ! conducted on the total u.  We keep normalBarotropicVelocity
      ! from the previous timestep.
      ! Note that normalBaroclinicVelocity may now include a
      ! barotropic component, because the weights layerThickness
      ! have changed.  That is OK, because the barotropicForcing
      ! variable subtracts out the barotropic component from the
      ! baroclinic.

#ifdef MPAS_OPENACC
      !$acc parallel present(normalVelocityCur, normalBarotropicVelocityCur, &
      !$acc    sshCur, layerThicknessCur, normalBaroClinicVelocityCur, &
      !$acc    normalVelocityNew, normalBaroclinicVelocityNew, sshNew, &
      !$acc    layerThicknessNew, minLevelCell, maxLevelCell,          &
      !$acc    sshSubcycleCur,sshSubcycleNew)
#endif

#ifdef MPAS_OPENACC
      !$acc loop collapse(2)
#else
      !$omp parallel
      !$omp do schedule(runtime) private(k)
#endif
      do iEdge = 1,nEdgesAll
      do k = 1,nVertLevels

         normalBaroclinicVelocityCur(k,iEdge) = &
                                 normalVelocityCur(k,iEdge) - &
                         normalBarotropicVelocityCur(iEdge)

         normalVelocityNew(k,iEdge) = normalVelocityCur(k,iEdge)

         normalBaroclinicVelocityNew(k,iEdge) = &
         normalBaroclinicVelocityCur(k,iEdge)
      end do
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

#ifdef MPAS_OPENACC
      !$acc loop gang
#else
      !$omp do schedule(runtime) private(k)
#endif
      do iCell = 1, nCellsAll
         sshNew(iCell) = sshCur(iCell)
         sshSubcycleCur(iCell) = sshCur(iCell)
         sshSubcycleNew(iCell) = sshCur(iCell)
#ifdef MPAS_OPENACC
         !$acc loop vector
#endif
         do k = 1,nVertLevels
            layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
         end do
      end do
#ifdef MPAS_OPENACC
      !$acc end parallel
#else
      !$omp end do
      !$omp end parallel
#endif

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr))
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr%memberName,&
                                     tracersGroupCur, 1)
            call mpas_pool_get_array(tracersPool, groupItr%memberName,&
                                     tracersGroupNew, 2)

            if ( associated(tracersGroupCur) .and. &
                 associated(tracersGroupNew) ) then

#ifdef MPAS_OPENACC
               !$acc enter data create(tracersGroupNew)
               !$acc enter data copyin(tracersGroupCur)
               !$acc parallel loop gang present(minLevelCell, maxLevelCell, &
               !$acc    tracersGroupNew, tracersGroupCur)
#else
               !$omp parallel
               !$omp do schedule(runtime) private(k)
#endif
               do iCell = 1, nCellsAll
               do k = minLevelCell(iCell), maxLevelCell(iCell)
#ifdef MPAS_OPENACC
               !$acc loop vector
#endif
                  do iTracer = 1,size(tracersGroupCur,1)
                     tracersGroupNew(iTracer,k,iCell) = &
                     tracersGroupCur(iTracer,k,iCell)
                  end do
               end do
               end do
#ifdef MPAS_OPENACC
               !$acc exit data copyout(tracersGroupNew)
               !$acc exit data delete(tracersGroupCur)
#else
               !$omp end do
               !$omp end parallel
#endif
            end if
         end if
      end do

      if (associated(highFreqThicknessNew)) then
#ifdef MPAS_OPENACC
         !$acc parallel loop gang &
         !$acc    present(highFreqThicknessCur, highFreqThicknessNew)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
         !$acc loop vector
#endif
         do k = 1,nVertLevels
            highFreqThicknessNew(k,iCell) = &
            highFreqThicknessCur(k,iCell)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      if (associated(lowFreqDivergenceNew)) then
#ifdef MPAS_OPENACC
         !$acc parallel loop gang &
         !$acc    present(lowFreqDivergenceCur, lowFreqDivergenceNew)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
         !$acc loop vector
#endif
         do k = 1,nVertLevels
            lowFreqDivergenceNew(k,iCell) = &
            lowFreqDivergenceCur(k,iCell)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      endif

      call mpas_timer_stop("si prep")

      call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)
      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)

      if (landIcePressureOn) then
         call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
         call mpas_pool_get_array(forcingPool, 'landIceDraft', landIceDraft)
      endif

#ifdef MPAS_OPENACC
      !$acc exit data delete(normalVelocityCur, normalBarotropicVelocityCur, &
      !$acc    sshCur, layerThicknessCur)
      !$acc exit data copyout(normalBaroClinicVelocityCur, normalVelocityNew, &
      !$acc    normalBaroclinicVelocityNew, sshNew, layerThicknessNew,        &
      !$acc    sshSubcycleCur,sshSubcycleNew)
      if (associated(highFreqThicknessNew)) then
         !$acc exit data delete(highFreqThicknessCur)
         !$acc exit data copyout(highFreqThicknessNew)
      endif
      if (associated(lowFreqDivergenceNew)) then
         !$acc exit data delete(lowFreqDivergenceCur)
         !$acc exit data copyout(lowFreqDivergenceNew)
      endif

      if (config_use_tidal_potential_forcing) then
         !$acc enter data copyin(pgf_sal)
         !$acc enter data create(sshSubcycleCurWithTides,tidalPotentialEta)
      endif
#endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN large outer timestep iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do splitImplicitStep = 1, numTSIterations

         if (config_disable_thick_all_tend .and. &
             config_disable_vel_all_tend .and. &
             config_disable_tr_all_tend) then
            exit ! don't compute in loop meant to update velocity,
                 ! thickness, and tracers
         end if

         call mpas_timer_start('si loop')

         stage1_tend_time = min(splitImplicitStep,2)

         call mpas_pool_get_array(statePool, 'normalVelocity', &
                           normalVelocityCur, stage1_tend_time)

#ifdef MPAS_OPENACC
         !$acc enter data copyin(layerThicknessCur, normalVelocityCur, sshCur)
         !$acc update device(layerThickEdgeFlux)
         if (config_use_freq_filtered_thickness .or. associated(highFreqThicknessNew) ) then
            !$acc enter data create(highFreqThicknessNew)
            !$acc enter data copyin(highFreqThicknessCur, highFreqThicknessTend)
         endif
#endif

         ! ---  update halos for diagnostic ocean boundary layer depth
         if (config_use_cvmix_kpp) then
            call mpas_timer_start("si halo diag obd")
            call mpas_dmpar_field_halo_exch(domain, 'boundaryLayerDepth')
            call mpas_timer_stop("si halo diag obd")
         end if

         ! ---  update halos for diagnostic variables
         call mpas_timer_start("si halo diag")

         call mpas_dmpar_field_halo_exch(domain, &
                                 'normalizedRelativeVorticityEdge')
         if (config_mom_del4 > 0.0_RKIND) then
           call mpas_dmpar_field_halo_exch(domain, 'divergence')
           call mpas_dmpar_field_halo_exch(domain, 'relativeVorticity')
         end if
         call mpas_timer_stop("si halo diag")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  Stage 1: Baroclinic velocity (3D) prediction, explicit with
         !           long timestep
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (config_use_freq_filtered_thickness) then

            call mpas_timer_start("si freq-filtered-thick computations")

            call ocn_tend_freq_filtered_thickness(tendPool, statePool, &
                                                  stage1_tend_time)

            call mpas_timer_stop("si freq-filtered-thick computations")

            call mpas_timer_start("si freq-filtered-thick halo update")

            call mpas_dmpar_field_halo_exch(domain,'tendHighFreqThickness')
            call mpas_dmpar_field_halo_exch(domain,'tendLowFreqDivergence')

            call mpas_timer_stop("si freq-filtered-thick halo update")

#ifdef MPAS_OPENACC
            !$acc parallel loop gang present(minLevelCell, maxLevelCell, &
            !$acc    highFreqThicknessNew, highFreqThicknessCur, highFreqThicknessTend)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k)
#endif
            do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
            !$acc loop vector
#endif
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               ! this is h^{hf}_{n+1}
               highFreqThicknessNew(k,iCell) = &
               highFreqThicknessCur(k,iCell) + dt* &
               highFreqThicknessTend(k,iCell)
            end do
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

         endif ! freq filtered thickness

         ! compute velocity tendencies, T(u*,w*,p*)
         call mpas_timer_start("si bcl vel")
         call mpas_timer_start('si bcl vel tend')

         ! compute vertAleTransportTop.  Use u (rather than &
         ! normalTransportVelocity) for momentum advection.
         ! Use the most recent time level available.

         if (associated(highFreqThicknessNew)) then
            call ocn_vert_transport_velocity_top( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalVelocityCur, sshCur, dt, &
                 vertAleTransportTop, err, highFreqThicknessNew)
         else
            call ocn_vert_transport_velocity_top( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalVelocityCur, sshCur, dt, &
                 vertAleTransportTop, err)
         endif

#ifdef MPAS_OPENACC
         !$acc exit data delete(layerThicknessCur, normalVelocityCur, sshCur)
         !$acc update host(vertAleTransportTop)
         if (config_use_freq_filtered_thickness .or. associated(highFreqThicknessNew) ) then
            !$acc exit data copyout(highFreqThicknessNew)
            !$acc exit data delete(highFreqThicknessCur, highFreqThicknessTend)
         endif
#endif

         call ocn_tend_vel(domain, tendPool, statePool, forcingPool, &
                           stage1_tend_time, domain % dminfo, dt)

         call mpas_timer_stop('si bcl vel tend')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do j=1,numClinicIterations(splitImplicitStep)

            call mpas_timer_start('bcl iters on linear Coriolis')

            ! Put f*normalBaroclinicVelocity^{perp} in
            ! normalVelocityNew as a work variable
            call ocn_fuperp(statePool, meshPool, 2)

            !TODO: look at loop optimizations here
            ! Only need to loop over owned cells, since there is a halo
            ! exchange immediately after this computation.

            allocate(uTemp(nVertLevels,nEdgesOwned))

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(k, cell1, cell2, uTemp, &
            !$omp         normalThicknessFluxSum, thicknessSum)
            do iEdge = 1, nEdgesOwned
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               ! could put this after with uTemp(maxleveledgetop+1:nvertlevels)=0
               uTemp(:,iEdge) = 0.0_RKIND
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  ! normalBaroclinicVelocityNew =
                  ! normalBaroclinicVelocityOld +
                  !                dt*(-f*normalBaroclinicVelocityPerp
                  !                    + T(u*,w*,p*) + g*grad(SSH*) )
                  ! Here uNew is a work variable containing
                  !  -fEdge(iEdge)*normalBaroclinicVelocityPerp(k,iEdge)
                  uTemp(k,iEdge) = normalBaroclinicVelocityCur(k,iEdge) &
                           + dt * (normalVelocityTend(k,iEdge) &
                           + normalVelocityNew(k,iEdge) &
                           + gravity * &
                             (sshNew(cell2) - sshNew(cell1)) &
                             /dcEdge(iEdge) )
               enddo ! vertical

               ! thicknessSum is initialized outside the loop because
               ! on land boundaries maxLevelEdgeTop=0, but we want to
               ! initialize thicknessSum with a nonzero value to avoid
               ! a NaN.
               normalThicknessFluxSum =                          &
                    layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge) &
                  * uTemp(minLevelEdgeBot(iEdge),iEdge)
               thicknessSum = &
                    layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)

               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                  normalThicknessFluxSum = normalThicknessFluxSum + &
                                 layerThickEdgeFlux(k,iEdge)*uTemp(k,iEdge)
                  thicknessSum = thicknessSum + &
                                 layerThickEdgeFlux(k,iEdge)
               enddo
               barotropicForcing(iEdge) = &
                              normalThicknessFluxSum/thicknessSum/dt

               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  ! These two steps are together here:
                  !{\bf u}'_{k,n+1} =
                  !    {\bf u}'_{k,n} - \Delta t {\overline {\bf G}}
                  !{\bf u}'_{k,n+1/2} = \frac{1}{2}\left(
                  !        {\bf u}^{'}_{k,n} +{\bf u}'_{k,n+1}\right)
                  ! so that normalBaroclinicVelocityNew is at time n+1/2
                  normalBaroclinicVelocityNew(k,iEdge) = 0.5_RKIND*( &
                  normalBaroclinicVelocityCur(k,iEdge) + uTemp(k,iEdge) - &
                         dt * barotropicForcing(iEdge))
               enddo

            enddo ! iEdge
            !$omp end do
            !$omp end parallel

            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(cell1, cell2, k, thicknessSum)
            do iEdge = 1, nEdgesHalo(config_num_halos+1)
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               thicknessSum = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)
               do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                  thicknessSum = thicknessSum + &
                                 layerThickEdgeFlux(k,iEdge)
               enddo
               bottomDepthEdge(iEdge) = thicknessSum &
                  - 0.5_RKIND*(sshNew(cell1) + sshNew(cell2))
            enddo ! iEdge
            !$omp end do
            !$omp end parallel

            deallocate(uTemp)

            call mpas_timer_start("si halo normalBaroclinicVelocity")
            call mpas_dmpar_field_halo_exch(domain, &
                              'normalBaroclinicVelocity', timeLevel=2)
            call mpas_timer_stop("si halo normalBaroclinicVelocity")

            call mpas_timer_stop('bcl iters on linear Coriolis')

         end do  ! do j=1,config_n_bcl_iter

         call mpas_timer_start('si halo barotropicForcing')
         call mpas_dmpar_field_halo_exch(domain, 'barotropicForcing')
         call mpas_timer_stop('si halo barotropicForcing')

         call mpas_timer_stop("si bcl vel")

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! END baroclinic iterations on linear Coriolis term
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Stage 2: Barotropic velocity (2D) prediction, split-implicit
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! The split-implicit barotropic mode solver
         !    - uses the preconditioned communication-avoiding
         !      (single-global synchronization) BiCGStab method 
         !      as a linear iterative solver
         !      (Kang et al., in preparation)
         !                    
         !    - uses the s-step conjugate gradient method (Chronopoulos
         !      and Gear, 1989) when 'config_btr_si_partition_match_mode'
         !      is turned on
         !                    
         !    - solves the nonlinear barotropic system but deals with
         !      it as the linear system by introducing the outer and                     
         !      inner solver iteration (the solver matrix is updated                    
         !      for every time step and inner solver iteration.)
         !                    
         !    - works with 'ras', 'block_jacobi', 'jacobi', 'none'
         !      preconditioners
         !                    
         !    - consists of seven substages
         !         :Stage 2.1. Initialization and preparation of vars
         !                     before outer iterations
         !          Stage 2.2. Computation of the initial residual
         !          Stage 2.3. The outer solver iteration
         !          Stage 2.4. The barotropic velocity update
         !          Stage 2.5. Re-computation of the initial residual
         !                     using lagged values
         !          Stage 2.6. Main (inner) solver iterations to 
         !                     obtain SSH at time (n+1) 
         !          Stage 2.7. The barotropic velocity update
         !
         !    - recommends to set 'config_n_ts_iter=2'
         !         :Time line of variables
         !
         !          Start of large outer time step iteration loop
         !             splitImplicitStep=1
         !             |  SSH and BtrVel at time (n)
         !             |  |  -> Advancing the barotropic system
         !             |  |     -> SSH and BtrVel at time (n+1)
         !             |  |        -> Averaging between time (n) 
         !             |  |                         and time (n+1)
         !             |  SSH and BtrVel at time (n+1/2)
         !             |
         !             splitImplicitStep=2
         !                SSH and BtrVel at time (n+1/2)
         !                |  -> Advancing the barotropic system
         !                |     -> SSH and BtrVel at time (n+3/2)
         !                |        -> Averaging between time (n+1/2)
         !                |                         and time (n+3/2)
         !                SSH and BtrVel at time (n+1)
         !
         !
         ! Reference: Kang et al. (2021): A scalable semi-implicit 
         !            barotropic mode solver for the MPAS-Ocean, JAMES
         !
         !
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Stage 2.1 : Preparation of variables before outer iterations
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPAS_OPENACC
         !$acc enter data copyin(sshCur,sshNew,                       &
         !$acc                   sshSubcycleCur,sshSubcycleNew,       &
         !$acc                   normalBarotropicVelocitySubcycleCur, &
         !$acc                   normalBarotropicVelocitySubcycleNew, &
         !$acc                   normalBarotropicVelocityCur,         &
         !$acc                   normalBarotropicVelocityNew,         &
         !$acc                   bottomDepthEdge)
         !$acc update device(barotropicForcing)
#endif

         call mpas_timer_start("si btr vel")

         cellHaloComputeCounter = config_num_halos
         edgeHaloComputeCounter = config_num_halos + 1

         nCells = nCellsHalo(cellHaloComputeCounter-1)
         nEdges = nEdgesHalo(edgeHaloComputeCounter-1)

         if (config_filter_btr_mode) then
#ifdef MPAS_OPENACC
            !$acc parallel loop present(barotropicForcing)
#else
            !$omp parallel
            !$omp do schedule(runtime)
#endif
            do iEdge = 1, nEdges
               barotropicForcing(iEdge) = 0.0_RKIND
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
         endif


         !-------------------------------------------------------------!
         ! BEGIN Large barotropic system iteration loop
         !-------------------------------------------------------------!
         do siLargeIter = 1, nSiLargeIter
         !-------------------------------------------------------------!

            ! Initialize variables for barotropic subcycling
            call mpas_timer_start('btr vel si init')
   
            if ( splitImplicitStep == 2 .or. siLargeIter > 1 ) then
               ! Here sshSubcycleNew and normalBarotropicVelocityCur
               ! are at time n+1/2.
   
#ifdef MPAS_OPENACC
               !$acc parallel present(sshSubcycleNew,sshCur,        &
               !$acc    sshSubcycleCur,normalBarotropicVelocityCur, &
               !$acc    normalBarotropicVelocitySubcycleCur)
               !$acc loop gang
#else
               !$omp parallel
               !$omp do schedule(runtime)
#endif
               do iCell = 1,nCellsAll
                  sshSubcycleNew(iCell) = &
                     0.5_RKIND*( sshSubcycleNew(iCell) &
                                +sshCur(iCell) )
                  sshSubcycleCur(iCell) = sshSubcycleNew(iCell)
               end do
#ifndef MPAS_OPENACC
               !$omp end do
#endif
   

#ifdef MPAS_OPENACC
               !$acc loop gang
#else
               !$omp do schedule(runtime)
#endif
               do iEdge = 1,nEdgesAll
                  normalBarotropicVelocityCur(iEdge) = &
                     normalBarotropicVelocitySubcycleCur(iEdge)
               end do
#ifdef MPAS_OPENACC
               !$acc end parallel
#else
               !$omp end do
               !$omp end parallel
#endif
            endif ! splitImplicitStep
   
#ifdef MPAS_OPENACC
            !$acc parallel loop gang &
            !$acc    present(nEdgesOnEdge,weightsOnEdge,        &
            !$acc            normalBarotropicVelocityCur,fEdge, &
            !$acc            barotropicCoriolisTerm,            &
            !$acc            edgesOnEdge)                       &
            !$acc    private(CoriolisTerm,i,eoe)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(CoriolisTerm, i, eoe)
#endif
            do iEdge = 1, nEdgesAll
               ! Compute the barotropic Coriolis term, -f*uPerp
               CoriolisTerm = 0.d0
               !$acc loop vector reduction(+:CoriolisTerm)
               do i = 1, nEdgesOnEdge(iEdge)
                  eoe = edgesOnEdge(i,iEdge)
                  CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge) &
                               * normalBarotropicVelocityCur(eoe)      &
                               * fEdge(eoe)
               end do ! i
                  barotropicCoriolisTerm(iEdge) = CoriolisTerm
            end do ! iEdge
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
   
            ! Subtract tidal potential from ssh, if needed
            ! Subtract the tidal potential from the current
            !    subcycle ssh and store and a work array.
            ! Then point sshSubcycleCur to the work array so the
            ! tidal potential terms are included in the grad
            ! operator inside the edge loop.
            if (config_use_tidal_potential_forcing) then
   
               call mpas_pool_get_array(forcingPool,            &
                                     'sshSubcycleCurWithTides', &
                                      sshSubcycleCurWithTides)
               call mpas_pool_get_array(forcingPool,        &
                                       'tidalPotentialEta', &
                                        tidalPotentialEta)
   
#ifdef MPAS_OPENACC
               !$acc update device(sshSubcycleCurWithTides,tidalPotentialEta)
 
               !$acc parallel loop gang  &
               !$acc    present(sshSubcycleCurWithTides,sshSubcycleCur, &
               !$acc       tidalPotentialEta,pgf_sal,sshSubcycleCur)
#else
               !$omp parallel
               !$omp do schedule(runtime)
#endif
               do iCell = 1, nCellsAll
                  sshSubcycleCurWithTides(iCell) = &
                           sshSubcycleCur(iCell) - &
                        tidalPotentialEta(iCell) - &
                     pgf_sal_on * pgf_sal(iCell) - &
                     (1.0_RKIND - pgf_sal_on) * self_attraction_and_loading_beta* &
                           sshSubcycleCur(iCell)
               end do
#ifndef MPAS_OPENACC
               !$omp end do
               !$omp end parallel
#endif
   
               call mpas_pool_get_array(forcingPool,            &
                                     'sshSubcycleCurWithTides', &
                                      sshSubcycleNew)
   
#ifdef MPAS_OPENACC
               !$acc update device(sshSubcycleNew)

               !$acc parallel loop gang  &
               !$acc    present(sshSubcycleCur,sshSubcycleNew)
#else
               !$omp parallel
               !$omp do schedule(runtime)
#endif
               do iCell = 1,nCellsAll
                  sshSubcycleCur(iCell) = sshSubcycleNew(iCell)
               end do
#ifndef MPAS_OPENACC
               !$omp end do
               !$omp end parallel
#endif

            endif !config_use_tidal_potential_forcing
   
            call mpas_timer_stop('btr vel si init')
   
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Stage 2.2 : Compution of the initial residual
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            call mpas_timer_start("si btr residual")
   
#ifdef MPAS_OPENACC
            !$acc parallel loop gang async                              &
            !$acc    present(nEdgesOnCell,edgesOnCell,cellsOnEdge,      &
            !$acc       sshSubcycleCur,bottomDepthEdge,dcEdge,          &
            !$acc       normalBarotropicVelocityCur,maxLevelEdgeTop,    &
            !$acc       barotropicCoriolisTerm,barotropicForcing,       &
            !$acc       edgeSignOnCell,dvEdge,areaCell,SIvec_r0,        &
            !$acc       SIvec_r00,SIvec_s0,SIvec_z0)                    &
            !$acc    private(sshTendb1,sshTendb2,sshTendAx,iEdge,cell1, &
            !$acc       cell2,sshEdgeCur,thicknessSumCur,sshDiffCur,    &
            !$acc       fluxb1,fluxb2,fluxAx,sshCurArea,i)
#else
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(sshTendb1,sshTendb2,sshTendAx,iEdge, &
            !$omp         cell1,cell2,sshEdgeCur,thicknessSumCur, &
            !$omp         sshDiffCur,fluxb1,fluxb2,fluxAx,sshCurArea)
#endif
            do iCell = 1, nPrecVec
               sshTendb1 = 0.0_RKIND
               sshTendb2 = 0.0_RKIND
               sshTendAx = 0.0_RKIND
   
               !$acc loop vector                                  &
               !$acc    reduction(+:sshTendb1,sshTendb2,sshTendAx)
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
                  if (maxLevelEdgeTop(iEdge).eq.0) cycle
   
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
   
                  ! Interpolation sshEdge
                  sshEdgeCur = 0.5_RKIND * ( sshSubcycleCur(cell1)  &
                                           + sshSubcycleCur(cell2) )
   
                  thicknessSumCur = si_ismf * sshEdgeCur    &
                                  + bottomDepthEdge(iEdge)
   
                  ! nabla (ssh^0)
                  sshDiffCur = (  sshSubcycleCur(cell2)   &
                                - sshSubcycleCur(cell1) ) &
                               /  dcEdge(iEdge)
   
                  fluxb1 = thicknessSumCur &
                         * normalBarotropicVelocityCur(iEdge)
                  fluxb2 = thicknessSumCur                    &
                         * (alpha2*gravity*sshDiffCur         &
                         + (-barotropicCoriolisTerm(iEdge)    &
                            -barotropicForcing(iEdge)))
                  fluxAx = thicknessSumCur * sshDiffCur
   
                  sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) &
                                        * fluxb1 * dvEdge(iEdge)
                  sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) &
                                        * fluxb2 * dvEdge(iEdge)
                  sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) &
                                        * fluxAx * dvEdge(iEdge)
               end do ! i
   
               sshTendb1  = R1_alpha1s_g_dt  * sshTendb1
               sshTendb2  = R1_alpha1_g      * sshTendb2
               sshCurArea = R1_alpha1s_g_dts * sshSubcycleCur(iCell) &
                                             * areaCell(iCell)
   
               SIvec_r0(iCell) = (-sshCurArea - sshTendb1 + sshTendb2) &
                                -(-sshCurArea - sshTendAx)
               SIvec_r00(iCell) = SIvec_r0(iCell)
               SIvec_s0(iCell)  = 0.0_RKIND
               SIvec_z0(iCell)  = 0.0_RKIND
            end do ! iCell
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
   
            ! Preconditioning --------------------------------------------!

            call si_precond(SIvec_r0,SIvec_rh0)
            !$acc wait
   
            call mpas_timer_start("si 1st halo")
#ifdef USE_GPU_AWARE_MPI
            call si_halo_exch(domain,SIvec_rh0)
#else
            call mpas_dmpar_field_halo_exch(domain, 'SIvec_rh0')
            !$acc update device(SIvec_rh0) async
#endif
            call mpas_timer_stop("si 1st halo")
   
            ! SpMV -------------------------------------------------------!
   
            call si_matvec_mul(SIvec_rh0,bottomDepthEdge, &
                               SIvec_w0,sshSubcycleCur)

            if ( si_algorithm == 'sbicg' ) then !!!

            ! Preconditioning --------------------------------------------!

               call si_precond(SIvec_w0,SIvec_wh0)
               !$acc wait

               call mpas_timer_start("si 1st halo")
#ifdef USE_GPU_AWARE_MPI
               call si_halo_exch(domain,SIvec_wh0)
#else
               call mpas_dmpar_field_halo_exch(domain, 'SIvec_wh0')
               !$acc update device(SIvec_wh0) async
#endif
               call mpas_timer_stop("si 1st halo")
      
               ! SpMV -------------------------------------------------------!
      
               call si_matvec_mul(SIvec_wh0,bottomDepthEdge, &
                                  SIvec_t0,sshSubcycleCur)
   
               ! Reduction --------------------------------------------------!
   
               call mpas_timer_start("si 1st global reduction")
               call si_global_reduction(SIcst_allreduce_local2, &
                                   SIcst_allreduce_global2,     &
                                   globalReprodSum2fld1,        &
                                   globalReprodSum2fld2,        &
                                   si_algorithm,2,domain,       &
                                   SIvec_r00,SIvec_r0,SIvec_w0)
               call mpas_timer_stop("si 1st global reduction")
   
               !$acc serial async present(SIcst_r00r0,SIcst_r00w0,    &
               !$acc           SIcst_alpha0,SIcst_beta0,SIcst_omega0, &
               !$acc           SIcst_alpha1,SIcst_beta1,SIcst_rho0,   &
               !$acc           SIcst_allreduce_global2)
               SIcst_r00r0  = SIcst_allreduce_global2(1)
               SIcst_r00w0  = SIcst_allreduce_global2(2)
               SIcst_rho0   = SIcst_r00r0
               SIcst_alpha0 = SIcst_rho0 / SIcst_r00w0
               SIcst_beta0  = 0.0_RKIND
               SIcst_omega0 = 0.0_RKIND

               SIcst_alpha1 = SIcst_alpha0
               SIcst_beta1  = SIcst_beta0
               SIcst_r00r0  = 0.0_RKIND
               SIcst_r00w0  = 0.0_RKIND
               !$acc end serial
               !$acc wait

            else if ( si_algorithm == 'scg' ) then !!!
           
               ! Reduction --------------------------------------------------!

               call mpas_timer_start("si 1st global reduction")
               call si_global_reduction(SIcst_allreduce_local2, &
                                   SIcst_allreduce_global2,     &
                                   globalReprodSum2fld1,        &
                                   globalReprodSum2fld2,        &
                                   si_algorithm,2,domain,       &
                                   SIvec_r0,SIvec_w0,SIvec_rh0)
               call mpas_timer_stop("si 1st global reduction")
   
               !$acc serial present(SIcst_r0rh0,SIcst_w0rh0,SIcst_alpha0, &
               !$acc           SIcst_beta0,SIcst_gamma0,                  &
               !$acc           SIcst_allreduce_global2)
               SIcst_r0rh0  = SIcst_allreduce_global2(1)
               SIcst_w0rh0  = SIcst_allreduce_global2(2)
               SIcst_alpha0 = SIcst_r0rh0 / SIcst_w0rh0
               SIcst_beta0  = 0.0_RKIND
               SIcst_gamma0 = SIcst_r0rh0
               SIcst_r0rh0  = 0.0_RKIND
               SIcst_w0rh0  = 0.0_RKIND
               !$acc end serial

            end if !!!
   
            call mpas_timer_stop("si btr residual")
   
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Stage 2.3 : Outer iterations 
            !             - lagged values are sufficiently up to date
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            call mpas_timer_start("si btr iteration")

            if ( si_algorithm == 'sbicg' ) then
               call si_solver_sbicg(domain,                                &
                    sshSubcycleNew,sshSubcycleCur,bottomDepthEdge,         &
                    tolerance_outer)
            else if ( si_algorithm == 'scg' ) then
               call si_solver_scg(domain,                                  &
                    sshSubcycleNew,sshSubcycleCur,bottomDepthEdge,         &
                    tolerance_outer)
            endif

            !   boundary update on sshNew
            call mpas_timer_start("si halo iter")
            !$acc update host(sshSubcycleNew)
            call mpas_dmpar_field_halo_exch(domain, 'sshSubcycle', &
                                            timeLevel=2)
            !$acc update device(sshSubcycleNew) async
            call mpas_timer_stop("si halo iter")
   
            call mpas_timer_stop("si btr iteration")
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Stage 2.4 : The barotropic velocity update
            !             using the lagged SSH
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            call mpas_timer_start("si btr vel update")

#ifdef MPAS_OPENACC
            !$acc parallel loop gang async &
            !$acc    present(sshNew,sshSubcycleNew)
#else
            !$omp parallel
            !$omp do schedule(runtime)
#endif
            do iCell = 1,nCellsAll
               ! Use of sshNew to save the lagged SSH
               sshNew(iCell) = sshSubcycleNew(iCell)
            end do
            !$omp end do
            !$omp end parallel
   
#ifdef MPAS_OPENACC
            !$acc parallel loop async &
            !$acc    present(edgeMask,cellsOnEdge,sshNew,dcEdge,    &
            !$acc       normalBarotropicVelocityNew,sshSubcycleCur, &
            !$acc       normalBarotropicVelocitySubcycleNew,        &
            !$acc       normalBarotropicVelocityCur,                &
            !$acc       barotropicCoriolisTerm,barotropicForcing,   &
            !$acc       nEdgesHalo)                                 &
            !$acc    private(temp_mask,cell1,cell2)
#else
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(temp_mask,cell1,cell2)
#endif
            do iEdge = 1,nEdgesHalo( edgeHaloComputeCounter-1 )
               temp_mask = edgeMask(1, iEdge)
   
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
                      
               ! Use of normBtrVelNew as a temp array for sshEdge
               normalBarotropicVelocityNew(iEdge) =       &
                          ( 0.5*sshNew(cell2)          &
                           +0.5*sshSubcycleCur(cell2)) &
                         -( 0.5*sshNew(cell1)          &
                           +0.5*sshSubcycleCur(cell1)) 
    
               normalBarotropicVelocitySubcycleNew(iEdge)              &
                  = temp_mask                                          &
                  * (normalBarotropicVelocityCur(iEdge)                &
                  - dt_si * (-barotropicCoriolisTerm(iEdge) + gravity     &
                          *normalBarotropicVelocityNew(iEdge)          &
                          /dcEdge(iEdge) - barotropicForcing(iEdge)))
            end do ! iEdge
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
   
#ifdef MPAS_OPENACC
            !$acc parallel loop gang async &
            !$acc    present(nEdgesOnEdge,edgesOnEdge,weightsOnEdge, &
            !$acc        normalBarotropicVelocitySubcycleNew,        &
            !$acc        normalBarotropicVelocityCur,fEdge,          &
            !$acc        barotropicCoriolisTerm,nEdgesHalo)          &
            !$acc    private(CoriolisTerm,eoe,i)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(CoriolisTerm, i, eoe)
#endif
            do iEdge = 1,nEdgesHalo(2)
               ! Compute the barotropic Coriolis term, -f*uPerp
               CoriolisTerm = 0.0_RKIND
               !$acc loop vector reduction(+:CoriolisTerm)
               do i = 1, nEdgesOnEdge(iEdge)
                  eoe = edgesOnEdge(i,iEdge)
                  CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge)       &
                               * 0.5_RKIND                                   &
                               * ( normalBarotropicVelocitySubcycleNew(eoe)  &
                                  +normalBarotropicVelocityCur(eoe)        ) &
                               * fEdge(eoe)
               end do
               barotropicCoriolisTerm(iEdge) = CoriolisTerm
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

#ifdef MPAS_OPENACC
            !$acc parallel loop gang async &
            !$acc    present(edgeMask,cellsOnEdge,dcEdge,         &
            !$acc       normalBarotropicVelocityNew,              &
            !$acc       normalBarotropicVelocitySubcycleNew,      &
            !$acc       normalBarotropicVelocityCur,              &
            !$acc       barotropicCoriolisTerm,barotropicForcing, &
            !$acc       nEdgesHalo) &
            !$acc    private(temp_mask,cell1,cell2)
#else
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(temp_mask,cell1,cell2)
#endif
            do iEdge = 1,nEdgesHalo(2)
               temp_mask = edgeMask(1, iEdge)
   
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
               normalBarotropicVelocitySubcycleNew(iEdge)            &
                  = temp_mask                                        &
                  * (normalBarotropicVelocityCur(iEdge)              &
                  - dt_si * (-barotropicCoriolisTerm(iEdge) + gravity   &
                          *normalBarotropicVelocityNew(iEdge)        &
                          /dcEdge(iEdge) - barotropicForcing(iEdge)))
            end do ! iEdge
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
   

#ifdef MPAS_OPENACC
            !$acc parallel loop gang async &
            !$acc    present(nEdgesOnEdge,edgesOnEdge,weightsOnEdge, &
            !$acc        normalBarotropicVelocitySubcycleNew,        &
            !$acc        normalBarotropicVelocityCur,fEdge,          &
            !$acc        barotropicCoriolisTerm)                     &
            !$acc    private(CoriolisTerm,eoe,i)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(CoriolisTerm, i, eoe)
#endif
            do iEdge = 1,nEdgesOwned
               ! Compute the barotropic Coriolis term, -f*uPerp
               CoriolisTerm = 0.0_RKIND
               !$acc loop vector reduction(+:CoriolisTerm)
               do i = 1, nEdgesOnEdge(iEdge)
                  eoe = edgesOnEdge(i,iEdge)
                  CoriolisTerm = CoriolisTerm + weightsOnEdge(i,iEdge)      &
                               * 0.5_RKIND                                  &
                               * ( normalBarotropicVelocitySubcycleNew(eoe) &
                                  +normalBarotropicVelocityCur(eoe) )       &
                               * fEdge(eoe)
               end do
               barotropicCoriolisTerm(iEdge) = CoriolisTerm
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
   
            call mpas_timer_start("si halo btr coriolis")
            !$acc wait
            !$acc update host(barotropicCoriolisTerm)
            call mpas_dmpar_field_halo_exch(domain, 'barotropicCoriolisTerm')
            !$acc update device(barotropicCoriolisTerm) async
            call mpas_timer_stop("si halo btr coriolis")
   
            call mpas_timer_stop ("si btr vel update")
   
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Stage 2.5 : Recompution of the initial residual
            !             with lagged values
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            call mpas_timer_start("si btr residual")
   
            ! SpMV -------------------------------------------------------!
   
#ifdef MPAS_OPENACC
            !$acc parallel loop gang async &
            !$acc    present(nEdgesOnCell,edgesOnCell,cellsOnEdge,        &
            !$acc       sshSubcycleCur,sshNew,bottomDepthEdge,dcEdge,     &
            !$acc       normalBarotropicVelocityCur,maxLevelEdgeTop,      &
            !$acc       barotropicCoriolisTerm,barotropicForcing,         &
            !$acc       edgeSignOnCell,dvEdge,areaCell,SIvec_r0,          &
            !$acc       SIvec_r00,SIvec_s0,SIvec_z0)                      &
            !$acc    private(sshTendb1,sshTendb2,sshTendAx,iEdge,cell1,   &
            !$acc       cell2,sshEdgeCur,sshEdgeLag,sshEdgeMid,k,         &
            !$acc       thicknessSumCur,thicknessSumLag,thicknessSumMid,  &
            !$acc       sshDiffCur,sshDiffLag,                            &
            !$acc       fluxb1,fluxb2,fluxAx,sshCurArea,sshLagArea)
#else
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(sshTendb1,sshTendb2,sshTendAx,iEdge, &
            !$omp         cell1,cell2,sshEdgeCur,sshEdgeLag,sshEdgeMid, &
            !$omp         thicknessSumCur,thicknessSumMid, &
            !$omp         thicknessSumLag,sshDiffCur,sshDiffLag, &
            !$omp         fluxb1,fluxb2,fluxAx,sshCurArea,sshLagArea)
#endif
            do iCell = 1, nPrecVec
               sshTendb1 = 0.0_RKIND
               sshTendb2 = 0.0_RKIND
               sshTendAx = 0.0_RKIND
   
               !$acc loop vector reduction(+:sshTendb1,sshTendb2,sshTendAx)
               do i = 1, nEdgesOnCell(iCell)
                  iEdge = edgesOnCell(i, iCell)
                  if (maxLevelEdgeTop(iEdge).eq.0) cycle
   
                  cell1 = cellsOnEdge(1, iEdge)
                  cell2 = cellsOnEdge(2, iEdge)
   
                  ! Interpolation sshEdge
                  sshEdgeCur = 0.5_RKIND * (  sshSubcycleCur(cell1)   &
                                            + sshSubcycleCur(cell2) )
                  sshEdgeLag = 0.5_RKIND * (sshNew(cell1) + sshNew(cell2))
                  sshEdgeMid = alpha1 * sshEdgeLag + alpha2 * sshEdgeCur
   
                  ! method 1, matches method 0 without pbcs,
                  ! works with pbcs.
                  thicknessSumCur = si_ismf * sshEdgeCur    &
                                  + bottomDepthEdge(iEdge)
                  thicknessSumLag = si_ismf * sshEdgeLag    &
                                  + bottomDepthEdge(iEdge)
                  thicknessSumMid = si_ismf * sshEdgeMid    &
                                  + bottomDepthEdge(iEdge)
   
                  ! nabla (ssh^0)
                  sshDiffCur = ( sshSubcycleCur(cell2)   &
                                -sshSubcycleCur(cell1) ) &
                               / dcEdge(iEdge)
                  sshDiffLag = (sshNew(cell2)-sshNew(cell1)) &
                               / dcEdge(iEdge)
   
                  fluxb1 = thicknessSumMid &
                         * normalBarotropicVelocityCur(iEdge)
                  fluxb2 = thicknessSumLag &
                         * ( alpha2*gravity*sshDiffCur &
                            + (-barotropicCoriolisTerm(iEdge) &
                               -barotropicForcing(iEdge)) )
                  fluxAx = thicknessSumLag * sshDiffLag
   
                  sshTendb1 = sshTendb1 + edgeSignOnCell(i, iCell) &
                                        * fluxb1 * dvEdge(iEdge)
                  sshTendb2 = sshTendb2 + edgeSignOnCell(i, iCell) &
                                        * fluxb2 * dvEdge(iEdge)
                  sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) &
                                        * fluxAx * dvEdge(iEdge)
               end do ! i
   
               sshTendb1  = R1_alpha1s_g_dt  * sshTendb1
               sshTendb2  = R1_alpha1_g      * sshTendb2
               sshCurArea = R1_alpha1s_g_dts * sshSubcycleCur(iCell) &
                                             * areaCell(iCell)
               sshLagArea = R1_alpha1s_g_dts * sshNew(iCell) &
                                             * areaCell(iCell)
   
               SIvec_r0(iCell) = (-sshCurArea - sshTendb1 + sshTendb2) &
                                -(-sshLagArea - sshTendAx)
               SIvec_r00(iCell) = SIvec_r0(iCell)
               SIvec_s0(iCell)  = 0.0_RKIND
               SIvec_z0(iCell)  = 0.0_RKIND
            end do ! iCell
            !$omp end do
            !$omp end parallel
   
            ! Preconditioning --------------------------------------------!

            call si_precond(SIvec_r0,SIvec_rh0)
            !$acc wait
   
            call mpas_timer_start("si halo r0")
#ifdef USE_GPU_AWARE_MPI
            call si_halo_exch(domain,SIvec_rh0)
#else
            call mpas_dmpar_field_halo_exch(domain, 'SIvec_rh0')
            !$acc update device(SIvec_rh0) async
#endif
            call mpas_timer_stop("si halo r0")
   
            ! SpMV -------------------------------------------------------!

            call si_matvec_mul(SIvec_rh0,bottomDepthEdge,SIvec_w0,sshNew)

            if ( si_algorithm == 'sbicg' ) then !!!

            ! Preconditioning --------------------------------------------!

               call si_precond(SIvec_w0,SIvec_wh0)
               !$acc wait
      
               call mpas_timer_start("si halo r0")
#ifdef USE_GPU_AWARE_MPI
               call si_halo_exch(domain,SIvec_wh0)
#else
               call mpas_dmpar_field_halo_exch(domain, 'SIvec_wh0')
               !$acc update device(SIvec_wh0) async
#endif
               call mpas_timer_stop("si halo r0")
      
               ! SpMV -------------------------------------------------------!
      
               call si_matvec_mul(SIvec_wh0,bottomDepthEdge,SIvec_t0,sshNew)
   
               ! Reduction --------------------------------------------------!
   
               call si_global_reduction(SIcst_allreduce_local2, &
                                   SIcst_allreduce_global2,     &
                                   globalReprodSum2fld1,        &
                                   globalReprodSum2fld2,        &
                                   si_algorithm,2,domain,       &
                                   SIvec_r00,SIvec_r0,SIvec_w0)
   
               !$acc serial async present(SIcst_r00r0,SIcst_r00w0,    &
               !$acc           SIcst_alpha0,SIcst_beta0,SIcst_omega0, &
               !$acc           SIcst_alpha1,SIcst_beta1,SIcst_rho0,   &
               !$acc           SIcst_allreduce_global2)
               SIcst_r00r0 = SIcst_allreduce_global2(1)
               SIcst_r00w0 = SIcst_allreduce_global2(2)
               SIcst_rho0   = SIcst_r00r0
               SIcst_alpha0 = SIcst_rho0 / SIcst_r00w0
               SIcst_beta0  = 0.0_RKIND
               SIcst_omega0 = 0.0_RKIND

               SIcst_alpha1 = SIcst_alpha0
               SIcst_beta1  = SIcst_beta0
               SIcst_r00r0  = 0.0_RKIND
               SIcst_r00w0  = 0.0_RKIND
               !$acc end serial

            else if ( si_algorithm == 'scg' ) then !!!
           
               ! Reduction --------------------------------------------------!
               call si_global_reduction(SIcst_allreduce_local2, &
                                   SIcst_allreduce_global2,     &
                                   globalReprodSum2fld1,        &
                                   globalReprodSum2fld2,        &
                                   si_algorithm,2,domain,       &
                                   SIvec_r0,SIvec_w0,SIvec_rh0)
   
               !$acc serial present(SIcst_r0rh0,SIcst_w0rh0,SIcst_alpha0, &
               !$acc           SIcst_beta0,SIcst_gamma0,                  &
               !$acc           SIcst_allreduce_global2)
               SIcst_r0rh0 = SIcst_allreduce_global2(1)
               SIcst_w0rh0 = SIcst_allreduce_global2(2)
               SIcst_alpha0 = SIcst_r0rh0 / SIcst_w0rh0
               SIcst_beta0  = 0.0_RKIND
               SIcst_gamma0 = SIcst_r0rh0
               SIcst_r0rh0  = 0.0_RKIND
               SIcst_w0rh0  = 0.0_RKIND
               !$acc end serial
   
            end if !!!

            call mpas_timer_stop("si btr residual")
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Stage 2.6 : Main (inner) solver iterations
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            call mpas_timer_start("si btr iteration")

            if ( si_algorithm == 'sbicg' ) then
               call si_solver_sbicg(domain,                                &
                    sshSubcycleNew,sshNew,bottomDepthEdge,tolerance_inner)
            else if ( si_algorithm == 'scg' ) then
               call si_solver_scg(domain,                                  &
                    sshSubcycleNew,sshNew,bottomDepthEdge,tolerance_inner)
            endif

            ! boundary update on SSHnew
            call mpas_timer_start("si halo iter")
            !$acc update host(sshSubcycleNew)
            call mpas_dmpar_field_halo_exch(domain, 'sshSubcycle', timeLevel=2)
            !$acc update device(sshSubcycleNew)
            call mpas_timer_stop("si halo iter")
   
            call mpas_timer_stop("si btr iteration")
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Stage 2.7 : The barotropic velocity update
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            call mpas_timer_start("si btr vel update")
   
#ifdef MPAS_OPENACC
            !$acc parallel loop gang &
            !$acc    present(edgeMask,cellsOnEdge,sshNew,dcEdge,    &
            !$acc       normalBarotropicVelocityNew,sshSubcycleCur, &
            !$acc       normalBarotropicVelocitySubcycleNew,        &
            !$acc       normalBarotropicVelocityCur,                &
            !$acc       barotropicCoriolisTerm,barotropicForcing,   &
            !$acc       sshSubcycleNew,nEdgesHalo)                  &
            !$acc    private(temp_mask,cell1,cell2)
#else
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(temp_mask,cell1,cell2)
#endif
            do iEdge = 1, nEdgesHalo(2)
               temp_mask = edgeMask(1, iEdge)
   
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
   
               ! Use of normBtrVelNew as a temp array for sshEdge
               normalBarotropicVelocityNew(iEdge) =    &
                          ( 0.5*sshSubcycleNew(cell2)  &
                           +0.5*sshSubcycleCur(cell2)) &
                         -( 0.5*sshSubcycleNew(cell1)  &
                           +0.5*sshSubcycleCur(cell1)) 
   
               normalBarotropicVelocitySubcycleNew(iEdge)              &
                  = temp_mask                                          &
                  * (normalBarotropicVelocityCur(iEdge)                &
                  - dt_si * (-barotropicCoriolisTerm(iEdge) + gravity  &
                          *normalBarotropicVelocityNew(iEdge)          &
                          /dcEdge(iEdge) - barotropicForcing(iEdge)))
            end do ! iEdge
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

#ifdef MPAS_OPENACC
            !$acc parallel loop gang  &
            !$acc    present(nEdgesOnEdge,edgesOnEdge,weightsOnEdge, &
            !$acc        normalBarotropicVelocitySubcycleNew,        &
            !$acc        normalBarotropicVelocityCur,fEdge,          &
            !$acc        barotropicCoriolisTerm,nEdgesHalo)          &
            !$acc    private(CoriolisTerm,eoe,i)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(CoriolisTerm, i, eoe)
#endif
            do iEdge = 1,nEdgesHalo(1)
               ! Compute the barotropic Coriolis term, -f*uPerp
               CoriolisTerm = 0.0_RKIND
               !$acc loop vector reduction(+:CoriolisTerm)
               do i = 1, nEdgesOnEdge(iEdge)
                  eoe = edgesOnEdge(i,iEdge)
                  CoriolisTerm =  CoriolisTerm + weightsOnEdge(i,iEdge)     &
                               * 0.5_RKIND                                  &
                               * ( normalBarotropicVelocitySubcycleNew(eoe) &
                                  +normalBarotropicVelocityCur(eoe) )       &
                               * fEdge(eoe)
               end do
               barotropicCoriolisTerm(iEdge) = CoriolisTerm
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
   
#ifdef MPAS_OPENACC
            !$acc parallel loop gang                                &
            !$acc    present(edgeMask,cellsOnEdge,dcEdge,           &
            !$acc       normalBarotropicVelocityNew,                &
            !$acc       normalBarotropicVelocitySubcycleNew,        &
            !$acc       normalBarotropicVelocityCur,                &
            !$acc       barotropicCoriolisTerm,barotropicForcing)   &
            !$acc    private(temp_mask,cell1,cell2)
#else
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(temp_mask,cell1,cell2)
#endif
            do iEdge = 1, nEdgesOwned
               temp_mask = edgeMask(1, iEdge)
   
               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
   
               normalBarotropicVelocitySubcycleNew(iEdge)              &
                  = temp_mask                                          &
                  * (normalBarotropicVelocityCur(iEdge)                &
                  - dt_si * (-barotropicCoriolisTerm(iEdge) + gravity  &
                          *normalBarotropicVelocityNew(iEdge)          &
                          /dcEdge(iEdge) - barotropicForcing(iEdge)))
            end do ! iEdge
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
   
#ifdef MPAS_OPENACC
            !$acc parallel loop gang &
            !$acc    present(cellsOnEdge,sshSubcycleCur,sshSubcycleNew, &
            !$acc       normalBarotropicVelocitySubcycleCur,            &
            !$acc       normalBarotropicVelocitySubcycleNew,            &
            !$acc       normalBarotropicVelocityCur,dcEdge,             &
            !$acc       normalBarotropicVelocityNew,edgeMask,           &
            !$acc       bottomDepthEdge,barotropicThicknessFlux,        &
            !$acc       barotropicCoriolisTerm,barotropicForcing,       &
            !$acc       maxLevelEdgeTop)                                &
            !$acc    private(cell1,cell2,sshEdge,thicknessSum)
#else
            !$omp parallel
            !$omp do schedule(runtime) &
            !$omp private(cell1,cell2,sshEdge,thicknessSum)
#endif
            do iEdge = 1, nEdgesOwned
               if (maxLevelEdgeTop(iEdge).eq.0) cycle

               cell1 = cellsOnEdge(1,iEdge)
               cell2 = cellsOnEdge(2,iEdge)
   
               ! normBtrVelSubCur is at time n+0.5 if splitImplicitStep=1
               !                     at time n+1.0 if splitImplicitStep=2
               normalBarotropicVelocitySubcycleCur(iEdge)                &
                  = 0.5_RKIND*normalBarotropicVelocitySubcycleNew(iEdge) &
                  + 0.5_RKIND*normalBarotropicVelocityCur(iEdge)
   
                         ! 0.25 = 0.5 * 0.5
               sshEdge = 0.25_RKIND * (  sshSubcycleCur(cell1)   &
                                       + sshSubcycleCur(cell2) ) &
                       + 0.25_RKIND * (  sshSubcycleNew(cell1)   &
                                       + sshSubcycleNew(cell2) )
   
               thicknessSum = sshEdge + bottomDepthEdge(iEdge)
   
               barotropicThicknessFlux(iEdge) &
                  = 0.5*(  normalBarotropicVelocitySubcycleNew(iEdge) &
                         + normalBarotropicVelocityCur(iEdge) )       &
                         * thicknessSum
            end do ! iEdge
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
   
            ! boundary update on F
            call mpas_timer_start("si halo btr vel")
            !$acc update host(barotropicThicknessFlux, &
            !$acc             normalBarotropicVelocitySubcycleCur)
            call mpas_dmpar_exch_group_create(domain, finalBtrGroupName)
            call mpas_dmpar_exch_group_add_field(domain, &
                      finalBtrGroupName, 'barotropicThicknessFlux')
            call mpas_dmpar_exch_group_add_field(domain,             &
                      finalBtrGroupName,                             &
                                 'normalBarotropicVelocitySubcycle', &
                                                       timeLevel=1)
            call mpas_dmpar_exch_group_full_halo_exch(domain, &
                      finalBtrGroupName)
            call mpas_dmpar_exch_group_destroy(domain, finalBtrGroupName)
            !$acc update device(barotropicThicknessFlux, &
            !$acc               normalBarotropicVelocitySubcycleCur)
            call mpas_timer_stop("si halo btr vel")
   
#ifdef MPAS_OPENACC
            !$acc parallel loop gang &
            !$acc    present(normalBarotropicVelocityNew, &
            !$acc            normalBarotropicVelocitySubcycleCur)
#else
            !$omp parallel
            !$omp do schedule(runtime)
#endif
            do iEdge = 1, nEdgesAll
               normalBarotropicVelocityNew(iEdge)              &
                  = normalBarotropicVelocitySubcycleCur(iEdge) 
            end do ! iEdge
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

            call mpas_timer_stop("si btr vel update")

         !-------------------------------------------------------------!
         end do ! siLargeIter
         !-------------------------------------------------------------!
         ! END   Large barotropic system iteration loop
         !-------------------------------------------------------------!

         call mpas_timer_stop("si btr vel")

         !$acc exit data copyout(sshNew,normalBarotropicVelocityNew)  &
         !$acc           delete (sshCur,                              &
         !$acc                   sshSubcycleCur,sshSubcycleNew,       &
         !$acc                   normalBarotropicVelocitySubcycleCur, &
         !$acc                   normalBarotropicVelocitySubcycleNew, &
         !$acc                   normalBarotropicVelocityCur,         &
         !$acc                   bottomDepthEdge)

         !$acc update host(barotropicThicknessFlux,barotropicForcing)

         ! Check that you can compute SSH using the total sum or the
         ! individual increments over the barotropic subcycles.
         ! efficiency: This next block of code is really a check for
         ! debugging, and can be removed later.
         call mpas_timer_start('btr si ssh verif')

         ! Correction velocity
         !   normalVelocityCorrection = (Flux - Sum(h u*))/H
         ! or, for the full latex version:
         !{\bf u}^{corr} = \left( {\overline {\bf F}}
         !  - \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}
         ! {\bf u}_k^{avg} \right)
         ! \left/ \sum_{k=1}^{N^{edge}} h_{k,*}^{edge}   \right.

         !Compute uTemp
         ! velocity for normalVelocityCorrection is
         ! normalBarotropicVelocity +
         ! normalBaroclinicVelocity + uBolus
         nEdges = nEdgesHalo(config_num_halos-1 )
         allocate(uTemp(nVertLevels,nEdges))

         !$omp parallel
         !$omp do schedule(runtime)
         do iEdge = 1, nEdges
            uTemp(:,iEdge) = normalBarotropicVelocityNew(iEdge) &
                         + normalBaroclinicVelocityNew(:,iEdge)
         end do
         !$omp end do
         !$omp end parallel

         call ocn_GM_add_to_transport_vel(uTemp, nEdges, nVertLevels)
         call ocn_MLE_add_to_transport_vel(uTemp, nEdges)

         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp private(k)
         do iEdge = 1, nEdges
            do k = 1, nVertLevels
               normalTransportVelocity(k,iEdge) =  &
                        normalBarotropicVelocityNew(iEdge)   + &
                        normalBaroclinicVelocityNew(k,iEdge)
            end do
         end do
         !$omp end do
         !$omp end parallel

         ! add GM and submesoscale contributions if requested
         call ocn_GM_add_to_transport_vel(normalTransportVelocity, nEdges, &
                                          nVertLevels)
         call ocn_MLE_add_to_transport_vel(normalTransportVelocity, nEdges)

         !$omp parallel
         !$omp do schedule(runtime) &
         !$omp private(k, normalThicknessFluxSum, &
         !$omp         thicknessSum, normalVelocityCorrection)
         do iEdge = 1, nEdges

            ! thicknessSum is initialized outside the loop because
            ! on land boundaries maxLevelEdgeTop=0, but I want to
            ! initialize thicknessSum with a nonzero value to avoid
            ! a NaN.
            normalThicknessFluxSum  &
               = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge) &
               * uTemp(minLevelEdgeBot(iEdge),iEdge)
            thicknessSum &
               = layerThickEdgeFlux(minLevelEdgeBot(iEdge),iEdge)

            do k = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
               normalThicknessFluxSum = normalThicknessFluxSum + &
                                        layerThickEdgeFlux(k,iEdge)* &
                                        uTemp(k,iEdge)
               thicknessSum = thicknessSum + &
                              layerThickEdgeFlux(k,iEdge)
            enddo

            normalVelocityCorrection = useVelocityCorrection* &
                          ((barotropicThicknessFlux(iEdge) -  &
                            normalThicknessFluxSum)/thicknessSum)

            do k = 1, nVertLevels

               ! normalTransportVelocity = normalBarotropicVelocity
               !                         + normalBaroclinicVelocity
               !                         + normalGMBolusVelocity
               !                         + normalMLEvelocity 
               !                         + normalVelocityCorrection
               ! This is u used in advective terms for layerThickness
               ! and tracers in tendency calls in stage 3.
               normalTransportVelocity(k,iEdge) = edgeMask(k,iEdge) &
                        *(normalTransportVelocity(k,iEdge) + & 
                          normalVelocityCorrection )
            enddo

         end do ! iEdge
         !$omp end do
         !$omp end parallel

         deallocate(uTemp)

         call mpas_timer_stop('btr si ssh verif')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Stage 3: Tracer, density, pressure, vert velocity prediction
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! only compute tendencies for active tracers on last large iteration
         if (splitImplicitStep < numTSIterations) then
            activeTracersOnly = .true.
         else
            activeTracersOnly = .false.
         endif

         ! Thickness tendency computations and thickness halo updates
         ! are completed before tracer tendency computations to allow
         ! monotonic advection.

         call mpas_timer_start('si thick tend')

         ! compute vertAleTransportTop.  Use normalTransportVelocity
         ! for advection of layerThickness and tracers.
         ! Use time level 1 values of layerThickness and
         ! layerThickEdgeFlux because layerThickness has not yet
         ! been computed for time level 2.
         call mpas_timer_start('thick vert trans vel top')
#ifdef MPAS_OPENACC
         !$acc enter data copyin(layerThicknessCur, sshCur)
         !$acc update device(layerThickEdgeFlux, normalTransportVelocity)
#endif
         if (associated(highFreqThicknessNew)) then
#ifdef MPAS_OPENACC
            !$acc enter data copyin(highFreqThicknessNew)
#endif
            call ocn_vert_transport_velocity_top( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalTransportVelocity, sshCur, &
                 dt, vertAleTransportTop, err, highFreqThicknessNew)
#ifdef MPAS_OPENACC
            !$acc exit data delete(highFreqThicknessNew)
#endif
         else
            call ocn_vert_transport_velocity_top( &
                 verticalMeshPool, layerThicknessCur, &
                 layerThickEdgeFlux, normalTransportVelocity, sshCur, &
                 dt, vertAleTransportTop, err)
         endif
#ifdef MPAS_OPENACC
         !$acc exit data delete(layerThicknessCur, sshCur)
         !$acc update host(vertAleTransportTop, normalTransportVelocity)
#endif
         call mpas_timer_stop('thick vert trans vel top')

         call ocn_tend_thick(tendPool, forcingPool)

         call mpas_timer_stop('si thick tend')

         ! update halo for thickness tendencies
         call mpas_timer_start("si halo thickness")

         call mpas_dmpar_field_halo_exch(domain, 'tendLayerThickness')

         call mpas_timer_stop("si halo thickness")

         call mpas_timer_start('si tracer tend', .false.)

         call ocn_tend_tracer(tendPool, statePool, forcingPool, &
                              meshPool, swForcingPool, &
                              dt, activeTracersOnly, 2)

         call mpas_pool_get_array(tracersPool, 'activeTracers', &
                                                tracersGroupCur, 1)
         call mpas_pool_get_array(tracersPool, 'activeTracers', &
                                                tracersGroupNew, 2)
         call mpas_pool_get_array(statePool,   'normalVelocity', &
                                                normalVelocityCur, 1)
         call mpas_pool_get_array(tracersTendPool,'activeTracersTend', &
                                                   activeTracersTend)

#ifdef MPAS_OPENACC
         !$acc enter data copyin(layerThicknessNew, normalVelocityNew)
         !$acc enter data copyin(normalBarotropicVelocityNew, &
         !$acc                   normalBaroclinicVelocityNew)
         !$acc enter data copyin(activeTracersNew)
         !$acc update device(layerThickEdgeFlux)
         if (config_use_freq_filtered_thickness) then
            !$acc enter data copyin(lowFreqDivergenceNew)
            !$acc enter data copyin(lowFreqDivergenceCur)
            !$acc enter data copyin(lowFreqDivergenceTend)
         endif
         if (splitImplicitStep < numTSIterations) then
            !$acc update device (normalTransportVelocity)
            !$acc enter data copyin(atmosphericPressure, seaIcePressure)
            !$acc enter data copyin(sshNew)
            !$acc update device(tracersSurfaceValue)
            if ( associated(normalGMBolusVelocity) ) then
               !$acc update device (normalGMBolusVelocity)
            end if
            if ( associated(normalMLEVelocity) ) then
               !$acc update device (normalMLEvelocity)
            end if
            if ( associated(frazilSurfacePressure) ) then
               !$acc enter data copyin(frazilSurfacePressure)
            endif
            if (landIcePressureOn) then
               !$acc enter data copyin(landIcePressure)
               !$acc enter data copyin(landIceDraft)
            endif
            if (config_use_freq_filtered_thickness) then
               !$acc enter data copyin(highFreqThicknessNew)
               !$acc enter data copyin(highFreqThicknessCur)
            endif
         elseif (splitImplicitStep == numTSIterations) then
            !$acc enter data copyin(layerThicknessCur)
            !$acc enter data copyin(layerThicknessTend)
            !$acc enter data copyin(normalBaroclinicVelocityCur)
            if (config_compute_active_tracer_budgets) then
               !$acc update device(activeTracerHorizontalAdvectionEdgeFlux)
               !$acc update device(activeTracerHorizontalAdvectionTendency)
               !$acc update device(activeTracerVerticalAdvectionTendency)
               !$acc update device(activeTracerHorMixTendency)
               !$acc update device(activeTracerSurfaceFluxTendency)
               !$acc update device(temperatureShortWaveTendency)
               !$acc update device(activeTracerNonLocalTendency)
            endif
         endif
#endif

         call mpas_timer_stop('si tracer tend')

         ! update halo for tracer tendencies
         call mpas_timer_start("si halo tracers")

         call mpas_pool_begin_iteration(tracersTendPool)
         do while (mpas_pool_get_next_member(tracersTendPool, &
                                             groupItr) )
            if (groupItr%memberType == MPAS_POOL_FIELD ) then
               ! Only compute tendencies for active tracers if
               ! activeTracersOnly flag is true.
               if (.not. activeTracersOnly .or. &
                   trim(groupItr%memberName)=='activeTracersTend') then
                  call mpas_dmpar_field_halo_exch(domain, &
                                                  groupItr%memberName)
               end if
            end if
         end do

         call mpas_timer_stop("si halo tracers")

         call mpas_timer_start('si loop fini')

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !  If iterating, reset variables for next iteration
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (splitImplicitStep < numTSIterations) then

            ! Get indices for dynamic tracers (Includes T&S).
            call mpas_pool_get_dimension_scalar(tracersPool,'activeGRP_start',&
                                                             startIndex)
            call mpas_pool_get_dimension_scalar(tracersPool,'activeGRP_end', &
                                                             endIndex)

            ! Only need T & S for earlier iterations,
            ! then all the tracers needed the last time through.

!! Keeping this calc on the CPU for now
#ifdef MPAS_OPENACC
            !$acc update host(layerThicknessNew)
            !$acc update host(activeTracersNew)
            !$acc update host(tracersSurfaceValue)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(i, k, temp_h, temp)
#endif
            do iCell = 1, nCellsAll
            do k = minLevelCell(iCell), maxLevelCell(iCell)

               ! this is h_{n+1}
               temp_h = layerThicknessCur(k,iCell) + dt* &
               layerThicknessTend(k,iCell)

               ! this is h_{n+1/2}
               layerThicknessNew(k,iCell) = 0.5* &
               (layerThicknessCur(k,iCell) + temp_h)

               do i = startIndex, endIndex
                  ! This is Phi at n+1
                  temp = (tracersGroupCur(i,k,iCell)* &
                          layerThicknessCur(k,iCell) + dt* &
                          activeTracersTend(i,k,iCell))/temp_h

                  ! This is Phi at n+1/2
                  tracersGroupNew(i,k,iCell) = 0.5_RKIND* &
                            (tracersGroupCur(i,k,iCell) + temp)
               end do ! tracer index
            end do ! vertical
            end do ! iCell
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

#ifdef MPAS_OPENACC
            !$acc update device(layerThicknessNew)
            !$acc update device(activeTracersNew)
#endif

            if (config_use_freq_filtered_thickness) then

#ifdef MPAS_OPENACC
               !$acc parallel loop present(minLevelCell, maxLevelCell) &
               !$acc   present(highFreqThicknessNew) &
               !$acc   present(highFreqThicknessCur) &
               !$acc   present(lowFreqDivergenceNew) &
               !$acc   present(lowFreqDivergenceCur) &
               !$acc   present(lowFreqDivergenceTend)
#else
               !$omp parallel
               !$omp do schedule(runtime) private(k, temp)
#endif
               do iCell = 1, nCellsAll
               do k = minLevelCell(iCell), maxLevelCell(iCell)

                  ! h^{hf}_{n+1} was computed in Stage 1

                  ! this is h^{hf}_{n+1/2}
                  highFreqThicknessNew(k,iCell) = 0.5_RKIND* &
                             (highFreqThicknessCur(k,iCell) + &
                              highFreqThicknessNew(k,iCell))

                  ! this is D^{lf}_{n+1}
                  temp = lowFreqDivergenceCur(k,iCell) + dt* &
                         lowFreqDivergenceTend(k,iCell)

                  ! this is D^{lf}_{n+1/2}
                  lowFreqDivergenceNew(k,iCell) = 0.5_RKIND* &
                           (lowFreqDivergenceCur(k,iCell) + temp)
               end do
               end do
#ifndef MPAS_OPENACC
               !$omp end do
               !$omp end parallel
#endif
            end if

#ifdef MPAS_OPENACC
            !$acc parallel loop collapse(2) present(normalVelocityNew) &
            !$acc   present(normalBaroclinicVelocityNew, edgeMask) &
            !$acc   present(normalBarotropicVelocityNew)
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k)
#endif
            do iEdge = 1, nEdgesAll
            do k = 1, nVertLevels

               ! u = normalBarotropicVelocity + normalBaroclinicVelocity
               ! here normalBaroclinicVelocity is at time n+1/2
               ! This is u used in next iteration or step
               normalVelocityNew(k,iEdge) = edgeMask(k,iEdge)* &
                         ( normalBarotropicVelocityNew(iEdge) + &
                           normalBaroclinicVelocityNew(k,iEdge) )

            enddo
            end do ! iEdge
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

            ! Efficiency note: We really only need this to compute
            ! layerThickEdgeFlux, density, pressure, and SSH
            ! in this diagnostics solve.
            call ocn_diagnostic_solve(dt, statePool, forcingPool, &
                                      meshPool, verticalMeshPool, scratchPool, &
                                      tracersPool, 2, full=.false.)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! If large iteration complete, compute all variables at
         ! time n+1
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         elseif (splitImplicitStep == numTSIterations) then

#ifdef MPAS_OPENACC
            !$acc parallel present(minLevelCell, maxLevelCell) &
            !$acc   present(layerThicknessTend) &
            !$acc   present(layerThicknessCur) &
            !$acc   present(minLevelEdgeBot, maxLevelEdgeTop) &
            !$acc   present(layerThickEdgeFlux) &
            !$acc   present(activeTracerHorizontalAdvectionEdgeFlux) &
            !$acc   present(layerThicknessNew) &
            !$acc   present(activeTracerHorizontalAdvectionTendency) &
            !$acc   present(activeTracerVerticalAdvectionTendency) &
            !$acc   present(activeTracerHorMixTendency) &
            !$acc   present(activeTracerSurfaceFluxTendency) &
            !$acc   present(temperatureShortWaveTendency) &
            !$acc   present(activeTracerNonLocalTendency)
#endif

#ifdef MPAS_OPENACC
            !$acc loop gang
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k)
#endif
            do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
            !$acc loop vector
#endif
            do k = minLevelCell(iCell), maxLevelCell(iCell)
               ! this is h_{n+1}
               layerThicknessNew(k,iCell) = &
                                     layerThicknessCur(k,iCell) + &
                                  dt*layerThicknessTend(k,iCell)
            end do
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

            if (config_compute_active_tracer_budgets) then

#ifdef MPAS_OPENACC
               !$acc loop gang
#else
               !$omp parallel
               !$omp do schedule(runtime) private(k)
#endif
               do iEdge = 1, nEdgesAll
#ifdef MPAS_OPENACC
               !$acc loop seq
#endif
               do k= minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
#ifdef MPAS_OPENACC
               !$acc loop vector
#endif
                  do iTracer = 1, size(activeTracerHorizontalAdvectionEdgeFlux,1)
                     activeTracerHorizontalAdvectionEdgeFlux(iTracer,k,iEdge) = &
                     activeTracerHorizontalAdvectionEdgeFlux(iTracer,k,iEdge) / &
                          layerThickEdgeFlux(k,iEdge)
                  enddo
               enddo
               enddo
#ifndef MPAS_OPENACC
               !$omp end do
#endif

#ifdef MPAS_OPENACC
               !$acc loop gang
#else
               !$omp do schedule(runtime) private(k)
#endif
               do iCell = 1, nCellsAll
               do k= minLevelCell(iCell), maxLevelCell(iCell)
#ifdef MPAS_OPENACC
               !$acc loop vector
                  do iTracer = 1, size(activeTracerHorizontalAdvectionTendency,1)
                     activeTracerHorizontalAdvectionTendency(iTracer,k,iCell) = &
                     activeTracerHorizontalAdvectionTendency(iTracer,k,iCell) / &
                               layerThicknessNew(k,iCell)

                     activeTracerVerticalAdvectionTendency(iTracer,k,iCell) = &
                     activeTracerVerticalAdvectionTendency(iTracer,k,iCell) / &
                               layerThicknessNew(k,iCell)

                     activeTracerHorMixTendency(iTracer,k,iCell) = &
                     activeTracerHorMixTendency(iTracer,k,iCell) / &
                                layerThicknessNew(k,iCell)

                     activeTracerSurfaceFluxTendency(iTracer,k,iCell) = &
                     activeTracerSurfaceFluxTendency(iTracer,k,iCell) / &
                                layerThicknessNew(k,iCell)

                     activeTracerNonLocalTendency(iTracer,k,iCell) = &
                     activeTracerNonLocalTendency(iTracer,k,iCell) / &
                                layerThicknessNew(k,iCell)
                  end do
#else
                  activeTracerHorizontalAdvectionTendency(:,k,iCell) = &
                  activeTracerHorizontalAdvectionTendency(:,k,iCell) / &
                            layerThicknessNew(k,iCell)

                  activeTracerVerticalAdvectionTendency(:,k,iCell) = &
                  activeTracerVerticalAdvectionTendency(:,k,iCell) / &
                            layerThicknessNew(k,iCell)

                  activeTracerHorMixTendency(:,k,iCell) = &
                  activeTracerHorMixTendency(:,k,iCell) / &
                             layerThicknessNew(k,iCell)

                  activeTracerSurfaceFluxTendency(:,k,iCell) = &
                  activeTracerSurfaceFluxTendency(:,k,iCell) / &
                             layerThicknessNew(k,iCell)

                  activeTracerNonLocalTendency(:,k,iCell) = &
                  activeTracerNonLocalTendency(:,k,iCell) / &
                             layerThicknessNew(k,iCell)
#endif

                  temperatureShortWaveTendency(k,iCell) = &
                  temperatureShortWaveTendency(k,iCell) / &
                             layerThicknessNew(k,iCell)
               end do
               end do
#ifndef MPAS_OPENACC
               !$omp end do
               !$omp end parallel
#endif
            endif

#ifdef MPAS_OPENACC
            !$acc end parallel
#endif

! This tracer block is still computed on the CPU
#ifdef MPAS_OPENACC
            !$acc update host(layerThicknessNew)
            !$acc update host(layerThicknessCur, layerThicknessTend)
            !$acc update host(activeTracersNew)
            !$acc update host(tracersSurfaceValue)
#endif

            call mpas_pool_begin_iteration(tracersPool)
            do while (mpas_pool_get_next_member(tracersPool, &
                                                groupItr) )
               if (groupItr%memberType == MPAS_POOL_FIELD) then
                  configName = 'config_use_'//trim(groupItr%memberName)
                  call mpas_pool_get_config(domain%configs, &
                                configName, config_use_tracerGroup)

                  if ( config_use_tracerGroup ) then
                     call mpas_pool_get_array(tracersPool, &
                                              groupItr%memberName, &
                                              tracersGroupCur, 1)
                     call mpas_pool_get_array(tracersPool, &
                                              groupItr%memberName, &
                                              tracersGroupNew, 2)

                     modifiedGroupName = &
                             trim(groupItr % memberName) // 'Tend'
                     call mpas_pool_get_array(tracersTendPool, &
                                              modifiedGroupName, &
                                              tracersGroupTend)

                     !$omp parallel
                     !$omp do schedule(runtime) private(k)
                     do iCell = 1, nCellsAll
                     do k = minLevelCell(iCell), maxLevelCell(iCell)
                        tracersGroupNew(:,k,iCell) = &
                       (tracersGroupCur(:,k,iCell) * &
                        layerThicknessCur(k,iCell) + dt* &
                       tracersGroupTend(:,k,iCell))/ &
                        layerThicknessNew(k,iCell)
                     end do
                     end do
#ifndef MPAS_OPENACC
                     !$omp end do
                     !$omp end parallel
#endif

                     ! limit salinity in separate loop
                     if (trim(groupItr%memberName) == &
                         'activeTracers' ) then
                        !$omp parallel
                        !$omp do schedule(runtime) private(k)
                        do iCell = 1, nCellsAll
                        do k = minLevelCell(iCell), maxLevelCell(iCell)
                           tracersGroupNew(indexSalinity,k,iCell) = &
                           max(0.001_RKIND,  &
                           tracersGroupNew(indexSalinity,k,iCell))
                        end do
                        end do
#ifndef MPAS_OPENACC
                        !$omp end do
                        !$omp end parallel
#endif
                     end if

                     ! Reset debugTracers to fixed value at the surface
                     if (trim(groupItr % memberName) == &
                         'debugTracers' .and. &
                         config_reset_debugTracers_near_surface) then

                        !$omp parallel
                        !$omp do schedule(runtime) private(k, lat)
                        do iCell = 1, nCellsAll

                           ! Reset tracer1 to 2 in top n layers
                           do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                tracersGroupNew(1,k,iCell) = 2.0_RKIND
                           end do

                           ! Reset tracer2 to 2 in top n layers
                           ! in zonal bands, and 1 outside
                           lat = latCell(iCell)*180.0_RKIND/pi
                           if (     lat>-60.0.and.lat<-55.0 &
                                .or.lat>-40.0.and.lat<-35.0 &
                                .or.lat>- 2.5.and.lat<  2.5 &
                                .or.lat> 35.0.and.lat< 40.0 &
                                .or.lat> 55.0.and.lat< 60.0 ) then
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(2,k,iCell) = 2.0_RKIND
                              end do
                           else
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(2,k,iCell) = 1.0_RKIND
                              end do
                           end if

                           ! Reset tracer3 to 2 in top n layers
                           ! in zonal bands, and 1 outside
                           lat = latCell(iCell)*180.0_RKIND/pi
                           if (     lat>-55.0.and.lat<-50.0 &
                                .or.lat>-35.0.and.lat<-30.0 &
                                .or.lat>-15.0.and.lat<-10.0 &
                                .or.lat> 10.0.and.lat< 15.0 &
                                .or.lat> 30.0.and.lat< 35.0 &
                                .or.lat> 50.0.and.lat< 55.0 ) then
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(3,k,iCell) = 2.0_RKIND
                              end do
                           else
                              do k = minLevelCell(iCell), minLevelCell(iCell)+config_reset_debugTracers_top_nLayers-1
                                 tracersGroupNew(3,k,iCell) = 1.0_RKIND
                              end do
                           end if
                        end do ! cells
#ifndef MPAS_OPENACC
                        !$omp end do
                        !$omp end parallel
#endif
                     end if ! debug tracers
                  end if ! use tracer group
               end if ! tracer
            end do ! tracer group

#ifdef MPAS_OPENACC
            !$acc update device(layerThicknessNew)
            !$acc update device(activeTracersNew)
#endif

            if (config_use_freq_filtered_thickness) then
#ifdef MPAS_OPENACC
               !$acc parallel present(minLevelCell, maxLevelCell) &
               !$acc   present(lowFreqDivergenceNew) &
               !$acc   present(lowFreqDivergenceCur) &
               !$acc   present(lowFreqDivergenceTend)
#endif

#ifdef MPAS_OPENACC
               !$acc loop gang
#else
               !$omp parallel
               !$omp do schedule(runtime) private(k)
#endif
               do iCell = 1, nCellsAll
#ifdef MPAS_OPENACC
               !$acc loop vector
#endif
               do k = minLevelCell(iCell), maxLevelCell(iCell)

                  ! h^{hf}_{n+1} was computed in Stage 1

                  ! this is D^{lf}_{n+1}
                  lowFreqDivergenceNew(k,iCell) = &
                  lowFreqDivergenceCur(k,iCell) + dt* &
                  lowFreqDivergenceTend(k,iCell)
               end do
               end do
#ifndef MPAS_OPENACC
               !$omp end do
               !$omp end parallel
#endif

#ifdef MPAS_OPENACC
               !$acc end parallel
#endif

            end if

            ! Recompute final u to go on to next step.
            ! u_{n+1} = normalBarotropicVelocity_{n+1} +
            !           normalBaroclinicVelocity_{n+1}
            ! Right now normalBaroclinicVelocityNew is at time n+1/2,
            ! so back compute to get normalBaroclinicVelocity at
            ! time n+1 using normalBaroclinicVelocity_{n+1/2} =
            !            1/2*(normalBaroclinicVelocity_n + u_Bcl_{n+1})
            ! so the following lines are
            ! u_{n+1} = normalBarotropicVelocity_{n+1} +
            !         2*normalBaroclinicVelocity_{n+1/2} -
            !         normalBaroclinicVelocity_n
            ! note that normalBaroclinicVelocity is recomputed at the
            ! beginning of the next timestep due to Imp Vert mixing,
            ! so normalBaroclinicVelocity does not have to be
            ! recomputed here.

#ifdef MPAS_OPENACC
            !$acc parallel present(minLevelEdgeBot, maxLevelEdgeTop) &
            !$acc   present(normalVelocityNew) &
            !$acc   present(normalBarotropicVelocityNew) &
            !$acc   present(normalBaroclinicVelocityNew) &
            !$acc   present(normalBaroclinicVelocityCur)
#endif

#ifdef MPAS_OPENACC
            !$acc loop gang
#else
            !$omp parallel
            !$omp do schedule(runtime) private(k)
#endif
            do iEdge = 1, nEdgesAll
#ifdef MPAS_OPENACC
            !$acc loop vector
#endif
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               normalVelocityNew(k,iEdge) = &
                               normalBarotropicVelocityNew(iEdge) + &
                             2*normalBaroclinicVelocityNew(k,iEdge) - &
                               normalBaroclinicVelocityCur(k,iEdge)
            end do
            end do ! iEdges
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

#ifdef MPAS_OPENACC
            !$acc end parallel
#endif

         endif ! splitImplicitStep

#ifdef MPAS_OPENACC
         !$acc exit data copyout(layerThicknessNew)
         !$acc exit data copyout(normalVelocityNew)
         !$acc exit data delete(normalBarotropicVelocityNew, &
         !$acc                  normalBaroclinicVelocityNew)
         !$acc exit data copyout(activeTracersNew)
         !$acc update host(layerThickEdgeFlux)
         if (config_use_freq_filtered_thickness) then
            !$acc exit data copyout(lowFreqDivergenceNew)
            !$acc exit data delete(lowFreqDivergenceCur)
            !$acc exit data delete(lowFreqDivergenceTend)
         endif
         if (splitImplicitStep < numTSIterations) then
            if (config_use_gm) then
               !$acc update host(vertGMBolusVelocityTop)
            end if
            if (config_submesoscale_enable) then
               !$acc update host(vertMLEBolusVelocityTop)
            end if
            !$acc exit data delete (atmosphericPressure, seaIcePressure)
            !$acc exit data copyout(sshNew)
            !$acc update host(layerThickEdgeMean)
            !$acc update host(relativeVorticity, circulation)
            !$acc update host(vertTransportVelocityTop, &
            !$acc             relativeVorticityCell, &
            !$acc             divergence, &
            !$acc             kineticEnergyCell, &
            !$acc             tangentialVelocity, &
            !$acc             vertVelocityTop)
            !$acc update host(normRelVortEdge, normPlanetVortEdge, &
            !$acc             normalizedRelativeVorticityCell)
            !$acc update host (surfacePressure)
            !$acc update host(zMid, zTop)
            !$acc update host(tracersSurfaceValue)
            !$acc update host(normalVelocitySurfaceLayer)
            !$acc update host(density, potentialDensity, displacedDensity)
            !$acc update host(thermExpCoeff,  &
            !$acc&            salineContractCoeff)
            !$acc update host(montgomeryPotential, pressure)
            if (landIcePressureOn) then
               !$acc exit data delete(landIcePressure)
               !$acc exit data delete(landIceDraft)
            endif
            if (config_use_freq_filtered_thickness) then
               !$acc exit data copyout(highFreqThicknessNew)
               !$acc exit data delete(highFreqThicknessCur)
            endif
            if ( associated(frazilSurfacePressure) ) then
               !$acc exit data delete(frazilSurfacePressure)
            endif
         elseif (splitImplicitStep == numTSIterations) then
            !$acc exit data delete(layerThicknessCur, layerThicknessTend)
            !$acc exit data delete(normalBaroclinicVelocityCur)
            if (config_compute_active_tracer_budgets) then
               !$acc update host(activeTracerHorizontalAdvectionEdgeFlux)
               !$acc update host(activeTracerHorizontalAdvectionTendency)
               !$acc update host(activeTracerVerticalAdvectionTendency)
               !$acc update host(activeTracerHorMixTendency)
               !$acc update host(activeTracerSurfaceFluxTendency)
               !$acc update host(temperatureShortWaveTendency)
               !$acc update host(activeTracerNonLocalTendency)
            endif
         endif
#endif

         call mpas_timer_stop('si loop fini')
         call mpas_timer_stop('si loop')

      end do  ! outer timestep loop (splitImplicitStep)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END large iteration loop
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPAS_OPENACC
      if (config_use_tidal_potential_forcing) then
         !$acc exit data delete(pgf_sal,tidalPotentialEta)
         !$acc exit data copyout(sshSubcycleCurWithTides)
      endif
#endif

      call mpas_timer_start("si implicit vert mix")

      ! Call ocean diagnostic solve in preparation for vertical mixing.
      ! Note it is called again after vertical mixing, because u and
      ! tracers change. For Richardson vertical mixing, only density,
      ! layerThickEdgeFlux, and kineticEnergyCell need to be computed.
      ! For kpp, more variables may be needed.  Either way, this could
      ! be made more efficient by only computing what is needed for the
      ! implicit vmix routine that follows.
#ifdef MPAS_OPENACC
      !$acc enter data copyin(layerThicknessNew, normalVelocityNew)
      !$acc update device (normalTransportVelocity)
      if (config_use_gm) then
         !$acc update device (normalGMBolusVelocity)
      end if
      if (config_submesoscale_enable) then
         !$acc update device (normalMLEvelocity)
      end if
      !$acc enter data copyin(atmosphericPressure, seaIcePressure)
      !$acc enter data copyin(sshNew)
      !$acc enter data copyin(activeTracersNew)
      !$acc update device(tracersSurfaceValue)
      if ( associated(frazilSurfacePressure) ) then
         !$acc enter data copyin(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc enter data copyin(landIcePressure)
         !$acc enter data copyin(landIceDraft)
      endif
#endif
      call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, verticalMeshPool, &
                                scratchPool, tracersPool, 2)

      call mpas_dmpar_field_halo_exch(domain, 'surfaceFrictionVelocity')
      call mpas_dmpar_field_halo_exch(domain, 'tangentialVelocity')

#ifdef MPAS_OPENACC
      !$acc update host(layerThickEdgeFlux, layerThickEdgeMean)
      !$acc update host(relativeVorticity, circulation)
      if (config_use_gm) then
         !$acc update host (vertGMBolusVelocityTop)
      end if
      if (config_submesoscale_enable) then
         !$acc update host (vertMLEBolusVelocityTop)
      end if
      !$acc update host(vertTransportVelocityTop, &
      !$acc             relativeVorticityCell, &
      !$acc             divergence, &
      !$acc             kineticEnergyCell, &
      !$acc             tangentialVelocity, &
      !$acc             vertVelocityTop)
      !$acc update host(normRelVortEdge, normPlanetVortEdge, &
      !$acc             normalizedRelativeVorticityCell)
      !$acc update host (surfacePressure)
      !$acc update host(zMid, zTop)
      !$acc exit data copyout(sshNew)
      !$acc exit data delete(activeTracersNew)
      !$acc update host(tracersSurfaceValue)
      !$acc update host(normalVelocitySurfaceLayer)
      !$acc exit data delete (atmosphericPressure, seaIcePressure)
      !$acc update host(density, potentialDensity, displacedDensity)
      !$acc update host(thermExpCoeff,  &
      !$acc&            salineContractCoeff)
      !$acc update host(montgomeryPotential, pressure)
      !$acc update host(RiTopOfCell, &
      !$acc             BruntVaisalaFreqTop)
      !$acc update host(tracersSurfaceLayerValue, &
      !$acc             indexSurfaceLayerDepth, &
      !$acc             normalVelocitySurfaceLayer, &
      !$acc             sfcFlxAttCoeff, &
      !$acc             surfaceFluxAttenuationCoefficientRunoff)
      if ( associated(frazilSurfacePressure) ) then
         !$acc exit data delete(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc exit data delete(landIcePressure)
         !$acc exit data delete(landIceDraft)
      endif
      !$acc exit data delete(layerThicknessNew, normalVelocityNew)
#endif

      call ocn_vmix_implicit(dt, meshPool, statePool, forcingPool, &
                             scratchPool, err, 2)

      ! Update halo on u and tracers, which were just updated for
      ! implicit vertical mixing.  If not done, this leads to lack of
      ! volume conservation.  It is required because halo updates in
      ! stage 3 are only conducted on tendencies, not on the velocity
      ! and tracer fields.  So this update is required to communicate
      ! the change due to implicit vertical mixing across the boundary.

      call mpas_timer_start('si vmix halos')

      call mpas_timer_start('si vmix halos normalVelFld')
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', &
                                      timeLevel=2)
      call mpas_timer_stop('si vmix halos normalVelFld')

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if (groupItr%memberType == MPAS_POOL_FIELD) then

            ! Reset iAge to zero where mask == 0
            if (config_use_idealAgeTracers.and.trim(groupItr%memberName) == 'idealAgeTracers') then
                call mpas_pool_get_array(tracersPool, &
                                         groupItr%memberName, &
                                         tracersGroupNew, 2)
                call mpas_pool_get_subpool(forcingPool, &
                                           'tracersIdealAgeFields', &
                                            tracersIdealAgeFieldsPool)
                call mpas_pool_get_array(tracersIdealAgeFieldsPool, &
                                         'idealAgeTracersIdealAgeMask', &
                                          tracerGroupIdealAgeMask)

                !$omp parallel
                !$omp do schedule(runtime)
                do iCell = 1, nCellsOwned
                   tracersGroupNew(1,1,iCell) = tracerGroupIdealAgeMask(1,iCell)*tracersGroupNew(1,1,iCell)
                end do ! cells
                !$omp end do
                !$omp end parallel
             endif

            ! Halo update on all tracer groups
            call mpas_dmpar_field_halo_exch(domain, groupItr%memberName, timeLevel=2)

         end if
      end do

      call mpas_timer_stop('si vmix halos')

      call mpas_timer_stop("si implicit vert mix")

      if ( configVertAdvMethod == vertAdvRemap ) then

         call mpas_timer_start("vertical remap")

         call mpas_pool_get_array(statePool, 'layerThicknessLag', layerThicknessLagNew, 2)

         ! Store the current, Lagrangian location for computation of the
         ! Eulerian vertical velocity induced by remapping
         layerThicknessLagNew = layerThicknessNew

         ! Perform the remapping
         call ocn_remap_vert_state(block, err)

         call mpas_timer_stop("vertical remap")

      end if

      call mpas_timer_start('si fini')

#ifdef MPAS_OPENACC
      !$acc enter data copyin(layerThicknessNew, normalVelocityNew)
      if (config_prescribe_velocity) then
         !$acc enter data copyin(normalVelocityCur) 
      endif
      if (config_prescribe_thickness) then
         !$acc enter data copyin(layerThicknessCur)
      endif
      !$acc update device (normalTransportVelocity)
      if (config_use_gm) then
         !$acc update device (normalGMBolusVelocity)
      end if
      if (config_submesoscale_enable) then
         !$acc update device (normalMLEvelocity)
      end if
      !$acc enter data copyin(atmosphericPressure, seaIcePressure)
      !$acc enter data copyin(sshNew)
      !$acc enter data copyin(activeTracersNew)
      !$acc update device(tracersSurfaceValue)
      if ( associated(frazilSurfacePressure) ) then
         !$acc enter data copyin(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc enter data copyin(landIcePressure)
         !$acc enter data copyin(landIceDraft)
      endif
#endif

      if (config_prescribe_velocity) then
#ifdef MPAS_OPENACC
         !$acc parallel loop collapse(2) present(normalVelocityNew, normalVelocityCur)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iEdge = 1, nEdgesAll
         do k=1,nVertLevels
            normalVelocityNew(k,iEdge) = normalVelocityCur(k,iEdge)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      if (config_prescribe_thickness) then
#ifdef MPAS_OPENACC
         !$acc parallel loop collapse(2) present(layerThicknessNew, layerThicknessCur)
#else
         !$omp parallel
         !$omp do schedule(runtime) private(k)
#endif
         do iCell = 1, nCellsAll
         do k=1,nVertLevels
            layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
         end do
         end do
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
      end if

      call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, verticalMeshPool, &
                                scratchPool, tracersPool, 2)

      ! Update the effective desnity in land ice if we're coupling to
      ! land ice
      call ocn_effective_density_in_land_ice_update(forcingPool, &
                                                    statePool, err)

      call mpas_timer_start('si final mpas reconstruct', .false.)

#ifdef MPAS_OPENACC
      !$acc enter data create(velocityX, velocityY, velocityZ)
      !$acc enter data create(velocityZonal, velocityMeridional)
      !$acc enter data create(gradSSHX, gradSSHY, gradSSHZ)
      !$acc enter data create(gradSSHZonal, gradSSHMeridional)
      !$acc enter data create(surfaceVelocity, SSHGradient)
#endif

      call mpas_reconstruct_gpu(meshPool, normalVelocityNew,       &
                                velocityX, velocityY, velocityZ,   &
                                velocityZonal, velocityMeridional, &
                                includeHalos = .true.)

      call mpas_reconstruct_gpu(meshPool, gradSSH,               &
                                gradSSHX, gradSSHY, gradSSHZ,    &
                                gradSSHZonal, gradSSHMeridional, &
                                includeHalos = .true.)

      call mpas_timer_stop('si final mpas reconstruct')


#ifdef MPAS_OPENACC
      !$acc parallel loop present(surfaceVelocity, velocityZonal, velocityMeridional, &
      !$acc   SSHGradient, gradSSHZonal, gradSSHMeridional)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         surfaceVelocity(indexSurfaceVelocityZonal,iCell) = &
                                   velocityZonal(minLevelCell(iCell),iCell)
         surfaceVelocity(indexSurfaceVelocityMeridional,iCell) = &
                                   velocityMeridional(minLevelCell(iCell),iCell)

         SSHGradient(indexSSHGradientZonal,iCell) = gradSSHZonal(iCell)
         SSHGradient(indexSSHGradientMeridional,iCell) = &
                                               gradSSHMeridional(iCell)
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

#ifdef MPAS_OPENACC
      !$acc update host(layerThickEdgeFlux, layerThickEdgeMean)
      !$acc update host(relativeVorticity, circulation)
      if (config_use_gm) then
         !$acc update host(vertGMBolusVelocityTop)
      end if
      if (config_submesoscale_enable) then
         !$acc update host(vertMLEBolusVelocityTop)
      end if
      !$acc update host(vertTransportVelocityTop, &
      !$acc             relativeVorticityCell, &
      !$acc             divergence, &
      !$acc             kineticEnergyCell, &
      !$acc             tangentialVelocity, &
      !$acc             vertVelocityTop)
      !$acc update host(normRelVortEdge, normPlanetVortEdge, &
      !$acc             normalizedRelativeVorticityCell)
      !$acc update host (surfacePressure)
      !$acc update host(zMid, zTop)
      !$acc exit data copyout(sshNew)
      !$acc exit data delete(activeTracersNew)
      !$acc update host(tracersSurfaceValue)
      !$acc update host(normalVelocitySurfaceLayer)
      !$acc exit data delete (atmosphericPressure, seaIcePressure)
      if ( associated(frazilSurfacePressure) ) then
         !$acc exit data delete(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc exit data delete(landIcePressure)
         !$acc exit data delete(landIceDraft)
      endif
      !$acc exit data copyout(layerThicknessNew, normalVelocityNew)
      !$acc update host(density, potentialDensity, displacedDensity)
      !$acc update host(thermExpCoeff,  &
      !$acc&            salineContractCoeff)
      !$acc update host(montgomeryPotential, pressure)
      !$acc update host(RiTopOfCell, &
      !$acc             BruntVaisalaFreqTop)
      !$acc update host(tracersSurfaceLayerValue, &
      !$acc             indexSurfaceLayerDepth, &
      !$acc             normalVelocitySurfaceLayer, &
      !$acc             sfcFlxAttCoeff, &
      !$acc             surfaceFluxAttenuationCoefficientRunoff)
      if (config_prescribe_velocity) then
         !$acc exit data delete(normalVelocityCur)
      endif
      if (config_prescribe_thickness) then
         !$acc exit data delete(layerThicknessCur)
      endif
      !$acc exit data copyout(velocityX, velocityY, velocityZ)
      !$acc exit data copyout(gradSSHX, gradSSHY, gradSSHZ)
      !$acc exit data copyout(velocityZonal, velocityMeridional)
      !$acc exit data copyout(gradSSHZonal, gradSSHMeridional)
      !$acc exit data copyout(surfaceVelocity, SSHGradient)
#endif
      call ocn_time_average_coupled_accumulate(statePool,forcingPool,2)

      if (config_use_GM .or. config_submesoscale_enable) then
         call ocn_reconstruct_eddy_vectors(meshPool)
      end if

      if (trim(config_land_ice_flux_mode) == 'coupled') then
         call mpas_timer_start("si effective density halo")
         call mpas_pool_get_field(statePool, &
                                 'effectiveDensityInLandIce', &
                                  effectiveDensityField, 2)
         call mpas_dmpar_exch_halo_field(effectiveDensityField)
         call mpas_timer_stop("si effective density halo")
      end if
      deallocate(bottomDepthEdge)

      call mpas_timer_stop('si fini')
      call mpas_timer_stop("si timestep")

   end subroutine ocn_time_integrator_si!}}}

!***********************************************************************
!
!  routine ocn_time_integration_si_init
!
!> \brief   Initialize split-implicit time stepping within MPAS-Ocean
!> \author  Mark Petersen
!> \date    September 2011
!
!> \author  Hyun-Gyu Kang (ORNL, for the split-implicit code)
!> \date    September 2019
!> \details
!>  This routine initializes variables required for the split-implicit
!>  method of integrating the ocean model forward in time
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integration_si_init(domain, dt)!{{{
     
      include 'mpif.h'

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: dt !< [in] time step (sec)

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (block_type), pointer :: &
         block ! structure with subdomain data

      type (mpas_pool_type), pointer :: &
         statePool,         &! structure holding state variables
         meshPool,          &! structure holding mesh variables
         diagnosticsPool     ! structure holding diagnostic variables

      integer, pointer :: nVertLevels

      integer ::         &
         iCell, iEdge, k,&! loop indices for cell, edge and vertical
         kmax,           &! index of deepest active edge
         ierr,           &! local error flag
         cell1, cell2     ! neighbor cell indices across edge

      integer :: nCells, nEdges, ihh, imm, iss

      integer, dimension(:), pointer :: nCellsArray, nEdgesArray,     &
         minLevelEdgeBot, maxLevelEdgeTop, minLevelCell, maxLevelCell

      integer, dimension(:,:), pointer :: cellsOnEdge

      real (kind=RKIND) ::       &
         normalThicknessFluxSum, &! vertical sum of thick flux
         layerThicknessSum,      &! vertical sum of layer thickness
         layerThicknessEdge1,    &! layer thickness on edge
         area_mean,              &! RMS of areaCell
         local_num_cells,        &! number of cells for each core
         local_area_sum,         &! local sum of area
         total_area_sum,         &! total sum of area
         tmp1,tmp2

      real (kind=RKIND), dimension(:), pointer :: &
         areaCell,               &! area of each cell
         bottomDepth,            &! bottom depth
         refBottomDepth,         &! reference bottom depth
         normalBarotropicVelocity ! normal barotropic velocity

      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThickness,           &! layer thickness cell center
         normalBaroclinicVelocity, &! normal baroclinic velocity
         normalVelocity             ! normal velocity (total)

      real (kind=RKIND), dimension(:), pointer :: ssh

#ifdef USE_GPU_AWARE_MPI
      type (mpas_communication_list), pointer :: sendList, recvList
      type (mpas_communication_list), pointer :: commListPtr,commListPtr2
      type (mpas_exchange_list), pointer :: exchListPtr
      integer :: iHalo, nHaloLayers, threadNum, nAdded, i, j,iComm,iExch
#endif

#ifndef USE_LAPACK
      call mpas_log_write( &
      'MPAS was not compiled with LAPACK/BLAS. LAPACK required for SI' &
      , MPAS_LOG_CRIT)
#endif

      !*** Set mask for using velocity correction
      if (config_vel_correction) then
         useVelocityCorrection = 1.0_RKIND
      else
         useVelocityCorrection = 0.0_RKIND
      endif

      !*** Determine the time integration type and set associated masks

      select case (trim(config_time_integrator))

      case ('split_implicit')

         call mpas_log_write( &
         '***********************************************************')
         call mpas_log_write( &
         'The split-implicit time integration is configured')
         call mpas_log_write( &
         '***********************************************************')

      case default
         call mpas_log_write('Incorrect choice config_time_integrator',&
                             MPAS_LOG_CRIT)

      end select

      !*** set number of baroclinic iterations on each outer
      !*** time step iteration (number can be different on the
      !*** first and last time step iteration)

      numTSIterations = config_n_ts_iter

      allocate(numClinicIterations(numTSIterations))

      numClinicIterations    = config_n_bcl_iter_mid ! most iterations
      numClinicIterations(1) = config_n_bcl_iter_beg ! first iteration
      numClinicIterations(numTSIterations)=config_n_bcl_iter_end !last

      block => domain % blocklist

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', &
                                                        nVertLevels)

      call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', &
                                                        nCellsArray)
      call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', &
                                                        nEdgesArray)

      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThickness, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocity, 1)
      call mpas_pool_get_array(statePool, 'normalBarotropicVelocity', &
                                           normalBarotropicVelocity, 1)
      call mpas_pool_get_array(statePool, 'normalBaroclinicVelocity', &
                                           normalBaroclinicVelocity, 1)

      call mpas_pool_get_array(meshPool, 'refBottomDepth', &
                                          refBottomDepth)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', &
                                          cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', &
                                          minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', &
                                          maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      call mpas_pool_get_array(statePool, 'ssh', ssh, 2)

      ! Get fields to use 1d-real halo exchange routine
      call mpas_pool_get_field(diagnosticsPool, 'SIvec_zh1', SIvec_zh1_field,1)
      call mpas_pool_get_field(diagnosticsPool, 'SIvec_wh0', SIvec_wh0_field,1)
      call mpas_pool_get_field(diagnosticsPool, 'SIvec_rh1', SIvec_rh1_field,1)

      nCells = nCellsArray(config_num_halos)
      nEdges = nEdgesArray(config_num_halos+1)

      if ( .not. config_do_restart ) then

      ! Compute barotropic velocity at first timestep
      ! This is only done upon start-up.
         if (config_filter_btr_mode) then
            do iCell = 1, nCells
               layerThickness(minLevelCell(iCell),iCell) = &
                  refBottomDepth(minLevelCell(iCell))
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

      endif ! .not. config_do_restart

      nCells = nCellsArray(1)

      ! Compute the root mean square of areaCell 
      ! for the solver tolerance
      local_num_cells = nCells
      call mpas_dmpar_sum_real(domain % dminfo, local_num_cells, &
                               total_num_cells)

      local_area_sum = 0.0_RKIND
      do iCell = 1,nCells
        local_area_sum = local_area_sum + areaCell(iCell)**2.0
      end do

      call mpas_dmpar_sum_real(domain % dminfo,local_area_sum, &
                                               total_area_sum)

      area_mean = dsqrt(total_area_sum / total_num_cells)
      ncpus = domain % dminfo % nprocs
      mean_num_cells = total_num_cells/ncpus

      ! Tolerance for the outer iteration
      tolerance_outer   = 0.01_RKIND * area_mean

      ! Tolerance for the inner iteration
      if ( config_btr_si_partition_match_mode ) then
         ! Tolerance for the partition match mode 
         tolerance_inner  = 1.0e-8_RKIND * area_mean
      else
         tolerance_inner  = config_btr_si_tolerance * area_mean
      endif


      ! Detection of ISMF (Temporariliy implemented. 
      !                    This will be revised in next SI version)
      ! If ISMF is detected, the split-implicit barotropic mode solver
      ! will solve a 'quasi-linear' barotropic system for the more 
      ! stable solver convergence.
      do iCell = 1, nCells
         k = maxLevelCell(iCell)
         zTop(k:nVertLevels,iCell) = -bottomDepth(iCell)      &
                                     +layerThickness(k,iCell)

         do k = maxLevelCell(iCell)-1, minLevelCell(iCell), -1
            zTop(k,iCell) = zTop(k+1,iCell) + layerThickness(k  ,iCell)
         end do

         ! copy zTop(1,iCell) into sea-surface height array
         ssh(iCell) = zTop(minLevelCell(iCell),iCell)
      end do

      tmp1 = minval(ssh)
      call mpas_dmpar_min_real(domain % dminfo, tmp1,tmp2 )

      si_ismf = 1
      if ( tmp2 < -10.d0 ) then
         si_ismf = 0
      endif

      ! Impliciness parameters
      alpha1= 0.5_RKIND
      alpha2= 1.0_RKIND - alpha1

      ! Get time step size to compute coefficients for the SI solver
      dt_si = dt

      ! Determination of nSiLargeIter (the barotropic system 
      !                                large iteration loop)
      !    - Currently, it is set as 2 when numTSIterations is 1.
      !    - Also, nSiLargeIter can be controlled by user if needed.
      !    - Higher nSiLargeIter can make simulations more stable and
      !      accurate, but runtime for the barotropic system will 
      !      increase.
      !    - Do not set larger than 2. This will be improved in next
      !      version.
      nSiLargeIter = config_n_btr_si_large_iter

      if ( numTSIterations == 1 ) then
         nSiLargeIter = 2
      else
         nSiLargeIter = config_n_btr_si_large_iter
         dt_si = dt_si / real(nSiLargeIter)
      endif


      ! Computation of coefficients which will be used in the SI solver
      R1_alpha1s_g_dts = 1.0_RKIND/((alpha1**2.0_RKIND) &
                           * gravity * dt_si**2.0_RKIND)
      R1_alpha1s_g_dt  = 1.0_RKIND/((alpha1**2.0_RKIND) &
                           * gravity * dt_si)
      R1_alpha1_g      = 1.0_RKIND/(gravity*alpha1)

      ! Compute preconditioner for the split-implicit barotropic
      !  mode solver
      call mpas_log_write(' Building a preconditioning matrix: ')
      call mpas_timer_start("preconditioning matrix build")
      call ocn_time_integrator_si_preconditioner(domain, dt)
      call mpas_timer_stop("preconditioning matrix build")

      !--------------------------------------------------------------------

      if (config_use_self_attraction_loading) then
         ! set pgf_sal_on to 1
         pgf_sal_on = 1.0_RKIND
      else
         pgf_sal_on = 0.0_RKIND
      endif

      self_attraction_and_loading_beta = config_self_attraction_and_loading_beta


#ifdef USE_GPU_AWARE_MPI
      !--------------------------------------------------------------------
      ! Initialization for the si_halo_exch routine
      !    - Save sendList and recvList
      !    - Allocate send & recv buffer
      !    - Define dest & source list
      !--------------------------------------------------------------------

      call mpi_barrier(domain%dminfo%comm,ierr)

      threadNum = mpas_threading_get_thread_num()
      if ( threadNum == 0 ) then
         nHaloLayers = size( SIvec_rh1_field % sendList % halos)
         allocate(haloLayers(nHaloLayers))
         do iHalo = 1, nHaloLayers
            haloLayers(iHalo) = iHalo
         end do
      end if

      call mpas_dmpar_build_comm_lists(                                    &
               SIvec_rh1_field % sendList, SIvec_rh1_field % recvList,     &
               haloLayers, SIvec_rh1_field % dimsizes, sendList, recvList )

      !-----------------------------------------------------------------

      i = 0
      nMaxSendList = 0
      commListPtr => sendList
      do while(associated(commListPtr))
         i = i + 1
         if ( nMaxSendList < commListPtr%nList ) then
              nMaxSendList = commListPtr%nList
         endif
         commListPtr => commListPtr % next
      end do
      nSendCommList = i

      i = 0
      nMaxRecvList = 0
      commListPtr => recvList
      do while(associated(commListPtr))
         i = i + 1
         if ( nMaxRecvList < commListPtr%nList ) then
              nMaxRecvList = commListPtr%nList
         endif
         commListPtr => commListPtr % next
      end do
      nRecvCommList = i

      allocate(exchSend(nSendCommList))
      allocate(exchRecv(nRecvCommList))

      !-----------------------------------------------------------------
      ! Exchange list for SEND

      commListPtr => sendList 
      i = 0
      do while(associated(commListPtr))
         i = i + 1

         allocate(exchSend(i)%buffer(commListPtr%nList))
         allocate(exchSend(i)%bufferOffset(nHaloLayers))
         allocate(exchSend(i)%halo(nHaloLayers))
         allocate(exchSend(i)%nExch(nHaloLayers))

         exchSend(i)%procID = commListPtr%procID
         exchSend(i)%nList  = commListPtr%nList
      
         nAdded = 0
         do iHalo = 1,nHaloLayers
            exchSend(i)%bufferOffset(iHalo) = nAdded
            exchListPtr => SIvec_rh1_field % sendList %      &
                           halos(haloLayers(iHalo)) % exchList
            j = 0
            do while (associated(exchListPtr))
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  do k = 1,exchListPtr % nList
                     nAdded = nAdded + 1 
                  end do
               endif
               j = j+1
               exchListPtr => exchListPtr % next
            end do
 
            exchSend(i)%nExch(iHalo) = j

            allocate(exchSend(i)%halo(iHalo)%endPointID(j))
            allocate(exchSend(i)%halo(iHalo)%nList(j))
            allocate(exchSend(i)%halo(iHalo)%exch(j))

            nullify(exchListPtr)
            exchListPtr => SIvec_rh1_field % sendList %      &
                           halos(haloLayers(iHalo)) % exchList
            j = 0
            do while (associated(exchListPtr))
               j = j+1
               exchSend(i)%halo(iHalo)%endPointID(j) = exchListPtr % endPointID
               exchSend(i)%halo(iHalo)%nList(j) = exchListPtr % nList

               allocate(exchSend(i)%halo(iHalo)%exch(j)%destList(exchListPtr % nList))
               allocate(exchSend(i)%halo(iHalo)%exch(j)%srcList(exchListPtr % nList))

               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  do k = 1,exchListPtr % nList
                     exchSend(i)%halo(iHalo)%exch(j)%destList(k) = &
                        exchListPtr % destList(k) + exchSend(i)%bufferOffset(iHalo)
                     exchSend(i)%halo(iHalo)%exch(j)%srcList(k) = &
                        exchListPtr % srcList(k)
                  end do
               endif
               exchListPtr => exchListPtr % next
            end do

         end do ! iHalo
         commListPtr => commListPtr  % next
      end do

      !-----------------------------------------------------------------
      ! Exchange list for RECV

      commListPtr => recvList 
      i = 0
      do while(associated(commListPtr))
         i = i + 1

         allocate(exchRecv(i)%buffer(commListPtr%nList))
         allocate(exchRecv(i)%bufferOffset(nHaloLayers))
         allocate(exchRecv(i)%halo(nHaloLayers))
         allocate(exchRecv(i)%nExch(nHaloLayers))

         exchRecv(i)%procID = commListPtr%procID
         exchRecv(i)%nList  = commListPtr%nList

         nAdded = 0
         do iHalo = 1,nHaloLayers
            exchRecv(i)%bufferOffset(iHalo) = nAdded
            exchListPtr => SIvec_rh1_field % recvList %      &
                           halos(haloLayers(iHalo)) % exchList
            j = 0
            do while (associated(exchListPtr))
               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  do k = 1,exchListPtr % nList
                     nAdded = nAdded + 1 
                  end do
               endif
               j = j+1
               exchListPtr => exchListPtr % next
            end do
 
            exchRecv(i)%nExch(iHalo) = j

            allocate(exchRecv(i)%halo(iHalo)%endPointID(j))
            allocate(exchRecv(i)%halo(iHalo)%nList(j))
            allocate(exchRecv(i)%halo(iHalo)%exch(j))

            nullify(exchListPtr)
            exchListPtr => SIvec_rh1_field % recvList %      &
                           halos(haloLayers(iHalo)) % exchList
            j = 0
            do while (associated(exchListPtr))
               j = j+1
               exchRecv(i)%halo(iHalo)%endPointID(j) = exchListPtr % endPointID
               exchRecv(i)%halo(iHalo)%nList(j) = exchListPtr % nList

               allocate(exchRecv(i)%halo(iHalo)%exch(j)%destList(exchListPtr % nList))
               allocate(exchRecv(i)%halo(iHalo)%exch(j)%srcList(exchListPtr % nList))

               if ( exchListPtr % endPointID == commListPtr % procID ) then
                  do k = 1,exchListPtr % nList
                     exchRecv(i)%halo(iHalo)%exch(j)%destList(k) = &
                        exchListPtr % destList(k)
                     exchRecv(i)%halo(iHalo)%exch(j)%srcList(k) = &
                        exchListPtr % srcList(k) + exchRecv(i)%bufferOffset(iHalo)
                  end do
               endif
               exchListPtr => exchListPtr % next
            end do

         end do ! iHalo
         commListPtr => commListPtr  % next
      end do

      nullify(sendList,recvList) 
      !-----------------------------------------------------------------

      !$acc enter data create(sbuffer,rbuffer)

      !$acc enter data copyin(exchRecv)
      do iComm = 1,nRecvCommList
         !$acc enter data copyin(exchRecv(iComm)%buffer, &
         !$acc                   exchRecv(iComm)%nExch,  &
         !$acc                   exchRecv(iComm)%halo,   &
         !$acc                   exchRecv(iComm)%procID)
         do iHalo = 1,config_num_halos
            !$acc enter data copyin(exchRecv(iComm)%halo(iHalo)%nList,      &
            !$acc                   exchRecv(iComm)%halo(iHalo)%endPointID, &
            !$acc                   exchRecv(iComm)%halo(iHalo)%exch )
            do iExch = 1,exchRecv(iComm)%nExch(iHalo)
               !$acc enter data copyin(                                   &
               !$acc    exchRecv(iComm)%halo(iHalo)%exch(iExch)%destList, &
               !$acc    exchRecv(iComm)%halo(iHalo)%exch(iExch)%srcList )
            end do
         end do
      end do

      !$acc enter data copyin(exchSend)
      do iComm = 1,nSendCommList
         !$acc enter data copyin(exchSend(iComm)%buffer, &
         !$acc                   exchSend(iComm)%nExch,  &
         !$acc                   exchSend(iComm)%halo,   &
         !$acc                   exchSend(iComm)%procID)
         do iHalo = 1,config_num_halos
            !$acc enter data copyin(exchSend(iComm)%halo(iHalo)%nList,      &
            !$acc                   exchSend(iComm)%halo(iHalo)%endPointID, &
            !$acc                   exchSend(iComm)%halo(iHalo)%exch)
            do iExch = 1,exchSend(iComm)%nExch(iHalo)
               !$acc enter data copyin(                                   &
               !$acc    exchSend(iComm)%halo(iHalo)%exch(iExch)%destList, &
               !$acc    exchSend(iComm)%halo(iHalo)%exch(iExch)%srcList)
            end do
         end do
      end do
#endif

   end subroutine ocn_time_integration_si_init!}}}

!***********************************************************************
!
!  routine ocn_time_integration_si_preconditioner
!
!> \brief   Construct a Block-Jacobi preconditioning matrix
!> \author  Hyun-Gyu Kang (Oak Ridge National Laboratory)
!> \date    September 2019
!> \details
!>  This routine constructs preconditioners for the
!>  split-implicit time stepper.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_si_preconditioner(domain, dt)!{{{

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
         meshPool ! structure holding mesh variables

      integer,dimension(:),pointer :: &
         globalCellId ! global index of each cell

      integer ::         &
         nCells, nEdges, &! number of cells or edges (excl halos) 
         iCell, iEdge,   &! loop indices for cell and edge
         cell1, cell2,   &! neighbor cell indices across edge
         nCellsHalo1st,  &! number of cells within 1st halo layer
         nCellsHalo2nd,  &! number of cells within 2st halo layer
         nPrecMatPacked   ! number of matrix elements 
                          !    of an uppder digonal matrix

      integer ::  i, j, info, itmp1

      real (kind=RKIND) :: thicknessSum, fluxAx, temp1

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      block => domain % blocklist
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_array(meshPool, 'indexToCellID', globalCellId)

      nCells   = nCellsOwned
      nEdges   = nEdgesOwned
      nCellsHalo1st = nCellsHalo(1)
      nCellsHalo2nd = nCellsHalo(2)

      si_algorithm = 'sbicg' ! Default iterative solver algorithm

      call mpas_log_write('   config_btr_si_preconditioner: ' &
                           // trim(config_btr_si_preconditioner) )


      if ( config_btr_si_partition_match_mode ) then

         call mpas_log_write( &
         '       Partition-match (bit-for-bit) mode is turned on.')
         call mpas_log_write( &
         '       The preconditioner is configured as jacobi.')
         config_btr_si_preconditioner = 'jacobi'
         call mpas_log_write( &
         '       The bit-for-bit allreduce is enabled.')
         call mpas_log_write( &
         '       The solver algorithm is configured as scg.')
         call mpas_log_write( &
         '       (Default is sbicg.)')
         si_algorithm = 'scg'

         ! Allocate temporary arrays for reproducible summation
         allocate( globalReprodSum2fld1(nCellsOwned,2),  &
                   globalReprodSum2fld2(nCellsOwned,2),  &
                   globalReprodSum3fld1(nCellsOwned,3),  &
                   globalReprodSum3fld2(nCellsOwned,3) )

      else

         if ( trim(config_btr_si_preconditioner) == 'ras' .or. &
              trim(config_btr_si_preconditioner) == 'block_jacobi') then
#if defined USE_MAGMA
            if ( int(mean_num_cells) > 4800 ) then
               call mpas_log_write( &
               '      nCells per core is larger than 4800 (MAGMA enabled).')
#elif defined USE_CUBLAS
            if ( int(mean_num_cells) > 4800 ) then
               call mpas_log_write( &
               '      nCells per core is larger than 4800 (CUBLAS enabled).')
#elif defined MPAS_OPENACC
            if ( int(mean_num_cells) > 2400 ) then
               call mpas_log_write( &
               '      nCells per core is larger than 2400 (no MAGMA or CUBLAS).')
#else
            if ( int(mean_num_cells) > 1200 ) then
               call mpas_log_write( &
               '      nCells per core is larger than 1200.')
#endif
               call mpas_log_write( &
               '      Because of memory and computational efficiency,')
               call mpas_log_write( &
               '      the preconditioner is configured as jacobi.')
               config_btr_si_preconditioner = 'jacobi'
               call mpas_log_write( &
               '      The solver algorithm is configured as scg.')
               call mpas_log_write( &
               '      (Default is sbicg.)')
               si_algorithm = 'scg'
            endif
         endif

      endif

      ! Restricted Additive Schwarz preconditioner --------------------!
      if ( trim(config_btr_si_preconditioner) == 'ras' .or. &
           trim(config_btr_si_preconditioner) == 'block_jacobi' ) then

! MAGMA IF ----!
#if defined(USE_MAGMA) || defined(USE_CUBLAS)
! MAGMA IF ----!

         if ( trim(config_btr_si_preconditioner) == 'ras' ) then
            nPrecVec = nCellsHalo2nd ! length of preconditioning vector
         else
            nPrecVec = nCells        ! length of preconditioning vector
         endif

         allocate(prec_ivmat(1:nPrecVec,1:nPrecVec))
                  prec_ivmat(:,:) = 0.0_RKIND

         do iCell = 1, nPrecVec

            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)

               ! method 1, matches method 0 without pbcs,
               ! works with pbcs.
               thicknessSum = min(bottomDepth(cell1), &
                                  bottomDepth(cell2))
               fluxAx = edgeSignOnCell(i,iCell) * dvEdge(iEdge) &
                      * thicknessSum / dcEdge(iEdge)

               if ( globalCellId(cell1) > 0 .and. &
                    globalCellId(cell1) < total_num_cells+1 ) then

                  if ( cell1 >= iCell .and. cell1 <= nPrecVec) then
                     prec_ivmat(iCell,cell1) = &
                     prec_ivmat(iCell,cell1) + fluxAx
                  endif
               endif

               if ( globalCellId(cell2) > 0 .and. &
                    globalCellId(cell2) < total_num_cells+1 ) then

                  if ( cell2 >= iCell .and. cell2 <= nPrecVec) then
                     prec_ivmat(iCell,cell2) = &
                     prec_ivmat(iCell,cell2) - fluxAx
                  endif
               endif
            end do ! i

            prec_ivmat(iCell,iCell) = &
            prec_ivmat(iCell,iCell)   &
               - (4.0_RKIND/(gravity*dt**2.0)) * areaCell(iCell)
         end do ! iCell

         ! Inversion of the preconditioning matrix
         prec_ivmat = -prec_ivmat

         ! DPOTRF computes the Cholesky factorizaton of a symmetric positive
         ! defnitie matrix.
#ifdef USE_LAPACK
         call DPOTRF('U',nPrecVec,prec_ivmat,nPrecVec,info)
#endif

         if (info.lt.0) stop 'ERROR in inverse routine: Matrix is numerically singular!'
         if (info.gt.0) stop 'ERROR in inverse routine: Matrix is not positive definite!'

         ! DPOTRI computes the inverse of a symmetric positive definite matrix,
         ! using the Cholesky factorization computed by DPOTRF
#ifdef USE_LAPACK
         call DPOTRI('U',nPrecVec,prec_ivmat,nPrecVec,info)
#endif

         prec_ivmat = -prec_ivmat ! (-) sign to restore a matrix sign

         if (info.ne.0) stop 'ERROR in inverse routine: Matrix inversion failed!'

! MAGMA IF ----!
#else
! MAGMA IF ----!

         if ( trim(config_btr_si_preconditioner) == 'ras' ) then
            nPrecVec = nCellsHalo2nd ! length of preconditioning vector
         else
            nPrecVec = nCells        ! length of preconditioning vector
         endif

         nPrecMatPacked = (nPrecVec*(nPrecVec+1))/2
                                  ! Packed size of preconditionig matrix

         allocate(prec_ivmat(1:nPrecMatPacked,1))
                  prec_ivmat(:,:) = 0.0_RKIND

         do iCell = 1, nPrecVec

            do i = 1, nEdgesOnCell(iCell)   
               iEdge = edgesOnCell(i, iCell)
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)

               ! method 1, matches method 0 without pbcs,
               ! works with pbcs.
               thicknessSum = min(bottomDepth(cell1), &
                                  bottomDepth(cell2))
               fluxAx = edgeSignOnCell(i,iCell) * dvEdge(iEdge) &
                      * thicknessSum / dcEdge(iEdge)

               if ( globalCellId(cell1) > 0 .and. &
                    globalCellId(cell1) < total_num_cells+1 ) then

                  if ( cell1 >= iCell .and. cell1 <= nPrecVec) then
                     prec_ivmat(iCell+((cell1-1)*cell1)/2,1) = &
                     prec_ivmat(iCell+((cell1-1)*cell1)/2,1) + fluxAx
                  endif
               endif

               if ( globalCellId(cell2) > 0 .and. &
                    globalCellId(cell2) < total_num_cells+1 ) then

                  if ( cell2 >= iCell .and. cell2 <= nPrecVec) then
                     prec_ivmat(iCell+((cell2-1)*cell2)/2,1) = &
                     prec_ivmat(iCell+((cell2-1)*cell2)/2,1) - fluxAx
                  endif
               endif
            end do ! i

            prec_ivmat(iCell+((iCell-1)*iCell)/2,1) = &
            prec_ivmat(iCell+((iCell-1)*iCell)/2,1)   &
               - (4.0_RKIND/(gravity*dt**2.0)) * areaCell(iCell)
         end do ! iCell

         ! Inversiion
           ! 1. Cholesky factorization of a real symmetric
           !    positive definite matirx A
         prec_ivmat(:,1) = -prec_ivmat(:,1)
#ifdef USE_LAPACK
         call DPPTRF('U',nPrecVec,prec_ivmat(:,1),info)
#endif
           ! 2. Inversion of a real symmetric positive definite
           !    matrix A using the Cholesky factorization
#ifdef USE_LAPACK
         call DPPTRI('U',nPrecVec,prec_ivmat(:,1),info)
#endif
         prec_ivmat(:,1) = -prec_ivmat(:,1)

! MAGMA IF ----!
#endif
! MAGMA IF ----!

      ! Jacobi preconditioner -----------------------------------------!
      else if ( trim(config_btr_si_preconditioner) == 'jacobi' ) then

         nPrecVec = nCells ! length of preconditioning vector

         allocate(prec_ivmat(1:nPrecVec,1))
                  prec_ivmat(:,1) = 0.0_RKIND

         do iCell = 1, nPrecVec

            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               cell1 = cellsOnEdge(1, iEdge)
               cell2 = cellsOnEdge(2, iEdge)

               ! method 1, matches method 0 without pbcs,
               ! works with pbcs.
               thicknessSum =  min(bottomDepth(cell1), &
                                   bottomDepth(cell2))

               fluxAx = edgeSignOnCell(i,iCell) * dvEdge(iEdge) &
                      * thicknessSum / dcEdge(iEdge)

               if (cell1 == iCell) then
                  prec_ivmat(iCell,1) = prec_ivmat(iCell,1) + fluxAx
                                                        ! reversed sign
               elseif ( cell2 == iCell) then
                  prec_ivmat(iCell,1) = prec_ivmat(iCell,1) - fluxAx
                                                        ! reversed sign
               endif
            end do ! i

            temp1 = prec_ivmat(iCell,1) &
                  - (4.0_RKIND/(gravity*dt**2.0)) * areaCell(iCell)

            prec_ivmat(iCell,1) = 1.0_RKIND / temp1

         end do ! iCell

      ! No preconditioner ---------------------------------------------!
      else if ( trim(config_btr_si_preconditioner) == 'none' ) then

         nPrecVec = nCells ! length of preconditioning vector

         allocate(prec_ivmat(1,1))
         prec_ivmat(:,1) = 1.0_RKIND

      else

         call mpas_log_write( &
         'Incorrect choice for config_btr_si_preconditioner: ' &
         // trim(config_btr_si_preconditioner) //              &
         '   choices are: ras, block_jacobi, jacobi, none',    &
         MPAS_LOG_CRIT)

      endif ! config_btr_si_preconditioner


#ifdef MPAS_OPENACC

#ifdef USE_MAGMA
      ! Allocate and setup the preconditioning matrix and vectors
      !    on GPU. Construct the preconditioing matrix on CPU
      !    first, then copy to GPU.

      call mpas_log_write( "MAGMA --- malloc GPU memory prec mat" )
      info = magmaf_dmalloc( dprec_ivmat, nPrecVec*nPrecVec)
      info = magmaf_dmalloc( dxvec, nPrecVec )
      info = magmaf_dmalloc( dyvec, nPrecVec )

      if (dprec_ivmat == 0 .or. dxvec == 0 .or. dyvec == 0) then
          call mpas_log_write( "MAGMA dmalloc failed", MPAS_LOG_CRIT)
      endif

      call mpas_log_write( "MAGMA --- set matrices" )
      call magmaf_queue_create( 0, queue )
      call magmaf_dsetmatrix( nPrecVec, nPrecVec, prec_ivmat, nPrecVec, &
                              dprec_ivmat, nPrecVec, queue )

      !$acc parallel loop present(SIvec_r0)
      do iCell = 1,nPrecVec
          SIvec_r0(iCell) = float(iCell)
      end do

      !$acc host_data use_device(SIvec_r0)
      dxvec = loc(SIvec_r0(1:nPrecVec))
      !$acc end host_data

      call magmaf_dsymv('U', nPrecVec, 1.0d0, dprec_ivmat, &
                             nPrecVec, dxvec, 1, 0.0d0,    &
                                       dyvec, 1, queue)
    
      call magmaf_dgetvector( nPrecVec,dyvec, 1, &
                              SIvec_rh0(1:nPrecVec), 1, queue )
#endif

      !$acc enter data copyin(prec_ivmat,                            &
      !$acc    SIcst_allreduce_local2,SIcst_allreduce_global2,       &
      !$acc    SIcst_allreduce_local3,SIcst_allreduce_global3,       &
      !$acc    SIcst_allreduce_local9,SIcst_allreduce_global9,       &
      !$acc    SIcst_alpha0,SIcst_alpha1,SIcst_beta0,SIcst_beta1,    &
      !$acc    SIcst_omega0,SIcst_gamma0,SIcst_gamma1,resid,         &
      !$acc    SIcst_r0rh0,SIcst_w0rh0,SIcst_r0r0 ,SIcst_r00r0,      &
      !$acc    SIcst_r00s0,SIcst_r00z0,SIcst_q0y0 ,SIcst_y0y0 ,      &
      !$acc    SIcst_r00y0,SIcst_r00t0,SIcst_r00v0,SIcst_q0q0 ,      &
      !$acc    SIcst_r00w0,SIcst_r00q0,SIcst_rho0,SIcst_rho1)

      if ( config_btr_si_partition_match_mode ) then
      !$acc enter data copyin(                                       &
      !$acc               globalReprodSum2fld1,globalReprodSum2fld2, &
      !$acc               globalReprodSum3fld1,globalReprodSum3fld2 )
      endif

#ifdef USE_CUBLAS
! Dummy preconditioning
      call si_precond(SIvec_r0,SIvec_rh0)
#endif

#endif

      !$acc serial present(                                          &
      !$acc    SIcst_allreduce_local2,SIcst_allreduce_global2,       &
      !$acc    SIcst_allreduce_local3,SIcst_allreduce_global3,       &
      !$acc    SIcst_allreduce_local9,SIcst_allreduce_global9,       &
      !$acc    SIcst_alpha0,SIcst_alpha1,SIcst_beta0,SIcst_beta1,    &
      !$acc    SIcst_omega0,SIcst_gamma0,SIcst_gamma1,resid,         &
      !$acc    SIcst_r0rh0,SIcst_w0rh0,SIcst_r0r0 ,SIcst_r00r0,      &
      !$acc    SIcst_r00s0,SIcst_r00z0,SIcst_q0y0 ,SIcst_y0y0 ,      &
      !$acc    SIcst_r00y0,SIcst_r00t0,SIcst_r00v0,SIcst_q0q0 ,      &
      !$acc    SIcst_r00w0,SIcst_r00q0,SIcst_rho0,SIcst_rho1)
      SIcst_allreduce_local2  = 0.0_RKIND
      SIcst_allreduce_global2 = 0.0_RKIND
      SIcst_allreduce_local3  = 0.0_RKIND
      SIcst_allreduce_global3 = 0.0_RKIND
      SIcst_allreduce_local9  = 0.0_RKIND
      SIcst_allreduce_global9 = 0.0_RKIND

      SIcst_alpha0 = 0.0_RKIND ; SIcst_alpha1 = 0.0_RKIND
      SIcst_beta0  = 0.0_RKIND ; SIcst_beta1  = 0.0_RKIND
      SIcst_gamma0 = 0.0_RKIND ; SIcst_gamma1 = 0.0_RKIND
      SIcst_r0rh0  = 0.0_RKIND ; SIcst_w0rh0  = 0.0_RKIND
      SIcst_r0r0   = 0.0_RKIND ; SIcst_r00r0  = 0.0_RKIND
      SIcst_r00s0  = 0.0_RKIND ; SIcst_r00z0  = 0.0_RKIND
      SIcst_q0y0   = 0.0_RKIND ; SIcst_y0y0   = 0.0_RKIND
      SIcst_r00y0  = 0.0_RKIND ; SIcst_r00t0  = 0.0_RKIND
      SIcst_r00v0  = 0.0_RKIND ; SIcst_q0q0   = 0.0_RKIND
      SIcst_r00w0  = 0.0_RKIND ; SIcst_r00q0  = 0.0_RKIND
      SIcst_rho0   = 0.0_RKIND ; SIcst_rho1   = 0.0_RKIND
      !$acc end serial

   end subroutine ocn_time_integrator_si_preconditioner !}}}

!***********************************************************************
!
!  routine si_matvec_mul
!
!> \brief   Matrix-vector multiplication for iterative solvers
!
!> \author  Hyun-Gyu Kang
!
!> \details
!   This routine performs a matrix-vector multiplication for iterative
!   solvers.
!
!-----------------------------------------------------------------------

   subroutine si_matvec_mul(SIvec_in,bottomDepthEdge,SIvec_out,ssh) !{{{

      real(kind=RKIND),dimension(:),intent(in)  :: &
          SIvec_in,      &! Input SI vector
          ssh,           &! Intermediate ssh
          bottomDepthEdge ! Bottom Depth at Edge

      real(kind=RKIND),dimension(:),intent(out) :: &
          SIvec_out  ! Output SI vector

      real(kind=RKIND) :: sshTendAx, sshEdge, thicknessSum, sshDiff
      real(kind=RKIND) :: fluxAx, sshArea

      integer :: iEdge, cell1, cell2, i, iCell
 
      call mpas_timer_start("si matvec_mul")

#ifdef MPAS_OPENACC
      !$acc parallel loop gang async                          &
      !$acc    present(SIvec_in,SIvec_out,ssh,                &
      !$acc       nEdgesOnCell,edgesOnCell,cellsOnEdge,       &
      !$acc       bottomDepthEdge, dcEdge,edgeSignOnCell,     &
      !$acc       dvEdge,areaCell,maxLevelEdgeTop)            &
      !$acc    private(sshTendAx,i,iEdge,cell1,cell2,sshEdge, &
      !$acc       thicknessSum,fluxAx,sshArea,sshDiff)
#else
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(sshTendAx,iEdge,cell1,cell2,sshEdge, &
      !$omp         thicknessSum,sshDiff,fluxAx,sshArea )
#endif
      do iCell = 1, nPrecVec
         sshTendAx = 0.0_RKIND

         !$acc loop vector reduction(+:sshTendAx)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (maxLevelEdgeTop(iEdge).eq.0) cycle

            cell1 = cellsOnEdge(1, iEdge)
            cell2 = cellsOnEdge(2, iEdge)

            ! Interpolation sshEdge
            sshEdge = 0.5_RKIND * (ssh(cell1) + ssh(cell2))

            ! method 1, matches method 0 without pbcs,
            ! works with pbcs.
            thicknessSum = si_ismf * sshEdge + bottomDepthEdge(iEdge)

            ! nabla (ssh^0)
            sshDiff = (SIvec_in(cell2)-SIvec_in(cell1)) / dcEdge(iEdge)

            fluxAx = thicknessSum * sshDiff

            sshTendAx = sshTendAx + edgeSignOnCell(i, iCell) &
                      * fluxAx * dvEdge(iEdge)
         end do ! i

         sshArea = R1_alpha1s_g_dts * SIvec_in(iCell) * areaCell(iCell)

         SIvec_out(iCell) = -sshArea - sshTendAx

      end do ! iCell
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      call mpas_timer_stop("si matvec_mul")
      
   end subroutine si_matvec_mul

!***********************************************************************
!
!  routine si_precond
!
!> \brief   Preconditioning for iterative solvers
!
!> \author  Hyun-Gyu Kang
!
!> \details
!   This routine performs preconditioing for iterative solvers.
!
!-----------------------------------------------------------------------

   subroutine si_precond(SIvec_in,SIvec_out)

      real(kind=RKIND), dimension(:), intent(in)  :: SIvec_in
                                        ! Input SI vector

      real(kind=RKIND), dimension(:), intent(out) :: SIvec_out
                                        ! SI vector to be preconditioned

      integer :: iCell

      call mpas_timer_start("si preconditioning")

      ! Preconditioning --------------------------------------------!
      if ( trim(config_btr_si_preconditioner) == 'ras' .or. &
           trim(config_btr_si_preconditioner) == 'block_jacobi' ) then
         ! Domain decomposition preconditioning 
         !    : Use BLAS for the symmetric
         !      matrix-vector multiplication

!#######################################################################
#if defined(USE_MAGMA)

         ! MAGMA
         !$acc host_data use_device(SIvec_in)
         dxvec = loc(SIvec_in(1:nPrecVec))
         !$acc end host_data

         call magmaf_dsymv('U', nPrecVec, 1.0_RKIND, dprec_ivmat, &
                                nPrecVec, dxvec, 1, 0.0_RKIND,    &
                                          dyvec, 1, queue)

#ifndef USE_GPU_AWARE_MPI
         ! GPU -> CPU copy the out vector for halo exch
         call magmaf_dgetvector( nPrecVec,dyvec, 1, &
                                 SIvec_out(1:nPrecVec), 1, queue)
#else                                 
         ! GPU -> GPU copy the out vector for halo exch
         !$acc host_data use_device(SIvec_out)
         call magmaf_dgetvector( nPrecVec,dyvec, 1, &
                                 SIvec_out(1:nPrecVec), 1, queue)
         !$acc end host_data
#endif                                 

!#######################################################################
#elif defined(USE_CUBLAS)

         ! CUBLAS
         !$acc host_data use_device(prec_ivmat,SIvec_in,SIvec_out)
         call cublasDsymv('U', nPrecVec, 1.0_RKIND,prec_ivmat(:,:), &
                               nPrecVec, SIvec_in(1:nPrecVec),      &
                                      1, 0.0_RKIND,                 &
                                         SIvec_out(1:nPrecVec), 1)
         !$acc end host_data
#ifndef USE_GPU_AWARE_MPI
         !$acc update host(SIvec_out)
#endif

!#######################################################################
#elif defined(USE_LAPACK)

         ! LAPACK on CPU 
         !    (also if no MAGMA & CUBLAS for GPU ; Data staging)
         !$acc update host(SIvec_in)
         call DSPMV('U', nPrecVec, 1.0_RKIND, prec_ivmat(:,1), &
                         SIvec_in(1:nPrecVec) , 1, 0.0_RKIND,  &
                         SIvec_out(1:nPrecVec), 1)
#endif

      elseif ( trim(config_btr_si_preconditioner) == 'jacobi' ) then
         ! Jacobi preconditioning
#ifdef MPAS_OPENACC
         !$acc parallel loop gang                       &
         !$acc    present(SIvec_out,SIvec_in,prec_ivmat)
#else
         !$omp parallel
         !$omp do schedule(runtime)
#endif
         do iCell = 1,nPrecVec
            SIvec_out(iCell) = SIvec_in(iCell) &
                             * prec_ivmat(iCell,1)
         end do
#if defined(MPAS_OPENACC) && !defined(USE_GPU_AWARE_MPI)
         !$acc update host(SIvec_out)
#else
         !$omp end do
         !$omp end parallel
#endif

      elseif ( trim(config_btr_si_preconditioner) == 'none' ) then

         ! No preconditioning
#ifdef MPAS_OPENACC
         !$acc parallel loop gang present(SIvec_out,SIvec_in)
#else
         !$omp parallel
         !$omp do schedule(runtime)
#endif
         do iCell = 1,nPrecVec
            SIvec_out(iCell) = SIvec_in(iCell)
         end do
#if defined(MPAS_OPENACC) && !defined(USE_GPU_AWARE_MPI)
         !$acc update host(SIvec_out)
#else
         !$omp end do
         !$omp end parallel
#endif

      end if

      call mpas_timer_stop("si preconditioning")

   end subroutine si_precond !}}}

!***********************************************************************
!
!  routine si_global_reduction
!
!> \brief   Global reduction for iterative solvers
!
!> \author  Hyun-Gyu Kang
!
!> \details
!   This routine performs global reduction for iterative solvers.
!   It uses a reproducible global summation when
!   'config_btr_si_partition_match_mode' is turned on.
!
!-----------------------------------------------------------------------

   subroutine si_global_reduction(localProd,globalProd,       & !{{{
                             reprodSum1,reprodSum2,           &
                             si_algorithm,nfld,domain,        &
                             SIvec_1,SIvec_2,SIvec_3,SIvec_4, &
                             SIvec_5,SIvec_6,SIvec_7)  

      include 'mpif.h'

      type (domain_type), intent(in) :: &
         domain  !< [inout] model state to advance forward

      real(kind=RKIND), dimension(:), intent(in) :: &
         SIvec_1,SIvec_2 ! SI vectors
      real(kind=RKIND), dimension(:), intent(in), optional :: &
         SIvec_3,SIvec_4,SIvec_5,SIvec_6,SIvec_7 ! SI vectors

      character(*), intent(in) :: si_algorithm ! Solver algorithm

      integer, intent(in) :: nfld ! Number of reductions

      real(kind=RKIND), dimension(:,:), intent(out) :: &
         reprodSum1,reprodSum2 ! Array for reprod. global summation

      real(kind=RKIND), dimension(:), intent(out) :: &
         localProd,globalProd ! Array for global summation

      integer :: iCell,ierr,threadNum

      call mpas_timer_start("si global reduction")

      if ( si_algorithm == 'scg' .and. nfld == 2 ) then !!!
   
         if ( config_btr_si_partition_match_mode ) then
    
#ifdef MPAS_OPENACC
            !$acc parallel loop                      &
            !$acc    present(reprodSum1,reprodSum2,  &
            !$acc            SIvec_1,SIvec_2,SIvec_3)
#else
            !$omp parallel
            !$omp do schedule(runtime)
#endif
            ! Reproducible sum of multiple fields over products
            do iCell = 1,nCellsOwned
               reprodSum1(iCell,1) = SIvec_1(iCell) !r0
               reprodSum1(iCell,2) = SIvec_2(iCell) !w0
    
               reprodSum2(iCell,1) = SIvec_3(iCell) !rh0
               reprodSum2(iCell,2) = SIvec_3(iCell) !rh0
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
        
            !$acc update host(reprodSum1,reprodSum2)
            globalProd(:) =                   &
               mpas_global_sum_nfld(reprodSum1, &
                                    reprodSum2, &
                                    domain%dminfo%comm)
            !$acc update device(globalProd)

         else
    
            !$acc parallel loop gang                     &
            !$acc    present(SIcst_r0rh0,SIcst_w0rh0,    &
            !$acc            SIvec_1,SIvec_2,SIvec_3)    &
            !$acc    reduction(+:SIcst_r0rh0,SIcst_w0rh0)
            do iCell = 1, nCellsOwned
               SIcst_r0rh0 = SIcst_r0rh0 + SIvec_1(iCell) &
                                         * SIvec_3(iCell)
               SIcst_w0rh0 = SIcst_w0rh0 + SIvec_2(iCell) &
                                         * SIvec_3(iCell)
            end do ! iCell

            !$acc serial present(SIcst_r0rh0,SIcst_w0rh0,localProd)
            localProd(1) = SIcst_r0rh0
            localProd(2) = SIcst_w0rh0
            SIcst_r0rh0 = 0.0_RKIND
            SIcst_w0rh0 = 0.0_RKIND
            !$acc end serial

            ! Global sum across CPUs
#ifdef USE_GPU_AWARE_MPI
            !$acc host_data use_device(localProd,globalProd)
#else
            !$acc update host(localProd)
#endif
            call mpas_dmpar_sum_real_array(domain % dminfo, nfld, &
                                           localProd,globalProd)
#ifdef USE_GPU_AWARE_MPI
            !$acc end host_data
#else
            !$acc update device(globalProd)
#endif
         endif
   
      elseif ( si_algorithm == 'scg' .and. nfld == 3 ) then !!!
   
         if ( config_btr_si_partition_match_mode ) then
   
            ! Reproducible sum of multiple fields over products
#ifdef MPAS_OPENACC
            !$acc parallel loop gang                 &
            !$acc    present(reprodSum1,reprodSum2,  &
            !$acc            SIvec_1,SIvec_2,SIvec_3)
#else
            !$omp parallel
            !$omp do schedule(runtime)
#endif
            do iCell = 1,nCellsOwned
               reprodSum1(iCell,1) = SIvec_1(iCell) !r1
               reprodSum1(iCell,2) = SIvec_2(iCell) !w1
               reprodSum1(iCell,3) = SIvec_1(iCell) !r1
      
               reprodSum2(iCell,1) = SIvec_3(iCell) !rh1
               reprodSum2(iCell,2) = SIvec_3(iCell) !rh1
               reprodSum2(iCell,3) = SIvec_1(iCell)  !r1
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
      
            !$acc update host(reprodSum1,reprodSum2)
            globalProd(:) =                   &
               mpas_global_sum_nfld(reprodSum1, &
                                    reprodSum2, &
                                    domain%dminfo%comm)
            !$acc update device(globalProd)

         else

            !$acc parallel loop present(SIcst_r0rh0,SIcst_w0rh0,     &
            !$acc    SIcst_r0r0,SIvec_1,SIvec_2,SIvec_3)             &
            !$acc    reduction(+:SIcst_r0rh0,SIcst_w0rh0,SIcst_r0r0)
            do iCell = 1,nCellsOwned
               SIcst_r0rh0 = SIcst_r0rh0 + SIvec_1(iCell) &
                                         * SIvec_3(iCell) ! s1

               SIcst_w0rh0 = SIcst_w0rh0 + SIvec_2(iCell) &
                                         * SIvec_3(iCell) ! z1

               SIcst_r0r0  = SIcst_r0r0  + SIvec_1(iCell) &
                                         * SIvec_1(iCell) ! z1
            end do

            !$acc serial present(SIcst_r0rh0,SIcst_w0rh0,SIcst_r0r0, &
            !$acc                localProd)
            localProd(1) = SIcst_r0rh0
            localProd(2) = SIcst_w0rh0
            localProd(3) = SIcst_r0r0
            SIcst_r0rh0 = 0.0_RKIND
            SIcst_w0rh0 = 0.0_RKIND
            SIcst_r0r0  = 0.0_RKIND
            !$acc end serial

            ! Global sum across CPUs
#ifdef USE_GPU_AWARE_MPI
            !$acc host_data use_device(localProd,globalProd)
#else
            !$acc update host(localProd)
#endif
            call MPI_Allreduce(localProd,globalProd,nfld,MPI_REAL8, &
                               MPI_SUM,domain%dminfo%comm,ierr)
#ifdef USE_GPU_AWARE_MPI
            !$acc end host_data
#else
            !$acc update device(globalProd)
#endif

         endif
   
      elseif ( si_algorithm == 'sbicg' .and. nfld == 2 ) then !!!

         if ( config_btr_si_partition_match_mode ) then

            ! Reproducible sum of multiple fields over products
#ifdef MPAS_OPENACC
            !$acc parallel loop                      &
            !$acc    present(reprodSum1,reprodSum2,  &
            !$acc            SIvec_1,SIvec_2,SIvec_3)
#else
            !$omp parallel
            !$omp do schedule(runtime)
#endif
            do iCell = 1,nCellsOwned
               reprodSum1(iCell,1) = SIvec_1(iCell) !r00
               reprodSum1(iCell,2) = SIvec_1(iCell) !r00

               reprodSum2(iCell,1) = SIvec_2(iCell) !r0
               reprodSum2(iCell,2) = SIvec_3(iCell) !w0
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif

            !$acc update host(reprodSum1,reprodSum2)
            globalProd(:) =                                &
               mpas_global_sum_nfld(reprodSum1,reprodSum2, &
                                    domain%dminfo%comm)
            !$acc update device(globalProd)

         else

            !$acc parallel loop present(SIcst_r00r0,SIcst_r00w0, &
            !$acc    SIvec_1,SIvec_2,SIvec_3)                    &
            !$acc    reduction(+:SIcst_r00r0,SIcst_r00w0)
            do iCell = 1, nCellsOwned
               SIcst_r00r0 = SIcst_r00r0 + SIvec_1(iCell) &
                                         * SIvec_2(iCell)
               SIcst_r00w0 = SIcst_r00w0 + SIvec_1(iCell) &
                                         * SIvec_3(iCell)
            end do ! iCell

            !$acc serial present(SIcst_r00r0,SIcst_r00w0,localProd)   
            localProd(1) = SIcst_r00r0
            localProd(2) = SIcst_r00w0
            SIcst_r00r0 = 0.0_RKIND
            SIcst_r00w0 = 0.0_RKIND
            !$acc end serial

            ! Global sum across CPUs
            threadNum = mpas_threading_get_thread_num()
            if ( threadNum == 0 ) then
#ifdef USE_GPU_AWARE_MPI
            !$acc host_data use_device(localProd,globalProd)
#else
            !$acc update host(localProd)
#endif
            call MPI_Allreduce(localProd,globalProd,nfld,MPI_REAL8, &
                               MPI_SUM,domain%dminfo%comm,ierr)
#ifdef USE_GPU_AWARE_MPI
            !$acc end host_data
#else
            !$acc update device(globalProd)
#endif
            endif ! thread

         endif

      elseif ( si_algorithm == 'sbicg' .and. nfld == 9 ) then !!!

         if ( config_btr_si_partition_match_mode ) then
 
            ! Reproducible sum of multiple fields over products
             
#ifdef MPAS_OPENACC
            !$acc parallel loop                                  &
            !$acc    present(reprodSum1,reprodSum2,              &
            !$acc       SIvec_1,SIvec_2,SIvec_3,SIvec_4,SIvec_5, &
            !$acc       SIvec_6,SIvec_7)
#else
            !$omp parallel
            !$omp do schedule(runtime)
#endif
            do iCell = 1,nCellsOwned
               globalReprodSum9fld1(iCell,1) = SIvec_1(iCell) !r00
               globalReprodSum9fld1(iCell,2) = SIvec_1(iCell) !r00
               globalReprodSum9fld1(iCell,3) = SIvec_2(iCell) !q0
               globalReprodSum9fld1(iCell,4) = SIvec_3(iCell) !y0
               globalReprodSum9fld1(iCell,5) = SIvec_1(iCell) !r00
               globalReprodSum9fld1(iCell,6) = SIvec_1(iCell) !r00
               globalReprodSum9fld1(iCell,7) = SIvec_1(iCell) !r00
               globalReprodSum9fld1(iCell,8) = SIvec_1(iCell) !r00
               globalReprodSum9fld1(iCell,9) = SIvec_2(iCell) !q0
    
               globalReprodSum9fld2(iCell,1) = SIvec_4(iCell) !s1
               globalReprodSum9fld2(iCell,2) = SIvec_5(iCell) !z1
               globalReprodSum9fld2(iCell,3) = SIvec_3(iCell) !y0
               globalReprodSum9fld2(iCell,4) = SIvec_3(iCell) !y0
               globalReprodSum9fld2(iCell,5) = SIvec_2(iCell) !q0
               globalReprodSum9fld2(iCell,6) = SIvec_3(iCell) !y0
               globalReprodSum9fld2(iCell,7) = SIvec_6(iCell) !t0
               globalReprodSum9fld2(iCell,8) = SIvec_7(iCell) !v0
               globalReprodSum9fld2(iCell,9) = SIvec_2(iCell) !q0
            end do
#ifndef MPAS_OPENACC
            !$omp end do
            !$omp end parallel
#endif
    
            !$acc update host(reprodSum1,reprodSum2)
            globalProd(:) =                                &
               mpas_global_sum_nfld(reprodSum1,reprodSum2, &
                                    domain%dminfo%comm)
    
         else
 
            !$acc parallel loop async present(SIcst_q0q0,                      &
            !$acc       SIcst_r00s0,SIcst_r00z0,SIcst_q0y0 ,SIcst_y0y0,  &
            !$acc       SIcst_r00q0,SIcst_r00y0,SIcst_r00t0,SIcst_r00v0, &
            !$acc       SIvec_1,SIvec_2,SIvec_3,SIvec_4,SIvec_5,SIvec_6, &
            !$acc       SIvec_7)                                         &
            !$acc    reduction(+:SIcst_q0q0,                             &
            !$acc       SIcst_r00s0,SIcst_r00z0,SIcst_q0y0 ,SIcst_y0y0,  &
            !$acc       SIcst_r00q0,SIcst_r00y0,SIcst_r00t0,SIcst_r00v0)
            do iCell = 1,nCellsOwned
               SIcst_r00s0 = SIcst_r00s0 + SIvec_1(iCell) &
                                         * SIvec_4(iCell) ! s1
       
               SIcst_r00z0 = SIcst_r00z0 + SIvec_1(iCell) &
                                         * SIvec_5(iCell) ! z1
       
               SIcst_q0y0  = SIcst_q0y0  + SIvec_2(iCell)  &
                                         * SIvec_3(iCell)
       
               SIcst_y0y0  = SIcst_y0y0  + SIvec_3(iCell)  &
                                         * SIvec_3(iCell)
       
               SIcst_r00q0 = SIcst_r00q0 + SIvec_1(iCell) &
                                         * SIvec_2(iCell)
       
               SIcst_r00y0 = SIcst_r00y0 + SIvec_1(iCell) &
                                         * SIvec_3(iCell)
       
               SIcst_r00t0 = SIcst_r00t0 + SIvec_1(iCell) &
                                         * SIvec_6(iCell)
       
               SIcst_r00v0 = SIcst_r00v0 + SIvec_1(iCell) &
                                         * SIvec_7(iCell)
       
               SIcst_q0q0 = SIcst_q0q0   + SIvec_2(iCell) &
                                         * SIvec_2(iCell)
            end do
    
            !$acc serial async present(localProd,SIcst_q0q0,               &
            !$acc    SIcst_r00s0,SIcst_r00z0,SIcst_q0y0 ,SIcst_y0y0, &
            !$acc    SIcst_r00q0,SIcst_r00y0,SIcst_r00t0,SIcst_r00v0)
            localProd(1) = SIcst_r00s0
            localProd(2) = SIcst_r00z0
            localProd(3) = SIcst_q0y0
            localProd(4) = SIcst_y0y0
            localProd(5) = SIcst_r00q0
            localProd(6) = SIcst_r00y0
            localProd(7) = SIcst_r00t0
            localProd(8) = SIcst_r00v0
            localProd(9) = SIcst_q0q0
            SIcst_r00s0 = 0.0_RKIND
            SIcst_r00z0 = 0.0_RKIND
            SIcst_q0y0  = 0.0_RKIND
            SIcst_y0y0  = 0.0_RKIND
            SIcst_r00q0 = 0.0_RKIND
            SIcst_r00y0 = 0.0_RKIND
            SIcst_r00t0 = 0.0_RKIND
            SIcst_r00v0 = 0.0_RKIND
            SIcst_q0q0  = 0.0_RKIND
            !$acc end serial
            !$acc wait

            ! Global sum across CPUs
            if ( threadNum == 0 ) then
#ifdef USE_GPU_AWARE_MPI
            !$acc host_data use_device(localProd,globalProd)
#else
            !$acc update host(localProd)
#endif
            call mpas_dmpar_sum_real_array(domain % dminfo, nfld,      &
                                           localProd,globalProd)
#ifdef USE_GPU_AWARE_MPI
            !$acc end host_data
#else
            !$acc update device(globalProd) async
#endif
            endif ! thread

         endif

      endif !!! si_algorithm & nfld

      call mpas_timer_stop("si global reduction")
   
   end subroutine si_global_reduction !}}}

!***********************************************************************
!
!  routine si_solver_scg
!
!> \brief   Linear iterative solver:
!              The s-step conjugate gradient algorithm
!> \author  Hyun-Gyu Kang
!
!> \details
!   This routine solves a linear system iteratively using 
!   the s-step conjugate gradient algorithm 
!   (Chronopoulos and Gear, 1989).
!
!-----------------------------------------------------------------------

   subroutine si_solver_scg(domain,sshSolved,ssh,bottomDepthEdge, & !{{{
                            tolerance)

      type (domain_type), intent(inout) :: &
         domain  !< [inout] model state to advance forward

      real(kind=RKIND), dimension(:), pointer, intent(in) :: &
         ssh                                    ! Intermediate ssh

      real(kind=RKIND), dimension(:), pointer, intent(inout) :: &
         sshSolved                              ! ssh to be solved

      real(kind=RKIND), dimension(:), intent(inout) :: &
         bottomDepthEdge                        ! Bottome depth at Edge

      real(kind=RKIND), intent(in) :: tolerance ! Tolerance of solver

      integer :: iter, iCell

      iter = 0

      resid = (tolerance+100.0)**2.0

      !-------------------------------------------------------------!
      do while ( dsqrt(resid) > tolerance )
      !-------------------------------------------------------------!

         iter = iter + 1

#ifdef MPAS_OPENACC
         !$acc parallel loop gang                                   &
         !$acc    present(SIvec_z1,SIvec_rh0,SIcst_beta0,           &
         !$acc       SIvec_z0,SIvec_s1,SIvec_w0,SIvec_s0,sshSolved, &
         !$acc       SIcst_alpha0,SIvec_r1,SIvec_r0)
#else
         !$omp parallel
         !$omp do schedule(runtime)
#endif
         do iCell = 1, nPrecVec
            SIvec_z1(iCell) = SIvec_rh0(iCell) &
                            + SIcst_beta0 * SIvec_z0(iCell)

            SIvec_s1(iCell) = SIvec_w0(iCell) &
                            + SIcst_beta0 * SIvec_s0(iCell)

            sshSolved(iCell) = sshSolved(iCell) &
                            + SIcst_alpha0 * SIvec_z1(iCell)

            SIvec_r1(iCell) = SIvec_r0(iCell) &
                            - SIcst_alpha0 * SIvec_s1(iCell)
         end do ! iCell
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
  
         ! Preconditioning -----------------------------------------!
         call si_precond(SIvec_r1,SIvec_rh1)
  
         call mpas_timer_start("si halo iter")
         call mpas_dmpar_exch_halo_field(SIvec_rh1_field)
#ifndef USE_GPU_AWARE_MPI
         !$acc update device(SIvec_rh1)
#endif
         call mpas_timer_stop("si halo iter")
  
  
         ! SpMV ----------------------------------------------------!
         call si_matvec_mul(SIvec_rh1,bottomDepthEdge,SIvec_w1,ssh)
  
         ! Reduction -----------------------------------------------!

         call si_global_reduction(SIcst_allreduce_local3,&
                             SIcst_allreduce_global3,    &
                             globalReprodSum3fld1,       &
                             globalReprodSum3fld2,       &
                             si_algorithm,3,domain,      &
                             SIvec_r1,SIvec_w1,SIvec_rh1)

         !$acc serial present(SIcst_r0rh0,SIcst_w0rh0,SIcst_gamma1,  &
         !$acc    SIcst_allreduce_global3,SIcst_omega0,SIcst_beta1,  &
         !$acc    SIcst_alpha1,SIcst_gamma0,SIcst_alpha0,resid,      &
         !$acc    SIcst_beta0,SIcst_r0r0)
         SIcst_r0rh0 = SIcst_allreduce_global3(1)
         SIcst_w0rh0 = SIcst_allreduce_global3(2)
               resid = SIcst_allreduce_global3(3)
         SIcst_gamma1 = SIcst_r0rh0
         SIcst_omega0 = SIcst_w0rh0
         SIcst_beta1  = SIcst_gamma1 / SIcst_gamma0
         SIcst_alpha1 = SIcst_gamma1                            &
                      /(SIcst_omega0 - SIcst_beta1*SIcst_gamma1 &
                                     / SIcst_alpha0 )
         SIcst_beta0  = SIcst_beta1
         SIcst_alpha0 = SIcst_alpha1
         SIcst_gamma0 = SIcst_gamma1

         SIcst_r0rh0 = 0.0_RKIND
         SIcst_w0rh0 = 0.0_RKIND
         SIcst_r0r0  = 0.0_RKIND
         !$acc end serial
         !$acc update host(resid)

#ifdef MPAS_OPENACC
         !$acc parallel loop gang                              &
         !$acc    present(SIvec_r0,SIvec_r1,SIvec_s0,SIvec_s1, &
         !$acc       SIvec_w0,SIvec_w1,SIvec_z0,SIvec_z1,      &
         !$acc       SIvec_rh0,SIvec_rh1)
#else
         !$omp parallel
         !$omp do schedule(runtime)
#endif
         do iCell = 1,nPrecVec
            SIvec_r0(iCell)  = SIvec_r1(iCell)
            SIvec_s0(iCell)  = SIvec_s1(iCell)
            SIvec_w0(iCell)  = SIvec_w1(iCell)
            SIvec_z0(iCell)  = SIvec_z1(iCell)
            SIvec_rh0(iCell) = SIvec_rh1(iCell)
         end do ! iCell
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif

         if ( iter > int(mean_num_cells*5) ) then
            call mpas_log_write(&
            '******************************************************')
            call mpas_log_write(&
            'Iteration number exceeds Max. #iteration: PROGRAM STOP')
            call mpas_log_write(&
            'Current #Iteration = $i ', intArgs=(/ iter /) )
            call mpas_log_write(&
            'Max.    #Iteration = $i ', &
            intArgs=(/ int(mean_num_cells*5) /) )
            call mpas_log_write(&
            '******************************************************')
            call mpas_log_write('',MPAS_LOG_CRIT)
         endif

      !-------------------------------------------------------------!
      end do ! do iter
      !-------------------------------------------------------------!
     
   end subroutine si_solver_scg !}}}

!***********************************************************************
!
!  routine si_solver_sbicg
!
!> \brief   Linear iterative solver:
!              The single-reduction bi-conjugate gradient
!              stabilization algorithm
!
!> \author  Hyun-Gyu Kang
!
!> \details
!   This routine solves a linear system iteratively using 
!   the single-reduction bi-conjugate gradient stabilization algorithm.
!
!-----------------------------------------------------------------------

   subroutine si_solver_sbicg(domain,sshSolved,ssh,bottomDepthEdge, & !{{{
                 tolerance)

      type (domain_type), intent(inout) :: &
         domain  !< [inout] model state to advance forward

      real(kind=RKIND), dimension(:), pointer, intent(in) :: &
         ssh                                    ! Intermediate ssh

      real(kind=RKIND), dimension(:), pointer, intent(inout) :: &
         sshSolved                              ! ssh to be solved

      real(kind=RKIND), dimension(:), intent(inout) :: &
         bottomDepthEdge                        ! Bottome depth at Edge

      real(kind=RKIND), intent(in) :: tolerance ! Tolerance of solver

      integer :: iter, iCell

      iter = 0
      resid = (tolerance+100.0)**2.0

      !-------------------------------------------------------------!
      do while ( dsqrt(resid) > tolerance )
      !-------------------------------------------------------------!
 
         iter = iter + 1
 
#ifdef MPAS_OPENACC
         !$acc parallel loop gang async present(                    &
         !$acc       SIvec_ph1,SIvec_rh0,SIvec_ph0,                 &
         !$acc       SIvec_sh0,SIvec_s1,SIvec_w0,                   &
         !$acc       SIvec_s0,SIvec_z0,SIvec_sh1,SIvec_wh0,         &
         !$acc       SIvec_zh0,SIvec_z1,SIvec_t0,SIvec_v0,SIvec_q0, &
         !$acc       SIvec_r0,SIvec_qh0,SIvec_y0,                   &
         !$acc       SIcst_alpha1,SIcst_beta1,SIcst_omega0)
#else
         !$omp parallel
         !$omp do schedule(runtime)
#endif
         do iCell = 1, nPrecVec
            SIvec_ph1(iCell) =  SIvec_rh0(iCell) + SIcst_beta1      &
                             * (SIvec_ph0(iCell) - SIcst_omega0     &
                                                 * SIvec_sh0(iCell))

            SIvec_s1(iCell)  =  SIvec_w0(iCell)  + SIcst_beta1      &
                             * (SIvec_s0(iCell)  - SIcst_omega0     &
                                                 * SIvec_z0(iCell))

            SIvec_sh1(iCell) =  SIvec_wh0(iCell) + SIcst_beta1      &
                             * (SIvec_sh0(iCell) - SIcst_omega0     &
                                                 * SIvec_zh0(iCell))

            SIvec_z1(iCell)  =  SIvec_t0(iCell)  + SIcst_beta1      &
                             * (SIvec_z0(iCell)  - SIcst_omega0     &
                                                 * SIvec_v0(iCell))

            SIvec_q0(iCell)  =  SIvec_r0(iCell)  - SIcst_alpha1     &
                                                 * SIvec_s1(iCell)

            SIvec_qh0(iCell) =  SIvec_rh0(iCell) - SIcst_alpha1     &
                                                 * SIvec_sh1(iCell)

            SIvec_y0(iCell)  =  SIvec_w0(iCell)  - SIcst_alpha1     &
                                                 * SIvec_z1(iCell)
         end do ! iCell
#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
 
         ! Preconditioning -----------------------------------------!
         call si_precond(SIvec_z1,SIvec_zh1)
         !$acc wait
 
         call mpas_timer_start("si halo iter")
#ifdef USE_GPU_AWARE_MPI
         call si_halo_exch(domain,SIvec_zh1)
#else
         call mpas_dmpar_exch_halo_field(SIvec_zh1_field)
         !$acc update device(SIvec_zh1) async
#endif
         call mpas_timer_stop("si halo iter")
 
 
         ! SpMV ----------------------------------------------------!
         call si_matvec_mul(SIvec_zh1,bottomDepthEdge,SIvec_v0,ssh)

         ! Reduction -----------------------------------------------!
         call si_global_reduction(SIcst_allreduce_local9,          &
                             SIcst_allreduce_global9,              &
                             globalReprodSum9fld1,                 &
                             globalReprodSum9fld2,                 &
                             si_algorithm,9,domain,                &
                             SIvec_r00,SIvec_q0,SIvec_y0,SIvec_s1, &
                             SIvec_z1,SIvec_t0,SIvec_v0)
 
         !$acc serial async present(SIcst_rho1,SIcst_r00q0,          &
         !$acc    SIcst_gamma1,SIcst_r00y0,SIcst_r00t0,SIcst_omega0, &
         !$acc    SIcst_alpha0,SIcst_r00v0,SIcst_beta0,SIcst_rho0,   &
         !$acc    SIcst_r00s0,SIcst_r00z0,SIcst_q0q0,resid,          &
         !$acc    SIcst_q0y0,SIcst_y0y0,SIcst_alpha1,SIcst_beta1,    &
         !$acc    SIcst_allreduce_global9)
         SIcst_r00s0 = SIcst_allreduce_global9(1)
         SIcst_r00z0 = SIcst_allreduce_global9(2)
         SIcst_q0y0  = SIcst_allreduce_global9(3)
         SIcst_y0y0  = SIcst_allreduce_global9(4)
         SIcst_r00q0 = SIcst_allreduce_global9(5)
         SIcst_r00y0 = SIcst_allreduce_global9(6)
         SIcst_r00t0 = SIcst_allreduce_global9(7)
         SIcst_r00v0 = SIcst_allreduce_global9(8)
         SIcst_q0q0  = SIcst_allreduce_global9(9)

         SIcst_alpha0 = SIcst_alpha1
         SIcst_beta0  = SIcst_beta1

         SIcst_omega0 = SIcst_q0y0 / SIcst_y0y0

         resid = SIcst_q0q0 &
               - 2.0_RKIND * SIcst_omega0 * SIcst_q0y0  &
               + SIcst_omega0**2.0_RKIND  * SIcst_y0y0

         SIcst_rho1 = SIcst_r00q0 - SIcst_omega0 * SIcst_r00y0

         SIcst_gamma1 = SIcst_r00y0 - SIcst_omega0 * SIcst_r00t0  &
                                    + SIcst_omega0 * SIcst_alpha0 &
                                                   * SIcst_r00v0

         SIcst_beta1 = (SIcst_alpha0/SIcst_omega0) &
                     * (SIcst_rho1 / SIcst_rho0)

         SIcst_alpha1 = SIcst_rho1 &
                      / ( SIcst_gamma1 + SIcst_beta1 * SIcst_r00s0   &
                                       - SIcst_beta1 * SIcst_omega0  &
                                                     * SIcst_r00z0 )
         SIcst_rho0 = SIcst_rho1

         SIcst_r00s0 = 0.0_RKIND
         SIcst_r00z0 = 0.0_RKIND
         SIcst_q0y0  = 0.0_RKIND
         SIcst_y0y0  = 0.0_RKIND
         SIcst_r00q0 = 0.0_RKIND
         SIcst_r00y0 = 0.0_RKIND
         SIcst_r00t0 = 0.0_RKIND
         SIcst_r00v0 = 0.0_RKIND
         SIcst_q0q0  = 0.0_RKIND
         !$acc end serial
         !$acc update host(resid) async

#ifdef MPAS_OPENACC
         !$acc parallel loop gang present(sshSolved,          &
         !$acc       SIvec_qh0,SIvec_q0,SIvec_y0,SIvec_wh0,   &
         !$acc       SIvec_t0,SIvec_v0,SIvec_w0, SIvec_r0,    &
         !$acc       SIvec_s0,SIvec_s1,SIvec_z0,SIvec_z1,     &
         !$acc       SIvec_rh0,SIvec_sh0,SIvec_sh1,           &
         !$acc       SIvec_ph0,SIvec_ph1,SIvec_zh0,SIvec_zh1, &
         !$acc       SIcst_alpha0,SIcst_omega0)
#else
         !$omp parallel
         !$omp do schedule(runtime)
#endif
         do iCell = 1,nPrecVec
            sshSolved(iCell) = sshSolved(iCell)                &
                             + SIcst_alpha0 * SIvec_ph1(iCell) &
                             + SIcst_omega0 * SIvec_qh0(iCell)

            SIvec_r0(iCell)  = SIvec_q0(iCell) - SIcst_omega0       &
                             * SIvec_y0(iCell)

            SIvec_rh0(iCell) =  SIvec_qh0(iCell) - SIcst_omega0     &
                             * (SIvec_wh0(iCell) - SIcst_alpha0     &
                                                 * SIvec_zh1(iCell))

            SIvec_w0(iCell)  =   SIvec_y0(iCell) - SIcst_omega0     &
                             * ( SIvec_t0(iCell) - SIcst_alpha0     &
                                                 * SIvec_v0(iCell))

            SIvec_s0(iCell)  = SIvec_s1(iCell)
            SIvec_z0(iCell)  = SIvec_z1(iCell)
            SIvec_sh0(iCell) = SIvec_sh1(iCell)
            SIvec_ph0(iCell) = SIvec_ph1(iCell)
            SIvec_zh0(iCell) = SIvec_zh1(iCell)
         end do

#ifndef MPAS_OPENACC
         !$omp end do
         !$omp end parallel
#endif
         !$acc wait
 
         ! Preconditioning -----------------------------------------!
         call si_precond(SIvec_w0,SIvec_wh0) !wh1 -> wh0; w1->w0
 
         call mpas_timer_start("si halo iter")
#ifdef USE_GPU_AWARE_MPI
         call si_halo_exch(domain,SIvec_wh0)
#else
         call mpas_dmpar_exch_halo_field(SIvec_wh0_field)
         !$acc update device(SIvec_wh0) async
#endif
         call mpas_timer_stop("si halo iter")
 
         ! SpMV ----------------------------------------------------!
         call si_matvec_mul(SIvec_wh0,bottomDepthEdge,SIvec_t0,ssh) 
                                             !wh1 -> wh0; t1 -> t0
         !$acc wait
 
         if ( iter > int(mean_num_cells*5) ) then
            call mpas_log_write(&
            '******************************************************')
            call mpas_log_write(&
            'Iteration number exceeds Max. #iteration: PROGRAM STOP')
            call mpas_log_write(&
            'Current #Iteration = $i ', intArgs=(/ iter /) )
            call mpas_log_write(&
            'Max.    #Iteration = $i ', &
            intArgs=(/ int(mean_num_cells*5) /) )
            call mpas_log_write(&
            '******************************************************')
            call mpas_log_write('',MPAS_LOG_CRIT)
         endif

      !-------------------------------------------------------------!
      end do ! do iter
      !-------------------------------------------------------------!

   end subroutine si_solver_sbicg !}}}

!***********************************************************************
!
!  routine si_halo_exch
!
!> \brief  Halo exchange for the SI baroropic mode solver
!
!> \author  Hyun-Gyu Kang
!
!> \details
!   This routine performs local halo exchanges for the semi-implicit
!   barotropic mode solver.
!   This routine enables the GPU-aware MPI halo exchanges during
!   solver iterations.
!
!-----------------------------------------------------------------------
#ifdef USE_GPU_AWARE_MPI

   subroutine si_halo_exch(domain,field_array) !{{{

      include 'mpif.h'

      type (domain_type), intent(inout) :: &
         domain  !< [inout] model state to advance forward
   
      real(kind=RKIND),dimension(:),intent(inout) :: field_array

      type (mpas_exchange_list), pointer :: exchListPtr

      integer,dimension(mpi_status_size) :: stats

      integer :: i,iComm, iHalo, iExch, mpi_ierr, threadNum

      threadNum = mpas_threading_get_thread_num()
     
      if ( threadNum == 0 ) then

      do iComm = 1,nRecvCommList
         call mpas_timer_start("irecv")
         ! Initiate mpi_irecv calls ---------------------------------------
         rbuffer => exchRecv(iComm)%buffer
         !$acc host_data use_device(rbuffer)
         call MPI_Irecv(                                               &
                             rbuffer, exchRecv(iComm)%nList, MPI_REAL8,&
              exchRecv(iComm)%procID, exchRecv(iComm)%procID,          &
              domain%dminfo%comm, exchRecv(iComm)%reqID, mpi_ierr )
         !$acc end host_data
         call mpas_timer_stop("irecv")

         ! Copy data into buffer, and initiate mpi_isend calls ------------
         call mpas_timer_start("send loop")
         !$acc parallel loop present(exchSend,field_array)
         do iHalo = 1,config_num_halos
            do iExch = 1,exchSend(iComm)%nExch(iHalo)
               if (exchSend(iComm)%halo(iHalo)%endPointID(iExch) == &
                   exchSend(iComm)%procID) then
                   !$acc loop vector
                   do i = 1,exchSend(iComm)%halo(iHalo)%nList(iExch)
                      exchSend(iComm)%buffer( exchSend(iComm)%  &
                                                  halo(iHalo)%  &
                                                  exch(iExch)%  &
                                                 destList(i)) = &
                      field_array(exchSend(iComm)%halo(iHalo)   &
                                     %exch(iExch)%srcList(i))
                   end do ! i
               endif
            end do ! iExch
         end do ! iHalo
         call mpas_timer_stop("send loop")

         call mpas_timer_start("isend")
         sbuffer => exchSend(iComm)%buffer
         !$acc host_data use_device(sbuffer)
         call MPI_Isend(                                               &
                             sbuffer, exchSend(iComm)%nList, MPI_REAL8,&
              exchSend(iComm)%procID, domain%dminfo%my_proc_id,        &
              domain%dminfo%comm, exchSend(iComm)%reqID, mpi_ierr )
         !$acc end host_data
         call mpas_timer_stop("isend")

      end do ! iComm

      ! Wait for mpi_irecv to finish, and unpack data from buffer ------
      do iComm = 1,nRecvCommList
         call MPI_Wait(exchRecv(iComm)%reqID,stats,mpi_ierr)

         call mpas_timer_start("recv loop")
         !$acc parallel loop present(exchRecv,field_array)
         do iHalo = 1,config_num_halos
            do iExch = 1,exchRecv(iComm)%nExch(iHalo)
               if (exchRecv(iComm)%halo(iHalo)%endPointID(iExch) == &
                   exchRecv(iComm)%procID) then
                   !acc loop vector
                   do i = 1,exchRecv(iComm)%halo(iHalo)%nList(iExch)
                      field_array(exchRecv(iComm)%halo(iHalo)    &
                                     %exch(iExch)%destList(i)) = &
                      exchRecv(iComm)%buffer(exchRecv(iComm)%    &
                                                 halo(iHalo)%    &
                                                 exch(iExch)%    &
                                                 srcList(i))
                   end do ! i
               endif
            end do ! iExch
         end do ! iHalo
         call mpas_timer_stop("recv loop")
      end do ! iComm

      ! Wait for mpi_isend to finish ----------------------------------
      do iComm = 1,nSendCommList
         call MPI_Wait(exchSend(iComm)%reqID,MPI_STATUS_IGNORE,mpi_ierr)
      end do

      endif ! threadNum

   end subroutine si_halo_exch

#endif

!***********************************************************************
!
!  routine ocn_time_integrator_si_variable_destroy
!
!> \brief   Destroy SI solver-related variables
!
!> \author  Hyun-Gyu Kang
!
!> \details
!   This routine deallocates CPU arrays allocated in the SI solver.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_si_variable_destroy() !{{{

      integer :: info,iComm,iHalo,iExch

#ifdef MPAS_OPENACC

#ifdef USE_MAGMA
      info = magmaf_free( dprec_ivmat )
      info = magmaf_free( dxvec )
      info = magmaf_free( dyvec )
#endif

      !$acc exit data delete(prec_ivmat,                             &
      !$acc    SIcst_allreduce_local2,SIcst_allreduce_global2,       &
      !$acc    SIcst_allreduce_local3,SIcst_allreduce_global3,       &
      !$acc    SIcst_allreduce_local9,SIcst_allreduce_global9,       &
      !$acc    SIcst_alpha0,SIcst_alpha1,SIcst_beta0,SIcst_beta1,    &
      !$acc    SIcst_omega0,SIcst_gamma0,SIcst_gamma1,resid,         &
      !$acc    SIcst_r0rh0,SIcst_w0rh0,SIcst_r0r0 ,SIcst_r00r0,      &
      !$acc    SIcst_r00s0,SIcst_r00z0,SIcst_q0y0 ,SIcst_y0y0 ,      &
      !$acc    SIcst_r00y0,SIcst_r00t0,SIcst_r00v0,SIcst_q0q0 ,      &
      !$acc    SIcst_r00w0,SIcst_r00q0,SIcst_rho0,SIcst_rho1)

      if ( config_btr_si_partition_match_mode ) then
      !$acc exit data delete(                                        &
      !$acc               globalReprodSum2fld1,globalReprodSum2fld2, &
      !$acc               globalReprodSum3fld1,globalReprodSum3fld2 )
      endif

#ifdef USE_GPU_AWARE_MPI
      do iComm = 1,nRecvCommList
         do iHalo = 1,config_num_halos
            do iExch = 1,exchRecv(iComm)%nExch(iHalo)
               !$acc exit data delete(                                    &
               !$acc    exchRecv(iComm)%halo(iHalo)%exch(iExch)%destList, &
               !$acc    exchRecv(iComm)%halo(iHalo)%exch(iExch)%srcList )
            end do
            !$acc exit data delete(exchRecv(iComm)%halo(iHalo)%nList,       &
            !$acc                   exchRecv(iComm)%halo(iHalo)%endPointID, &
            !$acc                   exchRecv(iComm)%halo(iHalo)%exch )
         end do
         !$acc exit data delete(exchRecv(iComm)%buffer,  &
         !$acc                   exchRecv(iComm)%nExch,  &
         !$acc                   exchRecv(iComm)%halo,   &
         !$acc                   exchRecv(iComm)%procID)
      end do

      do iComm = 1,nSendCommList
         do iHalo = 1,config_num_halos
            do iExch = 1,exchSend(iComm)%nExch(iHalo)
               !$acc exit data delete(                                    &
               !$acc    exchSend(iComm)%halo(iHalo)%exch(iExch)%destList, &
               !$acc    exchSend(iComm)%halo(iHalo)%exch(iExch)%srcList)
            end do
            !$acc exit data delete(exchSend(iComm)%halo(iHalo)%nList,       &
            !$acc                   exchSend(iComm)%halo(iHalo)%endPointID, &
            !$acc                   exchSend(iComm)%halo(iHalo)%exch)
         end do
         !$acc exit data delete(exchSend(iComm)%buffer,  &
         !$acc                   exchSend(iComm)%nExch,  &
         !$acc                   exchSend(iComm)%halo,   &
         !$acc                   exchSend(iComm)%procID)
      end do

      !$acc exit data delete(exchSend,exchRecv,sbuffer,rbuffer)
#endif

#endif

      deallocate(prec_ivmat)

      if ( config_btr_si_partition_match_mode ) then
         deallocate( globalReprodSum2fld1,  &
                     globalReprodSum2fld2,  &
                     globalReprodSum3fld1,  &
                     globalReprodSum3fld2 )
      endif

   end subroutine ocn_time_integrator_si_variable_destroy !}}}

!***********************************************************************

end module ocn_time_integration_si

! vim: foldmethod=marker
