!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_velocity_external
!
!> \MPAS land-ice velocity driver for external dycores
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for interfacing with 
!>  external velocity solvers.  These currently are LifeV (L1L2, First order),
!>  Albany (First order), and PHG (Stokes).
!>
!
!-----------------------------------------------------------------------

module li_velocity_external

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timer
   use li_setup
   use, intrinsic :: iso_c_binding
   use mpas_constants, only: gravity
   use li_constants

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

   public :: li_velocity_external_init, &
             li_velocity_external_block_init, &
             li_velocity_external_solve, &
             li_velocity_external_finalize

   interface
      ! Note: Could add all interface routines to this interface...
      ! For now, just trying it with this new routine.
      subroutine velocity_solver_set_parameters(gravity, config_ice_density, config_ocean_density, config_sea_level, config_default_flowParamA, config_enhancementFactor, &
                         config_flowLawExponent, config_dynamic_thickness, li_mask_ValueDynamicIce, li_mask_ValueIce) bind(C, name="velocity_solver_set_parameters")
         use iso_c_binding, only: C_INT, C_DOUBLE
         
         INTEGER(C_INT) :: li_mask_ValueDynamicIce, li_mask_ValueIce
         REAL(C_DOUBLE) :: config_ice_density, config_ocean_density, config_sea_level, config_default_flowParamA, &
                           config_enhancementFactor, config_flowLawExponent, config_dynamic_thickness, gravity
      end subroutine velocity_solver_set_parameters

   end interface

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************



contains



!***********************************************************************
!
!  routine li_velocity_external_init
!
!> \brief   Initializes velocity solver
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine initializes the ice velocity solver in
!>  external velocity solvers.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_external_init(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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
      integer, pointer :: config_num_halos, config_number_of_blocks
      character (len=StrKIND), pointer :: config_velocity_solver
      integer :: err_tmp

      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_num_halos', config_num_halos)
      call mpas_pool_get_config(liConfigs, 'config_number_of_blocks', config_number_of_blocks)
      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)

      ! Check for configuration options that are incompatible with external velocity solver conventions
      if (config_num_halos < 2) then
         write(stderrUnit,*) "Error: External velocity solvers require that config_num_halos >= 2"
         err_tmp = 1
      endif
      err = ior(err,err_tmp)

      if (config_number_of_blocks /= 0) then
         write(stderrUnit,*) "Error: External velocity solvers require that config_number_of_blocks=0"
         err_tmp = 1
      endif
      err = ior(err,err_tmp)

      ! Check if we are on a sphere - not supported by external dycores
      if (domain % on_a_sphere) then
         write(stderrUnit,*) "ERROR: External velocity solvers cannot be run with a spherical mesh."
         err_tmp = 1
      endif
      err = ior(err,err_tmp)

      ! These calls are needed for setting up the external velocity solvers
#if defined(USE_EXTERNAL_L1L2) || defined(USE_EXTERNAL_FIRSTORDER) || defined(USE_EXTERNAL_STOKES)
      !call external first order solver to set the grid of the velocity solver
      write(stdoutUnit,*) "Initializing external velocity solver."
      flush(stdoutUnit)  ! Flush log files before handing control to C++
      flush(stderrUnit)

      call velocity_solver_init_mpi(domain % dminfo % comm)
#else
      err = 1
      write(stderrUnit,*) "Error: To run with an external velocity solver you must compile MPAS with one."
#endif


      if (config_velocity_solver == 'Stokes') then
#ifdef USE_EXTERNAL_STOKES
          call interface_phg_init(domain, err)
#else
          write(stderrUnit,*) "Error: External Stokes library needed to run Stokes dycore."
          err = 1
          return
#endif
      endif
      err = ior(err,err_tmp)


      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_velocity_external_init."
      endif

   !--------------------------------------------------------------------
   end subroutine li_velocity_external_init



!***********************************************************************
!
!  routine li_velocity_external_block_init
!
!> \brief   Initializes blocks for external velocity solver use
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine initializes each block of the ice velocity solver in the
!>  external velocity solver.
!>   Note: LifeV/Albany/PHG only support one block per processor, but this has (hopefully)
!>   been written to work if that were to change.  (That's why all these external dycore init
!>   calls are in li_velocity_external_block_init instead of li_velocity_external_init.)
!
!-----------------------------------------------------------------------

   subroutine li_velocity_external_block_init(block, err)

      use li_mask

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (block_type), intent(in) :: &
         block          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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
      type (mpas_pool_type), pointer :: meshPool
      integer, pointer :: nCells, nEdges, nVertices, nCellsSolve, nEdgesSolve, nVerticesSolve, nVertInterfaces, maxNEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnCell, verticesOnEdge, edgesOnCell
      integer, dimension(:), pointer :: indexToCellID, indexToEdgeID, indexToVertexID, nEdgesOnCell
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, areaTriangle
      real (kind=RKIND), pointer :: radius
      type (field1DInteger), pointer :: indexToCellIDField, indexToEdgeIDField, indexToVertexIDField
      real (kind=RKIND), pointer :: config_ice_density, config_ocean_density,  config_sea_level, config_default_flowParamA, &
                                    config_enhancementFactor, config_flowLawExponent, config_dynamic_thickness

      ! halo exchange arrays
      integer, dimension(:), pointer :: sendCellsArray, &
                                        recvCellsArray, &
                                        sendVerticesArray, &
                                        recvVerticesArray, &
                                        sendEdgesArray, &
                                        recvEdgesArray

      err = 0

      !extract data from domain
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
      call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertInterfaces', nVertInterfaces)
      call mpas_pool_get_dimension(meshPool, 'maxEdges', maxNEdgesOnCell)
      call mpas_pool_get_config(meshPool, 'sphere_radius', radius)

      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(meshPool, 'indexToEdgeID', indexToEdgeID)
      call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'areaTriangle', areaTriangle)

      call mpas_pool_get_field(meshPool, 'indexToCellID', indexToCellIDField)
      call mpas_pool_get_field(meshPool, 'indexToEdgeID', indexToEdgeIDField)
      call mpas_pool_get_field(meshPool, 'indexToVertexID', indexToVertexIDField)

      ! build send and receive arrays using exchange_list
      call array_from_exchange_list(indexToCellIDField, sendCellsArray, recvCellsArray)
      call array_from_exchange_list(indexToEdgeIDField, sendEdgesArray, recvEdgesArray)
      call array_from_exchange_list(indexToVertexIDField, sendVerticesArray, recvVerticesArray)


#if defined(USE_EXTERNAL_L1L2) || defined(USE_EXTERNAL_FIRSTORDER) || defined(USE_EXTERNAL_STOKES)
      ! These calls are needed for using any of the external velocity solvers

      !zCell is supposed to be zero when working on planar geometries (radius = 0)
      !nVertLevels should be equal to nVertLevelsSolve (no splitting of the domain in the vertical direction)
      call mpas_timer_start("velocity_solver_set_grid_data")
      write(stdoutUnit,*) "Initializing external velocity solver grid data."
      flush(stdoutUnit)  ! Flush log files before handing control to C++
      flush(stderrUnit)

      call velocity_solver_set_grid_data(nCells, nEdges, nVertices, nVertInterfaces, &
              nCellsSolve, nEdgesSolve, nVerticesSolve, maxNEdgesOnCell, radius, &
              cellsOnEdge, cellsOnVertex, verticesOnCell, verticesOnEdge, edgesOnCell, &
              nEdgesOnCell, indexToCellID, &
              xCell, yCell, zCell, xVertex, yVertex, zVertex, areaTriangle, &
              sendCellsArray, recvCellsArray, &
              sendEdgesArray, recvEdgesArray, &
              sendVerticesArray, recvVerticesArray)
      call mpas_timer_stop("velocity_solver_set_grid_data")
#else
      write(stderrUnit,*) "Error: To run with an external velocity solver you must compile MPAS with one."
      err = 1
#endif

      !these can be deallocated because they have been copied on the c++ side
      deallocate(sendCellsArray, &
                 recvCellsArray, &
                 sendVerticesArray, &
                 recvVerticesArray, &
                 sendEdgesArray, &
                 recvEdgesArray)

      ! Set physical parameters needed on the other side
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)
      call mpas_pool_get_config(liConfigs, 'config_ocean_density', config_ocean_density)
      call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
      call mpas_pool_get_config(liConfigs, 'config_default_flowParamA', config_default_flowParamA)
      call mpas_pool_get_config(liConfigs, 'config_enhancementFactor', config_enhancementFactor)
      call mpas_pool_get_config(liConfigs, 'config_flowLawExponent', config_flowLawExponent)
      call mpas_pool_get_config(liConfigs, 'config_dynamic_thickness', config_dynamic_thickness)
#if defined(USE_EXTERNAL_L1L2) || defined(USE_EXTERNAL_FIRSTORDER) || defined(USE_EXTERNAL_STOKES)
      call velocity_solver_set_parameters(gravity, config_ice_density, config_ocean_density, config_sea_level, config_default_flowParamA, config_enhancementFactor, &
                                          config_flowLawExponent, config_dynamic_thickness, li_mask_ValueAlbanyActive, li_mask_ValueIce)
#endif


      ! === error check
      if (err > 0) then
          write (stderrUnit,*) "An error has occurred in li_velocity_external_block_init."
      endif

   !--------------------------------------------------------------------
   end subroutine li_velocity_external_block_init



!***********************************************************************
!
!  routine li_velocity_external_solve
!
!> \brief   Interface to call external velocity solvers
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine calls external first-order velocity solvers and/or Stokes velocity solvers.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_external_solve(meshPool, geometryPool, thermalPool, scratchPool, velocityPool, err)

       use li_mask

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
          meshPool !< Input: mesh information

      type (mpas_pool_type), intent(in) :: &
          geometryPool !< Input: geometry information

      type (mpas_pool_type), intent(in) :: &
          thermalPool !< Input: thermal information

      type (mpas_pool_type), intent(in) :: &
          scratchPool !< Input: scratch information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         velocityPool          !< Input/Output: velocity information

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

      real (kind=RKIND), dimension(:), pointer :: &
         thickness, bedTopography, lowerSurface, upperSurface, layerThicknessFractions, betaSolve, sfcMassBal
      real (kind=RKIND), dimension(:,:), pointer :: &
         normalVelocity, uReconstructX, uReconstructY, uReconstructZ
      real (kind=RKIND), dimension(:,:), pointer :: temperature
      real (kind=RKIND), pointer :: deltat
      integer, dimension(:), pointer :: vertexMask, cellMask, edgeMask, floatingEdges
      integer, dimension(:,:), pointer :: dirichletVelocityMask
      character (len=StrKIND), pointer :: config_velocity_solver
      logical, pointer :: config_always_compute_fem_grid
      logical, pointer :: config_output_external_velocity_solver_data
      real (kind=RKIND), pointer :: config_ice_density
      integer, pointer :: anyDynamicVertexMaskChanged
      integer, pointer :: dirichletMaskChanged
      integer, pointer :: nEdges
      type (field2dReal), pointer :: dissipationVertexField
      real (kind=RKIND), dimension(:,:), pointer :: heatDissipation  ! on cells
      integer :: iEdge
      real(kind=RKIND), parameter :: secondsInYear = 365.0 * 24.0 * 3600.0  !< The value of seconds in a year assumed by external dycores

      err = 0

      ! configs
      call mpas_pool_get_config(liConfigs, 'config_velocity_solver', config_velocity_solver)
      call mpas_pool_get_config(liConfigs, 'config_always_compute_fem_grid', config_always_compute_fem_grid)
      call mpas_pool_get_config(liConfigs, 'config_output_external_velocity_solver_data', config_output_external_velocity_solver_data)
      call mpas_pool_get_config(liConfigs, 'config_ice_density', config_ice_density)

      ! Mesh variables
      call mpas_pool_get_array(meshPool, 'layerThicknessFractions', layerThicknessFractions)
      call mpas_pool_get_array(meshPool, 'deltat', deltat)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      ! Geometry variables
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'lowerSurface', lowerSurface)
      call mpas_pool_get_array(geometryPool, 'upperSurface', upperSurface)
      call mpas_pool_get_array(geometryPool, 'vertexMask', vertexMask, timeLevel = 1)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(geometryPool, 'sfcMassBal', sfcMassBal)

      ! Thermal variables
      call mpas_pool_get_array(thermalPool, 'temperature', temperature)
      call mpas_pool_get_array(thermalPool, 'heatDissipation', heatDissipation)

      ! Velocity variables
      call mpas_pool_get_array(velocityPool, 'normalVelocity', normalVelocity)
      call mpas_pool_get_array(velocityPool, 'uReconstructX', uReconstructX)
      call mpas_pool_get_array(velocityPool, 'uReconstructY', uReconstructY)
      call mpas_pool_get_array(velocityPool, 'uReconstructZ', uReconstructZ)
      call mpas_pool_get_array(velocityPool, 'betaSolve', betaSolve)
      call mpas_pool_get_array(velocityPool, 'anyDynamicVertexMaskChanged', anyDynamicVertexMaskChanged)
      call mpas_pool_get_array(velocityPool, 'dirichletMaskChanged', dirichletMaskChanged)
      call mpas_pool_get_array(velocityPool, 'dirichletVelocityMask', dirichletVelocityMask, timeLevel = 1)
      call mpas_pool_get_array(velocityPool, 'floatingEdges', floatingEdges)


      ! ==================================================================
      ! External dycore calls to be made only when vertex mask changes
      ! ==================================================================

      ! Note these functions will always be called on the first solve because we
      ! initialize vertexMask to garbage which sets anyDynamicVertexMaskChanged to 1.
      if ((anyDynamicVertexMaskChanged == 1) .or. (config_always_compute_fem_grid) .or. &
          (dirichletMaskChanged == 1) ) then
         write(stdoutUnit,*) "Generating new external velocity solver FEM grid."
         flush(stdoutUnit)  ! Flush log files before handing control to C++
         flush(stderrUnit)
         call generate_fem_grid(config_velocity_solver, vertexMask, cellMask, dirichletVelocityMask, &
              floatingEdges, layerThicknessFractions, lowerSurface, thickness, err)
      endif


      ! ==================================================================
      ! External dycore calls to be made every time step (solve velocity!)
      ! ==================================================================

      ! convert from m/s (used by MPAS) to m/yr (used by external dycores)
      normalVelocity = normalVelocity * secondsInYear  ! this is intent(out) by dycores, but setting anyway for consistency
      uReconstructX =  uReconstructX  * secondsInYear
      uReconstructY =  uReconstructY  * secondsInYear

      write(stdoutUnit,*) "Beginning velocity solve using external velocity solver."
      flush(stdoutUnit)  ! Flush log files before handing control to C++
      flush(stderrUnit)

      select case (config_velocity_solver)
      case ('L1L2') ! ===============================================
#ifdef USE_EXTERNAL_L1L2
          call mpas_timer_start("velocity_solver_solve_L1L2")
          call velocity_solver_solve_L1L2(lowerSurface, thickness, betaSolve, temperature, &
                uReconstructX, uReconstructY,  &  ! Dirichlet boundary values to apply where dirichletVelocityMask=1
                normalVelocity, uReconstructX, uReconstructY)  ! return values
          call mpas_timer_stop("velocity_solver_solve_L1L2")

          if (config_output_external_velocity_solver_data) then
             ! Optional calls to have LifeV output data files
             call mpas_timer_start("velocity_solver export")
             call velocity_solver_export_2d_data(lowerSurface, thickness, betaSolve)
             call velocity_solver_export_L1L2_velocity();
             call mpas_timer_stop("velocity_solver export")
          endif
#else
              write(stderrUnit,*) "Error: External LifeV library needed to run L1L2 dycore."
              err = 1
              return
#endif

      case ('FO') ! ===============================================
#ifdef USE_EXTERNAL_FIRSTORDER

          ! Allocate scratch var for dissipation on vertices
          ! (Since it is a 3d field, making it scratch)
          call mpas_pool_get_field(scratchPool, 'workLevelVertex', dissipationVertexField)
          call mpas_allocate_scratch_field(dissipationVertexField, .true.)


          call mpas_timer_start("velocity_solver_solve_FO")
          call velocity_solver_solve_FO(bedTopography, lowerSurface, thickness, &
                betaSolve, sfcMassBal, temperature, &
                uReconstructX, uReconstructY,  &  ! Dirichlet boundary values to apply where dirichletVelocityMask=1
                normalVelocity, dissipationVertexField % array, uReconstructX, uReconstructY, deltat)  ! return values
!	  call velocity_solver_estimate_SS_SMB(normalVelocity, mesh % sfcMassBal % array)  ! this was used only for some ice2sea experiments, and is not a general routine to use
          call mpas_timer_stop("velocity_solver_solve_FO")
          ! Now interpolate from vertices to cell centers
          call li_interpolate_vertex_to_cell_2d(meshPool, dissipationVertexField % array, heatDissipation)
          heatDissipation = heatDissipation / (config_ice_density * cp_ice)
          call mpas_deallocate_scratch_field(dissipationVertexField, .true.)

          if (config_output_external_velocity_solver_data) then
             call mpas_timer_start("velocity_solver export")
             call velocity_solver_export_FO_velocity()
             call mpas_timer_stop("velocity_solver export")
          endif
#else
              write(stderrUnit,*) "Error: External library needed to run FO dycore."
              err = 1
              return
#endif

      case ('Stokes') ! ===============================================
#ifdef USE_EXTERNAL_STOKES
          call mpas_timer_start("velocity_solver_solve_stokes")
          call velocity_solver_solve_stokes(lowerSurface, thickness, betaSolve, temperature, &
                uReconstructX, uReconstructY,  &  ! Dirichlet boundary values to apply where dirichletVelocityMask=1
                normalVelocity, uReconstructX, uReconstructY, uReconstructZ)  ! return values
          uReconstructZ = uReconstructZ / (365.0*24.0*3600.0)  ! convert from m/yr to m/s
          call mpas_timer_stop("velocity_solver_solve_stokes")
#else
          write(stderrUnit,*) "Error: External Stokes library needed to run stokes dycore."
          err = 1
          return
#endif
      end select
      write(stdoutUnit,*) "Completed velocity solve using external velocity solver."

      ! convert from m/yr (used by external dycores) to m/s (used by MPAS)
      normalVelocity = normalVelocity / secondsInYear
      uReconstructX =  uReconstructX  / secondsInYear
      uReconstructY =  uReconstructY  / secondsInYear

      ! The external solver will calculate normalVelocity for all edges,
      ! but some of those edges are not dynamically active according to MPASLI's
      ! Voronoi grid conventions.  This zeros velocity on those edges.
      ! (Note: the choice of edges to get reconstructed used to be controlled
      !  inside the interface, but as logic got more complicated with Dirichlet
      !  b.c., that became unwieldly.  Look in the mask routine to see the logic
      !  for which edges are dynamic.)
      do iEdge = 1, nEdges
         if (.not. li_mask_is_dynamic_ice(edgeMask(iEdge)) ) normalVelocity(:,iEdge) = 0.0d0
      end do


   !--------------------------------------------------------------------
   end subroutine li_velocity_external_solve



!***********************************************************************
!
!  routine li_velocity_external_finalize
!
!> \brief   Finalizes external velocity solvers
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine finalizes the ice velocity solver in the external libraries.
!
!-----------------------------------------------------------------------

   subroutine li_velocity_external_finalize(err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

#if defined(USE_EXTERNAL_L1L2) || defined(USE_EXTERNAL_FIRSTORDER) || defined(USE_EXTERNAL_STOKES)
      ! This call is needed for using any of the external velocity solvers
   !   call velocity_solver_finalize()
#else
      write(stderrUnit,*) "Error: To run with an external velocity solver you must compile MPAS with one."
      err = 1
      return
#endif

   !--------------------------------------------------------------------
   end subroutine li_velocity_external_finalize



!***********************************************************************
!  private subroutines
!***********************************************************************



!***********************************************************************
!
!  routine interface_stokes_init
!
!> \brief   Initializes stokes external velocity solver
!> \author Matt Hoffman
!> \date   3 October 2013
!> \details
!>  This routine initializes the ice velocity solver in the stokes
!>  external library (currently only PHG).
!
!-----------------------------------------------------------------------

   subroutine interface_stokes_init(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

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

#ifdef USE_EXTERNAL_STOKES
      ! This call is needed for using any of the PHG velocity solvers
      call phg_init(domain % dminfo % comm)
#else
      write(stderrUnit,*) "Error: External Stokes library needed to run stokes dycore."
      err = 1
      return
#endif

   !--------------------------------------------------------------------
   end subroutine interface_stokes_init



!***********************************************************************
!
!  routine generate_fem_grid
!
!> \brief   Calls to interface to set FEM grid
!> \author Matt Hoffman
!> \date   2 April 2015
!> \details
!>  This routine calls functions in the C interface that generate the FEM grid.
!
!-----------------------------------------------------------------------

   subroutine generate_fem_grid(config_velocity_solver, vertexMask, cellMask, dirichletVelocityMask, floatingEdges, &
              layerThicknessFractions, lowerSurface, thickness, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      character (len=StrKIND), pointer :: config_velocity_solver
      integer, pointer, dimension(:), intent(in) :: vertexMask, cellMask, floatingEdges
      integer, pointer, dimension(:,:), intent(in) :: dirichletVelocityMask
      real(kind=RKIND), pointer, dimension(:), intent(in) :: layerThicknessFractions, &
             lowerSurface, thickness

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      err = 0

#if defined(USE_EXTERNAL_L1L2) || defined(USE_EXTERNAL_FIRSTORDER) || defined(USE_EXTERNAL_STOKES)
          call mpas_timer_start("velocity_solver_compute_2d_grid")
          call velocity_solver_compute_2d_grid(vertexMask, cellMask, dirichletVelocityMask, floatingEdges)
          call mpas_timer_stop("velocity_solver_compute_2d_grid")
#else
          write(stderrUnit,*) "Error: To run with an external velocity solver you must compile MPAS with one."
          err = 1
          return
#endif

          select case (config_velocity_solver)
          case ('L1L2')  ! ===============================================
#ifdef USE_EXTERNAL_L1L2
              call mpas_timer_start("velocity_solver_init_L1L2")
              !call velocity_solver_extrude_3d_grid(layerThicknessFractions, lowerSurface, thickness)
              call velocity_solver_init_L1L2(layerThicknessFractions)
              call mpas_timer_stop("velocity_solver_init_L1L2")
#else
              write(stderrUnit,*) "Error: External LifeV library needed to run L1L2 dycore."
              err = 1
              return
#endif

          case ('FO') ! ===============================================
#ifdef USE_EXTERNAL_FIRSTORDER
              call mpas_timer_start("velocity_solver_extrude_3d_grid")
              call velocity_solver_extrude_3d_grid(layerThicknessFractions, lowerSurface, thickness)
              call mpas_timer_stop("velocity_solver_extrude_3d_grid")
              call mpas_timer_start("velocity_solver_init_FO")
              call velocity_solver_init_FO(layerThicknessFractions)
              call mpas_timer_stop("velocity_solver_init_FO")
#else
              write(stderrUnit,*) "Error: External library needed to run FO dycore."
              err = 1
              return
#endif

          case ('Stokes') ! ===============================================
#ifdef USE_EXTERNAL_STOKES
              call mpas_timer_start("velocity_solver_extrude_3d_grid")
              call velocity_solver_extrude_3d_grid(layerThicknessFractions, lowerSurface, thickness)
              call mpas_timer_stop("velocity_solver_extrude_3d_grid")
              call mpas_timer_start("velocity_solver_init_stokes")
              call velocity_solver_init_stokes(layerThicknessFractions)
              call mpas_timer_stop("velocity_solver_init_stokes")
#else
              write(stderrUnit,*) "Error: External Stokes library needed to run stokes dycore."
              err = 1
              return
#endif
          end select

   !--------------------------------------------------------------------
   end subroutine generate_fem_grid



!***********************************************************************
!
!  routine array_from_exchange_list
!
!> \brief   Converts the MPAS Exchange Lists to flat arrays for external use
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine converts the MPAS Exchange Lists (type mpas_multihalo_exchange_list)
!>  to flat arrays for use by external dycores.
!-----------------------------------------------------------------------

   subroutine array_from_exchange_list(field, sendArray, recvArray)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (field1DInteger), pointer, intent(in) :: field   !< Input: the field that holds the MPAS exchange lists.
      ! Any 1d integer fields will work, but it is suggested to use one of indexToCellID, indexToEdgeID, or IndexToVertexID

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      integer, dimension(:), pointer :: sendArray !< Input/Output: the flat array of elements to send, should be unallocated on input
      integer, dimension(:), pointer :: recvArray !< Input/Output: the flat array of elements to receive, should be unallocated on input

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
!      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (mpas_communication_list), pointer :: sendCommList, recvCommList  ! Communication lists that will be setup from mpas_multihalo_exchange_list's as an intermediate step to flat arrays
      type (mpas_communication_list), pointer :: commListPtr  ! A temporary comm list pointer
      integer :: nHaloLayers, iHalo
      integer, dimension(:), pointer :: haloLayers  ! an array of halo layers, needed to setup comm lists - we want all of them
      type (field1DInteger), pointer :: fieldCursor
      integer :: nAdded, bufferOffset, i
      type (mpas_exchange_list), pointer :: exchListPtr


      ! ========================================================================
      ! Step 1: Generate communication lists from the mpas_multihalo_exchange_list's
      !   (this step is written to be compatible with multiple blocks per processor,
      !    even though that is not supported for external dycores.)
      !   This is done because communication lists have the various halos collapsed
      !   into a single list.
      ! ========================================================================

      ! == First prepare communication lists
      nHaloLayers = size(field % sendList % halos)
      allocate(haloLayers(nHaloLayers))
      do iHalo = 1, nHaloLayers
         haloLayers(iHalo) = iHalo
      end do
      ! Built new send/receive communication lists that have procID & nList filled out.
      call mpas_dmpar_build_comm_lists(field % sendList, field % recvList, haloLayers, field % dimsizes, sendCommList, recvCommList)


      ! == Next populate the commLists' ibuffer field with the element indices to communicate

      ! NOTE:  Looping over the various block's via the field linked list is NOT needed
      !  because packing the communication list with indices will be garbage
      !  if there is more than one block per processor since the indices are block specific.
      !  External dycores currently only support one block per processor and
      !  this subroutine would need substantial modification to support more.
      !  However I am keeping the code that traverses blocks
      !  because this section is taken from mpas_dmpar_exch_halo_field1d_integer
      !  and retaining it makes comparison to that subroutine easier.  The only
      !  difference is the assignements to the ibuffers.
      !  A check for 1 block per proc is in li_velocity_external_init.


      ! Allocate space in send lists, and copy data into buffer
      commListPtr => sendCommList
      do while(associated(commListPtr))  ! Traverse all the processors to be sent to.
        allocate(commListPtr % ibuffer(commListPtr % nList))
        nullify(commListPtr % rbuffer)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0

          fieldCursor => field
          do while(associated(fieldCursor))  ! This is the linked list traversal that is NOT needed.
            exchListPtr => fieldCursor % sendList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if (exchListPtr % endPointID == commListPtr % procID) then
                do  i = 1, exchListPtr % nList
                  commListPtr % ibuffer(exchListPtr % destList(i) + bufferOffset) = exchListPtr % srcList(i)  ! local indices to go into the send communication list
                  nAdded = nAdded + 1

                end do
              end if

              exchListPtr => exchListPtr % next
            end do

            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do

        commListPtr => commListPtr % next
      end do


      ! Allocate space in recv lists, and copy data into buffer
      commListPtr => recvCommList
      do while(associated(commListPtr))  ! Traverse all the processors to receive from.
        allocate(commListPtr % ibuffer(commListPtr % nList))
        nullify(commListPtr % rbuffer)
        bufferOffset = 0
        do iHalo = 1, nHaloLayers
          nAdded = 0
          fieldCursor => field
          do while(associated(fieldCursor))  ! This is the linked list traversal that is NOT needed.
            exchListPtr => fieldCursor % recvList % halos(haloLayers(iHalo)) % exchList
            do while(associated(exchListPtr))
              if (exchListPtr % endPointID == commListPtr % procID) then
                do i = 1, exchListPtr % nList
                  commListPtr % ibuffer( exchListPtr % srcList(i) + bufferOffset ) = exchListPtr % destList(i)  ! buffer index to go into the receive communication list
                end do
                nAdded = max(nAdded, maxval(exchListPtr % srcList))
              end if
              exchListPtr => exchListPtr % next
            end do

            fieldCursor => fieldCursor % next
          end do
          bufferOffset = bufferOffset + nAdded
        end do
        commListPtr => commListPtr % next
      end do


      ! ========================================================================
      ! Step 2: Flatten the communication lists to flat arrays
      ! ========================================================================
      call fill_exchange_array(sendCommList, sendArray)
      call fill_exchange_array(recvCommList, recvArray)


      ! Clean up
      call mpas_dmpar_destroy_communication_list(sendCommList)
      call mpas_dmpar_destroy_communication_list(recvCommList)
      deallocate(haloLayers)

   end subroutine array_from_exchange_list
!***********************************************************************


!***********************************************************************
!
!  routine fill_exchange_array
!
!> \brief   Fills the flat array for external use with information from an MPAS communication list
!> \author Matt Hoffman
!> \date   15 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine converts the MPAS Communication Lists (type mpas_communication_list)
!>  to flat arrays for use by external dycores.  The arrays have this format:
!>
!>  Pos 1: total size of array
!>  For each processor to be communicated with:
!>    Pos 1: processor ID
!>    Pos 2: nList (number of elements in this processor's sub-list
!>    Pos 3 to 3+nList-1:  local indices of elements to be communicated (using 0-based indexing for C/C++)
!-----------------------------------------------------------------------
   subroutine fill_exchange_array(commList, commArray)
      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      type (mpas_communication_list), pointer, intent(in) :: &
         commList  !< Input: Communication list to be flattened

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      integer, dimension(:), pointer :: commArray !< Input/Output: the flat array of elements to communicate, should be unallocated on input

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
!      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      integer :: arraySize ! size of array to house the send or receive list (sendArray, recvArray)
      integer :: offset !  offset for adding metadata about each processor into the flat commArray
      integer :: i
      type (mpas_communication_list), pointer :: commListPtr  ! A temporary comm list pointer for traversing linked lists

      arraySize = 1 !in first position we will store the size of the array
      commListPtr => commList
      do while (associated(commListPtr))
         ! for each processor to be communicated with, we will store the procID, nList, and then the list of local indices to be communicated
         arraySize = arraySize + commListPtr % nlist + 2
         commListPtr => commListPtr % next
      end do

      allocate(commArray(arraySize))

      commArray(1) = arraySize
      offset = 2    ! we will store the procID, nList before the list of local indices
      commListPtr => commList
      do while (associated(commListPtr))
         commArray(offset) = commListPtr % procID  ! store procID
         offset = offset + 1
         commArray(offset) = commListPtr % nlist   ! store nList
         do i = 1 , commListPtr % nlist
            commArray(i+offset) = commListPtr % ibuffer(i) -1  ! add the list of elements to be communicated, switching to 0-based indexing for C/C++
         end do
         offset = offset + commListPtr % nlist + 1

         commListPtr => commListPtr % next
      end do

   end subroutine fill_exchange_array
!***********************************************************************


end module li_velocity_external

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
