! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_init
!
!> \brief MPAS-ocean initialization routines.
!> \author Mark Petersen
!> \date   December 2013
!> \details
!>  This module contains routines to initialize variables at the
!>    beginning of an MPAS-Ocean simulation, or when starting the
!>    ocean analysis core.
!
!-----------------------------------------------------------------------

module ocn_init_routines

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_dmpar
   use mpas_constants

   use mpas_rbf_interpolation
   use mpas_vector_operations
   use mpas_vector_reconstruction

   use ocn_diagnostics
   use ocn_diagnostics_variables
   use ocn_gm
   use ocn_constants
   use ocn_config
   use ocn_mesh

   use ocn_surface_land_ice_fluxes
   use ocn_forcing

   interface ocn_init_add_tau_metadata
      module procedure ocn_init_add_tau_metadata_real
      module procedure ocn_init_add_tau_metadata_int
      module procedure ocn_init_add_tau_metadata_logical
      module procedure ocn_init_add_tau_metadata_character
   end interface

   private

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

   public :: &
      ocn_init_routines_block, &
      ocn_init_metadata

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_routines_block
!
!> \brief   Initialize blocks within MPAS-Ocean core
!> \author  Doug Jacobsen, Mark Petersen, Todd Ringler
!> \date    September 2011
!> \details
!>  This routine calls all block-level initializations required to begin a
!>  simulation with MPAS-Ocean
!
!-----------------------------------------------------------------------

   subroutine ocn_init_routines_block(block, dt, err)!{{{

      type (block_type), intent(inout) :: block
      real (kind=RKIND), intent(in) :: dt
      integer, intent(out) :: err

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool, tracersPool
      type (mpas_pool_type), pointer :: forcingPool, scratchPool
      integer :: iEdge, iCell
      integer :: err1

      integer, dimension(:,:), pointer :: boundaryCellTmp
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity

      real (kind=RKIND), dimension(:,:,:), pointer :: tracersGroup

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels

      real (kind=RKIND), dimension(:), pointer :: &
        seaIcePressure, atmosphericPressure

      real (kind=RKIND), dimension(:), pointer :: &
        frazilSurfacePressure

      real (kind=RKIND), dimension(:), pointer :: ssh

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers

      type (mpas_pool_iterator_type) :: groupItr

      call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
      call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
      call mpas_pool_get_dimension(block % dimensions, 'nVertices', nVertices)
      call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)

      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'state', statePool)
      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)

      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

      call mpas_pool_get_array(meshPool, 'boundaryCell', boundaryCellTmp)

      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

      call mpas_pool_get_array(forcingPool, 'seaIcePressure', seaIcePressure)
      call mpas_pool_get_array(forcingPool, 'atmosphericPressure', atmosphericPressure)
      call mpas_pool_get_array(forcingPool, 'frazilSurfacePressure', frazilSurfacePressure)

      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)

      if (.not. config_do_restart) then
         do iCell=1,nCells
            boundaryLayerDepth(iCell) = layerThickness(minLevelCell(iCell), iCell) * 0.5_RKIND
            indMLD(iCell) = minLevelCell(iCell)
         end do
      end if

#ifdef MPAS_OPENACC
      !$acc enter data copyin(layerThickness, normalVelocity)
      !$acc update device (normalTransportVelocity, &
      !$acc                normalGMBolusVelocity)
      !$acc enter data copyin(atmosphericPressure, seaIcePressure)
      !$acc enter data copyin(ssh)
      !$acc enter data copyin(activeTracers)
      !$acc update device(tracersSurfaceValue)
      if ( associated(frazilSurfacePressure) ) then
         !$acc enter data copyin(frazilSurfacePressure)
      endif
      if (landIcePressureOn) then
         !$acc enter data copyin(landIcePressure)
         !$acc enter data copyin(landIceDraft)
      endif
#endif
      call ocn_diagnostic_solve(dt,  statePool, forcingPool, meshPool, scratchPool, tracersPool)
#ifdef MPAS_OPENACC
      !$acc update host(layerThickEdgeFlux, layerThickEdgeMean)
      !$acc update host(relativeVorticity, circulation)
      !$acc update host(vertTransportVelocityTop, &
      !$acc             vertGMBolusVelocityTop, &
      !$acc             relativeVorticityCell, &
      !$acc             divergence, &
      !$acc             kineticEnergyCell, &
      !$acc             tangentialVelocity, &
      !$acc             vertVelocityTop)
      !$acc update host(normRelVortEdge, normPlanetVortEdge, &
      !$acc             normalizedRelativeVorticityCell)
      !$acc update host (surfacePressure)
      !$acc update host(zMid, zTop)
      !$acc exit data copyout(ssh)
      !$acc exit data delete(activeTracers)
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
      !$acc exit data delete(layerThickness, normalVelocity)
#endif

      layerThickness(:, nCells+1) = 0.0_RKIND


      do iEdge=1, nEdges
         normalVelocity(maxLevelEdgeTop(iEdge)+1:maxLevelEdgeBot(iEdge), iEdge) = 0.0_RKIND

         normalVelocity(maxLevelEdgeBot(iEdge)+1:nVertLevels,iEdge) = -1.0e34_RKIND
         
         normalVelocity(1:minLevelEdgeTop(iEdge)-1,iEdge) = -1.0e34_RKIND
      end do

      call mpas_pool_begin_iteration(tracersPool)
      do while ( mpas_pool_get_next_member(tracersPool, groupItr) )
         if ( groupItr % memberType == MPAS_POOL_FIELD ) then
            call mpas_pool_get_array(tracersPool, groupItr % memberName, tracersGroup, 1)
            if ( associated(tracersGroup) ) then
               do iCell=1,nCells
                  tracersGroup(:, maxLevelCell(iCell)+1:nVertLevels,iCell) =  -1.0e34_RKIND
                  tracersGroup(:, 1:minLevelCell(iCell)-1,iCell) =  -1.0e34_RKIND
               end do
            end if
         end if
      end do

      ! ------------------------------------------------------------------
      ! Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------
      do iEdge = 1, nEdges
         normalTransportVelocity(:, iEdge) = normalVelocity(:, iEdge)
      end do

      ! ------------------------------------------------------------------
      ! End: Accumulating various parametrizations of the transport velocity
      ! ------------------------------------------------------------------

      call mpas_rbf_interp_initialize(meshPool)
      call mpas_initialize_tangent_vectors(meshPool, edgeTangentVectors)

      call mpas_init_reconstruct(meshPool, includeHalos=.true.)

      call mpas_reconstruct(meshPool, normalVelocity,        &
                       velocityX,            &
                       velocityY,            &
                       velocityZ,            &
                       velocityZonal,        &
                       velocityMeridional    &
                      )

      call mpas_pool_initialize_time_levels(statePool)

      ! compute land-ice fluxes for potential output at startup
      call ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, err1, 1)
      err = ior(err, err1)
      call ocn_surface_land_ice_fluxes_build_arrays(meshPool, &
                                                    forcingPool, scratchPool, statePool, dt, err1)
      err = ior(err, err1)



   end subroutine ocn_init_routines_block!}}}


!***********************************************************************
!
!  routine ocn_init_metadata
!
!> \brief   Initialize any metadata for this processor
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine sets up any metadata for this MPI task and it's associated threads.
!>  The meta data could be related to performance data, or information about
!>  all the blocks on this processor.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_metadata(domain)!{{{
      type (domain_type), intent(inout) :: domain

      type (block_type), pointer :: block

      character (len=StrKIND) :: metaDataName

      integer :: iHalo
      integer :: numBlocks
      integer, dimension(:), pointer :: nCellsArray, nEdgesArray, nVerticesArray

      numBlocks = 0

      block => domain % blocklist
      do while ( associated(block) )
         numBlocks = numBlocks + 1

         call mpas_pool_get_dimension(block % dimensions, 'nCellsArray', nCellsArray)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesArray', nEdgesArray)
         call mpas_pool_get_dimension(block % dimensions, 'nVerticesArray', nVerticesArray)

         call ocn_init_add_tau_metadata( 'numCellHalos', size(nCellsArray) )
         call ocn_init_add_tau_metadata( 'numEdgeHalos', size(nEdgesArray) )
         call ocn_init_add_tau_metadata( 'numVertexHalos', size(nVerticesArray) )

         do iHalo = 1, size(nCellsArray)
            write(metaDataName, '(a8, i2)' ) 'cellHalo', iHalo
            call ocn_init_add_tau_metadata( metaDataName, nCellsArray(iHalo) )
         end do

         do iHalo = 1, size(nEdgesArray)
            write(metaDataName, '(a8, i2)' ) 'edgeHalo', iHalo
            call ocn_init_add_tau_metadata( metaDataName, nEdgesArray(iHalo) )
         end do

         do iHalo = 1, size(nVerticesArray)
            write(metaDataName, '(a10, i2)' ) 'vertexHalo', iHalo
            call ocn_init_add_tau_metadata( metaDataName, nVerticesArray(iHalo) )
         end do

         block => block % next
      end do

      call ocn_init_add_tau_metadata( 'numBlocks', numBlocks )

   end subroutine ocn_init_metadata!}}}


!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_real
!
!> \brief   Add tau real metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds a real value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_real( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      real (kind=RKIND), intent(in) :: dataValue

      character (len=StrKIND) :: dataString

#ifdef MPAS_TAU
      dataString = ''
      write( dataString, * ) dataValue
      call tau_metadata( trim(dataName), trim(dataString) )
#endif

   end subroutine ocn_init_add_tau_metadata_real!}}}

!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_int
!
!> \brief   Add tau integer metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds an integer value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_int( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      integer, intent(in) :: dataValue

      character (len=StrKIND) :: dataString

#ifdef MPAS_TAU
      dataString = ''
      write( dataString, * ) dataValue
      call tau_metadata( trim(dataName), trim(dataString) )
#endif

   end subroutine ocn_init_add_tau_metadata_int!}}}

!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_logical
!
!> \brief   Add tau logical metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds a logical value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_logical( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      logical, intent(in) :: dataValue

      character (len=StrKIND) :: dataString

#ifdef MPAS_TAU
      dataString = ''
      write( dataString, * ) dataValue
      call tau_metadata( trim(dataName), trim(dataString) )
#endif

   end subroutine ocn_init_add_tau_metadata_logical!}}}

!***********************************************************************
!
!  routine ocn_init_add_tau_metadata_character
!
!> \brief   Add tau character metadata
!> \author  Doug Jacobsen
!> \date    08/05/2016
!> \details
!>  This routine adds a character value metadata for TAU to this task.
!
!-----------------------------------------------------------------------
   subroutine ocn_init_add_tau_metadata_character( dataName, dataValue )!{{{
      character (len=*), intent(in) :: dataName
      character (len=*), intent(in) :: dataValue

#ifdef MPAS_TAU
      call tau_metadata( trim(dataName), trim(dataValue) )
#endif

   end subroutine ocn_init_add_tau_metadata_character!}}}


!***********************************************************************

end module ocn_init_routines

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
