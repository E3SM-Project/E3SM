!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_mesh_pool
!
!> \brief
!> \date 2020
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_mesh_pool

#include "gpu_macros.inc"

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_log

   implicit none
   private

   integer, public :: &
      nCells,         &
      nVerticesSolve, &
      vertexDegree

   integer, public, dimension(:), pointer :: & 
      nEdgesOnCell,   &
      solveStress,    &
      solveVelocity

   integer, public, dimension(:,:), pointer :: &
      verticesOnCell, &
      cellsOnVertex,  &
      cellVerticesAtVertex

   real(kind=RKIND), public, dimension(:), pointer :: &
      areaTriangle,   &
      tanLatVertexRotatedOverRadius, &
      icePressure,    &
      uVelocity,      &
      vVelocity,      &
      stressDivergenceU, &
      stressDivergenceV

   real(kind=RKIND), public, dimension(:,:), pointer :: &
      stress11,       &
      stress12,       &
      stress22

   real(kind=RKIND), public, dimension(:,:,:), pointer :: &
      basisGradientU, &
      basisGradientV, &
      basisIntegralsU,&
      basisIntegralsV,&
      basisIntegralsMetric

   public ::                   &
      seaice_mesh_pool_create, &
      seaice_mesh_pool_update, &
      seaice_mesh_pool_destroy

!-----------------------------------------------------------------------

contains



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_mesh_pool_create
!
!-----------------------------------------------------------------------

  subroutine seaice_mesh_pool_create(&
       domain)!{{{

    type(domain_type) :: &
         domain

    integer :: &
         blockCount

    type(block_type), pointer :: &
         block

    type (mpas_pool_type), pointer :: &
         meshPool,               &
         velocitySolverPool,     &
         velocityVariationalPool

    integer, pointer ::   &
         nCellsTmp,         &
         nVerticesSolveTmp, &
         vertexDegreeTmp

    blockCount = 0
    block => domain % blocklist
    do while ( associated(block) )

       blockCount = blockCount + 1
       if (blockCount > 1) then
          call mpas_log_write('seaice_mesh_pool_create: more than one block is no longer supported', MPAS_LOG_CRIT)
       endif

       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_subpool(block % structs, 'velocity_solver', velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)

       ! convert mesh dimensions from pointers to scalars
       call mpas_pool_get_dimension(meshPool, 'nCells', nCellsTmp)
       call MPAS_pool_get_dimension(meshPool, "nVerticesSolve", nVerticesSolveTmp)
       call MPAS_pool_get_dimension(meshPool, "vertexDegree", vertexDegreeTmp)

       nCells         = nCellsTmp
       nVerticesSolve = nVerticesSolveTmp
       vertexDegree   = vertexDegreeTmp

       ! point to existing arrays
       call mpas_pool_get_array(meshPool, 'nEdgesOnCell',   nEdgesOnCell)
       call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_array(meshPool, 'cellsOnVertex',  cellsOnVertex)
       call mpas_pool_get_array(meshPool, 'areaTriangle',   areaTriangle)

       call MPAS_pool_get_array(velocitySolverPool, "solveStress",       solveStress)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity",     solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "icePressure",       icePressure)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity",         uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity",         vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)

       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientU",       basisGradientU)
       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientV",       basisGradientV)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsU",      basisIntegralsU)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsV",      basisIntegralsV)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsMetric", basisIntegralsMetric)
       call MPAS_pool_get_array(velocityVariationalPool, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
       call MPAS_pool_get_array(velocityVariationalPool, "cellVerticesAtVertex", cellVerticesAtVertex)
       call MPAS_pool_get_array(velocityVariationalPool, "stress11",             stress11)
       call MPAS_pool_get_array(velocityVariationalPool, "stress22",             stress22)
       call MPAS_pool_get_array(velocityVariationalPool, "stress12",             stress12)

#if defined(MPAS_OPENMP_OFFLOAD) || defined(MPAS_OPENACC)
!$GPU ENTER_DATA COPY_IN_LP             &
!$GPUC   nCells,                        &
!$GPUC   nVerticesSolve,                &
!$GPUC   vertexDegree,                  &
!$GPUC   nEdgesOnCell,                  &
!$GPUC   verticesOnCell,                &
!$GPUC   cellsOnVertex,                 &
!$GPUC   areaTriangle,                  &
!$GPUC   solveStress,                   &
!$GPUC   solveVelocity,                 &
!$GPUC   icePressure,                   &
!$GPUC   uVelocity,                     &
!$GPUC   vVelocity,                     &
!$GPUC   stressDivergenceU,             &
!$GPUC   stressDivergenceV,             &
!$GPUC   basisGradientU,                &
!$GPUC   basisGradientV,                &
!$GPUC   basisIntegralsU,               &
!$GPUC   basisIntegralsV,               &
!$GPUC   basisIntegralsMetric,          &
!$GPUC   tanLatVertexRotatedOverRadius, &
!$GPUC   cellVerticesAtVertex,          &
!$GPUC   stress11,                      &
!$GPUC   stress12,                      &
!$GPUC   stress22                       &
!$GPUF
#endif

       block => block % next
    end do

  end subroutine seaice_mesh_pool_create!}}}
!-----------------------------------------------------------------------



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_mesh_pool_destroy
!
!-----------------------------------------------------------------------

  subroutine seaice_mesh_pool_destroy(&
       err)!{{{

    integer, intent(out) :: &
      err                   ! returned error flag

    err = 0

#if defined(MPAS_OPENMP_OFFLOAD) || defined(MPAS_OPENACC)
    ! first delete on device
!$GPU EXIT_DATA COPY_DEL_LP             &
!$GPUC   nCells,                        &
!$GPUC   nVerticesSolve,                &
!$GPUC   vertexDegree,                  &
!$GPUC   nEdgesOnCell,                  &
!$GPUC   verticesOnCell,                &
!$GPUC   cellsOnVertex,                 &
!$GPUC   areaTriangle,                  &
!$GPUC   solveStress,                   &
!$GPUC   solveVelocity,                 &
!$GPUC   icePressure,                   &
!$GPUC   uVelocity,                     &
!$GPUC   vVelocity,                     &
!$GPUC   stressDivergenceU,             &
!$GPUC   stressDivergenceV,             &
!$GPUC   basisGradientU,                &
!$GPUC   basisGradientV,                &
!$GPUC   basisIntegralsU,               &
!$GPUC   basisIntegralsV,               &
!$GPUC   basisIntegralsMetric,          &
!$GPUC   tanLatVertexRotatedOverRadius, &
!$GPUC   cellVerticesAtVertex,          &
!$GPUC   stress11,                      &
!$GPUC   stress12,                      &
!$GPUC   stress22                       &
!$GPUF
#endif

    ! then nullify on host
    nullify(nEdgesOnCell,      &
         verticesOnCell,       &
         cellsOnVertex,        &
         areaTriangle,         &
         solveStress,          &
         solveVelocity,        &
         icePressure,          &
         uVelocity,            &
         vVelocity,            &
         stressDivergenceU,    &
         stressDivergenceV,    &
         basisGradientU,       &
         basisGradientV,       &
         basisIntegralsU,      &
         basisIntegralsV,      &
         basisIntegralsMetric, &
         tanLatVertexRotatedOverRadius, &
         cellVerticesAtVertex, &
         stress11,             &
         stress12,             &
         stress22              &
    )

  end subroutine seaice_mesh_pool_destroy!}}}
!-----------------------------------------------------------------------



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_mesh_pool_update
!
!-----------------------------------------------------------------------

  subroutine seaice_mesh_pool_update(&
       domain)!{{{

    type(domain_type) :: &
         domain

#if defined(MPAS_OPENMP_OFFLOAD) || defined(MPAS_OPENACC)
    ! update arrays on device
!$GPU UPDATE_D_LP       &
!$GPUC   solveStress,   &
!$GPUC   solveVelocity, &
!$GPUC   icePressure,   &
!$GPUC   uVelocity,     &
!$GPUC   vVelocity,     &
!$GPUC   stress11,      &
!$GPUC   stress12,      &
!$GPUC   stress22       &
!$GPUF
#endif

  end subroutine seaice_mesh_pool_update!}}}
!-----------------------------------------------------------------------



end module seaice_mesh_pool
