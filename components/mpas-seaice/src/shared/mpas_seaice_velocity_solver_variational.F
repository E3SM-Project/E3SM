!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver_variational

#include "gpu_macros.inc"

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer
  use mpas_log, only: mpas_log_write

  implicit none

  private
  save

  public :: &
       seaice_init_velocity_solver_variational, &
       seaice_strain_tensor_variational, &
       seaice_average_strains_on_vertex, &
       seaice_stress_tensor_variational, &
       seaice_stress_divergence_variational, &
       seaice_final_divergence_shear_variational

contains

!-----------------------------------------------------------------------
! initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_velocity_solver_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 24 October 2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_velocity_solver_variational(domain)!{{{

    use seaice_mesh, only: &
         seaice_calc_local_coords

    use seaice_velocity_solver_variational_shared, only: &
         seaice_calc_variational_metric_terms, &
         seaice_cell_vertices_at_vertex

    use seaice_velocity_solver_wachspress, only: &
         seaice_init_velocity_solver_wachspress

    use seaice_velocity_solver_pwl, only: &
         seaice_init_velocity_solver_pwl

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    character(len=strKIND), pointer :: &
         config_variational_basis

    call mpas_timer_start("variational calc_metric_terms")
    call seaice_calc_variational_metric_terms(domain)
    call mpas_timer_stop("variational calc_metric_terms")

    call mpas_timer_start("variational vertices_at_vertex")
    call seaice_cell_vertices_at_vertex(domain)
    call mpas_timer_stop("variational vertices_at_vertex")

    call MPAS_pool_get_config(domain % configs, "config_variational_basis", config_variational_basis)

    if (trim(config_variational_basis) == "wachspress") then

       call seaice_init_velocity_solver_wachspress(domain)

    else if (trim(config_variational_basis) == "pwl") then

       call seaice_init_velocity_solver_pwl(domain)

    else if (trim(config_variational_basis) == "none") then

       continue

    else

       call MPAS_log_write("Unknown variational basis type: "//trim(config_variational_basis), MPAS_LOG_CRIT)

    endif

    call mpas_timer_start("variational denominator")
    call variational_denominator(domain)
    call mpas_timer_stop("variational denominator")

    !call homogenize_variational_basis_field()

  end subroutine seaice_init_velocity_solver_variational

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  variational_denominator
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 30th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine variational_denominator(domain)

    type (domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         boundaryPool, &
         velocityVariationalPool

    character(len=strKIND), pointer :: &
         config_variational_denominator_type

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:), pointer :: &
         nEdgesOnCell, &
         interiorVertex

    real(kind=RKIND), dimension(:), pointer :: &
         areaTriangle

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisIntegralsMetric

    real(kind=RKIND), dimension(:), pointer :: &
         variationalDenominator

    integer :: &
         iVertex, &
         iSurroundingCell, &
         iStressVertex, &
         iCell, &
         iVelocityVertex

    call MPAS_pool_get_config(domain % configs, "config_variational_denominator_type", config_variational_denominator_type)
    if (trim(config_variational_denominator_type) == "alternate") then

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "boundary", boundaryPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

          call MPAS_pool_get_dimension(meshPool, "nVertices", nVertices)
          call MPAS_pool_get_dimension(meshPool, "vertexDegree", vertexDegree)

          call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)
          call MPAS_pool_get_array(meshPool, "areaTriangle", areaTriangle)
          call MPAS_pool_get_array(meshPool, "cellsOnVertex", cellsOnVertex)

          call MPAS_pool_get_array(boundaryPool, "interiorVertex", interiorVertex)

          call MPAS_pool_get_array(velocityVariationalPool, "cellVerticesAtVertex", cellVerticesAtVertex)
          call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsMetric", basisIntegralsMetric)
          call MPAS_pool_get_array(velocityVariationalPool, "variationalDenominator", variationalDenominator)

          do iVertex = 1, nVertices

             variationalDenominator(iVertex) = 0.0_RKIND

             ! loop over surrounding cells
             do iSurroundingCell = 1, vertexDegree

                ! get the cell number of this cell
                iCell = cellsOnVertex(iSurroundingCell, iVertex)

                ! get the vertexOnCell number of the iVertex velocity point from cell iCell
                iVelocityVertex = cellVerticesAtVertex(iSurroundingCell,iVertex)

                ! loop over the vertices of the surrounding cell
                do iStressVertex = 1, nEdgesOnCell(iCell)

                   variationalDenominator(iVertex) = variationalDenominator(iVertex) + &
                        basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell)

                enddo ! iStressVertex

             enddo ! iSurroundingCell

             ! inverse
             variationalDenominator(iVertex) = variationalDenominator(iVertex)

          enddo ! iVertex

          blockPtr => blockPtr % next
       enddo

    else if (trim(config_variational_denominator_type) == "original") then

       blockPtr => domain % blocklist
       do while (associated(blockPtr))

          call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
          call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

          call MPAS_pool_get_dimension(meshPool, "nVertices", nVertices)

          call MPAS_pool_get_array(meshPool, "areaTriangle", areaTriangle)

          call MPAS_pool_get_array(velocityVariationalPool, "variationalDenominator", variationalDenominator)

          do iVertex = 1, nVertices
             variationalDenominator(iVertex) = areaTriangle(iVertex)
          enddo ! iVertex

          blockPtr => blockPtr % next
       enddo

    else

       call MPAS_log_write("Unknown variational denominator type: "//trim(config_variational_denominator_type), MPAS_LOG_CRIT)

    endif

  end subroutine variational_denominator

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  homogenize_variational_basis_field
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 9th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine homogenize_variational_basis_field()

    use seaice_mesh_pool, only: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsMetric, &
         basisIntegralsU, &
         basisIntegralsV, &
         nCells

    integer :: &
         iCell

    integer, parameter :: &
         iCellHomogenize = 1111

    !call homogenize_cell(basisGradientU, iCellHomogenize)
    !call homogenize_cell(basisGradientV, iCellHomogenize)
    !call homogenize_cell(basisIntegralsMetric, iCellHomogenize)
    !call homogenize_cell(basisIntegralsU, iCellHomogenize)
    !call homogenize_cell(basisIntegralsV, iCellHomogenize)

    do iCell = 1, nCells

       basisGradientU(:,:,iCell)       = basisGradientU(:,:,iCellHomogenize)
       basisGradientV(:,:,iCell)       = basisGradientV(:,:,iCellHomogenize)
       basisIntegralsMetric(:,:,iCell) = basisIntegralsMetric(:,:,iCellHomogenize)
       basisIntegralsU(:,:,iCell)      = basisIntegralsU(:,:,iCellHomogenize)
       basisIntegralsV(:,:,iCell)      = basisIntegralsV(:,:,iCellHomogenize)

    enddo ! iCell

  end subroutine homogenize_variational_basis_field

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  homogenize_variational_basis_field
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 9th January 2021
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine homogenize_cell(&
       field, &
       iCell)

    real(kind=RKIND), dimension(:,:,:), intent(inout) :: &
         field

    integer, intent(in) :: &
         iCell

    integer :: &
         iVertexOnCell1, &
         iVertexOnCell2, &
         nCanonicalValues, &
         iCV

    real(kind=RKIND), dimension(100) :: &
         canonicalValues

    logical :: &
         lFound

    nCanonicalValues = 0

    do iVertexOnCell1 = 1, 6
       do iVertexOnCell2 = 1, 6

          lFound = .false.
          do iCV = 1, nCanonicalValues

             if (abs(abs(canonicalValues(iCV) - abs(field(iVertexOnCell1,iVertexOnCell2,iCell)))) < 1e-12) then
                lFound = .true.
                field(iVertexOnCell1,iVertexOnCell2,iCell) = canonicalValues(iCV) * sign(1.0_RKIND, field(iVertexOnCell1,iVertexOnCell2,iCell))
                exit
             endif

          enddo ! iCV

          if (.not. lFound) then
             nCanonicalValues = nCanonicalValues + 1
             canonicalValues(nCanonicalValues) = abs(field(iVertexOnCell1,iVertexOnCell2,iCell))
          endif

       enddo ! iVertexOnCell2
    enddo ! iVertexOnCell1

    !do iVertexOnCell1 = 1, 6
    !   do iVertexOnCell2 = 1, 6
    !      write(*,*) iVertexOnCell1, iVertexOnCell2, field(iVertexOnCell1,iVertexOnCell2,iCell)
    !   enddo ! iVertexOnCell2
    !enddo ! iVertexOnCell1

  end subroutine homogenize_cell

!-----------------------------------------------------------------------
! time step
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_strain_tensor_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_strain_tensor_variational(domain)!{{{

    use seaice_mesh_pool, only: &
         nCells, &
         verticesOnCell, &
         nEdgesOnCell, &
         basisGradientU, &
         basisGradientV, &
         solveStress, &
         tanLatVertexRotatedOverRadius, &
         uVelocity, &
         vVelocity

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         velocityVariationalPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12

    integer :: &
         iCell, &
         iGradientVertex, &
         iBasisVertex, &
         iVertex, &
         jVertex

    real(kind=RKIND) :: &
         strain11Tmp, &
         strain22Tmp, &
         strain12Tmp

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

       call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11)
       call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22)
       call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12)

       ! loop over cells
#ifdef MPAS_OPENMP_OFFLOAD
!$omp target teams distribute parallel do
#elif MPAS_OPENACC
!$acc parallel loop gang worker
#else
!$omp parallel do default(shared) private(iGradientVertex, iBasisVertex, iVertex, jVertex, &
!$omp&                                    strain11Tmp, strain22Tmp, strain12Tmp)
#endif
       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             ! loop over velocity points surrounding cell - location of stress and derivative
             do iGradientVertex = 1, nEdgesOnCell(iCell)

                strain11Tmp = 0.0_RKIND
                strain22Tmp = 0.0_RKIND
                strain12Tmp = 0.0_RKIND

                ! loop over basis functions
                do iBasisVertex = 1, nEdgesOnCell(iCell)

                   iVertex = verticesOnCell(iBasisVertex,iCell)

                   strain11Tmp = strain11Tmp + uVelocity(iVertex) * basisGradientU(iBasisVertex,iGradientVertex,iCell)
                   strain22Tmp = strain22Tmp + vVelocity(iVertex) * basisGradientV(iBasisVertex,iGradientVertex,iCell)
                   strain12Tmp = strain12Tmp + 0.5_RKIND * (&
                        uVelocity(iVertex) * basisGradientV(iBasisVertex,iGradientVertex,iCell) + &
                        vVelocity(iVertex) * basisGradientU(iBasisVertex,iGradientVertex,iCell))

                enddo ! iVertexOnCell

                ! metric terms
                jVertex = verticesOnCell(iGradientVertex,iCell)

                strain11(iGradientVertex,iCell) = strain11Tmp - vVelocity(jVertex) * tanLatVertexRotatedOverRadius(jVertex)
                strain12(iGradientVertex,iCell) = strain12Tmp + uVelocity(jVertex) * tanLatVertexRotatedOverRadius(jVertex) * 0.5_RKIND
                strain22(iGradientVertex,iCell) = strain22Tmp

             enddo ! jVertexOnCell

          endif ! solveStress

       enddo ! iCell

       blockPtr => blockPtr % next
    end do

  end subroutine seaice_strain_tensor_variational!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_average_strains_on_vertex
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_average_strains_on_vertex(domain)

    use seaice_mesh_pool, only: &
         nCells, &
         nVerticesSolve, &
         cellVerticesAtVertex, &
         cellsOnVertex, &
         vertexDegree

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocityVariationalPool

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND) :: &
         strain11avg, &
         strain22avg, &
         strain12avg, &
         denominator

    integer :: &
         iVertex, &
         iVertexDegree, &
         iCell, &
         iVertexOnCell

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

       call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11)
       call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22)
       call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12)

       do iVertex = 1, nVerticesSolve

          strain11avg = 0.0_RKIND
          strain22avg = 0.0_RKIND
          strain12avg = 0.0_RKIND
          denominator = 0.0_RKIND

          do iVertexDegree = 1, vertexDegree

             iCell = cellsOnVertex(iVertexDegree,iVertex)

             if (iCell <= nCells) then

                iVertexOnCell = cellVerticesAtVertex(iVertexDegree,iVertex)

                strain11avg = strain11avg + strain11(iVertexOnCell,iCell) * areaCell(iCell)
                strain22avg = strain22avg + strain22(iVertexOnCell,iCell) * areaCell(iCell)
                strain12avg = strain12avg + strain12(iVertexOnCell,iCell) * areaCell(iCell)
                denominator = denominator + areaCell(iCell)

             endif

          enddo ! iVertexDegree

          strain11avg = strain11avg / denominator
          strain22avg = strain22avg / denominator
          strain12avg = strain12avg / denominator

          do iVertexDegree = 1, vertexDegree

             iCell = cellsOnVertex(iVertexDegree,iVertex)

             if (iCell <= nCells) then

                iVertexOnCell = cellVerticesAtVertex(iVertexDegree,iVertex)

                strain11(iVertexOnCell,iCell) = strain11avg
                strain22(iVertexOnCell,iCell) = strain22avg
                strain12(iVertexOnCell,iCell) = strain12avg

             endif

          enddo ! iCellOnVertex

       enddo ! iVertex

       blockPtr => blockPtr % next
    end do

  end subroutine seaice_average_strains_on_vertex

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_stress_tensor_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_stress_tensor_variational(domain)!{{{

    use seaice_velocity_solver_constitutive_relation, only: &
         constitutiveRelationType, &
         EVP_CONSTITUTIVE_RELATION, &
         REVISED_EVP_CONSTITUTIVE_RELATION, &
         LINEAR_CONSTITUTIVE_RELATION, &
         seaice_evp_constitutive_relation, &
         seaice_evp_constitutive_relation_revised, &
         seaice_linear_constitutive_relation, &
         eccentricitySquared, puny, dampingTimescale

    use seaice_mesh_pool, only: &
         nCells, &
         nEdgesOnCell, &
         stress11, &
         stress22, &
         stress12, &
         icePressure, &
         solveStress

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocitySolverPool, &
         velocityVariationalPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         replacementPressure

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), pointer :: &
         elasticTimeStep

    integer :: &
         iCell, &
         iVertexOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    real(kind=RKIND) :: &
         strainDivergence,    &
         strainTension,       &
         strainShearing,      &
         stress1,             &
         stress2,             &
         Delta,               &
         pressureCoefficient, &
         denominator

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

       call MPAS_pool_get_array(velocitySolverPool, "elasticTimeStep", elasticTimeStep)

       call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11)
       call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22)
       call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12)
       call MPAS_pool_get_array(velocityVariationalPool, "replacementPressure", replacementPressure)

       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

       if (constitutiveRelationType == EVP_CONSTITUTIVE_RELATION) then

          denominator = 1.0_RKIND + (0.5_RKIND * elasticTimeStep) / dampingTimescale

#ifdef MPAS_OPENMP_OFFLOAD
!$omp target teams distribute parallel do
#elif MPAS_OPENACC
!$acc parallel loop gang worker
#else
!$omp parallel do default(shared) private(iVertexOnCell)
#endif
          do iCell = 1, nCells

             replacementPressure(:,iCell) = 0.0_RKIND

             if (solveStress(iCell) == 1) then

#if defined(MPAS_OPENMP_OFFLOAD) || defined(MPAS_OPENACC)
                ! inline call to seaice_evp_constitutive_relation for GPUs
                do iVertexOnCell = 1, nEdgesOnCell(iCell)

                   ! convert from stress11 to stress1 etc
                   strainDivergence = strain11(iVertexOnCell,iCell) + strain22(iVertexOnCell,iCell)
                   strainTension    = strain11(iVertexOnCell,iCell) - strain22(iVertexOnCell,iCell)
                   strainShearing   = strain12(iVertexOnCell,iCell) * 2.0_RKIND

                   stress1 = stress11(iVertexOnCell,iCell) + stress22(iVertexOnCell,iCell)
                   stress2 = stress11(iVertexOnCell,iCell) - stress22(iVertexOnCell,iCell)

                   ! perform the constituitive relation
                   Delta = sqrt(strainDivergence*strainDivergence + &
                               (strainTension*strainTension + strainShearing*strainShearing) / eccentricitySquared)

                   pressureCoefficient                      = icePressure(iCell) / max(Delta,puny)
                   replacementPressure(iVertexOnCell,iCell) = pressureCoefficient * Delta

                   pressureCoefficient = (pressureCoefficient * elasticTimeStep) / (2.0_RKIND * dampingTimescale)

                   stress1  = (stress1  +  pressureCoefficient                        * (strainDivergence - Delta)) / denominator
                   stress2  = (stress2  + (pressureCoefficient / eccentricitySquared) *  strainTension            ) / denominator
                   stress12(iVertexOnCell,iCell) = (stress12(iVertexOnCell,iCell) &
                        + (pressureCoefficient / eccentricitysquared) * strainShearing * 0.5_RKIND) / denominator

                   ! convert back
                   stress11(iVertexOnCell,iCell) = 0.5_RKIND * (stress1 + stress2)
                   stress22(iVertexOnCell,iCell) = 0.5_RKIND * (stress1 - stress2)

#else
                   !$omp simd
                do iVertexOnCell = 1, nEdgesOnCell(iCell)

                   call seaice_evp_constitutive_relation(&
                        stress11(iVertexOnCell,iCell), &
                        stress22(iVertexOnCell,iCell), &
                        stress12(iVertexOnCell,iCell), &
                        strain11(iVertexOnCell,iCell), &
                        strain22(iVertexOnCell,iCell), &
                        strain12(iVertexOnCell,iCell), &
                        icePressure(iCell), &
                        replacementPressure(iVertexOnCell,iCell), &
                        areaCell(iCell), &
                        elasticTimeStep)
#endif
                enddo ! iVertexOnCell

             endif ! solveStress

          enddo ! iCell

       else if (constitutiveRelationType == REVISED_EVP_CONSTITUTIVE_RELATION) then

#ifdef MPAS_OPENMP
!$omp parallel do default(shared) private(iVertexOnCell)
#endif
          do iCell = 1, nCells

             if (solveStress(iCell) == 1) then

                !$omp simd
                do iVertexOnCell = 1, nEdgesOnCell(iCell)

                   call seaice_evp_constitutive_relation_revised(&
                        stress11(iVertexOnCell,iCell), &
                        stress22(iVertexOnCell,iCell), &
                        stress12(iVertexOnCell,iCell), &
                        strain11(iVertexOnCell,iCell), &
                        strain22(iVertexOnCell,iCell), &
                        strain12(iVertexOnCell,iCell), &
                        icePressure(iCell), &
                        replacementPressure(iVertexOnCell,iCell), &
                        areaCell(iCell))

                enddo ! iVertexOnCell

             endif ! solveStress

          enddo ! iCell

       else if (constitutiveRelationType == LINEAR_CONSTITUTIVE_RELATION) then

#ifdef MPAS_OPENMP
!$omp parallel do default(shared) private(iCell, iVertexOnCell)
#endif
          do iCell = 1, nCells

             if (solveStress(iCell) == 1) then

                !$omp simd
                do iVertexOnCell = 1, nEdgesOnCell(iCell)

                   call seaice_linear_constitutive_relation(&
                        stress11(iVertexOnCell,iCell), &
                        stress22(iVertexOnCell,iCell), &
                        stress12(iVertexOnCell,iCell), &
                        strain11(iVertexOnCell,iCell), &
                        strain22(iVertexOnCell,iCell), &
                        strain12(iVertexOnCell,iCell))

                enddo ! iVertexOnCell

             endif ! solveStress

          enddo ! iCell

       endif ! constitutiveRelationType

       blockPtr => blockPtr % next
    end do

  end subroutine seaice_stress_tensor_variational!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_stress_tensor_variational_linear
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 3rd October 2019
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_stress_tensor_variational_linear(domain)!{{{

    use seaice_velocity_solver_constitutive_relation, only: &
         seaice_linear_constitutive_relation

    use seaice_mesh_pool, only: &
         nCells, &
         nEdgesOnCell, &
         stress11, &
         stress22, &
         stress12, &
         solveStress

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         velocityVariationalPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12

    integer :: &
         iCell, &
         iVertexOnCell

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

       call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11)
       call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22)
       call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12)

#ifdef MPAS_OPENMP
!$omp parallel do default(shared) private(iCell, iVertexOnCell)
#endif
       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             !$omp simd
             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                call seaice_linear_constitutive_relation(&
                     stress11(iVertexOnCell,iCell), &
                     stress22(iVertexOnCell,iCell), &
                     stress12(iVertexOnCell,iCell), &
                     strain11(iVertexOnCell,iCell), &
                     strain22(iVertexOnCell,iCell), &
                     strain12(iVertexOnCell,iCell))

             enddo ! iVertexOnCell

          endif ! solveStress

       enddo ! iCell

       blockPtr => blockPtr % next
    end do

  end subroutine seaice_stress_tensor_variational_linear!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_stress_divergence_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_stress_divergence_variational(domain)!{{{

    use seaice_mesh_pool, only: &
         nVerticesSolve, &
         cellsOnVertex, &
         nEdgesOnCell, &
         areaTriangle, &
         vertexDegree, &
         basisIntegralsMetric, &
         basisIntegralsU, &
         basisIntegralsV, &
         cellVerticesAtVertex, &
         solveVelocity, &
         stress11, &
         stress22, &
         stress12, &
         tanLatVertexRotatedOverRadius

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         velocityVariationalPool

    real(kind=RKIND), dimension(:), pointer :: &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), dimension(:), pointer :: &
         variationalDenominator

    real(kind=RKIND) :: &
         stressDivergenceUVertex, &
         stressDivergenceVVertex, &
         stressDivergenceUCell, &
         stressDivergenceVCell

    integer :: &
         iVertex, &
         iSurroundingCell, &
         iCell, &
         iStressVertex, &
         iVelocityVertex

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)

       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)

       call MPAS_pool_get_array(velocityVariationalPool, "variationalDenominator", variationalDenominator)

       ! loop over velocity positions
#ifdef MPAS_OPENMP_OFFLOAD
!$omp target teams distribute parallel do
#elif MPAS_OPENACC
!$acc parallel loop gang worker
#else
!$omp parallel do default(shared) private(stressDivergenceUVertex, stressDivergenceVVertex, &
!$omp&   iSurroundingCell, iCell, iVelocityVertex, stressDivergenceUCell, stressDivergenceVCell, iStressVertex)
#endif
       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             stressDivergenceUVertex = 0.0_RKIND
             stressDivergenceVVertex = 0.0_RKIND

             ! loop over surrounding cells
             do iSurroundingCell = 1, vertexDegree

                ! get the cell number of this cell
                iCell = cellsOnVertex(iSurroundingCell, iVertex)

                ! get the vertexOnCell number of the iVertex velocity point from cell iCell
                iVelocityVertex = cellVerticesAtVertex(iSurroundingCell,iVertex)

                stressDivergenceUCell = 0.0_RKIND
                stressDivergenceVCell = 0.0_RKIND

                ! loop over the vertices of the surrounding cell
                do iStressVertex = 1, nEdgesOnCell(iCell)

                   ! normal & metric terms
                   stressDivergenceUCell = stressDivergenceUCell - &
                        stress11(iStressVertex,iCell) * basisIntegralsU(iStressVertex,iVelocityVertex,iCell) - &
                        stress12(iStressVertex,iCell) * basisIntegralsV(iStressVertex,iVelocityVertex,iCell) - &
                        stress12(iStressVertex,iCell) * basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) * &
                        tanLatVertexRotatedOverRadius(iVertex)

                   stressDivergenceVCell = stressDivergenceVCell - &
                        stress22(iStressVertex,iCell) * basisIntegralsV(iStressVertex,iVelocityVertex,iCell) - &
                        stress12(iStressVertex,iCell) * basisIntegralsU(iStressVertex,iVelocityVertex,iCell) + &
                        stress11(iStressVertex,iCell) * basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) * &
                        tanLatVertexRotatedOverRadius(iVertex)

                enddo ! iStressVertex

                stressDivergenceUVertex = stressDivergenceUVertex + stressDivergenceUCell
                stressDivergenceVVertex = stressDivergenceVVertex + stressDivergenceVCell

             enddo ! iSurroundingCell

             stressDivergenceU(iVertex) = stressDivergenceUVertex / variationalDenominator(iVertex)
             stressDivergenceV(iVertex) = stressDivergenceVVertex / variationalDenominator(iVertex)

          endif ! solveVelocity

       enddo ! iVertex

       blockPtr => blockPtr % next
    end do

  end subroutine seaice_stress_divergence_variational!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  final_divergence_shear_variational
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date July 9th 2015
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_final_divergence_shear_variational(domain)

    use seaice_velocity_solver_constitutive_relation, only: &
         eccentricitySquared

    use seaice_mesh_pool, only: &
         nCells, &
         solveStress, &
         nEdgesOnCell

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         blockPtr

    type(MPAS_pool_type), pointer :: &
         velocityVariationalPool, &
         velocitySolverPool, &
         ridgingPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), dimension(:), pointer :: &
         divergence, &
         shear, &
         ridgeConvergence, &
         ridgeShear

    real(kind=RKIND), dimension(:), allocatable :: &
         DeltaAverage

    real(kind=RKIND) :: &
         strainDivergenceSum, &
         strainTensionSum, &
         strainShearingSum, &
         strainDivergence, &
         strainTension, &
         strainShearing, &
         Delta

    logical, pointer :: &
         config_use_column_physics

    integer :: &
         iCell, &
         iVertexOnCell

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_variational", velocityVariationalPool)
       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11)
       call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22)
       call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12)

       call MPAS_pool_get_array(velocitySolverPool, "divergence", divergence)
       call MPAS_pool_get_array(velocitySolverPool, "shear", shear)

       allocate(DeltaAverage(nCells))

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             strainDivergenceSum = 0.0_RKIND
             strainTensionSum    = 0.0_RKIND
             strainShearingSum   = 0.0_RKIND
             DeltaAverage(iCell) = 0.0_RKIND

             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                strainDivergence = strain11(iVertexOnCell,iCell) + strain22(iVertexOnCell,iCell)
                strainTension    = strain11(iVertexOnCell,iCell) - strain22(iVertexOnCell,iCell)
                strainShearing   = strain12(iVertexOnCell,iCell) * 2.0_RKIND

                Delta = sqrt(strainDivergence**2 + (strainTension**2 + strainShearing**2) / eccentricitySquared)

                strainDivergenceSum = strainDivergenceSum + strainDivergence
                strainTensionSum    = strainTensionSum    + strainTension
                strainShearingSum   = strainShearingSum   + strainShearing
                DeltaAverage(iCell) = DeltaAverage(iCell) + Delta

             enddo ! iVertexOnCell

             divergence(iCell)   = strainDivergenceSum                              / real(nEdgesOnCell(iCell),RKIND)
             shear(iCell)        = sqrt(strainTensionSum**2 + strainShearingSum**2) / real(nEdgesOnCell(iCell),RKIND)
             DeltaAverage(iCell) = DeltaAverage(iCell)                              / real(nEdgesOnCell(iCell),RKIND)

          else

             divergence(iCell)   = 0.0_RKIND
             shear(iCell)        = 0.0_RKIND

          endif

       enddo ! iCell

       ! ridging parameters
       call MPAS_pool_get_config(blockPtr % configs, "config_use_column_physics", config_use_column_physics)

       if (config_use_column_physics) then

          call MPAS_pool_get_subpool(blockPtr % structs, "ridging", ridgingPool)

          call MPAS_pool_get_array(ridgingPool, "ridgeConvergence", ridgeConvergence)
          call MPAS_pool_get_array(ridgingPool, "ridgeShear", ridgeShear)

          do iCell = 1, nCells

             if (solveStress(iCell) == 1) then

                ridgeConvergence(iCell) = -min(divergence(iCell),0.0_RKIND)
                ridgeShear(iCell)       = 0.5_RKIND * (DeltaAverage(iCell) - abs(divergence(iCell)))

             else

                ridgeConvergence(iCell) = 0.0_RKIND
                ridgeShear(iCell)       = 0.0_RKIND

             endif

          enddo ! iCell

       endif ! config_use_column_physics

       ! units - for comparison to CICE
       divergence = divergence * 100.0_RKIND * 86400.0_RKIND
       shear      = shear      * 100.0_RKIND * 86400.0_RKIND

       ! cleanup
       deallocate(DeltaAverage)

       blockPtr => blockPtr % next
    enddo

  end subroutine seaice_final_divergence_shear_variational

!-----------------------------------------------------------------------

end module seaice_velocity_solver_variational
