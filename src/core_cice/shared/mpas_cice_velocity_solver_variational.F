!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_variational
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver_variational

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timer

  implicit none

  private
  save

  public :: &
       cice_init_velocity_solver_variational, &
       cice_strain_tensor_variational, &
       cice_stress_divergence_variational, &
       cice_internal_stress_variational, &
       cice_final_divergence_shear_variational

contains

!-----------------------------------------------------------------------
! initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_velocity_solver_variational
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 24 October 2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_velocity_solver_variational(&
       mesh, &
       velocity_variational, &
       boundary, &
       rotateCartesianGrid, &
       includeMetricTerms, &
       variationalBasisType, &
       integrationType, &
       integrationOrder)!{{{

    use cice_velocity_solver_wachspress, only: &
         cice_init_velocity_solver_wachspress

    use cice_velocity_solver_pwl, only: &
         cice_init_velocity_solver_pwl

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    type(MPAS_pool_type), pointer :: &
         velocity_variational, & !< Input/Output:
         boundary                !< Input/Output:

    logical, intent(in) :: &
         rotateCartesianGrid, & !< Input:
         includeMetricTerms     !< Input:

    character(len=*), intent(in) :: &
         variationalBasisType, & !< Input:
         integrationType         !< Input:

    integer, intent(in) :: &
         integrationOrder !< Input:

    if (trim(variationalBasisType) == "wachspress") then

       call cice_init_velocity_solver_wachspress(&
            mesh, &
            velocity_variational, &
            boundary, &
            rotateCartesianGrid, &
            includeMetricTerms, &
            integrationType, &
            integrationOrder)

    else if (trim(variationalBasisType) == "pwl") then

       call cice_init_velocity_solver_pwl(&
            mesh, &
            velocity_variational, &
            boundary, &
            rotateCartesianGrid, &
            includeMetricTerms)

    endif

  end subroutine cice_init_velocity_solver_variational

!-----------------------------------------------------------------------
! time step
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_internal_stress_variational
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_internal_stress_variational(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type (MPAS_pool_type), pointer :: &
         meshPool, &
         velocityVariationalPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         icePressure, &
         replacementPressure, &
         stressDivergenceU, &
         stressDivergenceV

    real(kind=RKIND), pointer :: &
         elasticTimeStep

    logical, pointer :: &
         revisedEVP

    integer, dimension(:), pointer :: &
         solveStress, &
         solveVelocity

    integer, dimension(:,:), pointer :: &
         cellVerticesAtVertex

    real(kind=RKIND), dimension(:), pointer :: &
         tanLatVertexRotatedOverRadius

    real(kind=RKIND), dimension(:,:), pointer :: &
         strain11, &
         strain22, &
         strain12, &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         basisGradientU, &
         basisGradientV, &
         basisIntegralsU, &
         basisIntegralsV, &
         basisIntegralsMetric

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_revised_evp", revisedEVP)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "icePressure", icePressure)
       call MPAS_pool_get_array(velocitySolverPool, "replacementPressure", replacementPressure)
       call MPAS_pool_get_array(velocitySolverPool, "elasticTimeStep", elasticTimeStep)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)

       call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11)
       call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22)
       call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12)
       call MPAS_pool_get_array(velocityVariationalPool, "stress11", stress11)
       call MPAS_pool_get_array(velocityVariationalPool, "stress22", stress22)
       call MPAS_pool_get_array(velocityVariationalPool, "stress12", stress12)
       call MPAS_pool_get_array(velocityVariationalPool, "cellVerticesAtVertex", cellVerticesAtVertex)
       call MPAS_pool_get_array(velocityVariationalPool, "tanLatVertexRotatedOverRadius", tanLatVertexRotatedOverRadius)
       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientU", basisGradientU)
       call MPAS_pool_get_array(velocityVariationalPool, "basisGradientV", basisGradientV)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsU", basisIntegralsU)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsV", basisIntegralsV)
       call MPAS_pool_get_array(velocityVariationalPool, "basisIntegralsMetric", basisIntegralsMetric)

       call mpas_timer_start("Velocity solver strain tensor")
       call cice_strain_tensor_variational(&
            meshPool, &
            strain11, &
            strain22, &
            strain12, &
            uVelocity, &
            vVelocity, &
            basisGradientU, &
            basisGradientV, &
            tanLatVertexRotatedOverRadius, &
            solveStress)
       call mpas_timer_stop("Velocity solver strain tensor")

       call mpas_timer_start("Velocity solver stress tensor")
       call cice_stress_tensor_variational(&
            meshPool, &
            stress11, &
            stress22, &
            stress12, &
            strain11, &
            strain22, &
            strain12, &
            icePressure, &
            replacementPressure, &
            solveStress, &
            elasticTimeStep, &
            revisedEVP)
       call mpas_timer_stop("Velocity solver stress tensor")

       call mpas_timer_start("Velocity solver stress divergence")
       call cice_stress_divergence_variational(&
            meshPool, &
            stressDivergenceU, &
            stressDivergenceV, &
            stress11, &
            stress22, &
            stress12, &
            basisIntegralsU, &
            basisIntegralsV, &
            basisIntegralsMetric, &
            tanLatVertexRotatedOverRadius, &
            cellVerticesAtVertex, &
            solveVelocity)
       call mpas_timer_stop("Velocity solver stress divergence")

       block => block % next
    end do

  end subroutine cice_internal_stress_variational!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_strain_tensor_variational
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_strain_tensor_variational(&
       mesh, &
       strain11, &
       strain22, &
       strain12, &
       uVelocity, &
       vVelocity, &
       basisGradientU, &
       basisGradientV, &
       tanLatVertexRotatedOverRadius, &
       solveStress)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:), intent(out) :: &
         strain11, & !< Output:
         strain22, & !< Output:
         strain12    !< Output:

    real(kind=RKIND), dimension(:), intent(in) :: &
         uVelocity, & !< Input:
         vVelocity, & !< Input:
         tanLatVertexRotatedOverRadius !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         basisGradientU, & !< Input:
         basisGradientV    !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress !< Input:

    integer :: &
         iCell, &
         iGradientVertex, &
         iBasisVertex, &
         iVertex, &
         jVertex

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         verticesOnCell

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)

    ! loop over cells
    do iCell = 1, nCells

       if (solveStress(iCell) == 1) then

          strain11(:,iCell) = 0.0_RKIND
          strain22(:,iCell) = 0.0_RKIND
          strain12(:,iCell) = 0.0_RKIND

          ! loop over velocity points surrounding cell - location of stress and derivative
          do iGradientVertex = 1, nEdgesOnCell(iCell)

             ! loop over basis functions
             do iBasisVertex = 1, nEdgesOnCell(iCell)

                iVertex = verticesOnCell(iBasisVertex,iCell)

                strain11(iGradientVertex,iCell) = strain11(iGradientVertex,iCell) + &
                     uVelocity(iVertex) * basisGradientU(iBasisVertex,iGradientVertex,iCell)

                strain22(iGradientVertex,iCell) = strain22(iGradientVertex,iCell) + &
                     vVelocity(iVertex) * basisGradientV(iBasisVertex,iGradientVertex,iCell)

                strain12(iGradientVertex,iCell) = strain12(iGradientVertex,iCell) + 0.5_RKIND * (&
                     uVelocity(iVertex) * basisGradientV(iBasisVertex,iGradientVertex,iCell) + &
                     vVelocity(iVertex) * basisGradientU(iBasisVertex,iGradientVertex,iCell))

             enddo ! iVertexOnCell

             ! metric terms
             jVertex = verticesOnCell(iGradientVertex,iCell)

             strain11(iGradientVertex,iCell) = strain11(iGradientVertex,iCell) - &
                  vVelocity(jVertex) * tanLatVertexRotatedOverRadius(jVertex)

             strain12(iGradientVertex,iCell) = strain12(iGradientVertex,iCell) + &
                  uVelocity(jVertex) * tanLatVertexRotatedOverRadius(jVertex) * 0.5_RKIND

          enddo ! jVertexOnCell

       endif ! solveStress

    enddo ! iCell

  end subroutine cice_strain_tensor_variational!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_stress_tensor_variational
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_stress_tensor_variational(&
       mesh, &
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12, &
       icePressure, &
       replacementPressure, &
       solveStress, &
       dtElastic, &
       revisedEVP)!{{{

    use cice_velocity_solver_constitutive_relation, only: &
         cice_evp_constitutive_relation, &
         cice_evp_constitutive_relation_revised

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:,:), intent(inout) :: &
         stress11, & !< Input/Output:
         stress22, & !< Input/Output:
         stress12    !< Input/Output:

    real(kind=RKIND), dimension(:), intent(out) :: &
         replacementPressure !< Output:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         strain11, & !< Input:
         strain22, & !< Input:
         strain12    !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         icePressure !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress !< Input:

    real(kind=RKIND), intent(in) :: &
         dtElastic !< Input:

    logical, intent(in) :: &
         revisedEVP !< Input:

    real(kind=RKIND), dimension(:), allocatable :: &
         replacementPressureCell

    integer :: &
         iCell, &
         iVertexOnCell

    integer, pointer :: &
         nCells, &
         maxEdges

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         areaCell

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nCells", nCells)
    call MPAS_pool_get_dimension(mesh, "maxEdges", maxEdges)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "areaCell", areaCell)

    allocate(replacementPressureCell(maxEdges))

    if (.not. revisedEVP) then

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                call cice_evp_constitutive_relation(&
                     stress11(iVertexOnCell,iCell), &
                     stress22(iVertexOnCell,iCell), &
                     stress12(iVertexOnCell,iCell), &
                     strain11(iVertexOnCell,iCell), &
                     strain22(iVertexOnCell,iCell), &
                     strain12(iVertexOnCell,iCell), &
                     icePressure(iCell), &
                     replacementPressureCell(iVertexOnCell), &
                     areaCell(iCell), &
                     dtElastic)

             enddo ! iVertexOnCell

             !replacementPressure(iCell) = replacementPressureCell(3) ! for comparison to CICE quads
             replacementPressure(iCell) = sum(replacementPressureCell(1:nEdgesOnCell(iCell))) / real(nEdgesOnCell(iCell),RKIND)

          endif ! solveStress

       enddo ! iCell

    else

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             do iVertexOnCell = 1, nEdgesOnCell(iCell)

                call cice_evp_constitutive_relation_revised(&
                     stress11(iVertexOnCell,iCell), &
                     stress22(iVertexOnCell,iCell), &
                     stress12(iVertexOnCell,iCell), &
                     strain11(iVertexOnCell,iCell), &
                     strain22(iVertexOnCell,iCell), &
                     strain12(iVertexOnCell,iCell), &
                     icePressure(iCell), &
                     replacementPressureCell(iVertexOnCell), &
                     areaCell(iCell))

             enddo ! iVertexOnCell

             !replacementPressure(iCell) = replacementPressureCell(3) ! for comparison to CICE quads
             replacementPressure(iCell) = sum(replacementPressureCell(1:nEdgesOnCell(iCell))) / real(nEdgesOnCell(iCell),RKIND)

          endif ! solveStress

       enddo ! iCell

    endif

    ! clean up
    deallocate(replacementPressureCell)

  end subroutine cice_stress_tensor_variational!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_stress_divergence_variational
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_stress_divergence_variational(&
       mesh, &
       stressDivergenceU, &
       stressDivergenceV, &
       stress11, &
       stress22, &
       stress12, &
       basisIntegralsU, &
       basisIntegralsV, &
       basisIntegralsMetric, &
       tanLatVertexRotatedOverRadius, &
       cellVerticesAtVertex, &
       solveVelocity)!{{{

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         stressDivergenceU, & !< Output:
         stressDivergenceV    !< Output:

    real(kind=RKIND), dimension(:,:), intent(in) :: &
         stress11, & !< Input:
         stress22, & !< Input:
         stress12    !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         basisIntegralsU, &   !< Input:
         basisIntegralsV, &   !< Input:
         basisIntegralsMetric !< Input:

    real(kind=RKIND), dimension(:), intent(in) :: &
         tanLatVertexRotatedOverRadius !< Input:

    integer, dimension(:,:), intent(in) :: &
         cellVerticesAtVertex !< Input:

    integer, dimension(:), intent(in) :: &
         solveVelocity !< Input:

    real(kind=RKIND) :: &
         stressDivergenceUCell, &
         stressDivergenceVCell

    integer :: &
         iVertex, &
         iSurroundingCell, &
         iCell, &
         iStressVertex, &
         iVelocityVertex

    integer, pointer :: &
         nVertices, &
         vertexDegree

    integer, dimension(:), pointer :: &
         nEdgesOnCell

    integer, dimension(:,:), pointer :: &
         cellsOnVertex, &
         verticesOnCell

    real(kind=RKIND), dimension(:), pointer :: &
         areaTriangle

    ! init variables
    call MPAS_pool_get_dimension(mesh, "nVertices", nVertices)
    call MPAS_pool_get_dimension(mesh, "vertexDegree", vertexDegree)

    call MPAS_pool_get_array(mesh, "nEdgesOnCell", nEdgesOnCell)
    call MPAS_pool_get_array(mesh, "cellsOnVertex", cellsOnVertex)
    call MPAS_pool_get_array(mesh, "verticesOnCell", verticesOnCell)
    call MPAS_pool_get_array(mesh, "areaTriangle", areaTriangle)

    ! loop over velocity positions
    do iVertex = 1, nVertices

       if (solveVelocity(iVertex) == 1) then

          stressDivergenceU(iVertex) = 0.0_RKIND
          stressDivergenceV(iVertex) = 0.0_RKIND

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

                ! normal terms
                stressDivergenceUCell = stressDivergenceUCell - &
                     stress11(iStressVertex,iCell) * basisIntegralsU(iStressVertex,iVelocityVertex,iCell) - &
                     stress12(iStressVertex,iCell) * basisIntegralsV(iStressVertex,iVelocityVertex,iCell)

                stressDivergenceVCell = stressDivergenceVCell - &
                     stress22(iStressVertex,iCell) * basisIntegralsV(iStressVertex,iVelocityVertex,iCell) - &
                     stress12(iStressVertex,iCell) * basisIntegralsU(iStressVertex,iVelocityVertex,iCell)

                ! metric terms
                stressDivergenceUCell = stressDivergenceUCell - &
                     stress12(iStressVertex,iCell) * basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) * tanLatVertexRotatedOverRadius(iVertex)

                stressDivergenceVCell = stressDivergenceVCell + &
                     stress11(iStressVertex,iCell) * basisIntegralsMetric(iStressVertex,iVelocityVertex,iCell) * tanLatVertexRotatedOverRadius(iVertex)

             enddo ! iStressVertex

             stressDivergenceU(iVertex) = stressDivergenceU(iVertex) + stressDivergenceUCell
             stressDivergenceV(iVertex) = stressDivergenceV(iVertex) + stressDivergenceVCell

          enddo ! iSurroundingCell

          stressDivergenceU(iVertex) = stressDivergenceU(iVertex) / areaTriangle(iVertex)
          stressDivergenceV(iVertex) = stressDivergenceV(iVertex) / areaTriangle(iVertex)

       endif ! solveVelocity

    enddo ! iVertex

  end subroutine cice_stress_divergence_variational!}}}

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

  subroutine cice_final_divergence_shear_variational(block)

    use cice_velocity_solver_constitutive_relation, only: &
         eccentricitySquared

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocityVariationalPool, &
         velocitySolverPool, &
         ridgingPool

    integer, pointer :: &
         nCells

    integer, dimension(:), pointer :: &
         nEdgesOnCell, &
         solveStress

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
         config_use_column_package

    integer :: &
         iCell, &
         iVertexOnCell

    call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
    call MPAS_pool_get_subpool(block % structs, "velocity_variational", velocityVariationalPool)
    call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

    call MPAS_pool_get_dimension(meshPool, "nCells", nCells)

    call MPAS_pool_get_array(meshPool, "nEdgesOnCell", nEdgesOnCell)

    call MPAS_pool_get_array(velocityVariationalPool, "strain11", strain11)
    call MPAS_pool_get_array(velocityVariationalPool, "strain22", strain22)
    call MPAS_pool_get_array(velocityVariationalPool, "strain12", strain12)

    call MPAS_pool_get_array(velocitySolverPool, "divergence", divergence)
    call MPAS_pool_get_array(velocitySolverPool, "shear", shear)
    call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

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
    call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)

    if (config_use_column_package) then

       call MPAS_pool_get_subpool(block % structs, "ridging", ridgingPool)

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

    endif ! config_use_column_package

    ! units - for comparison to CICE
    divergence = divergence * 100.0_RKIND * 86400.0_RKIND
    shear      = shear      * 100.0_RKIND * 86400.0_RKIND

    ! cleanup
    deallocate(DeltaAverage)

  end subroutine cice_final_divergence_shear_variational

!-----------------------------------------------------------------------

end module cice_velocity_solver_variational
