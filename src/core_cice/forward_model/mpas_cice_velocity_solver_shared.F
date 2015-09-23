!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_velocity_solver_shared
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

module cice_velocity_solver_shared

  use mpas_derived_types
  use mpas_pool_routines

  implicit none

  private
  save

  public :: &
       cice_init_evp, &
       cice_ice_strength, &
       cice_evp_constitutive_relation, &
       cice_evp_constitutive_relation_revised, &       
       cice_air_stress, &
       cice_ocean_stress, &
       cice_ocean_stress_coefficient, &
       cice_ocean_stress_final, &
       cice_solve_velocity, &
       cice_solve_velocity_revised, &       
       cice_surface_tilt, &
       cice_coriolis_force_coefficient, &
       cice_principal_stresses

  real(kind=RKIND), parameter, public :: &
       ciceAreaMinimum = 0.001_RKIND, &
       ciceMassMinimum = 0.01_RKIND
  
  ! general EVP parameters
  real(kind=RKIND), parameter, private :: &
       sinOceanTurningAngle = 0.0_RKIND, & ! northern hemisphere
       cosOceanTurningAngle = 1.0_RKIND, & ! northern hemisphere
       eccentricity = 2.0_RKIND, &
       dampingTimescaleParameter = 0.36_RKIND, &
       puny = 1.0e-11_RKIND
  
  real(kind=RKIND), parameter, public :: &
       eccentricitySquared = eccentricity**2

  real(kind=RKIND), private :: &
       dampingTimescale, &
       evpDampingCriterion

  ! Bouillon et al. 2013 parameters
  real(kind=RKIND), parameter, private :: &
       dampingRatioDenominator = 0.86_RKIND, & ! Se > 0.5
       dampingRatio = 5.5e-3_RKIND  ! xi = Sv/Sc < 1

  real(kind=RKIND), private :: &
       numericalInertiaCoefficient ! brlx

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_evp
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_init_evp(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: block

    integer, pointer :: &
         config_elastic_subcycle_number

    real(kind=RKIND), pointer :: &
         config_dt

    type(MPAS_pool_type), pointer :: &
         velocitySolver, &
         mesh

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep, &
         elasticTimeStep

    real(kind=RKIND), dimension(:), pointer :: &
         dvEdge

    integer, pointer :: &
         nEdgesSolve
        
    real(kind=RKIND) :: &
         gamma, &
         dvEdgeMin, &
         dvEdgeMinGlobal

    call MPAS_pool_get_config(domain % configs, "config_elastic_subcycle_number", config_elastic_subcycle_number)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)

    ! general EVP
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolver)
       call MPAS_pool_get_array(velocitySolver, "dynamicsTimeStep", dynamicsTimeStep)
       call MPAS_pool_get_array(velocitySolver, "elasticTimeStep", elasticTimeStep)

       dynamicsTimeStep = config_dt
       
       elasticTimeStep = dynamicsTimeStep / real(config_elastic_subcycle_number,RKIND)
       
       dampingTimescale = dampingTimescaleParameter * dynamicsTimeStep
       
       evpDampingCriterion = (1230.0_RKIND * dampingTimescale) / elasticTimeStep**2 

       block => block % next
    enddo

    ! find the minimum edge length in the grid
    dvEdgeMin = 1.0e30_RKIND
    
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
       
       call MPAS_pool_get_dimension(mesh, "nEdgesSolve", nEdgesSolve)
       
       call MPAS_pool_get_array(mesh, "dvEdge", dvEdge)

       dvEdgeMin = min(dvEdgeMin, minval(dvEdge(1:nEdgesSolve)))

       block => block % next
    enddo

    call mpas_dmpar_min_real(domain % dminfo, dvEdgeMin, dvEdgeMinGlobal)

    !!!! Testing!
    !dvEdgeMinGlobal = 8558.2317072059941_RKIND
    
    ! Bouillon et al. 2013
    block => domain % blocklist
    do while (associated(block))
       
       gamma = 0.25_RKIND * 1.0e11_RKIND * config_dt
       numericalInertiaCoefficient = (2.0_RKIND * dampingRatioDenominator * dampingRatio * gamma) / dvEdgeMinGlobal**2
       
       block => block % next
    enddo

  end subroutine cice_init_evp

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_ice_strength
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_ice_strength(domain)

    use ice_colpkg, only: &
         colpkg_ice_strength

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         tracersAggregatePool, &
         icestatePool, &
         tracersPool

    real(kind=RKIND), dimension(:), pointer:: &
         icePressure, &
         iceAreaCell, &
         iceVolumeCell, &
         openWaterArea

    real(kind=RKIND), dimension(:,:,:), pointer :: &
         iceAreaCategory, &
         iceVolumeCategory

    integer, dimension(:), pointer :: &
         solveStress

    type(field1DReal), pointer :: &
         icePressureField

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics

    integer, pointer :: &
         nCellsSolve, &
         nCategories

    integer :: &
         iCell

    real(kind=RKIND), parameter :: &
         Pstar = 2.75e4_RKIND, & ! constant in Hibler strength formula 
         Cstar = 20.0_RKIND      ! constant in Hibler strength formula 

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_column_package", config_use_column_package)
       call MPAS_pool_get_config(block % configs, "config_use_column_vertical_thermodynamics", config_use_column_vertical_thermodynamics)

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nCategories", nCategories)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)

       call MPAS_pool_get_array(velocitySolverPool, "icePressure", icePressure)
       call MPAS_pool_get_array(velocitySolverPool, "solveStress", solveStress)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)

       call MPAS_pool_get_array(icestatePool, "openWaterArea", openWaterArea)

       call MPAS_pool_get_array(tracersPool, "iceAreaCategory", iceAreaCategory, 1)
       call MPAS_pool_get_array(tracersPool, "iceVolumeCategory", iceVolumeCategory, 1)

       if (.not. config_use_column_package .or. &
            (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then

          do iCell = 1, nCellsSolve

             if (solveStress(iCell) == 1) then

                icePressure(iCell) = Pstar * iceVolumeCell(iCell) * exp(-Cstar*(1.0_RKIND-iceAreaCell(iCell)))

             else

                icePressure(iCell) = 0.0_RKIND

             endif ! solveStress

          enddo ! iCell

       else

          do iCell = 1, nCellsSolve

             icePressure(iCell) = 0.0_RKIND

             if (solveStress(iCell) == 1) then

                ! this routine doesnt reset icePressure
                call colpkg_ice_strength(&
                     nCategories, &
                     iceAreaCell(iCell), &
                     iceVolumeCell(iCell), &
                     openWaterArea(iCell), &
                     iceAreaCategory(1,:,iCell), &
                     iceVolumeCategory(1,:,iCell), &
                     icePressure(iCell))

             endif ! solveStress

          enddo ! iCell

       endif

       block => block % next
    enddo

    ! halo exchange ice strength
    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)
    call MPAS_pool_get_field(velocitySolverPool, "icePressure", icePressureField)
    call MPAS_dmpar_exch_halo_field(icePressureField)

  end subroutine cice_ice_strength

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_evp_constitutive_relation
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_evp_constitutive_relation(&
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12, &
       icePressure, &
       replacementPressure, &
       areaCell, &
       dtElastic)

    real(kind=RKIND), intent(inout) :: &
         stress11, & !< Input/Output: 
         stress22, & !< Input/Output: 
         stress12    !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         strain11, & !< Input:
         strain22, & !< Input:
         strain12    !< Input:

    real(kind=RKIND), intent(in) :: &
         icePressure, & !< Input:
         dtElastic, &   !< Input:
         areaCell       !< Input:

    real(kind=RKIND), intent(out) :: &
         replacementPressure !< Output: 

    real(kind=RKIND) :: &
         strainDivergence,    &
         strainTension,       &
         strainShearing,      &
         stress1,             &
         stress2,             &
         Delta,               &
         pressureCoefficient, &
         denominator
    
    ! convert from stress11 to stress1 etc
    strainDivergence = strain11 + strain22
    strainTension    = strain11 - strain22
    strainShearing   = strain12 * 2.0_RKIND
    
    stress1 = stress11 + stress22
    stress2 = stress11 - stress22

    ! perform the constituitive relation
    Delta = sqrt(strainDivergence**2 + (strainTension**2 + strainShearing**2) / eccentricitySquared)
    
    pressureCoefficient = icePressure / max(Delta,puny)
    replacementPressure = pressureCoefficient * Delta

    pressureCoefficient = (pressureCoefficient * dtElastic) / (2.0_RKIND * dampingTimescale)

    denominator = 1.0_RKIND + (0.5_RKIND * dtElastic) / dampingTimescale

    stress1  = (stress1  +  pressureCoefficient                        * (strainDivergence - Delta))  / denominator
    stress2  = (stress2  + (pressureCoefficient / eccentricitySquared) *  strainTension             ) / denominator
    stress12 = (stress12 + (pressureCoefficient / eccentricitysquared) *  strainShearing * 0.5_RKIND) / denominator
    
    ! convert back
    stress11 = 0.5_RKIND * (stress1 + stress2)
    stress22 = 0.5_RKIND * (stress1 - stress2)

  end subroutine cice_evp_constitutive_relation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_evp_constitutive_relation_revised
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_evp_constitutive_relation_revised(&
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12, &
       icePressure, &
       replacementPressure, &
       areaCell, &
       dtElastic)

    real(kind=RKIND), intent(inout) :: &
         stress11, & !< Input/Output: 
         stress22, & !< Input/Output: 
         stress12    !< Input/Output: 

    real(kind=RKIND), intent(in) :: &
         strain11, & !< Input:
         strain22, & !< Input:
         strain12    !< Input:

    real(kind=RKIND), intent(in) :: &
         icePressure, & !< Input:
         dtElastic, &   !< Input:
         areaCell       !< Input:

    real(kind=RKIND), intent(out) :: &
         replacementPressure !< Output:
    
    real(kind=RKIND) :: &
         strainDivergence,    &
         strainTension,       &
         strainShearing,      &
         stress1,             &
         stress2,             &
         Delta,               &
         pressureCoefficient, &
         denominator
    
    ! convert from stress11 to stress1 etc
    strainDivergence = strain11 + strain22
    strainTension    = strain11 - strain22
    strainShearing   = strain12 * 2.0_RKIND
    
    stress1 = stress11 + stress22
    stress2 = stress11 - stress22

    ! perform the constituitive relation
    Delta = sqrt(strainDivergence**2 + (strainTension**2 + strainShearing**2) / eccentricitySquared)
    
    pressureCoefficient = icePressure / max(Delta,puny)
    replacementPressure = pressureCoefficient * Delta

    pressureCoefficient = (pressureCoefficient * 2.0_RKIND * dampingRatio) / dampingRatioDenominator

    denominator = 1.0_RKIND + (2.0_RKIND * dampingRatio) / dampingRatioDenominator

    stress1  = (stress1  +  pressureCoefficient                        * (strainDivergence - Delta))  / denominator
    stress2  = (stress2  + (pressureCoefficient / eccentricitySquared) *  strainTension             ) / denominator
    stress12 = (stress12 + (pressureCoefficient / eccentricitysquared) *  strainShearing * 0.5_RKIND) / denominator
    
    ! convert back
    stress11 = 0.5_RKIND * (stress1 + stress2)
    stress22 = 0.5_RKIND * (stress1 - stress2)

  end subroutine cice_evp_constitutive_relation_revised

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_linear_constitutive_relation
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_linear_constitutive_relation(&
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12)

    real(kind=RKIND), intent(out) :: &
         stress11, & !< Output: 
         stress22, & !< Output: 
         stress12    !< Output: 

    real(kind=RKIND), intent(in) :: &
         strain11, & !< Input:
         strain22, & !< Input:
         strain12    !< Input:

    real(kind=RKIND), parameter :: &
         lambda = 1.0_RKIND

    stress11 = lambda * strain11
    stress22 = lambda * strain22
    stress12 = lambda * strain12

  end subroutine cice_linear_constitutive_relation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_solve_velocity_Hunke13
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_solve_velocity(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool

    integer, dimension(:), pointer :: &
         solveVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         uVelocityInitial, &
         vVelocityInitial, &
         totalMassVertex, &
         totalMassVertexfVertex, &
         stressDivergenceU, &
         stressDivergenceV, &
         airStressVertexU, &
         airStressVertexV, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff

    real(kind=RKIND), pointer :: &
         elasticTimeStep

    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(velocitySolverPool, "elasticTimeStep", elasticTimeStep)
       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityInitial", uVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityInitial", vVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)

       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             ! U equation
             leftMatrix(1,1) =  totalMassVertex(iVertex) / elasticTimeStep + oceanStressCoeff(iVertex) * cosOceanTurningAngle
             leftMatrix(1,2) = -totalMassVertexfVertex(iVertex) - &
                  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))

             ! V equation
             leftMatrix(2,1) =  totalMassVertexfVertex(iVertex) + &
                  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))
             leftMatrix(2,2) =  totalMassVertex(iVertex) / elasticTimeStep  + oceanStressCoeff(iVertex) * cosOceanTurningAngle

             ! right hand side of matrix solve
             rightHandSide(1) = stressDivergenceU(iVertex) + airStressVertexU(iVertex) + surfaceTiltForceU(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressU(iVertex) + (totalMassVertex(iVertex) * uVelocity(iVertex)) / elasticTimeStep

             rightHandSide(2) = stressDivergenceV(iVertex) + airStressVertexV(iVertex) + surfaceTiltForceV(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressV(iVertex) + (totalMassVertex(iVertex) * vVelocity(iVertex)) / elasticTimeStep

             ! solve the equation
             solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

             uVelocity(iVertex) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
             vVelocity(iVertex) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

          else

             uVelocity(iVertex) = 0.0_RKIND
             vVelocity(iVertex) = 0.0_RKIND

          endif ! solveVelocity

       enddo ! iVertex

       block => block % next
    end do

  end subroutine cice_solve_velocity
  
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_solve_velocity_revised
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_solve_velocity_revised(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool

    integer, dimension(:), pointer :: &
         solveVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity, &
         uVelocityInitial, &
         vVelocityInitial, &
         totalMassVertex, &
         totalMassVertexfVertex, &
         stressDivergenceU, &
         stressDivergenceV, &
         airStressVertexU, &
         airStressVertexV, &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCoeff

    real(kind=RKIND), pointer :: &
         dt

    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_dt", dt)
 
       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

       call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "uVelocityInitial", uVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocityInitial", vVelocityInitial)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceU", stressDivergenceU)
       call MPAS_pool_get_array(velocitySolverPool, "stressDivergenceV", stressDivergenceV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)

       do iVertex = 1, nVerticesSolve

          if (solveVelocity(iVertex) == 1) then

             ! U equation
             leftMatrix(1,1) =  (numericalInertiaCoefficient + 1.0_RKIND) * (totalMassVertex(iVertex) / dt) &
                             +  oceanStressCoeff(iVertex) * cosOceanTurningAngle
             leftMatrix(1,2) = -totalMassVertexfVertex(iVertex) &
                             -  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))

             ! V equation
             leftMatrix(2,1) =  totalMassVertexfVertex(iVertex) &
                             +  oceanStressCoeff(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,totalMassVertexfVertex(iVertex))
             leftMatrix(2,2) =  (numericalInertiaCoefficient + 1.0_RKIND) * (totalMassVertex(iVertex) / dt) &
                             +  oceanStressCoeff(iVertex) * cosOceanTurningAngle

             ! right hand side of matrix solve
             rightHandSide(1) = stressDivergenceU(iVertex) + airStressVertexU(iVertex) + surfaceTiltForceU(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressU(iVertex) + &
                  (totalMassVertex(iVertex) * (numericalInertiaCoefficient * uVelocity(iVertex) + uVelocityInitial(iVertex))) / dt

             rightHandSide(2) = stressDivergenceV(iVertex) + airStressVertexV(iVertex) + surfaceTiltForceV(iVertex) + &
                  oceanStressCoeff(iVertex) * oceanStressV(iVertex) + &
                  (totalMassVertex(iVertex) * (numericalInertiaCoefficient * vVelocity(iVertex) + vVelocityInitial(iVertex))) / dt

             ! solve the equation
             solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

             uVelocity(iVertex) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
             vVelocity(iVertex) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

          else

             uVelocity(iVertex) = 0.0_RKIND
             vVelocity(iVertex) = 0.0_RKIND

          endif ! solveVelocity

       enddo ! iVertex

       block => block % next
    end do

  end subroutine cice_solve_velocity_revised

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_air_stress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 19th September 2015
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_air_stress(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    logical, pointer :: &
         config_use_column_package, &
         config_use_column_vertical_thermodynamics, &
         config_use_air_stress

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         airStressVertexU, &
         airStressVertexV, &
         airStressCellU, &
         airStressCellV

    type(field1DReal), pointer :: &
         airStressCellUField, &
         airStressCellVField!, &

    ! calculate the air stress
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_package", config_use_column_package)
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_column_vertical_thermodynamics", config_use_column_vertical_thermodynamics)

    if (.not. config_use_column_package .or. &
         (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then
       call constant_air_stress(domain)
    endif

    ! check for no air stress
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_air_stress", config_use_air_stress)
    if (.not. config_use_air_stress) then
       block => domain % blocklist
       do while (associated(block))
          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
          call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
          call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)
          airStressCellU = 0.0_RKIND
          airStressCellV = 0.0_RKIND
          block => block % next
       end do
    endif ! .not. config_use_air_stress

    ! halo exchange air stress
    call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)
    call MPAS_pool_get_field(velocitySolverPool, "airStressCellU", airStressCellUField)
    call MPAS_pool_get_field(velocitySolverPool, "airStressCellV", airStressCellVField)
    call MPAS_dmpar_exch_halo_field(airStressCellUField)
    call MPAS_dmpar_exch_halo_field(airStressCellVField)

    ! interpolate air stress
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexU", airStressVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressVertexV", airStressVertexV)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            airStressVertexU, &
            airStressCellU)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            airStressVertexV, &
            airStressCellV)

       block => block % next
    end do

  end subroutine cice_air_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  constant_air_stress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine constant_air_stress(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         atmosCouplingPool, &
         tracersAggregatePool

    real(kind=RKIND), dimension(:), pointer :: &
         airStressCellU, &
         airStressCellV, &
         uAirVelocity, &
         vAirVelocity, &
         airDensity, &
         iceAreaCell

    real(kind=RKIND) :: &
         windSpeed

    integer, pointer :: &
         nCellsSolve

    integer :: &
         iCell

    real(kind=RKIND), parameter :: &
         airStressCoeff = 0.0012_RKIND 

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

       call MPAS_pool_get_array(velocitySolverPool, "airStressCellU", airStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "airStressCellV", airStressCellV)

       call MPAS_pool_get_array(atmosCouplingPool, "uAirVelocity", uAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "vAirVelocity", vAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "airDensity", airDensity)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)

       do iCell = 1, nCellsSolve

          windSpeed = sqrt(uAirVelocity(iCell)**2 + vAirVelocity(iCell)**2)

          airStressCellU(iCell) = airDensity(iCell) * windSpeed * airStressCoeff * uAirVelocity(iCell) * iceAreaCell(iCell)
          airStressCellV(iCell) = airDensity(iCell) * windSpeed * airStressCoeff * vAirVelocity(iCell) * iceAreaCell(iCell)

       enddo ! iCell

       block => block % next
    end do

  end subroutine constant_air_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_ocean_stress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_ocean_stress(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         oceanStressU, &
         oceanStressV, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         fVertex

    integer, dimension(:), pointer :: &
         solveVelocity

    logical, pointer :: &
         configUseOceanStress

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_ocean_stress", configUseOceanStress)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)

       if (configUseOceanStress) then

          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)

          call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                oceanStressU(iVertex) = uOceanVelocityVertex(iVertex) * cosOceanTurningAngle - &
                                        vOceanVelocityVertex(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex))
                oceanStressV(iVertex) = uOceanVelocityVertex(iVertex) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex)) + &
                                        vOceanVelocityVertex(iVertex) * cosOceanTurningAngle

             else

                oceanStressU(iVertex) = 0.0_RKIND
                oceanStressV(iVertex) = 0.0_RKIND

             endif ! solvePoints

          enddo ! iVertex

       else

          ! no ocean stress
          oceanStressU = 0.0_RKIND
          oceanStressV = 0.0_RKIND

       endif

       block => block % next
    end do

  end subroutine cice_ocean_stress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_ocean_stress_coefficient
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_ocean_stress_coefficient(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         icestatePool

    real(kind=RKIND), dimension(:), pointer :: &
         oceanStressCoeff, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         uVelocity, &
         vVelocity, &
         iceAreaVertex

    integer, dimension(:), pointer :: &
         solveVelocity

    logical, pointer :: &
         configUseOceanStress

    real(kind=RKIND), parameter :: &
         densityWater      = 1026.0_RKIND, &
         iceOceanDragCoeff = 0.00536_RKIND

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_ocean_stress", configUseOceanStress)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)

       if (configUseOceanStress) then

          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

          call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)

          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                oceanStressCoeff(iVertex) = iceOceanDragCoeff * densityWater * iceAreaVertex(iVertex) * &
                     sqrt((uOceanVelocityVertex(iVertex) - uVelocity(iVertex))**2 + &
                          (vOceanVelocityVertex(iVertex) - vVelocity(iVertex))**2)

             else

                oceanStressCoeff(iVertex) = 0.0_RKIND

             endif

          enddo ! iVertex

       else

          ! no ocean stress
          oceanStressCoeff = 0.0_RKIND

       endif

       block => block % next
    end do

  end subroutine cice_ocean_stress_coefficient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_ocean_stress_final
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_ocean_stress_final(domain)

    use cice_mesh, only: &
         cice_interpolate_vertex_to_cell

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         boundaryPool, &
         icestatePool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaVertex, &
         fVertex, &
         uVelocity, &
         vVelocity, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         oceanStressU, &
         oceanStressV, &
         oceanStressCellU, &
         oceanStressCellV, &
         oceanStressCoeff

    integer, dimension(:), pointer :: &
         solveVelocity

    type(field1DReal), pointer :: &
         oceanStressUField, &
         oceanStressVField

    logical, pointer :: &
         configUseOceanStress

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    ! get ocean stress coefficient
    call cice_ocean_stress_coefficient(domain)

    ! use stress config
    call MPAS_pool_get_config(domain % blocklist % configs, "config_use_ocean_stress", configUseOceanStress)

    ! get ocean stress on vertices
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
       call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)

       if (configUseOceanStress) then

          call MPAS_pool_get_dimension(meshPool, "nVerticesSolve", nVerticesSolve)

          call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

          call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)

          call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCoeff", oceanStressCoeff)
          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)

          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                oceanStressU(iVertex) = oceanStressCoeff(iVertex) * ((uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * cosOceanTurningAngle &
                                                                   - (vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex)))
                oceanStressV(iVertex) = oceanStressCoeff(iVertex) * ((vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * cosOceanTurningAngle &
                                                                   + (uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * sinOceanTurningAngle * sign(1.0_RKIND,fVertex(iVertex)))

                oceanStressU(iVertex) = oceanStressU(iVertex) / iceAreaVertex(iVertex)
                oceanStressV(iVertex) = oceanStressV(iVertex) / iceAreaVertex(iVertex)

             else

                oceanStressU(iVertex) = 0.0_RKIND
                oceanStressV(iVertex) = 0.0_RKIND

             endif ! solveVelocity

          enddo ! iVertex

       else

          ! no ocean stress
          oceanStressU = 0.0_RKIND
          oceanStressV = 0.0_RKIND

       endif

       block => block % next
    end do

    ! get ocean stress on cells
    if (configUseOceanStress) then

       ! halo exchange ocean stress
       call MPAS_pool_get_subpool(domain % blocklist % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_field(velocitySolverPool, "oceanStressU", oceanStressUField)
       call MPAS_pool_get_field(velocitySolverPool, "oceanStressV", oceanStressVField)

       call MPAS_dmpar_exch_halo_field(oceanStressUField)
       call MPAS_dmpar_exch_halo_field(oceanStressVField)

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

          call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
          call MPAS_pool_get_subpool(block % structs, "boundary", boundaryPool)
          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
          call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)

          call MPAS_pool_get_array(icestatePool, "iceAreaVertex", iceAreaVertex)

          call MPAS_pool_get_array(velocitySolverPool, "oceanStressU", oceanStressU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressV", oceanStressV)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)

          ! interpolate ocean stress to cells
          call cice_interpolate_vertex_to_cell(meshPool, boundaryPool, oceanStressCellU, oceanStressU)
          call cice_interpolate_vertex_to_cell(meshPool, boundaryPool, oceanStressCellV, oceanStressV)

          ! multiply ocean stress back by ice area
          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                oceanStressU(iVertex) = oceanStressU(iVertex) * iceAreaVertex(iVertex)
                oceanStressV(iVertex) = oceanStressV(iVertex) * iceAreaVertex(iVertex)

             endif ! solveVelocity

          enddo ! iVertex

          block => block % next
       end do

    else

       block => domain % blocklist
       do while (associated(block))

          call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellU", oceanStressCellU)
          call MPAS_pool_get_array(velocitySolverPool, "oceanStressCellV", oceanStressCellV)

          oceanStressCellU = 0.0_RKIND
          oceanStressCellV = 0.0_RKIND

          block => block % next
       end do

    endif ! configUseOceanStress

  end subroutine cice_ocean_stress_final
  
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_surface_tilt
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_surface_tilt(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         velocitySolverPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer:: &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         totalMassVertex, &
         fVertex

    integer, dimension(:), pointer :: &
         solveVelocity

    logical, pointer :: &
         configUseSurfaceTilt

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_use_surface_tilt", configUseSurfaceTilt)

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceU", surfaceTiltForceU)
       call MPAS_pool_get_array(velocitySolverPool, "surfaceTiltForceV", surfaceTiltForceV)

       if (configUseSurfaceTilt) then

          call MPAS_pool_get_array(velocitySolverPool, "solveVelocity", solveVelocity)
          call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
          call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)

          call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)

          call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

          do iVertex = 1, nVerticesSolve

             if (solveVelocity(iVertex) == 1) then

                surfaceTiltForceU(iVertex) = -fVertex(iVertex) * totalMassVertex(iVertex) * vOceanVelocityVertex(iVertex)
                surfaceTiltForceV(iVertex) =  fVertex(iVertex) * totalMassVertex(iVertex) * uOceanVelocityVertex(iVertex)

             else

                surfaceTiltForceU(iVertex) = 0.0_RKIND
                surfaceTiltForceV(iVertex) = 0.0_RKIND

             endif ! solveVelocity

          enddo ! iVertex

       else

          ! no surface tilt
          surfaceTiltForceU = 0.0_RKIND
          surfaceTiltForceV = 0.0_RKIND

       endif

       block => block % next
    end do

  end subroutine cice_surface_tilt

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_coriolis_force_coefficient
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_coriolis_force_coefficient(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         icestatePool, &
         velocitySolverPool, &
         meshPool

    real(kind=RKIND), dimension(:), pointer :: &
         totalMassVertexfVertex, &
         totalMassVertex, &
         fVertex

    integer, pointer :: &
         nVerticesSolve

    integer :: &
         iVertex

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)

       call MPAS_pool_get_subpool(block % structs, "icestate", icestatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)

       call MPAS_pool_get_array(icestatePool, "totalMassVertex", totalMassVertex)
       call MPAS_pool_get_array(velocitySolverPool, "totalMassVertexfVertex", totalMassVertexfVertex)
       call MPAS_pool_get_array(meshPool, "fVertex", fVertex)

       do iVertex = 1, nVerticesSolve

          totalMassVertexfVertex(iVertex) = totalMassVertex(iVertex) * fVertex(iVertex)

       enddo ! iVertex

       block => block % next
    end do

  end subroutine cice_coriolis_force_coefficient

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_principal_stresses
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_principal_stresses(&
       nPoints, &
       principalStress1, &
       principalStress2, &
       stress11, &
       stress22, &
       stress12, &
       replacementPressure)

    integer, intent(in) :: &
         nPoints

    real(kind=RKIND), dimension(:), intent(out) :: &
         principalStress1, & !< Output: 
         principalStress2    !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &         !< Input:
         stress22, &         !< Input:
         stress12, &         !< Input:
         replacementPressure !< Input:

    real(kind=RKIND) :: &
         sqrtContents

    integer :: &
         iPoint

    do iPoint = 1, nPoints

       if (replacementPressure(iPoint) > puny) then

          sqrtContents = (stress11(iPoint) + stress22(iPoint))**2 - &
                         4.0_RKIND * stress11(iPoint) * stress22(iPoint) + &
                         4.0_RKIND * stress12(iPoint)**2

          principalStress1(iPoint) = 0.5_RKIND * (stress11(iPoint) + stress22(iPoint)) + 0.5_RKIND * sqrt(sqrtContents)
          principalStress2(iPoint) = 0.5_RKIND * (stress11(iPoint) + stress22(iPoint)) - 0.5_RKIND * sqrt(sqrtContents)

          principalStress1(iPoint) = principalStress1(iPoint) / replacementPressure(iPoint)
          principalStress2(iPoint) = principalStress2(iPoint) / replacementPressure(iPoint)

       else
          
          principalStress1(iPoint) = 1.0e30_RKIND
          principalStress2(iPoint) = 1.0e30_RKIND

       endif

    enddo ! iPoint

  end subroutine cice_principal_stresses

!-----------------------------------------------------------------------

end module cice_velocity_solver_shared
