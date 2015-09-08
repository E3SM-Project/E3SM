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

  subroutine cice_ice_strength(&
       nCells, &
       nCategories, &
       icePressure, &
       iceAreaCell, &
       iceVolumeCell, &
       openWaterArea, &
       iceAreaCategory, &
       iceVolumeCategory, &
       solveStress, &
       config_use_column_package, &
       config_use_column_vertical_thermodynamics)

    use ice_colpkg, only: &
         colpkg_ice_strength

    integer, intent(in) :: &
         nCells, &  !< Input:
         nCategories !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         icePressure !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         iceAreaCell, &   !< Input:
         iceVolumeCell, & !< Input:
         openWaterArea    !< Input:

    real(kind=RKIND), dimension(:,:,:), intent(in) :: &
         iceAreaCategory, & !< Input:
         iceVolumeCategory  !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress !< Input:

    logical, intent(in) :: &
         config_use_column_package, &              !< Input:
         config_use_column_vertical_thermodynamics !< Input:

    integer :: &
         iCell

    real(kind=RKIND), parameter :: &
         Pstar = 2.75e4_RKIND, & ! constant in Hibler strength formula 
         Cstar = 20.0_RKIND      ! constant in Hibler strength formula 

    if (.not. config_use_column_package .or. &
       (config_use_column_package .and. .not. config_use_column_vertical_thermodynamics)) then

       do iCell = 1, nCells

          if (solveStress(iCell) == 1) then

             icePressure(iCell) = Pstar * iceVolumeCell(iCell) * exp(-Cstar*(1.0_RKIND-iceAreaCell(iCell)))

          else

             icePressure(iCell) = 0.0_RKIND

          endif ! solveStress

       enddo ! iCell

    else

       do iCell = 1, nCells

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

  subroutine cice_solve_velocity(&
       nPoints, &
       solvePoints, &
       uVelocity, &
       vVelocity, &
       uVelocityInitial, &
       vVelocityInitial, &
       mass, &
       massCoriolisParameter, &
       stressDivergenceU, &
       stressDivergenceV, &
       airStressU, &
       airStressV, &
       surfaceTiltForceU, &
       surfaceTiltForceV, &
       oceanStressU, &
       oceanStressV, &
       oceanStressCoeff, &
       dtElastic, &
       mesh)

    integer, intent(in) :: &
         nPoints !< Input:

    integer, dimension(:), intent(in) :: &
         solvePoints !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, & !< Output: 
         vVelocity    !< Output:

    real(kind=RKIND), dimension(:), intent(inout) :: &
         uVelocityInitial, &      !< Input:
         vVelocityInitial, &      !< Input:
         mass, &                  !< Input:
         massCoriolisParameter, & !< Input:
         stressDivergenceU, &     !< Input:
         stressDivergenceV, &     !< Input:
         airStressU, &            !< Input:
         airStressV, &            !< Input:
         surfaceTiltForceU, &     !< Input:
         surfaceTiltForceV, &     !< Input:
         oceanStressU, &          !< Input:
         oceanStressV, &          !< Input: 
         oceanStressCoeff         !< Input:

    real(kind=RKIND), intent(in) :: &
         dtElastic !< Input:

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator

    integer :: &
         iPoint

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then
          
          ! U equation
          leftMatrix(1,1) =  mass(iPoint) / dtElastic      + oceanStressCoeff(iPoint) * cosOceanTurningAngle
          leftMatrix(1,2) = -massCoriolisParameter(iPoint) - oceanStressCoeff(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,massCoriolisParameter(iPoint))

          ! V equation
          leftMatrix(2,1) =  massCoriolisParameter(iPoint) + oceanStressCoeff(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,massCoriolisParameter(iPoint))
          leftMatrix(2,2) =  mass(iPoint) / dtElastic      + oceanStressCoeff(iPoint) * cosOceanTurningAngle

          ! right hand side of matrix solve
          rightHandSide(1) = stressDivergenceU(iPoint) + airStressU(iPoint) + surfaceTiltForceU(iPoint) + &
                             oceanStressCoeff(iPoint) * oceanStressU(iPoint) + (mass(iPoint) * uVelocity(iPoint)) / dtElastic

          rightHandSide(2) = stressDivergenceV(iPoint) + airStressV(iPoint) + surfaceTiltForceV(iPoint) + &
                             oceanStressCoeff(iPoint) * oceanStressV(iPoint) + (mass(iPoint) * vVelocity(iPoint)) / dtElastic

          ! solve the equation
          solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

          uVelocity(iPoint) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
          vVelocity(iPoint) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

       else

          uVelocity(iPoint) = 0.0_RKIND
          vVelocity(iPoint) = 0.0_RKIND

       endif ! solvePoints       

    enddo ! iPoint

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

  subroutine cice_solve_velocity_revised(&
       nPoints, &
       solvePoints, &
       uVelocity, &
       vVelocity, &
       uVelocityInitial, &
       vVelocityInitial, &
       mass, &
       massCoriolisParameter, &
       stressDivergenceU, &
       stressDivergenceV, &
       airStressU, &
       airStressV, &
       surfaceTiltForceU, &
       surfaceTiltForceV, &
       oceanStressU, &
       oceanStressV, &
       oceanStressCoeff, &
       dt, &
       mesh)

    integer, intent(in) :: &
         nPoints !< Input:

    integer, dimension(:), intent(in) :: &
         solvePoints !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, & !< Output: 
         vVelocity    !< Output:

    real(kind=RKIND), dimension(:), intent(inout) :: &
         uVelocityInitial, &      !< Input:
         vVelocityInitial, &      !< Input:
         mass, &                  !< Input:
         massCoriolisParameter, & !< Input:
         stressDivergenceU, &     !< Input:
         stressDivergenceV, &     !< Input:
         airStressU, &            !< Input:
         airStressV, &            !< Input:
         surfaceTiltForceU, &     !< Input:
         surfaceTiltForceV, &     !< Input:
         oceanStressU, &          !< Input:
         oceanStressV, &          !< Input: 
         oceanStressCoeff         !< Input:

    real(kind=RKIND), intent(in) :: &
         dt !< Input:

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh !< Input:

    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator

    integer :: &
         iPoint
 
    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then
          
          ! U equation
          leftMatrix(1,1) =  (numericalInertiaCoefficient + 1.0_RKIND) * (mass(iPoint) / dt) &
                          +  oceanStressCoeff(iPoint) * cosOceanTurningAngle
          leftMatrix(1,2) = -massCoriolisParameter(iPoint) &
                          -  oceanStressCoeff(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,massCoriolisParameter(iPoint))

          ! V equation
          leftMatrix(2,1) =  massCoriolisParameter(iPoint) &
                          +  oceanStressCoeff(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,massCoriolisParameter(iPoint))
          leftMatrix(2,2) =  (numericalInertiaCoefficient + 1.0_RKIND) * (mass(iPoint) / dt) &
                          +  oceanStressCoeff(iPoint) * cosOceanTurningAngle

          ! right hand side of matrix solve
          rightHandSide(1) = stressDivergenceU(iPoint) + airStressU(iPoint) + surfaceTiltForceU(iPoint) + &
                             oceanStressCoeff(iPoint) * oceanStressU(iPoint) + &
                             (mass(iPoint) * (numericalInertiaCoefficient * uVelocity(iPoint) + uVelocityInitial(iPoint))) / dt

          rightHandSide(2) = stressDivergenceV(iPoint) + airStressV(iPoint) + surfaceTiltForceV(iPoint) + &
                             oceanStressCoeff(iPoint) * oceanStressV(iPoint) + &
                             (mass(iPoint) * (numericalInertiaCoefficient * vVelocity(iPoint) + vVelocityInitial(iPoint))) / dt

          ! solve the equation
          solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

          uVelocity(iPoint) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
          vVelocity(iPoint) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator
          
       else

          uVelocity(iPoint) = 0.0_RKIND
          vVelocity(iPoint) = 0.0_RKIND

       endif ! solvePoints       

    enddo ! iPoint

  end subroutine cice_solve_velocity_revised
  
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_air_stress
!
!> \brief 
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>  
!
!-----------------------------------------------------------------------

  subroutine cice_air_stress(&
       nPoints, &
       airStressU, &
       airStressV, &
       uAirVelocity, &
       vAirVelocity, &
       airDensity, &
       iceArea)

    integer, intent(in) :: &
         nPoints !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         airStressU, & !< Output: 
         airStressV    !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         uAirVelocity, & !< Input:
         vAirVelocity, & !< Input:
         airDensity, &   !< Input:
         iceArea         !< Input:

    real(kind=RKIND) :: &
         windSpeed

    integer :: &
         iPoint

    real(kind=RKIND), parameter :: &
         airStressCoeff = 0.0012_RKIND 

    do iPoint = 1, nPoints

       windSpeed = sqrt(uAirVelocity(iPoint)**2 + vAirVelocity(iPoint)**2)

       airStressU(iPoint) = airDensity(iPoint) * windSpeed * airStressCoeff * uAirVelocity(iPoint) * iceArea(iPoint)
       airStressV(iPoint) = airDensity(iPoint) * windSpeed * airStressCoeff * vAirVelocity(iPoint) * iceArea(iPoint)

    enddo ! iPoint

  end subroutine cice_air_stress

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

  subroutine cice_ocean_stress(&
       nPoints, &
       oceanStressU, &
       oceanStressV, &
       uOceanVelocity, &
       vOceanVelocity, &
       coriolisParameter, &
       solvePoints)

    integer, intent(in) :: &
         nPoints

    real(kind=RKIND), dimension(:), intent(out) :: &
         oceanStressU, & !< Output: 
         oceanStressV    !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocity, & !< Input:
         vOceanVelocity, & !< Input:
         coriolisParameter !< Input:

    integer, dimension(:), intent(in) :: &
         solvePoints !< Input:

    integer :: &
         iPoint

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          oceanStressU(iPoint) = uOceanVelocity(iPoint) * cosOceanTurningAngle - &
                                 vOceanVelocity(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))
          oceanStressV(iPoint) = uOceanVelocity(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint)) + &
                                 vOceanVelocity(iPoint) * cosOceanTurningAngle

          !oceanStressU(iPoint) = 0.0_RKIND
          !oceanStressV(iPoint) = 0.0_RKIND
          
       else

          oceanStressU(iPoint) = 0.0_RKIND
          oceanStressV(iPoint) = 0.0_RKIND

       endif ! solvePoints

    enddo ! iPoint

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

  subroutine cice_ocean_stress_coefficient(&
       nPoints, &
       oceanStressCoeff, &
       uOceanVelocity, &
       vOceanVelocity, &
       uVelocity, &
       vVelocity, &
       iceArea, &
       solvePoints) 

    integer, intent(in) :: &
         nPoints !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         oceanStressCoeff !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocity, & !< Input:
         vOceanVelocity, & !< Input:
         uVelocity, &      !< Input:
         vVelocity, &      !< Input:
         iceArea           !< Input:

    integer, dimension(:), intent(in) :: &
         solvePoints !< Input:

    real(kind=RKIND), parameter :: &
         densityWater      = 1026.0_RKIND, &
         iceOceanDragCoeff = 0.00536_RKIND

    integer :: &
         iPoint

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          oceanStressCoeff(iPoint) = iceOceanDragCoeff * densityWater * iceArea(iPoint) * &
                                     sqrt((uOceanVelocity(iPoint) - uVelocity(iPoint))**2 + &
                                          (vOceanVelocity(iPoint) - vVelocity(iPoint))**2)

          !oceanStressCoeff(iPoint) = 0.0_RKIND

       else

          oceanStressCoeff(iPoint) = 0.0_RKIND

       endif

    enddo ! iPoint

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

  subroutine cice_ocean_stress_final(&
       mesh, &
       boundary, &
       nVerticesSolve, &
       iceAreaVertex, &
       coriolisParameter, &
       uVelocity, &
       vVelocity, &
       uOceanVelocityVertex, &
       vOceanVelocityVertex, &
       oceanStressU, &
       oceanStressV, &
       oceanStressCoeff, &
       oceanStressCellU, &
       oceanStressCellV, &
       solveVelocity)

    use cice_mesh, only: &
         cice_interpolate_vertex_to_cell

    type(MPAS_pool_type), pointer, intent(in) :: &
         mesh, &
         boundary

    integer, intent(in) :: &
         nVerticesSolve

    real(kind=RKIND), dimension(:), intent(in) :: &    
         iceAreaVertex, &        !< Input:
         coriolisParameter, &    !< Input:
         uVelocity, &            !< Input:
         vVelocity, &            !< Input:
         uOceanVelocityVertex, & !< Input:
         vOceanVelocityVertex    !< Input:
    
    real(kind=RKIND), dimension(:), intent(out) :: &
         oceanStressU, &     !< Output: 
         oceanStressV, &     !< Output:
         oceanStressCellU, & !< Output:
         oceanStressCellV, & !< Output: 
         oceanStressCoeff    !< Output: 

    integer, dimension(:), intent(in) :: &
         solveVelocity !< Input:

    real(kind=RKIND) :: &
         oceanIceRelativeVelocity

    real(kind=RKIND), dimension(:), allocatable :: &
         oceanStressUVertexTmp, &
         oceanStressVVertexTmp

    integer :: &
         iVertex

    call cice_ocean_stress_coefficient(&
         nVerticesSolve, &
         oceanStressCoeff, &
         uOceanVelocityVertex, &
         vOceanVelocityVertex, &
         uVelocity, &
         vVelocity, &
         iceAreaVertex, &
         solveVelocity)

    allocate(oceanStressUVertexTmp(nVerticesSolve))
    allocate(oceanStressVVertexTmp(nVerticesSolve))

    do iVertex = 1, nVerticesSolve

       if (solveVelocity(iVertex) == 1) then

          oceanStressU(iVertex) = oceanStressCoeff(iVertex) * ((uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * cosOceanTurningAngle &
                                                             - (vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iVertex)))
          oceanStressV(iVertex) = oceanStressCoeff(iVertex) * ((vOceanVelocityVertex(iVertex) - vVelocity(iVertex)) * cosOceanTurningAngle &
                                                             + (uOceanVelocityVertex(iVertex) - uVelocity(iVertex)) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iVertex)))

          oceanStressUVertexTmp(iVertex) = oceanStressU(iVertex) / iceAreaVertex(iVertex)
          oceanStressVVertexTmp(iVertex) = oceanStressV(iVertex) / iceAreaVertex(iVertex)

       else

          oceanStressU(iVertex) = 0.0_RKIND
          oceanStressV(iVertex) = 0.0_RKIND

       endif ! solveVelocity

    enddo ! iVertex

    ! interpolate stress to cells
    call cice_interpolate_vertex_to_cell(mesh, boundary, oceanStressCellU, oceanStressUVertexTmp)
    call cice_interpolate_vertex_to_cell(mesh, boundary, oceanStressCellV, oceanStressVVertexTmp)

    deallocate(oceanStressUVertexTmp)
    deallocate(oceanStressVVertexTmp)

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

  subroutine cice_surface_tilt(&
       nPoints, &
       surfaceTiltForceU, &
       surfaceTiltForceV, &
       uOceanVelocity, &
       vOceanVelocity, &
       totalMass, &
       coriolisParameter, &
       solvePoints)

    integer, intent(in) :: &
         nPoints

    real(kind=RKIND), dimension(:), intent(out) :: &
         surfaceTiltForceU, & !< Output: 
         surfaceTiltForceV    !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocity, & !< Input:
         vOceanVelocity, & !< Input:
         totalMass, &      !< Input:
         coriolisParameter !< Input:

    integer, dimension(:), intent(in) :: &
         solvePoints !< Input:

    integer :: &
         iPoint

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          surfaceTiltForceU(iPoint) = -coriolisParameter(iPoint) * totalMass(iPoint) * vOceanVelocity(iPoint)
          surfaceTiltForceV(iPoint) =  coriolisParameter(iPoint) * totalMass(iPoint) * uOceanVelocity(iPoint)

          !surfaceTiltForceU(iPoint) = 0.0_RKIND
          !surfaceTiltForceV(iPoint) = 0.0_RKIND

       else

          surfaceTiltForceU(iPoint) = 0.0_RKIND
          surfaceTiltForceV(iPoint) = 0.0_RKIND

       endif ! solvePoints

    enddo ! iPoint

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

  subroutine cice_coriolis_force_coefficient(&
       nPoints, &
       totalMassCoriolisParameter, &
       totalMass, &
       coriolisParameter)

    integer, intent(in) :: &
         nPoints

    real(kind=RKIND), dimension(:), intent(out) :: &
         totalMassCoriolisParameter !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         totalMass, &      !< Input:
         coriolisParameter !< Input:

    integer :: &
         iPoint

    do iPoint = 1, nPoints

       totalMassCoriolisParameter(iPoint) = totalMass(iPoint) * coriolisParameter(iPoint)

    enddo ! iPoint

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
