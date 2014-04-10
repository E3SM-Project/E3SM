module cice_dynamics_shared

  use mpas_grid_types

  implicit none

  private
  public :: &
       init_evp, &
       ice_strength, &
       evp_constitutive_relation, &
       linear_constitutive_relation, &
       air_stress, &
       ocean_stress, &
       ocean_stress_coefficient, &
       solve_velocity, &
       surface_tilt, &
       principal_stresses

  real(kind=RKIND), parameter, public :: &
       sinOceanTurningAngle = 0.0_RKIND, & ! northern hemisphere
       cosOceanTurningAngle = 1.0_RKIND, & ! northern hemisphere
       areaMinimum = 0.001_RKIND, &
       massMinimum = 0.01_RKIND
  
  real(kind=RKIND), parameter :: &
       eccentricity = 2.0_RKIND, &
       eccentricitySquared = eccentricity**2, &
       dampingTimescaleParameter = 0.36_RKIND, &
       puny = 1.0e-11_RKIND

  real(kind=RKIND), public :: &
       dtDynamics, &
       dtElastic, &
       dampingTimescale, &
       elapsedTime, &
       evpDampingCriterion

  integer, parameter, public :: &
       nElasticSubcycle = 120 ! will be from namelist

contains

  !-------------------------------------------------------------

  subroutine init_evp(dt)

    real(kind=RKIND), intent(in) :: dt

    dtDynamics = dt

    dtElastic = dtDynamics / real(nElasticSubcycle,RKIND)

    dampingTimescale = dampingTimescaleParameter * dtDynamics

    evpDampingCriterion = (1230.0_RKIND * dampingTimescale) / dtElastic**2 

    elapsedTime = 0.0_RKIND

  end subroutine init_evp

  !-------------------------------------------------------------

  subroutine ice_strength(mesh,        icePressure,   &
                          iceAreaCell, iceVolumeCell, &
                          solveStress)

    type(mesh_type), intent(in) :: mesh

    real(kind=RKIND), dimension(:), intent(out) :: &
         icePressure

    real(kind=RKIND), dimension(:), intent(in) :: &
         iceAreaCell, &
         iceVolumeCell

    integer, dimension(:), intent(in) :: &
         solveStress

    integer :: &
         iCell

    real(kind=RKIND), parameter :: &
         Pstar = 2.75e4_RKIND, & ! constant in Hibler strength formula 
         Cstar = 20.0_RKIND      ! constant in Hibler strength formula 

    do iCell = 1, mesh % nCells

       if (solveStress(iCell) == 1) then

          icePressure(iCell) = Pstar * iceVolumeCell(iCell) * exp(-Cstar*(1.0_RKIND-iceAreaCell(iCell)))

       else

          icePressure(iCell) = 0.0_RKIND

       endif ! solveStress

    enddo ! iCell

  end subroutine ice_strength

  !-------------------------------------------------------------

  subroutine evp_constitutive_relation(stress11,    stress22,    stress12, &
                                       strain11,    strain22,    strain12, &
                                       icePressure, replacementPressure,   &
                                       areaCell,    dtElastic,   evpDamping)

    real(kind=RKIND), intent(inout) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), intent(in) :: &
         strain11, &
         strain22, &
         strain12

    real(kind=RKIND), intent(in) :: &
         icePressure, & 
         dtElastic,   &
         areaCell

    real(kind=RKIND), intent(out) :: &
         replacementPressure

    logical, intent(in) :: &
         evpDamping

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
    Delta = sqrt(strainDivergence**2 + &
         (strainTension**2 + strainShearing**2) / eccentricitySquared)

    if (evpDamping) then
       pressureCoefficient = icePressure  / max(Delta,4.0_RKIND*puny)
       replacementPressure = pressureCoefficient * Delta
       pressureCoefficient = min(pressureCoefficient,evpDampingCriterion*areaCell)
    else
       pressureCoefficient = icePressure / max(Delta,puny)
       replacementPressure = pressureCoefficient * Delta
    endif

    pressureCoefficient = (pressureCoefficient * dtElastic) / (2.0_RKIND * dampingTimescale)

    denominator = 1.0_RKIND + (0.5_RKIND * dtElastic) / dampingTimescale

    stress1  = (stress1  +  pressureCoefficient                        * (strainDivergence - Delta))  / denominator

    stress2  = (stress2  + (pressureCoefficient / eccentricitySquared) *  strainTension             ) / denominator

    stress12 = (stress12 + (pressureCoefficient / eccentricitysquared) *  strainShearing * 0.5_RKIND) / denominator

    ! convert back
    stress11 = 0.5_RKIND * (stress1 + stress2)
    stress22 = 0.5_RKIND * (stress1 - stress2)

  end subroutine evp_constitutive_relation

  !-------------------------------------------------------------

  subroutine solve_velocity(nPoints,           solvePoints,       &
                            uVelocity,         vVelocity,         &
                            mass,              coriolisParameter, &
                            stressDivergenceU, stressDivergenceV, &
                            airStressU,        airStressV,        &
                            surfaceTiltForceU, surfaceTiltForceV, &
                            oceanStressU,      oceanStressV,      &
                            oceanStressCoeff,  dtElastic,mesh)

    integer, intent(in) :: &
         nPoints

    integer, dimension(:), intent(in) :: &
         solvePoints

    real(kind=RKIND), dimension(:), intent(out) :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:), intent(in) :: &
         mass,              &
         coriolisParameter, &
         stressDivergenceU, &
         stressDivergenceV, &
         airStressU,        &
         airStressV,        &
         surfaceTiltForceU, &
         surfaceTiltForceV, &
         oceanStressU,      &
         oceanStressV,      &
         oceanStressCoeff

    real(kind=RKIND), intent(in) :: &
         dtElastic
    type(mesh_type), intent(in) :: mesh
    real(kind=RKIND), dimension(2) :: &
         rightHandSide

    real(kind=RKIND), dimension(2,2) :: &
         leftMatrix

    real(kind=RKIND) :: &
         solutionDenominator

    integer :: &
         iPoint,i,j

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          ! U equation
          leftMatrix(1,1) =  mass(iPoint) / dtElastic                 + oceanStressCoeff(iPoint) * cosOceanTurningAngle
          leftMatrix(1,2) = -mass(iPoint) * coriolisParameter(iPoint) - oceanStressCoeff(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))

          ! V equation
          leftMatrix(2,1) =  mass(iPoint) * coriolisParameter(iPoint) + oceanStressCoeff(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))
          leftMatrix(2,2) =  mass(iPoint) / dtElastic                 + oceanStressCoeff(iPoint) * cosOceanTurningAngle

          ! right hand side of matrix solve
          rightHandSide(1) = stressDivergenceU(iPoint) + airStressU(iPoint) + surfaceTiltForceU(iPoint) + &
                             oceanStressCoeff(iPoint) * oceanStressU(iPoint) + (mass(iPoint) * uVelocity(iPoint)) / dtElastic

          rightHandSide(2) = stressDivergenceV(iPoint) + airStressV(iPoint) + surfaceTiltForceV(iPoint) + &
                             oceanStressCoeff(iPoint) * oceanStressV(iPoint) + (mass(iPoint) * vVelocity(iPoint)) / dtElastic

          ! solve the equation
          solutionDenominator = leftMatrix(1,1) * leftMatrix(2,2) - leftMatrix(1,2) * leftMatrix(2,1)

          uVelocity(iPoint) = (leftMatrix(2,2) * rightHandSide(1) - leftMatrix(1,2) * rightHandSide(2)) / solutionDenominator
          vVelocity(iPoint) = (leftMatrix(1,1) * rightHandSide(2) - leftMatrix(2,1) * rightHandSide(1)) / solutionDenominator

          write(70,*) mesh % POPindxiv % array(iPoint), mesh % POPindxjv % array(iPoint), &
               uVelocity(iPoint), coriolisParameter(iPoint)

       else

          uVelocity(iPoint) = 0.0_RKIND
          vVelocity(iPoint) = 0.0_RKIND

       endif ! solvePoints

    enddo ! iCell

  end subroutine solve_velocity

  !-------------------------------------------------------------

  subroutine air_stress(airStressU, airStressV, uAirVelocity, vAirVelocity, airDensity, iceArea)

    real(kind=RKIND), dimension(:), intent(out) :: &
         airStressU, &
         airStressV

    real(kind=RKIND), dimension(:), intent(in) :: &
         uAirVelocity, &
         vAirVelocity, &
         airDensity, &
         iceArea

    real(kind=RKIND) :: &
         windSpeed

    integer :: &
         nPoints, &
         iPoint

    real(kind=RKIND), parameter :: &
         airStressCoeff = 0.0012_RKIND 
    
    nPoints = size(airStressU,1)

    do iPoint = 1, nPoints

       windSpeed = sqrt(uAirVelocity(iPoint)**2 + vAirVelocity(iPoint)**2)

       airStressU(iPoint) = airDensity(iPoint) * windSpeed * airStressCoeff * uAirVelocity(iPoint) * iceArea(iPoint)
       airStressV(iPoint) = airDensity(iPoint) * windSpeed * airStressCoeff * vAirVelocity(iPoint) * iceArea(iPoint)

    enddo ! iPoint

  end subroutine air_stress

  !-------------------------------------------------------------

  subroutine ocean_stress(oceanStressU,      oceanStressV,   &
                          uOceanVelocity,    vOceanVelocity, &
                          coriolisParameter, solvePoints)

    real(kind=RKIND), dimension(:), intent(out) :: &
         oceanStressU, &
         oceanStressV

    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocity, &
         vOceanVelocity, &
         coriolisParameter

    integer, dimension(:), intent(in) :: &
         solvePoints

    integer :: &
         nPoints, &
         iPoint

    nPoints = size(oceanStressU,1)

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          oceanStressU(iPoint) = uOceanVelocity(iPoint) * cosOceanTurningAngle - &
                                 vOceanVelocity(iPoint) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))
          oceanStressV(iPoint) = uOceanVelocity(iPoint) * sinOceanTurningAngle + &
                                 vOceanVelocity(iPoint) * cosOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))
          
       else

          oceanStressU(iPoint) = 0.0_RKIND
          oceanStressV(iPoint) = 0.0_RKIND

       endif ! solvePoints

    enddo ! iPoint

  end subroutine ocean_stress

  !-------------------------------------------------------------

  subroutine ocean_stress_coefficient(oceanStressCoeff, &
                                      uOceanVelocity, vOceanVelocity, &
                                      uVelocity, vVelocity, &
                                      iceArea, solvePoints) 

    real(kind=RKIND), dimension(:), intent(out) :: &
         oceanStressCoeff

    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocity, &
         vOceanVelocity, &
         uVelocity, &
         vVelocity, &
         iceArea

    integer, dimension(:), intent(in) :: &
         solvePoints

    real(kind=RKIND), parameter :: &
         densityWater      = 1026.0_RKIND, &
         iceOceanDragCoeff = 0.00536_RKIND

    integer :: &
         nPoints, &
         iPoint

    nPoints = size(oceanStressCoeff,1)

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          oceanStressCoeff(iPoint) = iceOceanDragCoeff * densityWater * iceArea(iPoint) * &
                                     sqrt((uOceanVelocity(iPoint) - uVelocity(iPoint))**2 + &
                                          (vOceanVelocity(iPoint) - vVelocity(iPoint))**2)
       else

          oceanStressCoeff(iPoint) = 0.0_RKIND

       endif

    enddo ! iPoint

  end subroutine ocean_stress_coefficient
  
  !-------------------------------------------------------------

  subroutine surface_tilt(surfaceTiltForceU, surfaceTiltForceV, &
                          uOceanVelocity,    vOceanVelocity,    &
                          totalMass,         coriolisParameter, &
                          solvePoints)

    real(kind=RKIND), dimension(:), intent(out) :: &
         surfaceTiltForceU, &
         surfaceTiltForceV

    real(kind=RKIND), dimension(:), intent(in) :: &
         uOceanVelocity,  &
         vOceanVelocity,  &
         totalMass,       &
         coriolisParameter

    integer, dimension(:), intent(in) :: &
         solvePoints

    integer :: &
         nPoints, &
         iPoint

    nPoints = size(surfaceTiltForceU,1)

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          surfaceTiltForceU(iPoint) = -coriolisParameter(iPoint) * totalMass(iPoint) * vOceanVelocity(iPoint)
          surfaceTiltForceV(iPoint) =  coriolisParameter(iPoint) * totalMass(iPoint) * uOceanVelocity(iPoint)

       else

          surfaceTiltForceU(iPoint) = 0.0_RKIND
          surfaceTiltForceV(iPoint) = 0.0_RKIND

       endif ! solvePoints

    enddo ! iPoint

  end subroutine surface_tilt

  !-------------------------------------------------------------

  subroutine principal_stresses(principalStress1, principalStress2, &
                                stress11, stress22, stress12, &
                                replacementPressure)

    real(kind=RKIND), dimension(:), intent(out) :: &
         principalStress1, &
         principalStress2

    real(kind=RKIND), dimension(:), intent(in) :: &
         stress11, &
         stress22, &
         stress12, &
         replacementPressure

    real(kind=RKIND) :: &
         sqrtContents

    integer :: &
         iPoint, &
         nPoints

    nPoints = size(principalStress1)

    open(48,file='principal.txt')

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

       write(48,*) principalStress1(iPoint), principalStress2(iPoint)

    enddo ! iPoint

    close(48)

  end subroutine principal_stresses

  !-------------------------------------------------------------
  ! Testing
  !-------------------------------------------------------------

  subroutine linear_constitutive_relation(stress11, stress22, stress12, &
                                          strain11, strain22, strain12)

    use cice_testing, only: lambda

    real(kind=RKIND), intent(out) :: &
         stress11, &
         stress22, &
         stress12

    real(kind=RKIND), intent(in) :: &
         strain11, &
         strain22, &
         strain12

    stress11 = lambda * strain11
    stress22 = lambda * strain22
    stress12 = lambda * strain12

  end subroutine linear_constitutive_relation

  !-------------------------------------------------------------

end module cice_dynamics_shared
