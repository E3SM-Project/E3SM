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
       cice_linear_constitutive_relation, &
       cice_air_stress, &
       cice_ocean_stress, &
       cice_ocean_stress_coefficient, &
       cice_ocean_stress_final, &
       cice_solve_velocity, &
       cice_surface_tilt, &
       cice_coriolis_force_coefficient, &
       cice_principal_stresses

  real(kind=RKIND), parameter, public :: &
       ciceAreaMinimum = 0.001_RKIND, &
       ciceMassMinimum = 0.01_RKIND
  
  ! local 
  real(kind=RKIND), parameter, private :: &
       sinOceanTurningAngle = 0.0_RKIND, & ! northern hemisphere
       cosOceanTurningAngle = 1.0_RKIND, & ! northern hemisphere
       eccentricity = 2.0_RKIND, &
       eccentricitySquared = eccentricity**2, &
       dampingTimescaleParameter = 0.36_RKIND, &
       puny = 1.0e-11_RKIND

  real(kind=RKIND), private :: &
       dampingTimescale, &
       evpDampingCriterion

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
         velocitySolver

    real(kind=RKIND), pointer :: &
         dynamicsTimeStep, &
         elasticTimeStep

    call MPAS_pool_get_config(domain % configs, "config_elastic_subcycle_number", config_elastic_subcycle_number)
    call MPAS_pool_get_config(domain % configs, "config_dt", config_dt)

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
       nPoints, &
       icePressure, &
       iceArea, &
       iceVolume, &
       solveStress)

    integer, intent(in) :: &
         nPoints !< Input:

    real(kind=RKIND), dimension(:), intent(out) :: &
         icePressure !< Output: 

    real(kind=RKIND), dimension(:), intent(in) :: &
         iceArea, & !< Input:
         iceVolume  !< Input:

    integer, dimension(:), intent(in) :: &
         solveStress !< Input:

    integer :: &
         iPoint

    real(kind=RKIND), parameter :: &
         Pstar = 2.75e4_RKIND, & ! constant in Hibler strength formula 
         Cstar = 20.0_RKIND      ! constant in Hibler strength formula 

    do iPoint = 1, nPoints

       if (solveStress(iPoint) == 1) then

          icePressure(iPoint) = Pstar * iceVolume(iPoint) * exp(-Cstar*(1.0_RKIND-iceArea(iPoint)))

       else

          icePressure(iPoint) = 0.0_RKIND

       endif ! solveStress

    enddo ! iPoint

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
       dtElastic, &
       evpDamping)

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

    logical, intent(in) :: &
         evpDamping !< Input:

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

  end subroutine cice_evp_constitutive_relation

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
!  cice_solve_velocity
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

    real(kind=RKIND), dimension(:), intent(in) :: &
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
          oceanStressV(iPoint) = uOceanVelocity(iPoint) * sinOceanTurningAngle + &
                                 vOceanVelocity(iPoint) * cosOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint))
          
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
       nPoints, &
       iceArea, &
       coriolisParameter, &
       uVelocity, &
       vVelocity, &
       uOceanVelocity, &
       vOceanVelocity, &
       oceanStressU, &
       oceanStressV, &
       oceanStressCoeff, &
       solvePoints)

    integer, intent(in) :: &
         nPoints

    real(kind=RKIND), dimension(:), intent(in) :: &    
         iceArea, &           !< Input:
         coriolisParameter, & !< Input:
         uVelocity, &         !< Input:
         vVelocity, &         !< Input:
         uOceanVelocity, &    !< Input:
         vOceanVelocity       !< Input:
    
    real(kind=RKIND), dimension(:), intent(out) :: &
         oceanStressU, &  !< Output: 
         oceanStressV, &  !< Output: 
         oceanStressCoeff !< Output: 

    integer, dimension(:), intent(in) :: &
         solvePoints !< Input:

    real(kind=RKIND) :: &
         oceanIceRelativeVelocity

    integer :: &
         iPoint

    call cice_ocean_stress_coefficient(&
         nPoints, &
         oceanStressCoeff, &
         uOceanVelocity, &
         vOceanVelocity, &
         uVelocity, &
         vVelocity, &
         iceArea, &
         solvePoints)

    do iPoint = 1, nPoints

       if (solvePoints(iPoint) == 1) then

          oceanStressU(iPoint) = oceanStressCoeff(iPoint) * ((uOceanVelocity(iPoint) - uVelocity(iPoint)) * cosOceanTurningAngle &
                                                           - (vOceanVelocity(iPoint) - vVelocity(iPoint)) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint)))
          oceanStressV(iPoint) = oceanStressCoeff(iPoint) * ((vOceanVelocity(iPoint) - vVelocity(iPoint)) * cosOceanTurningAngle &
                                                           + (uOceanVelocity(iPoint) - uVelocity(iPoint)) * sinOceanTurningAngle * sign(1.0_RKIND,coriolisParameter(iPoint)))

       else

          oceanStressU(iPoint) = 0.0_RKIND
          oceanStressV(iPoint) = 0.0_RKIND

       endif ! solvePoints

    enddo ! iPoint

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

    !open(48,file='principal.txt')

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

       !write(48,*) principalStress1(iPoint), principalStress2(iPoint)

    enddo ! iPoint

    !close(48)

  end subroutine cice_principal_stresses

!-----------------------------------------------------------------------

end module cice_velocity_solver_shared
