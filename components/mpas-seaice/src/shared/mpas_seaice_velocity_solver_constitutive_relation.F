!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_velocity_solver_constitutive_relation
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_velocity_solver_constitutive_relation

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_dmpar

  implicit none

  private
  save

  public :: &
       seaice_init_evp, &
       seaice_evp_constitutive_relation, &
       seaice_evp_constitutive_relation_revised, &
       seaice_linear_constitutive_relation

  ! constitutive relation options
  integer, public :: &
       constitutiveRelationType

  integer, public :: &
       EVP_CONSTITUTIVE_RELATION = 1, &
       REVISED_EVP_CONSTITUTIVE_RELATION = 2, &
       LINEAR_CONSTITUTIVE_RELATION = 3, &
       NONE_CONSTITUTIVE_RELATION = 4

  ! general EVP parameters
  real(kind=RKIND), parameter, private :: &
       eccentricity = 2.0_RKIND, &
       dampingTimescaleParameter = 0.36_RKIND

  real(kind=RKIND), parameter, public :: &
       eccentricitySquared = eccentricity**2, &
       puny = 1.0e-11_RKIND

  real(kind=RKIND), private :: &
       evpDampingCriterion

  ! Bouillon et al. 2013 parameters
  real(kind=RKIND), parameter, private :: &
       dampingRatioDenominator = 0.86_RKIND, & ! Se > 0.5
       dampingRatio = 5.5e-3_RKIND  ! xi = Sv/Sc < 1

  real(kind=RKIND), public :: &
       dampingTimescale, &
       numericalInertiaCoefficient ! brlx

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_init_evp
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_init_evp(domain)

    type(domain_type) :: domain

    type(block_type), pointer :: blockPtr

    character(len=strKIND), pointer :: &
         config_constitutive_relation_type

    type(MPAS_pool_type), pointer :: &
         velocitySolverPool, &
         meshPool

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

    ! constitutive relation type
    call MPAS_pool_get_config(domain % configs, "config_constitutive_relation_type", config_constitutive_relation_type)
    if (trim(config_constitutive_relation_type) == 'evp') then
       constitutiveRelationType = EVP_CONSTITUTIVE_RELATION
    else if (trim(config_constitutive_relation_type) == 'evp_revised') then
       constitutiveRelationType = REVISED_EVP_CONSTITUTIVE_RELATION
    else if (trim(config_constitutive_relation_type) == 'linear') then
       constitutiveRelationType = LINEAR_CONSTITUTIVE_RELATION
    else if (trim(config_constitutive_relation_type) == 'none') then
       constitutiveRelationType = NONE_CONSTITUTIVE_RELATION
    else
       call MPAS_log_write("config_constitutive_relation_type unknown: "//trim(config_constitutive_relation_type), MPAS_LOG_CRIT)
    endif

    ! general EVP
    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_array(velocitySolverPool, "dynamicsTimeStep", dynamicsTimeStep)
       call MPAS_pool_get_array(velocitySolverPool, "elasticTimeStep", elasticTimeStep)

       dampingTimescale = dampingTimescaleParameter * dynamicsTimeStep

       evpDampingCriterion = (1230.0_RKIND * dampingTimescale) / elasticTimeStep**2

       blockPtr => blockPtr % next
    enddo

    ! find the minimum edge length in the grid
    dvEdgeMin = 1.0e30_RKIND

    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       call MPAS_pool_get_subpool(blockPtr % structs, "mesh", meshPool)

       call MPAS_pool_get_dimension(meshPool, "nEdgesSolve", nEdgesSolve)

       call MPAS_pool_get_array(meshPool, "dvEdge", dvEdge)

       dvEdgeMin = min(dvEdgeMin, minval(dvEdge(1:nEdgesSolve)))

       blockPtr => blockPtr % next
    enddo

    call mpas_dmpar_min_real(domain % dminfo, dvEdgeMin, dvEdgeMinGlobal)

    ! Bouillon et al. 2013
    blockPtr => domain % blocklist
    do while (associated(blockPtr))

       gamma = 0.25_RKIND * 1.0e11_RKIND * dynamicsTimeStep
       numericalInertiaCoefficient = (2.0_RKIND * dampingRatioDenominator * dampingRatio * gamma) / dvEdgeMinGlobal**2

       blockPtr => blockPtr % next
    enddo

  end subroutine seaice_init_evp

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_evp_constitutive_relation
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_evp_constitutive_relation(&
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
#if defined CPRINTEL && defined MPAS_OPENMP
!$omp declare simd(seaice_evp_constitutive_relation) uniform(dtElastic), notinbranch &
!$omp& linear(ref(stress11,stress22,stress12,strain11,strain22,strain12,replacementPressure,icePressure,areaCell))
#endif

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

  end subroutine seaice_evp_constitutive_relation

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_evp_constitutive_relation_revised
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_evp_constitutive_relation_revised(&
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12, &
       icePressure, &
       replacementPressure, &
       areaCell)
#if defined CPRINTEL && defined MPAS_OPENMP
!$omp declare simd(seaice_evp_constitutive_relation_revised) uniform(icePressure,areaCell), notinbranch &
!$omp& linear(ref(stress11,stress22,stress12,strain11,strain22,strain12,replacementPressure))
#endif

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

  end subroutine seaice_evp_constitutive_relation_revised

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_linear_constitutive_relation
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

  subroutine seaice_linear_constitutive_relation(&
       stress11, &
       stress22, &
       stress12, &
       strain11, &
       strain22, &
       strain12)
#if defined CPRINTEL && defined MPAS_OPENMP
!$omp declare simd(seaice_linear_constitutive_relation) notinbranch &
!$omp& linear(ref(stress11,stress22,stress12,strain11,strain22,strain12))
#endif

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

  end subroutine seaice_linear_constitutive_relation

!-----------------------------------------------------------------------

end module seaice_velocity_solver_constitutive_relation
