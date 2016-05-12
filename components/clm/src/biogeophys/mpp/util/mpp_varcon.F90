module mpp_varcon

#ifdef USE_PETSC_LIB

  !
  ! !PUBLIC TYPES:
  implicit none
  save

#include "finclude/petscsys.h"

  PetscReal :: grav        = 9.80616d0   ! gravity constant [m/s2]

  PetscReal :: cpliq       = 4.188d3     ! Specific heat of water [J/kg-K]
  PetscReal :: cpice       = 2.11727d3   ! Specific heat of ice [J/kg-K]

  PetscReal :: denh2o      = 1.000d3     ! density of liquid water [kg/m3]
  PetscReal :: denice      = 0.917d3     ! density of ice [kg/m3]

  PetscReal :: tkair       = 0.023d0     ! thermal conductivity of air   [W/m/K]
  PetscReal :: tkice       = 2.290d0     ! thermal conductivity of ice   [W/m/K]
  PetscReal :: tkwat       = 0.57d0      ! thermal conductivity of water [W/m/K]
  PetscReal :: thk_bedrock = 3.0d0       ! thermal conductivity of 'typical' saturated granitic rock 
                                         ! (Clauser and Huenges, 1995)(W/m/K)

  PetscReal :: tfrz        = 273.15d0    ! freezing temperature [K]

  PetscReal :: cnfac       = 0.5d0       ! Crank Nicholson factor between 0 and 1


  ! Initialize landunit type constants

  integer, parameter :: istsoil    = 1  !soil         landunit type (natural vegetation)
  integer, parameter :: istcrop    = 2  !crop         landunit type
  integer, parameter :: istice     = 3  !land ice     landunit type (glacier)
  integer, parameter :: istice_mec = 4  !land ice (multiple elevation classes) landunit type
  integer, parameter :: istdlak    = 5  !deep lake    landunit type (now used for all lakes)
  integer, parameter :: istwet     = 6  !wetland      landunit type (swamp, marsh, etc.)

#endif

end module mpp_varcon
