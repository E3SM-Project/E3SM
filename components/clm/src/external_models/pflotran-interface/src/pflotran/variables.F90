module Variables_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  ! variables definitions
  PetscInt, parameter, public :: X_COORDINATE =             1
  PetscInt, parameter, public :: Y_COORDINATE =             2
  PetscInt, parameter, public :: Z_COORDINATE =             3
  PetscInt, parameter, public :: TEMPERATURE =              4
  PetscInt, parameter, public :: LIQUID_PRESSURE =          5
  PetscInt, parameter, public :: LIQUID_SATURATION =        6
  PetscInt, parameter, public :: GAS_SATURATION =           7
  PetscInt, parameter, public :: LIQUID_DENSITY =           8
  PetscInt, parameter, public :: LIQUID_DENSITY_MOL =       9
  PetscInt, parameter, public :: GAS_DENSITY =             10
  PetscInt, parameter, public :: GAS_DENSITY_MOL =         11
  PetscInt, parameter, public :: LIQUID_ENERGY =           12
  PetscInt, parameter, public :: GAS_ENERGY =              13
  PetscInt, parameter, public :: LIQUID_VISCOSITY =        14
  PetscInt, parameter, public :: GAS_VISCOSITY =           15
  PetscInt, parameter, public :: LIQUID_MOBILITY =         16
  PetscInt, parameter, public :: GAS_MOBILITY =            17
  PetscInt, parameter, public :: LIQUID_MOLE_FRACTION =    18
  PetscInt, parameter, public :: GAS_MOLE_FRACTION =       19
  PetscInt, parameter, public :: POROSITY =                20
  PetscInt, parameter, public :: PERMEABILITY =            21
  PetscInt, parameter, public :: PERMEABILITY_X =          22
  PetscInt, parameter, public :: PERMEABILITY_Y =          23
  PetscInt, parameter, public :: PERMEABILITY_Z =          24
  PetscInt, parameter, public :: PHASE =                   25
  PetscInt, parameter, public :: MATERIAL_ID =             26

  PetscInt, parameter, public :: PRIMARY_MOLALITY =        27
  PetscInt, parameter, public :: SECONDARY_MOLALITY =      28
  PetscInt, parameter, public :: TOTAL_MOLALITY =          29
  PetscInt, parameter, public :: PRIMARY_MOLARITY =        30
  PetscInt, parameter, public :: SECONDARY_MOLARITY =      31
  PetscInt, parameter, public :: TOTAL_MOLARITY =          32
  PetscInt, parameter, public :: MINERAL_VOLUME_FRACTION = 33
  PetscInt, parameter, public :: MINERAL_RATE =            34
  PetscInt, parameter, public :: MINERAL_SURFACE_AREA =    35
  PetscInt, parameter, public :: MINERAL_SATURATION_INDEX =36
  PetscInt, parameter, public :: PH =                      37
  PetscInt, parameter, public :: IMMOBILE_SPECIES =        38
  PetscInt, parameter, public :: SURFACE_CMPLX =           39
  PetscInt, parameter, public :: SURFACE_CMPLX_FREE =      40
  PetscInt, parameter, public :: SURFACE_SITE_DENSITY =    41
  PetscInt, parameter, public :: KIN_SURFACE_CMPLX =       42
  PetscInt, parameter, public :: KIN_SURFACE_CMPLX_FREE =  43
  PetscInt, parameter, public :: PRIMARY_ACTIVITY_COEF =   44
  PetscInt, parameter, public :: SECONDARY_ACTIVITY_COEF = 45
  PetscInt, parameter, public :: SC_FUGA_COEFF =           46
  PetscInt, parameter, public :: PRIMARY_KD =              47
  PetscInt, parameter, public :: TOTAL_SORBED =            48
  PetscInt, parameter, public :: TOTAL_SORBED_MOBILE =     49
  PetscInt, parameter, public :: COLLOID_MOBILE =          50
  PetscInt, parameter, public :: COLLOID_IMMOBILE =        51
  PetscInt, parameter, public :: AGE =                     52
  PetscInt, parameter, public :: STATE =                   53
  PetscInt, parameter, public :: PROCESS_ID =              54
  PetscInt, parameter, public :: ICE_SATURATION =          55
  PetscInt, parameter, public :: TOTAL_BULK =              56
  PetscInt, parameter, public :: ICE_DENSITY =             57
  PetscInt, parameter, public :: GAS_PRESSURE =            58
  PetscInt, parameter, public :: SECONDARY_TEMPERATURE =   59
  PetscInt, parameter, public :: SECONDARY_CONCENTRATION = 60
  PetscInt, parameter, public :: SEC_MIN_VOLFRAC =         61

  PetscInt, parameter, public :: SURFACE_LIQUID_HEAD        = 62
  PetscInt, parameter, public :: SURFACE_LIQUID_TEMPERATURE = 63
  PetscInt, parameter, public :: SURFACE_LIQUID_DENSITY     = 64

  PetscInt, parameter, public :: EH = 65
  PetscInt, parameter, public :: PE = 66
  PetscInt, parameter, public :: O2 = 67

  PetscInt, parameter, public :: PERMEABILITY_XY =         68
  PetscInt, parameter, public :: PERMEABILITY_YZ =         69
  PetscInt, parameter, public :: PERMEABILITY_XZ =         70

  PetscInt, parameter, public :: GEOMECH_DISP_X          = 71
  PetscInt, parameter, public :: GEOMECH_DISP_Y          = 72
  PetscInt, parameter, public :: GEOMECH_DISP_Z          = 73
  PetscInt, parameter, public :: STRAIN_XX               = 74
  PetscInt, parameter, public :: STRAIN_YY               = 75
  PetscInt, parameter, public :: STRAIN_ZZ               = 76
  PetscInt, parameter, public :: STRAIN_XY               = 77
  PetscInt, parameter, public :: STRAIN_YZ               = 78
  PetscInt, parameter, public :: STRAIN_ZX               = 79
  PetscInt, parameter, public :: STRESS_XX               = 80
  PetscInt, parameter, public :: STRESS_YY               = 81
  PetscInt, parameter, public :: STRESS_ZZ               = 82
  PetscInt, parameter, public :: STRESS_XY               = 83
  PetscInt, parameter, public :: STRESS_YZ               = 84
  PetscInt, parameter, public :: STRESS_ZX               = 85
  PetscInt, parameter, public :: GEOMECH_MATERIAL_ID     = 86
  PetscInt, parameter, public :: GEOMECH_REL_DISP_X      = 87
  PetscInt, parameter, public :: GEOMECH_REL_DISP_Y      = 88
  PetscInt, parameter, public :: GEOMECH_REL_DISP_Z      = 89

  PetscInt, parameter, public :: VOLUME                  = 90
  PetscInt, parameter, public :: TORTUOSITY              = 91
  PetscInt, parameter, public :: AIR_PRESSURE            = 92
  PetscInt, parameter, public :: CAPILLARY_PRESSURE      = 93
  PetscInt, parameter, public :: VAPOR_PRESSURE          = 94
  PetscInt, parameter, public :: SATURATION_PRESSURE     = 95
  PetscInt, parameter, public :: MAXIMUM_PRESSURE        = 96

  PetscInt, parameter, public :: EFFECTIVE_POROSITY      = 97
  PetscInt, parameter, public :: LIQUID_HEAD             = 98

  PetscInt, parameter, public :: GAS_CONCENTRATION       = 99

  PetscInt, parameter, public :: SEC_MIN_RATE            = 100
  PetscInt, parameter, public :: SEC_MIN_SI              = 101
  PetscInt, parameter, public :: RESIDUAL                = 102
  PetscInt, parameter, public :: MINERAL_POROSITY        = 103

  PetscInt, parameter, public :: SOIL_COMPRESSIBILITY    = 104
  PetscInt, parameter, public :: SOIL_REFERENCE_PRESSURE = 105

  PetscInt, parameter, public :: LIQUID_MASS_FRACTION    = 106
  PetscInt, parameter, public :: GAS_MASS_FRACTION       = 107

!  PetscInt, parameter, public :: LIQUID_VELOCITY_CELL_CENT = 58
!  PetscInt, parameter, public :: LIQUID_VELOCITY_CELL_FACE = 59

  PetscInt, parameter, public :: OIL_PRESSURE            = 108
  PetscInt, parameter, public :: OIL_SATURATION          = 109
  PetscInt, parameter, public :: OIL_DENSITY             = 110
  PetscInt, parameter, public :: OIL_DENSITY_MOL         = 111
  PetscInt, parameter, public :: OIL_ENERGY              = 112
  PetscInt, parameter, public :: OIL_MOBILITY            = 113

  PetscInt, parameter, public :: NATURAL_ID              = 114
  PetscInt, parameter, public :: REACTION_AUXILIARY      = 115

  PetscInt, parameter, public :: OIL_VISCOSITY           = 116

  PetscInt, parameter, public :: OIL_MOLE_FRACTION       = 115
  PetscInt, parameter, public :: SOLVENT_SATURATION      = 116
  PetscInt, parameter, public :: CAPILLARY_PRESSURE_OGC  = 117
  PetscInt, parameter, public :: FRACTURE                = 118

  PetscInt, parameter, public :: GAS_PERMEABILITY        = 119
  PetscInt, parameter, public :: GAS_PERMEABILITY_X      = 120
  PetscInt, parameter, public :: GAS_PERMEABILITY_Y      = 121
  PetscInt, parameter, public :: GAS_PERMEABILITY_Z      = 122

  PetscInt, parameter, public :: BUBBLE_POINT            = 123

end module Variables_module
