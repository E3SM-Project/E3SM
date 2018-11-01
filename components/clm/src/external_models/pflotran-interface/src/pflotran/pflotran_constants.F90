module PFLOTRAN_Constants_module

! IMPORTANT NOTE: This module can have no dependencies on other modules!!!
 
  implicit none

  private

#include "petsc/finclude/petscsys.h"
  ! MUST INCREMENT THIS NUMBER EVERYTIME A CHECKPOINT FILE IS MODIFIED TO PREVENT
  ! COMPATIBILITY ISSUES - geh.
  PetscInt, parameter, public :: CHECKPOINT_REVISION_NUMBER = 1
  
  PetscInt, parameter, public :: MAXSTRINGLENGTH = 512
  PetscInt, parameter, public :: MAXWORDLENGTH = 32
  PetscInt, parameter, public :: OUT_UNIT = 15
  PetscInt, parameter, public :: OUTPUT_UNIT = 16
  PetscInt, parameter, public :: IN_UNIT = 17
  ! If you increase MAX_IN_UNIT, you MUST ensure that no other units #
  ! lie between IN_UNIT and MAX_IN_UNIT, as these units are reserved
  ! for embedded input files.
  PetscInt, parameter, public :: MAX_IN_UNIT = 25
  PetscInt, parameter, public :: IUNIT_TEMP = 86
  ! EKG_UNIT = 87
  PetscInt, parameter, public :: INPUT_RECORD_UNIT = 88
  PetscInt, parameter, public :: HHISTORY_LENGTH = 1000
  ! HHISTORY_LENGTH is the length of the array used to store the differencing
  ! values h.
  
  ! formula weights
  PetscReal, parameter, public :: FMWNACL = 58.44277d0
  PetscReal, parameter, public :: FMWH2O = 18.01534d0  ! kg/kmol h2o
  PetscReal, parameter, public :: FMWCO2 = 44.0098d0
  PetscReal, parameter, public :: FMWAIR = 28.96d0
  PetscReal, parameter, public :: FMWOIL = 142.D0 ! used as deafault value

  ! constants
  PetscReal, parameter, public :: DAYS_PER_YEAR = 365.d0
  PetscReal, parameter, public :: H2O_CRITICAL_TEMPERATURE = 647.3d0  ! K
  PetscReal, parameter, public :: H2O_CRITICAL_PRESSURE = 22.064d6 ! Pa

  ! conversion factors
  PetscReal, parameter, public :: LOG_TO_LN = 2.30258509299d0
  PetscReal, parameter, public :: LN_TO_LOG = 0.434294481904d0  
  
  ! constants
                             ! from http://physics.nist.gov/cgi-bin/cuu/Value?r
  PetscReal, parameter, public :: IDEAL_GAS_CONSTANT = 8.31446d0 ! J/mol-K
  PetscReal, parameter, public :: HEAT_OF_FUSION = 3.34d5  ! J/kg
  PetscReal, parameter, public :: PI = 3.14159265359d0
  PetscReal, parameter, public :: FARADAY = 96485.3365d0 ! C/mol
  PetscReal, parameter, public :: EARTH_GRAVITY = 9.8068d0 ! m/s^2
  
  PetscInt, parameter, public :: ZERO_INTEGER = 0
  PetscInt, parameter, public :: ONE_INTEGER = 1
  PetscInt, parameter, public :: TWO_INTEGER = 2
  PetscInt, parameter, public :: THREE_INTEGER = 3
  PetscInt, parameter, public :: FOUR_INTEGER = 4
  PetscInt, parameter, public :: FIVE_INTEGER = 5
  PetscInt, parameter, public :: SIX_INTEGER = 6
  PetscInt, parameter, public :: SEVEN_INTEGER = 7
  PetscInt, parameter, public :: EIGHT_INTEGER = 8
  PetscInt, parameter, public :: NINE_INTEGER = 9
  PetscInt, parameter, public :: TEN_INTEGER = 10
  PetscInt, parameter, public :: ELEVEN_INTEGER = 11
  PetscInt, parameter, public :: TWELVE_INTEGER = 12
  PetscInt, parameter, public :: NEG_ONE_INTEGER = -1
  
  PetscMPIInt, parameter, public :: ZERO_INTEGER_MPI = ZERO_INTEGER
  PetscMPIInt, parameter, public :: ONE_INTEGER_MPI = ONE_INTEGER
  PetscMPIInt, parameter, public :: TWO_INTEGER_MPI = TWO_INTEGER
  PetscMPIInt, parameter, public :: THREE_INTEGER_MPI = THREE_INTEGER
  PetscMPIInt, parameter, public :: FOUR_INTEGER_MPI = FOUR_INTEGER
  PetscMPIInt, parameter, public :: FIVE_INTEGER_MPI = FIVE_INTEGER
  PetscMPIInt, parameter, public :: SIX_INTEGER_MPI = SIX_INTEGER
  PetscMPIInt, parameter, public :: SEVEN_INTEGER_MPI = SEVEN_INTEGER
  PetscMPIInt, parameter, public :: TEN_INTEGER_MPI = TEN_INTEGER
  PetscMPIInt, parameter, public :: ELEVEN_INTEGER_MPI = ELEVEN_INTEGER
  PetscMPIInt, parameter, public :: TWELVE_INTEGER_MPI = TWELVE_INTEGER
  PetscMPIInt, parameter, public :: MAXSTRINGLENGTH_MPI = MAXSTRINGLENGTH
  
  PetscInt, parameter, public :: X_DIRECTION = 1
  PetscInt, parameter, public :: Y_DIRECTION = 2
  PetscInt, parameter, public :: Z_DIRECTION = 3
  PetscInt, parameter, public :: LOWER = 1
  PetscInt, parameter, public :: UPPER = 2
  
  PetscInt, parameter, public :: TIME_NULL = 0
  PetscInt, parameter, public :: TIME_T = 1
  PetscInt, parameter, public :: TIME_TpDT = 2
  
  PetscInt, parameter, public :: SORPTION_LINEAR = 1
  PetscInt, parameter, public :: SORPTION_LANGMUIR = 2
  PetscInt, parameter, public :: SORPTION_FREUNDLICH  = 3
  
  ! Classes
  PetscInt, parameter, public :: NULL_CLASS = 0
  PetscInt, parameter, public :: FLOW_CLASS = 1
  PetscInt, parameter, public :: TRANSPORT_CLASS = 2
  
  ! Macros that are used as 'dm_index' values.  --RTM
  PetscInt, parameter, public :: ONEDOF = 1
  PetscInt, parameter, public :: NPHASEDOF = 2
  PetscInt, parameter, public :: THREENPDOF = 3
  PetscInt, parameter, public :: NFLOWDOF = 4
  PetscInt, parameter, public :: NTRANDOF = 5
  PetscInt, parameter, public :: SURF_ONEDOF = 6
  
  PetscInt, parameter, public :: GLOBAL = 1
  PetscInt, parameter, public :: LOCAL = 2
  PetscInt, parameter, public :: NATURAL = 3
  
  PetscInt, parameter, public :: NULL_MODE = 0
  
  ! flow modes
  PetscInt, parameter, public :: TH_MODE = 2
  
  ! transport modes
  PetscInt, parameter, public :: EXPLICIT_ADVECTION = 1
  
  ! condition types
  PetscInt, parameter, public :: NULL_CONDITION = 0
  PetscInt, parameter, public :: DIRICHLET_BC = 1
  PetscInt, parameter, public :: NEUMANN_BC = 2
  PetscInt, parameter, public :: DIRICHLET_ZERO_GRADIENT_BC = 3
  PetscInt, parameter, public :: ZERO_GRADIENT_BC = 4
  PetscInt, parameter, public :: HYDROSTATIC_BC = 5
  PetscInt, parameter, public :: SEEPAGE_BC = 6
  PetscInt, parameter, public :: MASS_RATE_SS = 7
  PetscInt, parameter, public :: VOLUMETRIC_RATE_SS = 8
  PetscInt, parameter, public :: SCALED_MASS_RATE_SS = 9
  PetscInt, parameter, public :: SCALED_VOLUMETRIC_RATE_SS = 10
  PetscInt, parameter, public :: CONCENTRATION_SS = 11
  PetscInt, parameter, public :: EQUILIBRIUM_SS = 12
  PetscInt, parameter, public :: CONDUCTANCE_BC = 13
  PetscInt, parameter, public :: UNIT_GRADIENT_BC = 14
  PetscInt, parameter, public :: SATURATION_BC = 15
  PetscInt, parameter, public :: HET_VOL_RATE_SS = 16
  PetscInt, parameter, public :: HET_MASS_RATE_SS = 17
  PetscInt, parameter, public :: HET_DIRICHLET_BC = 18
  PetscInt, parameter, public :: ENERGY_RATE_SS = 19
  PetscInt, parameter, public :: SCALED_ENERGY_RATE_SS = 20
  PetscInt, parameter, public :: HET_ENERGY_RATE_SS = 21
  PetscInt, parameter, public :: HET_SURF_SEEPAGE_BC = 22
  PetscInt, parameter, public :: SPILLOVER_BC = 23
  PetscInt, parameter, public :: SURFACE_DIRICHLET = 33
  PetscInt, parameter, public :: SURFACE_ZERO_GRADHEIGHT = 34
  PetscInt, parameter, public :: SURFACE_SPILLOVER = 35
  PetscInt, parameter, public :: HET_SEEPAGE_BC = 36
  PetscInt, parameter, public :: HET_CONDUCTANCE_BC = 37
  
  ! source/sink scaling options
  PetscInt, parameter, public :: SCALE_BY_PERM = 1
  PetscInt, parameter, public :: SCALE_BY_NEIGHBOR_PERM = 2
  PetscInt, parameter, public :: SCALE_BY_VOLUME = 3
  
  ! connection types
  PetscInt, parameter, public :: INTERNAL_CONNECTION_TYPE = 1
  PetscInt, parameter, public :: BOUNDARY_CONNECTION_TYPE = 2
  PetscInt, parameter, public :: INITIAL_CONNECTION_TYPE = 3
  PetscInt, parameter, public :: SRC_SINK_CONNECTION_TYPE = 4
  
  ! dofs for each mode
  PetscInt, parameter, public :: TH_PRESSURE_DOF = 1
  PetscInt, parameter, public :: TH_TEMPERATURE_DOF = 2
  PetscInt, parameter, public :: TH_CONDUCTANCE_DOF = 3

  ! phase ids
  PetscInt, parameter, public :: LIQUID_PHASE = 1
  PetscInt, parameter, public :: GAS_PHASE = 2
  PetscInt, parameter, public :: SOLID_PHASE = 3
  
  ! approaches to coupling reactive transport
  PetscInt, parameter, public :: GLOBAL_IMPLICIT = 0
  PetscInt, parameter, public :: OPERATOR_SPLIT = 1
  
  ! ids of non-petsc arrays
  PetscInt, parameter, public :: MATERIAL_ID_ARRAY = 1
  PetscInt, parameter, public :: SATURATION_FUNCTION_ID_ARRAY = 2
  
  ! interpolation methods
  PetscInt, parameter, public :: INTERPOLATION_NULL = 0
  PetscInt, parameter, public :: INTERPOLATION_STEP = 1
  PetscInt, parameter, public :: INTERPOLATION_LINEAR = 2
  
  ! surface/subsurface flags
  PetscInt, parameter, public :: SUBSURFACE = 0
  PetscInt, parameter, public :: SURFACE    = 1
  
  PetscInt, parameter, public :: DECOUPLED     = 0
  PetscInt, parameter, public :: SEQ_COUPLED = 1
  PetscInt, parameter, public :: FULLY_COUPLED = 2
  
  PetscInt, parameter, public :: KINEMATIC_WAVE = 1
  PetscInt, parameter, public :: DIFFUSION_WAVE = 2
  
  PetscReal, parameter, public :: MIN_SURFACE_WATER_HEIGHT = 1.0d-14

  ! print secondary continuum variable ids
  PetscInt, parameter, public :: PRINT_SEC_TEMP =           0
  PetscInt, parameter, public :: PRINT_SEC_CONC =           1
  PetscInt, parameter, public :: PRINT_SEC_MIN_VOLFRAC =    2
  PetscInt, parameter, public :: PRINT_SEC_MIN_RATE =       3
  PetscInt, parameter, public :: PRINT_SEC_MIN_SI =         4
  
  PetscInt, parameter, public :: PROCEED = 0
  PetscInt, parameter, public :: DONE = 1
  PetscInt, parameter, public :: FAIL = 2

  ! Grid type
  PetscInt, parameter, public :: NULL_GRID = 0
  PetscInt, parameter, public :: STRUCTURED_GRID = 1
  PetscInt, parameter, public :: UNSTRUCTURED_GRID = 2
  PetscInt, parameter, public :: IMPLICIT_UNSTRUCTURED_GRID = 3
  PetscInt, parameter, public :: EXPLICIT_UNSTRUCTURED_GRID = 4
  PetscInt, parameter, public :: POLYHEDRA_UNSTRUCTURED_GRID = 5
  PetscInt, parameter, public :: ONE_DIM_GRID = 1
  PetscInt, parameter, public :: TWO_DIM_GRID = 2
  PetscInt, parameter, public :: THREE_DIM_GRID = 3
  PetscInt, parameter, public :: VERTEX_CENTERED_OUTPUT_MESH = 1
  PetscInt, parameter, public :: CELL_CENTERED_OUTPUT_MESH = 2

  ! Macros that are used as 'vscatter_index' values
  PetscInt, parameter, public :: SURF_TO_SUBSURF = 1
  PetscInt, parameter, public :: SUBSURF_TO_SURF = 2
  
  ! Ice/water/vapor partitioning model
  PetscInt, parameter, public :: PAINTER_EXPLICIT = 1
  PetscInt, parameter, public :: PAINTER_KARRA_IMPLICIT = 2
  PetscInt, parameter, public :: PAINTER_KARRA_EXPLICIT = 3
  PetscInt, parameter, public :: DALL_AMICO = 4
  PetscInt, parameter, public :: PAINTER_KARRA_EXPLICIT_NOCRYO = 5
  PetscInt, parameter, public :: PAINTER_KARRA_EXPLICIT_SMOOTH = 6

  ! Relative permeability averaging
  PetscInt, parameter, public :: UPWIND = 1
  PetscInt, parameter, public :: HARMONIC = 2
  PetscInt, parameter, public :: DYNAMIC_HARMONIC = 3

  ! uninitialized values
  PetscInt, parameter, public :: UNINITIALIZED_INTEGER = -999
  PetscReal, parameter, public :: UNINITIALIZED_DOUBLE = -999.d0

  ! global solver convergence criteria
  PetscInt, parameter, public :: CONVERGENCE_OFF = -999
  PetscInt, parameter, public :: CONVERGENCE_CUT_TIMESTEP = -1
  PetscInt, parameter, public :: CONVERGENCE_KEEP_ITERATING = 0
  PetscInt, parameter, public :: CONVERGENCE_FORCE_ITERATION = 1
  PetscInt, parameter, public :: CONVERGENCE_CONVERGED = 2
  
  ! Dummy value
  PetscReal, parameter, public :: DUMMY_VALUE = UNINITIALIZED_DOUBLE
  
  interface Uninitialized
    module procedure UninitializedInteger
    module procedure UninitializedDouble
  end interface
  
  interface Initialized
    module procedure InitializedInteger
    module procedure InitializedDouble
  end interface
  
  public :: Initialized, &
            Uninitialized, &
            UninitializedMessage
  
contains

! ************************************************************************** !

function InitializedInteger(value)
  ! 
  ! Tests whether a variable is initialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none
  
  PetscInt :: value
  PetscBool :: InitializedInteger
  
  InitializedInteger = .not.Uninitialized(value)
  
end function InitializedInteger


! ************************************************************************** !

function UninitializedInteger(value)
  ! 
  ! Tests whether a variable is uninitialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none
  
  PetscInt :: value
  PetscBool :: UninitializedInteger
  
  UninitializedInteger = (value == UNINITIALIZED_INTEGER)
  
end function UninitializedInteger

! ************************************************************************** !

function InitializedDouble(value)
  ! 
  ! Tests whether a variable is initialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none
  
  PetscReal :: value
  PetscBool :: InitializedDouble

  InitializedDouble = .not.Uninitialized(value)
  
end function InitializedDouble

! ************************************************************************** !

function UninitializedDouble(value)
  ! 
  ! Tests whether a variable is uninitialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none
  
  PetscReal :: value
  PetscBool :: UninitializedDouble

  UninitializedDouble = (dabs(value-UNINITIALIZED_DOUBLE) < 1.d-20)
  
end function UninitializedDouble

! ************************************************************************** !

function UninitializedMessage(variable_name,routine_name)
  ! 
  ! Tests whether a variable is uninitialized based orginally being set to
  ! the value UNINITIALIZED_INTEGER
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  implicit none
  
  character(len=*) :: variable_name
  character(len=*) :: routine_name
  
  character(len=MAXSTRINGLENGTH) :: UninitializedMessage
  
  if (len_trim(routine_name) > 1) then
    UninitializedMessage = trim(variable_name) // &
                           ' uninitialized in ' // &
                           trim(routine_name) // '.'
  else
    UninitializedMessage = trim(variable_name) // &
                           ' uninitialized.'
  endif
  
end function UninitializedMessage

end module PFLOTRAN_Constants_module
