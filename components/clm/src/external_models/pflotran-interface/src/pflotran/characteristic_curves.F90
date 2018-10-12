module Characteristic_Curves_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  PetscReal, parameter :: DEFAULT_PCMAX = 1.d9
  
  type :: polynomial_type
    PetscReal :: low
    PetscReal :: high
    PetscReal :: coefficients(4)
  end type polynomial_type
 
  ! Begin Saturation Functions ------------------------------------------------
  type :: sat_func_base_type
    type(polynomial_type), pointer :: sat_poly
    type(polynomial_type), pointer :: pres_poly
    PetscReal :: Sr
    PetscReal :: pcmax
#ifdef SMOOTHING2
    type(polynomial_type), pointer :: sat_poly2      ! dry-end of the curve  ! added by F.-M. Yuan (2017-03-10)
    type(polynomial_type), pointer :: pres_poly2     ! dry-end of the curve  ! added by F.-M. Yuan (2017-03-10)
#endif
  contains
    procedure, public :: Init => SFBaseInit
    procedure, public :: Verify => SFBaseVerify
    procedure, public :: Test => SFBaseTest
    procedure, public :: SetupPolynomials => SFBaseSetupPolynomials
    procedure, public :: CapillaryPressure => SFBaseCapillaryPressure
    procedure, public :: Saturation => SFBaseSaturation
    !added pc function if ice exists
    procedure, public :: IceCapillaryPressure => SF_Ice_CapillaryPressure
  end type sat_func_base_type
  ! Default
  type, public, extends(sat_func_base_type) :: sat_func_default_type
  contains
    procedure, public :: Verify => SFDefaultVerify
    procedure, public :: CapillaryPressure => SFDefaultCapillaryPressure
    procedure, public :: Saturation => SFDefaultSaturation
  end type sat_func_default_type  
  type, public, extends(sat_func_base_type) :: sat_func_constant_type
    PetscReal :: constant_capillary_pressure
    PetscReal :: constant_saturation
  contains
    procedure, public :: Verify => SFConstantVerify
    procedure, public :: CapillaryPressure => SFConstantCapillaryPressure
    procedure, public :: Saturation => SFConstantSaturation
  end type sat_func_constant_type  
  type, public, extends(sat_func_base_type) :: sat_func_VG_type
    PetscReal :: alpha
    PetscReal :: m
  contains
    procedure, public :: Init => SF_VG_Init
    procedure, public :: Verify => SF_VG_Verify
    procedure, public :: CapillaryPressure => SF_VG_CapillaryPressure
    procedure, public :: Saturation => SF_VG_Saturation
  end type sat_func_VG_type  
  type, public, extends(sat_func_base_type) :: sat_func_BC_type
    PetscReal :: alpha
    PetscReal :: lambda
  contains
    procedure, public :: Init => SF_BC_Init
    procedure, public :: Verify => SF_BC_Verify
    procedure, public :: SetupPolynomials => SF_BC_SetupPolynomials
    procedure, public :: CapillaryPressure => SF_BC_CapillaryPressure
    procedure, public :: Saturation => SF_BC_Saturation
  end type sat_func_BC_type  
  type, public, extends(sat_func_base_type) :: sat_func_Linear_type
    PetscReal :: alpha
  contains
    procedure, public :: Init => SF_Linear_Init
    procedure, public :: Verify => SF_Linear_Verify
    procedure, public :: CapillaryPressure => SF_Linear_CapillaryPressure
    procedure, public :: Saturation => SF_Linear_Saturation
  end type sat_func_Linear_type
  ! BRAGFLO KRP1 modified van Genuchten-Parker Model
  ! very similar to van Genuchten-Mualem, but includes Sgr
  type, public, extends(sat_func_base_type) :: sat_func_BF_KRP1_type
    PetscReal :: Srg
    PetscReal :: alpha
    PetscReal :: m
  contains
    procedure, public :: Init => SF_BF_KRP1_Init
    procedure, public :: Verify => SF_BF_KRP1_Verify
    procedure, public :: CapillaryPressure => SF_BF_KRP1_CapillaryPressure
    procedure, public :: Saturation => SF_BF_KRP1_Saturation
  end type sat_func_BF_KRP1_type
  ! BRAGFLO KRP5 modified Linear Model
  ! very similar to Linear-Burdine, but includes Sgr
  type, public, extends(sat_func_base_type) :: sat_func_BF_KRP5_type
    PetscReal :: alpha
    PetscReal :: Srg
  contains
    procedure, public :: Init => SF_BF_KRP5_Init
    procedure, public :: Verify => SF_BF_KRP5_Verify
    procedure, public :: CapillaryPressure => SF_BF_KRP5_CapillaryPressure
    procedure, public :: Saturation => SF_BF_KRP5_Saturation
  end type sat_func_BF_KRP5_type
  ! BRAGFLO KRP9 modified Brooks-Corey Model
  type, public, extends(sat_func_base_type) :: sat_func_BF_KRP9_type
  contains
    procedure, public :: Init => SF_BF_KRP9_Init
    procedure, public :: Verify => SF_BF_KRP9_Verify
    procedure, public :: CapillaryPressure => SF_BF_KRP9_CapillaryPressure
    procedure, public :: Saturation => SF_BF_KRP9_Saturation
  end type sat_func_BF_KRP9_type
  ! BRAGFLO KRP4 modified Brooks-Corey Model
  type, public, extends(sat_func_BC_type) :: sat_func_BF_KRP4_type
    PetscReal :: Srg
    PetscInt :: pcmax_flag
  contains
    procedure, public :: Verify => SF_BF_KRP4_Verify
    procedure, public :: CapillaryPressure => SF_BF_KRP4_CapillaryPressure
    procedure, public :: Saturation => SF_BF_KRP4_Saturation
  end type sat_func_BF_KRP4_type
  ! BRAGFLO KRP11 
  type, public, extends(sat_func_base_type) :: sat_func_BF_KRP11_type
  contains
    procedure, public :: Init => SF_BF_KRP11_Init
    procedure, public :: Verify => SF_BF_KRP11_Verify
    procedure, public :: CapillaryPressure => SF_BF_KRP11_CapillaryPressure
    procedure, public :: Saturation => SF_BF_KRP11_Saturation
  end type sat_func_BF_KRP11_type 
  ! BRAGFLO KRP12 modified Brooks-Corey Model
  type, public, extends(sat_func_BC_type) :: sat_func_BF_KRP12_type
    PetscReal :: Srg
    PetscReal :: socmin
    PetscReal :: soceffmin
  contains
    procedure, public :: Verify => SF_BF_KRP12_Verify
    procedure, public :: CapillaryPressure => SF_BF_KRP12_CapillaryPressure
  end type sat_func_BF_KRP12_type
  ! modified Kosugi Model (Malama & Kuhlman, 2015)
  type, public, extends(sat_func_base_type) :: sat_func_mK_type
    PetscReal :: sigmaz, muz
    PetscReal :: rmax, r0
    PetscInt :: nparam
  contains
    procedure, public :: Verify => SF_mK_Verify
    procedure, public :: CapillaryPressure => SF_mK_CapillaryPressure
    procedure, public :: Saturation => SF_mK_Saturation
  end type sat_func_mK_type
  ! End Saturation Functions --------------------------------------------------

  ! Begin Relative Permeability Functions -------------------------------------
  type :: rel_perm_func_base_type
    type(polynomial_type), pointer :: poly
#ifdef SMOOTHING2
    type(polynomial_type), pointer :: poly2     ! dry-end of the curve  ! added by F.-M. Yuan (2017-03-10)
#endif
    PetscReal :: Sr
  contains
    procedure, public :: Init => RPFBaseInit
    procedure, public :: Verify => RPFBaseVerify
    procedure, public :: Test => RPF_Base_Test
    procedure, public :: SetupPolynomials => RPFBaseSetupPolynomials
    procedure, public :: RelativePermeability => RPF_Base_RelPerm
  end type rel_perm_func_base_type
  ! Default
  type, public, extends(rel_perm_func_base_type) :: rel_perm_func_default_type
  contains
    procedure, public :: Verify => RPFDefaultVerify
    procedure, public :: RelativePermeability => RPF_DefaultRelPerm
  end type rel_perm_func_default_type
  ! Mualem-VG-liq
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_VG_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_Mualem_VG_Liq_Init
    procedure, public :: Verify => RPF_Mualem_VG_Liq_Verify
    procedure, public :: SetupPolynomials => RPF_Mualem_SetupPolynomials
    procedure, public :: RelativePermeability => RPF_Mualem_VG_Liq_RelPerm
  end type rpf_Mualem_VG_liq_type
  ! Mualem-VG-gas
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_VG_gas_type
    PetscReal :: m
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Mualem_VG_Gas_Init
    procedure, public :: Verify => RPF_Mualem_VG_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_VG_Gas_RelPerm
  end type rpf_Mualem_VG_gas_type
  ! since the TOUGH2_Corey relative permeability function (IRP=7 in 
  ! TOUGH2 manual) calculates relative perm as a function of the 
  ! Mualem-based  liquid relative permeability when Srg = 0., we extend 
  ! the rpf_Mualem_type to save code
  type, public, extends(rpf_Mualem_VG_liq_type) :: rpf_TOUGH2_IRP7_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_TOUGH2_IRP7_Gas_Init
    procedure, public :: Verify => RPF_TOUGH2_IRP7_Gas_Verify
    procedure, public :: RelativePermeability => RPF_TOUGH2_IRP7_Gas_RelPerm
  end type rpf_TOUGH2_IRP7_gas_type
  ! Burdine-BC
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_BC_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_Burdine_BC_Liq_Init
    procedure, public :: Verify => RPF_Burdine_BC_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_BC_Liq_RelPerm
    procedure, public :: SetupPolynomials => RPF_Burdine_BC_Liq_SetupPolynomials  ! added by F.-M. Yuan (2017-03-10)
  end type rpf_Burdine_BC_liq_type
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_BC_gas_type
    PetscReal :: lambda
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Burdine_BC_Gas_Init
    procedure, public :: Verify => RPF_Burdine_BC_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_BC_Gas_RelPerm
  end type rpf_Burdine_BC_gas_type
  ! Mualem-BC
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_BC_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_Mualem_BC_Liq_Init
    procedure, public :: Verify => RPF_Mualem_BC_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_BC_Liq_RelPerm
  end type rpf_MUALEM_BC_liq_type
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_BC_gas_type
    PetscReal :: lambda
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Mualem_BC_Gas_Init
    procedure, public :: Verify => RPF_Mualem_BC_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_BC_Gas_RelPerm
  end type rpf_Mualem_BC_gas_type
  ! Burdine-VG
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_VG_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_Burdine_VG_Liq_Init
    procedure, public :: Verify => RPF_Burdine_VG_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_VG_Liq_RelPerm
  end type rpf_Burdine_VG_liq_type
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_VG_gas_type
    PetscReal :: m
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Burdine_VG_Gas_Init
    procedure, public :: Verify => RPF_Burdine_VG_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_VG_Gas_RelPerm
  end type rpf_Burdine_VG_gas_type
  ! Mualem-Linear
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_Linear_liq_type
    PetscReal :: pcmax
    PetscReal :: alpha
  contains
    procedure, public :: Init => RPF_Mualem_Linear_Liq_Init
    procedure, public :: Verify => RPF_Mualem_Linear_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_Linear_Liq_RelPerm
  end type rpf_Mualem_Linear_liq_type
  type, public, extends(rpf_Mualem_Linear_liq_type) :: & 
                        rpf_Mualem_Linear_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Mualem_Linear_Gas_Init
    procedure, public :: Verify => RPF_Mualem_Linear_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_Linear_Gas_RelPerm
  end type rpf_Mualem_Linear_gas_type
  ! Burdine-Linear
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_Linear_liq_type
  contains
    procedure, public :: Init => RPF_Burdine_Linear_Liq_Init
    procedure, public :: Verify => RPF_Burdine_Linear_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_Linear_Liq_RelPerm
  end type rpf_Burdine_Linear_liq_type
  type, public, extends(rel_perm_func_base_type) :: & 
                        rpf_Burdine_Linear_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Burdine_Linear_Gas_Init
    procedure, public :: Verify => RPF_Burdine_Linear_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_Linear_Gas_RelPerm
  end type rpf_Burdine_Linear_gas_type
  ! BRAGFLO KRP1 modified van Genuchten Model
  ! relperm equations for KRP1 is identical to Mualem van Genuchten formulation
  type, public, extends(rpf_Mualem_VG_liq_type) :: rpf_BRAGFLO_KRP1_liq_type
  contains
  end type rpf_BRAGFLO_KRP1_liq_type
  type, public, extends(rpf_Mualem_VG_gas_type) :: rpf_BRAGFLO_KRP1_gas_type
  contains
  end type rpf_BRAGFLO_KRP1_gas_type
  ! Burdine-Linear
  type, public, extends(rel_perm_func_base_type) :: rpf_BRAGFLO_KRP5_liq_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_BRAGFLO_KRP5_Liq_Init
    procedure, public :: Verify => RPF_BRAGFLO_KRP5_Liq_Verify
    procedure, public :: RelativePermeability => RPF_BRAGFLO_KRP5_Liq_RelPerm
  end type rpf_BRAGFLO_KRP5_liq_type
  type, public, extends(rpf_Burdine_Linear_gas_type) :: rpf_BRAGFLO_KRP5_gas_type
  contains
  end type rpf_BRAGFLO_KRP5_gas_type
  ! BRAGFLO KRP9
  type, public, extends(rel_perm_func_base_type) :: rpf_BRAGFLO_KRP9_liq_type
  contains
    procedure, public :: Init => RPF_BRAGFLO_KRP9_Liq_Init
    procedure, public :: Verify => RPF_BRAGFLO_KRP9_Liq_Verify
    procedure, public :: RelativePermeability => RPF_BRAGFLO_KRP9_Liq_RelPerm
  end type rpf_BRAGFLO_KRP9_liq_type
  type, public, extends(rpf_BRAGFLO_KRP9_liq_type) :: & 
                        rpf_BRAGFLO_KRP9_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_BRAGFLO_KRP9_Gas_Init
    procedure, public :: Verify => RPF_BRAGFLO_KRP9_Gas_Verify
    procedure, public :: RelativePermeability => RPF_BRAGFLO_KRP9_Gas_RelPerm
  end type rpf_BRAGFLO_KRP9_gas_type
  ! BRAGFLO KRP4 modified Brooks-Corey Model
  ! relperm equations for KRP4 is identical to Burdine Brooks Corey
  ! formulation, but with different conditions in Gas RelPerm
  type, public, extends(rpf_Burdine_BC_liq_type) :: rpf_BRAGFLO_KRP4_liq_type
  contains
  end type rpf_BRAGFLO_KRP4_liq_type
  type, public, extends(rpf_Burdine_BC_gas_type) :: rpf_BRAGFLO_KRP4_gas_type
  contains
    procedure, public :: RelativePermeability => RPF_BRAGFLO_KRP4_Gas_RelPerm
  end type rpf_BRAGFLO_KRP4_gas_type
  ! BRAGFLO KRP11
  type, public, extends(rel_perm_func_base_type) :: rpf_BRAGFLO_KRP11_liq_type
    PetscReal :: tolc
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_BRAGFLO_KRP11_Liq_Init
    procedure, public :: Verify => RPF_BRAGFLO_KRP11_Liq_Verify
    procedure, public :: RelativePermeability => RPF_BRAGFLO_KRP11_Liq_RelPerm
  end type rpf_BRAGFLO_KRP11_liq_type
  type, public, extends(rpf_BRAGFLO_KRP11_liq_type) :: & 
                        rpf_BRAGFLO_KRP11_gas_type
  contains
    procedure, public :: RelativePermeability => RPF_BRAGFLO_KRP11_Gas_RelPerm
  end type rpf_BRAGFLO_KRP11_gas_type
  ! BRAGFLO KRP12 modified Brooks-Corey Model
  ! relperm equations for KRP12 is identical to Burdine Brooks Corey
  ! formulation, but with different conditions and truncations
  type, public, extends(rpf_Burdine_BC_liq_type) :: rpf_BRAGFLO_KRP12_liq_type
  contains
    procedure, public :: RelativePermeability => RPF_BRAGFLO_KRP12_Liq_RelPerm
  end type rpf_BRAGFLO_KRP12_liq_type
  type, public, extends(rpf_Burdine_BC_gas_type) :: rpf_BRAGFLO_KRP12_gas_type
  contains
    procedure, public :: RelativePermeability => RPF_BRAGFLO_KRP12_Gas_RelPerm
  end type rpf_BRAGFLO_KRP12_gas_type
  ! Oil relative permeability functions
  type, public, extends(rel_perm_func_base_type) :: rpf_TOUGH2_Linear_oil_type
    PetscReal :: Sro !
  contains
    procedure, public :: Init => RPF_TOUGH2_Linear_Oil_Init 
    procedure, public :: Verify => RPF_TOUGH2_Linear_Oil_Verify
    procedure, public :: RelativePermeability => RPF_TOUGH2_Linear_Oil_RelPerm
  end type rpf_TOUGH2_Linear_Oil_type
  type, public, extends(rel_perm_func_base_type) :: RPF_Mod_BC_type
    PetscReal :: m   !exponential coeff. 
    PetscReal :: Srg 
    PetscReal :: Sro
    PetscReal :: kr_max
  contains
    procedure, public :: Init => RPF_Mod_BC_Init 
    procedure, public :: Verify => RPF_Mod_BC_Verify
    procedure, public :: SetupPolynomials => RPF_Mod_BC_SetupPolynomials
  end type RPF_Mod_BC_type
  type, public, extends(RPF_Mod_BC_type) :: RPF_Mod_BC_liq_type
  contains
    procedure, public :: RelativePermeability => RPF_Mod_BC_Liq_RelPerm
  end type RPF_Mod_BC_liq_type
  type, public, extends(RPF_Mod_BC_type) :: RPF_Mod_BC_oil_type
  contains
    procedure, public :: RelativePermeability => RPF_Mod_BC_Oil_RelPerm
  end type RPF_Mod_BC_oil_type
  ! Constant: for running tests with a fixed relative permeability
  type, public, extends(rel_perm_func_base_type) :: rel_perm_func_constant_type
    PetscReal :: kr
  contains
    procedure, public :: Verify => RPFConstantVerify
    procedure, public :: RelativePermeability => RPF_ConstantRelPerm
  end type rel_perm_func_constant_type
  ! modified Kosugi (Malama & Kuhlman, 2015) for liquid and gas
  type, public, extends(rel_perm_func_base_type) :: rpf_mK_liq_type
    PetscReal :: sigmaz
  contains
    procedure, public :: Verify => RPF_mK_Liq_Verify
    procedure, public :: RelativePermeability => RPF_mK_Liq_RelPerm
  end type rpf_mK_liq_type
  type, public, extends(rel_perm_func_base_type) :: rpf_mK_gas_type
    PetscReal :: Srg
    PetscReal :: sigmaz
  contains
    procedure, public :: Verify => RPF_mK_Gas_Verify
    procedure, public :: RelativePermeability => RPF_mK_Gas_RelPerm
  end type rpf_mK_gas_type
  ! End Relative Permeability Functions ---------------------------------------
 
  type, public :: characteristic_curves_type
    character(len=MAXWORDLENGTH) :: name
    PetscBool :: print_me
    PetscBool :: test
    class(sat_func_base_type), pointer :: saturation_function
    class(rel_perm_func_base_type), pointer :: liq_rel_perm_function
    class(rel_perm_func_base_type), pointer :: gas_rel_perm_function
    class(rel_perm_func_base_type), pointer :: oil_rel_perm_function
    class(characteristic_curves_type), pointer :: next
  end type characteristic_curves_type
  
  type, public :: characteristic_curves_ptr_type
    class(characteristic_curves_type), pointer :: ptr
  end type characteristic_curves_ptr_type 
  
  public :: CharacteristicCurvesCreate, &
            CharacteristicCurvesRead, &
            CharacteristicCurvesAddToList, &
            CharCurvesConvertListToArray, &
            CharacteristicCurvesGetID, &
            CharCurvesGetGetResidualSats, &
            CharacteristicCurvesDestroy, &
  ! required to be public for unit tests - Heeho Park
            SF_Constant_Create, &
            SF_VG_Create, &
            SF_BC_Create, &
            SF_Linear_Create, &
            SF_BF_KRP1_Create, &
            SF_BF_KRP5_Create, &
            SF_BF_KRP9_Create, &
            SF_BF_KRP4_Create, &
            SF_BF_KRP11_Create, &
            SF_BF_KRP12_Create, &
            SF_mK_Create, &
            RPF_Mualem_VG_Liq_Create, &
            RPF_Mualem_VG_Gas_Create, &
            RPF_Burdine_BC_Liq_Create, &
            RPF_Burdine_BC_Gas_Create, &
            RPF_TOUGH2_IRP7_Gas_Create, &
            RPF_Mualem_BC_Liq_Create, &
            RPF_Mualem_BC_Gas_Create, &
            RPF_Burdine_VG_Liq_Create, &
            RPF_Burdine_VG_Gas_Create, &
            RPF_Mualem_Linear_Liq_Create, &
            RPF_Mualem_Linear_Gas_Create, &
            RPF_Burdine_Linear_Liq_Create, &
            RPF_Burdine_Linear_Gas_Create, &
            RPF_BRAGFLO_KRP1_Liq_Create, &
            RPF_BRAGFLO_KRP1_Gas_Create, &
            RPF_BRAGFLO_KRP5_Liq_Create, &
            RPF_BRAGFLO_KRP5_Gas_Create, &
            RPF_BRAGFLO_KRP9_Liq_Create, &
            RPF_BRAGFLO_KRP9_Gas_Create, &
            RPF_BRAGFLO_KRP4_Liq_Create, &
            RPF_BRAGFLO_KRP4_Gas_Create, &
            RPF_BRAGFLO_KRP11_Liq_Create, &
            RPF_BRAGFLO_KRP11_Gas_Create, &
            RPF_BRAGFLO_KRP12_Liq_Create, &
            RPF_BRAGFLO_KRP12_Gas_Create, &
            RPF_mK_Liq_Create, &
            RPF_mK_Gas_Create, &
            PolynomialCreate

contains

! ************************************************************************** !

! Begin Characteristic Curves
function CharacteristicCurvesCreate()
  ! 
  ! Creates a characteristic curve object that holds parameters and pointers
  ! to functions for calculating saturation, capillary pressure, relative
  ! permeability, etc.
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/23/14
  ! 

  implicit none

  class(characteristic_curves_type), pointer :: CharacteristicCurvesCreate
  
  class(characteristic_curves_type), pointer :: characteristic_curves
  
  allocate(characteristic_curves)
  characteristic_curves%name = ''
  characteristic_curves%print_me = PETSC_FALSE
  characteristic_curves%test = PETSC_FALSE
  nullify(characteristic_curves%saturation_function)
  nullify(characteristic_curves%liq_rel_perm_function)
  nullify(characteristic_curves%gas_rel_perm_function)
  nullify(characteristic_curves%oil_rel_perm_function)
  nullify(characteristic_curves%next)

  CharacteristicCurvesCreate => characteristic_curves

end function CharacteristicCurvesCreate

! ************************************************************************** !

subroutine CharacteristicCurvesRead(this,input,option)
  ! 
  ! Reads in contents of a saturation_function card
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/21/09
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(characteristic_curves_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word, phase_keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  class(rel_perm_func_base_type), pointer :: rel_perm_function_ptr

  input%ierr = 0
  error_string = 'CHARACTERISTIC_CURVES'  
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
      case('SATURATION_FUNCTION')
! replacing read word that is capable of database lookup
!        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputReadWordDbaseCompatible(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'saturation_function_type', &
                           error_string)
        call StringToUpper(word)
        select case(word)
          case('CONSTANT')
            this%saturation_function => SF_Constant_Create()
          case('VAN_GENUCHTEN')
            this%saturation_function => SF_VG_Create()
          case('BROOKS_COREY')
            this%saturation_function => SF_BC_Create()
          case('LINEAR')
            this%saturation_function => SF_Linear_Create()
          case('BRAGFLO_KRP1')
            this%saturation_function => SF_BF_KRP1_Create()
          case('BRAGFLO_KRP4')
            this%saturation_function => SF_BF_KRP4_Create()
          case('BRAGFLO_KRP5')
            this%saturation_function => SF_BF_KRP5_Create()
          case('BRAGFLO_KRP9')
            this%saturation_function => SF_BF_KRP9_Create()
          case('BRAGFLO_KRP11')
            this%saturation_function => SF_BF_KRP11_Create()
          case('BRAGFLO_KRP12')
            this%saturation_function => SF_BF_KRP12_Create()
          case('MODIFIED_KOSUGI')
            this%saturation_function => SF_mK_Create()
          case default
            call InputKeywordUnrecognized(word,'SATURATION_FUNCTION',option)
        end select
        call SaturationFunctionRead(this%saturation_function,input,option)
      case('PERMEABILITY_FUNCTION')
        nullify(rel_perm_function_ptr)
        phase_keyword = 'NONE'
! replacing read word that is capable of database lookup
!        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputReadWordDbaseCompatible(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'permeability_function_type', &
                           error_string)
        call StringToUpper(word)
        select case(word)
          case('MUALEM','MUALEM_VG_LIQ')
            rel_perm_function_ptr => RPF_Mualem_VG_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MUALEM_VG_GAS')
            rel_perm_function_ptr => RPF_Mualem_VG_Gas_Create()
            phase_keyword = 'GAS'
          case('BURDINE','BURDINE_BC_LIQ')
            rel_perm_function_ptr => RPF_Burdine_BC_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BURDINE_BC_GAS')
            rel_perm_function_ptr => RPF_Burdine_BC_Gas_Create()
            phase_keyword = 'GAS'
          case('TOUGH2_IRP7_LIQ')
            rel_perm_function_ptr => RPF_Mualem_VG_Liq_Create()
            phase_keyword = 'LIQUID'
          case('TOUGH2_IRP7_GAS')
            rel_perm_function_ptr => RPF_TOUGH2_IRP7_Gas_Create()
            phase_keyword = 'GAS'
          case('MUALEM_BC_LIQ')
            rel_perm_function_ptr => RPF_Mualem_BC_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MUALEM_BC_GAS')
            rel_perm_function_ptr => RPF_Mualem_BC_Gas_Create()
            phase_keyword = 'GAS'
          case('BURDINE_VG_LIQ')
            rel_perm_function_ptr => RPF_Burdine_VG_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BURDINE_VG_GAS')
            rel_perm_function_ptr => RPF_Burdine_VG_Gas_Create()
            phase_keyword = 'GAS'
          case('MUALEM_LINEAR_LIQ')
            rel_perm_function_ptr => RPF_Mualem_Linear_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MUALEM_LINEAR_GAS')
            rel_perm_function_ptr => RPF_Mualem_Linear_Gas_Create()
            phase_keyword = 'GAS'
          case('BURDINE_LINEAR_LIQ')
            rel_perm_function_ptr => RPF_Burdine_Linear_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BURDINE_LINEAR_GAS')
            rel_perm_function_ptr => RPF_Burdine_Linear_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP1_LIQ')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP1_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP1_GAS')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP1_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP4_LIQ')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP4_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP4_GAS')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP4_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP5_LIQ')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP5_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP5_GAS')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP5_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP9_LIQ')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP9_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP9_GAS')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP9_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP11_LIQ')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP11_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP11_GAS')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP11_Gas_Create()
            phase_keyword = 'GAS'
          case('BRAGFLO_KRP12_LIQ')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP12_Liq_Create()
            phase_keyword = 'LIQUID'
          case('BRAGFLO_KRP12_GAS')
            rel_perm_function_ptr => RPF_BRAGFLO_KRP12_Gas_Create()
            phase_keyword = 'GAS'
          case('MODIFIED_KOSUGI_LIQ')
            rel_perm_function_ptr => RPF_mK_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MODIFIED_KOSUGI_GAS')
            rel_perm_function_ptr => RPF_mK_Gas_Create()
            phase_keyword = 'GAS'
          case('TOUGH2_LINEAR_OIL')
            rel_perm_function_ptr => RPF_TOUGH2_Linear_Oil_Create()
            phase_keyword = 'OIL'
          case('MOD_BC_LIQ')
            rel_perm_function_ptr => RPF_Mod_BC_Liq_Create()
            phase_keyword = 'LIQUID'
          case('MOD_BC_OIL')
            rel_perm_function_ptr => RPF_Mod_BC_Oil_Create()
            phase_keyword = 'OIL'
          case('CONSTANT')
            rel_perm_function_ptr => RPF_Constant_Create()
            ! phase_keyword = 'NONE'
          case default
            call InputKeywordUnrecognized(word,'PERMEABILITY_FUNCTION',option)
        end select
        call PermeabilityFunctionRead(rel_perm_function_ptr,phase_keyword, &
                                      input,option)
        ! if PHASE is specified, have to align correct pointer
        select case(phase_keyword)
          case('GAS')
            this%gas_rel_perm_function => rel_perm_function_ptr
          case('LIQUID')
            this%liq_rel_perm_function => rel_perm_function_ptr
          case('OIL')
            this%oil_rel_perm_function => rel_perm_function_ptr 
            ! PO: gas_rel_perm_fucntion initiated oil_rel_perm_function
            ! to pass the verification in CharacteristicCurvesVerify
            ! in case gas_rel_perm_function is not defined in the input
            ! We should change CharacteristicCurvesVerify instead
             this%gas_rel_perm_function => rel_perm_function_ptr
          case('NONE')
            option%io_buffer = 'PHASE has not been set for &
                               &CHARACTERISTIC_CURVES,PERMEABILITY_FUNCTION. &
                               &This is most likely a development issue, and &
                               &not an input deck mistake. Please e-mail &
                               &pflotran-dev [at] googlegroups [dot] com.' 
            call printErrMsg(option)
          case default
            call InputKeywordUnrecognized(word, &
              'PERMEABILITY_FUNCTION,PHASE',option)
        end select
      case('TEST') 
        this%test = PETSC_TRUE
      case('DEFAULT')
        this%saturation_function => SF_Default_Create()
        this%liq_rel_perm_function => RPF_Default_Create()
        this%gas_rel_perm_function => this%liq_rel_perm_function
      case default
        call InputKeywordUnrecognized(keyword,'charateristic_curves',option)
    end select 
  enddo
  
  call CharacteristicCurvesVerify(this,option)

end subroutine CharacteristicCurvesRead

! ************************************************************************** !

subroutine SaturationFunctionRead(saturation_function,input,option)
  !
  ! Reads in contents of a SATURATION_FUNCTION block
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(sat_func_base_type) :: saturation_function
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: smooth

  input%ierr = 0
  smooth = PETSC_FALSE
  error_string = 'CHARACTERISTIC_CURVES,SATURATION_FUNCTION,'
  select type(sf => saturation_function)
    class is(sat_func_constant_type)
      error_string = trim(error_string) // 'CONSTANT'
    class is(sat_func_VG_type)
      error_string = trim(error_string) // 'VAN_GENUCHTEN'
    class is(sat_func_BC_type)
      error_string = trim(error_string) // 'BROOKS_COREY'
    class is(sat_func_Linear_type)
      error_string = trim(error_string) // 'LINEAR'
    class is(sat_func_BF_KRP1_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP1'
    class is(sat_func_BF_KRP4_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP4'
    class is(sat_func_BF_KRP5_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP5'
    class is(sat_func_BF_KRP9_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP9'
    class is(sat_func_BF_KRP11_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP11'
    class is(sat_func_BF_KRP12_type)
      error_string = trim(error_string) // 'BRAGFLO_KRP12'
    class is(sat_func_mK_type)
      error_string = trim(error_string) // 'MODIFIED_KOSUGI'
  end select
  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   

    ! base 
    found = PETSC_TRUE
    select case(keyword)
      case('LIQUID_RESIDUAL_SATURATION') 
        call InputReadDouble(input,option,saturation_function%Sr)
        call InputErrorMsg(input,option,'LIQUID_RESIDUAL_SATURATION', &
                           error_string)
      case('MAX_CAPILLARY_PRESSURE') 
        call InputReadDouble(input,option,saturation_function%pcmax)
        call InputErrorMsg(input,option,'MAX_CAPILLARY_PRESSURE', &
                            error_string)
      case('SMOOTH')
        smooth = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select
    
    if (found) cycle
    
    select type(sf => saturation_function)
      class is(sat_func_constant_type)
        select case(keyword)
          case('CONSTANT_CAPILLARY_PRESSURE') 
            call InputReadDouble(input,option,sf%constant_capillary_pressure)
            call InputErrorMsg(input,option,'constant capillary pressure', &
                               error_string)
          case('CONSTANT_SATURATION') 
            call InputReadDouble(input,option,sf%constant_saturation)
            call InputErrorMsg(input,option,'constant saturation', &
                                error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'constant saturation function',option)
        end select
      class is(sat_func_VG_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,sf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'van Genuchten saturation function',option)
        end select
      class is(sat_func_BC_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'Brooks-Corey saturation function',option)
        end select
      class is(sat_func_Linear_type)
        select case(keyword)
          case('ALPHA')
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'Linear saturation function',option)
        end select
      class is(sat_func_BF_KRP1_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,sf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,sf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'BRAGFLO_KRP1 saturation function',option)
        end select
      class is(sat_func_BF_KRP5_type)
        select case(keyword)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,sf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'BRAGFLO_KRP5 saturation function',option)
        end select
      class is(sat_func_BF_KRP4_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,sf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case('KPC') 
            call InputReadInt(input,option,sf%pcmax_flag)
            call InputErrorMsg(input,option,'pcmax_flag, KPC',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'BRAGFLO_KRP4 saturation function',option)
        end select
      class is(sat_func_BF_KRP9_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'BRAGFLO_KRP9 saturation function',option)
        end select
      class is(sat_func_BF_KRP11_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'BRAGFLO_KRP11 saturation function',option)
        end select
      class is(sat_func_BF_KRP12_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,sf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,sf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,sf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case('SOCMIN') 
            call InputReadDouble(input,option,sf%socmin)
            call InputErrorMsg(input,option,'socmin',error_string)
          case('SOCEFFMIN') 
            call InputReadDouble(input,option,sf%soceffmin)
            call InputErrorMsg(input,option,'soceffmin',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'BRAGFLO_KRP12 saturation function',option)
          end select
        class is(sat_func_mK_type)
          select case(keyword)
            case('SIGMAZ')
              call InputReadDouble(input,option,sf%sigmaz)
              call InputErrorMsg(input,option,'sigmaz',error_string)
            case('MUZ')
              call InputReadDouble(input,option,sf%muz)
              call InputErrorMsg(input,option,'muz',error_string)
            case('RMAX')
              call InputReadDouble(input,option,sf%rmax)
              call InputErrorMsg(input,option,'rmax',error_string)
            case('R0')
              call InputReadDouble(input,option,sf%r0)
              call InputErrorMsg(input,option,'r0',error_string)
            case('NPARAM')
              call InputReadInt(input,option,sf%nparam)
              call InputErrorMsg(input,option,'nparam',error_string)
            case default
              call InputKeywordUnrecognized(keyword, &
                   'MODIFIED_KOSUGI saturation function',option)
          end select
      class default
        option%io_buffer = 'Read routine not implemented for ' &
                           // trim(error_string) // '.'
        call printErrMsg(option)
    end select
  enddo
  
  if (smooth) then
    call saturation_function%SetupPolynomials(option,error_string)
  endif

  select type(sf => saturation_function)
    class is(sat_func_constant_type)
      option%io_buffer = 'Constant saturation function is being used.'
      call printWrnMsg(option)
    class is(sat_func_VG_type)
    class is(sat_func_BC_type)
      if (.not.smooth) then
        option%io_buffer = 'Brooks-Corey saturation function is being used &
          &without SMOOTH option.'
        call printWrnMsg(option)
      endif
    class is(sat_func_Linear_type)
    class is(sat_func_BF_KRP1_type)
    class is(sat_func_BF_KRP5_type)
    class is(sat_func_BF_KRP4_type)
      if (.not.smooth) then
        option%io_buffer = 'Brooks-Corey saturation function is being used &
          &without SMOOTH option.'
        call printWrnMsg(option)
      endif
    class is(sat_func_BF_KRP12_type)
      if (.not.smooth) then
        option%io_buffer = 'Brooks-Corey saturation function is being used &
          &without SMOOTH option.'
        call printWrnMsg(option)
      endif
  end select

end subroutine SaturationFunctionRead

! ************************************************************************** !

subroutine PermeabilityFunctionRead(permeability_function,phase_keyword, &
                                    input,option)
  !
  ! Reads in contents of a PERMEABILITY_FUNCTION block
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none
  
  class(rel_perm_func_base_type) :: permeability_function
  character(len=MAXWORDLENGTH) :: phase_keyword
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, new_phase_keyword
  character(len=MAXSTRINGLENGTH) :: error_string
  PetscBool :: found
  PetscBool :: smooth

  input%ierr = 0
  smooth = PETSC_FALSE
  new_phase_keyword = 'NONE'
  error_string = 'CHARACTERISTIC_CURVES,PERMEABILITY_FUNCTION,'
  select type(rpf => permeability_function)
    class is(rpf_Mualem_VG_liq_type)
      error_string = trim(error_string) // 'MUALEM_VG_LIQ'
    class is(rpf_Mualem_VG_gas_type)
      error_string = trim(error_string) // 'MUALEM_VG_GAS'
    class is(rpf_Burdine_BC_liq_type)
      error_string = trim(error_string) // 'BURDINE_BC_LIQ'
    class is(rpf_Burdine_BC_gas_type)
      error_string = trim(error_string) // 'BURDINE_BC_GAS'
    class is(rpf_TOUGH2_IRP7_gas_type)
      error_string = trim(error_string) // 'TOUGH2_IRP7_GAS'
    class is(rpf_Mualem_BC_liq_type)
      error_string = trim(error_string) // 'MUALEM_BC_LIQ'
    class is(rpf_Mualem_BC_gas_type)
      error_string = trim(error_string) // 'MUALEM_BC_GAS'
    class is(rpf_Burdine_VG_liq_type)
      error_string = trim(error_string) // 'BURDINE_VG_LIQ'
    class is(rpf_Burdine_VG_gas_type)
      error_string = trim(error_string) // 'BURDINE_VG_GAS'
    class is(rpf_Mualem_Linear_liq_type)
      error_string = trim(error_string) // 'MUALEM_Linear_LIQ'
    class is(rpf_Mualem_Linear_gas_type)
      error_string = trim(error_string) // 'MUALEM_Linear_GAS'
    class is(rpf_Burdine_Linear_liq_type)
      error_string = trim(error_string) // 'BURDINE_Linear_LIQ'
    class is(rpf_Burdine_Linear_gas_type)
      error_string = trim(error_string) // 'BURDINE_Linear_GAS'
    class is(rpf_BRAGFLO_KRP1_liq_type)
      error_string = trim(error_string) // 'MUALEM_BF_KRP1_LIQ'
    class is(rpf_BRAGFLO_KRP1_gas_type)
      error_string = trim(error_string) // 'MUALEM_BF_KRP1_GAS'
    class is(rpf_BRAGFLO_KRP4_liq_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP4_LIQ'
    class is(rpf_BRAGFLO_KRP4_gas_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP4_GAS'  
    class is(rpf_BRAGFLO_KRP5_liq_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP5_LIQ'
    class is(rpf_BRAGFLO_KRP5_gas_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP5_GAS'
    class is(rpf_BRAGFLO_KRP9_liq_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP9_LIQ'
    class is(rpf_BRAGFLO_KRP9_gas_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP9_GAS'
    class is(rpf_BRAGFLO_KRP11_liq_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP11_LIQ'
    class is(rpf_BRAGFLO_KRP11_gas_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP11_GAS'
    class is(rpf_BRAGFLO_KRP12_liq_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP12_LIQ'
    class is(rpf_BRAGFLO_KRP12_gas_type)
      error_string = trim(error_string) // 'BURDINE_BF_KRP12_GAS'
    class is(rpf_mK_liq_type)
      error_string = trim(error_string) // 'MODIFIED_KOSUGI_LIQ'
    class is(rpf_mK_gas_type)
      error_string = trim(error_string) // 'MODIFIED_KOSUGI_GAS'
    class is(rpf_TOUGH2_Linear_oil_type)
      error_string = trim(error_string) // 'TOUGH2_Linear_OIL'
    class is(rpf_mod_BC_liq_type)
      error_string = trim(error_string) // 'Mod_BC_LIQ'
    class is(rpf_mod_BC_oil_type)
      error_string = trim(error_string) // 'Mod_BC_OIL'
    class is(rel_perm_func_constant_type)
      error_string = trim(error_string) // 'CONSTANT'
  end select

  do
    call InputReadPflotranString(input,option)
    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(keyword)   

    ! base 
    found = PETSC_TRUE
    select case(keyword)
      case('LIQUID_RESIDUAL_SATURATION') 
        call InputReadDouble(input,option,permeability_function%Sr)
        call InputErrorMsg(input,option,'residual_saturation',error_string)
      case('PHASE')
        call InputReadWord(input,option,new_phase_keyword,PETSC_TRUE)
        call InputErrorMsg(input,option,'phase',error_string)
        call StringToUpper(phase_keyword) 
      case('SMOOTH')
        smooth = PETSC_TRUE
      case default
        found = PETSC_FALSE
    end select
    if (found) cycle

    select type(rpf => permeability_function)
      class is(rpf_Mualem_VG_liq_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem van Genuchten liquid relative permeability function', &
              option)
        end select
      class is(rpf_Mualem_VG_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem van Genuchten gas relative permeability function', &
              option)
        end select
      class is(rpf_Burdine_BC_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine Brooks-Corey liquid relative permeability function', &
              option)
        end select
      class is(rpf_Burdine_BC_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine Brooks-Corey gas relative permeability function', &
              option)
        end select
      class is(rpf_TOUGH2_IRP7_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                   'TOUGH2 IRP7 gas relative permeability function',option)
        end select
      class is(rpf_Mualem_BC_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem Brooks-Corey liquid relative permeability function', &
              option)
        end select
      class is(rpf_Mualem_BC_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem Brooks-Corey gas relative permeability function', &
              option)
        end select
      class is(rpf_Burdine_VG_liq_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine van Genuchten liquid relative permeability function', &
              option)
        end select
      class is(rpf_Burdine_VG_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine van Genuchten gas relative permeability function', &
              option)
        end select
      class is(rpf_Mualem_Linear_liq_type)
        select case(keyword)
          case('MAX_CAPILLARY_PRESSURE') 
            call InputReadDouble(input,option,rpf%pcmax)
            call InputErrorMsg(input,option,'max_capillary_pressure',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,rpf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem Linear liquid relative permeability function', &
              option)
        end select
      class is(rpf_Mualem_Linear_gas_type)
        select case(keyword)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case('MAX_CAPILLARY_PRESSURE') 
            call InputReadDouble(input,option,rpf%pcmax)
            call InputErrorMsg(input,option,'max_capillary_pressure',error_string)
          case('ALPHA') 
            call InputReadDouble(input,option,rpf%alpha)
            call InputErrorMsg(input,option,'alpha',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mualem Linear gas relative permeability function', &
              option)
        end select
      class is(rpf_Burdine_Linear_liq_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine Linear liquid relative permeability function', &
              option)
        end select
      class is(rpf_Burdine_Linear_gas_type)
        select case(keyword)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Burdine Linear gas relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP1_liq_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP1 liquid relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP1_gas_type)
        select case(keyword)
          case('M') 
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP1 gas relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP5_liq_type)
        select case(keyword)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP5 liquid relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP5_gas_type)
        select case(keyword)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP5 gas relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP9_liq_type)
        select case(keyword)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP9 liq relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP9_gas_type)
        select case(keyword)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP9 gas relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP4_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP4 liq relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP4_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP4 gas relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP11_liq_type)
        select case(keyword)
          case('TOLC') 
            call InputReadDouble(input,option,rpf%tolc)
            call InputErrorMsg(input,option,'tolc',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP11 liq relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP11_gas_type)
        select case(keyword)
          case('TOLC') 
            call InputReadDouble(input,option,rpf%tolc)
            call InputErrorMsg(input,option,'tolc',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP11 gas relative permeability function', &
              option)
        end select  
      class is(rpf_BRAGFLO_KRP12_liq_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP4 liq relative permeability function', &
              option)
        end select
      class is(rpf_BRAGFLO_KRP12_gas_type)
        select case(keyword)
          case('LAMBDA') 
            call InputReadDouble(input,option,rpf%lambda)
            call InputErrorMsg(input,option,'lambda',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'BRAGFLO KRP4 gas relative permeability function', &
              option)
        end select
      class is(rpf_mK_liq_type)
        select case(keyword)
          case('SIGMAZ')
            call InputReadDouble(input,option,rpf%sigmaz)
            call InputErrorMsg(input,option,'sigmaz',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                 'MODIFIED_KOSUGI liquid relative permeability '//&
                 &'function',option)
        end select
      class is(rpf_mK_gas_type)
        select case(keyword)
          case('SIGMAZ')
            call InputReadDouble(input,option,rpf%sigmaz)
            call InputErrorMsg(input,option,'sigmaz',error_string)
          case('GAS_RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
                 'MODIFIED_KOSUGI gas relative permeability '//&
                 &'function',option)
        end select
      class is(rpf_TOUGH2_Linear_oil_type)
        select case(keyword)
          case('OIL_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Sro)
            call InputErrorMsg(input,option,'Sro',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'TOUGH2 LINEAR oil relative permeability function', &
              option)
        end select
      class is(rpf_mod_BC_liq_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m - power',error_string)
          case('OIL_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Sro)
            call InputErrorMsg(input,option,'Sro',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case('LIQUID_MAX_REL_PERM') 
            call InputReadDouble(input,option,rpf%kr_max)
            call InputErrorMsg(input,option,'kr_max',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mod BC liq relative permeability function', &
              option)
        end select
      class is(rpf_mod_BC_oil_type)
        select case(keyword)
          case('M')
            call InputReadDouble(input,option,rpf%m)
            call InputErrorMsg(input,option,'m - power',error_string)
          case('OIL_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Sro)
            call InputErrorMsg(input,option,'Sro',error_string)
          case('GAS_RESIDUAL_SATURATION') 
            call InputReadDouble(input,option,rpf%Srg)
            call InputErrorMsg(input,option,'Srg',error_string)
          case('OIL_MAX_REL_PERM') 
            call InputReadDouble(input,option,rpf%kr_max)
            call InputErrorMsg(input,option,'kr_max',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Mod BC oil relative permeability function', &
              option)
        end select
      class is(rel_perm_func_constant_type)
        select case(keyword)
          case('RESIDUAL_SATURATION')
            call InputReadDouble(input,option,rpf%Sr)
            call InputErrorMsg(input,option,'Sr',error_string)
          case('RELATIVE_PERMEABILITY') 
            call InputReadDouble(input,option,rpf%kr)
            call InputErrorMsg(input,option,'kr',error_string)
          case default
            call InputKeywordUnrecognized(keyword, &
              'Constant relative permeability function', &
              option)
        end select
      class default
        option%io_buffer = 'Read routine not implemented for relative ' // &
                           'permeability function class.'
        call printErrMsg(option)
    end select
  enddo

  
  ! for functions that are not phase-specific, check if PHASE was given:
  if (StringCompare('NONE',phase_keyword)) then
    phase_keyword = new_phase_keyword
    ! a liq or gas phase should now be specified for the non-phase-specific
    ! functions, so check if it was:
    if (StringCompare('NONE',new_phase_keyword)) then
      ! entering means the new phase keyword was also NONE (the default), so
      ! throw an error and abort:
      option%io_buffer = 'PHASE is not specified for ' // trim(error_string) 
      call printErrMsg(option)
    endif
  endif
  
  ! liquid phase relative permeability function check:
  if (StringCompare('LIQUID',phase_keyword)) then
    if (StringCompare('GAS',new_phase_keyword)) then
      ! user is requesting a liquid relative perm func for a gas phase:
      option%io_buffer = 'A liquid-phase relative permeability function &
                         &is being requested for the gas phase under ' &
                         // trim(error_string) // '.'
      call printErrMsg(option)
    endif
  endif
  
  ! gas phase relative permeability function check:
  if (StringCompare('GAS',phase_keyword)) then
    if (StringCompare('LIQUID',new_phase_keyword)) then
      ! user is requesting a gas relative perm func for a liquid phase:
      option%io_buffer = 'A gas-phase relative permeability function &
                         &is being requested for the liquid phase under ' &
                         // trim(error_string) // '.'
      call printErrMsg(option)
    endif
  endif

  if (smooth) then
    call permeability_function%SetupPolynomials(option,error_string)
  endif
  
end subroutine PermeabilityFunctionRead

! ************************************************************************** !

subroutine CharacteristicCurvesAddToList(new_characteristic_curves,list)
  ! 
  ! Adds a characteristic curves object to linked list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/02/07
  ! 

  implicit none
  
  class(characteristic_curves_type), pointer :: new_characteristic_curves
  class(characteristic_curves_type), pointer :: list

  class(characteristic_curves_type), pointer :: cur_characteristic_curves
  
  if (associated(list)) then
    cur_characteristic_curves => list
    ! loop to end of list
    do
      if (.not.associated(cur_characteristic_curves%next)) exit
      cur_characteristic_curves => cur_characteristic_curves%next
    enddo
    cur_characteristic_curves%next => new_characteristic_curves
  else
    list => new_characteristic_curves
  endif
  
end subroutine CharacteristicCurvesAddToList

! ************************************************************************** !

subroutine CharCurvesConvertListToArray(list,array,option)
  ! 
  ! Creates an array of pointers to the characteristic curves objects in the 
  ! list
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/11/07
  ! 

  use String_module
  use Option_module
  
  implicit none
  
  class(characteristic_curves_type), pointer :: list
  type(characteristic_curves_ptr_type), pointer :: array(:)
  type(option_type) :: option
    
  class(characteristic_curves_type), pointer :: cur_characteristic_curves
  PetscInt :: count

  count = 0
  cur_characteristic_curves => list
  do 
    if (.not.associated(cur_characteristic_curves)) exit
    count = count + 1
    cur_characteristic_curves => cur_characteristic_curves%next
  enddo
  
  if (associated(array)) deallocate(array)
  allocate(array(count))
  
  count = 0
  cur_characteristic_curves => list
  do 
    if (.not.associated(cur_characteristic_curves)) exit
    count = count + 1
    array(count)%ptr => cur_characteristic_curves
    if (cur_characteristic_curves%test .and. &
        option%myrank == option%io_rank) then
      call CharacteristicCurvesTest(cur_characteristic_curves,option)
    endif
    cur_characteristic_curves => cur_characteristic_curves%next
  enddo

end subroutine CharCurvesConvertListToArray

! ************************************************************************** !

function CharCurvesGetGetResidualSats(characteristic_curves,option)
  ! 
  ! Returns the residual saturations associated with a characteristic curves
  ! object
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  ! 

  use Option_module
  
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option

  PetscReal :: CharCurvesGetGetResidualSats(option%nphase)

  CharCurvesGetGetResidualSats(1) = &
    characteristic_curves%liq_rel_perm_function%Sr
  if (option%nphase > 1) then
    select type(rpf=>characteristic_curves%gas_rel_perm_function)
      class is(rpf_Mualem_VG_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_Mualem_VG_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_Burdine_BC_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_Burdine_BC_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_Mualem_BC_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_Mualem_BC_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_Burdine_VG_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_Burdine_VG_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_TOUGH2_IRP7_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_Mualem_Linear_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_Mualem_Linear_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_Burdine_Linear_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_Burdine_Linear_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_BRAGFLO_KRP1_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_BRAGFLO_KRP1_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_BRAGFLO_KRP5_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_BRAGFLO_KRP5_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_BRAGFLO_KRP9_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_BRAGFLO_KRP9_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_BRAGFLO_KRP4_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_BRAGFLO_KRP4_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_BRAGFLO_KRP11_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_BRAGFLO_KRP11_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_mK_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_mK_gas_type)
        CharCurvesGetGetResidualSats(2) = rpf%Srg
      class is(rpf_TOUGH2_Linear_oil_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sro
      class is(rpf_mod_BC_liq_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rpf_mod_BC_oil_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sro
      class is(rel_perm_func_constant_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class is(rel_perm_func_default_type)
        CharCurvesGetGetResidualSats(2) = rpf%Sr
      class default
        option%io_buffer = 'Relative permeability class not supported in ' // &
          'CharCurvesGetGetResidualSats.'
        call printErrMsg(option)
    end select
     
  endif

end function CharCurvesGetGetResidualSats

! ************************************************************************** !

function CharacteristicCurvesGetID(characteristic_curves_array, &
                                   characteristic_curves_name, &
                                   material_property_name, option)
  ! 
  ! Returns the ID of the characteristic curves object named
  ! "characteristic_curves_name"
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/11
  ! 

  use Option_module
  use String_module
  
  type(characteristic_curves_ptr_type), pointer :: &
    characteristic_curves_array(:)
  character(len=MAXWORDLENGTH) :: characteristic_curves_name
  character(len=MAXWORDLENGTH) :: material_property_name
  type(option_type) :: option

  PetscInt :: CharacteristicCurvesGetID

  CharacteristicCurvesGetID = 0
  do CharacteristicCurvesGetID = 1, size(characteristic_curves_array)
    if (StringCompare(characteristic_curves_name, &
                      characteristic_curves_array( &
                        CharacteristicCurvesGetID)%ptr%name)) then
      return
    endif
  enddo
  option%io_buffer = 'Characteristic curves "' // &
           trim(characteristic_curves_name) // &
           '" in material property "' // &
           trim(material_property_name) // &
           '" not found among available characteristic curves.'
  call printErrMsg(option)    

end function CharacteristicCurvesGetID

! ************************************************************************** !

subroutine CharacteristicCurvesTest(characteristic_curves,option)
  ! 
  ! Outputs values of characteristic curves over a range of values
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  use Option_module

  implicit none
  
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: phase

  call characteristic_curves%saturation_function%Test( &
                                                 characteristic_curves%name, &
                                                 option)
  phase = 'liquid'
  call characteristic_curves%liq_rel_perm_function%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
              
  if ( associated(characteristic_curves%gas_rel_perm_function) ) then
    phase = 'gas'
    call characteristic_curves%gas_rel_perm_function%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  endif

  if ( associated(characteristic_curves%oil_rel_perm_function) ) then
    phase = 'oil'
    call characteristic_curves%oil_rel_perm_function%Test( &
                                                 characteristic_curves%name, &
                                                 phase,option)
  end if
  
end subroutine CharacteristicCurvesTest

! ************************************************************************** !

subroutine CharacteristicCurvesVerify(characteristic_curves,option)
  ! 
  ! Checks if required parameters have been set for each curve type.
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/29/14
  !
  use Option_module

  implicit none
  
  class(characteristic_curves_type) :: characteristic_curves
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  string = 'CHARACTERISTIC_CURVES(' // trim(characteristic_curves%name) // &
           '),'

  call characteristic_curves%saturation_function%Verify(string,option)
  
  if (associated(characteristic_curves%liq_rel_perm_function) ) then
    call characteristic_curves%liq_rel_perm_function%Verify(string,option)
  else
    option%io_buffer = 'A liquid phase relative permeability curve has &
                       &not been set under CHARACTERISTIC_CURVES "' // &
                       trim(characteristic_curves%name) // '". A &
                       &PERMEABILITY_FUNCTION block must be specified &
                       &for the liquid phase.'
    call printErrMsg(option)
  end if

  if (associated(characteristic_curves%gas_rel_perm_function) ) then
    call characteristic_curves%gas_rel_perm_function%Verify(string,option)
  end if

  if ( associated(characteristic_curves%oil_rel_perm_function) ) then  
    call characteristic_curves%oil_rel_perm_function%Verify(string,option) 
  end if
  
end subroutine CharacteristicCurvesVerify

! ************************************************************************** !

! Begin Base Routines
function PolynomialCreate()

  implicit none
  
  type(polynomial_type), pointer :: PolynomialCreate  

  allocate(PolynomialCreate)
  PolynomialCreate%low = 0.d0
  PolynomialCreate%high = 0.d0
  PolynomialCreate%coefficients(:) = 0.d0
  
end function PolynomialCreate

! ************************************************************************** !

subroutine SFBaseInit(this)

  implicit none
  
  class(sat_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%sat_poly)
  nullify(this%pres_poly)
#ifdef SMOOTHING2
  nullify(this%sat_poly2)
  nullify(this%pres_poly2)
#endif
  this%Sr = UNINITIALIZED_DOUBLE
  this%pcmax = DEFAULT_PCMAX
  
end subroutine SFBaseInit

! ************************************************************************** !

subroutine SFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  
  
  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call printErrMsg(option)
  endif
  
end subroutine SFBaseVerify

! ************************************************************************** !

subroutine RPFBaseInit(this)

  implicit none
  
  class(rel_perm_func_base_type) :: this  

  ! Cannot allocate here.  Allocation takes place in daughter class
  nullify(this%poly)
#ifdef SMOOTHING2
  nullify(this%poly2)
#endif
  this%Sr = UNINITIALIZED_DOUBLE
  
end subroutine RPFBaseInit

! ************************************************************************** !

subroutine RPFBaseVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  if (Uninitialized(this%Sr)) then
    option%io_buffer = UninitializedMessage('LIQUID_RESIDUAL_SATURATION', &
                                            name)
    call printErrMsg(option)
  endif
  
end subroutine RPFBaseVerify

! ************************************************************************** !

subroutine SFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing saturation functions

  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)
  
end subroutine SFBaseSetupPolynomials

! ************************************************************************** !

subroutine RPFBaseSetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing relative permeability functions

  use Option_module
  
  implicit none
  
  class(rel_perm_func_base_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  option%io_buffer = 'Smoothing not supported for ' // trim(error_string)
  call printErrMsg(option)
  
end subroutine RPFBaseSetupPolynomials

! ************************************************************************** !

subroutine SFBaseCapillaryPressure(this,liquid_saturation, &
                                     capillary_pressure,option)
  use Option_module
  
  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseCapillaryPressure must be extended.'
  call printErrMsg(option)
  
end subroutine SFBaseCapillaryPressure

! ************************************************************************** !

subroutine SFBaseSaturation(this,capillary_pressure,liquid_saturation, &
                              dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFBaseSaturation must be extended.'
  call printErrMsg(option)
  
end subroutine SFBaseSaturation

! ************************************************************************** !

subroutine SFBaseTest(this,cc_name,option)

  use Option_module

  implicit none
  
  class(sat_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt, parameter :: num_values = 101
  PetscReal :: pc, pc_increment
  PetscReal :: perturbation
  PetscReal :: capillary_pressure(num_values)
  PetscReal :: capillary_pressure_pert(num_values)
  PetscReal :: liquid_saturation(num_values)
  PetscReal :: liquid_saturation_pert(num_values)
  PetscReal :: dsat_dpres(num_values)
  PetscReal :: dsat_dpres_numerical(num_values)
  PetscReal :: dummy_real(num_values)
  PetscInt :: count, i

  ! calculate saturation as a function of capillary pressure
  ! start at 1 Pa up to maximum capillary pressure
  pc = 1.d0
  pc_increment = 1.d0
  perturbation = 1.d-6
  count = 0
  do
    if (pc > this%pcmax) exit
    count = count + 1
    call this%Saturation(pc,liquid_saturation(count),dsat_dpres(count),option)
    capillary_pressure(count) = pc
    ! calculate numerical derivative dsat_dpres_numerical
    capillary_pressure_pert(count) = pc + pc*perturbation
    call this%Saturation(capillary_pressure_pert(count), &
                         liquid_saturation_pert(count),dummy_real(count),option)
    dsat_dpres_numerical(count) = (liquid_saturation_pert(count) - &
         & liquid_saturation(count))/(pc*perturbation)*(-1.d0) ! dPc/dPres
    ! get next value for pc
    if (pc > 0.99d0*pc_increment*10.d0) pc_increment = pc_increment*10.d0
    pc = pc + pc_increment
  enddo

  write(string,*) cc_name
  string = trim(cc_name) // '_pc_sat.dat'
  open(unit=86,file=string)
  write(86,*) '"capillary pressure", "saturation", "dsat/dpres", &
              &"dsat/dpres_numerical"'
  do i = 1, count
    write(86,'(4es14.6)') capillary_pressure(i), liquid_saturation(i), &
                          dsat_dpres(i), dsat_dpres_numerical(i)
  enddo
  close(86)

 ! calculate capillary pressure as a function of saturation
  do i = 1, num_values
    liquid_saturation(i) = dble(i-1)*0.01d0
    call this%CapillaryPressure(liquid_saturation(i),capillary_pressure(i), &
                                option)
  enddo
  count = num_values

  write(string,*) cc_name
  string = trim(cc_name) // '_sat_pc.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "capillary pressure"'
  do i = 1, count
    write(86,'(2es14.6)') liquid_saturation(i), capillary_pressure(i)
  enddo
  close(86)

end subroutine SFBaseTest

! ************************************************************************** !

subroutine RPF_Base_RelPerm(this,liquid_saturation,relative_permeability, &
                            dkr_sat,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'RPF_Base_RelPerm must be extended.'
  call printErrMsg(option)
  
end subroutine RPF_Base_RelPerm

! ************************************************************************** !

subroutine RPF_Base_Test(this,cc_name,phase,option)

  use Option_module

  implicit none
  
  class(rel_perm_func_base_type) :: this
  character(len=MAXWORDLENGTH) :: cc_name
  character(len=MAXWORDLENGTH) :: phase
  type(option_type), intent(inout) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: i
  PetscInt, parameter :: num_values = 101
  PetscReal :: perturbation
  PetscReal :: liquid_saturation(num_values)
  PetscReal :: liquid_saturation_pert(num_values)
  PetscReal :: kr(num_values)
  PetscReal :: kr_pert(num_values)
  PetscReal :: dkr_dsat(num_values)
  PetscReal :: dkr_dsat_numerical(num_values)
  PetscReal :: dummy_real(num_values)
  
  perturbation = 1.d-6

  do i = 1, num_values
    liquid_saturation(i) = dble(i-1)*0.01d0
    call this%RelativePermeability(liquid_saturation(i),kr(i),dkr_dsat(i), &
                                   option)
    ! calculate numerical derivative dkr_dsat_numerical
    liquid_saturation_pert(i) = liquid_saturation(i) &
                                + liquid_saturation(i)*perturbation
    call this%RelativePermeability(liquid_saturation_pert(i),kr_pert(i), &
                                   dummy_real(i),option)
    dkr_dsat_numerical(i) = (kr_pert(i) - kr(i))/ &
                            (liquid_saturation(i)*perturbation)
  enddo

  write(string,*) cc_name
  string = trim(cc_name) // '_' // trim(phase) // '_rel_perm.dat'
  open(unit=86,file=string)
  write(86,*) '"saturation", "' // trim(phase) // ' relative permeability", "' &
              // trim(phase) // ' dkr/dsat", "' // trim(phase) // &
              ' dkr/dsat_numerical"'
  do i = 1, size(liquid_saturation)
    write(86,'(4es14.6)') liquid_saturation(i), kr(i), dkr_dsat(i), &
                          dkr_dsat_numerical(i)
  enddo
  close(86)

end subroutine RPF_Base_Test
! End Base Routines

! ************************************************************************** !

! Begin SF: Default
function SF_Default_Create()

  ! Creates the default saturation function object

  implicit none
  
  class(sat_func_default_type), pointer :: SF_Default_Create
  
  allocate(SF_Default_Create)
  call SFBaseInit(SF_Default_Create)
  SF_Default_Create%Sr = 0.d0
  
end function SF_Default_Create

! ************************************************************************** !

subroutine SFDefaultVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_default_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  option%io_buffer = 'A default Saturation Function has been chosen in ' // &
    trim(name) // '.'
  call printWrnMsg(option)
  
end subroutine SFDefaultVerify

! ************************************************************************** !

subroutine SFDefaultCapillaryPressure(this,liquid_saturation, &
                                      capillary_pressure,option)
  use Option_module
  
  implicit none
  
  class(sat_func_default_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation < 1.d0) then
    option%io_buffer = 'SFDefaultCapillaryPressure is a dummy routine used &
      &for saturated flow only.  The user must specify a valid &
      &SATURATION_FUNCTION.'
    call printErrMsgByRank(option)
  endif

end subroutine SFDefaultCapillaryPressure

! ************************************************************************** !

subroutine SFDefaultSaturation(this,capillary_pressure,liquid_saturation, &
                               dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_default_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  option%io_buffer = 'SFDefaultSaturation is a dummy routine used &
    &for saturated flow only.  The user must specify a valid &
    &SATURATION_FUNCTION.'
  call printErrMsgByRank(option)

end subroutine SFDefaultSaturation

! ************************************************************************** !

function RPF_Default_Create()

  ! Creates the default relative permeability function object

  implicit none
  
  class(rel_perm_func_default_type), pointer :: RPF_Default_Create
  
  allocate(RPF_Default_Create)
  call RPFBaseInit(RPF_Default_Create)
  RPF_Default_Create%Sr = 0.d0
  
end function RPF_Default_Create

! ************************************************************************** !

subroutine RPFDefaultVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rel_perm_func_default_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  option%io_buffer = 'A default Relative Permeability Function has been ' // &
    'chosen in ' // trim(name) // '.'
  call printWrnMsg(option)

end subroutine RPFDefaultVerify

! ************************************************************************** !

subroutine RPF_DefaultRelPerm(this,liquid_saturation,relative_permeability, &
                            dkr_sat,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_default_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation < 1.d0) then
    option%io_buffer = 'RPF_Default_RelPerm is a dummy routine used &
      &for saturated flow only.  The user must specify a valid &
      &PERMEABILITY_FUNCTION.'
    call printErrMsgByRank(option)
  endif
  relative_permeability = 1.d0
  
end subroutine RPF_DefaultRelPerm
! End Default Routines

! ************************************************************************** !

subroutine SF_Ice_CapillaryPressure(this, pres_l, tc, &
                                   ice_pc, dice_pc_dpres, dice_pc_dt, option)
  !
  ! Computes the ice capillary_pressure as a function of Pres_l, Tc
  ! Mainly from Painter et al. (2011), Painter and Karra (2014)
  !
  ! Author: Fengming Yuan
  !         Based on relevant saturation_functions by Satish K.
  !         Revised with freezing-thawing zone smoothing following ATS algorithm
  ! Date: 06/01/2016
  !
  use Option_module
  use EOS_Water_module
  use Utility_module

  implicit none

  class(sat_func_base_type) :: this
  PetscReal, intent(in) :: pres_l      ! liquid water pressure head (atm. P adjusted), in Pa
  PetscReal, intent(in) :: tc          ! in oC
  PetscReal, intent(out) :: ice_pc, dice_pc_dt, dice_pc_dpres     ! in Pa
  type(option_type) :: option

  PetscReal :: pcgl, pw, tk

  PetscReal, parameter :: beta = 2.33d0           ! dimensionless -- ratio of soil ice surf. tension
  PetscReal, parameter :: T0   = 273.15d0         ! freezing-point at standard pressure: in K
  PetscReal, parameter :: Lf   = HEAT_OF_FUSION   ! fusion heat (in J/kg)
  PetscReal :: gamma, alpha, dalpha_drhol
  PetscReal :: rhol, drhol_dp, drhol_dt
  PetscReal :: rhoi

  PetscReal :: Tf, dTf_dt, dTf_dp
  PetscReal :: tftheta, dtftheta_dt, dtftheta_dp

  PetscReal :: Hfunc, dHfunc, tempreal
  PetscReal :: deltaTf, xTf, a, b, c
  PetscReal :: beta2, T_star, T_star_th, T_star_min, T_star_max
  PetscReal :: dice_pc_dp
  PetscReal, parameter :: dpc_dpres = -1.d0
  PetscReal :: PCGL_MAX_FRZ = 1.d8   ! max. soil matrical potential (Pa) for liq. water exists under super-cooling

  PetscErrorCode :: ierr

  !---------------------
  !
  pcgl = max(0.d0, option%reference_pressure - pres_l)       ! always non-negative (0 = saturated)
  if (pcgl > abs(this%pcmax)) pcgl = this%pcmax

  PCGL_MAX_FRZ = this%pcmax+1.0d0


  Tk = tc + T0

  ! if ice module turns on, 2-phase saturation recalculated (liq. and ice) under common 'pw' and 'tc'
  pw = max(option%reference_pressure, pres_l)

  ! --------------------

  ! constant 'rhol' (liq. water)
  rhol     = 999.8d0            ! kg/m3: kmol/m3*kg/kmol
  drhol_dp = 0.d0
  drhol_dt = 0.d0

  ! constant 'rhoi' (for ice)
  rhoi    = 916.7d0             ! kg/m3 at 273.15K

  if (option%use_th_freezing) then

    gamma       = beta*Lf
    alpha       = gamma/T0*rhol
    dalpha_drhol= gamma/T0

    Tf     = T0 - 1.d0/alpha*min(PCGL_MAX_FRZ,pcgl)                               ! P.-K. Eq.(10), omiga=1/beta
    dTf_dt = pcgl/alpha/alpha*(dalpha_drhol*drhol_dt)
    dTf_dp = (pcgl*dalpha_drhol*drhol_dp - alpha)/alpha/alpha   ! dpcgl_dp = 1.0
    if(pcgl>PCGL_MAX_FRZ) then
      dTf_dt = 0.d0
      dTf_dp = 0.d0
    endif

    xTf = Tk - T0

    select case (option%ice_model)

      case (PAINTER_EXPLICIT)

        ! explicit model from Painter (Comp. Geosci, 2011)
        ice_pc = pcgl
        dice_pc_dt = 0.d0
        dice_pc_dp = 1.d0

        if(tc<0.d0) then

          ice_pc = rhoi*beta*Lf*(-tc)/T0
          dice_pc_dt = -rhoi*beta*Lf/T0
          dice_pc_dp = 0.d0

        endif

      case (PAINTER_KARRA_EXPLICIT, PAINTER_KARRA_EXPLICIT_SMOOTH)

        ! The following is a slightly-modified version from PKE in saturation_function module
        ! without smoothing of freezing-thawing zone

        ! explicit model from Painter & Karra, VJZ (2014)
        tftheta = xTf/T0                              ! P.-K. Eq.(18): theta: (Tk-T0)/T0, assuming Tf~T0 (ignored FP depression) in Eq. (12)
        dtftheta_dt = 1.0d0/T0
        dtftheta_dp = 0.d0

        ice_pc     = -gamma * rhol*tftheta           ! P.-K. Eq.(18), first term (i.e. ice only)
        tempreal   = rhol*dtftheta_dt+tftheta*drhol_dt
        dice_pc_dt = -gamma * tempreal
        tempreal   = rhol*dtftheta_dp+tftheta*drhol_dp
        dice_pc_dp = -gamma * tempreal

        ! Heaviside function to truncate Eq. (18) at 'Tf': Tf = T0-1/alpha*pcgl, i.e. -alpha*(Tf-T0)=pcgl, the left side IS exactly ice_pc @Tf
        Hfunc = sign(0.5d0, -(Tk-Tf))+0.5d0
        dice_pc_dt =  dice_pc_dt*Hfunc                          ! no need to adjust dice_pc_dt due to dpcgl_dt = 0
        dice_pc_dp =  dice_pc_dp*Hfunc + (1.d0-Hfunc)           ! dpcgl_dp = 1.0
        ice_pc     =  ice_pc*Hfunc + pcgl*(1.d0-Hfunc)

        ! -----------
        ! smoothing 'ice_pc' when Tk ranging within deltaTf of T0, from PKE's PCice to 0.0
        ! from ATS, authored by Scott Painter et al.
        deltaTf = 1.0d-50              ! half-width of smoothing zone (by default, nearly NO smoothing)
        if(option%frzthw_halfwidth /= UNINITIALIZED_DOUBLE) deltaTf = option%frzthw_halfwidth
        if (option%ice_model==PAINTER_KARRA_EXPLICIT_SMOOTH .and. deltaTf>1.0d-50) then

#if 0
          ! F.-M. Yuan (2017-03-14): OFF
          ! symmetrically smoothing over 'T0', so may create a gap from Tf ~ T0-deltaTf (Tf could be far away from T0)
          if (abs(xTf)<deltaTf) then
            a = deltaTf/4.0d0
            b = -0.5d0
            c = 0.25d0/deltaTf

            tempreal   = a+b*xTf+c*xTf*xTf
            ice_pc     = alpha*tempreal
            dice_pc_dt = alpha*(b+2.0d0*c*xTf) + &
                       dalpha_drhol*drhol_dt*tempreal        ! in case we need this later

            ! a note here:
            ! (1) if xTf=-deltaTf (negative over Tf0), tempreal = deltaTf/4.0 + 0.5*deltaTf + 0.25*deltaTf = 1.0 * deltaTf
            !                   so, ice_pc = alpha * deltaTf = -alpha * xTf;
            !
            !                     and, b+2*c*xTf = -1, so ice_pc_dt = -alpha
            ! (2) if xTf=deltaTf (positive over Tf0), tempreal = deltaTf/4.0 - 0.5*deltaTf + 0.25*deltaTf = 0
            !                   so, ice_pc = 0.
            !                     and, b+2*x*xTf = 0, ice_pc_dt = 0
            ! Then it means the smoothing_curve ends exactly with the original format
          endif

#else

          ! using the new Heaveside Smoothing function:
          ! ice_pc@Tf-deltaTf --> ice_pc = pcgl @T0+deltaTf  (note: @Tf, ice_pc = pcgl by PKE, so this smoothing actually span over a large ranges around Tf~T0)
          ! (note: Tf could be far way from T0, and 'ice_pc' above trucated at 'Tf'; So this smoothing will span over Tf-deltaTf ~ T0+deltaTf asymmetrically)
          call HFunctionSmooth(Tk, Tf-deltaTf, T0+deltaTf, Hfunc, dHfunc)
          dice_pc_dt = dice_pc_dt * Hfunc + (ice_pc - pcgl) * dHfunc
          dice_pc_dp = (dice_pc_dp - 1.0d0) * Hfunc + 1.0d0          ! dHfunc_dp = 0, dpcgl_dp = 1
          ice_pc = (ice_pc - pcgl) * Hfunc + pcgl                    ! do this after derivatives
#endif
        endif !(option%ice_model==PAINTER_KARRA_EXPLICIT_SMOOTH)
        ! -----------

      case (DALL_AMICO)

        ! Model from Dall'Amico (2010) and Dall' Amico et al. (2011)
        ! rewritten following 'saturation_function.F90:SatFuncComputeIceDallAmico()'
        ! NOTE: here only calculate Pc1 and its derivatives
        ice_pc = pcgl
        dice_pc_dt = 0.d0
        dice_pc_dp = 1.d0

        !
        T_star_th  = 5.d-1                       ! unit: Kevin
        beta2      = beta !1.d0                  ! NOT sure why changed to 1. in saturation_function.F90.

        T_star = T0-1.d0/beta2/Lf/rhol*pcgl      ! 'T_star' shall be Tf same as other ice_model, when using same CONSTANTS
        T_star_min = T_star - T_star_th
        T_star_max = T_star

        ! H function
        if (Tk<T_star_min) then
          Hfunc = 1.0d0
          dHfunc= 0.d0
        else if (Tk>T_star_max) then
          Hfunc = 0.0d0
          dHfunc= 0.d0
        else
          tempreal = (Tk-T_star_min)/(T_star_max-T_star_min)
          tempreal = 1.d0 - tempreal*tempreal
          Hfunc = tempreal * tempreal
          dHfunc= -4.0d0*tempreal &
                  *(Tk-T_star_min)/(T_star_max-T_star_min)/(T_star_max-T_star_min)
        endif

        !
        ice_pc     = pcgl - beta2 *  (Tk-T_star)/T_star *Lf*rhol * Hfunc  ! L2030 in saturation_function.F90

        ! saturation_function.F90: L2053 -
        ! dsl_dT = dS1 * (-beta*L_f*rhol_l*H/T_star-beta*theta*L_f*rhol_l*dH_dT)
        ! i.e., dS1* [(-beta*L_f*rhol_l)*(H/T_star+theta*dH_dT)], with theta=(T-T_star)/T_star
        !       for second term:
        !           -beta/T_star*L_f*rhol_l*(H+(T-T_star)*dH_dT)
        !               (because dT_star/dt=0)
        dice_pc_dt = -beta2/T_star*Lf*rhol*        &
                     (Hfunc + (Tk-T_star)*dHfunc)

        ! saturation_function.F90: L2050 -
        !  dsl_dpl = -dS1*(1-T*T_0/T_star/T_star*H), in which
        !   (1) '-dS1' is w.r.t. '-Pc1' actually, because in L2045 it was negatived once
        !   (2) (1-T*T_0/T_star/T_star*H) is 'dp_fh2o_dP' in L.2034.
        !dice_pc_dp = -(1.0d0 - Tk*T0/T_star/T_star*Hfunc)   ! '-1' is for reverse '_dP' to '_dpc'

        ! HOWever, if we do derivative here for 'ice_pc' above, 'ice_pc' shall be
        !   ice_pc = pcgl - beta2*Lf*rhol*Hfunc* (Tk/T_star-1)
        !  w.r.t. 'pc', the second term can be simplified as: -beta2*Lf*rhol*Hfunc*Tk * (1/T_star)
        !              in which, dT_star_dp = -1.d0/beta2/Lf/rhol, and
        !                        d(1/T_star)_dp = dT_star_dp/T_star/T_star
        ! i.e.,
        dice_pc_dp = -(1.0d0 + Tk/T_star/T_star*Hfunc)   ! '-1' is for reverse '_dP' to '_dpc'

      case default
        option%io_buffer = 'SF_Ice_CapillaryPressure: characteristic-curve now only support ice-model: ' // &
          'PAINTER_KARRA_EXPLICIT, or PAINTER_KARRA_EXPLICIT_SMOOTH, or ' // &
          'PAINTER_EXPLICIT, or ' // &
          'DALL_AMICO '

        call printErrMsg(option)

    end select ! select case (option%ice_model)

    dice_pc_dpres = dice_pc_dp*dpc_dpres  ! convert to w.r.t 'pressure' from 'pc'
  else

    option%io_buffer = 'SF_Ice_CapillaryPressure: Ice model is OFF.'
    call printMsg(option)

  endif ! 'option%use_th_freezing'

end subroutine SF_Ice_CapillaryPressure

! ************************************************************************** !

! Begin SF: Constant
function SF_Constant_Create()

  ! Creates the default saturation function object

  implicit none
  
  class(sat_func_constant_type), pointer :: SF_Constant_Create
  
  allocate(SF_Constant_Create)
  call SFBaseInit(SF_Constant_Create)
  ! set Sr to zero as it doesn't matter, but must be initialized
  SF_Constant_Create%Sr = 0.d0 
  SF_Constant_Create%constant_capillary_pressure = UNINITIALIZED_DOUBLE
  SF_Constant_Create%constant_saturation = UNINITIALIZED_DOUBLE
  
end function SF_Constant_Create

! ************************************************************************** !

subroutine SFConstantVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_constant_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  character(len=MAXSTRINGLENGTH) :: string  

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,CONSTANT'
  endif
  call SFBaseVerify(this,string,option)
  select case(option%iflowmode)
    case(RICHARDS_MODE,TH_MODE)
      if (Initialized(this%constant_capillary_pressure)) then
        option%io_buffer = 'CONSTANT_CAPILLARY_PRESSURE is not supported for &
          &Richards or TH flow modes as CONSTANT_SATURATION must be applied. &
          &See ' // trim(string) // '.'
        call printErrMsg(option)
      endif
      if (Uninitialized(this%constant_saturation)) then
        option%io_buffer = 'CONSTANT_SATURATION must be specified for ' // &
          trim(string) // '.'
        call printErrMsg(option)
      endif
    case default
  end select

end subroutine SFConstantVerify

! ************************************************************************** !

subroutine SFConstantCapillaryPressure(this,liquid_saturation, &
                                      capillary_pressure,option)
  use Option_module
  
  implicit none
  
  class(sat_func_constant_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  capillary_pressure = this%constant_capillary_pressure

end subroutine SFConstantCapillaryPressure

! ************************************************************************** !

subroutine SFConstantSaturation(this,capillary_pressure,liquid_saturation, &
                               dsat_dpres,option)
  use Option_module

  implicit none
  
  class(sat_func_constant_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option

  liquid_saturation = this%constant_saturation
  dsat_dpres = 0.d0

end subroutine SFConstantSaturation
! End Constant Routines

! ************************************************************************** !

! Begin SF: van Genuchten
function SF_VG_Create()

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_VG_type), pointer :: SF_VG_Create
  
  allocate(SF_VG_Create)
  call SF_VG_Create%Init()
  
end function SF_VG_Create

! ************************************************************************** !

subroutine SF_VG_Init(this)

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_VG_type) :: this

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE
  
end subroutine SF_VG_Init

! ************************************************************************** !

subroutine SF_VG_Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_VG_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,VAN_GENUCHTEN'
  endif
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif   
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif   

end subroutine SF_VG_Verify

! ************************************************************************** !

subroutine SF_VG_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  
  implicit none
  
  class(sat_func_VG_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: n
  PetscReal :: Se
  PetscReal :: one_plus_pc_alpha_n
  PetscReal :: pc_alpha_n
  PetscReal :: pc_alpha
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  n = 1.d0/(1.d0-this%m)
  Se = (liquid_saturation-this%Sr)/(1.d0-this%Sr)
  one_plus_pc_alpha_n = Se**(-1.d0/this%m)
  pc_alpha_n = one_plus_pc_alpha_n - 1.d0
  pc_alpha = pc_alpha_n**(1.d0/n)
  capillary_pressure = pc_alpha/this%alpha
#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
  endif
#endif

  capillary_pressure = min(capillary_pressure,this%pcmax)

  ! notes by fmy @Mar-30-2016: this function appears NOT called anywhere by this moment (????)
  !        @May-05-2016: it's used when doing regression-tests of 'toil_ims' and 'general/liquid_gas.in'

end subroutine SF_VG_CapillaryPressure

! ************************************************************************** !

subroutine SF_VG_Saturation(this,capillary_pressure,liquid_saturation, &
                            dsat_dpres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_VG_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal, parameter :: pc_alpha_n_epsilon = 1.d-15
  PetscReal :: n
  PetscReal :: pc_alpha
  PetscReal :: pc_alpha_n
  PetscReal :: one_plus_pc_alpha_n
  PetscReal :: Se
  PetscReal :: dSe_dpc
  PetscReal, parameter :: dpc_dpres = -1.d0
  
  dsat_dpres = 0.d0
  
  if (associated(this%pres_poly)) then
    if (capillary_pressure < this%pres_poly%low) then
      liquid_saturation = 1.d0
      return
    else if (capillary_pressure < this%pres_poly%high) then
      call CubicPolynomialEvaluate(this%pres_poly%coefficients, &
                                   capillary_pressure,Se,dSe_dpc)
      liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
      dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
      return
    endif
  endif

  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    n = 1.d0/(1.d0-this%m)
    pc_alpha = capillary_pressure*this%alpha
    pc_alpha_n = pc_alpha**n
    !geh:  This conditional does not catch potential cancelation in 
    !      the dkr_sat deriviative calculation.  Therefore, I am setting
    !      an epsilon here
    !   if (1.d0 + pc_alpha_n == 1.d0) then ! check for zero perturbation
    if (pc_alpha_n < pc_alpha_n_epsilon) then 
      liquid_saturation = 1.d0
      !switch_to_saturated = PETSC_TRUE
      return
    endif
    one_plus_pc_alpha_n = 1.d0+pc_alpha_n
    Se = one_plus_pc_alpha_n**(-this%m)
    dSe_dpc = -this%m*n*this%alpha*pc_alpha_n/ &
            (pc_alpha*one_plus_pc_alpha_n**(this%m+1.d0))
    liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
    dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
  endif
  
end subroutine SF_VG_Saturation
! End SF: van Genuchten

! ************************************************************************** !

! Begin SF: Brooks-Corey
function SF_BC_Create()

  ! Creates the Brooks Corey capillary pressure function object

  implicit none
  
  class(sat_func_BC_type), pointer :: SF_BC_Create
  
  allocate(SF_BC_Create)
  call SF_BC_Create%Init()
  
end function SF_BC_Create

! ************************************************************************** !

subroutine SF_BC_Init(this)

  use Option_module

  implicit none
  
  class(sat_func_BC_type) :: this
  type(option_type) :: option

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%lambda = UNINITIALIZED_DOUBLE
  
end subroutine SF_BC_Init

! ************************************************************************** !

subroutine SF_BC_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(sat_func_BC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BROOKS_COREY'
  endif  
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif 
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif 
  
end subroutine SF_BC_Verify

! ************************************************************************** !

subroutine SF_BC_SetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing Brooks-Corey saturation function

  use Option_module
  use Utility_module
  
  implicit none
  
  class(sat_func_BC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  PetscReal :: b(4)

#ifdef SMOOTHING2
  ! F.-M. Yuan (2017-03-10): newly smoothing approach using Heaveside Function Smoothed
  ! a note here: when setup, only needs starting/ending points for smoothing

  ! WET-end
  if (associated(this%pres_poly)) call PolynomialDestroy(this%pres_poly)
  this%pres_poly => PolynomialCreate()
  this%pres_poly%low  = 0.000d0/this%alpha  ! have to make Se=1, exactly@pc=0 (otherwise, by BC equation, this is 1/alpha)
  this%pres_poly%high = 1.005d0/this%alpha  ! this is very sensitive (better close to 1/alpha)

  ! DRY-end
  if (associated(this%pres_poly2)) call PolynomialDestroy(this%pres_poly2)
  this%pres_poly2 => PolynomialCreate()
  this%pres_poly2%high= min(this%pcmax, &
          (0.001d0**(-1.d0/this%lambda))/this%alpha)       ! Se = 0.001 (note that if Se=0, pc=0)
  this%pres_poly2%low = min(this%pres_poly2%high-1.0d0, &  ! this is the starting point for smoothing, so it must be less than %high end.
          (0.005d0**(-1.d0/this%lambda))/this%alpha)       ! Se = 0.005


#else
  ! polynomial fitting pc as a function of saturation
  ! 1.05 is essentially pc*alpha (i.e. pc = 1.05/alpha)
  if (associated(this%sat_poly)) call PolynomialDestroy(this%sat_poly)
  this%sat_poly => PolynomialCreate()
  this%sat_poly%low = 1.05d0**(-this%lambda)
  this%sat_poly%high = 1.d0
  
  b = 0.d0
  ! fill right hand side
  ! capillary pressure at 1
  b(1) = 1.05d0/this%alpha 
  ! capillary pressure at 2
  b(2) = 0.d0
  ! derivative of pressure at saturation_1
  ! pc = Se**(-1/lambda)/alpha
  ! dpc_dSe = -1/lambda*Se**(-1/lambda-1)/alpha
  b(3) = -1.d0/this%lambda* &
          this%sat_poly%low**(-1.d0/this%lambda-1.d0)/ &
          this%alpha

  call QuadraticPolynomialSetup(this%sat_poly%low,this%sat_poly%high,b(1:3), &
                                ! indicates derivative given at 1
                                PETSC_TRUE) 
      
  this%sat_poly%coefficients(1:3) = b(1:3)

  ! polynomial fitting saturation as a function of pc
  !geh: cannot invert the pressure/saturation relationship above
  !     since it can result in saturations > 1 with both
  !     quadratic and cubic polynomials
  ! fill matix with values
  if (associated(this%pres_poly)) call PolynomialDestroy(this%pres_poly)
  this%pres_poly => PolynomialCreate()
  this%pres_poly%low = 0.95/this%alpha
  this%pres_poly%high = 1.05/this%alpha
  
  b = 0.d0
  ! Se at 1
  b(1) = 1.d0
  ! Se at 2
  b(2) = (this%pres_poly%high*this%alpha)** &
          (-this%lambda)
  ! derivative of Se at 1
  b(3) = 0.d0 
  ! derivative of Se at 2
  b(4) = -this%lambda/this%pres_poly%high* &
            (this%pres_poly%high*this%alpha)** &
              (-this%lambda)

  call CubicPolynomialSetup(this%pres_poly%low,this%pres_poly%high,b)

  this%pres_poly%coefficients(1:4) = b(1:4)

#endif
end subroutine SF_BC_SetupPolynomials

! ************************************************************************** !

subroutine SF_BC_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation using the
  ! Brooks-Corey formulation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  !
  use Option_module
  use Utility_module
  
  implicit none
  
  class(sat_func_BC_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dummy_real
#ifdef SMOOTHING2
  PetscReal :: Hfunc, dHfunc, x1, x0, y
#endif
  
  if (liquid_saturation <= this%Sr) then   ! F.-M. Yuan (2017-03-14): This IS not mathmatically correct, although this function seems NOT called except for unit test
#ifndef SMOOTHING2
    capillary_pressure = this%pcmax
    return
#endif
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  Se = (liquid_saturation-this%Sr)/(1.d0-this%Sr)
  if (associated(this%sat_poly)) then
    if (Se > this%sat_poly%low) then
      call QuadraticPolynomialEvaluate(this%sat_poly%coefficients(1:3), &
                                       Se,capillary_pressure,dummy_real)
      return
    endif
  endif
  capillary_pressure = (Se**(-1.d0/this%lambda))/this%alpha
#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
  endif
#endif  

  capillary_pressure = min(capillary_pressure,this%pcmax)

#ifdef SMOOTHING2
  ! fmy: by BC function, mathemaatically @Sr, pc not necessarily equals to pcmax
  !      So, @pcmax, S is not necessarily equaled to Sr, and requiring smoothing
  x1 = this%pcmax/10.0d0
  y  = (x1*this%alpha)**(-this%lambda)  ! Se
  x1 = this%Sr+(1.d0-this%Sr)*y         ! S to start smoothing

  x0 = this%pcmax
  y  = (x0*this%alpha)**(-this%lambda)
  x0 = this%Sr+(1.d0-this%Sr)*y         ! S to end smoothing

  call HFunctionSmooth(liquid_saturation, x1, x0, Hfunc, dHfunc)

  capillary_pressure = capillary_pressure*Hfunc + this%pcmax * (1.d0-Hfunc)

#endif
  
end subroutine SF_BC_CapillaryPressure

! ************************************************************************** !

subroutine SF_BC_Saturation(this,capillary_pressure,liquid_saturation, &
                            dsat_dpres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
    
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_BC_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: Se
  PetscReal :: dSe_dpc
  PetscReal, parameter :: dpc_dpres = -1.d0
#ifdef SMOOTHING2
  PetscReal :: Hfunc, dHfunc
#endif

  
  dsat_dpres = 0.d0
  
#ifndef SMOOTHING2
  ! reference #1
  if (associated(this%pres_poly)) then
    if (capillary_pressure < this%pres_poly%low) then
      liquid_saturation = 1.d0
      return
    else if (capillary_pressure < this%pres_poly%high) then
      call CubicPolynomialEvaluate(this%pres_poly%coefficients, &
                                   capillary_pressure,Se,dSe_dpc)
      liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
      dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
      return
    endif
  else
    if (capillary_pressure < 1.d0/this%alpha) then
      liquid_saturation = 1.d0
      dsat_dpres = 0.d0
      return
    endif
  endif

  pc_alpha_neg_lambda = (capillary_pressure*this%alpha)**(-this%lambda)
  Se = pc_alpha_neg_lambda
  dSe_dpc = -this%lambda/capillary_pressure*pc_alpha_neg_lambda
  liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
  dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
  
#else

  ! F.-M. Yuan (2017-03-10): needs full range function
  if (capillary_pressure < 1.d0/this%alpha) then
    liquid_saturation = 1.d0
    dsat_dpres = 0.d0
  else
    pc_alpha_neg_lambda = (capillary_pressure*this%alpha)**(-this%lambda)
    Se = pc_alpha_neg_lambda
    dSe_dpc = -this%lambda/capillary_pressure*pc_alpha_neg_lambda
    liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
    dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
  endif

  ! F.-M. Yuan (2017-03-10): new smoothing approach using Heaveside Function Smoothed
  ! WET-end
  if (associated(this%pres_poly)) then
    if (capillary_pressure<=this%pres_poly%high) then
      call HFunctionSmooth(capillary_pressure, this%pres_poly%high, this%pres_poly%low, &
                           Hfunc, dHfunc)
      !sat = 1.0 + (sat-1.0) * Hfunc, in which '1.0' is the base sat @pc=low end
      dsat_dpres = dsat_dpres*Hfunc + (liquid_saturation-1.d0) * dHfunc
      liquid_saturation = 1.0d0 + (liquid_saturation-1.d0)*Hfunc
    endif
  endif

  ! DRY-end
  if (associated(this%pres_poly2)) then
    if (capillary_pressure>=this%pres_poly2%low) then
      call HFunctionSmooth(capillary_pressure, this%pres_poly2%low, this%pres_poly2%high, &
                           Hfunc, dHfunc)
      !sat = sr + (sat-sr) * Hfunc, in which 'sr' is the base sat @pc=high (high end)
      dsat_dpres = dsat_dpres * Hfunc + (liquid_saturation-this%Sr) * dHfunc
      liquid_saturation = this%Sr+(liquid_saturation-this%Sr) * Hfunc
    endif
  endif

#endif

end subroutine SF_BC_Saturation
! End SF: Brooks-Corey

! ************************************************************************** !

! Begin SF: Linear Model
function SF_Linear_Create()

  ! Creates the Linear capillary pressure function object

  implicit none
  
  class(sat_func_Linear_type), pointer :: SF_Linear_Create
  
  allocate(SF_Linear_Create)
  call SF_Linear_Create%Init()
  
end function SF_Linear_Create

! ************************************************************************** !

subroutine SF_Linear_Init(this)

  ! Creates the Linear capillary pressure function object

  implicit none
  
  class(sat_func_Linear_type) :: this

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  
end subroutine SF_Linear_Init

! ************************************************************************** !

subroutine SF_Linear_Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_Linear_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,LINEAR'
  endif
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif   

end subroutine SF_Linear_Verify

! ************************************************************************** !

subroutine SF_Linear_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation
  ! 
  !   

  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  !
  use Option_module
  
  implicit none
  
  class(sat_func_Linear_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  Se = (liquid_saturation-this%Sr)/(1.d0-this%Sr)
  capillary_pressure = (1.d0/this%alpha-this%pcmax)*Se + this%pcmax
#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
  endif
#endif  

  capillary_pressure = min(capillary_pressure,this%pcmax)
  
end subroutine SF_Linear_CapillaryPressure

! ************************************************************************** !

subroutine SF_Linear_Saturation(this,capillary_pressure,liquid_saturation, &
                            dsat_dpres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  ! 
  !   
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_Linear_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dSe_dpc
  PetscReal, parameter :: dpc_dpres = -1.d0
  
  dsat_dpres = 0.d0

  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    Se = (this%pcmax-capillary_pressure) / (this%pcmax-1.d0/this%alpha)
    dSe_dpc = -1.d0/(this%pcmax-1.d0/this%alpha)
    liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
    dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
  endif 

end subroutine SF_Linear_Saturation
! End SF: Linear Model

! ************************************************************************** !

! Begin SF: BRAGFLO KRP1 Model
function SF_BF_KRP1_Create()

  ! Creates the BRAGFLO KRP1 capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP1_type), pointer :: SF_BF_KRP1_Create
  
  allocate(SF_BF_KRP1_Create)
  call SF_BF_KRP1_Create%Init()
  
end function SF_BF_KRP1_Create

! ************************************************************************** !

subroutine SF_BF_KRP1_Init(this)

  ! Creates the BRAGFLO KRP1 capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP1_type) :: this

  call SFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%alpha = UNINITIALIZED_DOUBLE
  this%m = UNINITIALIZED_DOUBLE
  
end subroutine SF_BF_KRP1_Init

! ************************************************************************** !

subroutine SF_BF_KRP1_Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_BF_KRP1_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,VAN_GENUCHTEN'
  endif
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif   
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif   
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('Srg',string)
    call printErrMsg(option)
  endif   

end subroutine SF_BF_KRP1_Verify

! ************************************************************************** !

subroutine SF_BF_KRP1_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation using the
  ! van Genuchten formulation with non-zero gas residual saturation
  ! BRAGFLO UM 6.02 pg 41, 42; eq.120 - eq.124
  ! Modified according to KRP=1 option of BRAGFLO
  ! Explanation: residual gas saturation is in the denominator of effective 
  ! saturation
  !
  ! Author: Heeho Park
  ! Date: 11/17/16
  !
  use Option_module
  
  implicit none
  
  class(sat_func_BF_KRP1_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se2
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  Se2 = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  capillary_pressure = (1/this%alpha)*(Se2**(-1/this%m)-1)**(1-this%m)
#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
  endif
#endif
  capillary_pressure = min(capillary_pressure,this%pcmax)
  
end subroutine SF_BF_KRP1_CapillaryPressure

! ************************************************************************** !

subroutine SF_BF_KRP1_Saturation(this,capillary_pressure,liquid_saturation, &
                            dsat_dpres,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation using the
  ! van Genuchten formulation with non-zero gas residual saturation
  ! BRAGFLO UM 6.02 pg 41, 42; eq.120 - eq.124
  ! Modified according to KRP=1 option of BRAGFLO
  ! Explanation: residual gas saturation is in the denominator of effective 
  ! saturation
  !
  ! Author: Heeho Park
  ! Date: 11/17/16
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_BF_KRP1_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: n
  PetscReal :: pc_alpha
  PetscReal :: pc_alpha_n
  PetscReal :: one_plus_pc_alpha_n
  PetscReal :: Se
  PetscReal :: dSe_dpc
  
  dsat_dpres = 0.d0
  
  if (associated(this%pres_poly)) then
    if (capillary_pressure < this%pres_poly%low) then
      liquid_saturation = 1.d0
      return
    else if (capillary_pressure < this%pres_poly%high) then
      call CubicPolynomialEvaluate(this%pres_poly%coefficients, &
                                   capillary_pressure,Se,dSe_dpc)
      liquid_saturation = this%Sr + (1.d0-this%Sr)*Se
      dsat_dpres = -(1.d0-this%Sr)*dSe_dpc
      return
    endif
  endif

  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    n = 1.d0/(1.d0-this%m)
    pc_alpha = capillary_pressure*this%alpha
    pc_alpha_n = pc_alpha**n
    one_plus_pc_alpha_n = 1.d0+pc_alpha_n
    Se = one_plus_pc_alpha_n**(-this%m)
    dSe_dpc = -this%m*n*this%alpha*pc_alpha_n/ &
            (pc_alpha*one_plus_pc_alpha_n**(this%m+1.d0))
    liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se
    dsat_dpres = -(1.d0-this%Sr-this%Srg)*dSe_dpc
  endif
  
end subroutine SF_BF_KRP1_Saturation
! End SF: BRAGFLO KRP1 Model

! ************************************************************************** !

! Begin SF: BRAGFLO KRP5 Model
function SF_BF_KRP5_Create()

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP5_type), pointer :: SF_BF_KRP5_Create
  
  allocate(SF_BF_KRP5_Create)
  call SF_BF_KRP5_Create%Init()
  
end function SF_BF_KRP5_Create

! ************************************************************************** !

subroutine SF_BF_KRP5_Init(this)

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP5_type) :: this

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine SF_BF_KRP5_Init

! ************************************************************************** !

subroutine SF_BF_KRP5_Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_BF_KRP5_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,LINEAR'
  endif
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif   
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('Srg',string)
    call printErrMsg(option)
  endif   
end subroutine SF_BF_KRP5_Verify

! ************************************************************************** !

subroutine SF_BF_KRP5_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation
  ! From BRAGFLO 6.02 User Manual page 121 KRP=5.
  ! I also referenced the source code of BRAGFLO 6.02

  ! Author: Heeho Park
  ! Date: 11/18/16
  !
  use Option_module
  
  implicit none
  
  class(sat_func_BF_KRP5_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  Se = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  capillary_pressure = (1.d0/this%alpha-this%pcmax)*Se + this%pcmax

  capillary_pressure = min(capillary_pressure,this%pcmax)
  
end subroutine SF_BF_KRP5_CapillaryPressure

! ************************************************************************** !

subroutine SF_BF_KRP5_Saturation(this,capillary_pressure,liquid_saturation, &
                            dsat_dpres,option)
  ! 
  ! Computes the saturation as a function of capillary_pressure
  ! From BRAGFLO 6.02 User Manual page 121 KRP=5.
  ! I also referenced the source code of BRAGFLO 6.02

  ! Author: Heeho Park
  ! Date: 11/18/16
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_BF_KRP5_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dSe_dpc
  
  dsat_dpres = 0.d0

  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    Se = (this%pcmax-capillary_pressure) / (this%pcmax-1.d0/this%alpha)
    dSe_dpc = -1.d0/(this%pcmax-1.d0/this%alpha)
    liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se
    dsat_dpres = -(1.d0-this%Sr-this%Srg)*dSe_dpc
  endif 

end subroutine SF_BF_KRP5_Saturation
! End SF: BRAGFLO KRP5 Model
                            
! ************************************************************************** !

! Begin SF: BRAGFLO KRP9 Model
function SF_BF_KRP9_Create()

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP9_type), pointer :: SF_BF_KRP9_Create
  
  allocate(SF_BF_KRP9_Create)
  call SF_BF_KRP9_Create%Init()
  
end function SF_BF_KRP9_Create

! ************************************************************************** !

subroutine SF_BF_KRP9_Init(this)

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP9_type) :: this

  call SFBaseInit(this)
  
end subroutine SF_BF_KRP9_Init

! ************************************************************************** !

subroutine SF_BF_KRP9_Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_BF_KRP9_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP9'
  endif
  call SFBaseVerify(this,string,option)

end subroutine SF_BF_KRP9_Verify

! ************************************************************************** !

subroutine SF_BF_KRP9_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation
  ! based on experimental measurements and analyses done by Vauclin et al.
  ! as discussed by Moridis and Pruess, and the BRAGFLO V6.02 Requirements
  ! Document and Verification and Validation Plan, Sandia National Laboratories,
  ! Carlsbad, NM. ERMS #558659.  
  ! Moridis, G. J., and K. Pruess.  1992.  TOUGH Simulations of 
  ! Updegraff's Set of Fluid and Heat Flow Problems.  LBL-32611, ERMS# 138458. 
  ! Berkeley, CA:  Lawrence Berkeley Laboratory.
  ! Author: Heeho Park
  ! Date: 03/26/15
  !
  use Option_module
  
  implicit none
  
  class(sat_func_BF_KRP9_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal, parameter :: a = 3783.0145d0
  PetscReal, parameter :: b = 2.9d0
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = 0.d0
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  Se = (1.d0-liquid_saturation)/(liquid_saturation)
  capillary_pressure = a*Se**(1.d0/b)
#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
  endif
#endif

!  capillary_pressure = min(capillary_pressure,this%pcmax)
  
end subroutine SF_BF_KRP9_CapillaryPressure

! ************************************************************************** !

subroutine SF_BF_KRP9_Saturation(this,capillary_pressure,liquid_saturation, &
                            dsat_dpres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  !  
  ! Author: Heeho Park
  ! Date: 03/26/15
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_BF_KRP9_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dS_dSe
  PetscReal :: dSe_dpc
  PetscReal, parameter :: dpc_dpres = -1.d0
  PetscReal, parameter :: a = 3783.0145d0
  PetscReal, parameter :: b = 2.9d0
  
  dsat_dpres = 0.d0
  
  if (capillary_pressure <= 0.d0) then
    liquid_saturation = 1.d0
    return
  else
    Se = (capillary_pressure/a)**(b)
    liquid_saturation = 1.d0 / (Se+1.d0)
    ! Python analytical derivative (Jenn Frederick)
    dS_dSe = -1.d0/(Se + 1.d0)**2
    dSe_dpc = b*(capillary_pressure/a)**b/capillary_pressure
    dsat_dpres = dS_dSe*dSe_dpc*dpc_dpres
  endif 

end subroutine SF_BF_KRP9_Saturation
! End SF: BRAGFLO KRP9 Model

! ************************************************************************** !

! Begin SF: BRAGFLO KRP4 Model

function SF_BF_KRP4_Create()

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP4_type), pointer :: SF_BF_KRP4_Create
  
  allocate(SF_BF_KRP4_Create)
  call SF_BF_KRP4_Create%Init()
  
end function SF_BF_KRP4_Create

! ************************************************************************** !

subroutine SF_BF_KRP4_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(sat_func_BF_KRP4_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BROOKS_COREY'
  endif  
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif 
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif 
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('Srg',string)
    call printErrMsg(option)
  endif 
  if (Uninitialized(this%pcmax_flag)) then
    option%io_buffer = UninitializedMessage('KPC',string)
    call printErrMsg(option)
  endif 
  
end subroutine SF_BF_KRP4_Verify
  
! ************************************************************************** !

subroutine SF_BF_KRP4_CapillaryPressure(this,liquid_saturation, &
                                        capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation using the
  ! Brooks-Corey formulation
  ! 
  ! Modified according to KRP=4 option of BRAGFLO
  ! Explanation: residual gas saturation in the calculation of effective 
  ! saturation
  ! There is no usage of Pc Max unless KPC card is defined as 2. If KPC = 0,
  ! then there is no cut off in Pc Max
  ! Author: Heeho Park
  ! Date: 11/14/15
  !
  use Option_module
  use Utility_module
  
  implicit none
  
  class(sat_func_BF_KRP4_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dummy_real
  
  if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  if (this%alpha > 1.0d20) then
    capillary_pressure = 0.d0
    return
  endif

  Se = (liquid_saturation-this%Sr)/(1.d0-this%Sr-this%Srg)
  
  if (Se > 1.d0) then
     Se = 1.d0
  else if (Se < 0.d0) then
     Se = 0.d0
  endif

  if (associated(this%sat_poly)) then
    if (Se > this%sat_poly%low) then
      call QuadraticPolynomialEvaluate(this%sat_poly%coefficients(1:3), &
                                       Se,capillary_pressure,dummy_real)
      return
    endif
  endif

  capillary_pressure = (Se**(-1.d0/this%lambda))/this%alpha
  
  if (this%pcmax_flag == 2) then
    capillary_pressure = min(capillary_pressure,this%pcmax)
  endif

end subroutine SF_BF_KRP4_CapillaryPressure

! ************************************************************************** !

subroutine SF_BF_KRP4_Saturation(this,capillary_pressure,liquid_saturation, &
                                 dsat_dpres,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation using the
  ! Brooks-Corey formulation
  ! 
  ! Modified according to KRP=4 option of BRAGFLO
  ! Explanation: residual gas saturation in the calculation of effective 
  ! saturation
  ! There is no usage of Pc Max unless KPC card is defined as 2. If KPC = 0,
  ! then there is no cut off in Pc Max
  ! Author: Heeho Park
  ! Date: 11/14/15
  ! 
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_BF_KRP4_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  PetscReal :: pc_alpha_neg_lambda
  PetscReal :: Se
  PetscReal :: dSe_dpc
  PetscReal, parameter :: dpc_dpres = -1.d0
  
  dsat_dpres = 0.d0
  
  ! reference #1
  if (associated(this%pres_poly)) then
    if (capillary_pressure < this%pres_poly%low) then
      liquid_saturation = 1.d0
      return
    else if (capillary_pressure < this%pres_poly%high) then
      call CubicPolynomialEvaluate(this%pres_poly%coefficients, &
                                   capillary_pressure,Se,dSe_dpc)
      liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se
      dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
      return
    endif
  else
    if (capillary_pressure < 1.d0/this%alpha) then
      liquid_saturation = 1.d0
      dsat_dpres = 0.d0
      return
    endif
  endif

  pc_alpha_neg_lambda = (capillary_pressure*this%alpha)**(-this%lambda)
  Se = pc_alpha_neg_lambda
  dSe_dpc = -this%lambda/capillary_pressure*pc_alpha_neg_lambda
  liquid_saturation = this%Sr + (1.d0-this%Sr-this%Srg)*Se
  dsat_dpres = (1.d0-this%Sr)*dSe_dpc*dpc_dpres
  
end subroutine SF_BF_KRP4_Saturation

! End SF: BRAGFLO KRP4 Model

! ************************************************************************** !
! Begin SF: BRAGFLO KRP11 Model
function SF_BF_KRP11_Create()

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP11_type), pointer :: SF_BF_KRP11_Create
  
  allocate(SF_BF_KRP11_Create)
  call SF_BF_KRP11_Create%Init()
  
end function SF_BF_KRP11_Create

! ************************************************************************** !

subroutine SF_BF_KRP11_Init(this)

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP11_type) :: this

  call SFBaseInit(this)
  
end subroutine SF_BF_KRP11_Init

! ************************************************************************** !

subroutine SF_BF_KRP11_Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(sat_func_BF_KRP11_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BRAGFLO_KRP11'
  endif
  call SFBaseVerify(this,string,option)

end subroutine SF_BF_KRP11_Verify

! ************************************************************************** !

subroutine SF_BF_KRP11_CapillaryPressure(this,liquid_saturation, &
                                         capillary_pressure,option)
  ! 
  ! KRP=11 of BRAGFLO
  ! capillary pressure is 0 at all times
  ! Author: Heeho Park
  ! Date: 03/26/15
  !
  use Option_module
  
  implicit none
  
  class(sat_func_BF_KRP11_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option

  capillary_pressure = 0.0d0
  
end subroutine SF_BF_KRP11_CapillaryPressure

! ************************************************************************** !

subroutine SF_BF_KRP11_Saturation(this,capillary_pressure,liquid_saturation, &
                                  dsat_dpres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
  ! 
  !   
  ! Author: Heeho Park
  ! Date: 03/26/15
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(sat_func_BF_KRP11_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option
  
  dsat_dpres = 0.d0

  liquid_saturation = 1.d0

end subroutine SF_BF_KRP11_Saturation
! End SF: BRAGFLO KRP11 Model

! ************************************************************************** !

! Begin SF: BRAGFLO KRP12 Model

function SF_BF_KRP12_Create()

  ! Creates the van Genutchten capillary pressure function object

  implicit none
  
  class(sat_func_BF_KRP12_type), pointer :: SF_BF_KRP12_Create
  
  allocate(SF_BF_KRP12_Create)
  call SF_BF_KRP12_Create%Init()
  
end function SF_BF_KRP12_Create

! ************************************************************************** !

subroutine SF_BF_KRP12_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(sat_func_BF_KRP12_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 
  
  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,BROOKS_COREY'
  endif  
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%socmin)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif 
  if (Uninitialized(this%soceffmin)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif 
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('Srg',string)
    call printErrMsg(option)
  endif 
  
end subroutine SF_BF_KRP12_Verify
! ************************************************************************** !

subroutine SF_BF_KRP12_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  ! 
  ! Computes the capillary_pressure as a function of saturation using the
  ! Brooks-Corey formulation
  ! 
  ! Modified according to KRP=12 option of BRAGFLO
  ! Explanation: The relative permeabilities are unchanged from the 
  ! modified Brooks-Corey model but the capillary presssure is 
  ! calculated with modified saturation
  ! Author: Heeho Park
  ! Date: 11/14/15
  !
  use Option_module
  use Utility_module
  
  implicit none
  
  class(sat_func_BF_KRP12_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dummy_real
  PetscReal :: soczro
  
  soczro = this%socmin - this%soceffmin
  
  if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  if (this%alpha > 1.0d20) then
    capillary_pressure = 0.d0
    return
  endif

  Se = (liquid_saturation-soczro)/(1.d0-soczro)
  Se = max(min(Se,1.0d0),this%soceffmin)
  
  if (associated(this%sat_poly)) then
    if (Se > this%sat_poly%low) then
      call QuadraticPolynomialEvaluate(this%sat_poly%coefficients(1:3), &
                                       Se,capillary_pressure,dummy_real)
      return
    endif
  endif

  capillary_pressure = (Se**(-1.d0/this%lambda))/this%alpha

end subroutine SF_BF_KRP12_CapillaryPressure

! End SF: BRAGFLO KRP12 Model  

! ************************************************************************** !
function SF_mK_Create()

  ! Creates the modified Kosugi saturation function object

  implicit none

  class(sat_func_mK_type), pointer :: SF_mK_Create

  allocate(SF_mK_Create)
  call SF_mK_Create%Init()

end function SF_mK_Create
! ************************************************************************** !

subroutine SF_mK_Verify(this,name,option)

  use Option_module

  implicit none

  class(sat_func_mK_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'SATURATION_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'SATURATION_FUNCTION,MODIFIED_KOSUGI'
  endif
  call SFBaseVerify(this,string,option)
  if (Uninitialized(this%sigmaz)) then
    option%io_buffer = UninitializedMessage('SIGMAZ',string)
    call printErrMsg(option)
  endif
  if (Uninitialized(this%muz)) then
    option%io_buffer = UninitializedMessage('MUZ',string)
    call printErrMsg(option)
  endif
  if (Uninitialized(this%nparam)) then
    option%io_buffer = UninitializedMessage('NPARAM',string)
    call printErrMsg(option)
  endif
  if (Uninitialized(this%rmax)) then
    ! rmax is used for both nparam 3 and 4
    option%io_buffer = UninitializedMessage('RMAX',string)
    call printErrMsg(option)
  endif
  select case(this%nparam)
    case(4)
      ! r0 is only used for nparam 4
      if (Uninitialized(this%r0)) then
        option%io_buffer = UninitializedMessage('R0',string)
        call printErrMsg(option)
      endif
      if (this%r0 >= this%rmax) then
        option%io_buffer = trim(string) // ' requires RMAX > R0'
        call printErrMsg(option)
      end if
    case(3)
      continue ! rmax handled above
    case default
      option%io_buffer = 'invalid NPARAM value in' // &
        trim(string) // '. Only NPARAM=(3,4) supported.'
      call printErrMsg(option)
  end select

end subroutine SF_MK_Verify

! ************************************************************************** !

subroutine SF_mK_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,option)
  !
  ! Computes the capillary_pressure as a function of saturation
  ! for modified Kosugi model.
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498–502. http://dx.doi.org/10.1111/gwat.12220
  !
  ! Author: Kris Kuhlman
  ! Date: 2017
  !
  use Option_module
  use Utility_module, only : InverseNorm

  implicit none

  PetscReal, parameter :: KAPPA = 1.49D-1 !  water in glass tube
  PetscReal, parameter :: LNKAP = log(KAPPA)
  PetscReal, parameter :: UNIT_CONVERSION = 9.982D+2*9.81d0/1.0D+2
  
  class(sat_func_mK_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  type(option_type), intent(inout) :: option

  PetscReal :: Se
  PetscReal :: inverse, exparg
  PetscReal :: hc, hmaxinv

  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif

  Se = (liquid_saturation - this%Sr)/(1.d0 - this%Sr)
  inverse = -InverseNorm(Se)
  exparg = this%sigmaz*inverse + LNKAP - this%muz

  hc = KAPPA/this%rmax
  capillary_pressure = exp(exparg) + hc
  if (this%nparam == 4) then
    hmaxinv = this%r0/KAPPA
    capillary_pressure = 1.d0/(1.d0/capillary_pressure + hmaxinv)
  end if

  capillary_pressure = min(capillary_pressure*UNIT_CONVERSION,this%pcmax)

end subroutine SF_mK_CapillaryPressure

! ************************************************************************** !

subroutine SF_mK_Saturation(this,capillary_pressure,liquid_saturation, &
                            dsat_dpres,option)
  !
  ! Computes the saturation (and associated derivatives) as a function of
  ! capillary pressure for modified Kosugi model
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498–502. http://dx.doi.org/10.1111/gwat.12220
  !
  ! Author: Kris Kuhlman
  ! Date: 2017
  !
  use Option_module
  use Utility_module, only : InverseNorm

  implicit none

  ! gnu & intel extension and required in f2008
  intrinsic :: erfc

  PetscReal, parameter :: KAPPA = 1.49D-1 ! water in glass tube
  PetscReal, parameter :: LNKAP = log(KAPPA)
  PetscReal, parameter :: SQRT2 = sqrt(2.0d0)
  PetscReal, parameter :: SQRTPI = sqrt(4.0d0*atan(1.0d0))
  PetscReal, parameter :: UNIT_CONVERSION = 9.982D+2*9.81d0/1.0D+2
  
  class(sat_func_mK_type) :: this
  PetscReal, intent(in) :: capillary_pressure
  PetscReal, intent(out) :: liquid_saturation
  PetscReal, intent(out) :: dsat_dpres
  type(option_type), intent(inout) :: option

  PetscReal :: hc, hmax, cap_press_scaled
  PetscReal :: rt2sz
  PetscReal :: lnArg, erfcArg

  dsat_dpres = 0.0d0
  cap_press_scaled = capillary_pressure/UNIT_CONVERSION
  
  hc = KAPPA/this%rmax
  if (cap_press_scaled <= hc) then
    liquid_saturation = 1.d0
    return
  end if

  if (this%nparam == 3) then
    lnArg = cap_press_scaled - hc
  else ! nparam == 4 
    hmax = KAPPA/this%r0
    if (cap_press_scaled >= hmax) then
      liquid_saturation = this%Sr
      return
    end if
    lnArg = 1.d0/(1.d0/cap_press_scaled - 1.d0/hmax) - hc
  end if

  rt2sz = SQRT2*this%sigmaz
  erfcArg = (log(lnArg) - LNKAP + this%muz)/rt2sz
  liquid_saturation = this%Sr + (1.0d0-this%Sr)*5.0D-1*erfc(erfcArg)
  dsat_dpres = exp(-erfcArg**2)/(SQRTPI*rt2sz*lnArg)/UNIT_CONVERSION

end subroutine SF_mK_Saturation
! End SF: modified Kosugi

! ************************************************************************** !

! Begin RPF: Mualem, Van Genuchten (Liquid)
function RPF_Mualem_VG_Liq_Create()

  ! Creates the van Genutchten Mualem relative permeability function object

  implicit none
  
  class(rpf_Mualem_vg_liq_type), pointer :: RPF_Mualem_VG_Liq_Create
  
  allocate(RPF_Mualem_VG_Liq_Create)
  call RPF_Mualem_VG_Liq_Create%Init()
  
end function RPF_Mualem_VG_Liq_Create

! ************************************************************************** !

subroutine RPF_Mualem_VG_Liq_Init(this)

  ! Initializes the van Genutchten Mualem relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_VG_liq_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Mualem_VG_Liq_Init

! ************************************************************************** !

subroutine RPF_Mualem_VG_Liq_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Mualem_VG_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif   
  
end subroutine RPF_Mualem_VG_Liq_Verify

! ************************************************************************** !

subroutine RPF_Mualem_SetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing Mualem permeability function

  use Option_module
  use Utility_module
  
  implicit none
  
  class(rpf_Mualem_VG_liq_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  PetscReal :: b(4)
  PetscReal :: one_over_m, Se_one_over_m, m

  if (associated(this%poly)) call PolynomialDestroy(this%poly)
  this%poly => PolynomialCreate()
  ! fill matix with values
  this%poly%low = 0.99d0  ! just below saturated
  this%poly%high = 1.d0   ! saturated
  
  m = this%m
  one_over_m = 1.d0/m
  Se_one_over_m = this%poly%low**one_over_m
  b(1) = 1.d0
  b(2) = sqrt(this%poly%low)*(1.d0-(1.d0-Se_one_over_m)**m)**2.d0
  b(3) = 0.d0
  b(4) = 0.5d0*b(2)/this%poly%low+ &
          2.d0*this%poly%low**(one_over_m-0.5d0)* &
          (1.d0-Se_one_over_m)**(m-1.d0)* &
          (1.d0-(1.d0-Se_one_over_m)**m)
  
  call CubicPolynomialSetup(this%poly%high,this%poly%low,b)
  
  this%poly%coefficients(1:4) = b(1:4)
  
end subroutine RPF_Mualem_SetupPolynomials

! ************************************************************************** !

subroutine RPF_Mualem_VG_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Mualem_VG_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: one_over_m
  PetscReal :: Se_one_over_m
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  if (associated(this%poly)) then
    if (Se > this%poly%low) then
      call CubicPolynomialEvaluate(this%poly%coefficients, &
                                   Se,relative_permeability,dkr_Se)
      return
    endif
  endif
  
  one_over_m = 1.d0/this%m
  Se_one_over_m = Se**one_over_m
  relative_permeability = sqrt(Se)*(1.d0-(1.d0-Se_one_over_m)**this%m)**2.d0
  dkr_Se = 0.5d0*relative_permeability/Se+ &
            2.d0*Se**(one_over_m-0.5d0)* &
                (1.d0-Se_one_over_m)**(this%m-1.d0)* &
                (1.d0-(1.d0-Se_one_over_m)**this%m)
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPF_Mualem_VG_Liq_RelPerm
! End RPF: Mualem, Van Genuchten (Liquid)

! ************************************************************************** !

! Begin RPF: Mualem, Van Genuchten (Gas)
function RPF_Mualem_VG_Gas_Create()

  ! Creates the van Genutchten Mualem gas relative permeability function object

  implicit none
  
  class(rpf_Mualem_VG_gas_type), pointer :: RPF_Mualem_VG_Gas_Create
  
  allocate(RPF_Mualem_VG_Gas_Create)
  call RPF_Mualem_VG_Gas_Create%Init()
  
end function RPF_Mualem_VG_Gas_Create

! ************************************************************************** !

subroutine RPF_Mualem_VG_Gas_Init(this)

  ! Initializes the van Genutchten Mualem gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_VG_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Mualem_VG_Gas_Init

! ************************************************************************** !

subroutine RPF_Mualem_VG_Gas_Verify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rpf_Mualem_VG_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM_VG_GAS'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif 
  
end subroutine RPF_Mualem_VG_Gas_Verify

! ************************************************************************** !

subroutine RPF_Mualem_VG_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Mualem_VG_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  relative_permeability = sqrt(Seg)*(1.d0-Se**(1.d0/this%m))**(2.d0*this%m)
  ! Mathematica analytical solution (Heeho Park)
  dkr_Se = -(1.d0-Se**(1.d0/this%m))**(2.d0*this%m)/(2.d0*sqrt(Seg)) &
          - 2.d0*sqrt(Seg)*Se**(1.d0/this%m-1.d0) &
          * (1.d0-Se**(1.d0/this%m))**(2.d0*this%m-1.d0)
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPF_Mualem_VG_Gas_RelPerm
! End RPF: Mualem, Van Genuchten (Gas)

! ************************************************************************** !

! RPF: Tough2 IRP7 w/ VG-Mualem (Gas)
function RPF_TOUGH2_IRP7_Gas_Create()

  ! Creates the Brooks-Corey Burdine gas relative permeability function object

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type), pointer :: RPF_TOUGH2_IRP7_Gas_Create
  
  allocate(RPF_TOUGH2_IRP7_Gas_Create)
  call RPF_TOUGH2_IRP7_Gas_Create%Init()
  
end function RPF_TOUGH2_IRP7_Gas_Create

! ************************************************************************** !

subroutine RPF_TOUGH2_IRP7_Gas_Init(this)

  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_TOUGH2_IRP7_Gas_Init

! ************************************************************************** !

subroutine RPF_TOUGH2_IRP7_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_TOUGH2_IRP7_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,TOUGH2_IRP7_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_TOUGH2_IRP7_Gas_Verify

! ************************************************************************** !

subroutine RPF_TOUGH2_IRP7_Gas_RelPerm(this,liquid_saturation, &
                                       relative_permeability,dkr_sat,option)
  ! 
  ! TOUGH2 IRP(7) equations from Appendix G of TOUGH2 user manual
  !
  use Option_module
  
  implicit none

  class(rpf_TOUGH2_IRP7_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: liquid_relative_permeability
  PetscReal :: liquid_dkr_sat
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  dkr_sat = dkr_sat / 0.d0
  dkr_sat = dkr_sat * 0.d0
  
                 ! essentially zero
  if (this%Srg <= 0.d0) then
    call RPF_Mualem_VG_Liq_RelPerm(this,liquid_saturation, &
                            liquid_relative_permeability, &
                            liquid_dkr_sat,option)
    relative_permeability = 1.d0 - liquid_relative_permeability
    return
  endif  
  
  if ((1.d0 - liquid_saturation) <= this%Srg) then
    relative_permeability = 0.d0
  else
    Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
    Seg = 1.d0 - Se
    relative_permeability = Seg**2*(1.d0-Se*Se)
    ! Mathematica Analytical solution (Heeho Park)
    dkr_Se = -2.d0*Seg**2.d0*Se - 2.d0*Seg*(1.d0-Se**2.d0)
    dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
    dkr_sat = dkr_Se * dSe_sat
  endif
    
end subroutine RPF_TOUGH2_IRP7_Gas_RelPerm
! End RPF: Tough2 IRP7 w/ VG-Mualem (Gas)

! ************************************************************************** !

! Begin RPF: Burdine, Brooks-Corey (Liquid)
function RPF_Burdine_BC_Liq_Create()

  ! Creates the Brooks-Corey Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_BC_Liq_type), pointer :: RPF_Burdine_BC_Liq_Create
  
  allocate(RPF_Burdine_BC_Liq_Create)
  call RPF_Burdine_BC_Liq_Create%Init()
  
end function RPF_Burdine_BC_Liq_Create

! ************************************************************************** !

subroutine RPF_Burdine_BC_Liq_Init(this)

  ! Initializes the Brooks-Corey Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_BC_Liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Burdine_BC_Liq_Init

! ************************************************************************** !

subroutine RPF_Burdine_BC_Liq_Verify(this,name,option)

  ! Initializes the Brooks-Corey Burdine relative permeability function object

  use Option_module
  
  implicit none
  
  class(rpf_Burdine_BC_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE'
  endif    
  call RPFBaseVerify(this,name,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif
  
end subroutine RPF_Burdine_BC_Liq_Verify

! ************************************************************************** !

subroutine RPF_Burdine_BC_Liq_SetupPolynomials(this,option,error_string)

  ! Sets up smoothing approach for Burdine-BC permeability function
  ! F.-M. Yuan (2017-03-10): applying for newly-written HFunctionSmooth to do tail-smoothing

  use Option_module
  use Utility_module

  implicit none

  class(rpf_Burdine_BC_liq_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string

#ifdef SMOOTHING2
  ! F.-M. Yuan (2017-03-10): applying for newly-written HFunctionSmooth to do tail-smoothing
  ! Here, we only need to know the smoothing startng/ending points, if with option of SMOOTH keyword in input cards

  ! WET-end smoothing
  if (associated(this%poly)) call PolynomialDestroy(this%poly)
  this%poly => PolynomialCreate()
  this%poly%low  = 0.999d0 ! normalized effective saturation (Se)  (starting smoothing point)
  this%poly%high = 1.000d0 ! normalized effective saturation (Se)  (saturated flow above this point)

  ! DRY-end smoothing
  if (associated(this%poly2)) call PolynomialDestroy(this%poly2)
  this%poly2 => PolynomialCreate()
  this%poly2%high = 0.050d0 ! normalized effective saturation (Se) (starting smoothing point)
  this%poly2%low  = 0.001d0 ! normalized effective saturation (Se)
                            !   (note that by BC function, Se may not never be ZERO)

#endif

end subroutine RPF_Burdine_BC_Liq_SetupPolynomials


! ************************************************************************** !

subroutine RPF_Burdine_BC_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Burdine_BC_Liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: power
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

#ifdef SMOOTHING2
  PetscReal :: Hfunc, dHfunc
#endif


  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  ! reference #1
  power = 3.d0+2.d0/this%lambda
  relative_permeability = Se**power
  dkr_Se = power*relative_permeability/Se          
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat

  ! smoothing curve at both dry and wet endings (F.-M. Yuan, 2017-03-10)
#ifdef SMOOTHING2
  ! WET-end smoothing up to 1.0@Se=1.0 (ending high point)
  if (associated(this%poly)) then
    if (Se>=this%poly%low) then
      call HFunctionSmooth(Se, this%poly%low, this%poly%high, Hfunc, dHfunc)

      !rpf = 1.0 + (rpf-1.0) * Hfunc, in which '1.0' is the base RP@Se=1.0
      dkr_sat = dkr_sat*Hfunc + (relative_permeability-1.0d0)*dHfunc
      relative_permeability = 1.0d0+(relative_permeability-1.0d0)*Hfunc

    endif

  endif

  ! DRY-end smoothing down to 0.0@Se=ending low point
  if (associated(this%poly2)) then
    if (Se<=this%poly2%high) then
      call HFunctionSmooth(Se, this%poly2%high, this%poly2%low, Hfunc, dHfunc)

      !rpf = rpf * Hfunc ( or, rpf = 0.0 + (rpf - 0.0)*Hfunc)
      dkr_sat = dkr_sat*Hfunc + relative_permeability*dHfunc
      relative_permeability = relative_permeability*Hfunc

    endif

  endif
#endif

end subroutine RPF_Burdine_BC_Liq_RelPerm
! End RPF: Burdine, Brooks-Corey (Liquid)

! ************************************************************************** !

! Begin RPF: Burdine, Brooks-Corey (Gas)
function RPF_Burdine_BC_Gas_Create()

  ! Creates the Brooks-Corey Burdine gas relative permeability function
  ! object

  implicit none
  
  class(rpf_Burdine_BC_gas_type), pointer :: RPF_Burdine_BC_Gas_Create
  
  allocate(RPF_Burdine_BC_Gas_Create)
  call RPF_Burdine_BC_Gas_Create%Init()
  
end function RPF_Burdine_BC_Gas_Create

! ************************************************************************** !

subroutine RPF_Burdine_BC_Gas_Init(this)

  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Burdine_BC_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Burdine_BC_Gas_Init

! ************************************************************************** !

subroutine RPF_Burdine_BC_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Burdine_BC_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE_BC_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_Burdine_BC_Gas_Verify

! ************************************************************************** !

subroutine RPF_Burdine_BC_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Burdine_BC_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  ! reference #1
  relative_permeability = Seg*Seg*(1.d0-Se**(1.d0+2.d0/this%lambda))
  ! Mathematica analytical solution (Heeho Park)
  dkr_Se = -(1.d0+2.d0/this%lambda)*Seg**2.d0*Se**(2.d0/this%lambda) &
           - 2.d0*Seg*(1.d0-Se**(1.d0+2.d0/this%lambda))
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPF_Burdine_BC_Gas_RelPerm
! End RPF: Burdine, Brooks-Corey (Gas)

! ************************************************************************** !

! Begin RPF: Mualem, Brooks-Corey (Liq)
function RPF_Mualem_BC_Liq_Create()

  ! Creates the Brooks-Corey Mualem liquid relative permeability function object

  implicit none
  
  class(rpf_Mualem_BC_liq_type), pointer :: RPF_Mualem_BC_Liq_Create
  
  allocate(RPF_Mualem_BC_Liq_Create)
  call RPF_Mualem_BC_Liq_Create%Init()
  
end function RPF_Mualem_BC_Liq_Create

! ************************************************************************** !

subroutine RPF_Mualem_BC_Liq_Init(this)

  ! Initializes the Brooks-Corey Mualem liquid relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_BC_liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE

end subroutine RPF_Mualem_BC_Liq_Init

! ************************************************************************** !

subroutine RPF_Mualem_BC_Liq_Verify(this,name,option)

  ! Initializes the Brooks-Corey Mualem liquid relative permeability function object

  use Option_module
  
  implicit none
  
  class(rpf_Mualem_BC_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM'
  endif    
  call RPFBaseVerify(this,name,option)
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif
  
end subroutine RPF_Mualem_BC_Liq_Verify

! ************************************************************************** !

subroutine RPF_Mualem_BC_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Mualem_BC_Liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: power
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  ! reference #1
  power = 2.5d0+2.d0/this%lambda
  relative_permeability = Se**power
  dkr_Se = power*relative_permeability/Se          
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat 

end subroutine RPF_Mualem_BC_Liq_RelPerm
! End RPF: Mualem, Brooks-Corey (Liq)

! ************************************************************************** !

! Begin RPF: Mualem, Brooks-Corey (Gas)
function RPF_Mualem_BC_Gas_Create()

  ! Creates the Brooks-Corey Mualem gas relative permeability function object

  implicit none
  
  class(rpf_Mualem_BC_gas_type), pointer :: RPF_Mualem_BC_Gas_Create
  
  allocate(RPF_Mualem_BC_Gas_Create)
  call RPF_Mualem_BC_Gas_Create%Init()
  
end function RPF_Mualem_BC_Gas_Create

! ************************************************************************** !

subroutine RPF_Mualem_BC_Gas_Init(this)

  ! Initializes the Brooks-Corey Mualem gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_BC_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Mualem_BC_Gas_Init

! ************************************************************************** !

subroutine RPF_Mualem_BC_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Mualem_BC_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM_BC_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_Mualem_BC_Gas_Verify

! ************************************************************************** !

subroutine RPF_Mualem_BC_Gas_RelPerm(this,liquid_saturation, &
                                       relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_Mualem_BC_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  ! reference Table 2
  relative_permeability = sqrt(Seg)* &
                             (1.d0-Se**(1.d0+1.d0/this%lambda))**2.d0
  ! Mathematica analytical solution (Heeho Park)
  dkr_Se = -2.d0*(1.d0+1.d0/this%lambda)*sqrt(Seg)*Se**(1.d0/this%lambda) &
          * (1.d0-Se**(1.d0+1.d0/this%lambda)) &
          - (1.d0-Se**(1.d0+1.d0/this%lambda))**2.d0/(2.d0*sqrt(Seg))
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPF_Mualem_BC_Gas_RelPerm
! End RPF: Mualem, Brooks-Corey (Gas)

! ************************************************************************** !

! Begin RPF: Burdine, Van Genuchten (Liq)
function RPF_Burdine_VG_Liq_Create()

  ! Creates the van Genutchten Mualem relative permeability function object

  implicit none
  
  class(rpf_burdine_vg_liq_type), pointer :: RPF_Burdine_VG_Liq_Create
  
  allocate(RPF_Burdine_VG_Liq_Create)
  call RPF_Burdine_VG_Liq_Create%Init()
  
end function RPF_Burdine_VG_Liq_Create

! ************************************************************************** !

subroutine RPF_Burdine_VG_Liq_Init(this)

  ! Initializes the van Genutchten Mualem relative permeability function object

  implicit none
  
  class(rpf_Burdine_VG_liq_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Burdine_VG_Liq_Init

! ************************************************************************** !

subroutine RPF_Burdine_VG_Liq_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Burdine_VG_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('M',string)
    call printErrMsg(option)
  endif   
  
end subroutine RPF_Burdine_VG_Liq_Verify

! ************************************************************************** !

subroutine RPF_Burdine_VG_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Burdine_VG_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: one_over_m
  PetscReal :: Se_one_over_m
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  one_over_m = 1.d0/this%m
  Se_one_over_m = Se**one_over_m
  relative_permeability = Se*Se*(1.d0-(1.d0-Se_one_over_m)**this%m)
  dkr_Se = 2.d0*relative_permeability/Se + &
                 Se*Se_one_over_m*(1.d0-Se_one_over_m)**(this%m-1.d0)
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPF_Burdine_VG_Liq_RelPerm
! End RPF: Burdine, Van Genuchten (Liq)

! ************************************************************************** !

! Begin RPF: Burdine, Van Genuchten (Gas)
function RPF_Burdine_VG_Gas_Create()

  ! Creates the Brooks-Corey Burdine gas relative permeability function object

  implicit none
  
  class(rpf_Burdine_VG_gas_type), pointer :: RPF_Burdine_VG_Gas_Create
  
  allocate(RPF_Burdine_VG_Gas_Create)
  call RPF_Burdine_VG_Gas_Create%Init()
  
end function RPF_Burdine_VG_Gas_Create

! ************************************************************************** !

subroutine RPF_Burdine_VG_Gas_Init(this)

  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Burdine_VG_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Burdine_VG_Gas_Init

! ************************************************************************** !

subroutine RPF_Burdine_VG_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Burdine_VG_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE_VG_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_Burdine_VG_Gas_Verify

! ************************************************************************** !

subroutine RPF_Burdine_VG_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14

  use Option_module
  
  implicit none

  class(rpf_Burdine_VG_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0

  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  ! reference Table 2
  relative_permeability = Seg*Seg*(1.d0-Se**(1.d0/this%m))**this%m
  dkr_Se = -Seg**2.d0*Se**(1.d0/this%m-1.d0) &
          *(1.d0-Se**(1.d0/this%m))**(this%m-1.d0) &
          - 2.d0*Seg*(1.d0-Se**(1.d0/this%m))**this%m
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPF_Burdine_VG_Gas_RelPerm
! End RPF: Burdine, Van Genuchten (Gas)

! ************************************************************************** !

! Begin RPF: Mualem, Linear (Liquid)
function RPF_Mualem_Linear_Liq_Create()

  ! Creates the Linear Mualem relative permeability function object

  implicit none
  
  class(rpf_Mualem_linear_liq_type), pointer :: RPF_Mualem_Linear_Liq_Create
  
  allocate(RPF_Mualem_Linear_Liq_Create)
  call RPF_Mualem_Linear_Liq_Create%Init()
  
end function RPF_Mualem_Linear_Liq_Create

! ************************************************************************** !

subroutine RPF_Mualem_Linear_Liq_Init(this)

  ! Initializes the Linear Mualem relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_Linear_liq_type) :: this

  call RPFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%pcmax = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Mualem_Linear_Liq_Init

! ************************************************************************** !

subroutine RPF_Mualem_Linear_Liq_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Mualem_Linear_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif  
  if (Uninitialized(this%pcmax)) then
    option%io_buffer = UninitializedMessage('MAX_CAPILLARY_PRESSURE',string)
    call printErrMsg(option)
  endif
  
end subroutine RPF_Mualem_Linear_Liq_Verify

! ************************************************************************** !

subroutine RPF_Mualem_Linear_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  !   
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Mualem_Linear_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: one_over_alpha
  PetscReal :: pct_over_pcmax
  PetscReal :: pc_over_pcmax
  PetscReal :: pc_log_ratio
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  one_over_alpha = 1.d0/this%alpha
  pct_over_pcmax = one_over_alpha/this%pcmax
  pc_over_pcmax = 1.d0-(1.d0-pct_over_pcmax)*Se
  pc_log_ratio = log(pc_over_pcmax) / log(pct_over_pcmax)
  relative_permeability = (Se**0.5d0)*(pc_log_ratio**2.d0)
  ! ***used Mathematica to verify***
  ! In[3]:
  ! D[Se^(1/2)*(Log[1 - (1 - pctoverpcmax)*Se]/Log[pctoverpcmax])^2, Se]
  ! Out[3]:
  ! (2 (-1 + pctoverpcmax) Sqrt[Se]
  !  Log[1 - (1 - pctoverpcmax) Se])/((1 - (1 - pctoverpcmax) Se) Log[
  !  pctoverpcmax]^2) + Log[1 - (1 - pctoverpcmax) Se]^2/(
  ! 2 Sqrt[Se] Log[pctoverpcmax]^2)
  dkr_Se = 2.d0*(-1.d0+pct_over_pcmax)*sqrt(Se)* log(pc_over_pcmax) / &
    (pc_over_pcmax*log(pct_over_pcmax)**2.d0) + &
    log(pc_over_pcmax)**2.d0 / (2.d0*sqrt(Se)*log(pct_over_pcmax)**2.d0)
  dSe_sat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_Se * dSe_sat
         
end subroutine RPF_Mualem_Linear_Liq_RelPerm
! End RPF: Mualem, Linear (Liquid)

! ************************************************************************** !

! Begin RPF: Mualem, Linear (Gas)
function RPF_Mualem_Linear_Gas_Create()

  ! Creates the Linear Mualem gas relative permeability function object

  implicit none
  
  class(rpf_Mualem_Linear_gas_type), pointer :: RPF_Mualem_Linear_Gas_Create
  
  allocate(RPF_Mualem_Linear_Gas_Create)
  call RPF_Mualem_Linear_Gas_Create%Init()
  
end function RPF_Mualem_Linear_Gas_Create
 
! ************************************************************************** !

subroutine RPF_Mualem_Linear_Gas_Init(this)

  ! Initializes the Linear Mualem gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Mualem_Linear_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  this%alpha = UNINITIALIZED_DOUBLE
  this%pcmax = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Mualem_Linear_Gas_Init

! ************************************************************************** !

subroutine RPF_Mualem_Linear_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Mualem_Linear_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM_LINEAR_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  if (Uninitialized(this%alpha)) then
    option%io_buffer = UninitializedMessage('ALPHA',string)
    call printErrMsg(option)
  endif  
  if (Uninitialized(this%pcmax)) then
    option%io_buffer = UninitializedMessage('MAX_CAPILLARY_PRESSURE',string)
    call printErrMsg(option)
  endif
  
end subroutine RPF_Mualem_Linear_Gas_Verify

! ************************************************************************** !

subroutine RPF_Mualem_Linear_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14

  use Option_module
  
  implicit none

  class(rpf_Mualem_Linear_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: liquid_relative_permeability
  PetscReal :: liquid_dkr_sat  
  PetscReal :: dkr_dSe
  PetscReal :: dSe_dsat
  
  call RPF_Mualem_Linear_Liq_RelPerm(this,liquid_saturation, &
                                     liquid_relative_permeability, &
                                     liquid_dkr_sat,option)
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0

  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  ! reference Table 2
  relative_permeability = Seg**0.5d0 * &
                 (1.d0-sqrt(liquid_relative_permeability*Se**(-0.5d0)))**2.d0
  ! Python analytical derivative (Jenn Frederick)
  dkr_dSe = 0.5d0*1.d0/Se*sqrt(Se**(-0.5d0)*liquid_relative_permeability)* &
    sqrt(1.d0-Se)*(1.d0-sqrt(Se**(-0.5d0)*liquid_relative_permeability))**1.0 &
    - (1.d0-sqrt(Se**(-0.5d0)*liquid_relative_permeability))**2.d0 &
    /(2.d0*sqrt(1.d0-Se))
  !one_over_apcm = 1.d0/(1.d-7)/(1.d9)
  !dkr_dSe = -2.0*Se**0.5*sqrt(Se**(-0.5)*log(-Se*(-one_over_apcm + 1.0) + 1.0)/log(one_over_apcm))*sqrt(-Se + 1.0)*(-0.25*Se**(-1.5)*log(-Se*(-one_over_apcm + 1.0) + 1.0)/log(one_over_apcm) + Se**(-0.5)*(one_over_apcm - 1.0)/(2*(-Se*(-one_over_apcm + 1.0) + 1.0)*log(one_over_apcm)))*(-sqrt(Se**(-0.5)*log(-Se*(-one_over_apcm + 1.0) + 1.0)/log(one_over_apcm)) + 1.0)**1.0*log(one_over_apcm)/log(-Se*(-one_over_apcm + 1.0) + 1.0) - (-sqrt(Se**(-0.5)*log(-Se*(-one_over_apcm + 1.0) + 1.0)/log(one_over_apcm)) + 1.0)**2.0/(2*sqrt(-Se + 1.0))
  dSe_dsat = 1.d0/(1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_dSe*dSe_dsat

end subroutine RPF_Mualem_Linear_Gas_RelPerm
! End RPF: Mualem, Linear (Gas)

! ************************************************************************** !

! Begin RPF: Burdine, Linear (Liquid)
function RPF_Burdine_Linear_Liq_Create()

  ! Creates the Linear Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_linear_liq_type), pointer :: RPF_Burdine_Linear_Liq_Create
  
  allocate(RPF_Burdine_Linear_Liq_Create)
  call RPF_Burdine_Linear_Liq_Create%Init()
  
end function RPF_Burdine_Linear_Liq_Create

! ************************************************************************** !

subroutine RPF_Burdine_Linear_Liq_Init(this)

  ! Initializes the Linear Burdine relative permeability function 
  ! object

  implicit none
  
  class(rpf_Burdine_Linear_liq_type) :: this

  call RPFBaseInit(this)
  
end subroutine RPF_Burdine_Linear_Liq_Init

! ************************************************************************** !

subroutine RPF_Burdine_Linear_Liq_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Burdine_Linear_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE'
  endif  
  call RPFBaseVerify(this,string,option)
  
end subroutine RPF_Burdine_Linear_Liq_Verify

! ************************************************************************** !

subroutine RPF_Burdine_Linear_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !  
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  ! 
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Burdine_Linear_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  relative_permeability = Se
  dkr_sat = 1.d0 / (1.d0 - this%Sr)
  
end subroutine RPF_Burdine_Linear_Liq_RelPerm
! End RPF: Burdine, Linear (Liquid)

! ************************************************************************** !

! Begin Burdine Linear (Gas)
function RPF_Burdine_Linear_Gas_Create()

  ! Creates the Linear Burdine gas relative permeability function object

  implicit none
  
  class(rpf_Burdine_Linear_gas_type), pointer :: RPF_Burdine_Linear_Gas_Create
  
  allocate(RPF_Burdine_Linear_Gas_Create)
  call RPF_Burdine_Linear_Gas_Create%Init()
  
end function RPF_Burdine_Linear_Gas_Create

! ************************************************************************** !

subroutine RPF_Burdine_Linear_Gas_Init(this)

  ! Initializes the Linear Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_Burdine_Linear_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_Burdine_Linear_Gas_Init

! ************************************************************************** !

subroutine RPF_Burdine_Linear_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_Burdine_Linear_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BURDINE_LINEAR_GAS&
             &/BRAGFLO_ KRP5'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_Burdine_Linear_Gas_Verify

! ************************************************************************** !

subroutine RPF_Burdine_Linear_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  !
  !
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  !

  use Option_module
  
  implicit none

  class(rpf_Burdine_Linear_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif
  
  Seg = 1.d0 - Se
  relative_permeability = Seg
  dkr_Se = -1.d0
  dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  dkr_sat = dkr_Se * dSe_sat
  
end subroutine RPF_Burdine_Linear_Gas_RelPerm
! End Burdine Linear (Gas)

! ************************************************************************** !

! Begin RPF: BRAGFLO KRP1 (Liq)
function RPF_BRAGFLO_KRP1_Liq_Create()

  ! Creates the KRP1 or VG_Mualem liq relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP1_liq_type), pointer :: RPF_BRAGFLO_KRP1_Liq_Create
  
  allocate(RPF_BRAGFLO_KRP1_Liq_Create)
  call RPF_BRAGFLO_KRP1_Liq_Create%Init()
  
end function RPF_BRAGFLO_KRP1_Liq_Create
! End RPF: BRAGFLO KRP4 (Liq)

! ************************************************************************** !

! Begin RPF: BRAGFLO KRP1 (Gas)
function RPF_BRAGFLO_KRP1_Gas_Create()

  ! Creates the KRP1 or VG_Mualem liq relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP1_gas_type), pointer :: RPF_BRAGFLO_KRP1_Gas_Create
  
  allocate(RPF_BRAGFLO_KRP1_Gas_Create)
  call RPF_BRAGFLO_KRP1_Gas_Create%Init()
  
end function RPF_BRAGFLO_KRP1_Gas_Create
! End RPF: BRAGFLO KRP4 (Gas)

! ************************************************************************** !

! Begin RPF: BRAGFLO KRP5 (Liquid)
function RPF_BRAGFLO_KRP5_Liq_Create()

  ! Creates the BRAGFLO KRP5 relative permeability function object
  ! Linear Burdine with gas residual saturation.

  implicit none
  
  class(rpf_BRAGFLO_KRP5_liq_type), pointer :: RPF_BRAGFLO_KRP5_Liq_Create
  
  allocate(RPF_BRAGFLO_KRP5_Liq_Create)
  call RPF_BRAGFLO_KRP5_Liq_Create%Init()
  
end function RPF_BRAGFLO_KRP5_Liq_Create

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP5_Liq_Init(this)

  ! Initializes the BRAGFLO KRP5 relative permeability function 
  ! object

  implicit none
  
  class(rpf_BRAGFLO_KRP5_liq_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_BRAGFLO_KRP5_Liq_Init

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP5_Liq_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_BRAGFLO_KRP5_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP5'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
    
end subroutine RPF_BRAGFLO_KRP5_Liq_Verify

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP5_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  !
  ! Author: Heeho Park
  ! Date: 11/18/16
  !
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_BRAGFLO_KRP5_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif
  
  relative_permeability = Se
  dkr_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
  
end subroutine RPF_BRAGFLO_KRP5_Liq_RelPerm
! End RPF: BRAGFLO KRP5 (Liquid)

! ************************************************************************** !

! Begin BRAGFLO KRP5 (Gas)
function RPF_BRAGFLO_KRP5_Gas_Create()

  ! Creates the BRAGFLO KRP5 gas relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP5_gas_type), pointer :: RPF_BRAGFLO_KRP5_Gas_Create
  
  allocate(RPF_BRAGFLO_KRP5_Gas_Create)
  call RPF_BRAGFLO_KRP5_Gas_Create%Init()
  
end function RPF_BRAGFLO_KRP5_Gas_Create
! End RPF: BRAGFLO KRP5 (Liquid)

! ************************************************************************** !
                                     
! Begin RPF: BRAGFLO KRP9 (Liquid)
function RPF_BRAGFLO_KRP9_Liq_Create()

  ! Creates the Linear Burdine relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP9_liq_type), pointer :: RPF_BRAGFLO_KRP9_Liq_Create
  
  allocate(RPF_BRAGFLO_KRP9_Liq_Create)
  call RPF_BRAGFLO_KRP9_Liq_Create%Init()
  
end function RPF_BRAGFLO_KRP9_Liq_Create

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP9_Liq_Init(this)

  ! Initializes the Linear Burdine relative permeability function 
  ! object

  implicit none
  
  class(rpf_BRAGFLO_KRP9_liq_type) :: this

  call RPFBaseInit(this)
  
end subroutine RPF_BRAGFLO_KRP9_Liq_Init

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP9_Liq_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_BRAGFLO_KRP9_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP9'
  endif  
  call RPFBaseVerify(this,string,option)
  
end subroutine RPF_BRAGFLO_KRP9_Liq_Verify

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP9_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! based on experimental measurements and analyses done by Vauclin et al.
  ! as discussed by Moridis and Pruess. 
  ! 14.	Moridis, G. J., and K. Pruess.  1992.  TOUGH Simulations of 
  ! Updegraff\92s Set of Fluid and Heat Flow Problems.  LBL-32611, ERMS# 138458. 
  ! Berkeley, CA:  Lawrence Berkeley Laboratory.
  ! Author: Heeho Park
  ! Date: 03/26/15
  ! 
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_BRAGFLO_KRP9_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dkr_dSe
  PetscReal :: dSe_dsat
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (1.d0-liquid_saturation)/(liquid_saturation)
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 0.d0
    return
  endif
  
  relative_permeability = 1.d0/(1.d0+28.768353d0*Se**1.7241379d0)
  ! Python analytical derivative (Jenn Frederick)
  dkr_dSe = -49.6006077278787*Se**0.7241379/(28.768353*Se**1.7241379+1.0)**2.0
  dSe_dsat = -1.0*(1.0/liquid_saturation) &
             - (1.0-liquid_saturation)/liquid_saturation**2.0
  dkr_sat = dkr_dSe * dSe_dsat
  
end subroutine RPF_BRAGFLO_KRP9_Liq_RelPerm
! End RPF: BRAGFLO KRP9 (Liquid)

! ************************************************************************** !

! Begin BRAGFLO KRP9 (Gas)
function RPF_BRAGFLO_KRP9_Gas_Create()

  ! Creates the Linear Burdine gas relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP9_gas_type), pointer :: RPF_BRAGFLO_KRP9_Gas_Create
  
  allocate(RPF_BRAGFLO_KRP9_Gas_Create)
  call RPF_BRAGFLO_KRP9_Gas_Create%Init()
  
end function RPF_BRAGFLO_KRP9_Gas_Create

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP9_Gas_Init(this)

  ! Initializes the Linear Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_BRAGFLO_KRP9_gas_type) :: this

  call RPFBaseInit(this)
  
end subroutine RPF_BRAGFLO_KRP9_Gas_Init

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP9_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_BRAGFLO_KRP9_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP9_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  
end subroutine RPF_BRAGFLO_KRP9_Gas_Verify

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP9_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! based on experimental measurements and analyses done by Vauclin et al.
  ! as discussed by Moridis and Pruess.  
  ! 14.	Moridis, G. J., and K. Pruess.  1992.  TOUGH Simulations of 
  ! Updegraff\92s Set of Fluid and Heat Flow Problems.  LBL-32611, ERMS# 138458. 
  ! Berkeley, CA:  Lawrence Berkeley Laboratory.
  ! Author: Heeho Park
  ! Date: 03/26/15
  ! 

  use Option_module
  
  implicit none

  class(rpf_BRAGFLO_KRP9_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: liquid_relative_permeability
  PetscReal :: liquid_dkr_sat

  relative_permeability = 0.d0
  dkr_sat = 0.d0

  Se = (1.d0-liquid_saturation)/(liquid_saturation)
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 1.d0
    return
  endif
  
  call RPF_BRAGFLO_KRP9_Liq_RelPerm(this,liquid_saturation, &
                        liquid_relative_permeability, &
                        liquid_dkr_sat,option)
  
  relative_permeability = 1.d0 - liquid_relative_permeability
  dkr_sat = -1.d0 * liquid_dkr_sat
  
end subroutine RPF_BRAGFLO_KRP9_Gas_RelPerm
! End RPF: BRAGFLO KRP9 (Gas)

! ************************************************************************** !

! Begin RPF: BRAGFLO KRP4 (Liq)
function RPF_BRAGFLO_KRP4_Liq_Create()

  ! Creates the KRP4 or BC_Burdine liq relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP4_liq_type), pointer :: RPF_BRAGFLO_KRP4_Liq_Create
  
  allocate(RPF_BRAGFLO_KRP4_Liq_Create)
  call RPF_BRAGFLO_KRP4_Liq_Create%Init()
  
end function RPF_BRAGFLO_KRP4_Liq_Create
! End RPF: BRAGFLO KRP4 (Liq)

! ************************************************************************** !

! Begin RPF: BRAGFLO KRP4 (Gas)
function RPF_BRAGFLO_KRP4_Gas_Create()

  ! Creates the KRP4 or BC_Burdine gas relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP4_gas_type), pointer :: RPF_BRAGFLO_KRP4_Gas_Create
  
  allocate(RPF_BRAGFLO_KRP4_Gas_Create)
  call RPF_BRAGFLO_KRP4_Gas_Create%Init()
  
end function RPF_BRAGFLO_KRP4_Gas_Create
! End RPF: BRAGFLO KRP4 (Gas)

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP4_Gas_Init(this)

  ! Initializes the Brooks-Corey Burdine gas relative permeability function 
  ! object

  implicit none
  
  class(rpf_BRAGFLO_KRP4_gas_type) :: this

  call RPFBaseInit(this)
  this%Srg = UNINITIALIZED_DOUBLE
  
end subroutine RPF_BRAGFLO_KRP4_Gas_Init

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP4_Gas_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_BRAGFLO_KRP4_gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP4_GAS'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_BRAGFLO_KRP4_Gas_Verify

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP4_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! (1) Chen, J., J.W. Hopmans, M.E. Grismer (1999) "Parameter estimation of
  !     of two-fluid capillary pressure-saturation and permeability functions",
  !     Advances in Water Resources, Vol. 22, No. 5, pp 479-493,
  !     http://dx.doi.org/10.1016/S0309-1708(98)00025-6.
  !   
  ! Author: Glenn Hammond
  ! Date: 12/11/07, 09/23/14
  ! 
  use Option_module
  
  implicit none

  class(rpf_BRAGFLO_KRP4_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: gas_saturation
  PetscReal :: dkr_Se
  PetscReal :: dSe_sat
  
  gas_saturation = 1.0d0 - liquid_saturation
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0  

  if (gas_saturation <= this%Srg) then
    relative_permeability = 0.0d0
  else if (liquid_saturation > this%Sr) then
    Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
    Seg = 1.d0 - Se
    ! reference #1
    relative_permeability = Seg*Seg*(1.d0-Se**(1.d0+2.d0/this%lambda))
    ! Mathematica analytical solution (Heeho Park)
    dkr_Se = -(1.d0+2.d0/this%lambda)*Seg**2.d0*Se**(2.d0/this%lambda) &
             - 2.d0*Seg*(1.d0-Se**(1.d0+2.d0/this%lambda))
    dSe_sat = 1.d0 / (1.d0 - this%Sr - this%Srg)
    dkr_sat = dkr_Se * dSe_sat
  else
    relative_permeability = 1.0d0
  endif
  
end subroutine RPF_BRAGFLO_KRP4_Gas_RelPerm
! End RPF: Burdine, Brooks-Corey (Gas)

! ************************************************************************** !
  
! Begin RPF: BRAGFLO KRP11 (Liquid)
function RPF_BRAGFLO_KRP11_Liq_Create()

  ! Creates the Linear Burdine relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP11_liq_type), pointer :: RPF_BRAGFLO_KRP11_Liq_Create
  
  allocate(RPF_BRAGFLO_KRP11_Liq_Create)
  call RPF_BRAGFLO_KRP11_Liq_Create%Init()
  
end function RPF_BRAGFLO_KRP11_Liq_Create

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP11_Liq_Init(this)

  ! Initializes the Linear Burdine relative permeability function 
  ! object

  implicit none
  
  class(rpf_BRAGFLO_KRP11_liq_type) :: this

  call RPFBaseInit(this)
  
end subroutine RPF_BRAGFLO_KRP11_Liq_Init

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP11_Liq_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_BRAGFLO_KRP11_liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,BRAGFLO_KRP11'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif
  if (Uninitialized(this%tolc)) then
    option%io_buffer = UninitializedMessage('TOLC',string)
    call printErrMsg(option)
  endif

end subroutine RPF_BRAGFLO_KRP11_Liq_Verify

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP11_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! KRP = 11 BRAGFLO relative permeability model
  ! the relative permeabilities decrease from 1 to zero linearly between
  ! the residual saturations (brine and gas) and the residual saturation
  ! plus a tolerance.
  !
  ! Author: Heeho Park
  ! Date: 03/26/15
  ! 
  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_BRAGFLO_KRP11_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: gas_saturation
  PetscReal :: tol
  
  gas_saturation = 1.d0 - liquid_saturation
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  tol = this%tolc * (1 - this%Sr - this%Srg)
  
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else if (gas_saturation <= this%Srg) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else if (liquid_saturation <= this%Sr+tol) then
    relative_permeability = (liquid_saturation - this%Sr)/tol
    dkr_sat = 1.d0/tol
  else if (gas_saturation <= this%Srg+tol) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  endif
    
end subroutine RPF_BRAGFLO_KRP11_Liq_RelPerm
! End RPF: BRAGFLO KRP11 (Liquid)

! ************************************************************************** !

! Begin BRAGFLO KRP11 (Gas)
function RPF_BRAGFLO_KRP11_Gas_Create()

  ! Creates the Linear Burdine gas relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP11_gas_type), pointer :: RPF_BRAGFLO_KRP11_Gas_Create
  
  allocate(RPF_BRAGFLO_KRP11_Gas_Create)
  call RPF_BRAGFLO_KRP11_Gas_Create%Init()
  
end function RPF_BRAGFLO_KRP11_Gas_Create

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP11_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! KRP = 11 BRAGFLO relative permeability model
  ! the relative permeabilities decrease from 1 to zero linearly between
  ! the residual saturations (brine and gas) and the residual saturation
  ! plus a tolerance.
  !
  ! Author: Heeho Park
  ! Date: 03/26/15
  ! 

  use Option_module
  
  implicit none

  class(rpf_BRAGFLO_KRP11_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: gas_saturation
  PetscReal :: tol
  
  gas_saturation = 1.d0 - liquid_saturation
  
  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  tol = this%tolc * (1 - this%Sr - this%Srg)
  
  if (liquid_saturation <= this%Sr) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else if (gas_saturation <= this%Srg) then
    relative_permeability = 0.d0
    dkr_sat = 0.d0
  else if (liquid_saturation <= this%Sr+tol) then
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  else if (gas_saturation <= this%Srg+tol) then
    relative_permeability = (gas_saturation - this%Srg)/tol
    dkr_sat = -1.d0/tol
  else
    relative_permeability = 1.d0
    dkr_sat = 0.d0
  endif
  
  end subroutine RPF_BRAGFLO_KRP11_Gas_RelPerm
! End RPF: BRAGFLO KRP11 (Gas)

! ************************************************************************** !

! Begin RPF: BRAGFLO KRP12 (Liq)
function RPF_BRAGFLO_KRP12_Liq_Create()

  ! Creates the KRP12 or BC_Burdine liq relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP12_liq_type), pointer :: RPF_BRAGFLO_KRP12_Liq_Create
  
  allocate(RPF_BRAGFLO_KRP12_Liq_Create)
  call RPF_BRAGFLO_KRP12_Liq_Create%Init()
  
end function RPF_BRAGFLO_KRP12_Liq_Create
! End RPF: BRAGFLO KRP12 (Liq)

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP12_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! Modified Brooks-Corey model KRP = 12 in BRAGFLO
  !   
  ! Author: Heeho Park
  ! Date: 11/13/15
  ! 
  use Option_module
  
  implicit none

  class(rpf_BRAGFLO_KRP12_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: power
  PetscReal :: dkr_dSe
  PetscReal :: dSe_dsat

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  
  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr)
  Se = max(min(Se,1.0d0),0.0d0)
  
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  else if (Se < this%Sr) then
    Se = this%Sr
  endif
  
  ! reference #1
  power = 3.d0+2.d0/this%lambda
  relative_permeability = Se**power
  dkr_dSe = power*relative_permeability/Se    
  dSe_dsat = 1.d0 / (1.d0 - this%Sr)
  dkr_sat = dkr_dSe * dSe_dsat
  
end subroutine RPF_BRAGFLO_KRP12_Liq_RelPerm
! End RPF: Burdine, Brooks-Corey (Liquid)
  
! ************************************************************************** !

! Begin RPF: BRAGFLO KRP12 (Gas)
function RPF_BRAGFLO_KRP12_Gas_Create()

  ! Creates the KRP12 or BC_Burdine gas relative permeability function object

  implicit none
  
  class(rpf_BRAGFLO_KRP12_gas_type), pointer :: RPF_BRAGFLO_KRP12_Gas_Create
  
  allocate(RPF_BRAGFLO_KRP12_Gas_Create)
  call RPF_BRAGFLO_KRP12_Gas_Create%Init()
  
end function RPF_BRAGFLO_KRP12_Gas_Create
! End RPF: BRAGFLO KRP12 (Gas)

! ************************************************************************** !

subroutine RPF_BRAGFLO_KRP12_Gas_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  ! 
  ! Modified Brooks-Corey model KRP = 12 in BRAGFLO
  !   
  ! Author:  Heeho Park
  ! Date: 11/13/15
  ! 
  use Option_module
  
  implicit none

  class(rpf_BRAGFLO_KRP12_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: Seg
  PetscReal :: gas_saturation
  PetscReal :: dkr_dSe
  PetscReal :: dSe_dsat
  
  gas_saturation = 1.0d0 - liquid_saturation
  relative_permeability = 0.d0
  dkr_sat = 0.d0

  if (gas_saturation <= this%Srg) then
    relative_permeability = 0.0d0
  else if (liquid_saturation > this%Sr) then
    Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sr - this%Srg)
    Se = max(min(Se,1.0d0),0.0d0)
    Seg = 1.d0 - Se
    ! reference #1
    relative_permeability = Seg*Seg*(1.d0-Se**(1.d0+2.d0/this%lambda))
    ! Mathematica analytical solution (Heeho Park)
    dkr_dSe = -(1.d0+2.d0/this%lambda)*Seg**2.d0*Se**(2.d0/this%lambda) &
              - 2.d0*Seg*(1.d0-Se**(1.d0+2.d0/this%lambda))
    dSe_dsat = 1.d0/(1.d0 - this%Sr - this%Srg)  
    dkr_sat = dkr_dSe * dSe_dsat
  else
    relative_permeability = 1.0d0
  endif
  
end subroutine RPF_BRAGFLO_KRP12_Gas_RelPerm
! End RPF: Burdine, Brooks-Corey (Gas)
  
! ************************************************************************** !

! Begin RPF: modified Kosugi (Liq)
function RPF_mK_Liq_Create()

  ! Creates the modified Kosugi liq relative permeability function object

  implicit none

  class(rpf_mK_liq_type), pointer :: RPF_mK_Liq_Create

  allocate(RPF_mK_Liq_Create)
  call RPF_mK_Liq_Create%Init()

end function RPF_mK_Liq_Create
! End RPF:  modified Kosugi (Liq)

! ************************************************************************** !

subroutine RPF_mK_Liq_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_mK_Liq_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'LIQUID_RELATIVE_PERM') > 0) then
    string = name
  else
    string = trim(name) // 'LIQUID_RELATIVE_PERM,MODIFIED_KOSUGI'
  endif
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%sigmaz)) then
    option%io_buffer = UninitializedMessage('SIGMAZ',string)
    call printErrMsg(option)
  endif

end subroutine RPF_mK_Liq_Verify
! ************************************************************************** !


subroutine RPF_mK_Liq_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  !
  ! Computes the relative permeability (and associated derivatives) as a
  ! function of saturation for modified Kosugi model
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498–502. http://dx.doi.org/10.1111/gwat.12220
  !
  ! Author: Kris Kuhlman
  ! Date: 2017
  !
  use Option_module
  use Utility_module

  implicit none

  ! gnu & intel extension and required in f2008
  intrinsic :: erfc

  PetscReal, parameter :: SQRT2 = sqrt(2.0d0)

  class(rpf_mK_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option

  PetscReal :: Se, dkr_Se
  PetscReal :: InvSatRange
  PetscReal :: erfcArg, erfcRes
  PetscReal :: invErfcRes
  PetscReal :: sqrtSe, expArg

  relative_permeability = 0.d0
  dkr_sat = 0.d0

  InvSatRange = 1.0d0/(1.0d0 - this%Sr)
  Se = (liquid_saturation - this%Sr)*InvSatRange
  if (Se >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Se <= 0.d0) then
    relative_permeability = 0.d0
    return
  endif

  invErfcRes = InverseNorm(Se)
  erfcArg = (this%sigmaz - invErfcRes)/SQRT2
  erfcRes = erfc(erfcArg)
  sqrtSe = sqrt(Se)
  relative_permeability = sqrtSe*erfcRes*5.0D-1

  ! from Wolfram Alpha (x -> Se)
  ! (InverseErfc[x] -> -1/Sqrt[x] InverseNorm[x/2])
  !
  ! D[(Sqrt[x] Erfc[sigmaz/Sqrt[2] + InverseErfc[2 x]])/2, x] =
  ! E^(InverseErfc[2 x]^2 - (simgaz/Sqrt[2] + InverseErfc[2 x])^2) * ...
  ! Sqrt[x] + Erfc[sigmaz/Sqrt[2] + InverseErfc[2 x]]/(4 Sqrt[x])
  expArg = 5.0D-1*invErfcRes**2 - erfcArg**2
  dkr_Se = erfcres/(4.0D0*sqrtSe) + sqrtSe*exp(expArg)

  ! InvSatRange = dSe/dsat
  dkr_sat = dkr_Se * InvSatRange 

end subroutine RPF_mK_Liq_RelPerm
! End RPF: modified Kosugi (Liquid)

! ************************************************************************** !

! Begin RPF: modified Kosugi (Gas)
function RPF_mK_Gas_Create()

  ! Creates the modified Kosugi gas relative permeability function object

  implicit none

  class(rpf_mK_gas_type), pointer :: RPF_mK_Gas_Create

  allocate(RPF_mK_Gas_Create)
  call RPF_mK_Gas_Create%Init()

end function RPF_mK_Gas_Create
! End RPF:  modified Kosugi (Gas)

! ************************************************************************** !

subroutine RPF_mK_Gas_Verify(this,name,option)

  use Option_module

  implicit none

  class(RPF_mK_Gas_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option

  character(len=MAXSTRINGLENGTH) :: string

  if (index(name,'GAS_RELATIVE_PERM') > 0) then
    string = name
  else
    string = trim(name) // 'GAS_RELATIVE_PERM,MODIFIED_KOSUGI'
  endif
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%sigmaz)) then
    option%io_buffer = UninitializedMessage('SIGMAZ',string)
    call printErrMsg(option)
  endif
  if (Uninitialized(this%srg)) then
    option%io_buffer = UninitializedMessage('SRG',string)
    call printErrMsg(option)
  endif

end subroutine RPF_mK_Gas_Verify
! ************************************************************************** !

subroutine RPF_mK_Gas_RelPerm(this,liquid_saturation, &
                              relative_permeability,dkr_sat,option)
  !
  ! Computes the relative permeability (and associated derivatives) as a
  ! function of saturation for modified Kosugi model
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498–502. http://dx.doi.org/10.1111/gwat.12220
  !
  ! Author: Kris Kuhlman
  ! Date: 2017
  !
  use Option_module
  use Utility_module, only : InverseNorm

  implicit none

  ! gnu & intel extension and required in f2008
  intrinsic :: erfc

  PetscReal, parameter :: SQRT2 = sqrt(2.0d0)

  class(rpf_mK_gas_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option

  PetscReal :: Se, Seg, InvSatRange
  PetscReal :: dkr_Se, dSe_sat
  PetscReal :: erfcArg, erfcRes
  PetscReal :: invErfcRes
  PetscReal :: sqrtSe, expArg

  InvSatRange = 1.d0/(1.d0 - this%Sr - this%Srg)
  Se = (liquid_saturation - this%Sr)*InvSatRange

  relative_permeability = 0.d0
  dkr_sat = 0.d0
  if (Se >= 1.d0) then
    relative_permeability = 0.d0
    return
  else if (Se <=  0.d0) then
    relative_permeability = 1.d0
    return
  endif

  Seg = 1.d0 - Se

  invErfcRes = InverseNorm(Seg)
  erfcArg = (this%sigmaz - invErfcRes)/SQRT2
  erfcRes = erfc(erfcArg)
  sqrtSe = sqrt(Seg)
  relative_permeability = sqrtSe*erfcRes*5.0D-1

  ! from Wolfram Alpha (x -> Seg)
  ! (InverseErfc[x] -> -1/Sqrt[x] InverseNorm[x/2])
  !
  ! D[(Sqrt[x] Erfc[sigmaz/Sqrt[2] + InverseErfc[2 x]])/2, x] =
  ! E^(InverseErfc[2 x]^2 - (simgaz/Sqrt[2] + InverseErfc[2 x])^2) * ...
  ! Sqrt[x] + Erfc[sigmaz/Sqrt[2] + InverseErfc[2 x]]/(4 Sqrt[x])
  expArg = 5.0D-1*invErfcRes**2 - erfcArg**2
  dkr_Se = erfcres/(4.0D0*sqrtSe) + sqrtSe*exp(expArg)

  ! -1 = dSeg/dSe
  ! InvSatRange = dSe/dsat
  dkr_sat = -1.d0 * dkr_Se * InvSatRange 

end subroutine RPF_MK_Gas_RelPerm
! End RPF: modified Kosigi (Gas)

! ************************************************************************** !

! Begin RPF: TOUGH2, Linear (Oil) 
function RPF_TOUGH2_Linear_Oil_Create()

  ! Creates the TOUGH2 Linear oil relative permeability function object
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/19/2015

  class(rpf_TOUGH2_Linear_oil_type), pointer :: RPF_TOUGH2_Linear_Oil_Create

  allocate(RPF_TOUGH2_Linear_Oil_Create)
  call RPF_TOUGH2_Linear_Oil_Create%Init()

end function RPF_TOUGH2_Linear_Oil_Create

! ************************************************************************** !

subroutine RPF_TOUGH2_Linear_Oil_Init(this)

  ! Initializes the TOUGH2 Linear Oil relative permeability function 
  ! object
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/19/2015

  implicit none
  
  class(rpf_TOUGH2_Linear_oil_type) :: this

  call RPFBaseInit(this)
  this%Sro = UNINITIALIZED_DOUBLE
  
end subroutine RPF_TOUGH2_Linear_Oil_Init

! ************************************************************************** !

subroutine RPF_TOUGH2_Linear_Oil_Verify(this,name,option)

  use Option_module

  implicit none
  
  class(rpf_TOUGH2_Linear_oil_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,TOUGH2_LINEAR_OIL'
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Sro)) then
    option%io_buffer = UninitializedMessage('OIL_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  
  
end subroutine RPF_TOUGH2_Linear_Oil_Verify

! ************************************************************************** !

subroutine RPF_TOUGH2_Linear_Oil_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 10/19/2015


  use Option_module
  
  implicit none

  class(rpf_TOUGH2_Linear_oil_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: So
  PetscReal :: Seo
  
  ! initialize to derivative to NaN so that not mistakenly used.
  dkr_sat = 0.d0
  dkr_sat = dkr_sat / 0.d0
  dkr_sat = dkr_sat * 0.d0

  So = 1.d0 - liquid_saturation

  Seo = (So - this%Sro) / (1.d0 - this%Sro)

  if (Seo >= 1.d0) then
    relative_permeability = 1.d0
    return
  else if (Seo <=  0.d0) then
    relative_permeability = 0.d0
    return
  endif

  relative_permeability = Seo

end subroutine RPF_TOUGH2_Linear_Oil_RelPerm
! End RPF: TOUGH2, Linear (Oil)

! ************************************************************************** !

!Beginning RPF Modified Brooks-Corey for liq and oil phase (RPF_Mod_BC_Oil)

!  procedure, public :: Init => RPF_Mod_BC_Oil_Init 
!  procedure, public :: Verify => RPF_Mod_BC_Oil_Verify
!  procedure, public :: SetupPolynomials => RPF_Mod_BC_SetupPolynomials
!  procedure, public :: RelativePermeability => RPF_Mod_BC_Oil_RelPerm

function RPF_Mod_BC_Liq_Create()

  ! Creates the Modified BC Oil relative permeability function object
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/20/2016

  class(rpf_mod_BC_liq_type), pointer :: RPF_Mod_BC_Liq_Create

  allocate(RPF_Mod_BC_Liq_Create)
  call RPF_Mod_BC_Liq_Create%Init()

end function RPF_Mod_BC_Liq_Create

! ************************************************************************** !

function RPF_Mod_BC_Oil_Create()

  ! Creates the Modified BC Oil relative permeability function object
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/20/2016

  class(rpf_mod_BC_oil_type), pointer :: RPF_Mod_BC_Oil_Create

  allocate(RPF_Mod_BC_Oil_Create)
  call RPF_Mod_BC_Oil_Create%Init()

end function RPF_Mod_BC_Oil_Create

! ************************************************************************** !

!subroutine RPF_Mod_BC_Oil_Init(this)
subroutine RPF_Mod_BC_Init(this)

  ! Initializes the Modified BC Oil relative permeability function object 
  ! object
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/20/2016

  implicit none
  
  !class(rpf_mod_BC_oil_type) :: this
  class(rpf_mod_BC_type) :: this

  call RPFBaseInit(this)
  this%m = UNINITIALIZED_DOUBLE
  this%Srg = UNINITIALIZED_DOUBLE
  this%Sro = UNINITIALIZED_DOUBLE
  this%kr_max = 1.0d0
   
!end subroutine RPF_Mod_BC_Oil_Init
end subroutine RPF_Mod_BC_Init

! ************************************************************************** !

!subroutine RPF_Mod_BC_Oil_Verify(this,name,option)
subroutine RPF_Mod_BC_Verify(this,name,option)

  use Option_module

  implicit none
  
  !class(rpf_mod_BC_oil_type) :: this
  class(rpf_mod_BC_type) :: this
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option
  
  character(len=MAXSTRINGLENGTH) :: string 

  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    select type(rpf => this)
      class is(rpf_mod_BC_liq_type) 
        string = trim(name) // 'PERMEABILITY_FUNCTION,MOD_BC_LIQ'
      class is(rpf_mod_BC_oil_type)
        string = trim(name) // 'PERMEABILITY_FUNCTION,MOD_BC_OIL'
    end select
  endif    
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%Sro)) then
    option%io_buffer = UninitializedMessage('OIL_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif
  if (Uninitialized(this%Srg)) then
    option%io_buffer = UninitializedMessage('GAS_RESIDUAL_SATURATION',string)
    call printErrMsg(option)
  endif  

  if (Uninitialized(this%m)) then
    option%io_buffer = UninitializedMessage('POWER EXPONENT',string)
    call printErrMsg(option)
  endif
  
!end subroutine RPF_Mod_BC_Oil_Verify
end subroutine RPF_Mod_BC_Verify


! ************************************************************************** !

subroutine RPF_Mod_BC_SetupPolynomials(this,option,error_string)

  ! Sets up polynomials for smoothing Modified BC permeability function

  use Option_module
  use Utility_module
  
  implicit none
  
  class(rpf_mod_BC_type) :: this
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: error_string
  
  PetscReal :: b(4)

  PetscReal :: Se_ph_low

  if (associated(this%poly)) call PolynomialDestroy(this%poly)
  this%poly => PolynomialCreate()
  ! fill matix with values
  this%poly%low = 0.99d0  ! just below saturated
  !this%poly%low = 0.95d0  ! just below saturated 
  this%poly%high = 1.d0   ! saturated
  Se_ph_low = this%poly%low
  !select type(rpf => this)
  !  class is(rpf_mod_BC_liq_type) 
  !    Se_ph_low = ( this%poly%low - this%Sr ) / &
  !                (1.0 - this%Sro - this%Sr - this%Srg)
  !  class is(rpf_mod_BC_oil_type)
  !    Se_ph_low = ( this%poly%low - this%Sro ) / &
  !                (1.0 - this%Sro - this%Sr - this%Srg) 
  !end select 

  b(1) = this%kr_max
  b(2) = this%kr_max * (Se_ph_low ** this%m)
  b(3) = 0.d0
  b(4) = this%m * this%kr_max * Se_ph_low ** (this%m - 1.0 )
  
  call CubicPolynomialSetup(this%poly%high,this%poly%low,b)
  
  this%poly%coefficients(1:4) = b(1:4)
  
end subroutine RPF_Mod_BC_SetupPolynomials

! ************************************************************************** !

subroutine RPF_Mod_BC_Liq_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/21/2016

  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Mod_BC_liq_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dkr_Se
  
  ! initialize to derivative to NaN so that not mistakenly used.
  dkr_sat = 0.d0
  dkr_sat = dkr_sat / 0.d0
  dkr_sat = dkr_sat * 0.d0

  Se = (liquid_saturation - this%Sr) / (1.d0 - this%Sro - this%Sr - this%Srg )

  if (Se >= 1.d0) then
    relative_permeability = this%kr_max
    return
  else if (Se <=  0.d0) then
    relative_permeability = 0.d0
    return
  endif

  if (associated(this%poly)) then
    if (Se > this%poly%low) then
      call CubicPolynomialEvaluate(this%poly%coefficients, &
                                   Se,relative_permeability,dkr_Se)
      return
    endif
  endif

  relative_permeability = this%kr_max * (Se ** this%m)

end subroutine RPF_Mod_BC_Liq_RelPerm

! ************************************************************************** !


subroutine RPF_Mod_BC_Oil_RelPerm(this,liquid_saturation, &
                                     relative_permeability,dkr_sat,option)
  ! 
  ! Computes the relative permeability (and associated derivatives) as a 
  ! function of saturation
  !
  ! Author: Paolo Orsini (OGS)
  ! Date: 02/20/2016

  use Option_module
  use Utility_module
  
  implicit none

  class(rpf_Mod_BC_oil_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  PetscReal :: So
  PetscReal :: Seo
  PetscReal :: dkr_Se
  
  ! initialize to derivative to NaN so that not mistakenly used.
  dkr_sat = 0.d0
  dkr_sat = dkr_sat / 0.d0
  dkr_sat = dkr_sat * 0.d0

  So = 1.d0 - liquid_saturation

  Seo = (So - this%Sro) / (1.d0 - this%Sro - this%Sr - this%Srg ) 

  if (Seo >= 1.d0) then
    relative_permeability = this%kr_max
    return
  else if (Seo <=  0.d0) then
    relative_permeability = 0.d0
    return
  endif

  if (associated(this%poly)) then
    if (Seo > this%poly%low) then
      call CubicPolynomialEvaluate(this%poly%coefficients, &
                                   Seo,relative_permeability,dkr_Se)
      return
    endif
  endif

  relative_permeability = this%kr_max * (Seo ** this%m)

end subroutine RPF_Mod_BC_Oil_RelPerm

!End RPF: Modified Brooks-Corey for the oil phase (RPF_Mod_BC_Oil)

! ************************************************************************** !

function RPF_Constant_Create()

  ! Creates the constant relative permeability function object

  implicit none
  
  class(rel_perm_func_constant_type), pointer :: RPF_Constant_Create
  
  allocate(RPF_Constant_Create)
  call RPFBaseInit(RPF_Constant_Create)
  ! set Sr = 0. to avoid uninitialized failure
  RPF_Constant_Create%Sr = 0.d0
  RPF_Constant_Create%kr = 0.d0
  
end function RPF_Constant_Create

! ************************************************************************** !

subroutine RPFConstantVerify(this,name,option)

  use Option_module
  
  implicit none
  
  class(rel_perm_func_constant_type) :: this  
  character(len=MAXSTRINGLENGTH) :: name
  type(option_type) :: option  

  character(len=MAXSTRINGLENGTH) :: string
  
  if (index(name,'PERMEABILITY_FUNCTION') > 0) then
    string = name
  else
    string = trim(name) // 'PERMEABILITY_FUNCTION,CONSTANT'
  endif  
  call RPFBaseVerify(this,string,option)
  if (Uninitialized(this%kr)) then
    option%io_buffer = UninitializedMessage('RELATIVE_PERMEABILITY',string)
    call printErrMsg(option)
  endif   

end subroutine RPFConstantVerify

! ************************************************************************** !

subroutine RPF_ConstantRelPerm(this,liquid_saturation,relative_permeability, &
                            dkr_sat,option)
  use Option_module

  implicit none
  
  class(rel_perm_func_constant_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: relative_permeability
  PetscReal, intent(out) :: dkr_sat
  type(option_type), intent(inout) :: option
  
  relative_permeability = this%kr
  dkr_sat = 0.d0
  
end subroutine RPF_ConstantRelPerm

! ************************************************************************** !

subroutine PolynomialDestroy(poly)
  ! 
  ! Destroys a polynomial smoother
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  type(polynomial_type), pointer :: poly
  
  if (.not.associated(poly)) return
  
  deallocate(poly)
  nullify(poly)

end subroutine PolynomialDestroy

! ************************************************************************** !

subroutine SaturationFunctionDestroy(sf)
  ! 
  ! Destroys a saturuation function
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(sat_func_base_type), pointer :: sf
  
  if (.not.associated(sf)) return
  
  call PolynomialDestroy(sf%sat_poly)
  call PolynomialDestroy(sf%pres_poly)
#ifdef SMOOTHING2
  call PolynomialDestroy(sf%sat_poly2)
  call PolynomialDestroy(sf%pres_poly2)
#endif
  deallocate(sf)
  nullify(sf)

end subroutine SaturationFunctionDestroy

! ************************************************************************** !

subroutine PermeabilityFunctionDestroy(rpf)
  ! 
  ! Destroys a saturuation function
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(rel_perm_func_base_type), pointer :: rpf
  
  if (.not.associated(rpf)) return
  
  call PolynomialDestroy(rpf%poly)
#ifdef SMOOTHING2
  call PolynomialDestroy(rpf%poly2)
#endif
  deallocate(rpf)
  nullify(rpf)

end subroutine PermeabilityFunctionDestroy

! ************************************************************************** !

recursive subroutine CharacteristicCurvesDestroy(cc)
  ! 
  ! Destroys a characteristic curve
  ! 
  ! Author: Glenn Hammond
  ! Date: 09/24/14
  ! 

  implicit none
  
  class(characteristic_curves_type), pointer :: cc
  
  if (.not.associated(cc)) return
  
  call CharacteristicCurvesDestroy(cc%next)
  
  call SaturationFunctionDestroy(cc%saturation_function)

  ! the liquid and gas relative permeability pointers may pointer to the
  ! same address. if so, destroy one and nullify the other.
  if (associated(cc%liq_rel_perm_function,cc%gas_rel_perm_function)) then
    call PermeabilityFunctionDestroy(cc%liq_rel_perm_function)
    nullify(cc%gas_rel_perm_function)
  !PO how about avoiding xxx_rel_perm_function => aaa_rel_perm_function? 
  !   it should semplify code. It seems we do this only to pass verify 
  else if (associated(cc%oil_rel_perm_function,cc%gas_rel_perm_function)) then 
    call PermeabilityFunctionDestroy(cc%oil_rel_perm_function)
    nullify(cc%gas_rel_perm_function)
  else
    call PermeabilityFunctionDestroy(cc%liq_rel_perm_function)
    call PermeabilityFunctionDestroy(cc%gas_rel_perm_function)
    !call PermeabilityFunctionDestroy(cc%oil_rel_perm_function)
  endif

  deallocate(cc)
  nullify(cc)
  
end subroutine CharacteristicCurvesDestroy

! ************************************************************************** !

end module Characteristic_Curves_module
