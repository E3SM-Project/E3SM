module Characteristic_Curves_Common_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Characteristic_Curves_Base_module

  implicit none

  private
  
!-----------------------------------------------------------------------------
!-- Saturation Functions -----------------------------------------------------
!-----------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_default_type
  contains
    procedure, public :: Verify => SFDefaultVerify
    procedure, public :: CapillaryPressure => SFDefaultCapillaryPressure
    procedure, public :: Saturation => SFDefaultSaturation
  end type sat_func_default_type  
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_constant_type
    PetscReal :: constant_capillary_pressure
    PetscReal :: constant_saturation
  contains
    procedure, public :: Verify => SFConstantVerify
    procedure, public :: CapillaryPressure => SFConstantCapillaryPressure
    procedure, public :: Saturation => SFConstantSaturation
  end type sat_func_constant_type
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_VG_type
    PetscReal :: alpha
    PetscReal :: m
  contains
    procedure, public :: Init => SF_VG_Init
    procedure, public :: Verify => SF_VG_Verify
    procedure, public :: CapillaryPressure => SF_VG_CapillaryPressure
    procedure, public :: Saturation => SF_VG_Saturation
  end type sat_func_VG_type  
  !---------------------------------------------------------------------------
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
  !---------------------------------------------------------------------------
  type, public, extends(sat_func_base_type) :: sat_func_Linear_type
    PetscReal :: alpha
  contains
    procedure, public :: Init => SF_Linear_Init
    procedure, public :: Verify => SF_Linear_Verify
    procedure, public :: CapillaryPressure => SF_Linear_CapillaryPressure
    procedure, public :: Saturation => SF_Linear_Saturation
  end type sat_func_Linear_type
  type, public, extends(sat_func_base_type) :: sat_func_mK_type
    PetscReal :: sigmaz, muz
    PetscReal :: rmax, r0
    PetscInt :: nparam
  contains
    procedure, public :: Init => SF_mK_Init
    procedure, public :: Verify => SF_mK_Verify
    procedure, public :: CapillaryPressure => SF_mK_CapillaryPressure
    procedure, public :: Saturation => SF_mK_Saturation
  end type sat_func_mK_type
  
!-----------------------------------------------------------------------------
!-- Relative Permeability Functions ------------------------------------------
!-----------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rel_perm_func_default_type
  contains
    procedure, public :: Verify => RPFDefaultVerify
    procedure, public :: RelativePermeability => RPF_DefaultRelPerm
  end type rel_perm_func_default_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_VG_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_Mualem_VG_Liq_Init
    procedure, public :: Verify => RPF_Mualem_VG_Liq_Verify
    procedure, public :: SetupPolynomials => RPF_Mualem_SetupPolynomials
    procedure, public :: RelativePermeability => RPF_Mualem_VG_Liq_RelPerm
  end type rpf_Mualem_VG_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_VG_gas_type
    PetscReal :: m
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Mualem_VG_Gas_Init
    procedure, public :: Verify => RPF_Mualem_VG_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_VG_Gas_RelPerm
  end type rpf_Mualem_VG_gas_type

  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_BC_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_Burdine_BC_Liq_Init
    procedure, public :: Verify => RPF_Burdine_BC_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_BC_Liq_RelPerm
    procedure, public :: SetupPolynomials => RPF_Burdine_BC_Liq_SetupPolynomials  ! added by F.-M. Yuan (2017-03-10)
  end type rpf_Burdine_BC_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_BC_gas_type
    PetscReal :: lambda
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Burdine_BC_Gas_Init
    procedure, public :: Verify => RPF_Burdine_BC_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_BC_Gas_RelPerm
  end type rpf_Burdine_BC_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_BC_liq_type
    PetscReal :: lambda
  contains
    procedure, public :: Init => RPF_Mualem_BC_Liq_Init
    procedure, public :: Verify => RPF_Mualem_BC_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_BC_Liq_RelPerm
  end type rpf_MUALEM_BC_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_BC_gas_type
    PetscReal :: lambda
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Mualem_BC_Gas_Init
    procedure, public :: Verify => RPF_Mualem_BC_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_BC_Gas_RelPerm
  end type rpf_Mualem_BC_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_VG_liq_type
    PetscReal :: m
  contains
    procedure, public :: Init => RPF_Burdine_VG_Liq_Init
    procedure, public :: Verify => RPF_Burdine_VG_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_VG_Liq_RelPerm
  end type rpf_Burdine_VG_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_VG_gas_type
    PetscReal :: m
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Burdine_VG_Gas_Init
    procedure, public :: Verify => RPF_Burdine_VG_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_VG_Gas_RelPerm
  end type rpf_Burdine_VG_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Mualem_Linear_liq_type
    PetscReal :: pcmax
    PetscReal :: alpha
  contains
    procedure, public :: Init => RPF_Mualem_Linear_Liq_Init
    procedure, public :: Verify => RPF_Mualem_Linear_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_Linear_Liq_RelPerm
  end type rpf_Mualem_Linear_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rpf_Mualem_Linear_liq_type) :: & 
                        rpf_Mualem_Linear_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Mualem_Linear_Gas_Init
    procedure, public :: Verify => RPF_Mualem_Linear_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Mualem_Linear_Gas_RelPerm
  end type rpf_Mualem_Linear_gas_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_Burdine_Linear_liq_type
  contains
    procedure, public :: Init => RPF_Burdine_Linear_Liq_Init
    procedure, public :: Verify => RPF_Burdine_Linear_Liq_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_Linear_Liq_RelPerm
  end type rpf_Burdine_Linear_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: & 
                        rpf_Burdine_Linear_gas_type
    PetscReal :: Srg
  contains
    procedure, public :: Init => RPF_Burdine_Linear_Gas_Init
    procedure, public :: Verify => RPF_Burdine_Linear_Gas_Verify
    procedure, public :: RelativePermeability => RPF_Burdine_Linear_Gas_RelPerm
  end type rpf_Burdine_Linear_gas_type  
  !---------------------------------------------------------------------------
  ! Constant: for running tests with a fixed relative permeability
  type, public, extends(rel_perm_func_base_type) :: rel_perm_func_constant_type
    PetscReal :: kr
  contains
    procedure, public :: Verify => RPFConstantVerify
    procedure, public :: RelativePermeability => RPF_ConstantRelPerm
  end type rel_perm_func_constant_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_mK_liq_type
    PetscReal :: sigmaz
  contains
    procedure, public :: Init => RPF_mK_Liq_Init
    procedure, public :: Verify => RPF_mK_Liq_Verify
    procedure, public :: RelativePermeability => RPF_mK_Liq_RelPerm
  end type rpf_mK_liq_type
  !---------------------------------------------------------------------------
  type, public, extends(rel_perm_func_base_type) :: rpf_mK_gas_type
    PetscReal :: Srg
    PetscReal :: sigmaz
  contains
    procedure, public :: Verify => RPF_mK_Gas_Verify
    procedure, public :: RelativePermeability => RPF_mK_Gas_RelPerm
  end type rpf_mK_gas_type
  
  public :: &! standard char. curves:
            SF_Default_Create, &
            SF_Constant_Create, &
            SF_VG_Create, &
            SF_BC_Create, &
            SF_Linear_Create, &
            SF_mK_Create, &
            ! standard rel. perm. curves:
            RPF_Default_Create, &
            RPF_Constant_Create, &  
            RPF_Mualem_VG_Liq_Create, &
            RPF_Mualem_VG_Gas_Create, &
            RPF_Burdine_BC_Liq_Create, &
            RPF_Burdine_BC_Gas_Create, &
            RPF_Mualem_BC_Liq_Create, &
            RPF_Mualem_BC_Gas_Create, &
            RPF_Burdine_VG_Liq_Create, &
            RPF_Burdine_VG_Gas_Create, &
            RPF_Mualem_Linear_Liq_Create, &
            RPF_Mualem_Linear_Gas_Create, &
            RPF_Burdine_Linear_Liq_Create, &
            RPF_Burdine_Linear_Gas_Create, &
            RPF_mK_Liq_Create, &
            RPF_mK_Gas_Create, &
            RPF_Mualem_VG_Liq_RelPerm
  
contains

! ************************************************************************** !
! ************************************************************************** !

function SF_Default_Create()

  ! Creates the default saturation function object

  implicit none
  
  class(sat_func_default_type), pointer :: SF_Default_Create
  
  allocate(SF_Default_Create)
  call SFBaseInit(SF_Default_Create)
  SF_Default_Create%Sr = 0.d0
  
  SF_Default_Create%analytical_derivative_available = PETSC_TRUE
  
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
                                      capillary_pressure,dpc_dsatl,option)
  use Option_module
  
  implicit none
  
  class(sat_func_default_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  if (liquid_saturation < 1.d0) then
    option%io_buffer = 'SFDefaultCapillaryPressure is a dummy routine used &
      &for saturated flow only.  The user must specify a valid &
      &SATURATION_FUNCTION.'
    call printErrMsgByRank(option)
  endif

end subroutine SFDefaultCapillaryPressure

! ************************************************************************** !

subroutine SFDefaultSaturation(this,capillary_pressure, &
                               liquid_saturation,dsat_dpres,option)
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  SF_Constant_Create%analytical_derivative_available = PETSC_TRUE
  
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
    case(TH_MODE)
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
                                       capillary_pressure,dpc_dsatl,option)
  use Option_module
  
  implicit none
  
  class(sat_func_constant_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  dpc_dsatl = 0.d0
  capillary_pressure = this%constant_capillary_pressure

end subroutine SFConstantCapillaryPressure

! ************************************************************************** !

subroutine SFConstantSaturation(this,capillary_pressure, &
                                liquid_saturation,dsat_dpres,option)
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

! ************************************************************************** !
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

! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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
                                   capillary_pressure,dpc_dsatl,option)
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
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: n
  PetscReal :: Se

  PetscReal :: neg_one_over_m
  PetscReal :: one_over_n
  PetscReal :: dSe_dsatl
  PetscReal :: Se_sup_neg_one_over_m
  PetscReal :: Se_sup_neg_one_over_m_minus_one
  
  dpc_dsatl = 0.d0

  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  n = 1.d0/(1.d0-this%m)
  neg_one_over_m = -1.d0/this%m
  one_over_n = 1.d0/n
  dSe_dsatl = 1.d0 / (1.d0-this%Sr)
  Se = (liquid_saturation-this%Sr)*dSe_dsatl
  Se_sup_neg_one_over_m = Se**neg_one_over_m
  Se_sup_neg_one_over_m_minus_one = Se_sup_neg_one_over_m - 1.d0
  capillary_pressure = (Se_sup_neg_one_over_m_minus_one**one_over_n)/this%alpha
  dpc_dsatl = capillary_pressure/Se_sup_neg_one_over_m_minus_one * &
              one_over_n * neg_one_over_m * Se_sup_neg_one_over_m / Se * &
              dSe_dsatl

#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
    dpc_dsatl = dpc_dsatl*(1.d0-liquid_saturation)/0.001d0 +
                capillary_pressure*(-1.d0)/0.001d0
  endif
#endif

  if (capillary_pressure > this%pcmax) then
    capillary_pressure = this%pcmax
    dpc_dsatl = 0.d0
  endif
  
end subroutine SF_VG_CapillaryPressure

! ************************************************************************** !

subroutine SF_VG_Saturation(this,capillary_pressure, &
                            liquid_saturation,dsat_dpres,option)
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

! ************************************************************************** !
! ************************************************************************** !

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

  call SFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
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
          (0.001d0**(-1.d0/this%lambda))/this%alpha)       ! Se = 0.001 (note that if Se=0, pc=0 by this eq.)
  this%pres_poly2%low = min(this%pres_poly2%high-1.0d2, &  ! this is the starting point for smoothing, so it must be less than %high end.
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
                                   capillary_pressure,dpc_dsatl,option)
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
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
#ifdef SMOOTHING2
  PetscReal :: Hfunc, dHfunc, x1, x0, y
#endif
  PetscReal :: dSe_dsatl
  PetscReal :: dpc_dSe
  PetscReal :: neg_one_over_lambda

  dpc_dsatl = 0.d0
  
  if (liquid_saturation <= this%Sr) then   ! F.-M. Yuan (2017-03-14): This IS not mathmatically correct, although this function seems NOT called except for unit test
#ifndef SMOOTHING2
    capillary_pressure = this%pcmax
    return
#endif
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  dSe_dsatl = 1.d0 / (1.d0-this%Sr)
  Se = (liquid_saturation-this%Sr)*dSe_dsatl
#ifndef SMOOTHING2
  if (associated(this%sat_poly)) then
    if (Se > this%sat_poly%low) then
      call QuadraticPolynomialEvaluate(this%sat_poly%coefficients(1:3), &
                                       Se,capillary_pressure,dpc_dSe)
      dpc_dsatl = dpc_dSe*dSe_dsatl
      return
    endif
  endif
#endif
  neg_one_over_lambda = -1.d0/this%lambda
  capillary_pressure = (Se**neg_one_over_lambda)/this%alpha
  dpc_dsatl = capillary_pressure/Se*neg_one_over_lambda*dSe_dsatl

#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
    dpc_dsatl = dpc_satl*(1.d0-liquid_saturation)/0.001d0 + &
                capillary_pressure*(-1.d0/0.001d0)
  endif
#endif  

#ifndef SMOOTHING2
  if (capillary_pressure > this%pcmax) then
    capillary_pressure = this%pcmax
    dpc_dsatl = 0.d0
  endif

#else
  ! WET-end
  if (associated(this%pres_poly)) then
    call this%Saturation(this%pres_poly%high,x1,y,option)
    x1 = this%Sr+(1.d0-this%Sr)*x1           ! S to start smoothing, towards saturated (S=1)
    if (liquid_saturation>x1) then
      call this%Saturation(this%pres_poly%low,x0,y,option)
      x0 = this%Sr+(1.d0-this%Sr)*x0         ! S to end smoothing @pc=pres_poly%low ( should be 0, i.e. S=1)

      call HFunctionSmooth(liquid_saturation, x1, x0, Hfunc, dHfunc)

      dpc_dsatl = (capillary_pressure-this%pres_poly%low)*dHfunc + dpc_dsatl*Hfunc
      capillary_pressure = capillary_pressure*Hfunc + this%pres_poly%low * (1.d0-Hfunc)
      if(capillary_pressure<this%pres_poly%low) then
        capillary_pressure=this%pres_poly%low
        dpc_dsatl = 0.d0
      end if
    end if
  end if

  ! DRY-end
  ! fmy: by BC function, mathemaatically @Sr, pc not necessarily equals to pcmax
  !      So, @pcmax, S is not necessarily equaled to Sr, and requiring smoothing
  if (associated(this%pres_poly2)) then
    call this%Saturation(this%pres_poly2%low,x1,y,option)
    x1 = this%Sr+(1.d0-this%Sr)*x1         ! S to start smoothing, towards S=Sr
    if (liquid_saturation<x1) then
      call this%Saturation(this%pres_poly2%high,x0,y,option)
      x0 = this%Sr+(1.d0-this%Sr)*x0         ! S to end smoothing @pc=pres_poly2%high (note: if Se=0, mathmatically pc is NaN or inf)

      call HFunctionSmooth(liquid_saturation, x1, x0, Hfunc, dHfunc)

      dpc_dsatl = (capillary_pressure-this%pres_poly2%high)*dHfunc + dpc_dsatl*Hfunc
      capillary_pressure = capillary_pressure*Hfunc + this%pres_poly2%high * (1.d0-Hfunc)
      if(capillary_pressure>this%pres_poly2%high) then
        capillary_pressure=this%pres_poly2%high
        dpc_dsatl = 0.d0
      end if

    end if
  end if

#endif
  
end subroutine SF_BC_CapillaryPressure

! ************************************************************************** !

subroutine SF_BC_Saturation(this,capillary_pressure, &
                            liquid_saturation,dsat_dpres,option)
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
    if (capillary_pressure<this%pres_poly%high) then
      call HFunctionSmooth(capillary_pressure, this%pres_poly%high, this%pres_poly%low, &
                           Hfunc, dHfunc)
      !sat = 1.0 + (sat-1.0) * Hfunc, in which '1.0' is the base sat @pc=low end
      dsat_dpres = dsat_dpres*Hfunc + (liquid_saturation-1.d0) * dHfunc
      liquid_saturation = 1.0d0 + (liquid_saturation-1.d0)*Hfunc
      if(liquid_saturation>1.0d0) then
        liquid_saturation=1.0d0
        dsat_dpres = 0.d0
      endif

    endif
  endif

  ! DRY-end
  if (associated(this%pres_poly2)) then
    if (capillary_pressure>this%pres_poly2%low) then
      call HFunctionSmooth(capillary_pressure, this%pres_poly2%low, this%pres_poly2%high, &
                           Hfunc, dHfunc)
      !sat = sr + (sat-sr) * Hfunc, in which 'sr' is the base sat @pc=high (high end)
      dsat_dpres = dsat_dpres * Hfunc + (liquid_saturation-this%Sr) * dHfunc
      liquid_saturation = this%Sr+(liquid_saturation-this%Sr) * Hfunc
      if(liquid_saturation<this%Sr) then
        liquid_saturation=this%Sr
        dsat_dpres = 0.d0
      endif

    endif
  endif

#endif

end subroutine SF_BC_Saturation

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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
                                       capillary_pressure,dpc_dsatl,option)
  ! 
  ! Computes the capillary pressure as a function of saturation.
  !
  ! Author: Bwalya Malama, Heeho Park
  ! Date: 11/14/14
  !
  use Option_module
  use Material_Aux_class
  
  implicit none
  
  class(sat_func_Linear_type) :: this
  PetscReal, intent(in) :: liquid_saturation
  PetscReal, intent(out) :: capillary_pressure
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option
  
  PetscReal :: Se
  PetscReal :: dSe_dsatl
  PetscReal :: one_over_alpha_minus_pcmax

  dpc_dsatl = 0.d0
  
  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif
  
  dSe_dsatl = 1.d0/(1.d0-this%Sr)
  Se = (liquid_saturation-this%Sr)*dSe_dsatl
  one_over_alpha_minus_pcmax = 1.d0/this%alpha-this%pcmax
  capillary_pressure = one_over_alpha_minus_pcmax*Se + this%pcmax
  dpc_dsatl = one_over_alpha_minus_pcmax*dSe_dsatl

#if defined(MATCH_TOUGH2)
  if (liquid_saturation > 0.999d0) then
    capillary_pressure = capillary_pressure*(1.d0-liquid_saturation)/0.001d0
    dpc_dsatl = dpc_satl*(1.d0-liquid_saturation)/0.001d0 + &
                capillary_pressure*(-1.d0/0.001d0)
  endif
#endif  

  if (capillary_pressure > this%pcmax) then
    capillary_pressure = this%pcmax
    dpc_dsatl = 0.d0
  endif
  
end subroutine SF_Linear_CapillaryPressure

! ************************************************************************** !

subroutine SF_Linear_Saturation(this,capillary_pressure, &
                                liquid_saturation,dsat_dpres,option)
  ! 
  ! Computes the saturation (and associated derivatives) as a function of 
  ! capillary pressure
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

! ************************************************************************** !
! ************************************************************************** !

function SF_mK_Create()

  ! Creates the modified Kosugi saturation function object

  implicit none

  class(sat_func_mK_type), pointer :: SF_mK_Create

  allocate(SF_mK_Create)
  call SF_mK_Create%Init()

end function SF_mK_Create

! ************************************************************************** !

subroutine SF_mK_Init(this)

  ! Initializes modified Kosugi saturation function object

  implicit none
  
  class(sat_func_mK_type) :: this

  call SFBaseInit(this)
  this%sigmaz = UNINITIALIZED_DOUBLE
  this%muz = UNINITIALIZED_DOUBLE
  this%rmax = UNINITIALIZED_DOUBLE
  this%r0 = UNINITIALIZED_DOUBLE
  this%nparam = UNINITIALIZED_INTEGER
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine SF_mK_Init

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

end subroutine SF_mK_Verify

! ************************************************************************** !

subroutine SF_mK_CapillaryPressure(this,liquid_saturation, &
                                   capillary_pressure,dpc_dsatl,option)
  !
  ! Computes the capillary_pressure as a function of saturation
  ! for modified Kosugi model.
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498�502. http://dx.doi.org/10.1111/gwat.12220
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
  PetscReal, intent(out) :: dpc_dsatl
  type(option_type), intent(inout) :: option

  PetscReal :: Se
  PetscReal :: inverse, exparg
  PetscReal :: hc, hmaxinv
  PetscReal :: dinverse_dSe
  PetscReal :: dSe_dsatl, dexparg_dinverse, dpc_dexparg
  PetscReal :: one_over_pc
  PetscReal :: tempreal

  dpc_dsatl = 0.d0

  if (liquid_saturation <= this%Sr) then
    capillary_pressure = this%pcmax
    return
  else if (liquid_saturation >= 1.d0) then
    capillary_pressure = 0.d0
    return
  endif

  dSe_dsatl = 1.d0/(1.d0 - this%Sr)
  Se = (liquid_saturation - this%Sr)*dSe_dsatl
!  inverse = -InverseNorm(Se)
  call InverseNorm(Se,inverse,PETSC_TRUE,dinverse_dSe)
  inverse = -1.d0*inverse
  dinverse_dSe = -1.d0*dinverse_dSe
  exparg = this%sigmaz*inverse + LNKAP - this%muz
  dexparg_dinverse = this%sigmaz

  hc = KAPPA/this%rmax
  dpc_dexparg = exp(exparg)
  capillary_pressure = dpc_dexparg + hc
  dpc_dsatl = dpc_dexparg*dexparg_dinverse*dinverse_dSe*dSe_dsatl
  if (this%nparam == 4) then
    hmaxinv = this%r0/KAPPA
    one_over_pc = 1.d0/capillary_pressure
    tempreal = 1.d0/(one_over_pc + hmaxinv)
    capillary_pressure = tempreal
    dpc_dsatl = capillary_pressure*tempreal*one_over_pc*one_over_pc*dpc_dsatl
  end if

  capillary_pressure = capillary_pressure*UNIT_CONVERSION
  dpc_dsatl = dpc_dsatl*UNIT_CONVERSION
  if (capillary_pressure > this%pcmax) then
    capillary_pressure = this%pcmax
    dpc_dsatl = 0.d0
  endif

end subroutine SF_mK_CapillaryPressure

! ************************************************************************** !

subroutine SF_mK_Saturation(this,capillary_pressure, &
                            liquid_saturation,dsat_dpres,option)
  !
  ! Computes the saturation (and associated derivatives) as a function of
  ! capillary pressure for modified Kosugi model
  !
  ! Malama, B. & K.L. Kuhlman, 2015. Unsaturated Hydraulic Conductivity
  ! Models Based on Truncated Lognormal Pore-size Distributions, Groundwater,
  ! 53(3):498�502. http://dx.doi.org/10.1111/gwat.12220
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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
    string = trim(name) // 'PERMEABILITY_FUNCTION,MUALEM_VG_LIQ'
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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

! ************************************************************************** !
! ************************************************************************** !

function RPF_Burdine_BC_Liq_Create()

  ! Creates the Brooks-Corey Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_BC_liq_type), pointer :: RPF_Burdine_BC_Liq_Create
  
  allocate(RPF_Burdine_BC_Liq_Create)
  call RPF_Burdine_BC_Liq_Create%Init()
  
end function RPF_Burdine_BC_Liq_Create

! ************************************************************************** !

subroutine RPF_Burdine_BC_Liq_Init(this)

  ! Initializes the Brooks-Corey Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_BC_liq_type) :: this

  call RPFBaseInit(this)
  this%lambda = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
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

  class(rpf_Burdine_BC_liq_type) :: this
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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
  if (Uninitialized(this%lambda)) then
    option%io_buffer = UninitializedMessage('LAMBDA',string)
    call printErrMsg(option)
  endif
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE

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

  class(rpf_Mualem_BC_liq_type) :: this
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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

! ************************************************************************** !
! ************************************************************************** !

function RPF_Mualem_Linear_Liq_Create()

  ! Creates the Linear Mualem relative permeability function object

  implicit none
  
  class(rpf_Mualem_linear_liq_type), pointer :: RPF_Mualem_Linear_Liq_Create
  
  allocate(RPF_Mualem_Linear_Liq_Create)
  call RPF_Mualem_Linear_Liq_Create%Init()
  
end function RPF_Mualem_Linear_Liq_Create

! ************************************************************************** !

subroutine RPF_Mualem_Linear_Liq_Init(this)

  ! Initializes the Linear Mualem relative permeability function object

  implicit none
  
  class(rpf_Mualem_Linear_liq_type) :: this

  call RPFBaseInit(this)
  this%alpha = UNINITIALIZED_DOUBLE
  this%pcmax = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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

! ************************************************************************** !
! ************************************************************************** !

function RPF_Burdine_Linear_Liq_Create()

  ! Creates the Linear Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_Linear_liq_type), pointer :: RPF_Burdine_Linear_Liq_Create
  
  allocate(RPF_Burdine_Linear_Liq_Create)
  call RPF_Burdine_Linear_Liq_Create%Init()
  
end function RPF_Burdine_Linear_Liq_Create

! ************************************************************************** !

subroutine RPF_Burdine_Linear_Liq_Init(this)

  ! Initializes the Linear Burdine relative permeability function object

  implicit none
  
  class(rpf_Burdine_Linear_liq_type) :: this

  call RPFBaseInit(this)
  
  this%analytical_derivative_available = PETSC_TRUE
  
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

! ************************************************************************** !
! ************************************************************************** !

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
  
  this%analytical_derivative_available = PETSC_TRUE
  
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

! ************************************************************************** !
! ************************************************************************** !

function RPF_mK_Liq_Create()

  ! Creates the modified Kosugi liq relative permeability function object

  implicit none

  class(rpf_mK_liq_type), pointer :: RPF_mK_Liq_Create

  allocate(RPF_mK_Liq_Create)
  call RPF_mK_Liq_Create%Init()

end function RPF_mK_Liq_Create

! ************************************************************************** !

subroutine RPF_mK_Liq_Init(this)

  ! Initializes modified Kosugi saturation function object

  implicit none
  
  class(rpf_mK_liq_type) :: this

  call RPFBaseInit(this)
  this%sigmaz = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_mK_Liq_Init

! ************************************************************************** !

subroutine RPF_mK_Liq_Verify(this,name,option)

  use Option_module

  implicit none

  class(rpf_mK_liq_type) :: this
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
  ! 53(3):498-502. http://dx.doi.org/10.1111/gwat.12220
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
  PetscReal :: dinvErfcRes_dSe

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

!  invErfcRes = InverseNorm(Se)
  call InverseNorm(Se,invErfcRes,PETSC_TRUE,dinvErfcRes_dSe)
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

! ************************************************************************** !
! ************************************************************************** !

function RPF_mK_Gas_Create()

  ! Creates the modified Kosugi gas relative permeability function object

  implicit none

  class(rpf_mK_gas_type), pointer :: RPF_mK_Gas_Create

  allocate(RPF_mK_Gas_Create)
  call RPF_mK_Gas_Create%Init()

end function RPF_mK_Gas_Create

! ************************************************************************** !

subroutine RPF_mK_Gas_Init(this)

  ! Initializes modified Kosugi saturation function object

  implicit none
  
  class(rpf_mK_gas_type) :: this

  call RPFBaseInit(this)
  this%sigmaz = UNINITIALIZED_DOUBLE
  
  this%analytical_derivative_available = PETSC_TRUE
  
end subroutine RPF_mK_Gas_Init

! ************************************************************************** !

subroutine RPF_mK_Gas_Verify(this,name,option)

  use Option_module

  implicit none

  class(rpf_mK_gas_type) :: this
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
  ! 53(3):498-502. http://dx.doi.org/10.1111/gwat.12220
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
  PetscReal :: dkr_Se
  PetscReal :: erfcArg, erfcRes
  PetscReal :: invErfcRes
  PetscReal :: sqrtSe, expArg
  PetscReal :: dinvErfcRes_dSeg

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

!  invErfcRes = InverseNorm(Seg)
  call InverseNorm(Seg,invErfcRes,PETSC_TRUE,dinvErfcRes_dSeg)
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

end subroutine RPF_mK_Gas_RelPerm

 
end module Characteristic_Curves_Common_module
