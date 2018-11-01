module EOS_Oil_module

#include "petsc/finclude/petscsys.h"
  use petscsys

  use PFLOTRAN_Constants_module
  use EOSData_module

  implicit none

  private


  ! module variables
  PetscReal :: fmw_oil           !kg/Kmol
  PetscReal :: constant_density  !kg/m3
  PetscReal :: constant_enthalpy
  PetscReal :: constant_viscosity
  PetscReal :: constant_sp_heat
  ! quadratic viscosity
  PetscReal :: quad_vis0
  PetscReal :: quad_vis_ref_pres(2)
  PetscReal :: quad_vis_ref_temp(2)
  PetscReal :: quad_vis_pres_coef(2)
  PetscReal :: quad_vis_temp_coef(2)
  !Reference Oil Density (e.g. used as surface density)
  PetscReal :: reference_density_kg
  ! parameters for linear density
  PetscReal :: compress_coeff      ! [kg/m3/Pa]
  PetscReal :: th_expansion_coeff  ! [kg/m3/°C]
  PetscReal :: den_linear_den0     ! [kg/m3]
  PetscReal :: den_linear_ref_pres ! [Pa]
  PetscReal :: den_linear_ref_temp ! [°C]

  ! quadratic enthalpy
  PetscReal :: quad_ent_ref_temp(2)
  PetscReal :: quad_ent_temp_coef(2)


  ! EOS databases
  class(eos_database_type), pointer :: eos_dbase
  class(eos_database_type), pointer :: eos_den_dbase
  class(eos_database_type), pointer :: eos_ent_dbase
  class(eos_database_type), pointer :: eos_vis_dbase

  ! PVT tables - eos_tables
  class(eos_table_type), pointer :: pvt_table

  ! when adding a new eos_database, remember to add it to EOSOilDBaseDestroy()

  ! In order to support generic EOS subroutines, we need the following:
  ! 1. An interface declaration that defines the argument list (best to have
  !    "Dummy" appended.
  ! 2. A procedure pointer that is initially set to null.  This pointer is
  !    pointed to the appropriate subroutine later on (e.g. EOSOilInit())
  ! 3. An interface for derivative/non-derivative versions

  ! procedure pointers

  procedure(EOSOilViscosityDummy), pointer :: EOSOilViscosityPtr => null()
  procedure(EOSOilDensityDummy), pointer :: EOSOilDensityPtr => null()
  procedure(EOSOilEnthalpyDummy), pointer :: EOSOilEnthalpyPtr => null()
  procedure(EOSOilDensityEnergyDummy), pointer :: &
    EOSOilDensityEnergyPtr => null()
  Procedure(EOSOilRSDummy), pointer :: EOSOilRSPtr => null()
  procedure(EOSOilCompressibilityDummy), pointer :: &
    EOSOilCompressibilityPtr => null()
  procedure(EOSOilViscosibilityDummy), pointer :: &
    EOSOilViscosibilityPtr => null()

  abstract interface
    subroutine EOSOilViscosityDummy(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr, &
                                    table_idxs)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! oil pressure [Pa]
      PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Vis     ! oil viscosity
      PetscReal, intent(out) :: dVis_dT ! derivative oil viscosity wrt temperature
      PetscReal, intent(out) :: dVis_dP ! derivative oil viscosity wrt Pressure
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSOilViscosityDummy
    subroutine EOSOilDensityDummy(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr, &
                                  table_idxs)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSOilDensityDummy
    subroutine EOSOilEnthalpyDummy(T,P,deriv,H,dH_dT,dH_dP,ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
      PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
      PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSOilEnthalpyDummy
    subroutine EOSOilDensityEnergyDummy(T,P,deriv,Rho,dRho_dT,dRho_dP, &
                                        H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr, &
                                        table_idxs)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
      PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
      PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
      PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
      PetscReal, intent(out) :: U       ! internal energy [J/kmol]
      PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
      PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSOilDensityEnergyDummy
    subroutine EOSOilRSDummy(T,P,deriv,RS,dRS_dT,dRS_dP,ierr,table_idxs)
      implicit none
      PetscReal, intent(in) :: T       ! temperature [C]
      PetscReal, intent(in) :: P       ! oil pressure [Pa]
      PetscBool, intent(in) :: deriv   ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: RS     ! Vg[sm3]/Vo[sm3] per unit of res. volume
      PetscReal, intent(out) :: dRS_dT ! derivative RS wrt temperature
      PetscReal, intent(out) :: dRS_dP ! derivative RS wrt Pressure
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSOilRSDummy
    subroutine EOSOilCompressibilityDummy(T,P,deriv,Co,dCo_dT,dCo_dP,ierr, &
                                          table_idxs)
      implicit none
      PetscReal, intent(in) :: T       ! temperature [C]
      PetscReal, intent(in) :: P       ! oil pressure [Pa]
      PetscBool, intent(in) :: deriv   ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Co     ! oil compressibility [1/Pa]
      PetscReal, intent(out) :: dCo_dT ! derivative Co wrt temperature
      PetscReal, intent(out) :: dCo_dP ! derivative Co wrt Pressure
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSOilCompressibilityDummy
    subroutine EOSOilViscosibilityDummy(T,P,deriv,Cvis,dCvis_dT,dCvis_dP,ierr, &
                                          table_idxs)
      implicit none
      PetscReal, intent(in) :: T       ! temperature [C]
      PetscReal, intent(in) :: P       ! oil pressure [Pa]
      PetscBool, intent(in) :: deriv   ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Cvis     ! oil viscosibility [1/Pa]
      PetscReal, intent(out) :: dCvis_dT ! derivative Cvis wrt temperature
      PetscReal, intent(out) :: dCvis_dP ! derivative Cvis wrt Pressure
      PetscErrorCode, intent(out) :: ierr
      PetscInt, pointer, optional, intent(inout) :: table_idxs(:)
    end subroutine EOSOilViscosibilityDummy
  end interface

  ! interfaces for derivative/non-derivative versions that are visible outside
  ! the module.
  interface EOSOilViscosity
    procedure EOSOilViscosityNoDerive
    procedure EOSOilViscosityDerive
  end interface
  interface EOSOilDensity
    procedure EOSOilDensityNoDerive
    procedure EOSOilDensityDerive
  end interface
  interface EOSOilEnthalpy
    procedure EOSOilEnthalpyNoDerive
    procedure EOSOilEnthalpyDerive
  end interface
  interface EOSOilDensityEnergy
    procedure EOSOilDenEnergyNoDerive
    procedure EOSOilDenEnergyDerive
  end interface
  interface EOSOilRS
    procedure EOSOilRSNoDerive
    procedure EOSOilRSDerive
  end interface
  interface EOSOilCompressibility
    procedure EOSOilCompressibilityNoDerive
    procedure EOSOilCompressibilityDerive
  end interface
  interface EOSOilViscosibility
    procedure EOSOilViscosibilityNoDerive
    procedure EOSOilViscosibilityDerive
  end interface

  public :: EOSOilInit, &
            EOSOilVerify, &
            EOSOilViscosity, &
            EOSOilDensity, &
            EOSOilEnthalpy, &
            EOSOilDensityEnergy, &
            EOSOilRS, &
            EOSOilCompressibility, &
            EOSOilViscosibility, &
            EOSOilInputRecord

  public :: EOSOilSetFMWConstant, &
            EOSOilGetFMW, &
            EOSOilSetReferenceDensity, &
            EOSOilGetReferenceDensity, &
            EOSOilSetViscosityConstant, &
            EOSOilSetViscosityQuad, &
            EOSOilSetVisQuadRefVis, &
            EOSOilSetVisQuadRefPres, &
            EOSOilSetVisQuadRefTemp, &
            EOSOilSetVisQuadPresCoef, &
            EOSOilSetVisQuadTempCoef, &
            EOSOilSetVisDBase, &
            EOSOilSetVisLinLogInterp, &
            EOSOilSetDensityConstant, &
            EOSOilSetDensityLinear, &
            EOSOilSetDensityInverseLinear, &
            EOSOilSetDenLinearRefDen, &
            EOSOilSetDenLinearComprCoef, &
            EOSOilSetDenLinearExpanCoef, &
            EOSOilSetDenLinearRefPres, &
            EOSOilSetDenLinearRefTemp, &
            EOSOilSetDenDBase, &
            EOSOilSetEnthalpyConstant, &
            EOSOilSetEnthalpyLinearTemp, &
            EOSOilSetEnthalpyQuadraticTemp, &
            EOSOilSetEntQuadRefTemp, &
            EOSOilSetEntQuadTempCoef, &
            EOSOilSetEntDBase, &
            EOSOilSetEOSDBase, &
            EOSOilSetPVDO, &
            EOSOilSetPVCO, &
            EOSOilTableProcess, &
            EOSOilDBaseDestroy

contains

! ************************************************************************** !

subroutine EOSOilInit()

  implicit none

  constant_density = UNINITIALIZED_DOUBLE
  constant_viscosity = UNINITIALIZED_DOUBLE
  constant_enthalpy = UNINITIALIZED_DOUBLE

  quad_vis0 = UNINITIALIZED_DOUBLE
  quad_vis_ref_pres(1:2) = UNINITIALIZED_DOUBLE
  quad_vis_ref_temp(1:2) = UNINITIALIZED_DOUBLE
  quad_vis_pres_coef(1:2) = UNINITIALIZED_DOUBLE
  quad_vis_temp_coef(1:2) = UNINITIALIZED_DOUBLE

  reference_density_kg = UNINITIALIZED_DOUBLE

  compress_coeff = UNINITIALIZED_DOUBLE
  th_expansion_coeff = UNINITIALIZED_DOUBLE
  den_linear_den0 = UNINITIALIZED_DOUBLE
  den_linear_ref_pres = UNINITIALIZED_DOUBLE
  den_linear_ref_temp = UNINITIALIZED_DOUBLE

  quad_ent_ref_temp(1:2) = UNINITIALIZED_DOUBLE
  quad_ent_temp_coef(1:2) = UNINITIALIZED_DOUBLE


  fmw_oil = FMWOIL !default oil formula weight C10H22 (142 g/mol)

  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS

  nullify(EOSOilViscosityPtr)
  nullify(EOSOilDensityPtr)
  nullify(EOSOilEnthalpyPtr)

  nullify(eos_dbase)
  nullify(eos_den_dbase)
  nullify(eos_ent_dbase)
  nullify(eos_vis_dbase)

  ! could decide for a default model, but it only if there is one that
  ! does nto require input parameters
  !EOSOilViscosityPtr => EOSOilViscosityConstant
  !EOSOilDensityPtr => EOSOilDensityConstant
  !EOSOilEnthalpyPtr => EOSOilEnthalpyConstant

end subroutine EOSOilInit

! ************************************************************************** !

subroutine EOSOilVerify(ierr,error_string)
  !
  ! Author: Paolo Orsini
  !
  ! to do : add check on unitialized coeffiecinets for linear and quadratic
  !         functions (density, viscosity, enthalpy)
  !
  implicit none

  PetscErrorCode, intent(out) :: ierr
  character(len=MAXSTRINGLENGTH), intent(out) :: error_string

  ierr = 0

  error_string = ''
  if (.not.associated(EOSOilDensityPtr) ) then
    error_string = trim(error_string) // ' Oil Density model not defined'
    ierr = 1
    return
  end if
  if (.not.associated(EOSOilEnthalpyPtr) ) then
    error_string = trim(error_string) // ' Oil Enthalpy model not defined'
    ierr = 1
    return
  end if
  if (.not.associated(EOSOilViscosityPtr) ) then
    error_string = trim(error_string) // ' Oil Viscosity model not defined'
    ierr = 1
    return
  end if

  if ( associated(EOSOilDensityPtr,EOSOilDensityConstant).and. &
      Uninitialized(constant_density) &
     ) then
    error_string = trim(error_string) // &
    ' Oil Constant Density selcected without providing a value'
    ierr = 1
    return
  end if

  if ( associated(EOSOilDensityPtr,EOSOilDensityEOSDBase) ) then
    if ( .not.eos_dbase%EOSPropPresent(EOS_DENSITY) ) then
      error_string = trim(error_string) // &
      ' Oil Density to be interpolated from database = ' // &
      eos_dbase%file_name // &
      ' which does not have data for density'
      ierr = 1
      return
    end if
  end if

  if ( associated(EOSOilDensityPtr,EOSOilDensityDenDBase) ) then
    if ( .not.eos_den_dbase%EOSPropPresent(EOS_DENSITY) ) then
      error_string = trim(error_string) // &
      ' Oil Density to be interpoalted from database = ' // &
      eos_dbase%file_name // &
      ' which does not have data for density'
      ierr = 1
      return
    end if
  end if

  if ( associated(EOSOilEnthalpyPtr,EOSOilEnthalpyConstant).and. &
      Uninitialized(constant_enthalpy) &
     ) then
    error_string = trim(error_string) // &
    ' Oil Constant Enthalpy selcected without providing a value'
    ierr = 1
    return
  end if

  if ( associated(EOSOilEnthalpyPtr,EOSOilEnthalpyEOSDBase) ) then
    if ( .not.eos_dbase%EOSPropPresent(EOS_ENTHALPY) ) then
      error_string = trim(error_string) // &
      ' Oil Enthalpy to be interpolated from database = ' // &
      eos_dbase%file_name // &
      ' which does not have data for Enthalpy'
      ierr = 1
      return
    end if
  end if

  if ( associated(EOSOilEnthalpyPtr,EOSOilEnthalpyEntDBase) ) then
    if ( .not.eos_ent_dbase%EOSPropPresent(EOS_ENTHALPY) ) then
      error_string = trim(error_string) // &
      ' Oil Density to be interpolated from database = ' // &
      eos_dbase%file_name // &
      ' which does not have data for enthalpy'
      ierr = 1
      return
    end if
  end if

  if ( associated(EOSOilViscosityPtr,EOSOilViscosityConstant).and. &
      Uninitialized(constant_viscosity) &
     ) then
    error_string = trim(error_string) // &
    ' Oil Constant Viscosity selcected without providing a value'
    ierr = 1
    return
  end if

  if ( associated(EOSOilViscosityPtr,EOSOilViscosityEOSDBase) ) then
    if ( .not.eos_dbase%EOSPropPresent(EOS_VISCOSITY) ) then
      error_string = trim(error_string) // &
      ' Oil Enthalpy to be interpolated from database = ' // &
      eos_dbase%file_name // &
      ' which does not have data for Viscosity'
      ierr = 1
      return
    end if
  end if

  if ( associated(EOSOilViscosityPtr,EOSOilViscosityVisDBase) ) then
    if ( .not.eos_vis_dbase%EOSPropPresent(EOS_VISCOSITY) ) then
      error_string = trim(error_string) // &
      ' Oil Density to be interpolated from database = ' // &
      eos_dbase%file_name // &
      ' which does not have data for visosity'
      ierr = 1
      return
    end if
  end if

  if ( associated(pvt_table) ) then
    if(Uninitialized(reference_density_kg)) then
      error_string = trim(error_string) // &
      'A reference (e.g. Surface) density must be specified ' // &
      'using either REFERENCE_DENSITY, SURFACE_DENSITY or STANDARD_DENSITY '
      ierr = 1
    end if
  end if

end subroutine EOSOilVerify

! ************************************************************************** !

subroutine EOSOilSetFMWConstant(fmw_input)

  implicit none

  PetscReal :: fmw_input

  fmw_oil = fmw_input

end subroutine EOSOilSetFMWConstant

! ************************************************************************** !

function EOSOilGetFMW()

  implicit none

  PetscReal :: EOSOilGetFMW

  EOSOilGetFMW = fmw_oil

end function EOSOilGetFMW

! ************************************************************************** !

subroutine EOSOilSetReferenceDensity(input_ref_density)

  implicit none

  PetscReal :: input_ref_density

  reference_density_kg = input_ref_density

end subroutine EOSOilSetReferenceDensity

! ************************************************************************** !

function EOSOilGetReferenceDensity()

  implicit none

  PetscReal :: EOSOilGetReferenceDensity

  EOSOilGetReferenceDensity= reference_density_kg

end function EOSOilGetReferenceDensity

! ************************************************************************** !

subroutine EOSOilSetViscosityConstant(viscosity)

  implicit none

  PetscReal :: viscosity

  constant_viscosity = viscosity
  EOSOilViscosityPtr => EOSOilViscosityConstant

end subroutine EOSOilSetViscosityConstant

! ************************************************************************** !

subroutine EOSOilSetViscosityQuad()

  implicit none

  EOSOilViscosityPtr => EOSOilQuadViscosity

end subroutine EOSOilSetViscosityQuad

! ************************************************************************** !

subroutine EOSOilSetVisQuadRefVis(vis0)

  implicit none

  PetscReal :: vis0

  quad_vis0 = vis0

end subroutine EOSOilSetVisQuadRefVis

! ************************************************************************** !

subroutine EOSOilSetVisQuadRefPres(p1,p2)

  implicit none

  PetscReal :: p1, p2

  quad_vis_ref_pres(1) = p1
  quad_vis_ref_pres(2) = p2

end subroutine EOSOilSetVisQuadRefPres

! ************************************************************************** !

subroutine EOSOilSetVisQuadRefTemp(t1,t2)

  implicit none

  PetscReal :: t1, t2

  quad_vis_ref_temp(1) = t1
  quad_vis_ref_temp(2) = t2

end subroutine EOSOilSetVisQuadRefTemp

! ************************************************************************** !

subroutine EOSOilSetVisQuadPresCoef(a1,a2)

  implicit none

  PetscReal :: a1, a2

  quad_vis_pres_coef(1) = a1
  quad_vis_pres_coef(2) = a2

end subroutine EOSOilSetVisQuadPresCoef

! ************************************************************************** !

subroutine EOSOilSetVisQuadTempCoef(b1,b2)

  implicit none

  PetscReal :: b1, b2

  quad_vis_temp_coef(1) = b1
  quad_vis_temp_coef(2) = b2

end subroutine EOSOilSetVisQuadTempCoef

! ************************************************************************** !

subroutine EOSOilSetVisDBase(filename,option)

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: filename
  type(option_type) :: option

  eos_vis_dbase => EOSDatabaseCreate(filename,'oil_den_database')
  call eos_vis_dbase%Read(option)

  !set property function pointers
  EOSOilViscosityPtr => EOSOilViscosityVisDBase

end subroutine EOSOilSetVisDBase

! ************************************************************************** !

subroutine EOSOilSetVisLinLogInterp(option)

  use Option_module

  implicit none

  type(option_type) :: option

  if ( associated(eos_dbase) ) then
    call eos_dbase%SetupVarLinLogInterp(EOS_VISCOSITY,option)
  end if

  if ( associated(eos_vis_dbase) ) then
    call eos_vis_dbase%SetupVarLinLogInterp(EOS_VISCOSITY,option)
  end if
  
  if ( associated(pvt_table) ) then
    call pvt_table%SetupVarLinLogInterp(EOS_VISCOSITY,option)
  end if  

end subroutine EOSOilSetVisLinLogInterp

! ************************************************************************** !

! ************************************************************************** !

subroutine EOSOilSetDensityConstant(density)

  implicit none

  PetscReal :: density

  constant_density = density
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS
  EOSOilDensityPtr => EOSOilDensityConstant

end subroutine EOSOilSetDensityConstant

! ************************************************************************** !

subroutine EOSOilSetDensityLinear()

  implicit none

  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS
  EOSOilDensityPtr => EOSOilDensityLinear

end subroutine EOSOilSetDensityLinear

! ************************************************************************** !

subroutine EOSOilSetDensityInverseLinear()

  implicit none

  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS
  EOSOilDensityPtr => EOSOilDensityInverseLinear

end subroutine EOSOilSetDensityInverseLinear


! ************************************************************************** !

subroutine EOSOilSetDenLinearRefDen(den0)

  implicit none

  PetscReal :: den0

  den_linear_den0 = den0

end subroutine EOSOilSetDenLinearRefDen

subroutine EOSOilSetDenLinearRefPres(ref_pres)

  implicit none

  PetscReal :: ref_pres

  den_linear_ref_pres = ref_pres

end subroutine EOSOilSetDenLinearRefPres

! ************************************************************************** !

subroutine EOSOilSetDenLinearRefTemp(ref_temp)

  implicit none

  PetscReal :: ref_temp

  den_linear_ref_temp = ref_temp

end subroutine EOSOilSetDenLinearRefTemp

! ************************************************************************** !

subroutine EOSOilSetDenLinearComprCoef(compress_c)

  implicit none

  PetscReal :: compress_c

  compress_coeff = compress_c

end subroutine EOSOilSetDenLinearComprCoef

! ************************************************************************** !

subroutine EOSOilSetDenLinearExpanCoef(expansion_c)

  implicit none

  PetscReal :: expansion_c

  th_expansion_coeff = expansion_c

end subroutine EOSOilSetDenLinearExpanCoef

! ************************************************************************** !

subroutine EOSOilSetDenDBase(filename,option)

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: filename
  type(option_type) :: option

  eos_den_dbase => EOSDatabaseCreate(filename,'oil_den_database')
  call eos_den_dbase%Read(option)

  !set property function pointers
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS
  EOSOilDensityPtr => EOSOilDensityDenDBase

end subroutine EOSOilSetDenDBase

! ************************************************************************** !

subroutine EOSOilSetEnthalpyConstant(enthalpy)

  implicit none

  PetscReal :: enthalpy

  constant_enthalpy = enthalpy
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS
  EOSOilEnthalpyPtr => EOSOilEnthalpyConstant

end subroutine EOSOilSetEnthalpyConstant

! ************************************************************************** !

subroutine EOSOilSetEnthalpyLinearTemp(specific_heat)

  implicit none

  PetscReal :: specific_heat

  constant_sp_heat = specific_heat
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS
  EOSOilEnthalpyPtr => EOSOilEnthalpyLinearTemp

  !write(*,*) "I am in EOS oil linear set up"

end subroutine EOSOilSetEnthalpyLinearTemp

! ************************************************************************** !

subroutine EOSOilSetEnthalpyQuadraticTemp()

  implicit none

  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS
  EOSOilEnthalpyPtr => EOSOilEnthalpyQuadTemp


end subroutine EOSOilSetEnthalpyQuadraticTemp

! ************************************************************************** !

subroutine EOSOilSetEntQuadRefTemp(t1,t2)

  implicit none

  PetscReal :: t1, t2

  quad_ent_ref_temp(1) = t1
  quad_ent_ref_temp(2) = t2

end subroutine EOSOilSetEntQuadRefTemp

! ************************************************************************** !

subroutine EOSOilSetEntQuadTempCoef(c1,c2)

  implicit none

  PetscReal :: c1, c2

  quad_ent_temp_coef(1) = c1
  quad_ent_temp_coef(2) = c2

end subroutine EOSOilSetEntQuadTempCoef

! ************************************************************************** !

subroutine EOSOilSetEntDBase(filename,option)

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: filename
  type(option_type) :: option

  eos_ent_dbase => EOSDatabaseCreate(filename,'oil_ent_database')
  call eos_ent_dbase%Read(option)

  !set property function pointers
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS
  EOSOilEnthalpyPtr => EOSOilEnthalpyEntDBase


end subroutine EOSOilSetEntDBase

! ************************************************************************** !

subroutine EOSOilSetEOSDBase(filename,option)

  use Option_module

  implicit none

  character(len=MAXWORDLENGTH) :: filename
  type(option_type) :: option

  eos_dbase => EOSDatabaseCreate(filename,'oil_database')
  call eos_dbase%Read(option)

  !set property function pointers
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyS
  if (.not.associated(EOSOilDensityPtr))  &
            EOSOilDensityPtr => EOSOilDensityEOSDBase
  if (.not.associated(EOSOilEnthalpyPtr)) &
    EOSOilEnthalpyPtr => EOSOilEnthalpyEOSDBase
  if (.not.associated(EOSOilViscosityPtr)) &
    EOSOilViscosityPtr => EOSOilViscosityEOSDBase

end subroutine EOSOilSetEOSDBase

! ************************************************************************** !

subroutine EOSOilViscosityConstant(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr, &
                                   table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Vis     ! oil viscosity
  PetscReal, intent(out) :: dVis_dT ! derivative oil viscosity wrt temperature
  PetscReal, intent(out) :: dVis_dP ! derivative oil viscosity wrt Pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:) !pvt table indices

  Vis = constant_viscosity

  dVis_dT = 0.0d0
  dVis_dP = 0.0d0

end subroutine EOSOilViscosityConstant

! ************************************************************************** !

subroutine EOSOilQuadViscosity(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr, &
                               table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Vis     ! oil viscosity
  PetscReal, intent(out) :: dVis_dT ! derivative oil viscosity wrt temperature
  PetscReal, intent(out) :: dVis_dP ! derivative oil viscosity wrt Pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:) !pvt table indices


  Vis = quad_vis0 + &
        quad_vis_pres_coef(1) * ( P - quad_vis_ref_pres(1) ) + &
        quad_vis_pres_coef(2) * ( P - quad_vis_ref_pres(2) )**2.0d0 + &
        quad_vis_temp_coef(1) * ( T - quad_vis_ref_temp(1) ) + &
        quad_vis_temp_coef(2) * ( T - quad_vis_ref_temp(2) )**2.0d0

  if (deriv) then
    dVis_dP = quad_vis_pres_coef(1) + &
              2.0d0 * quad_vis_pres_coef(2) * ( P - quad_vis_ref_pres(2) )
    dVis_dT = quad_vis_temp_coef(1) + &
              2.0d0 * quad_vis_temp_coef(2) * ( T - quad_vis_ref_temp(2) )
  end if

end subroutine EOSOilQuadViscosity

! ************************************************************************** !

subroutine EOSOilViscosityEOSDBase(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr, &
                                   table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Vis     ! oil viscosity [Pa-s]
  PetscReal, intent(out) :: dVis_dT ! derivative oil viscosity wrt db temperature [Pa-s/C]
  PetscReal, intent(out) :: dVis_dP ! derivative oil viscosity wrt db Pressure [Pa-s/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:) !pvt table indices

  !ierr initialised in EOSEOSProp
  !call eos_dbase%EOSProp(T,P,EOS_VISCOSITY,Vis,ierr)
  
  call eos_dbase%EOSPropGrad(T,P,EOS_VISCOSITY,Vis,dVis_dT,dVis_dP,ierr) 

  ! initialize to derivative to NaN so that not mistakenly used.
  !dVis_dT = InitToNan()
  !dVis_dP = InitToNan()

  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilViscosityEOSDBase - Viscosity derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilViscosityEOSDBase

! ************************************************************************** !

subroutine EOSOilViscosityVisDBase(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr, &
                                   table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Vis     ! oil viscosity [Pa-s]
  PetscReal, intent(out) :: dVis_dT ! derivative oil viscosity wrt temperature [Pa-s/C]
  PetscReal, intent(out) :: dVis_dP ! derivative oil viscosity wrt Pressure [Pa-s/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  !ierr initialised in EOSEOSProp
  !call eos_vis_dbase%EOSProp(T,P,EOS_VISCOSITY,Vis,ierr)
  call eos_vis_dbase%EOSPropGrad(T,P,EOS_VISCOSITY,Vis,dVis_dT,dVis_dP,ierr)

  ! initialize to derivative to NaN so that not mistakenly used.
  ! dVis_dT = InitToNan()
  ! dVis_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilViscosityVisDBase - Viscosity derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilViscosityVisDBase

! ************************************************************************** !

subroutine EOSOilViscosityTable(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr, &
                                   table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Vis     ! oil viscosity [Pa-s]
  PetscReal, intent(out) :: dVis_dT ! derivative oil viscosity wrt table temperature [Pa-s/C]
  PetscReal, intent(out) :: dVis_dP ! derivative oil viscosity wrt table Pressure [Pa-s/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:) !pvt table indices

  !ierr initialised in EOSEOSProp
  !call pvt_table%EOSProp(T,P,EOS_VISCOSITY,Vis,table_idxs,ierr)
  
  call pvt_table%EOSPropGrad(T,P,EOS_VISCOSITY,Vis,dVis_dT,dVis_dP, &
                             table_idxs,ierr)

  ! initialize to derivative to NaN so that not mistakenly used.
  ! dVis_dT = InitToNan()
  ! dVis_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilViscosityTable - Viscosity derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilViscosityTable

! ************************************************************************** !

subroutine EOSOilViscosityNoDerive(T,P,Rho,Vis,ierr,table_idxs)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]
  PetscReal, intent(out) :: Vis     ! oil viscosity
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:) !pvt table indices

  PetscReal :: dum1, dum2

  call EOSOilViscosityPtr(T,P,Rho,PETSC_FALSE,Vis,dum1,dum2,ierr,table_idxs)

end subroutine EOSOilViscosityNoDerive

! ************************************************************************** !

subroutine EOSOilViscosityDerive(T,P,Rho,Vis,dVis_dT,dVis_dP,ierr,table_idxs)
  
  implicit none

  PetscReal, intent(in) :: T         ! temperature [C]
  PetscReal, intent(in) :: P         ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho       ! oil density [kmol/m3]
  PetscReal, intent(out) :: Vis      ! oil viscosity [Pa-s]
  PetscReal, intent(out) :: dVis_dT  ! oil visc derviv WRT Temp [Pa-s/C]
  PetscReal, intent(out) :: dVis_dP  ! oil visc derviv WRT Press [Pa-s/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:) !pvt table indices


  call EOSOilViscosityPtr(T,P,Rho,PETSC_TRUE,Vis,dVis_dT,dVis_dP, &
                          ierr,table_idxs)

end subroutine EOSOilViscosityDerive

! ************************************************************************** !

subroutine EOSOilDensityConstant(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr, &
                                 table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

        ! kg/m3 * kmol/kg  = kmol/m3
  Rho = constant_density / fmw_oil ! kmol/m^3

  dRho_dT = 0.d0
  dRho_dP = 0.d0

end subroutine EOSOilDensityConstant

! ************************************************************************** !

subroutine EOSOilDensityLinear(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr, &
                               table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  Rho = den_linear_den0 + &
        compress_coeff * (P - den_linear_ref_pres ) - & ! compression
        th_expansion_coeff * (T - den_linear_ref_temp )    ! expansion

  ! conversion to molar density
        ! kg/m3 * kmol/kg  = kmol/m3
  Rho = Rho / fmw_oil ! kmol/m^3

  if (deriv) then
#if 0
    dRho_dT = compress_coeff / fmw_oil
    dRho_dP = - th_expansion_coeff / fmw_oil
#endif
    dRho_dP = compress_coeff / fmw_oil
    dRho_dT = - th_expansion_coeff / fmw_oil
  end if

end subroutine EOSOilDensityLinear

! ************************************************************************** !

subroutine EOSOilDensityInverseLinear(T,P,deriv,Rho,dRho_dT,dRho_dP,ierr, &
                                      table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  Rho = den_linear_den0 / &
        ( 1.0d0 - compress_coeff * (P - den_linear_ref_pres ) ) / &  ! compression
        ( 1.0d0 + th_expansion_coeff * (T - den_linear_ref_temp ) )  ! expansion

  ! conversion to molar density
        ! kg/m3 * kmol/kg  = kmol/m3
  Rho = Rho / fmw_oil ! kmol/m^3

  if (deriv) then
    dRho_dT = - th_expansion_coeff * den_linear_den0 *  &
        ( 1.0d0 + th_expansion_coeff * (T - den_linear_ref_temp) )**(-2.0) / &
        ( 1.0d0 - compress_coeff * (P - den_linear_ref_pres ) ) / &
         fmw_oil
    dRho_dP = compress_coeff * den_linear_den0 * &
        ( 1.0d0 - compress_coeff * (P - den_linear_ref_pres ) )**(-2.0) / &
        ( 1.0d0 + th_expansion_coeff * (T - den_linear_ref_temp ) ) / &
         fmw_oil
  end if

end subroutine EOSOilDensityInverseLinear

! ************************************************************************** !

subroutine EOSOilDensityEOSDBase(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr, &
                                 table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt db temperature [kmol/C]
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt db pressure [kmol/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  !ierr initialised in EOSEOSProp
  !call eos_dbase%EOSProp(T,P,EOS_DENSITY,Rho,ierr)

  call eos_dbase%EOSPropGrad(T,P,EOS_DENSITY,Rho,dRho_dT,dRho_dP,ierr)
  
  !PO todo: conversion when loaidng database to do this operation only once
  ! conversion to molar density
        ! kg/m3 * kmol/kg  = kmol/m3
  Rho = Rho / fmw_oil ! kmol/m^3
  dRho_dT = dRho_dT / fmw_oil
  dRho_dP = dRho_dP / fmw_oil

  ! initialize to derivative to NaN so that not mistakenly used.
  ! dRho_dT = InitToNan()
  ! dRho_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilDensityEOSDBase - Den derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilDensityEOSDBase

! ************************************************************************** !

subroutine EOSOilDensityDenDBase(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr, &
                                 table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt db temperature [kmol/C]
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt db pressure [kmol/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  !ierr initialised in EOSEOSProp
  !call eos_den_dbase%EOSProp(T,P,EOS_DENSITY,Rho,ierr)

  call eos_den_dbase%EOSPropGrad(T,P,EOS_DENSITY,Rho,dRho_dT,dRho_dP,ierr)

  !PO todo: conversion when loaidng database to do this operation only once
  ! conversion to molar density
        ! kg/m3 * kmol/kg  = kmol/m3
  Rho = Rho / fmw_oil ! kmol/m^3
  dRho_dT = dRho_dT / fmw_oil
  dRho_dP = dRho_dP / fmw_oil

  ! initialize derivative to NaN so that not mistakenly used.
  ! dRho_dT = InitToNan()
  ! dRho_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilDensityDenDBase - Den derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilDensityDenDBase

! ************************************************************************** !

subroutine EOSOilDensityTable(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr, &
                                 table_idxs)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt table temperature [kmol/C]
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt table pressure [kmol/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  !ierr initialised in EOSProp
  !Rho from pvt table is already in kmol/m3
  !call pvt_table%EOSProp(T,P,EOS_DENSITY,Rho,table_idxs,ierr)
  
  !Rho from pvt table is already in kmol/m3
  call pvt_table%EOSPropGrad(T,P,EOS_DENSITY,Rho,dRho_dT,dRho_dP, &
                             table_idxs,ierr)

  ! initialize derivative to NaN so that not mistakenly used.
  ! dRho_dT = InitToNan()
  ! dRho_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilDensityTable - Den derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilDensityTable

! ************************************************************************** !

subroutine EOSOilDensityNoDerive(T,P,Rho,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dum1, dum2

  ! derivatives are so cheap, just compute them
  call EOSOilDensityPtr(T, P, PETSC_FALSE, Rho, dum1, dum2,ierr,table_idxs)

end subroutine EOSOilDensityNoDerive

! ************************************************************************** !

subroutine EOSOilDensityDerive(T,P,Rho,dRho_dT,dRho_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call EOSOilDensityPtr(T,P,PETSC_TRUE,Rho,dRho_dT,dRho_dP,ierr,table_idxs)

end subroutine EOSOilDensityDerive

! ************************************************************************** !

subroutine EOSOilEnthalpyConstant(T,P,deriv,H,dH_dT,dH_dP,ierr)
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  H = constant_enthalpy ! J/kmol

  dH_dP = 0.d0
  dH_dT = 0.d0

end subroutine EOSOilEnthalpyConstant

! ************************************************************************** !

subroutine EOSOilEnthalpyLinearTemp(T,P,deriv,H,dH_dT,dH_dP,ierr)
  !
  ! Author: Paolo Orsini
  ! Date: 12/11/15
  !
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  H = constant_sp_heat * T * fmw_oil ! J/(kg °C) °C * Kg/Kmol = J/Kmol

  dH_dT = UNINITIALIZED_DOUBLE
  dH_dP = UNINITIALIZED_DOUBLE

  if (deriv) then
    dH_dP = 0.d0
    dH_dT = constant_sp_heat * fmw_oil
  end if

end subroutine EOSOilEnthalpyLinearTemp

! ************************************************************************** !

subroutine EOSOilEnthalpyQuadTemp(T,P,deriv,H,dH_dT,dH_dP,ierr)

  ! Author: Paolo Orsini
  ! Date: 6/23/16
  !
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  !H = constant_sp_heat * T * fmw_oil ! J/(kg °C) °C * Kg/Kmol = J/Kmol


       ! [ J/(kg °C) * °C ] * [ J/(kg °C °C) °C °C ] * Kg/Kmol = J/Kmol
  H = ( quad_ent_temp_coef(1) * (T-quad_ent_ref_temp(1)) + &
        0.5d0 * quad_ent_temp_coef(2) * (T-quad_ent_ref_temp(2))**2.0d0 ) &
        * fmw_oil

  dH_dT = UNINITIALIZED_DOUBLE
  dH_dP = UNINITIALIZED_DOUBLE

  if (deriv) then
    dH_dP = 0.d0
    dH_dT = ( quad_ent_temp_coef(1) + &
              quad_ent_temp_coef(2) * (T-quad_ent_ref_temp(2)) ) * fmw_oil
  end if

end subroutine EOSOilEnthalpyQuadTemp

! ************************************************************************** !


subroutine EOSOilEnthalpyEOSDBase(T,P,deriv,H,dH_dT,dH_dP,ierr)
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt db temperature [J/kmol/C]
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt db pressure [J/kmol/Pa]
  PetscErrorCode, intent(out) :: ierr

  !ierr initialised in EOSEOSProp
  !call eos_dbase%EOSProp(T,P,EOS_ENTHALPY,H,ierr)

  call eos_dbase%EOSPropGrad(T,P,EOS_ENTHALPY,H,dH_dT,dH_dP,ierr)

  !PO todo: conversion when loaidng database to do this operation only once
  ! conversion to molar energy
  ! J/kg * kg/Kmol = J/Kmol
  H = H  * fmw_oil
  dH_dT = dH_dT * fmw_oil
  dH_dP = dH_dP * fmw_oil
  
  ! initialize derivative to NaN so that not mistakenly used.
  ! dH_dT = InitToNan()
  ! dH_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilEnthalpyEOSDBase - H derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilEnthalpyEOSDBase

! ************************************************************************** !

subroutine EOSOilEnthalpyEntDBase(T,P,deriv,H,dH_dT,dH_dP,ierr)
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt db temperature [J/kmol/T]
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt db pressure [J/kmol/Pa]
  PetscErrorCode, intent(out) :: ierr

  !ierr initialised in EOSEOSProp
  !call eos_ent_dbase%EOSProp(T,P,EOS_ENTHALPY,H,ierr)

  call eos_ent_dbase%EOSPropGrad(T,P,EOS_ENTHALPY,H,dH_dT,dH_dP,ierr)

  !PO todo: conversion when loaidng database to do this operation only once
  ! conversion to molar energy
  ! J/kg * kg/Kmol = J/Kmol
  H = H  * fmw_oil
  dH_dT = dH_dT * fmw_oil
  dH_dP = dH_dP * fmw_oil

  ! initialize derivative to NaN so that not mistakenly used.
  ! dH_dT = InitToNan()
  ! dH_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilEnthalpyEntDBase - H derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilEnthalpyEntDBase

! ************************************************************************** !

subroutine EOSOilEnthalpyNoDerive(T,P,H,ierr)
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dum1, dum2

  call EOSOilEnthalpyPtr(T,P,PETSC_FALSE,H,dum1,dum2,ierr)

end subroutine EOSOilEnthalpyNoDerive

! ************************************************************************** !

subroutine EOSOilEnthalpyDerive(T,P,H,dH_dT,dH_dP,ierr)
  
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! enthalpy deriv WRT Temperature [J/kmol/C]
  PetscReal, intent(out) :: dH_dP   ! enthalpy deriv WRT Pressure [J/kmol/Pa]
  PetscErrorCode, intent(out) :: ierr

  call EOSOilEnthalpyPtr(T,P,PETSC_TRUE,H,dH_dT,dH_dP,ierr)

end subroutine EOSOilEnthalpyDerive

! ************************************************************************** !

subroutine EOSOilDensityEnergyS(T,P,deriv,Rho,dRho_dT,dRho_dP, &
                                H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr,table_idxs)
  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! deriv. internal energy wrt temperature
  PetscReal, intent(out) :: dU_dP   ! deriv. internal energy wrt pressure
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: Rho2

  call EOSOilDensityPtr(T,P,deriv,Rho,dRho_dT,dRho_dP,ierr,table_idxs)
  call EOSOilEnthalpyPtr(T,P,deriv,H,dH_dT,dH_dP,ierr)
  U = H - P/Rho

  if (deriv) then
    Rho2 = Rho*Rho
    dU_dP = dH_dP - 1.d0/Rho + dRho_dP*P/Rho2 !! maybe add checks if dRho_dP or P tiny so as to avoid rounding errors
    dU_dT = dH_dT + dRho_dT*P/Rho2 !!dP_dT = 0 because P and T independent.

    !! down the line will need 2nd order derivs:
    !! ddU_dPP = ddH_dPP + dRho_dP/Rho2 + dRho_dP/Rho2 + ddRho_dPP*P/Rho2 - 2.d0*dRho_dP*dRho_dP*P/Rho2/Rho
    !! ddU_dTT = ddH_dTT + ddRho_dTT*P/Rho2 - 2.d0*dRho_dT*dRho_dT*P/Rho2/Rho
  else
    dU_dT = InitToNan()
    dU_dP = InitToNan()
  end if   
end subroutine EOSOilDensityEnergyS

! ************************************************************************** !

subroutine EOSOilDenEnergyNoDerive(T,P,Rho,H,U,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dum1, dum2, dum3, dum4, dum5, dum6

  call EOSOilDensityEnergyPtr(T,P,PETSC_FALSE,Rho,dum1,dum2, &
                              H,dum3,dum4,U,dum5,dum6,ierr,table_idxs)


end subroutine EOSOilDenEnergyNoDerive

! **************************************************************************** !

subroutine EOSOilDenEnergyDerive(T,P,Rho,dRho_dT,dRho_dP,H,dH_dT,dH_dP, &
                                    U,dU_dT,dU_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! oil density deriv WRT Temp [kmol/m^3/C]
  PetscReal, intent(out) :: dRho_dP ! oil density deriv WRT Press [kmol/m^3/Pa]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! enthalpy deriv WRT Temp [J/kmol/C]
  PetscReal, intent(out) :: dH_dP   ! enthalpy deriv WRT Press [J/kmol/Pa]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscReal, intent(out) :: dU_dT   ! int. energy deriv WRT Temp [J/kmol/C]
  PetscReal, intent(out) :: dU_dP   ! int. energy deriv WRT Press [J/kmol/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call EOSOilDensityEnergyPtr(T,P,PETSC_TRUE,Rho,dRho_dT,dRho_dP, &
                              H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr,table_idxs)

end subroutine EOSOilDenEnergyDerive

! **************************************************************************** !

subroutine EOSOilRSTable(T,P,deriv,RS,dRS_dT,dRS_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T       ! temperature [C]
  PetscReal, intent(in) :: P       ! oil pressure [Pa]
  PetscBool, intent(in) :: deriv   ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: RS     ! Vg[sm3]/Vo[sm3] per unit of res. volume
  PetscReal, intent(out) :: dRS_dT ! derivative RS wrt table temperature [sm3/sm3/C]
  PetscReal, intent(out) :: dRS_dP ! derivative RS wrt table Pressure [sm3/sm3/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  !ierr initialised in EOSProp
  !call pvt_table%EOSProp(T,P,EOS_RS,Rs,table_idxs,ierr)
  
  call pvt_table%EOSPropGrad(T,P,EOS_RS,Rs,dRS_dT,dRS_dP,table_idxs,ierr)

  ! initialize derivative to NaN so that not mistakenly used.
  ! dRS_dT = InitToNan()
  ! dRS_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilRSTable - RS derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilRSTable

! **************************************************************************** !

subroutine EOSOilRSNoDerive(T,P,RS,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T       ! temperature [C]
  PetscReal, intent(in) :: P       ! oil pressure [Pa]
  PetscReal, intent(out) :: RS     ! Vg[sm3]/Vo[sm3] per unit of res. volume
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dum1, dum2

  call EOSOilRSPtr(T,P,PETSC_FALSE,RS,dum1,dum2,ierr,table_idxs)

end subroutine EOSOilRSNoDerive

! **************************************************************************** !

subroutine EOSOilRSDerive(T,P,RS,dRS_dT,dRS_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T       ! temperature [C]
  PetscReal, intent(in) :: P       ! oil pressure [Pa]
  PetscReal, intent(out) :: RS     ! Vg[sm3]/Vo[sm3] per unit of res. volume
  PetscReal, intent(out) :: dRS_dT ! RS deriv WRT Temp [sm3/sm3/C]
  PetscReal, intent(out) :: dRS_dP ! RS deriv WRT Press [sm3/sm3/Pa]
  
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call EOSOilRSPtr(T,P,PETSC_TRUE,RS,dRS_dT,dRS_dP,ierr,table_idxs)

end subroutine EOSOilRSDerive

! **************************************************************************** !

subroutine EOSOilCompressibilityTable(T,P,deriv,Co,dCo_dT,dCo_dP,ierr, &
                                      table_idxs)
  implicit none
  PetscReal, intent(in) :: T       ! temperature [C]
  PetscReal, intent(in) :: P       ! oil pressure [Pa]
  PetscBool, intent(in) :: deriv   ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Co     ! oil compressibility [1/Pa]
  PetscReal, intent(out) :: dCo_dT ! derivative Co wrt table temperature [1/Pa] * [1/C]
  PetscReal, intent(out) :: dCo_dP ! derivative Co wrt table Pressure [1/Pa] * [1/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  !ierr initialised in EOSProp
  !call pvt_table%EOSProp(T,P,EOS_COMPRESSIBILITY,Co,table_idxs,ierr)

  call pvt_table%EOSPropGrad(T,P,EOS_COMPRESSIBILITY,Co,dCo_dT,dCo_dP, &
                             table_idxs,ierr) 
   
  ! ! initialize derivative to NaN so that not mistakenly used.
  ! dCo_dT = InitToNan()
  ! dCo_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilCompressibilityTable - Co derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilCompressibilityTable

! **************************************************************************** !

subroutine EOSOilCompressibilityNoDerive(T,P,Co,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T       ! temperature [C]
  PetscReal, intent(in) :: P       ! oil pressure [Pa]
  PetscReal, intent(out) :: Co     ! oil compressibility [1/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dum1, dum2

  call EOSOilCompressibilityPtr(T,P,PETSC_FALSE,Co,dum1,dum2,ierr,table_idxs)

end subroutine EOSOilCompressibilityNoDerive

! **************************************************************************** !

subroutine EOSOilCompressibilityDerive(T,P,Co,dCo_dT,dCo_dP,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T       ! temperature [C]
  PetscReal, intent(in) :: P       ! oil pressure [Pa]
  PetscReal, intent(out) :: Co     ! oil compressibility [1/Pa]
  PetscReal, intent(out) :: dCo_dT ! oil compress. derive WRT Temp [1/Pa/C]
  PetscReal, intent(out) :: dCo_dP ! oil compress. derive WRT Temp [1/Pa/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call EOSOilCompressibilityPtr(T,P,PETSC_TRUE,Co,dCo_dT,dCo_dP,ierr, &
                                table_idxs)

end subroutine EOSOilCompressibilityDerive

! **************************************************************************** !

subroutine EOSOilViscosibilityTable(T,P,deriv,Cvis,dCvis_dT,dCvis_dP,ierr, &
                                      table_idxs)
  implicit none
  PetscReal, intent(in) :: T       ! temperature [C]
  PetscReal, intent(in) :: P       ! oil pressure [Pa]
  PetscBool, intent(in) :: deriv   ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: Cvis     ! oil viscosibility [1/Pa]
  PetscReal, intent(out) :: dCvis_dT ! derivative Cvis wrt table temperature [1/Pa/C]
  PetscReal, intent(out) :: dCvis_dP ! derivative Cvis wrt table Pressure [1/Pa/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  !ierr initialised in EOSProp
  !call pvt_table%EOSProp(T,P,EOS_VISCOSIBILITY,Cvis,table_idxs,ierr)
 
  call pvt_table%EOSPropGrad(T,P,EOS_VISCOSIBILITY,Cvis,dCvis_dT,dCvis_dP, &
                        table_idxs,ierr)

  ! ! initialize derivative to NaN so that not mistakenly used.
  ! dCvis_dT = InitToNan()
  ! dCvis_dP = InitToNan()
  ! 
  ! if (deriv) then
  !   ! not yet implemented
  !   ierr = 99 !error 99 points out that deriv are asked but not available yet.
  !   print*, "EOSOilCompressibilityTable - Co derivatives not supported"
  !   stop
  ! end if

end subroutine EOSOilViscosibilityTable

! **************************************************************************** !

subroutine EOSOilViscosibilityNoDerive(T,P,Cvis,ierr,table_idxs)

  implicit none

  PetscReal, intent(in) :: T       ! temperature [C]
  PetscReal, intent(in) :: P       ! oil pressure [Pa]
  PetscReal, intent(out) :: Cvis   ! oil viscosibility [1/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: dum1, dum2

  call EOSOilViscosibilityPtr(T,P,PETSC_FALSE,Cvis,dum1,dum2,ierr,table_idxs)

end subroutine EOSOilViscosibilityNoDerive

! **************************************************************************** !

subroutine EOSOilViscosibilityDerive(T,P,Cvis,dCvis_dT,dCvis_dP,ierr, &
                                     table_idxs)

  implicit none

  PetscReal, intent(in) :: T       ! temperature [C]
  PetscReal, intent(in) :: P       ! oil pressure [Pa]
  PetscReal, intent(out) :: Cvis   ! oil viscosibility [1/Pa]
  PetscReal, intent(out) :: dCvis_dT ! oil visc. deriv WRT Temp [1/Pa/C]
  PetscReal, intent(out) :: dCvis_dP ! oil visc. deriv WRT Press [1/Pa/Pa]
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  call EOSOilViscosibilityPtr(T,P,PETSC_TRUE,Cvis,dCvis_dT,dCvis_dP, &
                              ierr,table_idxs)

end subroutine EOSOilViscosibilityDerive

! **************************************************************************** !


function InitToNan()

implicit none

PetscReal :: InitToNan

InitToNan = 0.0
InitToNan = 1.0/InitToNan
InitToNan = 0.0d0*InitToNan

return

end function InitToNan

! **************************************************************************** !

subroutine EOSOilSetPVDO(input,option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/18/17
  !
  ! Set up a PVDO table

  use Option_module
  use Input_Aux_module
  use Lookup_Table_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: db_var => null()
  character(len=MAXWORDLENGTH) :: internal_units, user_units
  PetscInt :: data_idx

  pvt_table => EOSTableCreate('PVDO',option)
  
  pvt_table%num_prop = 2
  
  ! units initially assing default values - overwritten by units specified 
  ! in the table input 
  internal_units = '' !assign default value by SetDefaultInternalUnits
  user_units = ''     !assign default value by SetMetricUnits
  
  !adding FVF 
  data_idx = 1 !position of FVF in the table (after pressure)
  db_var => CreateLookupTableVar(EOS_FVF,internal_units,user_units,data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !adding VISCOSITY 
  data_idx = 2 !position of FVF in the table (after pressure)
  db_var => CreateLookupTableVar(EOS_VISCOSITY,internal_units,user_units, &
                                 data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !set Default internal must be called before Set Metric
  call pvt_table%SetDefaultInternalUnits(option)
  call pvt_table%SetMetricUnits(option)

  call pvt_table%Read(input,option)

  call EOSTableAddToList(pvt_table,eos_table_list)

  EOSOilViscosityPtr => EOSOilViscosityTable
  EOSOilDensityPtr => EOSOilDensityTable

end subroutine EOSOilSetPVDO

! **************************************************************************** !

subroutine EOSOilSetPVCO(input,option)
  !
  ! Author: Paolo Orsini
  ! Date: 10/31/17
  !
  ! Set up a PVCO table

  use Option_module
  use Input_Aux_module
  use Lookup_Table_module

  implicit none

  type(input_type), pointer :: input
  type(option_type) :: option

  type(lookup_table_var_type), pointer :: db_var => null()
  character(len=MAXWORDLENGTH) :: internal_units, user_units
  PetscInt :: data_idx

  pvt_table => EOSTableCreate('PVCO',option)


  pvt_table%num_prop = 5
  
  ! units initially assing default values - overwritten by units specified 
  ! in the table input 
  internal_units = '' !assign default value by SetDefaultInternalUnits
  user_units = ''     !assign default value by SetMetricUnits
  
  !adding RS 
  data_idx = 1 !position in the table after pressure
  db_var => CreateLookupTableVar(EOS_RS,internal_units,user_units,data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !adding FVF
  data_idx = 2 !position in the table after pressure
  db_var => CreateLookupTableVar(EOS_FVF,internal_units,user_units,data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !adding VISCOSITY
  data_idx = 3 !position in the table after pressure
  db_var => CreateLookupTableVar(EOS_VISCOSITY,internal_units,user_units, &
                                 data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !adding COMPRESSIBILITY
  data_idx = 4 !position in the table after pressure
  db_var => CreateLookupTableVar(EOS_COMPRESSIBILITY,internal_units, &
                                 user_units,data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !adding VISCOSIBILITY
  data_idx = 5 !position in the table after pressure
  db_var => CreateLookupTableVar(EOS_VISCOSIBILITY,internal_units,user_units, &
                                 data_idx)
  call pvt_table%AddEOSProp(db_var,option)
  nullify(db_var)

  !set Default internal must be called before Set Metric
  call pvt_table%SetDefaultInternalUnits(option)
  call pvt_table%SetMetricUnits(option)

  call pvt_table%Read(input,option)

  call EOSTableAddToList(pvt_table,eos_table_list)

  EOSOilViscosityPtr => EOSOilViscosityTable
  EOSOilDensityPtr => EOSOilDensityTable
  EOSOilRSPtr => EOSOilRSTable
  EOSOilCompressibilityPtr => EOSOilCompressibilityTable
  EOSOilViscosibilityPtr => EOSOilViscosibilityTable

end subroutine EOSOilSetPVCO

! ************************************************************************** !

subroutine EOSOilTableProcess(option,FMW_gas,ref_den_gas_kg)
  !
  ! Author: Paolo Orsini
  ! Date: 10/28/17
  !
  ! Processes oil pvt table - once the entire input deck has been read

  use Option_module

  implicit none

  type(option_type) :: option
  PetscReal, intent(in) :: FMW_gas
  PetscReal, intent(in) :: ref_den_gas_kg

  if (.not.associated(pvt_table)) return

  select case(pvt_table%name)
    case("PVDO","PVCO")
      call pvt_table%ConvertFVFtoMolarDensity(fmw_oil,reference_density_kg)
  end select

  select case(pvt_table%name)
    case("PVCO") !will add PVTO when available
      call ConvertRSVoltoRSMolar(FMW_gas,ref_den_gas_kg)
  end select

end subroutine EOSOilTableProcess

! ************************************************************************** !

subroutine ConvertRSVoltoRSMolar(FMW_gas,ref_den_gas_kg)
  !
  ! Author: Paolo Orsini
  ! Date: 11/08/17
  !
  ! Convert volumetric RS to molar RS

  use Lookup_Table_module

  implicit none

  PetscReal, intent(in) :: FMW_gas
  PetscReal, intent(in) :: ref_den_gas_kg

  PetscInt :: data_idx
  PetscReal, pointer :: var_data(:,:) => null()
  type(lookup_table_var_ptr_type), pointer :: var_array(:) => null()  

  var_array => pvt_table%lookup_table_gen%var_array
  var_data => pvt_table%lookup_table_gen%var_data


  data_idx = var_array(EOS_RS)%ptr%data_idx
  !mol/mol = (kg/sm3 * kmol/kg)_gas / (kg/sm3 * kmol/kg)_oil * (sm3_g / sm3_o)
  var_data(data_idx,:) = &
           (ref_den_gas_kg / FMW_gas) / (reference_density_kg / fmw_oil) * &
            var_data(data_idx,:)
  
  !from this point on in the data_idx there is not FVF but EOS_DENSITY
  !modify variable pointig this
  var_array(EOS_RS)%ptr%internal_units = 'mol/mol'
  var_array(EOS_RS)%ptr%user_units = 'mol/mol'
  var_array(EOS_RS)%ptr%conversion_factor = 1.0

  nullify(var_array)
  nullify(var_data)

end subroutine ConvertRSVoltoRSMolar

! ************************************************************************** !

subroutine EOSOilInputRecord()
  !
  ! Prints ingested equation of state information to the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 05/04/2016
  !

  implicit none

  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'OIL'

  write(id,'(a)') 'EOSOilInputRecord not implemented: &
                  &Email jmfrede@sandia.gov for more information if using &
                  &OIL modes and EOS information is wanted.'

  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'

end subroutine EOSOilInputRecord

! ************************************************************************** !
subroutine EOSOilDBaseDestroy()

  implicit none

  !destroy databases
  call EOSDatabaseDestroy(eos_dbase)
  call EOSDatabaseDestroy(eos_den_dbase)
  call EOSDatabaseDestroy(eos_ent_dbase)
  call EOSDatabaseDestroy(eos_vis_dbase)

  !nullify EOS table pointer - pvt tables are deallocated when destroying
  !the table list
  nullify(pvt_table)

end subroutine EOSOilDBaseDestroy

! ************************************************************************** !

end module EOS_Oil_module
