module EOS_Oil_module
 
  use PFLOTRAN_Constants_module
  use EOSDatabase_module

  implicit none

  private
  
#include "petsc/finclude/petscsys.h"

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

  ! these should be define as astract interfaces, because there are no   
  ! precedures named as xxxDummy that have such interfaces
  interface
    subroutine EOSOilViscosityDummy(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! oil pressure [Pa]
      PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]  
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Vis     ! oil viscosity 
      PetscReal, intent(out) :: dVis_dT ! derivative oil viscosity wrt temperature
      PetscReal, intent(out) :: dVis_dP ! derivative oil viscosity wrt Pressure
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSOilViscosityDummy
    subroutine EOSOilDensityDummy(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr)
      implicit none
      PetscReal, intent(in) :: T        ! temperature [C]
      PetscReal, intent(in) :: P        ! pressure [Pa]
      PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
      PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
      PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
      PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
      PetscErrorCode, intent(out) :: ierr
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
                                        H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
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
    end subroutine EOSOilDensityEnergyDummy
  end interface 

  ! interfaces for derivative/non-derivative versions that are visible outside
  ! the module.
  interface EOSOilViscosity
    procedure EOSOilViscosityNoDerive
  ! procedure EOSOilViscosityDerive
  end interface
  interface EOSOilDensity
    procedure EOSOilDensityNoDerive
    procedure EOSOilDensityDerive
  end interface
  interface EOSOilEnthalpy
    procedure EOSOilEnthalpyNoDerive
   ! procedure EOSOilEnthalpyDerive
  end interface
  interface EOSOilDensityEnergy
    procedure EOSOilDenEnergyNoDerive
   ! procedure EOSOilDenEnergyDerive
  end interface  


  public :: EOSOilInit, &
            EOSOilVerify, &
            EOSOilViscosity, &
            EOSOilDensity, &
            EOSOilEnthalpy, & 
            EOSOilDensityEnergy, &
            EOSOilInputRecord

  public :: EOSOilSetFMWConstant, &
            EOSOilGetFMW, &
            EOSOilSetViscosityConstant, &
            EOSOilSetViscosityQuad, &
            EOSOilSetVisQuadRefVis, &
            EOSOilSetVisQuadRefPres, &
            EOSOilSetVisQuadRefTemp, &
            EOSOilSetVisQuadPresCoef, &
            EOSOilSetVisQuadTempCoef, &  
            EOSOilSetVisDBase, &
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

  compress_coeff = UNINITIALIZED_DOUBLE  
  th_expansion_coeff = UNINITIALIZED_DOUBLE 
  den_linear_den0 = UNINITIALIZED_DOUBLE
  den_linear_ref_pres = UNINITIALIZED_DOUBLE
  den_linear_ref_temp = UNINITIALIZED_DOUBLE

  quad_ent_ref_temp(1:2) = UNINITIALIZED_DOUBLE 
  quad_ent_temp_coef(1:2) = UNINITIALIZED_DOUBLE 


  fmw_oil = FMWOIL !default oil formula weight C10H22 (142 g/mol)

  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms

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

subroutine EOSOilSetDensityConstant(density)

  implicit none
  
  PetscReal :: density
  
  constant_density = density  
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms
  EOSOilDensityPtr => EOSOilDensityConstant
  
end subroutine EOSOilSetDensityConstant

! ************************************************************************** !

subroutine EOSOilSetDensityLinear()

  implicit none
  
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms
  EOSOilDensityPtr => EOSOilDensityLinear
  
end subroutine EOSOilSetDensityLinear

! ************************************************************************** !

subroutine EOSOilSetDensityInverseLinear()

  implicit none
  
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms
  EOSOilDensityPtr => EOSOilDensityInverseLinear
  
end subroutine EOSOilSetDensityInverseLinear


! ************************************************************************** !

subroutine EOSOilSetDenLinearRefDen(den0)

  implicit none
  
  PetscReal :: den0

  den_linear_den0 = den0 
  
end subroutine EOSOilSetDenLinearRefDen

! ************************************************************************** !

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
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms 
  EOSOilDensityPtr => EOSOilDensityDenDBase

end subroutine EOSOilSetDenDBase

! ************************************************************************** !

subroutine EOSOilSetEnthalpyConstant(enthalpy)

  implicit none
  
  PetscReal :: enthalpy
  
  constant_enthalpy = enthalpy  
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms
  EOSOilEnthalpyPtr => EOSOilEnthalpyConstant
  
end subroutine EOSOilSetEnthalpyConstant

! ************************************************************************** !

subroutine EOSOilSetEnthalpyLinearTemp(specific_heat)

  implicit none
  
  PetscReal :: specific_heat 
  
  constant_sp_heat = specific_heat  
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms
  EOSOilEnthalpyPtr => EOSOilEnthalpyLinearTemp

  !write(*,*) "I am in EOS oil linear set up"  

end subroutine EOSOilSetEnthalpyLinearTemp

! ************************************************************************** !

subroutine EOSOilSetEnthalpyQuadraticTemp()

  implicit none
  
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms
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
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms 
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
  EOSOilDensityEnergyPtr => EOSOilDensityEnergyTOilIms 
  if (.not.associated(EOSOilDensityPtr))  &
            EOSOilDensityPtr => EOSOilDensityEOSDBase
  if (.not.associated(EOSOilEnthalpyPtr)) &
    EOSOilEnthalpyPtr => EOSOilEnthalpyEOSDBase
  if (.not.associated(EOSOilViscosityPtr)) &
    EOSOilViscosityPtr => EOSOilViscosityEOSDBase

end subroutine EOSOilSetEOSDBase

! ************************************************************************** !

subroutine EOSOilViscosityConstant(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr)
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

  Vis = constant_viscosity

  dVis_dT = 0.0d0
  dVis_dP = 0.0d0
  
end subroutine EOSOilViscosityConstant

! ************************************************************************** !

subroutine EOSOilQuadViscosity(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr)
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

subroutine EOSOilViscosityEOSDBase(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr)
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

  !ierr initialised in EOSEOSProp 
  call eos_dbase%EOSProp(T,P,EOS_VISCOSITY,Vis,ierr)

  dVis_dT = 0.0d0
  dVis_dP = 0.0d0

  if (deriv) then
    ! not yet implemented
    ierr = 99 !error 99 points out that deriv are asked but not available yet. 
    print*, "EOSOilViscosityEOSDBase - Viscosity derivatives not supported"
    stop
  end if
  
end subroutine EOSOilViscosityEOSDBase

! ************************************************************************** !

subroutine EOSOilViscosityVisDBase(T,P,Rho,deriv,Vis,dVis_dT,dVis_dP,ierr)
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

  !ierr initialised in EOSEOSProp 
  call eos_vis_dbase%EOSProp(T,P,EOS_VISCOSITY,Vis,ierr)

  dVis_dT = 0.0d0
  dVis_dP = 0.0d0

  if (deriv) then
    ! not yet implemented
    ierr = 99 !error 99 points out that deriv are asked but not available yet. 
    print*, "EOSOilViscosityVisDBase - Viscosity derivatives not supported"
    stop
  end if
  
end subroutine EOSOilViscosityVisDBase

! ************************************************************************** !

subroutine EOSOilViscosityNoDerive(T,P,Rho,Vis,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! oil pressure [Pa]
  PetscReal, intent(in) :: Rho      ! oil density [kmol/m3]  
  PetscReal, intent(out) :: Vis     ! oil viscosity 
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dum1, dum2
  
  call EOSOilViscosityPtr(T,P,Rho,PETSC_FALSE,Vis,dum1,dum2,ierr)
  
end subroutine EOSOilViscosityNoDerive

! ************************************************************************** !

subroutine EOSOilDensityConstant(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr)
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
        ! kg/m3 * kmol/kg  = kmol/m3
  Rho = constant_density / fmw_oil ! kmol/m^3

  dRho_dT = 0.d0
  dRho_dP = 0.d0

end subroutine EOSOilDensityConstant

! ************************************************************************** !

subroutine EOSOilDensityLinear(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr)
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

  Rho = den_linear_den0 + &
        compress_coeff * (P - den_linear_ref_pres ) - & ! compression 
        th_expansion_coeff * (T - den_linear_ref_temp )    ! expansion

  ! conversion to molar density
        ! kg/m3 * kmol/kg  = kmol/m3
  Rho = Rho / fmw_oil ! kmol/m^3

  if (deriv) then
    dRho_dT = compress_coeff / fmw_oil
    dRho_dP = - th_expansion_coeff / fmw_oil
  end if

end subroutine EOSOilDensityLinear

! ************************************************************************** !

subroutine EOSOilDensityInverseLinear(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr)
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

subroutine EOSOilDensityEOSDBase(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr)
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

  !ierr initialised in EOSEOSProp 
  call eos_dbase%EOSProp(T,P,EOS_DENSITY,Rho,ierr)

  ! conversion to molar density
        ! kg/m3 * kmol/kg  = kmol/m3
  Rho = Rho / fmw_oil ! kmol/m^3

  if (deriv) then
    ! not yet implemented
    ierr = 99 !error 99 points out that deriv are asked but not available yet. 
    print*, "EOSOilDensityEOSDBase - Den derivatives not supported"
    stop
  end if

end subroutine EOSOilDensityEOSDBase

! ************************************************************************** !

subroutine EOSOilDensityDenDBase(T, P, deriv, Rho, dRho_dT, dRho_dP, ierr)
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

  !ierr initialised in EOSEOSProp 
  call eos_den_dbase%EOSProp(T,P,EOS_DENSITY,Rho,ierr)

  ! conversion to molar density
        ! kg/m3 * kmol/kg  = kmol/m3
  Rho = Rho / fmw_oil ! kmol/m^3

  if (deriv) then
    ! not yet implemented
    ierr = 99 !error 99 points out that deriv are asked but not available yet. 
    print*, "EOSOilDensityDenDBase - Den derivatives not supported"
    stop
  end if

end subroutine EOSOilDensityDenDBase

! ************************************************************************** !

subroutine EOSOilDensityNoDerive(T,P,Rho,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2
  
  ! derivatives are so cheap, just compute them
  call EOSOilDensityPtr(T, P, PETSC_FALSE, Rho, dum1, dum2, ierr)
  
end subroutine EOSOilDensityNoDerive

! ************************************************************************** !

subroutine EOSOilDensityDerive(T,P,Rho,dRho_dT,dRho_dP,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: dRho_dT ! derivative oil density wrt temperature
  PetscReal, intent(out) :: dRho_dP ! derivative oil density wrt pressure
  PetscErrorCode, intent(out) :: ierr
  
  call EOSOilDensityPtr(T,P,PETSC_TRUE,Rho,dRho_dT,dRho_dP,ierr)
                          
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
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  !ierr initialised in EOSEOSProp 
  call eos_dbase%EOSProp(T,P,EOS_ENTHALPY,H,ierr)


  ! conversion to molar energy
  ! J/kg * kg/Kmol = J/Kmol
  H = H  * fmw_oil  

  dH_dT = UNINITIALIZED_DOUBLE
  dH_dP = UNINITIALIZED_DOUBLE

  if (deriv) then
    ! not yet implemented
    ierr = 99 !error 99 points out that deriv are asked but not available yet. 
    print*, "EOSOilEnthalpyEOSDBase - H derivatives not supported"
    stop  
  end if

end subroutine EOSOilEnthalpyEOSDBase

! ************************************************************************** !

subroutine EOSOilEnthalpyEntDBase(T,P,deriv,H,dH_dT,dH_dP,ierr)
  implicit none
  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscBool, intent(in) :: deriv    ! indicate if derivatives are needed or not
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: dH_dT   ! derivative enthalpy wrt temperature
  PetscReal, intent(out) :: dH_dP   ! derivative enthalpy wrt pressure
  PetscErrorCode, intent(out) :: ierr

  !ierr initialised in EOSEOSProp 
  call eos_ent_dbase%EOSProp(T,P,EOS_ENTHALPY,H,ierr)


  ! conversion to molar energy
  ! J/kg * kg/Kmol = J/Kmol
  H = H  * fmw_oil  

  dH_dT = UNINITIALIZED_DOUBLE
  dH_dP = UNINITIALIZED_DOUBLE

  if (deriv) then
    ! not yet implemented
    ierr = 99 !error 99 points out that deriv are asked but not available yet. 
    print*, "EOSOilEnthalpyEntDBase - H derivatives not supported"
    stop  
  end if

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

subroutine EOSOilDensityEnergyTOilIms(T,P,deriv,Rho,dRho_dT,dRho_dP, &
                                    H,dH_dT,dH_dP,U,dU_dT,dU_dP,ierr)
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

  call EOSOilDensityPtr(T,P,deriv,Rho,dRho_dT,dRho_dP,ierr)
  call EOSOilEnthalpyPtr(T,P,deriv,H,dH_dT,dH_dP,ierr)  

  U = H - P/Rho

  dU_dT = UNINITIALIZED_DOUBLE
  dU_dP = UNINITIALIZED_DOUBLE

  if (deriv) then
    print*, "EOSOilDensityEnergyTOilIms - U derivatives not supported"
    stop  
  end if   


end subroutine EOSOilDensityEnergyTOilIms

! ************************************************************************** !

subroutine EOSOilDenEnergyNoDerive(T,P,Rho,H,U,ierr)

  implicit none

  PetscReal, intent(in) :: T        ! temperature [C]
  PetscReal, intent(in) :: P        ! pressure [Pa]
  PetscReal, intent(out) :: Rho     ! oil density [kmol/m^3]
  PetscReal, intent(out) :: H       ! enthalpy [J/kmol]
  PetscReal, intent(out) :: U       ! internal energy [J/kmol]
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dum1, dum2, dum3, dum4, dum5, dum6
 
  call EOSOilDensityEnergyPtr(T,P,PETSC_FALSE,Rho,dum1,dum2, &
                              H,dum3,dum4,U,dum5,dum6,ierr) 
 

end subroutine EOSOilDenEnergyNoDerive

! **************************************************************************** !

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
  
  call EOSDatabaseDestroy(eos_dbase)
  call EOSDatabaseDestroy(eos_den_dbase)
  call EOSDatabaseDestroy(eos_ent_dbase)
  call EOSDatabaseDestroy(eos_vis_dbase)

end subroutine EOSOilDBaseDestroy

! ************************************************************************** !

end module EOS_Oil_module
