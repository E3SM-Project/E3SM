module EOS_Water_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Geometry_module

  implicit none

  private
  
  ! module variables
  PetscReal :: constant_density
  PetscReal :: constant_enthalpy
  PetscReal :: constant_viscosity
  PetscReal :: constant_steam_density
  PetscReal :: constant_steam_enthalpy

  ! exponential
  PetscReal :: exponent_reference_density
  PetscReal :: exponent_reference_pressure
  PetscReal :: exponent_water_compressibility
  
  ! planes for planar eos
  type(plane_type) :: water_density_tp_plane
  type(plane_type) :: water_enthalpy_tp_plane
  type(plane_type) :: steam_density_tp_plane
  type(plane_type) :: steam_enthalpy_tp_plane
  

  ! quadratic
  PetscReal :: quadratic_reference_density
  PetscReal :: quadratic_reference_pressure
  PetscReal :: quadratic_wat_compressibility

  ! linear
  PetscReal :: linear_reference_density
  PetscReal :: linear_reference_pressure
  PetscReal :: linear_water_compressibility

  ! In order to support generic EOS subroutines, we need the following:
  ! 1. An interface declaration that defines the argument list (best to have 
  !    "Dummy" appended.
  ! 2. A procedure pointer that is initially set to null.  This pointer is
  !    pointed to the appropriate subroutine later on (e.g. EOSWaterInit())
  ! 3. An interface for derivative/non-derivative versions

  ! procedure pointer declarations
  ! standard versions
  procedure(EOSWaterViscosityDummy), pointer :: EOSWaterViscosityPtr => null()
  procedure(EOSWaterSatPressDummy), pointer :: &
    EOSWaterSaturationPressurePtr => null()
  procedure(EOSWaterDensityDummy), pointer :: EOSWaterDensityPtr => null()
  procedure(EOSWaterEnthalpyDummy), pointer :: EOSWaterEnthalpyPtr => null()
  procedure(EOSWaterSteamDenEnthDummy), pointer :: &
    EOSWaterSteamDensityEnthalpyPtr => null()
  procedure(EOSWaterDensityIceDummy), pointer :: &
    EOSWaterDensityIcePtr => null()
  ! extended versions
  procedure(EOSWaterViscosityExtDummy), pointer :: &
    EOSWaterViscosityExtPtr => null()
  procedure(EOSWaterDensityExtDummy), pointer :: &
    EOSWaterDensityExtPtr => null()
  procedure(EOSWaterEnthalpyExtDummy), pointer :: &
    EOSWaterEnthalpyExtPtr => null()
    
  ! interface blocks
  interface
  ! standard versions
    subroutine EOSWaterViscosityDummy(T, P, PS, dPS_dT, &
                                      calculate_derivatives, VW, &
                                      dVW_dT, dVW_dP, ierr)
      implicit none
      PetscReal, intent(in) :: T, P, PS, dPS_dT
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: VW
      PetscReal, intent(out) :: dVW_dT, dVW_dP
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterViscosityDummy
    subroutine EOSWaterSatPressDummy(T, calculate_derivatives, &
                                     PS, dPS_dT, ierr)
      implicit none
      PetscReal, intent(in) :: T
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: PS, dPS_dT
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterSatPressDummy
    subroutine EOSWaterDensityDummy(t,p,calculate_derivatives, &
                                    dw,dwmol,dwp,dwt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dw,dwmol,dwp,dwt
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterDensityDummy
    subroutine EOSWaterEnthalpyDummy(t,p,calculate_derivatives, &
                                     hw,hwp,hwt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: hw,hwp,hwt
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterEnthalpyDummy
    subroutine EOSWaterSteamDenEnthDummy(t,p,calculate_derivatives, &
                                         dg,dgmol,hg,dgp,dgt,hgp,hgt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dg,dgmol,dgp,dgt
      PetscReal, intent(out) :: hg,hgp,hgt
      PetscErrorCode, intent(out) :: ierr 
    end subroutine EOSWaterSteamDenEnthDummy
    subroutine EOSWaterDensityIceDummy(t,p,calculate_derivatives, &
                                       dw,dwp,dwt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dw,dwp,dwt
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterDensityIceDummy
    ! Extended versions
    subroutine EOSWaterViscosityExtDummy(T, P, PS, dPS_dT, aux, &
                                         calculate_derivatives, VW, &
                                         dVW_dT, dVW_dP, ierr)
      implicit none
      PetscReal, intent(in) :: T, P, PS, dPS_dT, aux(*)
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: VW
      PetscReal, intent(out) :: dVW_dT, dVW_dP
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterViscosityExtDummy
    subroutine EOSWaterDensityExtDummy(t,p,aux,calculate_derivatives, &
                                       dw,dwmol,dwp,dwt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscReal, intent(in) :: aux(*)
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: dw,dwmol,dwp,dwt
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterDensityExtDummy
    subroutine EOSWaterEnthalpyExtDummy(t,p,aux,calculate_derivatives, &
                                        hw,hwp,hwt,ierr)
      implicit none
      PetscReal, intent(in) :: t
      PetscReal, intent(in) :: p
      PetscReal, intent(in) :: aux(*)
      PetscBool, intent(in) :: calculate_derivatives
      PetscReal, intent(out) :: hw,hwp,hwt
      PetscErrorCode, intent(out) :: ierr
    end subroutine EOSWaterEnthalpyExtDummy
  end interface
  
  ! interfaces for derivative/non-derivative versions that are visible outside
  ! the module.
  ! standard versions
  interface EOSWaterViscosity
    procedure EOSWaterViscosityNoDerive
    procedure EOSWaterViscosityDerive
  end interface
  interface EOSWaterSaturationPressure
    procedure EOSWaterSatPresNoDerive
    procedure EOSWaterSatPresDerive
  end interface
  interface EOSWaterDensity
    procedure EOSWaterDensityNoDerive
    procedure EOSWaterDensityDerive
  end interface
  interface EOSWaterEnthalpy
    procedure EOSWaterEnthalpyNoDerive
    procedure EOSWaterEnthalpyDerive
  end interface
  interface EOSWaterSteamDensityEnthalpy
    procedure EOSWaterSteamDenEnthNoDerive
    procedure EOSWaterSteamDenEnthDerive
  end interface
  interface EOSWaterDensityIce
    procedure EOSWaterDensityIceNoDerive
    procedure EOSWaterDensityIceDerive
  end interface
  ! Extended versions
  interface EOSWaterViscosityExt
    procedure EOSWaterViscosityExtNoDerive
    procedure EOSWaterViscosityExtDerive
  end interface
  interface EOSWaterDensityExt
    procedure EOSWaterDensityExtNoDerive
    procedure EOSWaterDensityExtDerive
!geh: very useful for debuggin
!    procedure EOSWaterDensityExtNumericalDerive
  end interface
  interface EOSWaterEnthalpyExt
    procedure EOSWaterEnthalpyExtNoDerive
    procedure EOSWaterEnthalpyExtDerive
  end interface

  ! the "public" definition that makes subroutines visible outside.
  public :: EOSWaterInit, &
            EOSWaterVerify, &
            EOSWaterViscosity, &
            EOSWaterSaturationPressure, &
            EOSWaterDensity, &
            EOSWaterEnthalpy, &
            EOSWaterSteamDensityEnthalpy, &
            EOSWaterDuanMixture, &
            EOSWaterViscosityNaCl, &
            EOSWaterInternalEnergyIce, &
            EOSWaterDensityIcePainter, &
            EOSWaterSaturationTemperature, &
            EOSWaterDensityIce, &
            EOSWaterDensityTGDPB01, &
            EOSWaterViscosityExt, &
            EOSWaterDensityExt, &
            EOSWaterEnthalpyExt, &
            EOSWaterInputRecord

  public :: EOSWaterSetDensity, &
            EOSWaterSetEnthalpy, &
            EOSWaterSetViscosity, &
            EOSWaterSetSteamDensity, &
            EOSWaterSetSteamEnthalpy
            
            
  public :: TestEOSWaterBatzleAndWang, &
            EOSWaterTest, &
            EOSWaterSteamTest
 
  contains

! ************************************************************************** !

subroutine EOSWaterInit()

  implicit none
  
  constant_density = UNINITIALIZED_DOUBLE
  constant_viscosity = UNINITIALIZED_DOUBLE
  constant_enthalpy = UNINITIALIZED_DOUBLE
  constant_steam_density = UNINITIALIZED_DOUBLE
  constant_steam_enthalpy = UNINITIALIZED_DOUBLE
  exponent_reference_density = UNINITIALIZED_DOUBLE
  exponent_reference_pressure = UNINITIALIZED_DOUBLE
  exponent_water_compressibility = UNINITIALIZED_DOUBLE
  quadratic_reference_density = UNINITIALIZED_DOUBLE
  quadratic_reference_pressure = UNINITIALIZED_DOUBLE
  quadratic_wat_compressibility = UNINITIALIZED_DOUBLE
  
  ! standard versions  
  EOSWaterDensityPtr => EOSWaterDensityIFC67
  EOSWaterEnthalpyPtr => EOSWaterEnthalpyIFC67
  EOSWaterViscosityPtr => EOSWaterViscosity1
  EOSWaterSaturationPressurePtr => EOSWaterSaturationPressureIFC67
  EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDensityEnthalpyIFC67
  EOSWaterDensityIcePtr => EOSWaterDensityIcePainter
  
  ! extended versions
  EOSWaterViscosityExtPtr => EOSWaterViscosityKestinExt
  EOSWaterDensityExtPtr => EOSWaterDensityBatzleAndWangExt
  
end subroutine EOSWaterInit

! ************************************************************************** !

subroutine EOSWaterVerify(ierr,error_string)

  implicit none
  
  PetscErrorCode, intent(out) :: ierr
  character(len=MAXSTRINGLENGTH), intent(out) :: error_string
  
  ierr = 0
  error_string = ''
  if ((associated(EOSWaterDensityPtr,EOSWaterDensityIFC67) .and. &
        Initialized(constant_density)) .or. &
      (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyIFC67) .and. &
        Initialized(constant_enthalpy)) &
     ) then
    ierr = 1
  endif

  if (associated(EOSWaterDensityPtr,EOSWaterDensityConstant) .and. &
      Uninitialized(constant_density)) then
    error_string = trim(error_string) // &
      ' CONSTANT density not set.'
    ierr = 1
  endif
  
  if (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyConstant) .and. &
      Uninitialized(constant_enthalpy)) then
    error_string = trim(error_string) // &
      ' CONSTANT enthalpy not set.'
    ierr = 1
  endif
  
  if (associated(EOSWaterDensityPtr,EOSWaterDensityExponential) .and. &
      (Uninitialized(exponent_reference_density) .or. & 
       Uninitialized(exponent_reference_pressure) .or. &
       Uninitialized(exponent_water_compressibility))) then
    error_string = trim(error_string) // &
      ' Exponential parameters incorrect.'
    ierr = 1
  endif

  if (associated(EOSWaterDensityPtr,EOSWaterDensityLinear) .and. &
      (Uninitialized(linear_reference_density) .or. & 
       Uninitialized(linear_reference_pressure) .or. &
       Uninitialized(linear_water_compressibility))) then
    error_string = trim(error_string) // &
      ' Linear parameters incorrect.'
     ierr = 1
  endif
 
  if (associated(EOSWaterDensityPtr,EOSWaterDensityBRAGFLO) .and. &
      (Uninitialized(exponent_reference_density) .or. & 
       Uninitialized(exponent_reference_pressure) .or. &
       Uninitialized(exponent_water_compressibility))) then
    error_string = trim(error_string) // &
      ' BRAGFLO parameters incorrect.'
    ierr = 1
  endif

  if ((associated(EOSWaterViscosityPtr, &
                  EOSWaterViscosityConstant) .and. &
       Uninitialized(constant_viscosity)) .or. &
      (associated(EOSWaterViscosityPtr, &
                  EOSWaterViscosity1) .and. &
       Initialized(constant_viscosity))) then
    ierr = 1
  endif
  
  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDenEnthConstant) .and. &
      (Uninitialized(constant_steam_density) .or. &
       Uninitialized(constant_steam_enthalpy))) then
    if (Uninitialized(constant_steam_density)) then
      error_string = trim(error_string) // &
        ' CONSTANT steam density not set.'
    endif
    if (Uninitialized(constant_steam_enthalpy)) then
      error_string = trim(error_string) // &
        ' CONSTANT steam enthalpy not set.'
    endif
    ierr = 1
  endif
  
end subroutine EOSWaterVerify

! ************************************************************************** !

subroutine EOSWaterSetDensity(keyword,aux)

  implicit none
  
  character(len=*) :: keyword
  PetscReal, optional :: aux(*)
  
  select case(keyword)
    case('CONSTANT')
      constant_density = aux(1)  
      EOSWaterDensityPtr => EOSWaterDensityConstant
    case('DEFAULT','IFC67')
      EOSWaterDensityPtr => EOSWaterDensityIFC67
      EOSWaterDensityExtPtr => EOSWaterDensityBatzleAndWangExt
    case('EXPONENTIAL')
      exponent_reference_density = aux(1)
      exponent_reference_pressure = aux(2)
      exponent_water_compressibility = aux(3)  
      EOSWaterDensityPtr => EOSWaterDensityExponential
    case('LINEAR')
      linear_reference_density = aux(1)
      linear_reference_pressure = aux(2)
      linear_water_compressibility = aux(3)
      EOSWaterDensityPtr => EOSWaterDensityLinear
    case('BRAGFLO')
      exponent_reference_density = aux(1)
      exponent_reference_pressure = aux(2)
      exponent_water_compressibility = aux(3)  
      EOSWaterDensityPtr => EOSWaterDensityBRAGFLO
    case('QUADRATIC')
      if (Initialized(aux(1))) then
        quadratic_reference_density = 999.014d0 !kg/m3
      else
        quadratic_reference_density = aux(1)
      end if
      if (Initialized(aux(2))) then
        quadratic_reference_pressure = 1.0d5 !Pa
      else
        quadratic_reference_pressure = aux(2)
      end if
      if (Initialized(aux(3))) then
        quadratic_wat_compressibility = 3.94769306686405d-10 !1/Pa
      else
        quadratic_wat_compressibility = aux(3)
      end if
      EOSWaterDensityPtr => EOSWaterDensityQuadratic
    case('TRANGENSTEIN')
      EOSWaterDensityPtr => EOSWaterDensityTrangenstein
    case('PLANAR')
      EOSWaterDensityPtr => EOSWaterDensityTPPlanar
      call EOSWaterDensityTPPlanarSetup()
    case('TGDPB01')
      EOSWaterDensityPtr => EOSWaterDensityTGDPB01
    case('PAINTER')
      EOSWaterDensityPtr => EOSWaterDensityPainter
    case('BATZLE_AND_WANG')
      EOSWaterDensityPtr => EOSWaterDensityBatzleAndWang
      EOSWaterDensityExtPtr => EOSWaterDensityBatzleAndWangExt
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetDensity().'      
      stop
  end select
  
end subroutine EOSWaterSetDensity

! ************************************************************************** !

subroutine EOSWaterSetEnthalpy(keyword,aux)

  implicit none
  
  character(len=*) :: keyword
  PetscReal, optional :: aux(*)
  
  select case(keyword)
    case('CONSTANT')
      constant_enthalpy = aux(1)  
      EOSWaterEnthalpyPtr => EOSWaterEnthalpyConstant
    case('DEFAULT','IFC67')
      EOSWaterEnthalpyPtr => EOSWaterEnthalpyIFC67
    case('PLANAR')
      EOSWaterEnthalpyPtr => EOSWaterEnthalpyTPPlanar  
      call EOSWaterEnthalpyTPPlanarSetup()
    case('PAINTER')
      EOSWaterEnthalpyPtr => EOSWaterEnthalpyPainter      
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetEnthalpy().'
      stop
  end select
  
end subroutine EOSWaterSetEnthalpy

! ************************************************************************** !

subroutine EOSWaterSetViscosity(keyword,aux)

  implicit none
  
  character(len=*) :: keyword
  PetscReal, optional :: aux(*)
  
  select case(keyword)
    case('CONSTANT')
      constant_viscosity = aux(1)  
      EOSWaterViscosityPtr => EOSWaterViscosityConstant
      EOSWaterViscosityExtPtr => EOSWaterViscosityConstantExt
    case('DEFAULT')
      EOSWaterViscosityPtr => EOSWaterViscosity1
      EOSWaterViscosityExtPtr => EOSWaterViscosityKestinExt
    case('BATZLE_AND_WANG')
      EOSWaterViscosityPtr => EOSWaterViscosityBatzleAndWang
      EOSWaterViscosityExtPtr => EOSWaterViscosityBatzleAndWangExt
    case('GRABOWSKI')
      EOSWaterViscosityPtr => EOSWaterViscosityGrabowski
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetViscosity().'  
      stop 
  end select
  
end subroutine EOSWaterSetViscosity


! ************************************************************************** !

subroutine EOSWaterSetSteamDensity(keyword,aux)

  implicit none
  
  character(len=*) :: keyword
  PetscReal, optional :: aux(*)
  
  select case(keyword)
    case('CONSTANT')
      constant_steam_density = aux(1)  
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDenEnthConstant
    case('PLANAR')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDenEnthTPPlanar
      call EOSWaterSteamDenEnthTPPlanarSetup()
    case('IFC67')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDensityEnthalpyIFC67
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetSteamDensity().'
      stop
  end select
  
end subroutine EOSWaterSetSteamDensity

! ************************************************************************** !

subroutine EOSWaterSetSteamEnthalpy(keyword,aux)

  implicit none
  
  character(len=*) :: keyword
  PetscReal, optional :: aux(*)
  
  select case(keyword)
    case('CONSTANT')
      constant_steam_enthalpy = aux(1)  
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDenEnthConstant
    case('PLANAR')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDenEnthTPPlanar
      call EOSWaterSteamDenEnthTPPlanarSetup()
    case('DEFAULT','IFC67')
      EOSWaterSteamDensityEnthalpyPtr => EOSWaterSteamDensityEnthalpyIFC67
    case default
      print *, 'Unknown pointer type "' // trim(keyword) // &
        '" in EOSWaterSetSteamEnthalpy().'
      stop
  end select
  
end subroutine EOSWaterSetSteamEnthalpy

! ************************************************************************** !

subroutine EOSWaterViscosityNoDerive(T, P, PS, VW, ierr)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(out) :: VW ! water viscosity
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dPS_dT ! derivative of PS with respect to temp
  PetscReal :: dum1, dum2
  
  ierr = 0
  dPS_dT = 0.d0
  call EOSWaterViscosityPtr(T, P, PS, dPS_dT, PETSC_FALSE, VW, &
                            dum1, dum2, ierr)
  
end subroutine EOSWaterViscosityNoDerive

! ************************************************************************** !

subroutine EOSWaterViscosityDerive(T, P, PS, dPS_dT, VW, dVW_dT, &
                                   dVW_dP, ierr)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(out) :: ierr
  
  ierr = 0
  call EOSWaterViscosityPtr(T, P, PS, dPS_dT, PETSC_TRUE, VW, &
                            dVW_dT, dVW_dP, ierr)
  
end subroutine EOSWaterViscosityDerive

! ************************************************************************** !

subroutine EOSWaterSatPresNoDerive(T, PS, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscReal, intent(out) :: PS ! Saturation pres. and derivative
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dummy
  
  ierr = 0
  call EOSWaterSaturationPressurePtr(T, PETSC_FALSE, PS, dummy, ierr)
  
end subroutine EOSWaterSatPresNoDerive

! ************************************************************************** !

subroutine EOSWaterSatPresDerive(T, PS, dPS_dT, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(out) :: ierr
  
  ierr = 0
  call EOSWaterSaturationPressurePtr(T, PETSC_TRUE, PS, dPS_dT, ierr)
  
end subroutine EOSWaterSatPresDerive

! ************************************************************************** !

subroutine EOSWaterDensityNoDerive(t,p,dw,dwmol,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: dw,dwmol
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2
  
  ierr = 0
  call EOSWaterDensityPtr(t,p,PETSC_FALSE,dw,dwmol,dum1,dum2,ierr)
  
end subroutine EOSWaterDensityNoDerive

! ************************************************************************** !

subroutine EOSWaterDensityDerive(t,p,dw,dwmol,dwp,dwt,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr
  
  ierr = 0
  call EOSWaterDensityPtr(t,p,PETSC_TRUE,dw,dwmol,dwp,dwt,ierr)
  
end subroutine EOSWaterDensityDerive

! ************************************************************************** !

subroutine EOSWaterEnthalpyNoDerive(t,p,hw,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: hw
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2
  
  ierr = 0
  call EOSWaterEnthalpyPtr(t,p,PETSC_FALSE,hw,dum1,dum2,ierr)
  
end subroutine EOSWaterEnthalpyNoDerive

! ************************************************************************** !

subroutine EOSWaterEnthalpyDerive(t,p,hw,hwp,hwt,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(out) :: ierr
  
  ierr = 0
  call EOSWaterEnthalpyPtr(t,p,PETSC_TRUE,hw,hwp,hwt,ierr)
  
end subroutine EOSWaterEnthalpyDerive

! ************************************************************************** !

subroutine EOSWaterViscosityExtNoDerive(T, P, PS, aux, VW, ierr)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: VW ! water viscosity
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dPS_dT ! derivative of PS with respect to temp
  PetscReal :: dum1, dum2
  
  ierr = 0
  dPS_dT = 0.d0
  call EOSWaterViscosityExtPtr(T, P, PS, dPS_dT, aux, PETSC_FALSE, VW, &
                               dum1, dum2, ierr)
  
end subroutine EOSWaterViscosityExtNoDerive

! ************************************************************************** !

subroutine EOSWaterViscosityExtDerive(T, P, PS, dPS_dT, aux, VW, dVW_dT, &
                                      dVW_dP, ierr)

  implicit none

  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(out) :: ierr
  
  ierr = 0
  call EOSWaterViscosityExtPtr(T, P, PS, dPS_dT, aux, PETSC_TRUE, VW, &
                               dVW_dT, dVW_dP, ierr)
  
end subroutine EOSWaterViscosityExtDerive

! ************************************************************************** !

subroutine EOSWaterDensityExtNoDerive(t,p,aux,dw,dwmol,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: dw,dwmol
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2
  
  ierr = 0
  call EOSWaterDensityExtPtr(t,p,aux,PETSC_FALSE,dw,dwmol,dum1,dum2,ierr)
  
end subroutine EOSWaterDensityExtNoDerive

! ************************************************************************** !

subroutine EOSWaterDensityExtDerive(t,p,aux,dw,dwmol,dwp,dwt,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr
  
  ierr = 0
  call EOSWaterDensityExtPtr(t,p,aux,PETSC_TRUE,dw,dwmol,dwp,dwt,ierr)
  
end subroutine EOSWaterDensityExtDerive

! ************************************************************************** !

subroutine EOSWaterEnthalpyExtNoDerive(t,p,aux,hw,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: aux(*)  
  PetscReal, intent(out) :: hw
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2
  
  ierr = 0
  call EOSWaterEnthalpyExtPtr(t,p,aux,PETSC_FALSE,hw,dum1,dum2,ierr)
  
end subroutine EOSWaterEnthalpyExtNoDerive

! ************************************************************************** !

subroutine EOSWaterEnthalpyExtDerive(t,p,aux,hw,hwp,hwt,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(in) :: aux(*)  
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(out) :: ierr
  
  ierr = 0
  call EOSWaterEnthalpyExtPtr(t,p,aux,PETSC_TRUE,hw,hwp,hwt,ierr)
  
end subroutine EOSWaterEnthalpyExtDerive

! ************************************************************************** !

subroutine EOSWaterViscosity1(T, P, PS, dPS_dT, calculate_derivatives, &
                              VW, dVW_dT, dVW_dP, ierr)

! Calculates the viscosity of water and derivatives as a function of 
! temperature, pressure, and saturation pressure.

  implicit none
    
  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: EX, PHI, AM, pwr, aln10
  
  EX  = 247.8d0/(T+133.15d0)
  PHI = 1.0467d0*(T-31.85d0)
  !geh: added max(P,PS)-PS in place of P-PS
  !geh: here P should be the maximum of Pl and Pg
  AM  = 1.d0+PHI*(max(P,PS)-PS)*1.d-11
  pwr = 10.d0**EX
  VW = 1.d-7*AM*241.4d0*pwr
    
  if (calculate_derivatives) then
    aln10 = log(10.d0)
    dVW_dT = VW/AM*1.d-11* &
            ! dAM_PHI_dT       dAM_PS_dT
            (1.0467d0*(P-PS) - PHI*dPS_dT) - &
            ! dpwr_EX_dT
            VW*aln10*247.8d0/(T+133.15d0)**2
    dVW_dP = VW/AM*PHI*1.d-11
  else
    dVW_dT = UNINITIALIZED_DOUBLE
    dVW_dP = UNINITIALIZED_DOUBLE
  endif
 
end subroutine EOSWaterViscosity1

! ************************************************************************** !

subroutine EOSWaterViscosityConstant(T, P, PS, dPS_dT, &
                                     calculate_derivatives, &
                                      VW, dVW_dT, dVW_dP, ierr)

! Calculates the viscosity of water and derivatives as a function of 
! temperature, pressure, and saturation pressure.

  implicit none
    
  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(out) :: ierr
  
  VW = constant_viscosity
  
  dVW_dT = 0.d0
  dVW_dP = 0.d0
  
end subroutine EOSWaterViscosityConstant

! ************************************************************************** !

subroutine EOSWaterViscosityConstantExt(T, P, PS, dPS_dT, aux,&
                                        calculate_derivatives, &
                                        VW, dVW_dT, dVW_dP, ierr)

! Calculates the viscosity of water and derivatives as a function of 
! temperature, pressure, and saturation pressure.

  implicit none
    
  PetscReal, intent(in) :: T, P, PS ! temperature, pressure, saturation_press
  PetscReal, intent(in) :: dPS_dT ! derivative of PS with respect to temp
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW ! water viscosity
  PetscReal, intent(out) :: dVW_dT, dVW_dP ! derivatives
  PetscErrorCode, intent(out) :: ierr
  
  VW = constant_viscosity
  
  dVW_dT = 0.d0
  dVW_dP = 0.d0
  
end subroutine EOSWaterViscosityConstantExt

! ************************************************************************** !

subroutine EOSWaterSaturationPressureIFC67(T, calculate_derivatives, &
                                           PS, dPS_dT, ierr)

  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal, save, dimension(9) :: A(9)
  PetscReal :: TC, SC, PCAP, E1, E2
  PetscReal :: one_m_tc, one_m_tc_sq, E2_bottom
  PetscReal :: dTC_dT, dSC_dTC, dE1_dTC, dE2_dTC, dPC_dSC, dPC_dTC
    
  DATA A/ &
    -7.691234564d0,-2.608023696d1,-1.681706546d2,6.423285504d1, &
    -1.189646225d2,4.167117320d0,2.097506760E1,1.d9,6.d0/

  if (T .GT. 500.d0) then
    ierr = 1
    return
  end if
  TC = (T+273.15d0)/H2O_CRITICAL_TEMPERATURE
  one_m_tc = 1.d0-TC
  one_m_tc_sq = one_m_tc*one_m_tc
  SC = A(1)*one_m_tc+A(2)*one_m_tc_sq+A(3)*one_m_tc**3.d0+ &
        A(4)*one_m_tc**4.d0+A(5)*one_m_tc**5.d0
  E1 = TC*(1.d0+A(6)*one_m_tc+A(7)*one_m_tc_sq)
  E2_bottom = A(8)*one_m_tc_sq+A(9)
  E2 = one_m_tc/E2_bottom
  PCAP = EXP(SC/E1-E2)
   
  PS = PCAP*H2O_CRITICAL_PRESSURE
    
  if (calculate_derivatives) then
    dTC_dT = 1.d0/H2O_CRITICAL_TEMPERATURE
    dSC_dTC = -A(1)-2.d0*A(2)*one_m_tc-3.d0*A(3)*one_m_tc_sq- &
              4.d0*A(4)*one_m_tc**3.-5.d0*A(5)*one_m_tc**4.
    dE1_dTC = (1.d0+A(6)*one_m_tc+A(7)*one_m_tc_sq)+ &
              TC*(-A(6)-2.d0*A(7)*one_m_tc)
    dE2_dTC = -1.d0/E2_bottom+one_m_tc/(E2_bottom*E2_bottom)*2.d0*one_m_tc
    dPC_dTC = (-SC/(E1*E1)*dE1_dTC-dE2_dTC)*PCAP
    dPC_dSC = 1.d0/E1*PCAP
    dPS_dT = (dPC_dSC*dSC_dTC+dPC_dTC)*dTC_dT*H2O_CRITICAL_PRESSURE
  else
    dPS_dT = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterSaturationPressureIFC67

! ************************************************************************** !

subroutine EOSWaterSatPresWagnerPruss(T, calculate_derivatives, &
                                      PS, dPS_dT, ierr)
  ! 
  ! Calculates the saturation pressure of water as a function of temperature
  ! based on Wagner W. and A. Pruss (1993) "International Equations for the 
  ! Saturation Properties of Ordinary Water Substance. Revised According to 
  ! the International Temperature Scale of 1990. Addendum to J. Phys. Chem. 
  ! Ref. Data 16, 893 (1987)".
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/16
  ! 
  implicit none

  PetscReal, intent(in) :: T ! temperature
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: PS, dPS_dT ! Saturation pres. and derivative
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal, parameter :: a(6) = [-7.85951783d0,1.84408259d0,-11.7866497d0, &
                                  22.6807411d0,-15.9618719d0,1.80122502d0]
  PetscReal, parameter :: Tc = 647.096d0 ! [K]
  PetscReal, parameter :: Pc = 22.064d6  ! [Pa]
  PetscReal :: tau
  PetscReal :: T_over_Tc
  PetscReal :: polynomial
  PetscReal :: dpolynomial_dtau, dtau_dT
  
  T_over_Tc = (T+293.15d0)/Tc
  tau = 1.d0 - T_over_Tc
  polynomial = a(1)*tau+a(2)*tau**1.5d0+a(3)*tau**3.d0+ &
               a(4)*tau**3.5d0+a(5)*tau**4.d0+a(6)*tau**7.5d0
  PS = Pc*exp(polynomial/T_over_Tc)
  if (calculate_derivatives) then
    dtau_dT = -1.d0/Tc
    dpolynomial_dtau = a(1)+1.5d0*a(2)*sqrt(tau)+3.d0*a(3)*tau*tau+ &
                       3.5d0*a(4)*tau**2.5d0+4.d0*a(5)*tau**3.d0+ &
                       7.5d0*a(6)*tau**6.5d0
    dPS_dT = PS*(dpolynomial_dtau/T_over_Tc*dtau_dT+polynomial*Tc)
  else
    dPS_dT = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterSatPresWagnerPruss

! ************************************************************************** !

subroutine EOSWaterDensityIFC67(t,p,calculate_derivatives,dw,dwmol, &
                                dwp,dwt,ierr)

!  This subroutine calculates water and steam-gas mixture properties.
!  The water and steam properties are valid in the range of:
!
!            0 < p < 165.4 * 10^5 pascals (165.4 bars)
!            0 < t < 350 centigrade (623.15 Kelvin)
!
!  The properties cover densities, enthalpies, internal energies,
!  and partial derivatives of these quanties with respect to
!  pressure and temperature.
!
!  For saturated fluid, it will also calculate water saturation
!  temperature for the specified pressure using Newton-Raphson and
!  the derivative dts/dp (=tsp) or Ps for a given temperature.
!
!  Ref.: International Formulation Committee of the Sixth International
!       Conference on Properties of Steam (1967).

  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr
  
  PetscInt :: i
    
  PetscReal, save :: aa(0:22)
  PetscReal, save :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12
  
  PetscReal :: beta,beta2x,beta4,theta,utheta,theta2x,theta18,theta20
  PetscReal :: xx,yy,zz
  PetscReal :: u0,u1,u2,u3,u4,u5,u6,u7,u8,u9
  PetscReal :: tempreal
!  PetscReal :: v0_1, v1_1, v2_1, v3_1, v4_1
!  PetscReal :: v1_2, v2_2, v3_2, v4_2, v20_2, v40_2
!  PetscReal :: v1_3, v2_3, v3_3, v4_3
!  PetscReal :: v1_4, v2_4, v3_4
!  PetscReal :: v1_5, v2_5
!  PetscReal :: v1_6
!  PetscReal :: term1,term2,term2t,term3,term3t,term3p,term4,term4t,term4p, &
!               term5,term5t,term5p,term6,term6t,term6p,term7,term7t,term7p
!  PetscReal :: dv2t,dv2p,dv3t
  PetscReal :: vr,ypt,yptt,zpt,zpp,vrpt,vrpp,cnv
  PetscReal :: tc1,pc1,vc1,utc1,upc1,vc1mol
  PetscReal :: d2z_dp2    ! 2nd derivative of z w.r.t. pressure
  PetscReal :: d2vr_dp2   ! 2nd derivative of vr w.r.t. pressure
  PetscReal :: d2wmol_dp2 ! 2nd derivative of density w.r.t. pressure
  PetscReal, parameter :: zero = 0.d0
  PetscReal, parameter :: one = 1.d0
  PetscReal, parameter :: two = 2.d0
  PetscReal, parameter :: three = 3.d0
  PetscReal, parameter :: four = 4.d0
  PetscReal, parameter :: five = 5.d0
  PetscReal, parameter :: six = 6.d0
  PetscReal, parameter :: seven = 7.d0
  PetscReal, parameter :: eight = 8.d0
  PetscReal, parameter :: nine = 9.d0
  PetscReal, parameter :: ten = 10.d0
  
  data aa/ &
!-----data aa0,aa1,aa2,aa3/
        6.824687741d03,-5.422063673d02,-2.096666205d04, 3.941286787d04, &
!-----data aa4,aa5,aa6,aa7/
        -6.733277739d04, 9.902381028d04,-1.093911774d05, 8.590841667d04, &
!-----data aa8,aa9,aa10,aa11/
        -4.511168742d04, 1.418138926d04,-2.017271113d03, 7.982692717d00, &
!-----data aa12,aa13,aa14,aa15/
        -2.616571843d-2, 1.522411790d-3, 2.284279054d-2, 2.421647003d02, &
!-----data aa16,aa17,aa18,aa19/
        1.269716088d-10,2.074838328d-7, 2.174020350d-8, 1.105710498d-9, &
!-----data aa20,aa21,aa22/    
        1.293441934d01, 1.308119072d-5, 6.047626338d-14/

  data a1,a2,a3,a4/ &
  8.438375405d-1, 5.362162162d-4, 1.720000000d00, 7.342278489d-2/
  data a5,a6,a7,a8/ &
  4.975858870d-2, 6.537154300d-1, 1.150000000d-6, 1.510800000d-5/
  data a9,a10,a11,a12/ &
  1.418800000d-1, 7.002753165d00, 2.995284926d-4, 2.040000000d-1/
    
  tc1 = H2O_CRITICAL_TEMPERATURE    ! K
  pc1 = H2O_CRITICAL_PRESSURE       ! Pa 
  vc1 = 0.00317d0  ! m^3/kg
  utc1 = one/tc1   ! 1/C
  upc1 = one/pc1   ! 1/Pa
  vc1mol = vc1*FMWH2O ! m^3/kmol
    
  theta = (t+273.15d0)*utc1
  theta2x = theta*theta
  theta18 = theta**18.d0
  theta20 = theta18*theta2x
    
  beta = p*upc1
  beta2x = beta*beta
  beta4  = beta2x*beta2x
    
  yy = one-a1*theta2x-a2*theta**(-6.d0)
  xx = a3*yy*yy-two*(a4*theta-a5*beta)
    
!   Note: xx may become negative near the critical point-pcl.
  if (xx.gt.zero) then
    xx = sqrt(xx)
  else
    write(*,*) 'Warning: negative term in density (eos_water.F90:&
      &EOSWaterDensityIFC67):'
    write(*,*) 't= ',t,' p= ',p,' xx= ',xx
    ierr = 1
    xx = 1.d-6               !set arbitrarily
  end if
  zz = yy + xx                                     
  u0 = -five/17.d0
  u1 = aa(11)*a5*zz**u0
  u2 = one/(a8+theta**11.d0)
  u3 = aa(17)+(two*aa(18)+three*aa(19)*beta)*beta
  u4 = one/(a7+theta18*theta)
  u5 = (a10+beta)**(-4.d0)
  u6 = a11-three*u5
  u7 = aa(20)*theta18*(a9+theta2x)
  u8 = aa(15)*(a6-theta)**9.d0
    
  vr = u1+aa(12)+theta*(aa(13)+aa(14)*theta)+u8*(a6-theta) &
        +aa(16)*u4-u2*u3-u6*u7+(three*aa(21)*(a12-theta) &
        +four*aa(22)*beta/theta20)*beta2x
    
  dwmol = one/(vr*vc1mol) ! kmol/m^3
  dw = one/(vr*vc1) ! kg/m^3

  ! ypt used for enthalpy even if derivative not calculated
  ypt = six*a2*theta**(-7.d0)-two*a1*theta
  
  !---calculate derivatives for water density
  if (calculate_derivatives) then
    zpt = ypt+(a3*yy*ypt-a4)/xx
    zpp = a5/xx
    u9 = u0*u1/zz
    vrpt = u9*zpt+aa(13)+two*aa(14)*theta-ten*u8 &
        -19.d0*aa(16)*u4*u4*theta18+11.d0*u2*u2*u3*theta**10.d0 &
        -aa(20)*u6*(18.d0*a9*theta18+20.d0*theta20)/theta &
        -(three*aa(21)+80.d0*aa(22)*beta/(theta20*theta))*beta2x
    
    vrpp = u9*zpp-u2*(two*aa(18)+six*aa(19)*beta)-12.d0*u7*u5/ &
        (a10+beta)+(six*aa(21)*(a12-theta)+12.d0*aa(22)*beta/ &
        theta20)*beta
    
    cnv = -one/(vc1mol*vr*vr)
    dwt = cnv*vrpt*utc1 ! kmol/m^3/C
    dwp = cnv*vrpp*upc1 ! kmol/m^3/Pa

    ! 2nd derivative w.r.t pressure
    d2z_dp2 = -a5*a5/xx/xx/xx

    d2vr_dp2 = u9*(u0-one)/zz*zpp*zpp + &
               u0*u1/zz*d2z_dp2 + &
               six*u2*aa(19) + &
               60.d0*u7*u5/(a10+beta)/(a10+beta) + &
               six*(a12 - theta) + &
               24.d0*aa(22)*beta/theta20

    d2wmol_dp2 = -cnv*upc1*upc1*(2.d0/vr*vrpp*vrpp -  d2vr_dp2) ! kmol/m^3/Pa^2

  else
    dwt = UNINITIALIZED_DOUBLE
    dwp = UNINITIALIZED_DOUBLE
    d2wmol_dp2 = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityIFC67

! ************************************************************************** !

subroutine EOSWaterEnthalpyIFC67(t,p,calculate_derivatives,hw, &
                                 hwp,hwt,ierr)

!  This subroutine calculates water and steam-gas mixture properties.
!  The water and steam properties are valid in the range of:
!
!            0 < p < 165.4 * 10^5 pascals (165.4 bars)
!            0 < t < 350 centigrade (623.15 Kelvin)
!
!  The properties cover densities, enthalpies, internal energies,
!  and partial derivatives of these quanties with respect to
!  pressure and temperature.
!
!  For saturated fluid, it will also calculate water saturation
!  temperature for the specified pressure using Newton-Raphson and
!  the derivative dts/dp (=tsp) or Ps for a given temperature.
!
!  Ref.: International Formulation Committee of the Sixth International
!       Conference on Properties of Steam (1967).

  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(out) :: ierr
  
  PetscInt :: i
    
  PetscReal, save :: aa(0:22)
  PetscReal, save :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12
  
  PetscReal :: beta,beta2x,beta4,theta,utheta,theta2x,theta18,theta20
  PetscReal :: xx,yy,zz
  PetscReal :: u0,u1!,u2,u3,u4,u5,u6,u7,u8,u9
  PetscReal :: tempreal
  PetscReal :: v0_1, v1_1, v2_1, v3_1, v4_1
  PetscReal :: v1_2, v2_2, v3_2, v4_2, v20_2, v40_2
  PetscReal :: v1_3, v2_3, v3_3, v4_3
  PetscReal :: v1_4, v2_4, v3_4
  PetscReal :: v1_5, v2_5
  PetscReal :: v1_6
  PetscReal :: term1,term2,term2t,term3,term3t,term3p,term4,term4t,term4p, &
               term5,term5t,term5p,term6,term6t,term6p,term7,term7t,term7p
  PetscReal :: dv2t,dv2p,dv3t
!  PetscReal :: vr,ypt,yptt,zpt,zpp,vrpt,vrpp,cnv
  PetscReal :: ypt,yptt,zpt,zpp
  PetscReal :: tc1,pc1,vc1,utc1,upc1,vc1mol
  PetscReal, parameter :: zero = 0.d0
  PetscReal, parameter :: one = 1.d0
  PetscReal, parameter :: two = 2.d0
  PetscReal, parameter :: three = 3.d0
  PetscReal, parameter :: four = 4.d0
  PetscReal, parameter :: five = 5.d0
  PetscReal, parameter :: six = 6.d0
  PetscReal, parameter :: seven = 7.d0
  PetscReal, parameter :: eight = 8.d0
  PetscReal, parameter :: nine = 9.d0
  PetscReal, parameter :: ten = 10.d0
  
  data aa/ &
!-----data aa0,aa1,aa2,aa3/
        6.824687741d03,-5.422063673d02,-2.096666205d04, 3.941286787d04, &
!-----data aa4,aa5,aa6,aa7/
        -6.733277739d04, 9.902381028d04,-1.093911774d05, 8.590841667d04, &
!-----data aa8,aa9,aa10,aa11/
        -4.511168742d04, 1.418138926d04,-2.017271113d03, 7.982692717d00, &
!-----data aa12,aa13,aa14,aa15/
        -2.616571843d-2, 1.522411790d-3, 2.284279054d-2, 2.421647003d02, &
!-----data aa16,aa17,aa18,aa19/
        1.269716088d-10,2.074838328d-7, 2.174020350d-8, 1.105710498d-9, &
!-----data aa20,aa21,aa22/    
        1.293441934d01, 1.308119072d-5, 6.047626338d-14/

  data a1,a2,a3,a4/ &
  8.438375405d-1, 5.362162162d-4, 1.720000000d00, 7.342278489d-2/
  data a5,a6,a7,a8/ &
  4.975858870d-2, 6.537154300d-1, 1.150000000d-6, 1.510800000d-5/
  data a9,a10,a11,a12/ &
  1.418800000d-1, 7.002753165d00, 2.995284926d-4, 2.040000000d-1/
    
  tc1 = H2O_CRITICAL_TEMPERATURE    ! K
  pc1 = H2O_CRITICAL_PRESSURE     ! Pa 
  vc1 = 0.00317d0  ! m^3/kg
  utc1 = one/tc1   ! 1/C
  upc1 = one/pc1   ! 1/Pa
  vc1mol = vc1*FMWH2O ! m^3/kmol
    
  theta = (t+273.15d0)*utc1
  theta2x = theta*theta
  theta18 = theta**18.d0
  theta20 = theta18*theta2x
    
  beta = p*upc1
  beta2x = beta*beta
  beta4  = beta2x*beta2x
    
  yy = one-a1*theta2x-a2*theta**(-6.d0)
  xx = a3*yy*yy-two*(a4*theta-a5*beta)
!   Note: xx may become negative near the critical point-pcl.
  if (xx.gt.zero) then
    xx = sqrt(xx)
  else
    write(*,*) 'Warning: negative term in density (eos_water.F90:&
      &EOSWaterEnthalpyIFC67):'
    write(*,*) 't= ',t,' p= ',p,' xx= ',xx
    ierr = 1
    xx = 1.d-6               !set arbitrarily
  end if
  zz = yy + xx                                     
  u0 = -five/17.d0
  u1 = aa(11)*a5*zz**u0
#if 0  
  u2 = one/(a8+theta**11)
  u3 = aa(17)+(two*aa(18)+three*aa(19)*beta)*beta
  u4 = one/(a7+theta18*theta)
  u5 = (a10+beta)**(-4)
  u6 = a11-three*u5
  u7 = aa(20)*theta18*(a9+theta2x)
  u8 = aa(15)*(a6-theta)**9
    
  vr = u1+aa(12)+theta*(aa(13)+aa(14)*theta)+u8*(a6-theta) &
        +aa(16)*u4-u2*u3-u6*u7+(three*aa(21)*(a12-theta) &
        +four*aa(22)*beta/theta20)*beta2x
    
  dwmol = one/(vr*vc1mol) ! kmol/m^3
  dw = one/(vr*vc1) ! kg/m^3

  !---calculate derivatives for water density
  if (calculate_derivatives) then
    u9 = u0*u1/zz
    vrpt = u9*zpt+aa(13)+two*aa(14)*theta-ten*u8 &
        -19.d0*aa(16)*u4*u4*theta18+11.d0*u2*u2*u3*theta**10.d0 &
        -aa(20)*u6*(18.d0*a9*theta18+20.d0*theta20)/theta &
        -(three*aa(21)+80.d0*aa(22)*beta/(theta20*theta))*beta2x
    
    vrpp = u9*zpp-u2*(two*aa(18)+six*aa(19)*beta)-12.d0*u7*u5/ &
        (a10+beta)+(six*aa(21)*(a12-theta)+12.d0*aa(22)*beta/ &
        theta20)*beta
    
    cnv = -one/(vc1mol*vr*vr)
    dwt = cnv*vrpt*utc1 ! kmol/m^3/C
    dwp = cnv*vrpp*upc1 ! kmol/m^3/Pa
  else
    dwt = UNINITIALIZED_DOUBLE
    dwp = UNINITIALIZED_DOUBLE
  endif
#endif  

  ! ypt used for enthalpy even if derivative not calculated
  ypt = six*a2*theta**(-7.d0)-two*a1*theta
  

!---compute enthalpy internal energy and derivatives for water
  utheta = one/theta
  term1 = aa(0)*theta
  term2 = -aa(1)
  ! term2t is part of the derivative calc., but left here to avoid
  ! recomputing the expensive do loop
  term2t = zero
  do i = 3,10
    tempreal = dfloat(i-2)*aa(i)*theta**(i-1)
    term2t = term2t+tempreal*utheta*dfloat(i-1)
    term2 = term2+tempreal                            
  end do
    
  ! "v" section 1
  v0_1 = u1/a5
  v2_1 = 17.d0*(zz/29.d0-yy/12.d0)+five*theta*ypt/12.d0
  v3_1 = a4*theta-(a3-one)*theta*yy*ypt
  v1_1 = zz*v2_1+v3_1
  term3 = v0_1*v1_1
  
  ! block 1 removed from here

  ! "v" section 2
  v1_2 = nine*theta+a6
  v20_2 = (a6-theta)
  v2_2 = v20_2**9.d0
  v3_2 = a7+20.d0*theta**19.d0
  v40_2 = a7+theta**19.d0
  v4_2 = one/(v40_2*v40_2)
  ! term4p is a derivative, but left due to dependency in term4
  term4p = aa(12)-aa(14)*theta2x+aa(15)*v1_2*v2_2+aa(16)*v3_2*v4_2
  term4 = term4p*beta
  
  ! block 2 removed from here
    
  ! "v" section 3
  v1_3 = beta*(aa(17)+aa(18)*beta+aa(19)*beta2x)
  v2_3 = 12.d0*theta**11.d0+a8
  v4_3 = one/(a8+theta**11.d0)
  v3_3 = v4_3*v4_3
  term5 = v1_3*v2_3*v3_3
  
  ! block 3 removed from here
    
  ! "v" section 4
  v1_4 = (a10+beta)**(-3.d0)+a11*beta
  v3_4 = (17.d0*a9+19.d0*theta2x)
  v2_4 = aa(20)*theta18*v3_4       
  term6 = v1_4*v2_4

  ! block 4 removed from here
    
  ! "v" section 5
  v1_5 = 21.d0*aa(22)/theta20*beta4
  v2_5 = aa(21)*a12*beta2x*beta
  term7 = v1_5+v2_5  
    
  ! "v" section 6
  v1_6 = pc1*vc1mol
  hw = (term1-term2+term3+term4-term5+term6+term7)*v1_6
    
  if (calculate_derivatives) then

    zpt = ypt+(a3*yy*ypt-a4)/xx
    zpp = a5/xx
  
    ! block 1
    yptt = -two*a1-42.d0*a2/theta**8.d0
    dv2t = 17.d0*(zpt/29.d0-ypt/12.d0)+five/12.d0*(ypt+theta*yptt) 
    dv3t = a4-(a3-one)*(theta*yy*yptt+yy*ypt+theta*ypt*ypt)
    dv2p = 17.d0*zpp/29.d0
    v4_1 = five*v1_1/(17.d0*zz)       
    term3t = v0_1*(zz*dv2t+(v2_1-v4_1)*zpt+dv3t)
    term3p = v0_1*(zz*dv2p+(v2_1-v4_1)*zpp)
  
    ! block 2
    term4t = (-two*aa(14)*theta+nine*aa(15)*(v2_2-v1_2*v2_2/v20_2) &
             +38.d0*theta18*aa(16)*(ten*v4_2-v3_2*v4_2/v40_2))*beta

    ! block 3
    term5p = v3_3*v2_3*(aa(17)+two*aa(18)*beta+three*aa(19)*beta2x)
    term5t = v1_3*(132.d0*v3_3*theta**10.d0-22.d0*v2_3*v3_3*v4_3*theta**10.d0)
  
    ! block 4
    term6p = v2_4*(a11-three*(a10+beta)**(-4.d0))
    term6t = v1_4*aa(20)*theta18*(18.d0*v3_4*utheta+38.d0*theta)
    
    ! block 5
    term7p = beta2x*(three*aa(21)*a12+84.d0*aa(22)*beta/theta20)
    term7t = -420.d0*aa(22)*beta4/(theta20*theta)

    hwp = (term3p+term4p-term5p+term6p+term7p)*vc1mol
    hwt = (aa(0)-term2t+term3t+term4t-term5t+term6t+term7t)*v1_6*utc1
  else
    hwp = UNINITIALIZED_DOUBLE
    hwt = UNINITIALIZED_DOUBLE
  endif
    
end subroutine EOSWaterEnthalpyIFC67

! ************************************************************************** !

subroutine EOSWaterDensityConstant(t,p,calculate_derivatives,dw,dwmol, &
                                   dwp,dwt,ierr)
  implicit none
  
  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr
  
  dw = constant_density ! kg/m^3
  dwmol = dw/FMWH2O ! kmol/m^3
  
  dwp = 0.d0
  dwt = 0.d0

end subroutine EOSWaterDensityConstant

! ************************************************************************** !

subroutine EOSWaterEnthalpyConstant(t,p,calculate_derivatives, &
                                    hw,hwp,hwt,ierr)
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(out) :: ierr
  
  hw = constant_enthalpy ! J/kmol
  
  hwp = 0.d0
  hwt = 0.d0
  
end subroutine EOSWaterEnthalpyConstant

! ************************************************************************** !

subroutine EOSWaterDensityExponential(t,p,calculate_derivatives, &
                                      dw,dwmol,dwp,dwt,ierr)
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr
  
  ! kg/m^3
  dw = exponent_reference_density*exp(exponent_water_compressibility* &
                                     (p-exponent_reference_pressure))
  dwmol = dw/FMWH2O ! kmol/m^3
  
  if (calculate_derivatives) then
    dwp = dwmol*exponent_water_compressibility !kmol/m^3/Pa
  else
    dwp = UNINITIALIZED_DOUBLE
  endif
  dwt = 0.d0

end subroutine EOSWaterDensityExponential

! ************************************************************************** !

subroutine EOSWaterDensityLinear(t,p,calculate_derivatives, &
                                      dw,dwmol,dwp,dwt,ierr)
  !
  ! Water density linear model
  !
  ! Author: Satish Karra 
  ! Date: 06/19/17

  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr

  ! kg/m^3
  dw = linear_reference_density*(1.d0 + & 
         linear_water_compressibility*(p-linear_reference_pressure))
  
  dwmol = dw/FMWH2O ! kmol/m^3
  
  if (calculate_derivatives) then
    dwp = linear_reference_density*linear_water_compressibility/FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
  endif
  dwt = 0.d0

end subroutine EOSWaterDensityLinear

! ************************************************************************** !

subroutine EOSWaterDensityBRAGFLO(t,p,calculate_derivatives, &
                                  dw,dwmol,dwp,dwt,ierr)
  !
  ! Water density based on formulation in BRAGFLO.  The BRAGFLO user manual
  ! is incorrect as it does not include the truncation in the code (see
  ! subroutine DENO near line 6848 of Bragflo.f
  !
  ! Author: Glenn Hammond
  ! Date: 06/02/17
                                  
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade (ignored)
  PetscReal, intent(in) :: p   ! Pressure in Pascal
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: p_adjust
  
  ! kg/m^3
  p_adjust = max(p,0.d0)
  dw = exponent_reference_density* &
         exp(min(80.d0,exponent_water_compressibility* &
                       (p_adjust-exponent_reference_pressure)))
  
  dwmol = dw/FMWH2O ! kmol/m^3
  
  if (calculate_derivatives) then
    print *, 'Analytical derivatives in EOSWaterDensityBRAGFLO() &
      &currently not possible due to discontinuities in the formulation.'
    stop
    !dwp = dwmol*exponent_water_compressibility !kmol/m^3/Pa
  else
    dwp = UNINITIALIZED_DOUBLE
  endif
  dwt = 0.d0

end subroutine EOSWaterDensityBRAGFLO

! ************************************************************************** !

subroutine EOSWaterDensityQuadratic(t,p,calculate_derivatives, &
                                      dw,dwmol,dwp,dwt,ierr)
  !
  ! Water density quadratic model
  !
  ! Author: Paolo Orsini
  ! Date: 02/10/17

  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: X_pr

  ! [-] =                         1/Pa *  Pa                        
  X_pr = quadratic_wat_compressibility * (p - quadratic_reference_pressure )

  ! kg/m^3
  dw = quadratic_reference_density * (1.0 + X_pr + (X_pr**2)/2.0d0 )

  dwmol = dw/FMWH2O ! kmol/m^3
  
  if (calculate_derivatives) then
          ! kg/m^3 * 1/Pa * [-] / (kg/kmol) = kmol/m^3/Pa
    dwp = quadratic_reference_density * quadratic_wat_compressibility * &
          (1.0 + X_pr) / FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
  endif
  dwt = 0.d0

end subroutine EOSWaterDensityQuadratic

! ************************************************************************** !

subroutine EOSWaterDensityTrangenstein(t,p,calculate_derivatives, &
                                       dw,dwmol,dwp,dwt,ierr)
  !
  ! Water density model taken from Trangenstein
  ! Trangenstein, J. A. “Analysis of a Model and Sequential Numerical Method
  ! for Thermal Reservoir Simulation,” The Mathematics of Oil Recovery, 
  ! King, P.R. (ed.), Oxford, UK, Clarendon Press (1992), 359.
  !
  ! Author: Paolo Orsini
  ! Date: 02/22/17

  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr

  PetscReal, parameter :: a0 = 9.99839520d+02
  PetscReal, parameter :: a1 = 1.69551760d+01
  PetscReal, parameter :: a2 = -7.98700000d-03
  PetscReal, parameter :: a3 = -4.61704610d-05
  PetscReal, parameter :: a4 = 1.05563020d-07
  PetscReal, parameter :: a5 = -2.80543530d-10
  PetscReal, parameter :: a6 = 1.68798500d-02
  PetscReal, parameter :: a7 = 10.20

  PetscReal, parameter :: cpw = 4.00d-005  !1/atm
  ! conversion 1/atm -> 1/Pa         
  ! cpw_Pa = cpw * 1.01325*1.0d-1

  PetscReal :: t_numer, d_t_numer_dt, t_denom, d_t_denom_dt,  p_exponent, d_p_exponent_dp
  
  dw = (a0 + a1*t + a2 * t**2.0 + a3 * t**3.0 + a4 * t**4.0 + a5 * t**5.0 ) / &
       ( 1.0 + a6 * t ) * &           !Pa -> MPa
       dexp( cpw / (1.01325*1.0d-1) * (p * 1.0d-6 - a7 ) )
               !1/atm -> 1/MPa    

  dwmol = dw/FMWH2O ! kmol/m^3
  
  if (calculate_derivatives) then

    !!! DS - FIRST ATTEMPT at these derivatives
    t_numer = (a0 + a1*t + a2 * t**2.0 + a3 * t**3.0 + a4 * t**4.0 + a5 * t**5.0 ) 
    d_t_numer_dt = (a1 + 2.0* a2 * t +  3.0 * a3 * t**2.0 + 4.0 * a4 * t**3.0 + 5.0 * a5 * t**4.0 ) 

    t_denom = ( 1.0 + a6 * t )
    d_t_denom_dt = a6 

    p_exponent = cpw / (1.01325*1.0d-1) * (p * 1.0d-6 - a7 )
    d_p_exponent_dp = cpw / (1.01325*1.0d-1) * 1.0d-6 
    
    dwt = (d_t_numer_dt / t_denom) - (t_numer * d_t_denom_dt / t_denom / t_denom  )
    dwt = dwt * dexp(p_exponent)

    dwp = dw * d_p_exponent_dp 

    !! derivatives are undestood to be in mol units so
    dwt = dwt/FMWH2O
    dwp = dwp/FMWH2O

    !PO TODO add derivatives
    !print *, 'Derivatives not set up in EOSWaterDensityTrangenstein.'
    !stop
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityTrangenstein

! ************************************************************************** !

subroutine EOSWaterViscosityGrabowski(T, P, PS, dPS_dT, &
                                      calculate_derivatives, VW, &
                                      dVW_dT, dVW_dP, ierr)
  ! 
  ! Grabowski, J. W. and Rubin, B. “A Preliminary Numerical Simulation 
  ! Study of In-situ Combustion in a Cold Lake Oil Sands Reservoir” 
  ! Journal of Canadian Petroleum Technology, (1981) 20, No. 2, 79-89.
  !
  ! Author: Paolo Orsini
  ! Date: 02/26/17
  ! 
  implicit none
  PetscReal, intent(in) :: T       ! C
  PetscReal, intent(in) :: P       ! Pa
  PetscReal, intent(in) :: PS      ! Pa
  PetscReal, intent(in) :: dPS_dT  ! Pa/C
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW     ! Pa-s
  PetscReal, intent(out) :: dVW_dT, dVW_dP
  PetscErrorCode, intent(out) :: ierr
  
  ! convert from centipoise to Pa-s (1 cP = 1.d-3 Pa-s)
  PetscReal, parameter :: centipoise_to_Pa_s = 1.d-3
  PetscReal, parameter :: a = 2.1850
  PetscReal, parameter :: b = 0.04012
  PetscReal, parameter :: c = 5.1547d-6

  PetscReal :: t_F ! temperature in F
  PetscReal :: denom, d_denom_d_t_F, d_VW_d_t_F, d_t_F_d_T
  
  t_F = T*(9.0/5.0)+32.0

  VW = a / (-1.0  + b * t_F + c * t_F**2.0 )
  VW = VW * centipoise_to_Pa_s
       
  if (calculate_derivatives) then

    !!! DS - FIRST ATTEMPT at these derivatives
    denom = (-1.0  + b * t_F + c * t_F**2.0 )
    d_denom_d_t_F = ( b + 2.0 * c * t_F )
    !! simple chain rule for this:
    d_VW_d_t_F = -1.0 * a * d_denom_d_t_F / denom / denom

    d_t_F_d_T = 9.0/5.0
    !! and again simple chain rule:
    dVW_dT = d_t_F_d_T * d_VW_d_t_F
    !! scaling needed:
    dVW_dT = dVW_dT * centipoise_to_Pa_s


    
    dVW_dP = 0.d0

    !PO TODO add derivatives
    !print *, 'Derivatives not set up in EOSWaterViscosityGrabowski'
    !stop
    !dVW_dT =
  else
    dVW_dP = UNINITIALIZED_DOUBLE
    dVW_dT = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterViscosityGrabowski

! ************************************************************************** !

subroutine EOSWaterDensityTPPlanarSetup()
  ! 
  ! Setups up plane for density of water as a function of temperature and 
  ! pressure using a simple plane equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/31/17
  ! 
  implicit none
  
  PetscReal, parameter :: p0 = 1.d6
  PetscReal, parameter :: t0 = 30.d0
  PetscReal, parameter :: drho_dp = 4.42d-7
  PetscReal, parameter :: drho_dT = -0.312d0
  PetscReal, parameter :: rho_reference = 996.d0
  
  PetscReal :: dw
  
  call GeomComputePlaneWithGradients(water_density_tp_plane,p0,t0, &
                                     rho_reference,drho_dp,drho_dT)

end subroutine EOSWaterDensityTPPlanarSetup

! ************************************************************************** !

subroutine EOSWaterDensityTPPlanar(t,p,calculate_derivatives, &
                                   dw,dwmol,dwp,dwt,ierr)
  ! 
  ! Calculates the density of water as a function of temperature and pressure
  ! using a simple plane equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/04/16
  ! 
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr
  
  ! kg/m^3
  dw = GeometryGetPlaneZIntercept(water_density_tp_plane,p,t)
  dwmol = dw/FMWH2O ! kmol/m^3
  
  if (calculate_derivatives) then
    call GeomGetPlaneGradientinXandY(water_density_tp_plane,dwp,dwt)
    dwp = dwp/FMWH2O
    dwt = dwt/FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityTPPlanar

! ************************************************************************** !

subroutine EOSWaterEnthalpyTPPlanarSetup()
  ! 
  ! Setups up plane for enthalpy of water as a function of temperature and 
  ! pressure using a simple plane equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/31/17
  ! 
  implicit none
  
  PetscReal, parameter :: p0 = 1.d6
  PetscReal, parameter :: t0 = 30.d0
  PetscReal, parameter :: dh_dp = 1.65d-2
  PetscReal, parameter :: dh_dT = 7.5d4 
  PetscReal, parameter :: h_reference = 2.27d6 ! J/kmol
  
  call GeomComputePlaneWithGradients(water_enthalpy_tp_plane,p0,t0, &
                                     h_reference,dh_dp,dh_dT)  

end subroutine EOSWaterEnthalpyTPPlanarSetup

! ************************************************************************** !

subroutine EOSWaterEnthalpyTPPlanar(t,p,calculate_derivatives, &
                                    hw,hwp,hwt,ierr)
  ! 
  ! Calculates the enthalpy of water as a function of temperature and pressure
  ! using a simple plane equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/04/16
  ! 
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: hw,hwp,hwt
  PetscErrorCode, intent(out) :: ierr
  
  ! kg/m^3
  hw = GeometryGetPlaneZIntercept(water_enthalpy_tp_plane,p,t)
  
  if (calculate_derivatives) then
    call GeomGetPlaneGradientinXandY(water_enthalpy_tp_plane,hwp,hwt)
  else
    hwp = UNINITIALIZED_DOUBLE
    hwt = UNINITIALIZED_DOUBLE
  endif
  
end subroutine EOSWaterEnthalpyTPPlanar

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthNoDerive(t,p,dg,dgmol,hg,ierr)

  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: p   ! Vapor Pressure in Pascals.
  PetscReal, intent(out) :: dg,dgmol
  PetscReal, intent(out) :: hg
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal :: dum1, dum2, dum3, dum4
  
  call EOSWaterSteamDensityEnthalpyPtr(t,p,PETSC_FALSE,dg,dgmol,hg, &
                                       dum1,dum2,dum3,dum4,ierr)
  
end subroutine EOSWaterSteamDenEnthNoDerive

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthDerive(t,pv,dg,dgmol,hg, &
                                      dgp,dgt,hgp,hgt,ierr)
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: pv  ! Vapor Pressure in Pascals.
  PetscReal, intent(out) :: dg,dgmol,dgp,dgt
  PetscReal, intent(out) :: hg,hgp,hgt
  PetscErrorCode, intent(out) :: ierr
  
  call EOSWaterSteamDensityEnthalpyPtr(t,pv,PETSC_TRUE,dg,dgmol,hg,dgp, &
                                       dgt,hgp,hgt,ierr)
  
end subroutine EOSWaterSteamDenEnthDerive

! ************************************************************************** !

subroutine EOSWaterSteamDensityEnthalpyIFC67(t,pv,calculate_derivatives, &
                                             dg,dgmol,hg, &
                                             dgp,dgt,hgp,hgt,ierr)
! t/C  p/Pa dgmol/(mol/m^3)  h/J/kmol
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: pv  ! Vapor Pressure in Pascals.
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dg,dgmol,dgp,dgt
  PetscReal, intent(out) :: hg,hgp,hgt
  PetscErrorCode, intent(out) :: ierr
  
  PetscInt, save :: n(8),ll(8),x(8,2),z(8,3)
  PetscReal, save :: xi1,xl0,xl1,xl2
  PetscReal, save :: b(8,2),bb(0:9,0:6)
  PetscReal :: sumbx(8),sumbxt(8)
  
  PetscInt :: i,j
  PetscReal, save :: delp,delt
  PetscReal :: beta,betap,betal,betalp,betalp1,betal1,ubeta,theta, &
            thetap,utheta
  PetscReal :: xx,xxp
  PetscReal :: f1,f2,fim1,fim2,sum,sumt,sum1,sum1t,sum1p, &
            sum2,sum2t
  PetscReal :: u1,u1t,u1p,u2,u2p,u2t,u3,u3p,u3t,v1,v1t
  PetscReal :: term,term1,term1t,term1p,term2,term2t,term2p, &
            term3,term3t,term3p,term4,term4t,term4p,term5,term5t,term5p
  PetscReal :: hr,hrp,hrt,hrpt,hrpp
  PetscReal :: vr,vrpt,vrpp
  PetscReal :: tc1,pc1,vc1,utc1,upc1,vc1mol
  PetscReal, parameter :: zero = 0.d0
  PetscReal, parameter :: one = 1.d0
  PetscReal, parameter :: two = 2.d0
  PetscReal, parameter :: three = 3.d0
  PetscReal, parameter :: four = 4.d0
  PetscReal, parameter :: five = 5.d0
  PetscReal, parameter :: six = 6.d0
  PetscReal, parameter :: seven = 7.d0
  PetscReal, parameter :: eight = 8.d0
  PetscReal, parameter :: nine = 9.d0
  PetscReal, parameter :: ten = 10.d0
 
  data delt,delp/1.d-6,1.d-6/
    
  data n/2,3,2,2,3,2,2,2/
  data ll/0,0,0,0,0,1,1,2/
  data x/0,0,0,0,0,14,19,54, &
          0,0,0,0,0, 0, 0,27/
  data z/13,18,18,25,32,12,24,24, &
          3, 2,10,14,28,11,18,14, &
          0, 1, 0, 0,24, 0, 0, 0/
    
  data b/7.6333333333d-1,0.d0,0.d0,0.d0,0.d0,4.006073948d-1, &
          8.636081627d-2,-8.532322921d-1, &
          0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,3.460208861d-1/

!---sub-region 2
  data bb/ &
!---bb(0-9,0)
          1.683599274d1,  8*0.d0        , 1.936587558d2, &
!---bb(0-9,1)
          2.856067796d01, 6.670375918d-2, 8.390104328d-2, &
          4.520918904d-1,-5.975336707d-1, 5.958051609d-1, &
          1.190610271d-1, 1.683998803d-1, 6.552390126d-3, &
          -1.388522425d+3, &
!---bb(0-9,2)
          -5.438923329d01, 1.388983801d0 , 2.614670893d-2, &
          1.069036614d-1,-8.847535804d-2,-5.159303373d-1, &
          -9.867174132d-2,-5.809438001d-2, 5.710218649d-4, &
          4.126607219d03, &
!---bb(0-9,3)
          4.330662834d-1, 0.d0          ,-3.373439453d-2, &
          0.d0          , 0.d0          , 2.075021122d-1, &
          0.d0          , 0.d0          , 0.d0          , &
          -6.508211677d03, &
!---bb(0-9,4)
          -6.547711697d-1, 8*0.d0        , 5.745984054d03, &
!---bb(0-9,5)
          8.565182058d-2, 8*0.d0        ,-2.693088365d03, &
!---bb(0-9,6)
                          9*0.d0        , 5.235718623d02/
    
  data xi1/4.260321148d0/
  data xl0,xl1,xl2/15.74373327d0,-34.17061978d0,19.31380707d0/
    
  tc1 = H2O_CRITICAL_TEMPERATURE       ! K
  pc1 = H2O_CRITICAL_PRESSURE        ! Pa
  vc1 = 0.00317d0     ! m^3/kg
  utc1 = one/tc1
  upc1 = one/pc1
  vc1mol = vc1*FMWH2O ! m^3/kmol

  theta  = (t+273.15d0)*utc1
  beta   = pv*upc1
  ubeta  = one/beta
  utheta = one/theta
  xx = exp(b(1,1)*(one-theta))
    
!---compute steam density and derivatives
    
  term1 = xi1*theta*ubeta
  term1t = xi1*ubeta
  term1p = -term1*ubeta        
    
  do i = 1,8
    sum = zero
    sumt = zero
    do j = 1,n(i)
      u1 = bb(i,j)*xx**z(i,j)     
      sumt = sumt-b(1,1)*z(i,j)*u1                     
      sum = sum + u1                     
    end do
    sumbx(i) = sum
    sumbxt(i) = sumt
  end do
    
  term2  = sumbx(1)+beta*(two*sumbx(2)+beta*(three*sumbx(3) &
            +beta*(four*sumbx(4)+five*beta*sumbx(5))))
  term2t = sumbxt(1)+beta*(two*sumbxt(2)+beta*(three*sumbxt(3) &
            +beta*(four*sumbxt(4)+five*beta*sumbxt(5))))
  term2p = two*sumbx(2)+beta*(six*sumbx(3)+beta*(12.d0*sumbx(4) &
            +20.d0*beta*sumbx(5)))
    
  term3  = zero
  term3p = zero
  term3t = zero
    
  do i = 6,8
    fim1 = dfloat(i-1)
    fim2 = fim1-one    
    
    sum2 = zero
    sum2t = zero
    do j = 1,ll(i)
      u1 = b(i,j)*xx**x(i,j)
      sum2t = sum2t-b(1,1)*x(i,j)*u1
      sum2 = sum2 + u1
    end do
    
    f1 = fim2*beta**(1-i)
    f2 = beta**(2-i)
    u1 = one/(f2+sum2)**2
    term = u1*f1*sumbx(i)
    term3 = term3+term      
    u2 = two*term/(f2+sum2)
    term3t = term3t+u1*(f1*sumbxt(i)+u2*sum2t)    
    term3p = term3p-u1*(f1*fim1*sumbx(i)+u2*fim2*f2)*ubeta
  end do
    
  term4 = bb(9,0)
  term4t = zero  
  do i = 1,6
    u1 = bb(9,i)*xx**i 
    term4t = term4t-b(1,1)*dfloat(i)*u1
    term4 = term4 + u1             
  end do
    
  betal = xl0+theta*(xl1+xl2*theta)
  betalp = xl1+ two*xl2*theta
  u1 = 11.d0*(beta/betal)**10
  term4t = u1*(term4t-10.d0*betalp*term4/betal)
  term4 = term4*u1                  
  term4p = 10.d0*term4*ubeta             
    
  vr = term1-term2-term3+term4
  vrpt = term1t-term2t-term3t+term4t
  vrpp = term1p-term2p-term3p+term4p
    
  u1 = -one/(vc1mol*vr)
  dgmol = -u1
  dg = one/(vc1*vr)
  
  if (calculate_derivatives) then
    u1 = u1/vr
    dgt = u1*vrpt*utc1
    dgp = u1*vrpp*upc1
  else
    dgt = UNINITIALIZED_DOUBLE
    dgp = UNINITIALIZED_DOUBLE
  endif
    
!---compute steam enthalpy, internal energy, and derivatives
  thetap = theta+delt 
  betap  = beta +delp 
  xxp    = exp(b(1,1)*(one-thetap))
    
  term1  = bb(0,0)*theta
  term1t = bb(0,0)*thetap
    
  term2  = zero
  term2t = zero
  do i = 1,5
    u1 = bb(0,i)*(dfloat(i)-two)
    term2t = term2t+u1*thetap**(i-1)
    term2  = term2+u1*theta**(i-1)
  end do
    
  term3  = zero
  term3t = zero
  term3p = zero
  u1 = one
  u1p = one
  do i = 1,5
    u1 = u1*beta
    u1p = u1p*betap
    
    sum = zero
    sumt = zero
    do j = 1,n(i)
      sumt = sumt+bb(i,j)*(one+z(i,j)*b(1,1)*thetap)*xxp**z(i,j)
      sum  = sum+bb(i,j)*(one+z(i,j)*b(1,1)*theta)*xx**z(i,j)
    end do
    
    term3t = term3t+u1*sumt
    term3p = term3p+u1p*sum
    term3  = term3+u1*sum
  end do
    
  term4  = zero
  term4t = zero
  term4p = zero
  do i = 6,8
    
    sum1  = zero
    sum2  = zero
    sum1t = zero
    sum2t = zero
    
    do j = 1,ll(i)
      u1 = b(i,j)*xxp**x(i,j)
      sum1t = sum1t+x(i,j)*u1               
      sum2t = sum2t+u1                 
      u1 = b(i,j)*xx**x(i,j)
      sum1 = sum1+x(i,j)*u1               
      sum2 = sum2+u1                 
    end do
    
    u1 = one/(beta**(2-i)+sum2)
    u2 = one-b(1,1)*theta*sum1*u1                  
    u3 = b(1,1)*theta
    
    u1t = one/(beta**(2-i)+sum2t)
    u2t = one-b(1,1)*thetap*sum1t*u1t                 
    u3t = b(1,1)*thetap
    
    u1p = one/(betap**(2-i)+sum2)
    u2p = one-b(1,1)*theta*sum1*u1p                   
    u3p = u3          
    
    sum1 = zero
    sum1t = zero
    sum1p = zero
    do j = 1,n(i)
      sum1t = sum1t + bb(i,j)*xxp**z(i,j)*(u2t+z(i,j)*u3t)
      sum1p = sum1p + bb(i,j)*xx **z(i,j)*(u2p+z(i,j)*u3p)
      sum1  = sum1  + bb(i,j)*xx **z(i,j)*(u2 +z(i,j)*u3)
    end do
    
    term4t = term4t+sum1t*u1t
    term4p = term4p+sum1p*u1p
    term4  = term4+sum1*u1
  end do
    
  u1 = ten*betalp/betal
  term5 = (one+theta*u1)*bb(9,0)
    
  betal1  = xl0+thetap*(xl1+xl2*thetap)
  betalp1 = xl1+two*xl2*thetap
  u1t     = ten*betalp1/betal1
  term5t  = (one+thetap*u1t)*bb(9,0)
    
  do i = 1,6
    v1     = one+theta*(u1+dfloat(i)*b(1,1))
    v1t    = one+thetap*(u1t+dfloat(i)*b(1,1))
    term5t = v1t*bb(9,i)*xxp**i+term5t                        
    term5  =  v1*bb(9,i)*xx **i+term5                         
  end do
    
  term5  = term5*beta*(beta/betal)**10
  term5t = term5t*beta*(beta/betal1)**10
  term5p = term5*(betap/beta)**11        
    
  hr   = term1 -term2 -term3 -term4 +term5
  hrt  = term1t-term2t-term3t-term4t+term5t
  hrp  = term1 -term2 -term3p-term4p+term5p
  hrpt = (hrt-hr)/delt
  hrpp = (hrp-hr)/delp
    
  v1 = pc1*vc1mol  ! Pa = (nRT/V) = J/m^3 : J/m^3 * m^3/kmol = J/kmol
  hg = hr*v1       ! J/kmol
  
  if (calculate_derivatives) then
    hgt = hrpt*v1*utc1
    hgp = hrpp*vc1mol
  else
    hgt = UNINITIALIZED_DOUBLE
    hgp = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterSteamDensityEnthalpyIFC67

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthConstant(t,pv,calculate_derivatives, &
                                        dg,dgmol,hg,dgp,dgt,hgp,hgt,ierr)
! t/C  p/Pa dgmol/(mol/m^3)  h/J/kmol
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: pv  ! Vapor Pressure in Pascals.
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dg,dgmol,dgp,dgt
  PetscReal, intent(out) :: hg,hgp,hgt
  PetscErrorCode, intent(out) :: ierr
  
  dg = constant_steam_density
  dgmol = dg/FMWH2O
  hg = constant_steam_enthalpy
  
  if (calculate_derivatives) then
    dgt = 0.d0
    dgp = 0.d0
    hgt = 0.d0
    hgp = 0.d0
  else
    dgt = UNINITIALIZED_DOUBLE
    dgp = UNINITIALIZED_DOUBLE
    hgt = UNINITIALIZED_DOUBLE
    hgp = UNINITIALIZED_DOUBLE
  endif
  
end subroutine EOSWaterSteamDenEnthConstant

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthTPPlanarSetup()
  ! 
  ! Setups up plane for density adn enthalpy of steam as a function of 
  ! temperature and pressure using a simple plane equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/31/17
  ! 
  implicit none
  
  PetscReal, parameter :: p0 = 1.d5
  PetscReal, parameter :: t0 = 30.d0
  PetscReal, parameter :: drho_dp = 7.8d-6
  PetscReal, parameter :: drho_dT = -3.54d-3
  PetscReal, parameter :: rho_reference = 0.846d0 ! kg/m^3
  PetscReal, parameter :: dh_dp = -5.2d0
  PetscReal, parameter :: dh_dT = 4.d4 
  PetscReal, parameter :: h_reference = 4.54d7 ! J/kmol
  
  call GeomComputePlaneWithGradients(steam_density_tp_plane,p0,t0, &
                                     rho_reference,drho_dp,drho_dT)  
  call GeomComputePlaneWithGradients(steam_enthalpy_tp_plane,p0,t0, &
                                     h_reference,dh_dp,dh_dT)  

end subroutine EOSWaterSteamDenEnthTPPlanarSetup

! ************************************************************************** !

subroutine EOSWaterSteamDenEnthTPPlanar(t,pv,calculate_derivatives, &
                                        dv,dvmol,hv,dvp,dvt,hvp,hvt,ierr)
! t/C  p/Pa dgmol/(mol/m^3)  h/J/kmol
                                          ! 
  ! Calculates the enthalpy of water as a function of temperature and pressure
  ! using a simple plane equation
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/04/16
  ! 
  implicit none
  
  PetscReal, intent(in) :: t   ! Temperature in centigrade.
  PetscReal, intent(in) :: pv  ! Vapor Pressure in Pascals.
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dv,dvmol,dvp,dvt ! kmol/m^3
  PetscReal, intent(out) :: hv,hvp,hvt       ! J/kmol
  PetscErrorCode, intent(out) :: ierr
  
  ! FMWAIR = 28.96d0

  ! kg/m^3
  dv = GeometryGetPlaneZIntercept(steam_density_tp_plane,pv,t)
  dvmol = dv/FMWH2O ! kmol/m^3
  hv = GeometryGetPlaneZIntercept(steam_enthalpy_tp_plane,pv,t)
  
  if (calculate_derivatives) then
    call GeomGetPlaneGradientinXandY(steam_density_tp_plane,dvp,dvt)
    dvp = dvp/FMWH2O
    dvt = dvt/FMWH2O
    call GeomGetPlaneGradientinXandY(steam_enthalpy_tp_plane,hvp,hvt)
  else
    dvp = UNINITIALIZED_DOUBLE
    dvt = UNINITIALIZED_DOUBLE
    hvp = UNINITIALIZED_DOUBLE
    hvt = UNINITIALIZED_DOUBLE    
  endif  
  
end subroutine EOSWaterSteamDenEnthTPPlanar

! ************************************************************************** !

subroutine EOSWaterDuanMixture(t,p,xmol,y_nacl,avgmw,dw_kg,denmix)

! Duan et al. (2008) Energy and Fuels, v 22, 1666-1674.

  implicit none

  PetscReal :: t,tk,p,xco2,xmol,x1,y_nacl,vphi_a1,vphi_a2,vphi,denmix,pw_kg,dw_kg,avgmw

  PetscReal :: fmwh2o = 18.01534d0
  PetscReal :: fmwco2 = 44.0098d0
  PetscReal :: fmwnacl = 58.44277d0
  PetscReal :: dummy
  PetscErrorCode :: ierr

  !duan mixing **************************
  tk = t + 273.15D0; xco2 = xmol;
  call EOSWaterDensity(t,p,pw_kg,dummy,ierr)
  x1 = 1.D0-xco2;
  vphi_a1 = (0.3838402D-3*tk - 0.5595385D0)*tk + 0.30429268D3 + &
            (-0.72044305D5 + 0.63003388D7/tk)/tk;  
  vphi_a2 = (-0.57709332D-5*tk + 0.82764653D-2)*tk - 0.43813556D1 + &
            (0.10144907D4 - 0.86777045D5/tk)/tk;  
  vphi = (1.D0 + vphi_a1 + vphi_a2*p*1.D-6)*(fmwh2o*1.D-3/pw_kg); 
  vphi = x1*((1.D0 - y_nacl)*fmwh2o + y_nacl*fmwnacl)*1.D-3/dw_kg + xco2*vphi;
  denmix = (x1*((1.D0 - y_nacl)*fmwh2o + y_nacl*fmwnacl) + xco2*fmwco2)*1.D-3/vphi;
  denmix = denmix/avgmw
  
end subroutine EOSWaterDuanMixture

! ************************************************************************** !

subroutine EOSWaterViscosityNaCl (t,p_Pa,xnacl,visnacl)

  !viscosity: Kestin et al. (1981)

  implicit none

  PetscReal, intent(in) :: t        ! [C]
  PetscReal, intent(in) :: p_Pa     ! [Pa]
  PetscReal, intent(in) :: xnacl    ! [-]
  PetscReal, intent(out) :: visnacl ! [Pa-s]


  PetscReal, save :: a1,a2,a3,b1,b2,b3,c1,c2,c3,c4,wnacl
  PetscReal :: ak,bk,ck
  PetscReal :: beta,betap,betas,betaw
  PetscReal :: tt,mnacl,fac,mu0,ms
  PetscReal :: p_GPa

  data a1,a2,a3 / 3.324d-2, 3.624d-3, -1.879d-4 /
  data b1,b2,b3 / -3.96d-2, 1.02d-2, -7.02d-4 /
  data c1,c2,c3,c4 / 1.2378d0, -1.303d-3, 3.06d-6, 2.55d-8 /

  data wnacl / 58.44277d-3 / ! (kg/mol NaCl)

  !convert pressure to GPa
  p_GPa = p_Pa*1.d-9

  mnacl = xnacl/(1.d0-xnacl)/wnacl

  tt = 20.d0-t
  ck = (c1 + (c2 + (c3+c4*tt)*tt)*tt)*tt/(96.d0+t)
  ak = (a1 + (a2 + a3*mnacl)*mnacl)*mnacl
  bk = (b1 + (b2 + b3*mnacl)*mnacl)*mnacl

  ms = 6.044d0 + (2.8d-3 + 3.6d-5*t)*t
  fac = mnacl/ms
  betaw = -1.297d0 + (5.74d-2 + (-6.97d-4 + (4.47d-6 - 1.05d-8*t)*t)*t)*t
  betas = 0.545d0 + 2.8d-3 * t - betaw
  betap = (2.5d0 + (-2.d0 + 0.5d0*fac)*fac)*fac
  beta = betas*betap + betaw

  mu0 = 1001.74d-6 * 10.d0**(ak + ck*(bk + 1.d0))

  visnacl = mu0*(1.d0 + beta*p_GPa)

end subroutine EOSWaterViscosityNaCl

! ************************************************************************** !

subroutine EOSWaterViscosityKestinExt(T, P, PS, dPS_dT, aux, &
                                      calculate_derivatives, VW, &
                                      dVW_dT, dVW_dP, ierr)

  !viscosity: Kestin et al. (1981)

  implicit none

  PetscReal, intent(in) :: T   ! C
  PetscReal, intent(in) :: P   ! Pa
  PetscReal, intent(in) :: PS  ! Pa
  PetscReal, intent(in) :: dPS_dT
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW ! Pa-s
  PetscReal, intent(out) :: dVW_dT, dVW_dP
  PetscErrorCode, intent(out) :: ierr

  PetscReal, save :: a1,a2,a3,b1,b2,b3,c1,c2,c3,c4,wnacl
  PetscReal :: xnacl    ! [-]
  PetscReal :: ak,bk,ck
  PetscReal :: beta,betap,betas,betaw
  PetscReal :: tt,mnacl,fac,mu0,ms
  PetscReal :: p_GPa

  data a1,a2,a3 / 3.324d-2, 3.624d-3, -1.879d-4 /
  data b1,b2,b3 / -3.96d-2, 1.02d-2, -7.02d-4 /
  data c1,c2,c3,c4 / 1.2378d0, -1.303d-3, 3.06d-6, 2.55d-8 /

  data wnacl / 58.44277d-3 / ! (kg/mol NaCl)

  if (calculate_derivatives) then
    print *, 'Derivatives not set up in EOSWaterViscosityKestinExt().'
    stop
  endif

  !convert pressure to GPa
  p_GPa = P*1.d-9
  xnacl = aux(1)

  mnacl = xnacl/(1.d0-xnacl)/wnacl

  tt = 20.d0-t
  ck = (c1 + (c2 + (c3+c4*tt)*tt)*tt)*tt/(96.d0+t)
  ak = (a1 + (a2 + a3*mnacl)*mnacl)*mnacl
  bk = (b1 + (b2 + b3*mnacl)*mnacl)*mnacl

  ms = 6.044d0 + (2.8d-3 + 3.6d-5*t)*t
  fac = mnacl/ms
  betaw = -1.297d0 + (5.74d-2 + (-6.97d-4 + (4.47d-6 - 1.05d-8*t)*t)*t)*t
  betas = 0.545d0 + 2.8d-3 * t - betaw
  betap = (2.5d0 + (-2.d0 + 0.5d0*fac)*fac)*fac
  beta = betas*betap + betaw

  mu0 = 1001.74d-6 * 10.d0**(ak + ck*(bk + 1.d0))

  VW = mu0*(1.d0 + beta*p_GPa)

end subroutine EOSWaterViscosityKestinExt

! ************************************************************************** !

subroutine EOSWaterSaturationTemperature(ps,ts_guess,ts,t_ps,ierr)

  !  This function calculates saturation temperature for a given Ps c  
  !  Ref.: International Formulation Committee of the Sixth International
  !       Conference on Properties of Steam (1967).

  !    ps  = saturation pressure (pascals)
  !    ts  = saturation temperature (deg. C)
  !    tsp = estimated ts on entry and dT/dps on return

  implicit none

  PetscReal, intent(in) :: ps
  PetscReal, intent(in) :: ts_guess
  PetscReal, intent(out) :: ts
  PetscReal, intent(out) :: t_ps
  PetscErrorCode :: ierr
  

  PetscReal, parameter :: epsilon = 1.d-10
  PetscReal, parameter :: tc1 = H2O_CRITICAL_TEMPERATURE
  PetscReal, parameter :: pc1 = H2O_CRITICAL_PRESSURE
      
  PetscReal :: theta, beta, u1, err
  PetscReal :: t1num, t1nump
  PetscReal :: t1, t1den, t1denp, term1, term1p, t2, term2, term2p
  PetscReal :: f, fp
  PetscReal :: kn(9)
  PetscInt :: itr

  data kn / -7.691234564d0,-2.608023696d1,-1.681706546d2, &
            6.423285504d1,-1.189646225d2, 4.167117320d0, &
            2.097506760d1, 1.d9         , 6.d0/

!geh  if (ipvtab.eq.0 .or. tsp.gt.369.d0) then

!-------newton-raphson iteration for calculating ts by analytical funcs

!-------on entry, ts_guess = estimated ts
  theta = (ts_guess+273.15d0)/tc1
  beta  = ps/pc1

  u1  = 1.d0-theta
  itr = 0

  do         
    itr = itr + 1

    t1num  = u1*(kn(1)+u1*(kn(2)+u1*(kn(3)+u1*(kn(4)+u1*kn(5)))))
    t1nump = u1*(2.d0*kn(2)+u1*(3.d0*kn(3)+u1*(4.d0*kn(4)+5.d0* &
             u1*kn(5))))+kn(1)
    t1     = 1.d0+u1*(kn(6)+kn(7)*u1)
    t1den  = 1.d0/(theta*t1)
    t1denp = theta*(kn(6)+2.d0*kn(7)*u1)-t1
    term1  = t1num*t1den
    term1p = (t1nump-term1*t1denp)*t1den

    t2     = 1.d0/(kn(8)*u1*u1+kn(9))
    term2  = u1*t2
    term2p = t2*(1.d0-2.d0*kn(8)*u1*u1*t2)
    f      = exp(term1-term2)-beta
    fp     = (f+beta)*(term1p-term2p)
    err    = f/fp
    u1     = u1-err
    theta  = 1.d0-u1
    
    if (dabs(err) <= epsilon .or. itr >= 20) exit

  enddo

  ts = theta*tc1-273.15d0

!-------Note-(dbeta/dtheta) = -fp ; tsp = dT/dps
  t_ps = -tc1/(pc1*fp)

end subroutine EOSWaterSaturationTemperature

! ************************************************************************** !

subroutine EOSWaterDensityIcePainter(T, P, calculate_derivatives, &
                                     den_ice, dden_ice_dT, dden_ice_dP, ierr)
  ! Subroutine to calculate the density of ice at given temperature
  ! and pressure
  ! T is in deg C, P is in Pa, density is in kmol/m3
  implicit none
  
  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscBool, intent(in) :: calculate_derivatives

  PetscReal, intent(out) :: den_ice
  PetscReal, intent(out) :: dden_ice_dT
  PetscReal, intent(out) :: dden_ice_dP
  PetscErrorCode, intent(out) :: ierr

  PetscReal, parameter :: P_ref = 1.d5
  PetscReal, parameter :: alpha = 3.3d-10
  PetscReal, parameter :: beta = 1.53d-4

  den_ice = 5.09424d1*(1.d0 + alpha*(P - P_ref) - beta*(T)) !in Kmol/m3
  if (calculate_derivatives) then
    dden_ice_dT = 5.09424d1*(-beta)
    dden_ice_dP = 5.09424d1*alpha
  else
    dden_ice_dT = UNINITIALIZED_DOUBLE
    dden_ice_dP = UNINITIALIZED_DOUBLE
  endif
  
end subroutine EOSWaterDensityIcePainter

! ************************************************************************** !

subroutine EOSWaterInternalEnergyIce(T, u_ice, du_ice_dT)
  ! Subroutine to calculate the internal energy of ice at given temperature and 
  ! pressure
  ! T is in deg C, internal energy is in J/mol
  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(out) :: u_ice
  PetscReal, intent(out) :: du_ice_dT
  PetscErrorCode :: ierr
  
  PetscReal, parameter :: a = -10.6644d0
  PetscReal, parameter :: b = 0.1698d0
  PetscReal, parameter :: c = 198148.d0
  PetscReal, parameter :: T_ref = 273.15d0

  ! from Maier-Kelly type fit (integrated tref to t)
  ! in J/mol

  u_ice = a*(T) + b/2.d0*((T + T_ref)**(2.d0) - T_ref**(2.d0)) + &
          c*(1.d0/T_ref - 1.d0/(T + T_ref))
  u_ice = u_ice - HEAT_OF_FUSION*FMWH2O*1.d-3   ! kJ/kmol
  du_ice_dT = a + b*(T + T_ref) + c/((T + T_ref)**(2.d0)) !kJ/kmol/K
  
end subroutine EOSWaterInternalEnergyIce

! ************************************************************************** !

subroutine EOSWaterDensityPainter(t,p,calculate_derivatives,dw,dwmol, &
                                  dwp,dwt,ierr)

! wateos_simple: Simple water equation of state from Scott Painter
! Author: Satish Karra, LANL
! Date: 02/1/12
! T in C, P in Pa
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr

  PetscReal, parameter :: a = 999.915d0
  PetscReal, parameter :: b = 0.0416516d0
  PetscReal, parameter :: c = -0.0100836d0
  PetscReal, parameter :: d = 0.000206355
  PetscReal, parameter :: alpha = 5.0d-10     ! in Pa^(-1)
  PetscReal, parameter :: T_ref = 273.15d0    ! in K
  PetscReal, parameter :: P_ref = 1.0d5       ! in Pa

  PetscReal :: den_w_one_bar, T_K
  PetscReal :: u_J_kg, h_J_kg
  PetscReal :: du_dt

  ! Density of water
  T_K = T + T_ref    ! convert to Kelvin
  den_w_one_bar = a + b*(T_K - T_ref) + c*(T_K - T_ref)**(2.d0) + &
                  d*(T_K - T_ref)**(3.d0)
  dw = den_w_one_bar*(1 + alpha*(P - P_ref))
  dwmol = dw/FMWH2O     ! in mol

  ! Internal energy
  u_J_kg = 4.217*1.0d3*(T_K - T_ref)    ! in J/kg
  h_J_kg = u_J_kg + P/dw    ! in J/kg

  if (calculate_derivatives) then
    ! Derivatives of density
    dwp = 1/FMWH2O*den_w_one_bar*alpha    ! in Kmol/Pa
    dwt = 1/FMWH2O*(1 + alpha*(P - P_ref))*(b + 2.d0*c*(T_K - T_ref) + &
                              3.d0*d*(T_K - T_ref)**(2.d0))      ! in Kmol/K
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityPainter

! ************************************************************************** !

subroutine EOSWaterEnthalpyPainter(T, P, calculate_derivatives, &
                                   h_J_kmol, dh_dp, dh_dt, ierr)

! wateos_simple: Simple water equation of state from Scott Painter
! Author: Satish Karra, LANL
! Date: 02/1/12
! T in C, P in Pa
  implicit none

  PetscReal, intent(in) :: T
  PetscReal, intent(in) :: P
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: h_J_kmol
  PetscReal :: den_water_kg, den_water_kmol
  PetscReal :: dden_water_dp, dden_water_dt
  PetscReal, intent(out) :: dh_dp, dh_dt

  PetscErrorCode, intent(out) :: ierr

  PetscReal, parameter :: a = 999.915d0
  PetscReal, parameter :: b = 0.0416516d0
  PetscReal, parameter :: c = -0.0100836d0
  PetscReal, parameter :: d = 0.000206355
  PetscReal, parameter :: alpha = 5.0d-10     ! in Pa^(-1)
  PetscReal, parameter :: T_ref = 273.15d0    ! in K
  PetscReal, parameter :: P_ref = 1.0d5       ! in Pa

  PetscReal :: den_w_one_bar, T_K
  PetscReal :: u_J_kg, h_J_kg
  PetscReal :: du_dt

  ! Density of water
  T_K = T + T_ref    ! convert to Kelvin
  den_w_one_bar = a + b*(T_K - T_ref) + c*(T_K - T_ref)**(2.d0) + &
                  d*(T_K - T_ref)**(3.d0)
  den_water_kg = den_w_one_bar*(1 + alpha*(P - P_ref))
  den_water_kmol = den_water_kg/FMWH2O     ! in mol

  ! Internal energy
  u_J_kg = 4.217*1.0d3*(T_K - T_ref)    ! in J/kg
  h_J_kg = u_J_kg + P/den_water_kg    ! in J/kg
  h_J_kmol = h_J_kg*FMWH2O     ! in J/kmol

  if (calculate_derivatives) then
    ! Derivatives of density
    dden_water_dp = 1/FMWH2O*den_w_one_bar*alpha    ! in Kmol/Pa
    dden_water_dt = 1/FMWH2O*(1 + alpha*(P - P_ref))*(b + 2.d0*c*(T_K - T_ref) + &
                              3.d0*d*(T_K - T_ref)**(2.d0))      ! in Kmol/K

    ! Derivatives of enthalpy
    dh_dp = FMWH2O/den_water_kg   ! in J/kmol/Pa
    du_dt = 4.217*1.d3                  ! in J/kg/K
    dh_dt = FMWH2O*(du_dt + P*(-1.d0/den_water_kg**(2.d0))* &
                    dden_water_dt*FMWH2O)    ! in MJ/kmol/K
  else
    dden_water_dp = UNINITIALIZED_DOUBLE
    dden_water_dp = UNINITIALIZED_DOUBLE
    dh_dp = UNINITIALIZED_DOUBLE
    du_dt = UNINITIALIZED_DOUBLE
    dh_dt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterEnthalpyPainter

! ************************************************************************** !

subroutine EOSWaterDensityIceNoDerive(t,p,dw,ierr)

  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: dw
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dum1, dum2

  call EOSWaterDensityIcePtr(t,p,PETSC_FALSE,dw,dum1,dum2,ierr)

end subroutine EOSWaterDensityIceNoDerive

! ************************************************************************** !

subroutine EOSWaterDensityIceDerive(t,p,dw,dwp,dwt,ierr)
  implicit none

  PetscReal, intent(in) :: t
  PetscReal, intent(in) :: p
  PetscReal, intent(out) :: dw,dwp,dwt
  PetscErrorCode, intent(out) :: ierr

  call EOSWaterDensityIcePtr(t,p,PETSC_TRUE,dw,dwp,dwt,ierr)

end subroutine EOSWaterDensityIceDerive

! ************************************************************************** !

subroutine EOSWaterDensityTGDPB01(t, p, calculate_derivatives, &
                                  dw, dwmol, dwp, dwt, ierr)

  ! 
  ! Tanaka M. , G. Girard, R. Davis, A. Peuto, and N. Bignell. 2001.
  ! Recommended table for the density of water between 0 °C
  ! and 40 °C based on recent experimental reports. Metrologia,
  ! 38:301-309 [doi:10.1088/0026-1394/38/4/3].
  ! 
  ! Author: Gautam Bisht, LBNL
  ! Date: 24/06/15
  ! 
  implicit none

  PetscReal, intent(in) :: t   ! Temperature in centigrade
  PetscReal, intent(in) :: p   ! Pressure in Pascals
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw,dwmol,dwp,dwt
  PetscErrorCode, intent(out) :: ierr

  PetscReal,parameter :: a1 = -3.983035d0     ! [degC]
  PetscReal,parameter :: a2 = 301.797d0       ! [degC]
  PetscReal,parameter :: a3 = 522528.9d0      ! [degC^{2}]
  PetscReal,parameter :: a4 = 69.34881d0      ! [degC]
  PetscReal,parameter :: a5 = 999.974950d0    ! [kg m^{-3}]
  PetscReal,parameter :: k0 = 50.74d-11       ! [Pa^{-1}]
  PetscReal,parameter :: k1 = -0.326d-11      ! [Pa^{-1} degC^{-1}]
  PetscReal,parameter :: k2 = 0.00416d-11     ! [Pa^{-1} degC^{-2}]
  PetscReal,parameter :: p0 = 101325.d0       ! [Pa]
  PetscReal :: t_c
  PetscReal :: dent
  PetscReal :: kappa
  PetscReal :: ddent_dt
  PetscReal :: ddent_dt_1
  PetscReal :: ddent_dt_2
  PetscReal :: ddent_dt_3
  PetscReal :: ddent_dp
  PetscReal :: dkappa_dp
  PetscReal :: dkappa_dt
  PetscReal :: dden_dt

  ! Density of water as function of temperature
  dent = a5*(1.d0 - ((t + a1)**2.d0)*(t + a2)/a3/(t + a4))

  ! Compressibility of water
  kappa = (1.d0 + (k0 + k1*t + k2*t**2.d0)*(p - p0))

  ! Density of water
  dw    = dent*kappa ! [kg m^{-3}]
  dwmol = dw/FMWH2O  ! [kmol m^{-3}]

  if (calculate_derivatives) then
    ! Derivative
    ddent_dp = 0.d0
    ddent_dt_1 = -((t + a1)**2.d0)/a3/(t + a4)
    ddent_dt_2 = -2.d0*(t + a1)*(t + a2)/a3/(t + a4)
    ddent_dt_3 =  ((t + a1)**2.d0)*(t + a2)/a3/((t + a4)**2.d0)
    ddent_dt   = a5*(ddent_dt_1 + ddent_dt_2 + ddent_dt_3)

    dkappa_dp = (k0 + k1*t + k2*t**2.d0)
    dkappa_dt = (k1 + 2.d0*k2*t)*(p - p0)

    dwt = (ddent_dt*kappa + dent*dkappa_dt)/FMWH2O ! [kmol m^{-3} degC^{-1}]
    dwp = (ddent_dp*kappa + dent*dkappa_dp)/FMWH2O ! [kmol m^{-3} Pa^{-1}]
  else
    dwt = UNINITIALIZED_DOUBLE
    dwp = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityTGDPB01

! ************************************************************************** !

subroutine EOSWaterDensityBatzleAndWang(tin, pin, calculate_derivatives, &
                                        dw, dwmol, dwp, dwt, ierr)

  ! 
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Equation 27a
  !
  ! Author: Glenn Hammond
  ! Date: 02/08/16
  ! 
  implicit none

  PetscReal, intent(in) :: tin   ! Temperature in centigrade
  PetscReal, intent(in) :: pin   ! Pressure in Pascal
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw ! kg/m^3
  PetscReal, intent(out) :: dwmol ! kmol/m^3
  PetscReal, intent(out) :: dwp ! kmol/m^3-Pa
  PetscReal, intent(out) :: dwt ! kmol/m^3-C
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal, parameter :: g_cm3_to_kg_m3 = 1.d3
  PetscReal, parameter :: Pa_to_MPa = 1.d-6
  PetscReal :: t ! temperature in Celcius
  PetscReal :: t_sq, t_cub
  PetscReal :: p_MPa, p_MPa_sq
  
  t = tin
  t_sq = t*t
  t_cub = t*t_sq
  p_MPa = pin*Pa_to_MPa
  p_MPa_sq = p_MPa*p_MPa
  
  ! temperature is in C and pressure in MPa
  ! g/cm^3
  ! Eq. 27a
  dw = 1.d0 + 1.d-6*(-80.d0*t - 3.3d0*t_sq + 1.75d-3*t_cub + &
                     489.d0*p_MPa - 2.d0*t*p_MPa + 1.6d-2*t_sq*p_MPa - &
                     1.3d-5*t_cub*p_MPa - 3.33d-1*p_MPa_sq - 2.d-3*t*p_MPa_sq)
  ! convert from g/cm^3 to kg/m^3
  dw = dw * g_cm3_to_kg_m3
  dwmol = dw/FMWH2O ! kmol/m^3 
  
  if (calculate_derivatives) then
    dwp = 1.d-6*(489.d0 - 2.d0*t + 1.6d-2*t_sq - 1.3d-5*t_cub - &
                 6.66d-1*p_MPa - 4.d-3*t*p_MPa) *  &
          g_cm3_to_kg_m3/FMWH2O
    ! convert from kmol/m^3-MPa to kmol/m^3-Pa
    dwp = dwp*Pa_to_MPa
    dwt = 1.d-6*(-80.d0 - 6.6d0*t + 5.25d-3*t_sq - 2.d0*p_MPa + &
                 3.2d-2*t*p_MPa - 3.9d-5*t_sq*p_MPa - 2.d-3*p_MPa_sq) * &
          g_cm3_to_kg_m3/FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityBatzleAndWang

! ************************************************************************** !

subroutine EOSWaterDensityBatzleAndWangExt(tin, pin, aux, &
                                           calculate_derivatives, &
                                           dw, dwmol, dwp, dwt, ierr)

  ! 
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Equation 27b
  !
  ! Author: Glenn Hammond
  ! Date: 02/08/16
  ! 
  implicit none

  PetscReal, intent(in) :: tin   ! Temperature in centigrade
  PetscReal, intent(in) :: pin   ! Pressure in Pascal
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: dw ! kg/m^3
  PetscReal, intent(out) :: dwmol ! kmol/m^3
  PetscReal, intent(out) :: dwp ! kmol/m^3-Pa
  PetscReal, intent(out) :: dwt ! kmol/m^3-C
  PetscErrorCode, intent(out) :: ierr
  
  PetscReal, parameter :: g_cm3_to_kg_m3 = 1.d3
  PetscReal, parameter :: Pa_to_MPa = 1.d-6
  PetscReal :: t_C ! temperature in Celcius
  PetscReal :: p_MPa
  PetscReal :: s
  
  t_C = tin
  p_MPa = pin*Pa_to_MPa
  s = aux(1)
  
  call EOSWaterDensityPtr(tin, pin, calculate_derivatives, &
                          dw, dwmol, dwp, dwt, ierr)
  
  ! temperature is in C and pressure in MPa
  ! kg/m^3
  ! Eq. 27b
  dw = dw + &
       s*(0.668d0 + 0.44d0*s + &
          1.d-6*(300.d0*p_MPa - 2400.d0*p_MPa*s + &
                 t_C*(80.d0 + 3.d0*t_C - 3300.d0*s - 13.d0*p_MPa + &
                      47.d0*p_Mpa*s))) * &
       g_cm3_to_kg_m3
  dwmol = dw/FMWH2O ! kmol/m^3 
  
  if (calculate_derivatives) then
        ! v - this dwp is in the correct units of kmol/m^3-Pa
    dwp = dwp + &
          s*(1.d-6*(300.d0 - 2400.d0*s + t_C*(-13.d0 + 47.d0*s))) * &
                                 ! v - convert from kmol/m^3-MPa to kmol/m^3-Pa
          g_cm3_to_kg_m3/FMWH2O*Pa_to_MPa
    dwt = dwt + &
          s*(1.d-6*(80.d0 + 6.d0*t_C - 3300.d0*s - 13.d0*p_MPa + &
                    47.d0*p_Mpa*s)) * &
          g_cm3_to_kg_m3/FMWH2O
  else
    dwp = UNINITIALIZED_DOUBLE
    dwt = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterDensityBatzleAndWangExt

! ************************************************************************** !

subroutine EOSWaterViscosityBatzleAndWang(T, P, PS, dPS_dT, &
                                          calculate_derivatives, VW, &
                                          dVW_dT, dVW_dP, ierr)
  ! 
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Equation 32
  !
  ! Author: Glenn Hammond
  ! Date: 02/08/16
  ! 
  implicit none
  PetscReal, intent(in) :: T       ! C
  PetscReal, intent(in) :: P       ! Pa
  PetscReal, intent(in) :: PS      ! Pa
  PetscReal, intent(in) :: dPS_dT  ! Pa/C
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW     ! Pa-s
  PetscReal, intent(out) :: dVW_dT, dVW_dP
  PetscErrorCode, intent(out) :: ierr
  
  ! convert from centipoise to Pa-s (1 cP = 1.d-3 Pa-s)
  PetscReal, parameter :: centipoise_to_Pa_s = 1.d-3
  PetscReal :: t_C
  PetscReal :: exponential_term
  PetscReal :: temperature_term
  
  t_C = T
  
  ! this is Eq. 32 without all the salt terms.
  ! -0.057138d0 = -1.d0*0.42d0*(-0.17d0)**2.d0+0.045d0  
  exponential_term = -0.057138d0*t_C**0.8d0
  temperature_term = 1.65d0*exp(exponential_term)
  VW = 0.1d0 + temperature_term
  VW = VW * centipoise_to_Pa_s
       
  if (calculate_derivatives) then
    dVW_dP = 0.d0
    dVW_dT = 0.8d0*temperature_term*exponential_term/t_C*centipoise_to_Pa_s
  else
    dVW_dP = UNINITIALIZED_DOUBLE
    dVW_dT = UNINITIALIZED_DOUBLE
  endif

end subroutine EOSWaterViscosityBatzleAndWang

! ************************************************************************** !

subroutine EOSWaterViscosityBatzleAndWangExt(T, P, PS, dPS_dT, aux, &
                                             calculate_derivatives, VW, &
                                             dVW_dT, dVW_dP, ierr)
  ! 
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Equation 32
  !
  ! Author: Glenn Hammond
  ! Date: 02/08/16
  ! 
  implicit none

  PetscReal, intent(in) :: T, P, PS, dPS_dT
  PetscReal, intent(in) :: aux(*)
  PetscBool, intent(in) :: calculate_derivatives
  PetscReal, intent(out) :: VW
  PetscReal, intent(out) :: dVW_dT, dVW_dP
  PetscErrorCode, intent(out) :: ierr
  
  ! convert from centipoise to Pa-s (1 cP = 1.d-3 Pa-s)
  PetscReal, parameter :: centipoise_to_Pa_s = 1.d-3
  PetscReal :: t_C
  PetscReal :: s
  PetscReal :: exponential_term
  PetscReal :: temperature_term
  
  s = aux(1)
  t_C = T
  
  exponential_term = -1.d0*(0.42d0*(s**0.8d0-0.17d0)**2.d0 + 0.045d0)* &
                     t_C**0.8d0
  temperature_term = (1.65d0 + 91.9d0*s**3.d0)*exp(exponential_term)
  VW = 0.1d0 + 0.333d0*s + temperature_term
  VW = VW * centipoise_to_Pa_s
       
  if (calculate_derivatives) then
    dVW_dP = 0.d0
    dVW_dT = 0.8d0*temperature_term*exponential_term/t_C*centipoise_to_Pa_s
  else
    dVW_dP = UNINITIALIZED_DOUBLE
    dVW_dT = UNINITIALIZED_DOUBLE
  endif
  
end subroutine EOSWaterViscosityBatzleAndWangExt

! ************************************************************************** !

subroutine TestEOSWaterBatzleAndWang()
  ! 
  ! From Batlze M. and Z. Wang (1992) Seismic properties of fluids, Geophysics,
  ! Vol. 57, No. 11, Pg. 1396-1408.
  !
  ! Code for testing
  !
  ! Author: Glenn Hammond
  ! Date: 02/15/16
  ! 
  PetscReal :: p, t, dw, dwmol, dwp, dwt, ps, dps_dt, vw, dvw_dt, dvw_dp
  PetscReal :: t2, dw2, p2, vw2
  PetscReal :: aux(1)
  PetscErrorCode :: ierr
  
  p = 0.101325d6
  t = 20.d0
  ps = -999.d0
  dps_dt = -999.d0
  aux(1) = 0.14d0
  
  call EOSWaterSetDensity('BATZLE_AND_WANG')
  call EOSWaterSetViscosity('BATZLE_AND_WANG')
  
  call EOSWaterDensityBatzleAndWang(t,p, PETSC_TRUE, &
                                    dw, dwmol, dwp, dwt, ierr)
  print *, 'dw:    ', dw
  print *, 'dwmol: ', dwmol
  print *, 'dwp:   ', dwp
  print *, 'dwt:   ', dwt

  t2 = t + 1.d-6*t
  call EOSWaterDensityBatzleAndWang(t2,p, PETSC_TRUE, &
                                    dw2, dwmol, dwp, dwt, ierr)
  print *
  print *, 'dw2:   ', dw2
  print *, 'dw(t): ', (dw2-dw)/(t2-t)
  
  p2 = p + 1.d-6*p
  call EOSWaterDensityBatzleAndWang(t,p2, PETSC_TRUE, &
                                    dw2, dwmol, dwp, dwt, ierr)
  print *
  print *, 'dw2:   ', dw2
  print *, 'dw(p): ', (dw2-dw)/(p2-p)
  
  
  call EOSWaterDensityBatzleAndWangExt(t,p,aux, PETSC_TRUE, &
                                    dw, dwmol, dwp, dwt, ierr)
  print *, 'Density-Ext'
  print *, 'dw:    ', dw
  print *, 'dwmol: ', dwmol
  print *, 'dwp:   ', dwp
  print *, 'dwt:   ', dwt

  call EOSWaterViscosityBatzleAndWang(t, p, PS, dPS_dT, &
                                      PETSC_TRUE, vw, &
                                      dvw_dt, dvw_dp, ierr)  
  print *
  print *, 'vw:      ', vw
  print *, 'dvw_dp:  ', dvw_dp
  print *, 'dvw_dt:  ', dvw_dt
  
  
  call EOSWaterViscosityBatzleAndWangExt(t, p, PS, dPS_dT, aux, &
                                         PETSC_TRUE, vw, &
                                         dvw_dt, dvw_dp, ierr) 
  print *, 'Ext-'
  print *, 'vw:      ', vw
  print *, 'dvw_dp:  ', dvw_dp
  print *, 'dvw(t)t:  ', dvw_dt
  
  call EOSWaterViscosityBatzleAndWangExt(t2, p, PS, dPS_dT, aux, &
                                         PETSC_TRUE, vw2, &
                                         dvw_dt, dvw_dp, ierr) 
  
  print *, 'Ext-numerical'
  print *, 'vw:      ', vw2
  print *, 'dvw(t)t:  ', (vw2-vw)/(t2-t)
  
  call EOSWaterViscosityBatzleAndWangExt(t, p, PS, dPS_dT, aux, &
                                         PETSC_TRUE, vw, &
                                         dvw_dt, dvw_dp, ierr) 
  print *, 'Ext-S'
  print *, 'vw:      ', vw
  print *, 'dvw_dp:  ', dvw_dp
  print *, 'dvw(t)t:  ', dvw_dt
  
end subroutine TestEOSWaterBatzleAndWang

! ************************************************************************** !

subroutine EOSWaterDensityExtNumericalDerive(t,p,aux,dw,dwmol,dwp,dwt,ierr)

  implicit none

  PetscReal, intent(in) :: t     ! Temperature in centigrade
  PetscReal, intent(in) :: p     ! Pressure in Pascal
  PetscReal, intent(in) :: aux(*)
  PetscReal, intent(out) :: dw ! kg/m^3
  PetscReal, intent(out) :: dwmol ! kmol/m^3
  PetscReal, intent(out) :: dwp ! kmol/m^3-Pa
  PetscReal, intent(out) :: dwt ! kmol/m^3-C
  PetscErrorCode, intent(out) :: ierr

  PetscReal :: dwp_analytical, dwt_analytical
  PetscReal :: dwp_numerical, dwt_numerical
  PetscReal :: dw_no_salt, dwp_no_salt, dwt_no_salt
  PetscReal :: dum1, dum2, dum3
  PetscReal :: dwmol_tpert, dwmol_ppert
  PetscReal :: tpert, ppert, t_plus_tpert, p_plus_ppert
  PetscReal :: salinity(1)
  PetscReal :: pert_tol = 1.d-6
  
  tpert = t*pert_tol
  ppert = p*pert_tol
  t_plus_tpert = t + tpert
  p_plus_ppert = p + ppert
  salinity(1) = aux(1)

#if 0 
  ! test against non-extended version
  call EOSWaterDensityPtr(t,p,PETSC_TRUE,dw,dwmol,dwp_analytical, &
                          dwt_analytical,ierr)
  dwp = dwp_analytical
  dwt = dwt_analytical
#else
  call EOSWaterDensityExtPtr(t,p,salinity,PETSC_TRUE, &
                             dw,dwmol,dwp_analytical,dwt_analytical,ierr)
  dwp = dwp_analytical
  dwt = dwt_analytical
  call EOSWaterDensityExtPtr(t_plus_tpert,p,salinity,PETSC_FALSE, &
                             dum1,dwmol_tpert,dum2,dum3,ierr)
  call EOSWaterDensityExtPtr(t,p_plus_ppert,salinity,PETSC_FALSE, &
                             dum1,dwmol_ppert,dum2,dum3,ierr)

  dwp_numerical = (dwmol_ppert-dwmol)/ppert
  dwt_numerical = (dwmol_tpert-dwmol)/tpert

  if (.not.PETSC_TRUE) then
    dwp = dwp_numerical
    dwt = dwt_numerical
  else
    dwp = dwp_analytical
    dwt = dwt_analytical
  endif

  if (dabs((dwp_numerical-dwp_analytical)/dwp_numerical) > 1.d-4) then
    print *, p, t, salinity(1), dw, dwmol, dwp, dwp_analytical, dwp_numerical
  endif
#endif
  
end subroutine EOSWaterDensityExtNumericalDerive

! **************************************************************************** !

subroutine EOSWaterInputRecord()
  ! 
  ! Prints ingested equation of state information to the input record file.
  ! 
  ! Author: Jenn Frederick
  ! Date: 05/04/2016
  ! 
  
  implicit none
  
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: id = INPUT_RECORD_UNIT

  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'WATER'
  
  ! water density [kg/m^3]
  if (associated(EOSWaterDensityPtr,EOSWaterDensityConstant)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(word1,*) constant_density
    write(id,'(a)') 'constant, ' // adjustl(trim(word1)) // ' kg/m^3'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityExponential)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'exponential'
    write(id,'(a29)',advance='no') 'exp. ref. density: '
    write(word1,*) exponent_reference_density
    write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    write(id,'(a29)',advance='no') 'exp. ref. pressure: '
    write(word1,*) exponent_reference_pressure
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    write(id,'(a29)',advance='no') 'exp. water compressibility: '
    write(word1,*) exponent_water_compressibility
    write(id,'(a)') adjustl(trim(word1)) // ' 1/Pa'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityLinear)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'linear'
    write(id,'(a29)',advance='no') 'linear ref. density: '
    write(word1,*) linear_reference_density
    write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    write(id,'(a29)',advance='no') 'linear ref. pressure: '
    write(word1,*) linear_reference_pressure
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    write(id,'(a29)',advance='no') 'linear water compressibility: '
    write(word1,*) linear_water_compressibility
    write(id,'(a)') adjustl(trim(word1)) // ' 1/Pa'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityBRAGFLO)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'BRAGFLO'
    write(id,'(a29)',advance='no') 'exp. ref. density: '
    write(word1,*) exponent_reference_density
    write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    write(id,'(a29)',advance='no') 'exp. ref. pressure: '
    write(word1,*) exponent_reference_pressure
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    write(id,'(a29)',advance='no') 'exp. water compressibility: '
    write(word1,*) exponent_water_compressibility
    write(id,'(a)') adjustl(trim(word1)) // ' 1/Pa'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityIFC67)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'default, IFC67'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityTGDPB01)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'TGDPB01'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityPainter)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'PAINTER'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityBatzleAndWang) .and. &
      associated(EOSWaterDensityExtPtr,EOSWaterDensityBatzleAndWangExt)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'Batzle and Wang'
  endif
  if (associated(EOSWaterDensityPtr,EOSWaterDensityQuadratic)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'Quadratic'
    write(id,'(a29)',advance='no') 'quad. ref. density: '
    write(word1,*) quadratic_reference_density
    write(id,'(a)') adjustl(trim(word1)) // ' kg/m^3'
    write(id,'(a29)',advance='no') 'quad. ref. pressure: '
    write(word1,*) quadratic_reference_pressure
    write(id,'(a)') adjustl(trim(word1)) // ' Pa'
    write(id,'(a29)',advance='no') 'quad. water compressibility: '
    write(word1,*) quadratic_wat_compressibility
    write(id,'(a)') adjustl(trim(word1)) // ' 1/Pa'
  end if
  if (associated(EOSWaterDensityPtr,EOSWaterDensityTrangenstein)) then
    write(id,'(a29)',advance='no') 'water density: '
    write(id,'(a)') 'Trangenstein'
  end if
  
  ! water viscosity [Pa-s]
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosityConstant)) then
    write(id,'(a29)',advance='no') 'water viscosity: '
    write(word1,*) constant_viscosity
    write(id,'(a)') 'constant, ' // trim(word1) // ' Pa-sec'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosity1)) then
    write(id,'(a29)',advance='no') 'water viscosity: '
    write(id,'(a)') 'default'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosityBatzleAndWang)) then
    write(id,'(a29)',advance='no') 'water viscosity: '
    write(id,'(a)') 'Batzle and Wang'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosityGrabowski)) then
    write(id,'(a29)',advance='no') 'water viscosity: '
    write(id,'(a)') 'Grabowski Batzle'
  endif
  
  ! water enthalpy [J/kmol]
  if (associated(EOSWaterViscosityPtr,EOSWaterEnthalpyConstant)) then
    write(id,'(a29)',advance='no') 'water enthalpy: '
    write(word1,*) constant_enthalpy
    write(id,'(a)') 'constant, ' // trim(word1) // ' J/kmol'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterEnthalpyIFC67)) then
    write(id,'(a29)',advance='no') 'water enthalpy: '
    write(id,'(a)') 'default, IFC67'
  endif
  if (associated(EOSWaterViscosityPtr,EOSWaterEnthalpyPainter)) then
    write(id,'(a29)',advance='no') 'water enthalpy: '
    write(id,'(a)') 'PAINTER'
  endif
  
  ! steam density
  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDenEnthConstant)) then
    write(id,'(a29)',advance='no') 'steam density: '
    write(word1,*) constant_steam_density
    write(id,'(a)') 'constant, ' // trim(word1) // ' kg/m^3'
  endif
  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDensityEnthalpyIFC67)) then
    write(id,'(a29)',advance='no') 'steam density: '
    write(id,'(a)') 'default, IFC67'
  endif
  
  ! steam enthalpy [J/kmol]
  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDenEnthConstant)) then
    write(id,'(a29)',advance='no') 'steam enthalpy: '
    write(word1,*) constant_steam_enthalpy
    write(id,'(a)') 'constant, ' // trim(word1) // ' J/kmol'
  endif
  if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                 EOSWaterSteamDensityEnthalpyIFC67)) then
    write(id,'(a29)',advance='no') 'steam enthalpy: '
    write(id,'(a)') 'default, IFC67'
  endif
  
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  
end subroutine EOSWaterInputRecord

! ************************************************************************** !

subroutine EOSWaterTest(temp_low,temp_high,pres_low,pres_high, &
                        ntemp,npres,uniform_temp,uniform_pres,filename)

  implicit none

  PetscReal :: temp_low
  PetscReal :: temp_high
  PetscReal :: pres_low
  PetscReal :: pres_high
  PetscInt :: npres
  PetscInt :: ntemp
  PetscBool :: uniform_temp
  PetscBool :: uniform_pres
  character(len=MAXWORDLENGTH) :: filename

  PetscReal, allocatable :: temp(:)
  PetscReal, allocatable :: pres(:)
  PetscReal, allocatable :: density_kg(:,:)
  PetscReal, allocatable :: enthalpy(:,:)
  PetscReal, allocatable :: viscosity(:,:)
  PetscReal, allocatable :: saturation_pressure_array(:)
  PetscReal :: dum1, dum2, dum3, dum4
  PetscInt :: itemp, ipres
  PetscReal :: ln_low, ln_high
  PetscReal :: saturation_pressure
  PetscReal :: NaN
  character(len=MAXWORDLENGTH) :: eos_density_name
  character(len=MAXWORDLENGTH) :: eos_enthalpy_name
  character(len=MAXWORDLENGTH) :: eos_viscosity_name
  character(len=MAXWORDLENGTH) :: eos_saturation_pressure_name
  character(len=MAXSTRINGLENGTH) :: header, string

  PetscErrorCode :: ierr

  NaN = 0.d0
  NaN = 1.d0/NaN
  NaN = 0.d0*NaN

  allocate(temp(ntemp))
  temp = UNINITIALIZED_DOUBLE
  allocate(pres(ntemp))
  pres = UNINITIALIZED_DOUBLE
  allocate(density_kg(npres,ntemp))
  density_kg = UNINITIALIZED_DOUBLE
  allocate(viscosity(npres,ntemp))
  viscosity = UNINITIALIZED_DOUBLE
  allocate(enthalpy(npres,ntemp))
  enthalpy = UNINITIALIZED_DOUBLE
  allocate(saturation_pressure_array(ntemp))
  saturation_pressure_array = UNINITIALIZED_DOUBLE

  if (uniform_pres) then
    do ipres = 1, npres
      pres(ipres) = &
        (pres_high-pres_low)/max(dble(npres-1),1.d0) * (ipres-1) + pres_low
    enddo
  else
    ln_high = log(pres_high)
    ln_low = log(pres_low)
    do ipres = 1, npres
      pres(ipres) = &
        exp((ln_high-ln_low)/max(dble(npres-1),1.d0) * (ipres-1) + ln_low)
    enddo
  endif

  if (uniform_temp) then
    do itemp = 1, ntemp
      temp(itemp) = &
        (temp_high-temp_low)/max(dble(ntemp-1),1.d0) * (itemp-1) + temp_low
    enddo
  else
    ln_high = log(temp_high)
    ln_low = log(temp_low)
    do itemp = 1, ntemp
      temp(itemp) = &
        exp((ln_high-ln_low)/max(dble(ntemp-1),1.d0) * (itemp-1) + ln_low)
    enddo
  endif

  ! density
  if (associated(EOSWaterDensityPtr,EOSWaterDensityConstant)) then
    eos_density_name = 'Constant'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityExponential)) then
    eos_density_name = 'Exponential'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityLinear)) then
    eos_density_name = 'Linear'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityBRAGFLO)) then
    eos_density_name = 'BRAGFLO'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityIFC67)) then
    eos_density_name = 'IFC67'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityTGDPB01)) then
    eos_density_name = 'TGDPB01'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityPainter)) then
    eos_density_name = 'Painter'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityBatzleAndWang)) then
    eos_density_name = 'Batzle and Wang'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityQuadratic)) then
    eos_density_name = 'Quadratic'
  else if (associated(EOSWaterDensityPtr,EOSWaterDensityTrangenstein)) then
    eos_density_name = 'Trangenstein'
  else 
    eos_density_name = 'Unknown'
  endif

  ! enthalpy
  if (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyConstant)) then
    eos_enthalpy_name = 'Constant'
  else if (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyIFC67)) then
    eos_enthalpy_name = 'IFC67'
  else if (associated(EOSWaterEnthalpyPtr,EOSWaterEnthalpyPainter)) then
    eos_enthalpy_name = 'Painter'
  else
    eos_enthalpy_name = 'Unknown'
  endif

  ! viscosity
  if (associated(EOSWaterViscosityPtr,EOSWaterViscosityConstant)) then
    eos_viscosity_name = 'Constant'
  else if (associated(EOSWaterViscosityPtr,EOSWaterViscosity1)) then
    eos_viscosity_name = 'Default'
  else if (associated(EOSWaterViscosityPtr,EOSWaterViscosityBatzleAndWang)) then
    eos_viscosity_name = 'Batzle and Wang'
  else if (associated(EOSWaterViscosityPtr,EOSWaterViscosityGrabowski)) then
    eos_viscosity_name = 'Grabowski'
  else
    eos_viscosity_name = 'Unknown'
  endif

  ! saturation pressure
  if (associated(EOSWaterSaturationPressurePtr, &
                 EOSWaterSaturationPressureIFC67)) then
    eos_saturation_pressure_name = 'IFC67'
  else
    eos_saturation_pressure_name = 'Unknown'
  endif

  do itemp = 1, ntemp
    do ipres = 1, npres
      ! EOSWaterSaturationPressurePtr() must come before call to 
      ! EOSWaterViscosityPtr()
      call EOSWaterSaturationPressurePtr(temp(itemp),PETSC_FALSE, &
                                         saturation_pressure,dum1,ierr)
      if (ipres == 1) &
        saturation_pressure_array(itemp) = saturation_pressure
      call EOSWaterDensityPtr(temp(itemp),pres(ipres),PETSC_FALSE, &
                              density_kg(ipres,itemp), &
                              dum1,dum2,dum3,ierr)
      call EOSWaterEnthalpyPtr(temp(itemp),pres(ipres),PETSC_FALSE, &
                               enthalpy(ipres,itemp),dum1,dum2,ierr)
      call EOSWaterViscosityPtr(temp(itemp),pres(ipres),saturation_pressure, &
                                dum1,PETSC_FALSE,viscosity(ipres,itemp), &
                                dum2,dum3,ierr)
    enddo
  enddo

100 format(100es16.8)
  if (len_trim(filename) == 0) then
    string = 'eos_water_test.txt'
  else
    string = filename
  endif
  open(unit=IUNIT_TEMP,file=string)
  header = 'T[C], P[Pa], &
    &Density (' // trim(eos_density_name) // ') [kg/m^3], &
    &Enthalpy (' // trim(eos_enthalpy_name) // ') [J/kmol], &
    &Viscosity (' // trim(eos_viscosity_name) // ') [Pa-s], &
    &Saturation Pressure (' // trim(eos_saturation_pressure_name) // ') [Pa]'
  write(IUNIT_TEMP,'(a)') trim(header)
  write(IUNIT_TEMP,'(100i9)') ntemp, npres
  do itemp = 1, ntemp
    do ipres = 1, npres
      write(IUNIT_TEMP,100) temp(itemp), pres(ipres), &
            density_kg(ipres,itemp), enthalpy(ipres,itemp), &
            viscosity(ipres,itemp), saturation_pressure_array(itemp)
    enddo
  enddo
  close(IUNIT_TEMP)

  deallocate(temp)
  deallocate(pres)
  deallocate(density_kg)
  deallocate(enthalpy)
  deallocate(viscosity)
  deallocate(saturation_pressure_array)

end subroutine EOSWaterTest

! ************************************************************************** !

subroutine EOSWaterSteamTest(temp_low,temp_high,pres_low,pres_high, &
                             ntemp,npres,uniform_temp,uniform_pres,filename)

  implicit none

  PetscReal :: temp_low
  PetscReal :: temp_high
  PetscReal :: pres_low
  PetscReal :: pres_high
  PetscInt :: npres
  PetscInt :: ntemp
  PetscBool :: uniform_temp
  PetscBool :: uniform_pres
  character(len=MAXWORDLENGTH) :: filename

  PetscReal, allocatable :: temp(:)
  PetscReal, allocatable :: pres(:)
  PetscReal, allocatable :: density_kg(:,:)
  PetscReal, allocatable :: density_kg_psat(:)
  PetscReal, allocatable :: enthalpy(:,:)
  PetscReal, allocatable :: enthalpy_psat(:)
  PetscReal, allocatable :: viscosity(:,:)
  PetscReal, allocatable :: saturation_pressure_array(:)
  PetscReal :: dum1, dum2, dum3, dum4, dum5
  PetscInt :: itemp, ipres
  PetscReal :: ln_low, ln_high
  PetscReal :: saturation_pressure
  PetscReal :: NaN
  character(len=MAXWORDLENGTH) :: eos_density_name
  character(len=MAXWORDLENGTH) :: eos_enthalpy_name
  character(len=MAXWORDLENGTH) :: eos_saturation_pressure_name
  character(len=MAXSTRINGLENGTH) :: header, string

  PetscErrorCode :: ierr

  NaN = 0.d0
  NaN = 1.d0/NaN
  NaN = 0.d0*NaN

  allocate(temp(ntemp))
  temp = UNINITIALIZED_DOUBLE
  allocate(pres(ntemp))
  pres = UNINITIALIZED_DOUBLE
  allocate(density_kg(npres,ntemp))
  density_kg = UNINITIALIZED_DOUBLE
  allocate(density_kg_psat(ntemp))
  density_kg_psat = UNINITIALIZED_DOUBLE
  allocate(viscosity(npres,ntemp))
  viscosity = UNINITIALIZED_DOUBLE
  allocate(enthalpy(npres,ntemp))
  enthalpy = UNINITIALIZED_DOUBLE
  allocate(enthalpy_psat(ntemp))
  enthalpy_psat = UNINITIALIZED_DOUBLE
  allocate(saturation_pressure_array(ntemp))
  saturation_pressure_array = UNINITIALIZED_DOUBLE

  if (uniform_pres) then
    do ipres = 1, npres
      pres(ipres) = (pres_high-pres_low)/dble(npres-1) * (ipres-1) + pres_low
    enddo
  else
    ln_high = log(pres_high)
    ln_low = log(pres_low)
    do ipres = 1, npres
      pres(ipres) = exp((ln_high-ln_low)/dble(npres-1) * (ipres-1) + ln_low)
    enddo
  endif

  if (uniform_temp) then
    do itemp = 1, ntemp
      temp(itemp) = (temp_high-temp_low)/dble(ntemp-1) * (itemp-1) + temp_low
    enddo
  else
    ln_high = log(temp_high)
    ln_low = log(temp_low)
    do itemp = 1, ntemp
      temp(itemp) = exp((ln_high-ln_low)/dble(ntemp-1) * (itemp-1) + ln_low)
    enddo
  endif

  ! density
  if (associated(EOSWaterSteamDensityEnthalpyPtr,EOSWaterSteamDenEnthConstant)) then
    eos_density_name = 'Constant'
  else if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                      EOSWaterSteamDensityEnthalpyIFC67)) then
    eos_density_name = 'IFC67'
  else 
    eos_density_name = 'Unknown'
  endif

  ! enthalpy
  if (associated(EOSWaterSteamDensityEnthalpyPtr,EOSWaterSteamDenEnthConstant)) then
    eos_enthalpy_name = 'Constant'
  else if (associated(EOSWaterSteamDensityEnthalpyPtr, &
                      EOSWaterSteamDensityEnthalpyIFC67)) then
    eos_enthalpy_name = 'IFC67'
  else
    eos_enthalpy_name = 'Unknown'
  endif

  ! saturation pressure
  if (associated(EOSWaterSaturationPressurePtr, &
                 EOSWaterSaturationPressureIFC67)) then
    eos_saturation_pressure_name = 'IFC67'
  else
    eos_saturation_pressure_name = 'Unknown'
  endif

  do itemp = 1, ntemp
    do ipres = 1, npres
      ! EOSWaterSaturationPressurePtr() must come before call to 
      ! EOSWaterViscosityPtr()
      call EOSWaterSaturationPressurePtr(temp(itemp),PETSC_FALSE, &
                                         saturation_pressure,dum1,ierr)
      if (ipres == 1) then
        saturation_pressure_array(itemp) = saturation_pressure
        call EOSWaterSteamDensityEnthalpyPtr(temp(itemp),saturation_pressure, &
                                             PETSC_FALSE, &
                                             density_kg_psat(itemp),dum1, &
                                             enthalpy_psat(itemp), &
                                             dum2,dum3,dum4,dum5,ierr)
      endif
      call EOSWaterSteamDensityEnthalpyPtr(temp(itemp),pres(ipres),PETSC_FALSE, &
                                           density_kg(ipres,itemp),dum1, &
                                           enthalpy(ipres,itemp), &
                                           dum2,dum3,dum4,dum5,ierr)
    enddo
  enddo

100 format(100es16.8)
  if (len_trim(filename) == 0) then
    string = 'eos_water_steam_test.txt'
  else
    string = filename
  endif
  open(unit=IUNIT_TEMP,file=string)
  header = 'T[C], P[Pa], &
    &Steam Density (' // trim(eos_density_name) // ') [kg/m^3], &
    &Steam Enthalpy (' // trim(eos_enthalpy_name) // ') [J/kmol]'
  write(IUNIT_TEMP,'(a)') trim(header)
  write(IUNIT_TEMP,'(100i9)') ntemp, npres
  do itemp = 1, ntemp
    do ipres = 1, npres
      write(IUNIT_TEMP,100) temp(itemp), pres(ipres), &
            density_kg(ipres,itemp), enthalpy(ipres,itemp)
    enddo
  enddo
  close(IUNIT_TEMP)

  if (len_trim(filename) == 0) then
    string = 'eos_water_steam_psat_test.txt'
  else
    itemp = index(filename,'.')
    string = filename
    string(itemp:itemp+4) = '_psat'
    string(itemp+5:) = filename(itemp:)
  endif
  open(unit=IUNIT_TEMP,file=string)
  header = 'T[C], &
    &Water Saturation Pressure (' // trim(eos_saturation_pressure_name) // ') [Pa], &
    &Steam Density @ Saturation Pressure (' // trim(eos_density_name) // ') [kg/m^3], &
    &Steam Enthalpy @ Saturation Pressure (' // trim(eos_enthalpy_name) // ') [J/kmol]'
  write(IUNIT_TEMP,'(a)') trim(header)
  write(IUNIT_TEMP,'(100i9)') ntemp
  do itemp = 1, ntemp
    write(IUNIT_TEMP,100) temp(itemp), &
          saturation_pressure_array(itemp), &
          density_kg_psat(itemp), enthalpy_psat(itemp)
  enddo
  close(IUNIT_TEMP)

  deallocate(temp)
  deallocate(pres)
  deallocate(density_kg)
  deallocate(enthalpy)
  deallocate(density_kg_psat)
  deallocate(enthalpy_psat)
  deallocate(saturation_pressure_array)

end subroutine EOSWaterSteamTest

! ************************************************************************** !

end module EOS_Water_module
