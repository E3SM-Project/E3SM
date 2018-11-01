module Auxiliary_module
  
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Global_Aux_module
  use TH_Aux_module
  use Reactive_Transport_Aux_module
  use Material_Aux_class
  use Secondary_Continuum_Aux_module
  use InlineSurface_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

  type, public :: auxiliary_type 
    type(global_type), pointer :: Global
    type(reactive_transport_type), pointer :: RT
    type(th_type), pointer :: TH

    type(material_type), pointer :: Material
    type(sc_heat_type), pointer :: SC_heat
    type(sc_rt_type), pointer :: SC_RT
    type(inlinesurface_type), pointer :: InlineSurface
  end type auxiliary_type
  
  public :: AuxInit, &
            AuxDestroy

contains

! ************************************************************************** !

subroutine AuxInit(aux)
  ! 
  ! Nullifies pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/09/08
  ! 

  implicit none
  
  type(auxiliary_type) :: aux
  
  nullify(aux%Global)
  nullify(aux%RT)
  nullify(aux%TH)
  
  nullify(aux%Material)
  nullify(aux%SC_heat)
  nullify(aux%SC_RT)
  nullify(aux%InlineSurface)

end subroutine AuxInit

! ************************************************************************** !

subroutine AuxDestroy(aux)
  ! 
  ! Deallocates any allocated pointers in auxiliary object
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/09/08
  ! 

  implicit none
  
  type(auxiliary_type) :: aux
  
  call GlobalAuxDestroy(aux%Global)
  call RTAuxDestroy(aux%RT)
  call THAuxDestroy(aux%TH)

  call MaterialAuxDestroy(aux%Material)
  call SecondaryAuxHeatDestroy(aux%SC_heat)
  call SecondaryAuxRTDestroy(aux%SC_RT)
  call InlineSurfaceAuxDestroy(aux%InlineSurface)
  
  nullify(aux%Global)
  nullify(aux%RT)

  nullify(aux%Material)
  nullify(aux%SC_Heat)
  nullify(aux%SC_RT)
  nullify(aux%InlineSurface)

end subroutine AuxDestroy

end module Auxiliary_module
