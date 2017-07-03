module Auxiliary_module
  
  use Global_Aux_module
  use TH_Aux_module
  use Richards_Aux_module
  use Reactive_Transport_Aux_module
  use Mphase_Aux_module
  use Immis_Aux_module
  use Miscible_Aux_module
  use Flash2_Aux_module
  use General_Aux_module
  !use TOilIms_Aux_module
  use Material_Aux_class
  use Secondary_Continuum_Aux_module
  use InlineSurface_Aux_module
  
  !use PM_Base_Aux_module  !new auxvar data structure
  use PM_TOilIms_Aux_module  !new auxvar data structure  

  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  type, public :: auxiliary_type 
    type(global_type), pointer :: Global
    type(reactive_transport_type), pointer :: RT
    type(th_type), pointer :: TH
    type(richards_type), pointer :: Richards
    type(mphase_type), pointer :: Mphase
    type(immis_type), pointer :: Immis
    type(miscible_type), pointer :: Miscible
    type(flash2_type), pointer :: Flash2
    type(general_type), pointer :: General
    type(material_type), pointer :: Material
    type(sc_heat_type), pointer :: SC_heat
    type(sc_rt_type), pointer :: SC_RT
    !type(toil_ims_type), pointer :: TOil_ims replaced by new auxvar data structure
    class(pm_toil_ims_aux_type), pointer :: TOil_ims
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
  nullify(aux%Richards)
  
  nullify(aux%Mphase)
  nullify(aux%Immis)
  nullify(aux%Flash2)
  nullify(aux%Miscible)
  nullify(aux%General)
  nullify(aux%TOil_ims)
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
  call RichardsAuxDestroy(aux%Richards)
  call MphaseAuxDestroy(aux%Mphase)
  call MiscibleAuxDestroy(aux%Miscible)
  call GeneralAuxDestroy(aux%General)
  call TOilImsAuxDestroy(aux%TOil_ims) 
  !deallocate(aux%TOil_ims); nullify(aux%TOil_ims) 
  call MaterialAuxDestroy(aux%Material)
  call SecondaryAuxHeatDestroy(aux%SC_heat)
  call SecondaryAuxRTDestroy(aux%SC_RT)
  call InlineSurfaceAuxDestroy(aux%InlineSurface)
  
  nullify(aux%Global)
  nullify(aux%RT)
  nullify(aux%Richards)
  nullify(aux%Mphase)
  nullify(aux%Immis)
  nullify(aux%Miscible)
  nullify(aux%General)
  nullify(aux%TOil_ims) 
  nullify(aux%Material)
  nullify(aux%SC_Heat)
  nullify(aux%SC_RT)
  nullify(aux%InlineSurface)

end subroutine AuxDestroy

end module Auxiliary_module
