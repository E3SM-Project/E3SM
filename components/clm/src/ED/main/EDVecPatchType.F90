module EDVecPatchType

  !----------------------------------------------------
  ! define the ED pft structure
  !----------------------------------------------------

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type 

  type, public :: EDVecPatch_type

     integer , pointer :: ED_patch(:)     ! Is there a patch of vegetation defined for this p?  
     integer , pointer :: ED_bareground(:)! Is this patch the designated bare ground fraction?  
     real(r8), pointer :: wtED(:)         ! What is the weight of each patch in the ED natural vegetation area? 

   contains

     procedure, public :: Init
     
  end type EDVecPatch_type

  type(EDVecPatch_type), target :: EDpft  ! Vector ED patch data structure
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !ARGUMENTS:
    class(EDVecPatch_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! LOCAL VARAIBLES:
    integer  :: begp, endp
    !------------------------------------------------------------------------

    begp = bounds%begp; endp=bounds%endp

    allocate(this%ED_patch      (begp:endp)) 
    allocate(this%ED_bareground (begp:endp)) 
    allocate(this%wtED          (begp:endp)) 

  end subroutine Init

end module EDVecPatchType
