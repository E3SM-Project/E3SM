module Reaction_Surface_Complexation_Aux_module
  
  use Reaction_Database_Aux_module

  use PFLOTRAN_Constants_module

  implicit none
  
  private 

#include "petsc/finclude/petscsys.h"
  
  PetscInt, parameter, public :: SRFCMPLX_RXN_NULL = 0
  PetscInt, parameter, public :: SRFCMPLX_RXN_EQUILIBRIUM = 1
  PetscInt, parameter, public :: SRFCMPLX_RXN_MULTIRATE_KINETIC = 2
  PetscInt, parameter, public :: SRFCMPLX_RXN_KINETIC = 3

  ! surface complexation surface types
  PetscInt, parameter, public :: NULL_SURFACE = 0
  PetscInt, parameter, public :: COLLOID_SURFACE = 1
  PetscInt, parameter, public :: MINERAL_SURFACE = 2
  PetscInt, parameter, public :: ROCK_SURFACE = 3

  type, public :: surface_complex_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXWORDLENGTH) :: free_site_name
    PetscReal :: free_site_stoich
    PetscReal :: Z
    PetscReal :: forward_rate
    PetscReal :: backward_rate
    PetscBool :: print_me
    ! pointer that can be used to index the master list
    type(surface_complex_type), pointer :: ptr
    type(database_rxn_type), pointer :: dbaserxn
    type(surface_complex_type), pointer :: next
  end type surface_complex_type

  type, public :: surface_complexation_rxn_type
    PetscInt :: id
    PetscInt :: itype
    PetscInt :: free_site_id
    character(len=MAXWORDLENGTH) :: free_site_name
    PetscBool :: free_site_print_me
    PetscBool :: site_density_print_me
    PetscInt :: surface_itype
    PetscInt :: mineral_id
    character(len=MAXWORDLENGTH) :: surface_name
    PetscReal :: site_density ! site density in mol/m^3 bulk
    PetscReal, pointer :: rates(:)
    PetscReal, pointer :: site_fractions(:)
    PetscReal :: kinmr_scale_factor
    type(surface_complex_type), pointer :: complex_list
    type (surface_complexation_rxn_type), pointer :: next
  end type surface_complexation_rxn_type  
  
  type, public :: srfcplx_constraint_type
    ! Any changes here must be incorporated within ReactionProcessConstraint()
    ! where constraints are reordered
    character(len=MAXWORDLENGTH), pointer :: names(:)
    PetscReal, pointer :: constraint_conc(:)
    ! the term basis below indicates that this quantity is calculated
    ! internally in ReactionEquilibrateConstraint.
    PetscReal, pointer :: basis_free_site_conc(:)
  end type srfcplx_constraint_type


  type, public :: surface_complexation_type

    type(surface_complexation_rxn_type), pointer :: rxn_list
    type(surface_complex_type), pointer :: complex_list

    ! surface complexes
    PetscInt :: nsrfcplx
    character(len=MAXWORDLENGTH), pointer :: srfcplx_names(:)
    PetscBool, pointer :: srfcplx_print(:)
    PetscInt, pointer :: srfcplxspecid(:,:)
    PetscReal, pointer :: srfcplxstoich(:,:)
    PetscInt, pointer :: srfcplxh2oid(:)
    PetscReal, pointer :: srfcplxh2ostoich(:)
    PetscReal, pointer :: srfcplx_free_site_stoich(:)
    PetscReal, pointer :: srfcplx_logK(:)
    PetscReal, pointer :: srfcplx_logKcoef(:,:)
    PetscReal, pointer :: srfcplx_Z(:)  ! valence
    
    ! surface complexation reaction (general members)
    PetscInt :: nsrfcplxrxn
    character(len=MAXWORDLENGTH), pointer :: srfcplxrxn_site_names(:)
    PetscBool, pointer :: srfcplxrxn_site_print(:)
    PetscBool, pointer :: srfcplxrxn_site_density_print(:)
    PetscInt, pointer :: srfcplxrxn_to_surf(:)
    PetscInt, pointer :: srfcplxrxn_surf_type(:)
    PetscInt, pointer :: srfcplxrxn_to_complex(:,:)
    ! site density in 
    ! (1) mol/m^3 bulk 
    ! (2) mol/m^3 mineral, which * mineral volume fraction = mol/m^3 bulk
    ! (3) mol/kg rock, which * 1-porosity = mol/m^3 bulk
    PetscReal, pointer :: srfcplxrxn_site_density(:) 
    ! this flag indicates that the stoichiometry for free sites in one of the
    ! reactions is not 1, and thus we must use nonlinear iteration to solve
    PetscBool, pointer :: srfcplxrxn_stoich_flag(:)

    ! equilibrium
    PetscInt :: neqsrfcplx
    PetscInt :: neqsrfcplxrxn 
    PetscInt, pointer :: eqsrfcplxrxn_to_srfcplxrxn(:)
!geh: will not use for now.  allocate eqsrfcplx_conc to size of total #
!     of surface complexes
    !PetscInt, pointer :: srfcplx_to_eqsrfcplx(:)
    
    ! kinetic
    PetscInt :: nkinsrfcplx
    PetscInt :: nkinsrfcplxrxn
    PetscInt, pointer :: kinsrfcplxrxn_to_srfcplxrxn(:)
    PetscInt, pointer :: kinsrfcplx_to_name(:,:)
    PetscReal, pointer :: kinsrfcplx_forward_rate(:,:)
    PetscReal, pointer :: kinsrfcplx_backward_rate(:,:)  

    ! multirate kinetic surface complexation
    PetscInt :: nkinmrsrfcplx
    PetscInt :: nkinmrsrfcplxrxn 
    PetscInt, pointer :: kinmrsrfcplxrxn_to_srfcplxrxn(:)
    ! the zeroth entry of kinmr_nrate will hold the max nrate
    PetscInt, pointer :: kinmr_nrate(:)
    PetscReal, pointer :: kinmr_rate(:,:)
    PetscReal, pointer :: kinmr_frac(:,:)
    
  end type surface_complexation_type

  public :: SurfaceComplexationCreate, &
            SurfaceComplexationRxnCreate, &
            SurfaceComplexCreate, &
            SurfaceComplexConstraintCreate, &
            SrfCplxGetSrfCplxCountInRxnType, &
            SrfCplxMapMasterSrfCplxToRxn, &
            SurfaceComplexationDestroy, &
            SurfaceComplexationRxnDestroy, &
            SurfaceComplexDestroy, &
            SurfaceComplexConstraintDestroy
             
contains

! ************************************************************************** !

function SurfaceComplexationCreate()
  ! 
  ! Allocate and initialize surface complexation
  ! object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  use Option_module

  implicit none
  
  type(surface_complexation_type), pointer :: SurfaceComplexationCreate
  
  type(surface_complexation_type), pointer :: surface_complexation

  allocate(surface_complexation)  

  nullify(surface_complexation%rxn_list)
  nullify(surface_complexation%complex_list)

  ! surface complexes
  surface_complexation%nsrfcplx = 0
  
  nullify(surface_complexation%srfcplx_names)
  nullify(surface_complexation%srfcplx_print)

  nullify(surface_complexation%srfcplxspecid)
  nullify(surface_complexation%srfcplxstoich)
  nullify(surface_complexation%srfcplxh2oid)
  nullify(surface_complexation%srfcplxh2ostoich)
  nullify(surface_complexation%srfcplx_free_site_stoich)
  nullify(surface_complexation%srfcplx_logK)
  nullify(surface_complexation%srfcplx_logKcoef)
  nullify(surface_complexation%srfcplx_Z)

  ! surface complexation reaction (general members)
  surface_complexation%nsrfcplxrxn = 0
  nullify(surface_complexation%srfcplxrxn_site_names)
  nullify(surface_complexation%srfcplxrxn_site_print)
  nullify(surface_complexation%srfcplxrxn_site_density_print)

  nullify(surface_complexation%srfcplxrxn_to_surf)
  nullify(surface_complexation%srfcplxrxn_surf_type)
  nullify(surface_complexation%srfcplxrxn_to_complex)
  nullify(surface_complexation%srfcplxrxn_site_density)
  nullify(surface_complexation%srfcplxrxn_stoich_flag) 
  
  ! equilibrium
  surface_complexation%neqsrfcplx = 0
  surface_complexation%neqsrfcplxrxn = 0
  nullify(surface_complexation%eqsrfcplxrxn_to_srfcplxrxn) 
  
  ! kinetic   
  surface_complexation%nkinsrfcplx = 0
  surface_complexation%nkinsrfcplxrxn = 0
  nullify(surface_complexation%kinsrfcplxrxn_to_srfcplxrxn) 
  nullify(surface_complexation%kinsrfcplx_to_name) 
  nullify(surface_complexation%kinsrfcplx_forward_rate)
  nullify(surface_complexation%kinsrfcplx_backward_rate)

  ! multirate kinetic surface complexation
  surface_complexation%nkinmrsrfcplx = 0
  surface_complexation%nkinmrsrfcplxrxn = 0
  nullify(surface_complexation%kinmrsrfcplxrxn_to_srfcplxrxn) 
  nullify(surface_complexation%kinmr_nrate)
  nullify(surface_complexation%kinmr_rate)
  nullify(surface_complexation%kinmr_frac)

  SurfaceComplexationCreate => surface_complexation
  
end function SurfaceComplexationCreate

! ************************************************************************** !

function SurfaceComplexationRxnCreate()
  ! 
  ! Allocate and initialize a surface complexation
  ! reaction
  ! 
  ! Author: Peter Lichtner
  ! Date: 10/21/08
  ! 

  implicit none
    
  type(surface_complexation_rxn_type), pointer :: SurfaceComplexationRxnCreate

  type(surface_complexation_rxn_type), pointer :: srfcplxrxn
  
  allocate(srfcplxrxn)
  srfcplxrxn%free_site_id = 0
  srfcplxrxn%itype = SRFCMPLX_RXN_NULL
  srfcplxrxn%free_site_name = ''
  srfcplxrxn%free_site_print_me = PETSC_FALSE
  srfcplxrxn%site_density_print_me = PETSC_FALSE

  srfcplxrxn%surface_itype = NULL_SURFACE
  srfcplxrxn%mineral_id = 0
  srfcplxrxn%surface_name = ''
  srfcplxrxn%site_density = 0.d0
  srfcplxrxn%kinmr_scale_factor = 1.d0
  nullify(srfcplxrxn%rates)
  nullify(srfcplxrxn%site_fractions)
  
  nullify(srfcplxrxn%complex_list)
  nullify(srfcplxrxn%next)
  
  SurfaceComplexationRxnCreate => srfcplxrxn
  
end function SurfaceComplexationRxnCreate

! ************************************************************************** !

function SurfaceComplexCreate()
  ! 
  ! Allocate and initialize a surface complex reaction
  ! 
  ! Author: Peter Lichtner
  ! Date: 10/21/08
  ! 

  implicit none
    
  type(surface_complex_type), pointer :: SurfaceComplexCreate

  type(surface_complex_type), pointer :: srfcplx
  
  allocate(srfcplx)
  srfcplx%id = 0
  srfcplx%name = ''
  srfcplx%free_site_name = ''
  srfcplx%Z = 0.d0
  srfcplx%free_site_stoich = 0.d0
  srfcplx%forward_rate = 0.d0
  ! default is UNINITIALIZED_INTEGER in case the only the forward rate is defined.  In that case
  ! the backward rate will be calculated as a function of the forward rate and
  ! the equilibrium coefficient (logK).
  srfcplx%backward_rate = UNINITIALIZED_DOUBLE
  srfcplx%print_me = PETSC_FALSE
  nullify(srfcplx%ptr)
  nullify(srfcplx%dbaserxn)
  nullify(srfcplx%next)
  
  SurfaceComplexCreate => srfcplx
  
end function SurfaceComplexCreate

! ************************************************************************** !

function SurfaceComplexConstraintCreate(surface_complexation,option)
  ! 
  ! Creates a surface complex constraint object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/09
  ! 

  use Option_module
  
  implicit none
  
  type(surface_complexation_type) :: surface_complexation
  type(option_type) :: option
  type(srfcplx_constraint_type), pointer :: SurfaceComplexConstraintCreate

  type(srfcplx_constraint_type), pointer :: constraint  

  allocate(constraint)
  allocate(constraint%names(surface_complexation%nkinsrfcplx))
  constraint%names = ''
  allocate(constraint%constraint_conc(surface_complexation%nkinsrfcplx))
  constraint%constraint_conc = 0.d0
  allocate(constraint%basis_free_site_conc( &
                           surface_complexation%nkinsrfcplxrxn))
  constraint%basis_free_site_conc = 0.d0

  SurfaceComplexConstraintCreate => constraint

end function SurfaceComplexConstraintCreate

! ************************************************************************** !

function SrfCplxGetSrfCplxCountInRxnType(surface_complexation,rxn_type)
  ! 
  ! Deallocates a surface complexation reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/22/12
  ! 

  implicit none
  
  type(surface_complexation_type) :: surface_complexation
  PetscInt :: rxn_type
  
  type(surface_complexation_rxn_type), pointer :: cur_srfcplx_rxn
  type(surface_complex_type), pointer :: cur_srfcplx
  
  PetscInt :: SrfCplxGetSrfCplxCountInRxnType
  
  SrfCplxGetSrfCplxCountInRxnType = 0

  ! to determine the number of unique equilibrium surface complexes,
  ! we negate the ids of the complexes in the master list as a flag
  ! to avoid duplicate counts.  We then traverse the list of 
  ! complexes in the rxn and see how many ids are negated
  cur_srfcplx => surface_complexation%complex_list
  do
    if (.not.associated(cur_srfcplx)) exit
    cur_srfcplx%id = -abs(cur_srfcplx%id)
    cur_srfcplx => cur_srfcplx%next
  enddo   
  cur_srfcplx_rxn => surface_complexation%rxn_list
  do
    if (.not.associated(cur_srfcplx_rxn)) exit
    if (cur_srfcplx_rxn%itype == rxn_type) then
      cur_srfcplx => cur_srfcplx_rxn%complex_list
      do
        if (.not.associated(cur_srfcplx)) exit
        ! recall that complexes in rxns point to complexes in master list
        ! through a ptr member of the derived type
        if (cur_srfcplx%ptr%id < 0) then
          cur_srfcplx%ptr%id = abs(cur_srfcplx%ptr%id)
          SrfCplxGetSrfCplxCountInRxnType = SrfCplxGetSrfCplxCountInRxnType + 1
        endif
        cur_srfcplx => cur_srfcplx%next
      enddo
    endif
    cur_srfcplx_rxn => cur_srfcplx_rxn%next
  enddo
  ! unflag complexes (see comment above)
  cur_srfcplx => surface_complexation%complex_list
  do
    if (.not.associated(cur_srfcplx)) exit
    cur_srfcplx%id = abs(cur_srfcplx%id)
    cur_srfcplx => cur_srfcplx%next
  enddo   

end function SrfCplxGetSrfCplxCountInRxnType

! ************************************************************************** !

subroutine SrfCplxMapMasterSrfCplxToRxn(surface_complexation,rxn_type)
  ! 
  ! Maps surface complexes from the master list
  ! to a compressed rxn list (e.g. for array)
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/22/12
  ! 

  implicit none
  
  type(surface_complexation_type) :: surface_complexation
  PetscInt :: rxn_type
  
  type(surface_complexation_rxn_type), pointer :: cur_srfcplx_rxn
  type(surface_complex_type), pointer :: cur_srfcplx
  PetscInt, allocatable :: srfcplx_to_rxnsrfcplx(:)
  PetscInt :: isrfcplx, isrfcplx_in_rxn
  
  isrfcplx = SrfCplxGetSrfCplxCountInRxnType(surface_complexation,rxn_type)
  allocate(srfcplx_to_rxnsrfcplx(isrfcplx))

  ! flag the master list
  cur_srfcplx => surface_complexation%complex_list
  isrfcplx = 0
  do
    if (.not.associated(cur_srfcplx)) exit
    ! don't need the ptr here since it is the master list
    isrfcplx = isrfcplx + 1
    cur_srfcplx%id = -abs(cur_srfcplx%id)
    cur_srfcplx => cur_srfcplx%next
  enddo
  allocate(srfcplx_to_rxnsrfcplx(isrfcplx))
  srfcplx_to_rxnsrfcplx = 0
  
  ! determine which surface complexes in master list are include in 
  ! reaction by flagging them
  cur_srfcplx_rxn => surface_complexation%rxn_list
  do
    if (.not.associated(cur_srfcplx_rxn)) exit
    if (cur_srfcplx_rxn%itype == rxn_type) then
      cur_srfcplx => cur_srfcplx_rxn%complex_list
      do
        if (.not.associated(cur_srfcplx)) exit
        cur_srfcplx%ptr%id = abs(cur_srfcplx%ptr%id)
        cur_srfcplx => cur_srfcplx%next
      enddo
    endif
    cur_srfcplx_rxn => cur_srfcplx_rxn%next
  enddo
  
  ! load flagged names in master list into the list of names in the
  ! same order as master list
  ! flag the master list
  cur_srfcplx => surface_complexation%complex_list
  isrfcplx = 0
  isrfcplx_in_rxn = 0
  do
    if (.not.associated(cur_srfcplx)) exit
    ! don't need the ptr here since it is the master list
    isrfcplx = isrfcplx + 1
    if (cur_srfcplx%id > 0) then
      isrfcplx_in_rxn = isrfcplx_in_rxn + 1
      srfcplx_to_rxnsrfcplx(isrfcplx) = isrfcplx_in_rxn
    endif
    ! clear all flags
    cur_srfcplx%ptr%id = abs(cur_srfcplx%ptr%id)
    cur_srfcplx => cur_srfcplx%next
  enddo   
  
  select case(rxn_type)
    case(SRFCMPLX_RXN_EQUILIBRIUM)
!      if (.not.associated(surface_complexation%srfcplx_to_eqsrfcplx)) then
!        allocate(surface_complexation% &
!                   srfcplx_to_eqsrfcplx(size(srfcplx_to_rxnsrfcplx)))
!      endif
!      surface_complexation%srfcplx_to_eqsrfcplx = srfcplx_to_rxnsrfcplx
    case(SRFCMPLX_RXN_MULTIRATE_KINETIC)
    case(SRFCMPLX_RXN_KINETIC)
  end select
  
  deallocate(srfcplx_to_rxnsrfcplx)
  
end subroutine SrfCplxMapMasterSrfCplxToRxn

! ************************************************************************** !

subroutine SurfaceComplexationRxnDestroy(srfcplxrxn)
  ! 
  ! Deallocates a surface complexation reaction
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/08
  ! 

  use Utility_module

  implicit none
    
  type(surface_complexation_rxn_type), pointer :: srfcplxrxn

  type(surface_complex_type), pointer :: cur_srfcplx, prev_srfcplx
  
  if (.not.associated(srfcplxrxn)) return
  
  cur_srfcplx => srfcplxrxn%complex_list
  do
    if (.not.associated(cur_srfcplx)) exit
    prev_srfcplx => cur_srfcplx
    cur_srfcplx => cur_srfcplx%next
    call SurfaceComplexDestroy(prev_srfcplx)
    nullify(prev_srfcplx)
  enddo
  
  call DeallocateArray(srfcplxrxn%rates)
  call DeallocateArray(srfcplxrxn%site_fractions)
  
  deallocate(srfcplxrxn)  
  nullify(srfcplxrxn)

end subroutine SurfaceComplexationRxnDestroy

! ************************************************************************** !

subroutine SurfaceComplexDestroy(srfcplx)
  ! 
  ! Deallocates a surface complex
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/21/08
  ! 

  implicit none
    
  type(surface_complex_type), pointer :: srfcplx

  if (.not.associated(srfcplx)) return
  
  if (associated(srfcplx%dbaserxn)) &
    call DatabaseRxnDestroy(srfcplx%dbaserxn)
  nullify(srfcplx%dbaserxn)
  nullify(srfcplx%next)

  deallocate(srfcplx)  
  nullify(srfcplx)

end subroutine SurfaceComplexDestroy

! ************************************************************************** !

subroutine SurfaceComplexConstraintDestroy(constraint)
  ! 
  ! Destroys a surface complex constraint
  ! object
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/21/09
  ! 

  use Utility_module, only : DeallocateArray
  implicit none
  
  type(srfcplx_constraint_type), pointer :: constraint
  
  if (.not.associated(constraint)) return
  
  call DeallocateArray(constraint%names)
  call DeallocateArray(constraint%constraint_conc)
  call DeallocateArray(constraint%basis_free_site_conc)
  
  deallocate(constraint)
  nullify(constraint)

end subroutine SurfaceComplexConstraintDestroy

! ************************************************************************** !

subroutine SurfaceComplexationDestroy(surface_complexation)
  ! 
  ! Deallocates a reaction object
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/12
  ! 

  use Utility_module, only: DeallocateArray 
  
  implicit none

  type(surface_complexation_type), pointer :: surface_complexation
  
  type(surface_complexation_rxn_type), pointer :: cur_srfcplxrxn, prev_srfcplxrxn
  type(surface_complex_type), pointer :: cur_srfcplx, prev_srfcplx

  if (.not.associated(surface_complexation)) return
  
  ! surface complexation reactions
  cur_srfcplxrxn => surface_complexation%rxn_list
  do
    if (.not.associated(cur_srfcplxrxn)) exit
    prev_srfcplxrxn => cur_srfcplxrxn
    cur_srfcplxrxn => cur_srfcplxrxn%next
    call SurfaceComplexationRxnDestroy(prev_srfcplxrxn)
  enddo    
  nullify(surface_complexation%rxn_list)
  
  ! surface complexes
  cur_srfcplx => surface_complexation%complex_list
  do
    if (.not.associated(cur_srfcplx)) exit
    prev_srfcplx => cur_srfcplx
    cur_srfcplx => cur_srfcplx%next
    call SurfaceComplexDestroy(prev_srfcplx)
    nullify(prev_srfcplx)
  enddo

  ! surface complexes
  call DeallocateArray(surface_complexation%srfcplx_names)
  call DeallocateArray(surface_complexation%srfcplx_print)
  call DeallocateArray(surface_complexation%srfcplxspecid)
  call DeallocateArray(surface_complexation%srfcplxstoich)
  call DeallocateArray(surface_complexation%srfcplxh2oid)
  call DeallocateArray(surface_complexation%srfcplxh2ostoich)
  call DeallocateArray(surface_complexation%srfcplx_free_site_stoich)
  call DeallocateArray(surface_complexation%srfcplx_logK)
  call DeallocateArray(surface_complexation%srfcplx_logKcoef)
  call DeallocateArray(surface_complexation%srfcplx_Z)
    
  ! surface complexation reaction (general members)
  call DeallocateArray(surface_complexation%srfcplxrxn_site_names)
  call DeallocateArray(surface_complexation%srfcplxrxn_site_print)
  call DeallocateArray(surface_complexation%srfcplxrxn_site_density_print)
  call DeallocateArray(surface_complexation%srfcplxrxn_to_surf)
  call DeallocateArray(surface_complexation%srfcplxrxn_surf_type)
  call DeallocateArray(surface_complexation%srfcplxrxn_to_complex)
  call DeallocateArray(surface_complexation%srfcplxrxn_site_density)
  call DeallocateArray(surface_complexation%srfcplxrxn_stoich_flag)

  ! equilibrium
  call DeallocateArray(surface_complexation%eqsrfcplxrxn_to_srfcplxrxn)
    
  ! kinetic
  call DeallocateArray(surface_complexation%kinsrfcplxrxn_to_srfcplxrxn)
  call DeallocateArray(surface_complexation%kinsrfcplx_to_name)
  call DeallocateArray(surface_complexation%kinsrfcplx_forward_rate)
  call DeallocateArray(surface_complexation%kinsrfcplx_backward_rate)

  ! multirate kinetic surface complexation
  call DeallocateArray(surface_complexation%kinmrsrfcplxrxn_to_srfcplxrxn)
  call DeallocateArray(surface_complexation%kinmr_nrate)
  call DeallocateArray(surface_complexation%kinmr_rate)
  call DeallocateArray(surface_complexation%kinmr_frac)
  
  deallocate(surface_complexation)
  nullify(surface_complexation)

end subroutine SurfaceComplexationDestroy

end module Reaction_Surface_Complexation_Aux_module
