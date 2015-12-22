#ifdef USE_PETSC_LIB


module ConditionType

  ! !USES:
  use ConnectionSetType  , only : connection_set_type
  use clm_varctl         , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"

  type, public :: condition_type

     character (len=256)                :: name                   ! name for condition
     character (len=256)                :: units                  ! units

     PetscInt                           :: id                     ! identifier of condition within the list
     PetscInt                           :: itype                  ! identifier for type of condition
     PetscInt                           :: region_itype           ! identifier for region
     PetscInt                           :: ncells
     PetscReal, pointer                 :: value(:)               ! Magnitude of the condition

     PetscInt                           :: list_id_of_other_goveq ! List ID of other governing equation
     PetscBool                          :: swap_order             ! FALSE(default): "upwind cell  " = BC; "downwind cell" = Internal cell
                                                                  ! TRUE          : "downwind cell" = BC; "upwind cell  " = Internal cell

     type(connection_set_type), pointer :: conn_set               ! Applicable to BC condition type
     type(condition_type), pointer      :: next                   ! Pointer to next condition

  end type condition_type

  type, public :: condition_list_type
     PetscInt                         :: num_condition_list
     type(condition_type), pointer    :: first
     type(condition_type), pointer    :: last
  end type condition_list_type

  public :: ConditionNew
  public :: ConditionDestroy
  public :: ConditionPrintInfo
  public :: ConditionListInit
  public :: ConditionListAddCondition
  public :: ConditionListClean
  public :: ConditionListPrintInfo
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  function ConditionNew()
    !
    ! !DESCRIPTION:
    ! Return a new condition
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_type),pointer :: ConditionNew
    type(condition_type),pointer :: cond

    allocate(cond)

    cond%name                   = ""
    cond%units                  = ""
    cond%id                     = -1
    cond%itype                  = -1
    cond%region_itype           = -1
    cond%ncells                 = 0

    cond%list_id_of_other_goveq = -1
    cond%swap_order             = PETSC_FALSE

    nullify(cond%value    )
    nullify(cond%conn_set )
    nullify(cond%next     )

    ConditionNew => cond

  end function ConditionNew

  !------------------------------------------------------------------------
  subroutine ConditionListInit(list)
    !
    ! !DESCRIPTION:
    ! Initialize an existing condition list
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_list_type) :: list

    list%num_condition_list   = 0

    nullify(list%first)
    nullify(list%last )

  end subroutine ConditionListInit

  !------------------------------------------------------------------------
  subroutine ConditionPrintInfo(cur_cond)
    !
    ! !DESCRIPTION:
    ! Prints information about the condition
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_type) :: cur_cond

    write(iulog,*) '      Condition_name : ',trim(cur_cond%name)
    write(iulog,*) '      Condition_units: ',trim(cur_cond%units)

  end subroutine ConditionPrintInfo

  !------------------------------------------------------------------------
  subroutine ConditionListPrintInfo(list)
    !
    ! !DESCRIPTION:
    ! Prints information about all conditions present in the condition list
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_list_type), intent(in) :: list
    !
    type(condition_type), pointer         :: cur_cond

    cur_cond => list%first
    if (.not.associated(cur_cond)) return

    write(iulog,*) '    Condition-List_num : ',list%num_condition_list

    do
       if (.not.associated(cur_cond)) exit
       call ConditionPrintInfo(cur_cond)
       write(iulog,*) '      Number_of_Conn : ',cur_cond%conn_set%num_connections
       write(iulog,*)''
       cur_cond => cur_cond%next
    enddo
    write(iulog,*)''

  end subroutine ConditionListPrintInfo


  !------------------------------------------------------------------------
  subroutine ConditionListAddCondition(list, new_cond)
    !
    ! !DESCRIPTION:
    ! Add a condition to a condition list
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_list_type)    :: list
    type(condition_type),pointer :: new_cond

    list%num_condition_list = list%num_condition_list + 1
    new_cond%id             = list%num_condition_list

    if (.not.associated(list%first)) then
       list%first => new_cond
    endif

    if (associated(list%last)) then
       list%last%next => new_cond
    endif

    list%last => new_cond

  end subroutine ConditionListAddCondition


  !------------------------------------------------------------------------
  subroutine ConditionListClean(list)
    !
    ! !DESCRIPTION:
    ! Release all allocated memory
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_list_type), intent(inout) :: list
    !
    ! !LOCAL VARIABLES:
    type(condition_type), pointer            :: cur_cond, prev_cond
    
    cur_cond => list%first
    do
       if (.not.associated(cur_cond)) exit
       prev_cond => cur_cond
       cur_cond => cur_cond%next
       call ConditionDestroy(prev_cond)
    enddo

    call ConditionListInit(list)

  end subroutine ConditionListClean


  !------------------------------------------------------------------------
  subroutine ConditionDestroy(cond)
    !
    ! !DESCRIPTION:
    ! Release all allocated memory
    !
    ! !USES:
    use ConnectionSetType  , only : ConnectionSetDestroy
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_type) :: cond

    cond%name           = ""
    cond%units          = ""
    cond%id             = -1
    cond%itype          = -1
    cond%region_itype   = -1
    cond%ncells         = 0

    if (associated(cond%value)) deallocate(cond%value)
    nullify(cond%value)

    call ConnectionSetDestroy(cond%conn_set)
    nullify(cond%conn_set)

  end subroutine ConditionDestroy

end module ConditionType
#endif
