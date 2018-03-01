module ConditionType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use ConnectionSetType  , only : connection_set_type
  use mpp_varctl         , only : iulog
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: condition_type

     character (len=256)                :: name                       ! name for condition
     character (len=256)                :: units                      ! units

     PetscInt                           :: id                         ! identifier of condition within the list
     PetscInt                           :: itype                      ! identifier for type of condition
     PetscInt                           :: region_itype               ! identifier for region
     PetscInt                           :: ncells
     PetscReal, pointer                 :: value(:)                   ! Magnitude of the condition

     PetscInt                           :: num_other_goveqs           ! Number of other governing equations
     PetscInt, pointer                  :: list_id_of_other_goveqs(:) ! List ID of other governing equations
     PetscInt, pointer                  :: itype_of_other_goveqs(:)   ! Type of other governing equations (e.g. GE_THERM_SSW_TBASED, GE_THERM_SNOW_TBASED, etc)
     PetscBool                          :: swap_order                 ! FALSE(default): "upwind cell  " = BC; "downwind cell" = Internal cell
                                                                      ! TRUE          : "downwind cell" = BC; "upwind cell  " = Internal cell
     PetscBool, pointer                 :: swap_order_of_other_goveqs(:)
     PetscBool, pointer                 :: coupled_via_intauxvar_with_other_goveqns(:)

     type(connection_set_type), pointer :: conn_set                   ! Applicable to BC condition type
     type(condition_type), pointer      :: next                       ! Pointer to next condition

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
  public :: CondListGetNumCondsExcptCondItype
  public :: CondListGetNumCellsForCondsExcptCondItype
  public :: CondListGetConnIDDnForCondsExcptCondItype
  public :: CondListGetCondNamesExcptCondItype
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

    cond%num_other_goveqs       = 0
    cond%swap_order             = PETSC_FALSE

    nullify(cond%value                                   )
    nullify(cond%list_id_of_other_goveqs                 )
    nullify(cond%itype_of_other_goveqs                   )
    nullify(cond%conn_set                                )
    nullify(cond%next                                    )
    nullify(cond%swap_order_of_other_goveqs              )
    nullify(cond%coupled_via_intauxvar_with_other_goveqns)

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
  subroutine CondListGetNumCondsExcptCondItype( &
       list, cond_itype, num_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of conditions excluding conditions of
    ! itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_list_type) , intent(in)  :: list
    PetscInt                  , intent(in)  :: cond_itype
    PetscInt                  , intent(out) :: num_conds
    !
    ! !LOCAL VARIABLES:
    type(condition_type), pointer            :: cur_cond

    num_conds = 0
    cur_cond => list%first
    do
       if (.not.associated(cur_cond)) exit
       if (cur_cond%itype /= cond_itype) then
          num_conds = num_conds + 1
       endif
       cur_cond => cur_cond%next
    enddo

  end subroutine CondListGetNumCondsExcptCondItype

  !------------------------------------------------------------------------
  subroutine CondListGetNumCellsForCondsExcptCondItype( &
       list, cond_itype, num_cells_for_conds)
    !
    ! !DESCRIPTION:
    ! Returns the total number of cells associated with all conditions excluding
    ! conditions of itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_list_type) , intent(in)           :: list
    PetscInt                  , intent(in)           :: cond_itype
    PetscInt                  , intent(out), pointer :: num_cells_for_conds(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                         :: num_conds
    type(condition_type)      , pointer              :: cur_cond

    call CondListGetNumCondsExcptCondItype(list, cond_itype, num_conds)

    if (num_conds == 0) then

       nullify(num_cells_for_conds)

    else

       allocate(num_cells_for_conds(num_conds))

       num_conds = 0
       cur_cond => list%first
       do
          if (.not.associated(cur_cond)) exit
          if (cur_cond%itype /= cond_itype) then
             num_conds = num_conds + 1
             num_cells_for_conds(num_conds) = cur_cond%ncells
          endif
          cur_cond => cur_cond%next
       enddo

    end if

  end subroutine CondListGetNumCellsForCondsExcptCondItype

  !------------------------------------------------------------------------
  subroutine CondListGetConnIDDnForCondsExcptCondItype( &
       list, cond_itype, num_cells, cell_id_dn)
    !
    ! !DESCRIPTION:
    ! Returns the downwind cell associated with all conditions excluding
    ! conditions of itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_list_type) , intent(in)           :: list
    PetscInt                  , intent(in)           :: cond_itype
    PetscInt                  , intent(out)          :: num_cells
    PetscInt                  , intent(out), pointer :: cell_id_dn(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                                         :: icond
    PetscInt                                         :: iconn
    PetscInt                                         :: count
    PetscInt                                         :: num_conds
    PetscInt                  , pointer              :: num_cells_for_conds(:)
    type(condition_type)      , pointer              :: cur_cond
    type(connection_set_type) , pointer              :: cur_conn_set

    call CondListGetNumCondsExcptCondItype(list, cond_itype, num_conds)
    call CondListGetNumCellsForCondsExcptCondItype(list, cond_itype, num_cells_for_conds)

    num_cells = 0
    do icond = 1, num_conds
       num_cells = num_cells + num_cells_for_conds(icond)
    end do

    if (num_cells == 0) then

       nullify(cell_id_dn)

    else

       allocate(cell_id_dn(num_cells))

       count = 0
       cur_cond => list%first
       do
          if (.not.associated(cur_cond)) exit

          if (cur_cond%itype /= cond_itype) then
             cur_conn_set => cur_cond%conn_set
             do iconn = 1, cur_conn_set%num_connections
                count = count + 1
                cell_id_dn(count) = cur_conn_set%conn(iconn)%GetIDDn()
             end do
          endif
          cur_cond => cur_cond%next
       enddo

    end if

  end subroutine CondListGetConnIDDnForCondsExcptCondItype

  !------------------------------------------------------------------------
  subroutine CondListGetCondNamesExcptCondItype( &
       list, cond_itype, cond_names)
    !
    ! !DESCRIPTION:
    ! Returns the name of all conditions except condition of
    ! itype = cond_type_to_exclude
    !
    implicit none
    !
    ! !ARGUMENTS
    type(condition_list_type) , intent(in) :: list
    PetscInt                  , intent(in) :: cond_itype
    character (len=256)       , pointer    :: cond_names(:)
    !
    ! !LOCAL VARIABLES:
    PetscInt                               :: num_conds
    type(condition_type)     , pointer     :: cur_cond

    call CondListGetNumCondsExcptCondItype( &
         list, cond_itype, num_conds)

    if (num_conds == 0) then

       nullify(cond_names)

    else

       allocate(cond_names(num_conds))

       num_conds = 0
       cur_cond => list%first
       do
          if (.not.associated(cur_cond)) exit
          if (cur_cond%itype /= cond_itype) then
             num_conds = num_conds + 1
             cond_names(num_conds) = trim(cur_cond%name)
          endif
          cur_cond => cur_cond%next
       enddo

    end if

  end subroutine CondListGetCondNamesExcptCondItype

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

#endif

end module ConditionType
