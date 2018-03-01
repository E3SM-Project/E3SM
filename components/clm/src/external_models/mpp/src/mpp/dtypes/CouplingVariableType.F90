module CouplingVariableType

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use mpp_varctl         , only : iulog
  use petscsys
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: coupling_variable_type
     PetscInt                      :: variable_type                     ! id of the coupling variable (e.g. VAR_PRESSURE)
     PetscInt                      :: num_cells                         ! number of cells associated with the coupling variable
     PetscInt                      :: rank_of_coupling_goveqn           ! rank of the coupling governing equation (GE) in the SoE list
     PetscBool                     :: variable_is_bc_in_coupling_goveqn ! variable in the coupling goveqn is a boundary (=TRUE)
     PetscInt                      :: offset_of_bc_in_current_goveqn    ! offset of the boundary auxvars in the current GE
     PetscInt                      :: rank_of_bc_in_current_goveqn      ! rank of the boundary condition within the boundary
                                                                        ! condition list of the current GE
     PetscInt                      :: rank_of_bc_in_coupling_goveqn     ! rank of the boundary condition within the boundary
                                                                        ! condition list of the coupling GE provided
                                                                        ! variable_is_bc_in_coupling_goveqn = TRUE 
     type (coupling_variable_type), pointer :: next                              ! pointer to the next coupling variable
  end type coupling_variable_type

  type, public :: coupling_variable_list_type
     PetscInt                              :: num_coupling_vars
     type(coupling_variable_type), pointer :: first
     type(coupling_variable_type), pointer :: last
  end type coupling_variable_list_type

  public :: CouplingVariableCreate
  public :: CouplingVariablePrintInfo
  public :: CouplingVariableDestroy
  public :: CouplingVariableListCreate
  public :: CouplingVariableListDestroy
  public :: CouplingVariableListAddCouplingVar
  public :: CouplingVariableListPrintInfo

  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  function CouplingVariableCreate()
    !
    ! !DESCRIPTION:
    ! Return a new coupling var
    !
    implicit none
    !
    ! !ARGUMENTS
    type(coupling_variable_type), pointer :: CouplingVariableCreate
    type(coupling_variable_type), pointer :: cpl_var

    allocate(cpl_var)

    cpl_var%variable_type                     = 0
    cpl_var%num_cells                         = 0
    cpl_var%rank_of_coupling_goveqn           = 0
    cpl_var%variable_is_bc_in_coupling_goveqn = PETSC_FALSE
    cpl_var%offset_of_bc_in_current_goveqn    = 0
    cpl_var%rank_of_bc_in_current_goveqn      = 0
    cpl_var%rank_of_bc_in_coupling_goveqn     = 0

    nullify(cpl_var%next)

    CouplingVariableCreate => cpl_var

  end function CouplingVariableCreate

  !------------------------------------------------------------------------
  subroutine CouplingVariablePrintInfo(cpl_var)
    !
    ! !DESCRIPTION:
    ! Prints information about the coupler variable
    !
    implicit none
    !
    ! !ARGUMENTS
    type(coupling_variable_type) :: cpl_var

    write(iulog,*) '      Coupler_variable_type: ',cpl_var%variable_type

  end subroutine CouplingVariablePrintInfo

  !------------------------------------------------------------------------
  subroutine CouplingVariableDestroy(cpl_var)
    !
    ! !DESCRIPTION:
    ! Resets coupling variable
    !
    implicit none
    !
    ! !ARGUMENTS
    type(coupling_variable_type) :: cpl_var

    cpl_var%variable_type                     = 0
    cpl_var%num_cells                         = 0
    cpl_var%rank_of_coupling_goveqn           = 0
    cpl_var%variable_is_bc_in_coupling_goveqn = PETSC_FALSE
    cpl_var%offset_of_bc_in_current_goveqn    = 0
    cpl_var%rank_of_bc_in_current_goveqn      = 0
    cpl_var%rank_of_bc_in_coupling_goveqn     = 0

  end subroutine CouplingVariableDestroy

  !------------------------------------------------------------------------
  subroutine CouplingVariableListCreate(list)
    !
    ! !DESCRIPTION:
    ! Initialize an existing condition list
    !
    implicit none
    !
    ! !ARGUMENTS
    type(coupling_variable_list_type) :: list

    list%num_coupling_vars   = 0

    nullify(list%first)
    nullify(list%last )

  end subroutine CouplingVariableListCreate

  !------------------------------------------------------------------------
  subroutine CouplingVariableListAddCouplingVar(list, new_cpl_var)
    !
    ! !DESCRIPTION:
    ! Add a condition to a condition list
    !
    implicit none
    !
    ! !ARGUMENTS
    type(coupling_variable_list_type)     :: list
    type(coupling_variable_type), pointer :: new_cpl_var

    list%num_coupling_vars = list%num_coupling_vars + 1

    if (.not.associated(list%first)) then
       list%first => new_cpl_var
    endif

    if (associated(list%last)) then
       list%last%next => new_cpl_var
    endif

    list%last => new_cpl_var

  end subroutine CouplingVariableListAddCouplingVar

  !------------------------------------------------------------------------
  subroutine CouplingVariableListPrintInfo(list)
    !
    ! !DESCRIPTION:
    ! Prints information about all coupling variables present in the list
    !
    implicit none
    !
    ! !ARGUMENTS
    type(coupling_variable_list_type), intent(in) :: list
    !
    type(coupling_variable_type), pointer         :: cur_cpl_var

    cur_cpl_var => list%first
    if (.not.associated(cur_cpl_var)) then
       write(iulog,*) '    Coupling-Variable-List_num : ',list%num_coupling_vars
       return
    else
       write(iulog,*) '    Coupling-Variable-List_num : ',list%num_coupling_vars
    end if

    do
       if (.not.associated(cur_cpl_var)) exit
       call CouplingVariablePrintInfo(cur_cpl_var)
       write(iulog,*)''
       cur_cpl_var => cur_cpl_var%next
    enddo
    write(iulog,*)''

  end subroutine CouplingVariableListPrintInfo

  !------------------------------------------------------------------------
  subroutine CouplingVariableListDestroy(list)
    !
    ! !DESCRIPTION:
    ! Release all allocated memory
    !
    implicit none
    !
    ! !ARGUMENTS
    type(coupling_variable_list_type), intent(inout) :: list
    !
    ! !LOCAL VARIABLES:
    type(coupling_variable_type), pointer            :: cur_cpl_var, prev_cpl_var
    
    cur_cpl_var => list%first
    do
       if (.not.associated(cur_cpl_var)) exit
       prev_cpl_var => cur_cpl_var
       cur_cpl_var  => cur_cpl_var%next
       call CouplingVariableDestroy(prev_cpl_var)
    enddo

    call CouplingVariableListCreate(list)

  end subroutine CouplingVariableListDestroy

#endif USE_PETSC_LIB

end module CouplingVariableType
