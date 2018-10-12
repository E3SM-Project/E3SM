module Coupler_module
 
#include "petsc/finclude/petscsys.h"
  use petscsys
  use Condition_module
  use Connection_module
  use Region_module
 
  use PFLOTRAN_Constants_module

  implicit none

  private
 
      
  ! coupler types
  PetscInt, parameter, public :: INITIAL_COUPLER_TYPE = 1
  PetscInt, parameter, public :: BOUNDARY_COUPLER_TYPE = 2
  PetscInt, parameter, public :: SRC_SINK_COUPLER_TYPE = 3
  PetscInt, parameter, public :: COUPLER_IPHASE_INDEX = 1

  type, public :: coupler_type
    PetscInt :: id                                      ! id of coupler
    character(len=MAXWORDLENGTH) :: name                ! name of coupler
    PetscInt :: itype                                   ! integer defining type
    character(len=MAXWORDLENGTH) :: ctype               ! character string defining type
    character(len=MAXWORDLENGTH) :: flow_condition_name ! character string defining name of condition to be applied
    character(len=MAXWORDLENGTH) :: tran_condition_name ! character string defining name of condition to be applied
    character(len=MAXWORDLENGTH) :: region_name         ! character string defining name of region to be applied
    character(len=MAXWORDLENGTH) :: well_spec_name      ! character string defining name of well_spec to be applied 
    PetscInt :: iflow_condition                         ! id of condition in condition array/list
    PetscInt :: itran_condition                         ! id of condition in condition array/list
    PetscInt :: iregion                                 ! id of region in region array/list
    PetscInt :: iface                                   ! for structured grids only
    PetscInt, pointer :: flow_aux_mapping(:)            ! maps flow_aux_real_var to primarhy dof
    PetscInt, pointer :: flow_bc_type(:)                ! id of boundary condition type
    PetscInt, pointer :: flow_aux_int_var(:,:)          ! auxiliary array for integer value
    PetscReal, pointer :: flow_aux_real_var(:,:)        ! auxiliary array for real values
    type(flow_condition_type), pointer :: flow_condition     ! pointer to condition in condition array/list
    type(tran_condition_type), pointer :: tran_condition     ! pointer to condition in condition array/list
    type(region_type), pointer :: region                ! pointer to region in region array/list
    type(connection_set_type), pointer :: connection_set ! pointer to an array/list of connections
    PetscInt :: numfaces_set
    type(coupler_type), pointer :: next                 ! pointer to next coupler
  end type coupler_type
  
  type, public :: coupler_ptr_type
    type(coupler_type), pointer :: ptr
  end type coupler_ptr_type
    
  type, public :: coupler_list_type
    PetscInt :: num_couplers
    type(coupler_type), pointer :: first
    type(coupler_type), pointer :: last
    type(coupler_ptr_type), pointer :: array(:)    
  end type coupler_list_type
  
  public :: CouplerCreate, &
            CouplerDestroy, &
            CouplerInitList, &
            CouplerAddToList, &
            CouplerRead, &
            CouplerDestroyList, &
            CouplerGetNumConnectionsInList, &
            CouplerListComputeConnections, &
            CouplerGetPtrFromList
  
  interface CouplerCreate
    module procedure CouplerCreate1
    module procedure CouplerCreate2
    module procedure CouplerCreateFromCoupler
  end interface
    
contains

! ************************************************************************** !

function CouplerCreate1()
  ! 
  ! CouplerCreate: Creates a coupler
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none

  type(coupler_type), pointer :: CouplerCreate1
  
  type(coupler_type), pointer :: coupler
  
  allocate(coupler)
  coupler%id = 0
  coupler%name = ''
  coupler%itype = BOUNDARY_COUPLER_TYPE
  coupler%ctype = "boundary"
  coupler%flow_condition_name = ""
  coupler%tran_condition_name = ""
  coupler%region_name = ""
  coupler%well_spec_name = "" 
  coupler%iflow_condition = 0
  coupler%itran_condition = 0
  coupler%iregion = 0
  coupler%iface = 0
  nullify(coupler%flow_aux_mapping)
  nullify(coupler%flow_bc_type)
  nullify(coupler%flow_aux_int_var)
  nullify(coupler%flow_aux_real_var)
  nullify(coupler%flow_condition)
  nullify(coupler%tran_condition)
  nullify(coupler%region)
  nullify(coupler%connection_set)
  nullify(coupler%next)
  
  CouplerCreate1 => coupler

end function CouplerCreate1

! ************************************************************************** !

function CouplerCreate2(itype)
  ! 
  ! Creates a coupler
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none

  PetscInt :: itype
  
  type(coupler_type), pointer :: CouplerCreate2
  
  type(coupler_type), pointer :: coupler
  
  coupler => CouplerCreate1()
  coupler%itype = itype
  select case(itype)
    case(INITIAL_COUPLER_TYPE)
      coupler%ctype = 'initial'
    case(BOUNDARY_COUPLER_TYPE)
      coupler%ctype = 'boundary'
    case(SRC_SINK_COUPLER_TYPE)
      coupler%ctype = 'source_sink'
  end select

  CouplerCreate2 => coupler

end function CouplerCreate2

! ************************************************************************** !

function CouplerCreateFromCoupler(coupler)
  ! 
  ! Creates a coupler
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none
  
  type(coupler_type), pointer :: coupler
  
  type(coupler_type), pointer :: CouplerCreateFromCoupler
  type(coupler_type), pointer :: new_coupler

  new_coupler => CouplerCreate1()

  new_coupler%id = coupler%id
  new_coupler%name = coupler%name
  new_coupler%itype = coupler%itype
  new_coupler%ctype = coupler%ctype
  new_coupler%flow_condition_name = coupler%flow_condition_name
  new_coupler%tran_condition_name = coupler%tran_condition_name
  new_coupler%region_name = coupler%region_name
  new_coupler%well_spec_name = coupler%well_spec_name
  new_coupler%iflow_condition = coupler%iflow_condition
  new_coupler%itran_condition = coupler%itran_condition
  new_coupler%iregion = coupler%iregion
  new_coupler%iface = coupler%iface

  ! these must remain null  
  nullify(coupler%flow_condition)
  nullify(coupler%tran_condition)
  nullify(coupler%region)
  nullify(coupler%flow_aux_mapping)
  nullify(coupler%flow_bc_type)
  nullify(coupler%flow_aux_int_var)
  nullify(coupler%flow_aux_real_var)
  nullify(coupler%connection_set)
  nullify(coupler%next)

  CouplerCreateFromCoupler => new_coupler

end function CouplerCreateFromCoupler

! ************************************************************************** !

subroutine CouplerInitList(list)
  ! 
  ! Initializes a coupler list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none

  type(coupler_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_couplers = 0

end subroutine CouplerInitList

! ************************************************************************** !

subroutine CouplerRead(coupler,input,option)
  ! 
  ! Reads a coupler from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(coupler_type) :: coupler
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','COUPLER')   
    call StringToUpper(word)      
    
    select case(trim(word))
    
      case('REGION','SURF_REGION')
        call InputReadWord(input,option,coupler%region_name,PETSC_TRUE)
      case('FLOW_CONDITION','SURF_FLOW_CONDITION')
        call InputReadWord(input,option,coupler%flow_condition_name,PETSC_TRUE)
      case('TRANSPORT_CONDITION')
        call InputReadWord(input,option,coupler%tran_condition_name,PETSC_TRUE)
      case default
        call InputKeywordUnrecognized(word,'coupler ',option)
    end select 
  
  enddo  

end subroutine CouplerRead

! ************************************************************************** !

subroutine CouplerAddToList(new_coupler,list)
  ! 
  ! Adds a new coupler to a coupler list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  type(coupler_type), pointer :: new_coupler
  type(coupler_list_type) :: list
  
  list%num_couplers = list%num_couplers + 1
  new_coupler%id = list%num_couplers
  if (.not.associated(list%first)) list%first => new_coupler
  if (associated(list%last)) list%last%next => new_coupler
  list%last => new_coupler
  
end subroutine CouplerAddToList

! ************************************************************************** !

subroutine CouplerListComputeConnections(grid,option,coupler_list)
  ! 
  ! computes connectivity for a list of couplers
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/20/08
  ! 

  use Option_module
  use Grid_module
  
  implicit none
 
  type(grid_type) :: grid
  type(option_type) :: option
  type(coupler_list_type), pointer :: coupler_list
  
  type(coupler_type), pointer :: coupler
  PetscInt :: offset
  
  if (.not.associated(coupler_list)) return
  
  offset = 0
  coupler => coupler_list%first
  do
    if (.not.associated(coupler)) exit 
    call CouplerComputeConnections(grid,option,coupler)
    if (associated(coupler%connection_set)) then
      coupler%connection_set%offset = offset
      offset = offset + coupler%connection_set%num_connections
    endif
    coupler => coupler%next
  enddo

end subroutine CouplerListComputeConnections

! ************************************************************************** !

subroutine CouplerComputeConnections(grid,option,coupler)
  ! 
  ! computes connectivity coupler to a grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/20/08
  ! 

  use Connection_module
  use Option_module
  use Region_module
  use Grid_module
  use Dataset_Base_class
  use Dataset_Gridded_HDF5_class
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Explicit_module, only : UGridExplicitSetBoundaryConnect, &
                                           UGridExplicitSetConnections
  
  implicit none
 
  type(grid_type) :: grid
  type(option_type) :: option
  type(coupler_type), pointer :: coupler_list
  
  PetscInt :: iconn
  PetscInt :: cell_id_local, cell_id_ghosted
  PetscInt :: connection_itype
  PetscInt :: iface
  type(connection_set_type), pointer :: connection_set
  type(region_type), pointer :: region
  type(coupler_type), pointer :: coupler
  PetscBool :: nullify_connection_set
  PetscErrorCode :: ierr

  if (.not.associated(coupler)) return
  
  nullify_connection_set = PETSC_FALSE
  select case(coupler%itype)
    case(INITIAL_COUPLER_TYPE)
      if (associated(coupler%flow_condition)) then
        if (associated(coupler%flow_condition%pressure)) then
          if (coupler%flow_condition%pressure%itype /= HYDROSTATIC_BC .and. &
              coupler%flow_condition%pressure%itype /= SEEPAGE_BC .and. &
              coupler%flow_condition%pressure%itype /= CONDUCTANCE_BC) then
            select type(selector => coupler%flow_condition%pressure%dataset)
              class is(dataset_gridded_hdf5_type)
              class default
                nullify_connection_set = PETSC_TRUE
            end select
          endif
        else if (associated(coupler%flow_condition%concentration)) then
          ! need to calculate connection set
        endif
        !geh: this is a workaround for defining temperature with a gridded
        !     dataset.  still need to set up the connections.
        if (associated(coupler%flow_condition%temperature)) then
          select type(selector => coupler%flow_condition%temperature%dataset)
            class is(dataset_gridded_hdf5_type)
              nullify_connection_set = PETSC_FALSE
          end select
        endif
      else
        nullify_connection_set = PETSC_TRUE
      endif
      connection_itype = INITIAL_CONNECTION_TYPE
    case(SRC_SINK_COUPLER_TYPE)
      connection_itype = SRC_SINK_CONNECTION_TYPE
    case(BOUNDARY_COUPLER_TYPE)
      connection_itype = BOUNDARY_CONNECTION_TYPE
  end select
  
  if (nullify_connection_set) then
    nullify(coupler%connection_set)
    return
  endif
  
  region => coupler%region

  select case(grid%itype)
    case(EXPLICIT_UNSTRUCTURED_GRID)
      if (associated(region%explicit_faceset)) then
        connection_set => &
          UGridExplicitSetBoundaryConnect(grid%unstructured_grid% &
                                            explicit_grid, &
                                          region%cell_ids, &
                                     region%explicit_faceset%face_centroids, &
                                     region%explicit_faceset%face_areas, &
                                     region%name,option)
      else
        connection_set => &
          UGridExplicitSetConnections(grid%unstructured_grid% &
                                        explicit_grid, &
                                      region%cell_ids, &
                                      connection_itype,option)
      endif
    case default
      connection_set => ConnectionCreate(region%num_cells,connection_itype)
    
      ! if using higher order advection, allocate associated arrays
      if (option%itranmode == EXPLICIT_ADVECTION .and. &
          option%transport%tvd_flux_limiter /= 1 .and. &  ! 1 = upwind
          connection_set%itype == BOUNDARY_CONNECTION_TYPE) then
        ! connections%id_up2 should remain null as it will not be used
        allocate(connection_set%id_dn2(size(connection_set%id_dn)))
        connection_set%id_dn2 = 0
      endif  

      iface = coupler%iface
      do iconn = 1,region%num_cells
    
        cell_id_local = region%cell_ids(iconn)
        if (associated(region%faces)) iface = region%faces(iconn)
    
        connection_set%id_dn(iconn) = cell_id_local

        call GridPopulateConnection(grid,connection_set,iface,iconn, &
                                    cell_id_local,option)
      enddo
  end select

  coupler%connection_set => connection_set
  nullify(connection_set)
 
end subroutine CouplerComputeConnections

! ************************************************************************** !

function CouplerGetNumConnectionsInList(list)
  ! 
  ! Returns the number of connections associated
  ! with all couplers in the list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/19/07
  ! 

  implicit none
  
  type(coupler_list_type) :: list
  
  PetscInt :: CouplerGetNumConnectionsInList
  type(coupler_type), pointer :: coupler
  
  CouplerGetNumConnectionsInList = 0
  coupler => list%first
  
  do
    if (.not.associated(coupler)) exit
    CouplerGetNumConnectionsInList = CouplerGetNumConnectionsInList + &
                                     coupler%connection_set%num_connections
    coupler => coupler%next
  enddo

end function CouplerGetNumConnectionsInList

! ************************************************************************** !

function CouplerGetPtrFromList(coupler_name,coupler_list)
  ! 
  ! Returns a pointer to the coupler matching
  ! coupler_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  use String_module

  implicit none
  
  type(coupler_type), pointer :: CouplerGetPtrFromList
  character(len=MAXWORDLENGTH) :: coupler_name
  PetscInt :: length
  type(coupler_list_type) :: coupler_list

  type(coupler_type), pointer :: coupler
    
  nullify(CouplerGetPtrFromList)

  coupler => coupler_list%first
  do 
    if (.not.associated(coupler)) exit
    length = len_trim(coupler_name)
    if (length == len_trim(coupler%name) .and. &
        StringCompare(coupler%name,coupler_name,length)) then
      CouplerGetPtrFromList => coupler
      return
    endif
    coupler => coupler%next
  enddo
  
end function CouplerGetPtrFromList

! ************************************************************************** !

subroutine CouplerDestroyList(coupler_list)
  ! 
  ! Deallocates a list of couplers
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  type(coupler_list_type), pointer :: coupler_list
  
  type(coupler_type), pointer :: coupler, prev_coupler
  
  if (.not.associated(coupler_list)) return
  
  coupler => coupler_list%first
  do 
    if (.not.associated(coupler)) exit
    prev_coupler => coupler
    coupler => coupler%next
    call CouplerDestroy(prev_coupler)
  enddo
  
  coupler_list%num_couplers = 0
  nullify(coupler_list%first)
  nullify(coupler_list%last)
  if (associated(coupler_list%array)) deallocate(coupler_list%array)
  nullify(coupler_list%array)
  
  deallocate(coupler_list)
  nullify(coupler_list)

end subroutine CouplerDestroyList

! ************************************************************************** !

subroutine CouplerDestroy(coupler)
  ! 
  ! Destroys a coupler
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 
  use Utility_module, only : DeallocateArray
  implicit none
  
  type(coupler_type), pointer :: coupler
  
  if (.not.associated(coupler)) return

  ! since the below are simply pointers to objects in list that have already
  ! or will be deallocated from the list, nullify instead of destroying
  
  nullify(coupler%flow_condition)     ! since these are simply pointers to 
  nullify(coupler%tran_condition)     ! since these are simply pointers to 
  nullify(coupler%region)        ! conditoins in list, nullify

  call DeallocateArray(coupler%flow_aux_mapping)
  call DeallocateArray(coupler%flow_bc_type)
  call DeallocateArray(coupler%flow_aux_int_var)
  call DeallocateArray(coupler%flow_aux_real_var)

  call ConnectionDestroy(coupler%connection_set)
  nullify(coupler%connection_set)

  deallocate(coupler)
  nullify(coupler)

end subroutine CouplerDestroy

end module Coupler_module
