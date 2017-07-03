module Geomechanics_Coupler_module
 
  use Geomechanics_Condition_module
  use Geomechanics_Region_module
  use PFLOTRAN_Constants_module
 
  implicit none

  private
 
#include "petsc/finclude/petscsys.h"

  ! coupler types
  ! SK: Note that there is no initial coupler since we solve 
  ! a quasi-static problem for geomechanics (when coupled to flow, otherwise
  ! it is a steady state problem)
  PetscInt, parameter, public :: GM_BOUNDARY_COUPLER_TYPE = 1
  PetscInt, parameter, public :: GM_SRC_SINK_COUPLER_TYPE = 2

   type, public :: geomech_coupler_type
    PetscInt :: id                           ! id of coupler
    character(len=MAXWORDLENGTH) :: name                         ! name of coupler
    PetscInt :: itype                        ! integer defining type
    character(len=MAXWORDLENGTH) :: ctype                        ! character string defining type
    character(len=MAXWORDLENGTH) :: geomech_condition_name       ! character string defining name of condition to be applied
    character(len=MAXWORDLENGTH) :: region_name                  ! character string defining name of region to be applied
    PetscInt :: igeomech_condition           ! id of condition in condition array/list
    PetscInt :: iregion                      ! id of region in region array/list
    PetscInt, pointer :: geomech_aux_int_var(:,:)     ! auxiliary array for integer value
    PetscReal, pointer :: geomech_aux_real_var(:,:)    ! auxiliary array for real values
    type(geomech_condition_type), pointer :: geomech_condition            ! pointer to condition in condition array/list
    type(gm_region_type), pointer :: region                       ! pointer to region in region array/list
    type(geomech_coupler_type), pointer :: next                         ! pointer to next coupler
  end type geomech_coupler_type
  
  type, public :: geomech_coupler_ptr_type
    type(geomech_coupler_type), pointer :: ptr
  end type geomech_coupler_ptr_type
    
  type, public :: geomech_coupler_list_type
    PetscInt :: num_couplers
    type(geomech_coupler_type), pointer :: first
    type(geomech_coupler_type), pointer :: last
    type(geomech_coupler_type), pointer :: array(:)    
  end type geomech_coupler_list_type
  
  public :: GeomechCouplerCreate, &
            GeomechCouplerDestroy, &
            GeomechCouplerInitList, &
            GeomechCouplerAddToList, &
            GeomechCouplerRead, &
            GeomechCouplerDestroyList, &
            GeomechCouplerGetPtrFromList

  interface GeomechCouplerCreate
    module procedure GeomechCouplerCreate1
    module procedure GeomechCouplerCreate2
    module procedure GeomechCouplerCreateFromGeomechCoupler
  end interface
    
contains

! ************************************************************************** !

function GeomechCouplerCreate1()
  ! 
  ! GeomechCouplerCreate: Creates a coupler
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  implicit none

  type(geomech_coupler_type), pointer :: GeomechCouplerCreate1
  
  type(geomech_coupler_type), pointer :: coupler
  
  allocate(coupler)
  coupler%id = 0
  coupler%name = ''
  coupler%itype = GM_BOUNDARY_COUPLER_TYPE
  coupler%ctype = "boundary"
  coupler%geomech_condition_name = ""
  coupler%region_name = ""
  coupler%igeomech_condition = 0
  coupler%iregion = 0
  nullify(coupler%geomech_aux_int_var)
  nullify(coupler%geomech_aux_real_var)
  nullify(coupler%geomech_condition)
  nullify(coupler%region)
  nullify(coupler%next)
  
  GeomechCouplerCreate1 => coupler

end function GeomechCouplerCreate1

! ************************************************************************** !

function GeomechCouplerCreate2(itype)
  ! 
  ! Creates a coupler
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  implicit none

  PetscInt :: itype
  
  type(geomech_coupler_type), pointer :: GeomechCouplerCreate2
  
  type(geomech_coupler_type), pointer :: coupler
  
  coupler => GeomechCouplerCreate1()
  coupler%itype = itype
  select case(itype)
    case(GM_BOUNDARY_COUPLER_TYPE)
      coupler%ctype = 'boundary'
    case(GM_SRC_SINK_COUPLER_TYPE)
      coupler%ctype = 'source_sink'
  end select

  GeomechCouplerCreate2 => coupler

end function GeomechCouplerCreate2

! ************************************************************************** !

function GeomechCouplerCreateFromGeomechCoupler(coupler)
  ! 
  ! Creates a coupler
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  implicit none
  
  type(geomech_coupler_type), pointer :: coupler
  
  type(geomech_coupler_type), pointer :: GeomechCouplerCreateFromGeomechCoupler
  type(geomech_coupler_type), pointer :: new_coupler

  new_coupler => GeomechCouplerCreate1()

  new_coupler%id = coupler%id
  new_coupler%name = coupler%name
  new_coupler%itype = coupler%itype
  new_coupler%ctype = coupler%ctype
  new_coupler%geomech_condition_name = coupler%geomech_condition_name
  new_coupler%region_name = coupler%region_name
  new_coupler%igeomech_condition = coupler%igeomech_condition
  new_coupler%iregion = coupler%iregion

  ! these must remain null  
  nullify(coupler%geomech_condition)
  nullify(coupler%region)
  nullify(coupler%geomech_aux_int_var)
  nullify(coupler%geomech_aux_real_var)
  nullify(coupler%next)

  GeomechCouplerCreateFromGeomechCoupler => new_coupler

end function GeomechCouplerCreateFromGeomechCoupler

! ************************************************************************** !

subroutine GeomechCouplerInitList(list)
  ! 
  ! Initializes a coupler list
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  implicit none

  type(geomech_coupler_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_couplers = 0

end subroutine GeomechCouplerInitList

! ************************************************************************** !

subroutine GeomechCouplerRead(coupler,input,option)
  ! 
  ! Reads a coupler from the input file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(geomech_coupler_type) :: coupler
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS COUPLER')   
    call StringToUpper(word)      
    
    select case(trim(word))
    
      case('GEOMECHANICS_REGION')
        call InputReadWord(input,option,coupler%region_name,PETSC_TRUE)
      case('GEOMECHANICS_CONDITION')
        call InputReadWord(input,option,coupler%geomech_condition_name, &
                           PETSC_TRUE)
      case default
        call InputKeywordUnrecognized(word, &
                     'geomechanics coupler',option)
    end select 
  
  enddo  

end subroutine GeomechCouplerRead

! ************************************************************************** !

subroutine GeomechCouplerAddToList(new_coupler,list)
  ! 
  ! Adds a new coupler to a coupler list
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  implicit none
  
  type(geomech_coupler_type), pointer :: new_coupler
  type(geomech_coupler_list_type) :: list
  
  list%num_couplers = list%num_couplers + 1
  new_coupler%id = list%num_couplers
  if (.not.associated(list%first)) list%first => new_coupler
  if (associated(list%last)) list%last%next => new_coupler
  list%last => new_coupler
  
end subroutine GeomechCouplerAddToList

! ************************************************************************** !

function GeomechCouplerGetPtrFromList(coupler_name,coupler_list)
  ! 
  ! Returns a pointer to the geomech coupler
  ! matching coupler_name
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  use String_module

  implicit none
  
  type(geomech_coupler_type), pointer :: GeomechCouplerGetPtrFromList
  character(len=MAXWORDLENGTH) :: coupler_name
  PetscInt :: length
  type(geomech_coupler_list_type) :: coupler_list

  type(geomech_coupler_type), pointer :: coupler
    
  nullify(GeomechCouplerGetPtrFromList)

  coupler => coupler_list%first
  do 
    if (.not.associated(coupler)) exit
    length = len_trim(coupler_name)
    if (length == len_trim(coupler%name) .and. &
        StringCompare(coupler%name,coupler_name,length)) then
      GeomechCouplerGetPtrFromList => coupler
      return
    endif
    coupler => coupler%next
  enddo
  
end function GeomechCouplerGetPtrFromList

! ************************************************************************** !

subroutine GeomechCouplerDestroyList(coupler_list)
  ! 
  ! Deallocates a list of geomech couplers
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  implicit none
  
  type(geomech_coupler_list_type), pointer :: coupler_list
  
  type(geomech_coupler_type), pointer :: coupler, prev_coupler
  
  if (.not.associated(coupler_list)) return
  
  coupler => coupler_list%first
  do 
    if (.not.associated(coupler)) exit
    prev_coupler => coupler
    coupler => coupler%next
    call GeomechCouplerDestroy(prev_coupler)
  enddo
  
  coupler_list%num_couplers = 0
  nullify(coupler_list%first)
  nullify(coupler_list%last)
  if (associated(coupler_list%array)) deallocate(coupler_list%array)
  nullify(coupler_list%array)
  
  deallocate(coupler_list)
  nullify(coupler_list)

end subroutine GeomechCouplerDestroyList

! ************************************************************************** !

subroutine GeomechCouplerDestroy(coupler)
  ! 
  ! Destroys a coupler
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/13/13
  ! 

  implicit none
  
  type(geomech_coupler_type), pointer :: coupler
  
  if (.not.associated(coupler)) return
  
  ! since the below are simply pointers to objects in list that have already
  ! or will be deallocated from the list, nullify instead of destroying
  
  nullify(coupler%geomech_condition)     ! since these are simply pointers to 
  nullify(coupler%region)                ! conditions in list, nullify

  if (associated(coupler%geomech_aux_int_var)) &
    deallocate(coupler%geomech_aux_int_var)
  nullify(coupler%geomech_aux_int_var)
  if (associated(coupler%geomech_aux_real_var)) &
    deallocate(coupler%geomech_aux_real_var)
  nullify(coupler%geomech_aux_real_var)
   
  deallocate(coupler)
  nullify(coupler)

end subroutine GeomechCouplerDestroy

  
end module Geomechanics_Coupler_module
