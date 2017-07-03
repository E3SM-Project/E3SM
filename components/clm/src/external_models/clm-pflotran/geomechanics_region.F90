module Geomechanics_Region_module
 
  use Geometry_module
  use PFLOTRAN_Constants_module
  
  implicit none

  private

#include "petsc/finclude/petscsys.h"
 
  type, public :: gm_region_type
    PetscInt :: id
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXSTRINGLENGTH) :: filename
    type(point3d_type), pointer :: coordinates(:)
    PetscInt :: num_verts 
    PetscInt, pointer :: vertex_ids(:) 
    type(gm_region_type), pointer :: next
  end type gm_region_type
  
  type, public :: gm_region_ptr_type
    type(gm_region_type), pointer :: ptr
  end type gm_region_ptr_type
  
  type, public :: gm_region_list_type
    PetscInt :: num_regions
    type(gm_region_type), pointer :: first
    type(gm_region_type), pointer :: last
    type(gm_region_type), pointer :: array(:)
  end type gm_region_list_type

  interface GeomechRegionCreate
    module procedure GeomechRegionCreateWithList
    module procedure GeomechRegionCreateWithNothing
    module procedure GeomechRegionCreateWithGeomechRegion
  end interface GeomechRegionCreate
  
  interface GeomechRegionReadFromFile
    module procedure GeomechRegionReadFromFileId
    module procedure GeomechRegionReadFromFilename
  end interface GeomechRegionReadFromFile
  
   public :: GeomechRegionCreate, &
             GeomechRegionDestroy, &
             GeomechRegionAddToList, &
             GeomechRegionReadFromFile, &
             GeomechRegionDestroyList, &
             GeomechRegionRead, &
             GeomechRegionInitList, &
             GeomechRegionGetPtrFromList
             
 contains

! ************************************************************************** !

function GeomechRegionCreateWithNothing()
  ! 
  ! Creates a region with no arguments for
  ! geomechanics
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  implicit none
  
  type(gm_region_type), pointer :: GeomechRegionCreateWithNothing
  
  type(gm_region_type), pointer :: region
  
  allocate(region)
  region%id = 0
  region%name = ""
  region%filename = ""
  region%num_verts = 0
  nullify(region%coordinates)
  nullify(region%vertex_ids)
  nullify(region%next)
  
  GeomechRegionCreateWithNothing => region

end function GeomechRegionCreateWithNothing

! ************************************************************************** !

function GeomechRegionCreateWithList(list)
  ! 
  ! Creates a region from a list of vertices
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  implicit none
  
  PetscInt :: list(:)
  
  type(gm_region_type), pointer :: GeomechRegionCreateWithList
  
  type(gm_region_type), pointer :: region

  region => GeomechRegionCreateWithNothing()
  region%num_verts = size(list)
  allocate(region%vertex_ids(region%num_verts))
  region%vertex_ids = list
  
  GeomechRegionCreateWithList => region

end function GeomechRegionCreateWithList

! ************************************************************************** !

function GeomechRegionCreateWithGeomechRegion(region)
  ! 
  ! Creates a copy of a region
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  use Grid_Unstructured_Cell_module

  implicit none
  
  type(gm_region_type), pointer :: GeomechRegionCreateWithGeomechRegion
  type(gm_region_type), pointer :: region
  
  type(gm_region_type), pointer :: new_region
  PetscInt :: icount, temp_int
  
  new_region => GeomechRegionCreateWithNothing()
  
  new_region%id = region%id
  new_region%name = region%name
  new_region%filename = region%filename
  new_region%num_verts = region%num_verts
  if (associated(region%coordinates)) then
    call GeometryCopyCoordinates(region%coordinates, &
                                 new_region%coordinates)
  endif
  if (associated(region%vertex_ids)) then
    allocate(new_region%vertex_ids(new_region%num_verts))
    new_region%vertex_ids(1:new_region%num_verts) = &
    region%vertex_ids(1:new_region%num_verts)
  endif
 
  GeomechRegionCreateWithGeomechRegion => new_region
  
end function GeomechRegionCreateWithGeomechRegion

! ************************************************************************** !

subroutine GeomechRegionInitList(list)
  ! 
  ! Initializes a region list
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  implicit none

  type(gm_region_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_regions = 0

end subroutine GeomechRegionInitList

! ************************************************************************** !

subroutine GeomechRegionAddToList(new_region,list)
  ! 
  ! Adds a new region to a region list
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  implicit none
  
  type(gm_region_type), pointer :: new_region
  type(gm_region_list_type) :: list
  
  list%num_regions = list%num_regions + 1
  new_region%id = list%num_regions
  if (.not.associated(list%first)) list%first => new_region
  if (associated(list%last)) list%last%next => new_region
  list%last => new_region
  
end subroutine GeomechRegionAddToList

! ************************************************************************** !

subroutine GeomechRegionRead(region,input,option)
  ! 
  ! Reads a region from the input file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  use Input_Aux_module
  use String_module
  use Option_module
  
  implicit none
  
  type(option_type) :: option
  type(gm_region_type) :: region
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word
 
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS_REGION')
    call StringToUpper(keyword)   

    select case(trim(keyword))
      case('COORDINATE')
        allocate(region%coordinates(1))
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x) 
        if (InputError(input)) then
          input%ierr = 0
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'GEOMECHANICS_REGION')
          call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x)
        endif
        call InputErrorMsg(input,option,'x-coordinate','GEOMECHANICS_REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%y)
        call InputErrorMsg(input,option,'y-coordinate','GEOMECHANICS_REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%z)
        call InputErrorMsg(input,option,'z-coordinate','GEOMECHANICS_REGION')
      case('COORDINATES')
        call GeometryReadCoordinates(input,option,region%name, &
                                     region%coordinates)
      case('FILE')
        call InputReadNChars(input,option,region%filename, &
            MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','GEOMECHANICS_REGION')
        call GeomechRegionReadFromFilename(region,option,region%filename)
      case('LIST')
        option%io_buffer = 'GEOMECHANICS_REGION LIST currently not implemented'
        call printErrMsg(option)
      case default
        call InputKeywordUnrecognized(keyword,'GEOMECHANICS_REGION',option)
    end select
  enddo
 
end subroutine GeomechRegionRead

! ************************************************************************** !

subroutine GeomechRegionReadFromFilename(region,option,filename)
  ! 
  ! Reads a list of vertex ids from a file named
  ! filename
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  use Input_Aux_module
  use Option_module
  use Utility_module
  
  implicit none
  
  type(gm_region_type) :: region
  type(option_type) :: option
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: filename
  
  input => InputCreate(IUNIT_TEMP,filename,option)
  call GeomechRegionReadFromFileId(region,input,option)          
  call InputDestroy(input)         

end subroutine GeomechRegionReadFromFilename

! ************************************************************************** !

subroutine GeomechRegionReadFromFileId(region,input,option)
  ! 
  ! Reads a list of vertex ids from an open file
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  use Input_Aux_module
  use Option_module
  use Utility_module
  use Logging_module
  use Grid_Unstructured_Cell_module
  
  implicit none
  
  type(gm_region_type) :: region
  type(option_type) :: option
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=1) :: backslash
  character(len=MAXSTRINGLENGTH) :: string, string1

  PetscInt, pointer :: temp_int_array(:)
  PetscInt, pointer :: vertex_ids(:)
  PetscInt :: max_size
  PetscInt :: count
  PetscInt :: temp_int
  PetscInt :: input_data_type
  PetscInt :: ii
  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: remainder
  PetscErrorCode :: ierr

  max_size = 1000
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  allocate(temp_int_array(max_size))
  allocate(vertex_ids(max_size))
  
  temp_int_array = 0
  vertex_ids = 0
  
  count = 0
  call InputReadPflotranString(input, option)
  do 
    call InputReadInt(input, option, temp_int)
    if (InputError(input)) exit
    count = count + 1
    temp_int_array(count) = temp_int
  enddo

  if (count == 1) then
    !
    ! Input data contains only cell ids
    !
    vertex_ids(1) = temp_int_array(1)
    count = 1

    ! Read the data
    do
      call InputReadPflotranString(input, option)
      if (InputError(input)) exit
      call InputReadInt(input, option, temp_int)
      if (.not.InputError(input)) then
        count = count + 1
        vertex_ids(count) = temp_int
      endif
      if (count+1 > max_size) then ! resize temporary array
        call reallocateIntArray(vertex_ids, max_size)
      endif
    enddo

    ! Depending on processor rank, save only a portion of data
    region%num_verts = count/option%mycommsize
      remainder = count - region%num_verts*option%mycommsize
    if (option%myrank < remainder) region%num_verts = region%num_verts + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_verts, istart, ONE_INTEGER_MPI, MPIU_INTEGER, &
                    MPI_SUM, option%mycomm, ierr)
    call MPI_Scan(region%num_verts, iend, ONE_INTEGER_MPI, MPIU_INTEGER, &
                   MPI_SUM, option%mycomm, ierr)

    ! Allocate memory and save the data
    region%num_verts = iend - istart
    allocate(region%vertex_ids(region%num_verts))
    region%vertex_ids(1:region%num_verts) = vertex_ids(istart+1:iend)
    deallocate(vertex_ids)
  else
   option%io_buffer = 'Provide one vertex_id per line, GEOMECHANICS_REGION.'
   call printErrMsg(option) 
  endif
  
  deallocate(temp_int_array)
  
#ifdef GEOMECH_DEBUG
  write(string,*) option%myrank
  write(string1,*) region%name
  string = 'geomech_region_' // trim(adjustl(string1)) // '_vertex_ids' &
    // trim(adjustl(string)) // '.out'
  open(unit=86,file=trim(string))
  do ii = 1, region%num_verts
    write(86,'(i5)') region%vertex_ids(ii)
  enddo  
  close(86)
#endif    
    
end subroutine GeomechRegionReadFromFileId

! ************************************************************************** !

function GeomechRegionGetPtrFromList(region_name,region_list)
  ! 
  ! Returns a pointer to the region matching region_name
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  use String_module

  implicit none
  
  type(gm_region_type), pointer :: GeomechRegionGetPtrFromList
  character(len=MAXWORDLENGTH) :: region_name
  PetscInt :: length
  type(gm_region_list_type) :: region_list

  type(gm_region_type), pointer :: region
    
  nullify(GeomechRegionGetPtrFromList)
  region => region_list%first
  
  do 
    if (.not.associated(region)) exit
    length = len_trim(region_name)
    if (length == len_trim(region%name) .and. &
        StringCompare(region%name,region_name,length)) then
      GeomechRegionGetPtrFromList => region
      return
    endif
    region => region%next
  enddo
  
end function GeomechRegionGetPtrFromList

! ************************************************************************** !

subroutine GeomechRegionDestroyList(region_list)
  ! 
  ! Deallocates a list of regions
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  implicit none
  
  type(gm_region_list_type), pointer :: region_list
  
  type(gm_region_type), pointer :: region, prev_region
  
  if (.not.associated(region_list)) return
  
  region => region_list%first
  do 
    if (.not.associated(region)) exit
    prev_region => region
    region => region%next
    call GeomechRegionDestroy(prev_region)
  enddo
  
  region_list%num_regions = 0
  nullify(region_list%first)
  nullify(region_list%last)
  if (associated(region_list%array)) deallocate(region_list%array)
  nullify(region_list%array)
  
  deallocate(region_list)
  nullify(region_list)

end subroutine GeomechRegionDestroyList

! ************************************************************************** !

subroutine GeomechRegionDestroy(region)
  ! 
  ! Deallocates a region
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 06/06/13
  ! 

  implicit none
  
  type(gm_region_type), pointer :: region
  
  if (.not.associated(region)) return
  
  if (associated(region%vertex_ids)) deallocate(region%vertex_ids)
  nullify(region%vertex_ids)
  if (associated(region%coordinates)) deallocate(region%coordinates)
  nullify(region%coordinates)
  
  nullify(region%next)

  deallocate(region)
  nullify(region)
  
end subroutine GeomechRegionDestroy

end module Geomechanics_Region_module
