module Region_module
 
  use Geometry_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"

  PetscInt, parameter, public :: DEFINED_BY_BLOCK = 1
  PetscInt, parameter, public :: DEFINED_BY_COORD = 2
  PetscInt, parameter, public :: DEFINED_BY_CELL_IDS = 3
  PetscInt, parameter, public :: DEFINED_BY_CELL_AND_FACE_IDS = 4
  PetscInt, parameter, public :: DEFINED_BY_VERTEX_IDS = 5
  PetscInt, parameter, public :: DEFINED_BY_SIDESET_UGRID = 6
  PetscInt, parameter, public :: DEFINED_BY_FACE_UGRID_EXP = 7
  PetscInt, parameter, public :: DEFINED_BY_POLY_BOUNDARY_FACE = 8
  PetscInt, parameter, public :: DEFINED_BY_POLY_CELL_CENTER = 9
  PetscInt, parameter, public :: DEFINED_BY_CARTESIAN_BOUNDARY = 10

  type, public :: block_type        
    PetscInt :: i1,i2,j1,j2,k1,k2    
    type(block_type), pointer :: next
  end type block_type
 
  type, public :: region_type
    PetscInt :: id
    PetscInt :: def_type
    character(len=MAXWORDLENGTH) :: name
    character(len=MAXSTRINGLENGTH) :: filename
    PetscInt :: i1,i2,j1,j2,k1,k2
    type(point3d_type), pointer :: coordinates(:)
    PetscInt :: iface
    PetscInt :: num_cells
    PetscInt, pointer :: cell_ids(:)
    PetscInt, pointer :: faces(:)
    !TODO(geh): Tear anything to do with structured/unstructured grids other
    !           than cell id ane face id out of region.
    PetscInt, pointer :: vertex_ids(:,:) ! For Unstructured mesh
    PetscInt :: num_verts              ! For Unstructured mesh
    type(region_sideset_type), pointer :: sideset
    type(region_explicit_face_type), pointer :: explicit_faceset
    type(polygonal_volume_type), pointer :: polygonal_volume
    type(region_type), pointer :: next
  end type region_type
  
  type, public :: region_ptr_type
    type(region_type), pointer :: ptr
  end type region_ptr_type
  
  type, public :: region_list_type
    PetscInt :: num_regions
    type(region_type), pointer :: first
    type(region_type), pointer :: last
    type(region_type), pointer :: array(:)
  end type region_list_type

  type, public :: region_sideset_type
    PetscInt :: nfaces
    PetscInt, pointer :: face_vertices(:,:)
  end type region_sideset_type
  
  type, public :: region_explicit_face_type
    type(point3d_type), pointer :: face_centroids(:)
    PetscReal, pointer :: face_areas(:)
  end type region_explicit_face_type
  
  interface RegionCreate
    module procedure RegionCreateWithBlock
    module procedure RegionCreateWithList
    module procedure RegionCreateWithNothing
    module procedure RegionCreateWithRegion    
  end interface RegionCreate
  
  interface RegionReadFromFile
    module procedure RegionReadFromFileId
    module procedure RegionReadFromFilename
    module procedure RegionReadSideSet
    module procedure RegionReadExplicitFaceSet
  end interface RegionReadFromFile
  
  public :: RegionCreate, &
            RegionRead, &
            RegionReadFromFile, &
            RegionInitList, &
            RegionAddToList, &
            RegionGetPtrFromList, & 
            RegionDestroyList, &
            RegionReadSideSet, &
            RegionCreateSideset, &
            RegionInputRecord, &
            RegionDestroy
  
contains

! ************************************************************************** !

function RegionCreateWithNothing()
  ! 
  ! Creates a region with no arguments
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none
  
  type(region_type), pointer :: RegionCreateWithNothing
  
  type(region_type), pointer :: region
  
  allocate(region)
  region%id = 0
  region%def_type = 0
  region%name = ""
  region%filename = ""
  region%i1 = 0
  region%i2 = 0
  region%j1 = 0
  region%j2 = 0
  region%k1 = 0
  region%k2 = 0
  region%iface = 0
  region%num_cells = 0
  ! By default it is assumed that the region is applicable to strucutred grid,
  ! unless explicitly stated in pflotran input file
  region%num_verts = 0
  nullify(region%coordinates)
  nullify(region%cell_ids)
  nullify(region%faces)
  nullify(region%vertex_ids)
  nullify(region%sideset)
  nullify(region%explicit_faceset)
  nullify(region%polygonal_volume)
  nullify(region%next)
  
  RegionCreateWithNothing => region

end function RegionCreateWithNothing

! ************************************************************************** !

function RegionCreateSideset()
  ! 
  ! Creates a sideset
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/11
  ! 

  implicit none
  
  type(region_sideset_type), pointer :: RegionCreateSideset
  
  type(region_sideset_type), pointer :: sideset
  
  allocate(sideset)
  sideset%nfaces = 0
  nullify(sideset%face_vertices)
  
  RegionCreateSideset => sideset

end function RegionCreateSideset

! ************************************************************************** !

function RegionCreateExplicitFaceSet()
  ! 
  ! Creates a sideset
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/11
  ! 

  implicit none
  
  type(region_explicit_face_type), pointer :: RegionCreateExplicitFaceSet
  
  type(region_explicit_face_type), pointer :: explicit_faceset
  
  allocate(explicit_faceset)
  nullify(explicit_faceset%face_centroids)
  nullify(explicit_faceset%face_areas)
  
  RegionCreateExplicitFaceSet => explicit_faceset

end function RegionCreateExplicitFaceSet

! ************************************************************************** !

function RegionCreateWithBlock(i1,i2,j1,j2,k1,k2)
  ! 
  ! Creates a region with i,j,k indices for arguments
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none
  
  PetscInt :: i1, i2, j1, j2, k1, k2
  
  type(region_type), pointer :: RegionCreateWithBlock

  type(region_type), pointer :: region
  
  region => RegionCreateWithNothing()
  region%i1 = i1
  region%i2 = i2
  region%j1 = j2
  region%j2 = j2
  region%k1 = k1
  region%k2 = k2
  region%num_cells = (abs(i2-i1)+1)*(abs(j2-j1)+1)* &
                                    (abs(k2-k1)+1)
                                    
  RegionCreateWithBlock => region                                    

end function RegionCreateWithBlock

! ************************************************************************** !

function RegionCreateWithList(list)
  ! 
  ! RegionCreate: Creates a region from a list of cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none
  
  PetscInt :: list(:)
  
  type(region_type), pointer :: RegionCreateWithList
  
  type(region_type), pointer :: region

  region => RegionCreateWithNothing()
  region%num_cells = size(list)
  allocate(region%cell_ids(region%num_cells))
  region%cell_ids = list
  
  RegionCreateWithList => region

end function RegionCreateWithList

! ************************************************************************** !

function RegionCreateWithRegion(region)
  ! 
  ! Creates a copy of a region
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/22/08
  ! 

  use Grid_Unstructured_Cell_module

  implicit none
  
  type(region_type), pointer :: RegionCreateWithRegion
  type(region_type), pointer :: region
  
  type(region_type), pointer :: new_region
  PetscInt :: icount, temp_int
  
  new_region => RegionCreateWithNothing()
  
  new_region%id = region%id
  new_region%def_type = region%def_type
  new_region%name = region%name
  new_region%filename = region%filename
  new_region%i1 = region%i1
  new_region%i2 = region%i2
  new_region%j1 = region%j1
  new_region%j2 = region%j2
  new_region%k1 = region%k1
  new_region%k2 = region%k2
  new_region%iface = region%iface
  new_region%num_cells = region%num_cells
  new_region%num_verts = region%num_verts
  if (associated(region%coordinates)) then
    call GeometryCopyCoordinates(region%coordinates, &
                                 new_region%coordinates)
  endif
  if (associated(region%cell_ids)) then
    allocate(new_region%cell_ids(new_region%num_cells))
    new_region%cell_ids(1:new_region%num_cells) = &
      region%cell_ids(1:region%num_cells)
  endif
  if (associated(region%faces)) then
    allocate(new_region%faces(new_region%num_cells))
    new_region%faces(1:new_region%num_cells) = &
      region%faces(1:region%num_cells)
  endif
  if (associated(region%vertex_ids)) then
    allocate(new_region%vertex_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts))
    new_region%vertex_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts) = &
    region%vertex_ids(0:MAX_VERT_PER_FACE,1:new_region%num_verts)
  endif
  if (associated(region%sideset)) then
    new_region%sideset => RegionCreateSideSet()
    new_region%sideset%nfaces = region%sideset%nfaces
    allocate(new_region%sideset%face_vertices( &
               size(region%sideset%face_vertices,1), &
               size(region%sideset%face_vertices,2)))
    new_region%sideset%face_vertices = region%sideset%face_vertices
  endif
  if (associated(region%explicit_faceset)) then
    new_region%explicit_faceset => RegionCreateExplicitFaceSet()
    allocate(new_region%explicit_faceset%face_centroids( &
               size(region%explicit_faceset%face_centroids)))
    new_region%explicit_faceset%face_centroids = &
      region%explicit_faceset%face_centroids
    do icount = 1, size(region%explicit_faceset%face_centroids)
      new_region%explicit_faceset%face_centroids(icount)%x = &
        region%explicit_faceset%face_centroids(icount)%x
      new_region%explicit_faceset%face_centroids(icount)%y = &
        region%explicit_faceset%face_centroids(icount)%y
      new_region%explicit_faceset%face_centroids(icount)%z = &
        region%explicit_faceset%face_centroids(icount)%z
    enddo
    allocate(new_region%explicit_faceset%face_areas( &
               size(region%explicit_faceset%face_areas)))
    new_region%explicit_faceset%face_areas = &
      region%explicit_faceset%face_areas
  endif
  if (associated(region%polygonal_volume)) then
    new_region%polygonal_volume => GeometryCreatePolygonalVolume()
    if (associated(region%polygonal_volume%xy_coordinates)) then
      call GeometryCopyCoordinates(region%polygonal_volume%xy_coordinates, &
                                   new_region%polygonal_volume%xy_coordinates)
    endif
    if (associated(region%polygonal_volume%xz_coordinates)) then
      call GeometryCopyCoordinates(region%polygonal_volume%xz_coordinates, &
                                   new_region%polygonal_volume%xz_coordinates)
    endif
    if (associated(region%polygonal_volume%yz_coordinates)) then
      call GeometryCopyCoordinates(region%polygonal_volume%yz_coordinates, &
                                   new_region%polygonal_volume%yz_coordinates)
    endif
  endif
  
  RegionCreateWithRegion => new_region
  
end function RegionCreateWithRegion

! ************************************************************************** !

subroutine RegionInitList(list)
  ! 
  ! Initializes a region list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  ! 

  implicit none

  type(region_list_type) :: list
  
  nullify(list%first)
  nullify(list%last)
  nullify(list%array)
  list%num_regions = 0

end subroutine RegionInitList

! ************************************************************************** !

subroutine RegionAddToList(new_region,list)
  ! 
  ! Adds a new region to a region list
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  ! 

  implicit none
  
  type(region_type), pointer :: new_region
  type(region_list_type) :: list
  
  list%num_regions = list%num_regions + 1
  new_region%id = list%num_regions
  if (.not.associated(list%first)) list%first => new_region
  if (associated(list%last)) list%last%next => new_region
  list%last => new_region
  
end subroutine RegionAddToList

! ************************************************************************** !

subroutine RegionRead(region,input,option)
  ! 
  ! Reads a region from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/20/08
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Input_Aux_module
  use String_module
  use Option_module
  use Grid_Structured_module
  
  implicit none
  
  type(option_type) :: option
  type(region_type) :: region
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: keyword, word

  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit
    
    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','REGION')
    call StringToUpper(keyword)   

    select case(trim(keyword))
    
      case('BLOCK')
        region%def_type = DEFINED_BY_BLOCK
        call InputReadInt(input,option,region%i1)
        if (InputError(input)) then
          input%ierr = 0
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'REGION')
          call InputReadInt(input,option,region%i1) 
        endif
        call InputErrorMsg(input,option,'i1','REGION')
        call InputReadInt(input,option,region%i2)
        call InputErrorMsg(input,option,'i2','REGION')
        call InputReadInt(input,option,region%j1)
        call InputErrorMsg(input,option,'j1','REGION')
        call InputReadInt(input,option,region%j2)
        call InputErrorMsg(input,option,'j2','REGION')
        call InputReadInt(input,option,region%k1)
        call InputErrorMsg(input,option,'k1','REGION')
        call InputReadInt(input,option,region%k2)
        call InputErrorMsg(input,option,'k2','REGION')
      case('CARTESIAN_BOUNDARY')
        region%def_type = DEFINED_BY_CARTESIAN_BOUNDARY
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'cartesian boundary face','REGION')
        call StringToUpper(word)
        select case(word)
          case('WEST')
            region%iface = WEST_FACE
          case('EAST')
            region%iface = EAST_FACE
          case('NORTH')
            region%iface = NORTH_FACE
          case('SOUTH')
            region%iface = SOUTH_FACE
          case('BOTTOM')
            region%iface = BOTTOM_FACE
          case('TOP')
            region%iface = TOP_FACE
          case default
            option%io_buffer = 'Cartesian boundary face "' // trim(word) // &
              '" not recognized.'
            call printErrMsg(option)
        end select
      case('COORDINATE')
        region%def_type = DEFINED_BY_COORD
        allocate(region%coordinates(1))
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x) 
        if (InputError(input)) then
          input%ierr = 0
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,'REGION')
          call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%x)
        endif
        call InputErrorMsg(input,option,'x-coordinate','REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%y)
        call InputErrorMsg(input,option,'y-coordinate','REGION')
        call InputReadDouble(input,option,region%coordinates(ONE_INTEGER)%z)
        call InputErrorMsg(input,option,'z-coordinate','REGION')
      case('COORDINATES')
        region%def_type = DEFINED_BY_COORD
        call GeometryReadCoordinates(input,option,region%name, &
                                     region%coordinates)
      case('POLYGON')
        if (.not.associated(region%polygonal_volume)) then
          region%polygonal_volume => GeometryCreatePolygonalVolume()
        endif
        do
          call InputReadPflotranString(input,option)
          if (InputError(input)) exit
          if (InputCheckExit(input,option)) exit
          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'keyword','REGION')
          call StringToUpper(word)   
          select case(trim(word))
            case('TYPE')
              call InputReadWord(input,option,word,PETSC_TRUE)
              call InputErrorMsg(input,option,'polygon type','REGION')
              call StringToUpper(word)
              select case(word)
                case('BOUNDARY_FACES_IN_VOLUME')
                  region%def_type = DEFINED_BY_POLY_BOUNDARY_FACE
                case('CELL_CENTERS_IN_VOLUME')
                  region%def_type = DEFINED_BY_POLY_CELL_CENTER
                case default
                  option%io_buffer = 'REGION->POLYGON->"' // trim(word) // &
                    '" not recognized.'
                  call printErrMsg(option)
              end select
            case('XY')
              call GeometryReadCoordinates(input,option,region%name, &
                                         region%polygonal_volume%xy_coordinates)
            case('XZ')
              call GeometryReadCoordinates(input,option,region%name, &
                                         region%polygonal_volume%xz_coordinates)
            case('YZ')
              call GeometryReadCoordinates(input,option,region%name, &
                                         region%polygonal_volume%yz_coordinates)
            case default
              option%io_buffer = 'Keyword not recognized for REGION POLYGON.'
              call printErrMsg(option)
          end select
        enddo
      case('FILE')
        call InputReadNChars(input,option,region%filename, &
                             MAXSTRINGLENGTH,PETSC_TRUE)
        call InputErrorMsg(input,option,'filename','REGION')
      case('LIST')
        option%io_buffer = 'REGION LIST currently not implemented'
        call printErrMsg(option)
      case('FACE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'face','REGION')
        call StringToUpper(word)
        select case(word)
          case('WEST')
            region%iface = WEST_FACE
          case('EAST')
            region%iface = EAST_FACE
          case('NORTH')
            region%iface = NORTH_FACE
          case('SOUTH')
            region%iface = SOUTH_FACE
          case('BOTTOM')
            region%iface = BOTTOM_FACE
          case('TOP')
            region%iface = TOP_FACE
          case default
            option%io_buffer = 'FACE "' // trim(word) // &
              '" not recognized.'
            call printErrMsg(option)
        end select
      case default
        call InputKeywordUnrecognized(keyword,'REGION',option)
    end select
  enddo
 
end subroutine RegionRead

! ************************************************************************** !

subroutine RegionReadFromFilename(region,option,filename)
  ! 
  ! Reads a list of cells from a file named filename
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  ! 

  use Input_Aux_module
  use Option_module
  use Utility_module
  
  implicit none
  
  type(region_type) :: region
  type(option_type) :: option
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: filename
  
  input => InputCreate(IUNIT_TEMP,filename,option)
  call RegionReadFromFileId(region,input,option)          
  call InputDestroy(input)         

end subroutine RegionReadFromFilename

! ************************************************************************** !

subroutine RegionReadFromFileId(region,input,option)
  ! 
  ! Reads a list of cells from an open file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/29/07
  ! 

#include <petsc/finclude/petscsys.h>
  use petscsys
  use Input_Aux_module
  use Option_module
  use Utility_module
  use Logging_module
  use Grid_Unstructured_Cell_module
  
  implicit none
  
  type(region_type) :: region
  type(option_type) :: option
  type(input_type), pointer :: input
  
  character(len=MAXWORDLENGTH) :: word
  character(len=1) :: backslash

  PetscInt, pointer :: temp_int_array(:)
  PetscInt, pointer :: cell_ids_p(:)
  PetscInt, pointer :: face_ids_p(:)
  PetscInt, pointer :: vert_id_0_p(:)
  PetscInt, pointer :: vert_id_1_p(:)
  PetscInt, pointer :: vert_id_2_p(:)
  PetscInt, pointer :: vert_id_3_p(:)
  PetscInt, pointer :: vert_id_4_p(:)
  PetscInt :: max_size
  PetscInt :: count
  PetscInt :: temp_int
  PetscInt :: input_data_type
  PetscInt :: ii
  PetscInt :: istart
  PetscInt :: iend
  PetscInt :: remainder
  PetscErrorCode :: ierr

  PetscInt, parameter :: CELL_IDS_ONLY = 1
  PetscInt, parameter :: CELL_IDS_WITH_FACE_IDS = 2
  PetscInt, parameter :: VERTEX_IDS = 3

  call PetscLogEventBegin(logging%event_region_read_ascii,ierr);CHKERRQ(ierr)
  
  !TODO(geh): clean and optimize this subroutine
  
  max_size = 1000
  backslash = achar(92)  ! 92 = "\" Some compilers choke on \" thinking it
                          ! is a double quote as in c/c++
  
  allocate(temp_int_array(max_size))
  allocate(cell_ids_p(max_size))
  allocate(face_ids_p(max_size))
  allocate(vert_id_0_p(max_size))
  allocate(vert_id_0_p(max_size))
  allocate(vert_id_1_p(max_size))
  allocate(vert_id_2_p(max_size))
  allocate(vert_id_3_p(max_size))
  allocate(vert_id_4_p(max_size))
  
  temp_int_array = 0
  cell_ids_p = 0
  face_ids_p = 0
  vert_id_0_p = 0
  vert_id_1_p = -1
  vert_id_2_p = -1
  vert_id_3_p = -1
  vert_id_4_p = -1
  
  
  count = 0

  ! Determine if region definition in the input data is one of the following:
  !  1) Contains cell ids only : Only ONE entry per line
  !  2) Contains cell ids and face ids: TWO entries per line
  !  3) Contains vertex ids that make up the face: MORE than two entries per
  !     line
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
    input_data_type = CELL_IDS_ONLY
    cell_ids_p(1) = temp_int_array(1)
    count = 1
    region%def_type = DEFINED_BY_CELL_IDS

    ! Read the data
    do
      call InputReadPflotranString(input, option)
      if (InputError(input)) exit
      call InputReadInt(input, option, temp_int)
      if (.not.InputError(input)) then
        count = count + 1
        cell_ids_p(count) = temp_int
      endif
      if (count+1 > max_size) then ! resize temporary array
        call reallocateIntArray(cell_ids_p, max_size)
      endif
    enddo

    ! Depending on processor rank, save only a portion of data
    region%num_cells = count/option%mycommsize
      remainder = count - region%num_cells*option%mycommsize
    if (option%myrank < remainder) region%num_cells = region%num_cells + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_cells, istart, ONE_INTEGER_MPI, MPIU_INTEGER, &
                    MPI_SUM, option%mycomm, ierr)
    call MPI_Scan(region%num_cells, iend, ONE_INTEGER_MPI, MPIU_INTEGER, &
                   MPI_SUM, option%mycomm, ierr)

    ! Allocate memory and save the data
    region%num_cells = iend - istart
    allocate(region%cell_ids(region%num_cells))
    region%cell_ids(1:region%num_cells) = cell_ids_p(istart+1:iend)
    deallocate(cell_ids_p)

  else if (count == 2) then
    !
    ! Input data contains cell ids + face ids
    !
    input_data_type = CELL_IDS_WITH_FACE_IDS
    cell_ids_p(1) = temp_int_array(1)
    face_ids_p(1) = temp_int_array(2)
    count = 1 ! reset the counter to represent the num of rows read
    region%def_type = DEFINED_BY_CELL_AND_FACE_IDS

    ! Read the data
    do
      call InputReadPflotranString(input, option)
      if (InputError(input)) exit
      call InputReadInt(input, option, temp_int)
      if (InputError(input)) exit
      count = count + 1
      cell_ids_p(count) = temp_int

      call InputReadInt(input,option,temp_int)
      if (InputError(input)) then
        option%io_buffer = 'ERROR while reading the region "' // &
          trim(region%name) // '" from file'
        call printErrMsg(option)
      endif
      face_ids_p(count) = temp_int
      if (count+1 > max_size) then ! resize temporary array
        call reallocateIntArray(cell_ids_p, max_size)
        call reallocateIntArray(face_ids_p, max_size)
      endif
    enddo

    ! Depending on processor rank, save only a portion of data
    region%num_cells = count/option%mycommsize
      remainder = count - region%num_cells*option%mycommsize
    if (option%myrank < remainder) region%num_cells = region%num_cells + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_cells,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(region%num_cells,iend,ONE_INTEGER_MPI,MPIU_INTEGER, &
                   MPI_SUM,option%mycomm,ierr)

    ! Allocate memory and save the data
    allocate(region%cell_ids(region%num_cells))
    allocate(region%faces(region%num_cells))
    region%cell_ids(1:region%num_cells) = cell_ids_p(istart + 1:iend)
    region%faces(1:region%num_cells) = face_ids_p(istart + 1:iend)
    deallocate(cell_ids_p)
    deallocate(face_ids_p)

  else
    !
    ! Input data contains vertices
    !
    input_data_type = VERTEX_IDS
    vert_id_0_p(1) = temp_int_array(1)
    vert_id_1_p(1) = temp_int_array(2)
    vert_id_2_p(1) = temp_int_array(3)
    vert_id_3_p(1) = temp_int_array(4)
    if (vert_id_0_p(1) == 4 ) vert_id_4_p(1) = temp_int_array(5)
    count = 1 ! reset the counter to represent the num of rows read
    region%def_type = DEFINED_BY_VERTEX_IDS

    ! Read the data
    do
      call InputReadPflotranString(input,option)
      if (InputError(input)) exit
      call InputReadInt(input,option,temp_int)
      if (InputError(input)) exit
      count = count + 1
      vert_id_0_p(count) = temp_int

      vert_id_4_p(count) = UNINITIALIZED_INTEGER
      do ii = 1, vert_id_0_p(count)
        call InputReadInt(input,option,temp_int)
        if (InputError(input)) then
          option%io_buffer = 'ERROR while reading the region "' // &
            trim(region%name) // '" from file'
          call printErrMsg(option)
        endif

        select case(ii)
          case(1)
            vert_id_1_p(count) = temp_int
          case(2)
            vert_id_2_p(count) = temp_int
          case(3)
            vert_id_3_p(count) = temp_int
          case(4)
            vert_id_4_p(count) = temp_int
        end select

        if (count+1 > max_size) then ! resize temporary array
          call reallocateIntArray(vert_id_0_p,max_size)
          call reallocateIntArray(vert_id_1_p,max_size)
          call reallocateIntArray(vert_id_2_p,max_size)
          call reallocateIntArray(vert_id_3_p,max_size)
          call reallocateIntArray(vert_id_4_p,max_size)
        endif
      enddo
    enddo

    ! Depending on processor rank, save only a portion of data
    region%num_verts = count/option%mycommsize
      remainder = count - region%num_verts*option%mycommsize
    if (option%myrank < remainder) region%num_verts = region%num_verts + 1
    istart = 0
    iend   = 0
    call MPI_Exscan(region%num_verts,istart,ONE_INTEGER_MPI,MPIU_INTEGER, &
                    MPI_SUM,option%mycomm,ierr)
    call MPI_Scan(region%num_verts,iend,ONE_INTEGER_MPI,MPIU_INTEGER, &
                   MPI_SUM,option%mycomm,ierr)

    ! Allocate memory and save the data
    region%num_verts = iend - istart
    allocate(region%vertex_ids(0:MAX_VERT_PER_FACE,1:region%num_verts))
    region%vertex_ids(0,1:region%num_verts) = vert_id_0_p(istart + 1: iend)
    region%vertex_ids(1,1:region%num_verts) = vert_id_1_p(istart + 1: iend)
    region%vertex_ids(2,1:region%num_verts) = vert_id_2_p(istart + 1: iend)
    region%vertex_ids(3,1:region%num_verts) = vert_id_3_p(istart + 1: iend)
    region%vertex_ids(4,1:region%num_verts) = vert_id_4_p(istart + 1: iend)
    deallocate(vert_id_0_p)
    deallocate(vert_id_1_p)
    deallocate(vert_id_2_p)
    deallocate(vert_id_3_p)
    deallocate(vert_id_4_p)

  endif
  
#if 0  
  count = 1
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    call InputReadInt(input,option,temp_int)
    if (.not.InputError(input)) then
      count = count + 1
      temp_int_array(count) = temp_int
      write(*,*) count,temp_int
    endif
    if (count+1 > max_size) then ! resize temporary array
      call reallocateIntArray(temp_int_array,max_size)
    endif
  enddo

  if (count > 0) then
    region%num_cells = count
    allocate(region%cell_ids(count))
    region%cell_ids(1:count) = temp_int_array(1:count)
  else
    region%num_cells = 0
    nullify(region%cell_ids)
  endif
#endif
  deallocate(temp_int_array) 

  call PetscLogEventEnd(logging%event_region_read_ascii,ierr);CHKERRQ(ierr)

end subroutine RegionReadFromFileId

! ************************************************************************** !

subroutine RegionReadSideSet(sideset,filename,option)
  ! 
  ! Reads an unstructured grid sideset
  ! 
  ! Author: Glenn Hammond
  ! Date: 12/19/11
  ! 

#include <petsc/finclude/petscsys.h>
  use petscsys
  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  type(region_sideset_type) :: sideset
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string, hint
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: num_faces_local_save
  PetscInt :: num_faces_local
  PetscInt :: num_to_read
  PetscInt, parameter :: max_nvert_per_face = 4
  PetscInt, allocatable :: temp_int_array(:,:)

  PetscInt :: iface, ivertex, irank, num_vertices
  PetscInt :: remainder
  PetscErrorCode :: ierr
  PetscMPIInt :: status_mpi(MPI_STATUS_SIZE)
  PetscMPIInt :: int_mpi
  PetscInt :: fileid
  
  fileid = 86
  input => InputCreate(fileid,filename,option)

! Format of sideset file
! type: T=triangle, Q=quadrilateral
! vertn(Q) = 4
! vertn(T) = 3
! -----------------------------------------------------------------
! num_faces  (integer)
! type vert1 vert2 ... vertn  ! for face 1 (integers)
! type vert1 vert2 ... vertn  ! for face 2
! ...
! ...
! type vert1 vert2 ... vertn  ! for face num_faces
! -----------------------------------------------------------------

  hint = 'Unstructured Sideset'

  call InputReadPflotranString(input,option)
  string = 'unstructured sideset'
  call InputReadStringErrorMsg(input,option,hint)  

  ! read num_faces
  call InputReadInt(input,option,sideset%nfaces)
  call InputErrorMsg(input,option,'number of faces',hint)

  ! divide faces across ranks
  num_faces_local = sideset%nfaces/option%mycommsize 
  num_faces_local_save = num_faces_local
  remainder = sideset%nfaces - num_faces_local*option%mycommsize
  if (option%myrank < remainder) num_faces_local = &
                                 num_faces_local + 1

  ! allocate array to store vertices for each faces
  allocate(sideset%face_vertices(max_nvert_per_face, &
                                 num_faces_local))
  sideset%face_vertices = UNINITIALIZED_INTEGER

  ! for now, read all faces from ASCII file through io_rank and communicate
  ! to other ranks
  if (option%myrank == option%io_rank) then
    allocate(temp_int_array(max_nvert_per_face, &
                            num_faces_local_save+1))
    ! read for other processors
    do irank = 0, option%mycommsize-1
      temp_int_array = UNINITIALIZED_INTEGER
      num_to_read = num_faces_local_save
      if (irank < remainder) num_to_read = num_to_read + 1

      do iface = 1, num_to_read
        ! read in the vertices defining the cell face
        call InputReadPflotranString(input,option)
        call InputReadStringErrorMsg(input,option,hint)  
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'face type',hint)
        call StringToUpper(word)
        select case(word)
          case('Q')
            num_vertices = 4
          case('T')
            num_vertices = 3
          case('L')
            num_vertices = 2
          case default
            option%io_buffer = 'Unknown face type: ' // trim(word)
        end select
        do ivertex = 1, num_vertices
          call InputReadInt(input,option,temp_int_array(ivertex,iface))
          call InputErrorMsg(input,option,'vertex id',hint)
        enddo
      enddo
      ! if the faces reside on io_rank
      if (irank == option%io_rank) then
#if UGRID_DEBUG
        write(string,*) num_faces_local
        string = trim(adjustl(string)) // ' faces stored on p0'
        print *, trim(string)
#endif
        sideset%nfaces = num_faces_local
        sideset%face_vertices(:,1:num_faces_local) = &
          temp_int_array(:,1:num_faces_local)
      else
        ! otherwise communicate to other ranks
#if UGRID_DEBUG
        write(string,*) num_to_read
        write(word,*) irank
        string = trim(adjustl(string)) // ' faces sent from p0 to p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
        int_mpi = num_to_read*max_nvert_per_face
        call MPI_Send(temp_int_array,int_mpi,MPIU_INTEGER,irank, &
                      num_to_read,option%mycomm,ierr)
      endif
    enddo
    deallocate(temp_int_array)
  else
    ! other ranks post the recv
#if UGRID_DEBUG
        write(string,*) num_faces_local
        write(word,*) option%myrank
        string = trim(adjustl(string)) // ' faces received from p0 at p' // &
                 trim(adjustl(word))
        print *, trim(string)
#endif
    sideset%nfaces = num_faces_local
    int_mpi = num_faces_local*max_nvert_per_face
    call MPI_Recv(sideset%face_vertices,int_mpi, &
                  MPIU_INTEGER,option%io_rank, &
                  MPI_ANY_TAG,option%mycomm,status_mpi,ierr)
  endif

!  unstructured_grid%nlmax = num_faces_local
!  unstructured_grid%num_vertices_local = num_vertices_local

  call InputDestroy(input)

end subroutine RegionReadSideSet

! ************************************************************************** !

subroutine RegionReadExplicitFaceSet(explicit_faceset,cell_ids,filename,option)
  ! 
  ! Reads an unstructured grid explicit region
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/18/12
  ! 
#include <petsc/finclude/petscsys.h>
  use petscsys
  use Input_Aux_module
  use Option_module
  use String_module
  
  implicit none
  
  type(region_explicit_face_type), pointer :: explicit_faceset
  PetscInt, pointer :: cell_ids(:)
  character(len=MAXSTRINGLENGTH) :: filename
  type(option_type) :: option
  
  type(input_type), pointer :: input
  character(len=MAXSTRINGLENGTH) :: string, hint
  character(len=MAXWORDLENGTH) :: word
  PetscInt :: fileid
  
  PetscInt :: num_connections
  PetscInt :: iconn
  
  explicit_faceset => RegionCreateExplicitFaceSet()
  
  fileid = 86
  input => InputCreate(fileid,filename,option)
  
! Format of explicit unstructured grid file
! id_ = integer
! x_, y_, z_, area_ = real
! definitions
! id_ = id of grid cell
! x_ = x coordinate of cell face
! y_ = y coordinate of cell face
! z_ = z coordinate of cell face
! area_ = area of grid cell face
! -----------------------------------------------------------------
! CONNECTIONS <integer>   integer = # connections (M)
! id_1 x_1 y_1 z_1 area_1
! id_2 x_2 y_2 z_2 area_2
! ...
! ...
! id_M x_M y_M z_M area_M
! -----------------------------------------------------------------

  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit

    call InputReadWord(input,option,word,PETSC_FALSE)
    call StringToUpper(word)
    hint = trim(word)
  
    select case(word)
      case('CONNECTIONS')
        hint = 'Explicit Unstructured Grid CONNECTIONS in file: ' // &
          trim(adjustl(filename))
        call InputReadInt(input,option,num_connections)
        call InputErrorMsg(input,option,'number of connections',hint)
        
        allocate(cell_ids(num_connections))
        cell_ids = 0
        allocate(explicit_faceset%face_areas(num_connections))
        explicit_faceset%face_areas = 0    
        allocate(explicit_faceset%face_centroids(num_connections))
        do iconn = 1, num_connections
          explicit_faceset%face_centroids(iconn)%x = 0.d0
          explicit_faceset%face_centroids(iconn)%y = 0.d0
          explicit_faceset%face_centroids(iconn)%z = 0.d0
        enddo
        do iconn = 1, num_connections
          call InputReadPflotranString(input,option)
          call InputReadStringErrorMsg(input,option,hint)  
          call InputReadInt(input,option,cell_ids(iconn))
          call InputErrorMsg(input,option,'cell id',hint)
          call InputReadDouble(input,option, &
                               explicit_faceset%face_centroids(iconn)%x)
          call InputErrorMsg(input,option,'face x coordinate',hint)
          call InputReadDouble(input,option, &
                               explicit_faceset%face_centroids(iconn)%y)
          call InputErrorMsg(input,option,'face y coordinate',hint)
          call InputReadDouble(input,option, &
                               explicit_faceset%face_centroids(iconn)%z)
          call InputErrorMsg(input,option,'face z coordinate',hint)
          call InputReadDouble(input,option, &
                               explicit_faceset%face_areas(iconn))
          call InputErrorMsg(input,option,'face area',hint)
        enddo
      case default
        call InputKeywordUnrecognized(word, &
               'REGION (explicit unstructured grid)',option)
    end select
  enddo

  call InputDestroy(input)

end subroutine RegionReadExplicitFaceSet

! ************************************************************************** !

function RegionGetPtrFromList(region_name,region_list)
  ! 
  ! Returns a pointer to the region matching region_name
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  use String_module

  implicit none
  
  type(region_type), pointer :: RegionGetPtrFromList
  character(len=MAXWORDLENGTH) :: region_name
  PetscInt :: length
  type(region_list_type) :: region_list

  type(region_type), pointer :: region
    
  nullify(RegionGetPtrFromList)
  region => region_list%first
  
  do 
    if (.not.associated(region)) exit
    length = len_trim(region_name)
    if (length == len_trim(region%name) .and. &
        StringCompare(region%name,region_name,length)) then
      RegionGetPtrFromList => region
      return
    endif
    region => region%next
  enddo
  
end function RegionGetPtrFromList

! **************************************************************************** !

subroutine RegionInputRecord(region_list)
  ! 
  ! Prints ingested region information to the input record file
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/30/2016
  ! 
  use Grid_Structured_module

  implicit none

  type(region_list_type), pointer :: region_list
  
  type(region_type), pointer :: cur_region
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: k
  PetscInt :: id = INPUT_RECORD_UNIT
  character(len=10) :: sFormat, iFormat
  
  sFormat = '(ES14.7)'
  iFormat = '(I10)'

  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'REGIONS'
  
  cur_region => region_list%first
  do
    if (.not.associated(cur_region)) exit
    write(id,'(a29)',advance='no') 'region: '
    write(id,'(a)') adjustl(trim(cur_region%name))
    if (len_trim(cur_region%filename) > 0) then
      write(id,'(a29)',advance='no') 'from file: '
      write(id,'(a)') adjustl(trim(cur_region%filename)) 
    endif
    
    select case (cur_region%def_type)
    !--------------------------------
      case (DEFINED_BY_BLOCK)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'BLOCK'
        write(id,'(a29)',advance='no') 'I indices: '
        write(word1,iFormat) cur_region%i1
        write(word2,iFormat) cur_region%i2
        write(id,'(a)') adjustl(trim(word1)) // ' ' // adjustl(trim(word2))
        write(id,'(a29)',advance='no') 'J indices: '
        write(word1,iFormat) cur_region%j1
        write(word2,iFormat) cur_region%j2
        write(id,'(a)') adjustl(trim(word1)) // ' ' // adjustl(trim(word2))
        write(id,'(a29)',advance='no') 'K indices: '
        write(word1,iFormat) cur_region%k1
        write(word2,iFormat) cur_region%k2
        write(id,'(a)') adjustl(trim(word1)) // ' ' // adjustl(trim(word2))
    !--------------------------------
      case (DEFINED_BY_CARTESIAN_BOUNDARY)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'CARTESIAN BOUNDARY'
    !--------------------------------
      case (DEFINED_BY_COORD)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'COORDINATE(S)'
        write(id,'(a29)',advance='no') 'X coordinate(s): '
        string = ''
        do k = 1,size(cur_region%coordinates)
         write(word1,sFormat) cur_region%coordinates(k)%x
         string = adjustl(trim(string)) // ' ' // adjustl(trim(word1))
        enddo
        write(id,'(a)') adjustl(trim(string)) // ' m'
        write(id,'(a29)',advance='no') 'Y coordinate(s): '
        string = ''
        do k = 1,size(cur_region%coordinates)
          write(word1,sFormat) cur_region%coordinates(k)%y
          string = adjustl(trim(string)) // ' ' // adjustl(trim(word1))
        enddo
        write(id,'(a)') adjustl(trim(string)) // ' m'
        write(id,'(a29)',advance='no') 'Z coordinate(s): '
        string = ''
        do k = 1,size(cur_region%coordinates)
          write(word1,sFormat) cur_region%coordinates(k)%z
          string = adjustl(trim(string)) // ' ' // adjustl(trim(word1))
        enddo
        write(id,'(a)') adjustl(trim(string)) // ' m'
    !--------------------------------
      case (DEFINED_BY_CELL_AND_FACE_IDS)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'CELL AND FACE IDS'
    !--------------------------------
      case (DEFINED_BY_CELL_IDS)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'CELL IDS'
    !--------------------------------
      case (DEFINED_BY_VERTEX_IDS)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'VERTEX IDS'
    !--------------------------------
      case (DEFINED_BY_FACE_UGRID_EXP)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'FACE UNSTRUCTURED GRID EXPLICIT'
    !--------------------------------
      case (DEFINED_BY_POLY_BOUNDARY_FACE)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'POLYGON BOUNDARY FACES IN VOLUME'
    !--------------------------------
      case (DEFINED_BY_POLY_CELL_CENTER)
        write(id,'(a29)',advance='no') 'defined by: '
        write(id,'(a)') 'POLYGON CELL CENTERS IN VOLUME'
    !--------------------------------
    end select
    
    if (cur_region%iface /= 0) then
      write(id,'(a29)',advance='no') 'face: '
      select case (cur_region%iface)
        case (WEST_FACE)
          write(id,'(a)') 'west'
        case (EAST_FACE)
          write(id,'(a)') 'east'
        case (NORTH_FACE)
          write(id,'(a)') 'north'
        case (SOUTH_FACE)
          write(id,'(a)') 'south'
        case (BOTTOM_FACE)
          write(id,'(a)') 'bottom'
        case (TOP_FACE)
          write(id,'(a)') 'top'
      end select
    endif
    
    write(id,'(a29)') '---------------------------: '
    cur_region => cur_region%next
  enddo
  
end subroutine RegionInputRecord

! **************************************************************************** !

subroutine RegionDestroySideset(sideset)
  ! 
  ! Deallocates a unstructured grid side set
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/09
  ! 

  implicit none
  
  type(region_sideset_type), pointer :: sideset
  
  if (.not.associated(sideset)) return
  
  if (associated(sideset%face_vertices)) deallocate(sideset%face_vertices)
  nullify(sideset%face_vertices)
  
  deallocate(sideset)
  nullify(sideset)
  
end subroutine RegionDestroySideset

! ************************************************************************** !

subroutine RegionDestroyExplicitFaceSet(explicit_faceset)
  ! 
  ! Deallocates a unstructured grid explicit grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/18/12
  ! 

  use Utility_module, only : DeallocateArray

  implicit none
  
  type(region_explicit_face_type), pointer :: explicit_faceset
  
  if (.not.associated(explicit_faceset)) return
  
  if (associated(explicit_faceset%face_centroids)) &
    deallocate(explicit_faceset%face_centroids)
  nullify(explicit_faceset%face_centroids)
  call DeallocateArray(explicit_faceset%face_areas)
  
  deallocate(explicit_faceset)
  nullify(explicit_faceset)
  
end subroutine RegionDestroyExplicitFaceSet

! ************************************************************************** !

subroutine RegionDestroyList(region_list)
  ! 
  ! Deallocates a list of regions
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  type(region_list_type), pointer :: region_list
  
  type(region_type), pointer :: region, prev_region
  
  if (.not.associated(region_list)) return
  
  region => region_list%first
  do 
    if (.not.associated(region)) exit
    prev_region => region
    region => region%next
    call RegionDestroy(prev_region)
  enddo
  
  region_list%num_regions = 0
  nullify(region_list%first)
  nullify(region_list%last)
  if (associated(region_list%array)) deallocate(region_list%array)
  nullify(region_list%array)
  
  deallocate(region_list)
  nullify(region_list)

end subroutine RegionDestroyList

! ************************************************************************** !

subroutine RegionDestroy(region)
  ! 
  ! Deallocates a region
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  type(region_type), pointer :: region
  
  if (.not.associated(region)) return
  
  call DeallocateArray(region%cell_ids)
  call DeallocateArray(region%faces)
  if (associated(region%coordinates)) deallocate(region%coordinates)
  nullify(region%coordinates)
  call RegionDestroySideset(region%sideset)
  call RegionDestroyExplicitFaceSet(region%explicit_faceset)
  call GeometryDestroyPolygonalVolume(region%polygonal_volume)
  
  if (associated(region%vertex_ids)) deallocate(region%vertex_ids)
  nullify(region%vertex_ids)
  
  nullify(region%next)

  deallocate(region)
  nullify(region)

end subroutine RegionDestroy

end module Region_module
