module Grid_Structured_module

#include "petsc/finclude/petscsys.h"
  use petscsys
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none
 
  private
 
! structured grid faces
  PetscInt, parameter, public :: NULL_FACE = 0
  PetscInt, parameter, public :: WEST_FACE = 1
  PetscInt, parameter, public :: EAST_FACE = 2
  PetscInt, parameter, public :: SOUTH_FACE = 3
  PetscInt, parameter, public :: NORTH_FACE = 4
  PetscInt, parameter, public :: BOTTOM_FACE = 5
  PetscInt, parameter, public :: TOP_FACE = 6

  PetscInt, parameter, public :: CARTESIAN_GRID = 3
  PetscInt, parameter, public :: CYLINDRICAL_GRID = 4
  PetscInt, parameter, public :: SPHERICAL_GRID = 5

  type, public :: grid_structured_type

    character(len=MAXWORDLENGTH) :: ctype
    PetscInt :: itype  ! type of grid (e.g. structured, unstructured, etc.)
    PetscInt :: global_offset
    PetscInt :: nx, ny, nz    ! Global domain dimensions of the grid.
    PetscInt :: nxy, nmax     ! nx * ny, nx * ny * nz
    PetscInt :: npx, npy, npz ! Processor partition in each direction.
    PetscInt :: npx_final, npy_final, npz_final ! actual decomposition
    PetscInt :: nlx, nly, nlz ! Local grid dimension w/o ghost nodes.
    PetscInt :: ngx, ngy, ngz ! Local grid dimension with ghost nodes.
    PetscInt :: lxs, lys, lzs 
      ! Global indices of non-ghosted corner (starting) of local domain.
    PetscInt :: gxs, gys, gzs
      ! Global indices of ghosted starting corner of local domain.
    PetscInt :: lxe, lye, lze, gxe, gye, gze
      ! Global indices of non-ghosted/ghosted ending corner of local domain.
    PetscInt :: nlxy, nlxz, nlyz
    PetscInt :: ngxy, ngxz, ngyz
    
    PetscInt :: istart, jstart, kstart, iend, jend, kend
      ! istart gives the ghosted local x-index of the non-ghosted starting 
      ! (lower left)corner. iend gives the local x-index of the non-ghosted 
      ! ending corner. jstart, jend correspond to y-index, kstart, kend to 
      ! z-index.  These are all zero-based indexing    

    PetscInt :: nlmax  ! Total number of non-ghosted nodes in local domain.
    PetscInt :: ngmax  ! Number of ghosted & non-ghosted nodes in local domain.
    PetscInt :: nlmax_faces  ! Total number of non-ghosted faces in local domain.
    PetscInt :: ngmax_faces  ! Number of ghosted & non-ghosted faces in local domain.

    PetscReal :: local_origin(3) ! local origin of non-ghosted grid
    PetscReal :: bounds(3,2)

    ! grid spacing for each direction for global domain
    PetscReal, pointer :: dx_global(:), dy_global(:), dz_global(:)
    ! grid spacing for each direction for local, ghosted domain
    PetscReal, pointer :: dxg_local(:), dyg_local(:), dzg_local(:)
    
    PetscBool :: invert_z_axis
    
    PetscReal, pointer :: dx(:), dy(:), dz(:)  ! ghosted grid spacings for each grid cell
    
    PetscInt, pointer :: cell_neighbors(:,:)
    
  end type grid_structured_type

  public :: StructGridCreate, &
            StructGridDestroy, &
            StructGridCreateDM, &
            StructGridComputeLocalBounds, &
            StructGridComputeInternConnect, &
            StructGridCreateVecFromDM, &
            StructGridMapIndices, &
            StructGridComputeSpacing, &
            StructGridComputeCoord, &
            StructGridReadDXYZ, &
            StructGridComputeVolumes, &
            StructGridPopulateConnection, &
            StructGridGetIJKFromCoordinate, &
            StructGridGetIJKFromLocalID, &
            StructGridGetIJKFromGhostedID, &
            StructGridGetGhostedNeighbors, &
            StructGridCreateTVDGhosts, &
            StructGridGetGhostedNeighborsCorners
            
contains

! ************************************************************************** !

function StructGridCreate()
  ! 
  ! Creates a structured grid object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/07
  ! 

  implicit none
  
  type(grid_structured_type), pointer :: StructGridCreate

  type(grid_structured_type), pointer :: structured_grid

  allocate(structured_grid)
  
  structured_grid%ctype = ''
  structured_grid%itype = 0
  structured_grid%global_offset = 0
  structured_grid%nx = 0
  structured_grid%ny = 0
  structured_grid%nz = 0
  structured_grid%npx = PETSC_DECIDE
  structured_grid%npy = PETSC_DECIDE
  structured_grid%npz = PETSC_DECIDE
  
  structured_grid%npx_final = 0
  structured_grid%npy_final = 0
  structured_grid%npz_final = 0

  structured_grid%nx = 0
  structured_grid%ny = 0
  structured_grid%nz = 0
  
  structured_grid%nxy = 0
  structured_grid%nmax = 0
  
  structured_grid%nlx = 0
  structured_grid%nly = 0
  structured_grid%nlz = 0
  structured_grid%nlxy = 0
  structured_grid%nlxz = 0
  structured_grid%nlyz = 0
  structured_grid%nlmax = 0
  
  structured_grid%ngx = 0
  structured_grid%ngy = 0
  structured_grid%ngz = 0
  structured_grid%ngxy = 0
  structured_grid%ngxz = 0
  structured_grid%ngyz = 0
  structured_grid%ngmax = 0

  structured_grid%lxs = 0
  structured_grid%lys = 0
  structured_grid%lzs = 0

  structured_grid%gxs = 0
  structured_grid%gys = 0
  structured_grid%gzs = 0

  structured_grid%lxe = 0
  structured_grid%lye = 0
  structured_grid%lze = 0

  structured_grid%gxe = 0
  structured_grid%gye = 0
  structured_grid%gze = 0

  structured_grid%istart = 0
  structured_grid%jstart = 0
  structured_grid%kstart = 0
  structured_grid%iend = 0
  structured_grid%jend = 0
  structured_grid%kend = 0
  
  nullify(structured_grid%dx_global)
  nullify(structured_grid%dy_global)
  nullify(structured_grid%dz_global)

  nullify(structured_grid%dxg_local)
  nullify(structured_grid%dyg_local)
  nullify(structured_grid%dzg_local)

  nullify(structured_grid%dx)
  nullify(structured_grid%dy)
  nullify(structured_grid%dz)
  
  nullify(structured_grid%cell_neighbors)
 
  structured_grid%local_origin = -1.d20
  structured_grid%bounds = -1.d20
  
  structured_grid%invert_z_axis = PETSC_FALSE
  
  StructGridCreate => structured_grid
  
end function StructGridCreate

! ************************************************************************** !

subroutine StructGridCreateDM(structured_grid,da,ndof,stencil_width, &
                              stencil_type,option)
  ! 
  ! StructGridCreateDMs: Creates structured distributed, parallel meshes/grids
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/22/07
  ! 
#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Option_module
#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif
        
  implicit none

  type(option_type) :: option
  type(grid_structured_type) :: structured_grid
  DM :: da
  PetscInt :: ndof
  PetscInt :: stencil_width
  PetscEnum :: stencil_type

  PetscErrorCode :: ierr

#ifdef CLM_PFLOTRAN
  PetscInt :: i, ncell
#endif

  !-----------------------------------------------------------------------
  ! Generate the DM object that will manage communication.
  !-----------------------------------------------------------------------
  ! This code is for the DMDACreate3D() interface in PETSc versions >= 3.2 --RTM
#ifndef CLM_PFLOTRAN
  call DMDACreate3D(option%mycomm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, &
                    DM_BOUNDARY_NONE,stencil_type, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                    structured_grid%npx,structured_grid%npy, &
                    structured_grid%npz,ndof,stencil_width, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    da,ierr);CHKERRQ(ierr)

#else

  if (option%mapping_files) then
    call DMDACreate3D(option%mycomm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, &
                    DM_BOUNDARY_NONE,stencil_type, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                    structured_grid%npx,structured_grid%npy, &
                    structured_grid%npz,ndof,stencil_width, &
                    PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                    da,ierr);CHKERRQ(ierr)

  else
    ! when coupling with CLM with NO mapping files,
    ! CLM domain is decomposed along X and Y directions (round-robin approach)
    ! in order to match both with each other, it's better to do direct (explicit) decompose as following
    ! (Fengming Yuan @ornl, 2016 Feb.)
    ! some checking
    if (structured_grid%npx /= size(clm_pf_idata%clm_lx)) then
      option%io_buffer = 'clm domain X-direction decompose incorrect '
      call printErrMsg(option)
    end if
    if (structured_grid%npy /= size(clm_pf_idata%clm_ly)) then
      option%io_buffer = 'clm domain Y-direction decompose incorrect '
      call printErrMsg(option)
    end if
    if (structured_grid%npz /= size(clm_pf_idata%clm_lz)) then
      option%io_buffer = 'clm domain Z-direction decompose incorrect '
      call printErrMsg(option)
    end if

    ncell = 0
    do i=1, structured_grid%npx
      ncell = ncell + clm_pf_idata%clm_lx(i)
    end do
    if (structured_grid%nx /= ncell) then
      option%io_buffer = 'clm domain decomposed X-direction cell no. sum NOT matches with grid%nx'
      call printErrMsg(option)
    end if

    ncell = 0
    do i=1, structured_grid%npy
      ncell = ncell + clm_pf_idata%clm_ly(i)
    end do
    if (structured_grid%ny /= ncell) then
      option%io_buffer = 'clm domain decomposed Y-direction cell no. sum NOT matches with grid%ny'
      call printErrMsg(option)
    end if

    ncell = 0
    do i=1, structured_grid%npz
      ncell = ncell + clm_pf_idata%clm_lz(i)
    end do
    if (structured_grid%nz /= ncell) then
      option%io_buffer = 'clm domain decomposed Z-direction cell no. sum NOT matches with grid%nz'
      call printErrMsg(option)
    end if

    call DMDACreate3D(option%mycomm,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, &
                    DM_BOUNDARY_NONE,stencil_type, &
                    structured_grid%nx,structured_grid%ny,structured_grid%nz, &
                    structured_grid%npx,structured_grid%npy,structured_grid%npz, &
                    ndof,stencil_width, &
                    clm_pf_idata%clm_lx, &
                    clm_pf_idata%clm_ly, &
                    clm_pf_idata%clm_lz, &
                    da,ierr);CHKERRQ(ierr)

  endif
#endif

  call DMSetFromOptions(da,ierr);CHKERRQ(ierr)
  call DMSetup(da,ierr);CHKERRQ(ierr)
  call DMDAGetInfo(da,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                   PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                   structured_grid%npx_final,structured_grid%npy_final, &
                   structured_grid%npz_final, &
                   PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                   PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER, &
                   PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)

end subroutine StructGridCreateDM

! ************************************************************************** !

subroutine StructGridComputeLocalBounds(structured_grid,da,option)
  ! 
  ! Computes the corners for the local portion
  ! of the structured grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/13/08
  ! 

#include "petsc/finclude/petscdm.h"
  use petscdm
  use Option_module
  
  implicit none

  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  DM :: da

  PetscErrorCode :: ierr

  ! get corner information
  call DMDAGetCorners(da, structured_grid%lxs, &
      structured_grid%lys, structured_grid%lzs, structured_grid%nlx, &
      structured_grid%nly, structured_grid%nlz, ierr);CHKERRQ(ierr)
     
  structured_grid%lxe = structured_grid%lxs + structured_grid%nlx
  structured_grid%lye = structured_grid%lys + structured_grid%nly
  structured_grid%lze = structured_grid%lzs + structured_grid%nlz
  structured_grid%nlxy = structured_grid%nlx * structured_grid%nly
  structured_grid%nlxz = structured_grid%nlx * structured_grid%nlz
  structured_grid%nlyz = structured_grid%nly * structured_grid%nlz
  structured_grid%nlmax = structured_grid%nlx * structured_grid%nly * structured_grid%nlz
     
  ! get ghosted corner information
  call DMDAGetGhostCorners(da, structured_grid%gxs, &
      structured_grid%gys, structured_grid%gzs, structured_grid%ngx, &
      structured_grid%ngy, structured_grid%ngz, ierr);CHKERRQ(ierr)
     
  structured_grid%gxe = structured_grid%gxs + structured_grid%ngx
  structured_grid%gye = structured_grid%gys + structured_grid%ngy
  structured_grid%gze = structured_grid%gzs + structured_grid%ngz
  structured_grid%ngxy = structured_grid%ngx * structured_grid%ngy
  structured_grid%ngxz = structured_grid%ngx * structured_grid%ngz
  structured_grid%ngyz = structured_grid%ngy * structured_grid%ngz
  structured_grid%ngmax = structured_grid%ngx * structured_grid%ngy * structured_grid%ngz
  
  structured_grid%global_offset = 0
  call MPI_Exscan(structured_grid%nlmax,structured_grid%global_offset, &
                  ONE_INTEGER_MPI,MPIU_INTEGER,MPI_SUM,option%mycomm,ierr)  
   
end subroutine StructGridComputeLocalBounds

! ************************************************************************** !

subroutine StructGridCreateVecFromDM(da,vector,vector_type)
  ! 
  ! Creates a PETSc vector from a DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/08
  ! 
#include "petsc/finclude/petscdm.h"
  use petscdm

  implicit none

  DM :: da
  Vec :: vector
  PetscInt :: vector_type
  
  PetscErrorCode :: ierr

  select case (vector_type)
    case(GLOBAL)
      call DMCreateGlobalVector(da,vector,ierr);CHKERRQ(ierr)
    case(LOCAL)
      call DMCreateLocalVector(da,vector,ierr);CHKERRQ(ierr)
    case(NATURAL)
      call DMDACreateNaturalVector(da,vector,ierr);CHKERRQ(ierr)
  end select

end subroutine StructGridCreateVecFromDM

! ************************************************************************** !

subroutine StructGridReadDXYZ(structured_grid,input,option)
  ! 
  ! Reads structured grid spacing from input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  use Utility_module, only : UtilityReadArray
  use Option_module
  use Input_Aux_module
  
  implicit none
  
  type(grid_structured_type) :: structured_grid
  type(input_type), pointer :: input
  type(option_type) :: option
  character(len=MAXSTRINGLENGTH) :: string
  
  PetscInt :: i

  ! the arrays are allocated within UtilityReadArray
  string = 'DXYZ : X'
  call UtilityReadArray(structured_grid%dx_global, &
                        structured_grid%nx,string,input,option)
  string = 'DXYZ : Y'
  call UtilityReadArray(structured_grid%dy_global, &
                        structured_grid%ny,string,input,option)
  string = 'DXYZ : Z'
  call UtilityReadArray(structured_grid%dz_global, &
                        structured_grid%nz,string,input,option)
    
  if (OptionPrintToFile(option)) then
    write(option%fid_out,'(/," *DXYZ ")')
    write(option%fid_out,'("  dx  ",/,(1p10e12.4))') &
      (structured_grid%dx_global(i),i=1,structured_grid%nx)
    write(option%fid_out,'("  dy  ",/,(1p10e12.4))') &
      (structured_grid%dy_global(i),i=1,structured_grid%ny)
    write(option%fid_out,'("  dz  ",/,(1p10e12.4))') &
      (structured_grid%dz_global(i),i=1,structured_grid%nz)
  endif

end subroutine StructGridReadDXYZ

! ************************************************************************** !

subroutine StructGridComputeSpacing(structured_grid,origin_global,option)
  ! 
  ! Computes structured grid spacing
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/26/07
  ! 

  use Option_module
  
  implicit none
  
  type(grid_structured_type) :: structured_grid
  PetscReal :: origin_global(3)
  type(option_type) :: option
  
  PetscInt :: i, j, k, ghosted_id
  PetscReal :: tempreal
  PetscErrorCode :: ierr

  allocate(structured_grid%dxg_local(structured_grid%ngx))
  structured_grid%dxg_local = 0.d0
  allocate(structured_grid%dyg_local(structured_grid%ngy))
  structured_grid%dyg_local = 0.d0
  allocate(structured_grid%dzg_local(structured_grid%ngz))
  structured_grid%dzg_local = 0.d0
  
  if (.not.associated(structured_grid%dx_global)) then
    ! indicates that the grid spacings still need to be computed
    if (structured_grid%bounds(1,1) < -1.d19) then 
      option%io_buffer = 'Bounds have not been set for grid and DXYZ ' // & 
        'does not exist'
      call printErrMsg(option)
    endif
    allocate(structured_grid%dx_global(structured_grid%nx))
    allocate(structured_grid%dy_global(structured_grid%ny))
    allocate(structured_grid%dz_global(structured_grid%nz))
      
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        structured_grid%dx_global = &
          (structured_grid%bounds(X_DIRECTION,UPPER)- &
           structured_grid%bounds(X_DIRECTION,LOWER)) / &
          dble(structured_grid%nx)
        structured_grid%dy_global = &
          (structured_grid%bounds(Y_DIRECTION,UPPER)- &
           structured_grid%bounds(Y_DIRECTION,LOWER)) / &
          dble(structured_grid%ny)
        structured_grid%dz_global = &
          (structured_grid%bounds(Z_DIRECTION,UPPER)- &
           structured_grid%bounds(Z_DIRECTION,LOWER)) / &
          dble(structured_grid%nz)
      case(CYLINDRICAL_GRID)
        structured_grid%dx_global = &
          (structured_grid%bounds(X_DIRECTION,UPPER)- &
           structured_grid%bounds(X_DIRECTION,LOWER)) / &
          dble(structured_grid%nx)
        structured_grid%dy_global = 1.d0
        structured_grid%dz_global = &
          (structured_grid%bounds(Z_DIRECTION,UPPER)- &
           structured_grid%bounds(Z_DIRECTION,LOWER)) / &
          dble(structured_grid%nz)
      case(SPHERICAL_GRID)
        structured_grid%dx_global = &
          (structured_grid%bounds(X_DIRECTION,UPPER)- &
           structured_grid%bounds(X_DIRECTION,LOWER)) / &
          dble(structured_grid%nx)
        structured_grid%dy_global = 1.d0
        structured_grid%dz_global = 1.d0
    end select
  else
    ! x-direction
    if (structured_grid%itype == CARTESIAN_GRID .or. &
        structured_grid%itype == CYLINDRICAL_GRID .or. &
        structured_grid%itype == SPHERICAL_GRID) then
      tempreal = origin_global(X_DIRECTION)
      structured_grid%bounds(X_DIRECTION,LOWER) = tempreal
      do i = 1, structured_grid%nx
        tempreal = tempreal + structured_grid%dx_global(i)
      enddo
      structured_grid%bounds(X_DIRECTION,UPPER) = tempreal
    endif
    ! y-direction
    if (structured_grid%itype == CARTESIAN_GRID) then
      tempreal = origin_global(Y_DIRECTION)
      structured_grid%bounds(Y_DIRECTION,LOWER) = tempreal
      do j = 1, structured_grid%ny
        tempreal = tempreal + structured_grid%dy_global(j)
      enddo
      structured_grid%bounds(Y_DIRECTION,UPPER) = tempreal
    else
      structured_grid%bounds(Y_DIRECTION,LOWER) = 0.d0
      structured_grid%bounds(Y_DIRECTION,UPPER) = 1.d0
    endif
    ! z-direction    
    if (structured_grid%itype == CARTESIAN_GRID .or. &
        structured_grid%itype == CYLINDRICAL_GRID) then
      tempreal = origin_global(Z_DIRECTION)
      structured_grid%bounds(Z_DIRECTION,LOWER) = tempreal
      do k = 1, structured_grid%nz
        tempreal = tempreal + structured_grid%dz_global(k)
      enddo
      structured_grid%bounds(Z_DIRECTION,UPPER) = tempreal
    else
      structured_grid%bounds(Y_DIRECTION,LOWER) = 0.d0
      structured_grid%bounds(Z_DIRECTION,UPPER) = 1.d0
    endif
  endif

  option%io_buffer = 'Domain Bounds (x y z):'
  call printMsg(option)
  write(option%io_buffer,'(2x,3es18.10)') structured_grid%bounds(:,LOWER)
  call printMsg(option)
  write(option%io_buffer,'(2x,3es18.10)') structured_grid%bounds(:,UPPER)
  call printMsg(option)

  structured_grid%dxg_local(1:structured_grid%ngx) = &
    structured_grid%dx_global(structured_grid%gxs+1:structured_grid%gxe)
  structured_grid%dyg_local(1:structured_grid%ngy) = &
    structured_grid%dy_global(structured_grid%gys+1:structured_grid%gye)
  structured_grid%dzg_local(1:structured_grid%ngz) = &
    structured_grid%dz_global(structured_grid%gzs+1:structured_grid%gze)
        
  allocate(structured_grid%dx(structured_grid%ngmax))
  structured_grid%dx = 0.d0
  allocate(structured_grid%dy(structured_grid%ngmax))
  structured_grid%dy = 0.d0
  allocate(structured_grid%dz(structured_grid%ngmax))
  structured_grid%dz = 0.d0
 
  do k = 1, structured_grid%ngz
    do j = 1, structured_grid%ngy
      do i = 1, structured_grid%ngx
        ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
        structured_grid%dx(ghosted_id) = structured_grid%dxg_local(i)
        structured_grid%dy(ghosted_id) = structured_grid%dyg_local(j)
        structured_grid%dz(ghosted_id) = structured_grid%dzg_local(k)
      enddo
    enddo
  enddo
  
end subroutine StructGridComputeSpacing

! ************************************************************************** !

subroutine StructGridComputeCoord(structured_grid,option,origin_global, &
                                      grid_x,grid_y,grid_z, &
                                      x_min,x_max,y_min,y_max,z_min,z_max)
  ! 
  ! Computes structured coordinates in x,y,z
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

  use Option_module
  
implicit none
  
  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  PetscReal :: origin_global(3)
  PetscReal :: grid_x(:), grid_y(:), grid_z(:)
  PetscReal :: x_min, x_max, y_min, y_max, z_min, z_max

! PetscErrorCode :: ierr
  PetscInt :: i, j, k, ghosted_id
  PetscReal :: x, y, z
  PetscInt :: prevnode

  x_min = origin_global(X_DIRECTION)
  y_min = origin_global(Y_DIRECTION)
  z_min = origin_global(Z_DIRECTION)
    
  do i=1,structured_grid%lxs
    x_min = x_min + structured_grid%dx_global(i)
  enddo
  do j=1,structured_grid%lys
    y_min = y_min + structured_grid%dy_global(j)
  enddo
  do k=1,structured_grid%lzs
    z_min = z_min + structured_grid%dz_global(k)
  enddo
   
  ! set min and max bounds of domain in coordinate directions
  structured_grid%local_origin(X_DIRECTION) = x_min
  structured_grid%local_origin(Y_DIRECTION) = y_min
  structured_grid%local_origin(Z_DIRECTION) = z_min
  x_max = x_min
  y_max = y_min
  z_max = z_min
  
  ! there is an issue with cumulative round off error where the summed 
  ! maximum extend may not match the global bound on the domain. this is 
  ! an issue if a region (e.g. point is located on the upper global bound.
  ! to present this, we set the maximum extend (at the global upper 
  ! boundary) to the non-summed value
  if (structured_grid%lxe < structured_grid%nx) then
    do i=structured_grid%istart,structured_grid%iend
      x_max = x_max + structured_grid%dxg_local(i+1)
    enddo
  else 
    x_max = structured_grid%bounds(X_DIRECTION,UPPER) 
  endif
  if (structured_grid%lye < structured_grid%ny) then
    do j=structured_grid%jstart,structured_grid%jend
      y_max = y_max + structured_grid%dyg_local(j+1)
    enddo
  else 
    y_max = structured_grid%bounds(Y_DIRECTION,UPPER) 
  endif
  if (structured_grid%lze < structured_grid%nz) then
    do k=structured_grid%kstart,structured_grid%kend
      z_max = z_max + structured_grid%dzg_local(k+1)
    enddo
  else 
    z_max = structured_grid%bounds(Z_DIRECTION,UPPER) 
  endif

! fill in grid cell coordinates
  ghosted_id = 0
  if (structured_grid%kstart > 0) then
    z = structured_grid%local_origin(Z_DIRECTION) - &
        0.5d0*structured_grid%dzg_local(1)
  else
    z = structured_grid%local_origin(Z_DIRECTION) + &
        0.5d0*structured_grid%dzg_local(1)
  endif
  do k=1, structured_grid%ngz
    if (structured_grid%jstart > 0) then
      y = structured_grid%local_origin(Y_DIRECTION) - &
          0.5d0*structured_grid%dyg_local(1)
    else
      y = structured_grid%local_origin(Y_DIRECTION) + &
          0.5d0*structured_grid%dyg_local(1)
    endif
    do j=1, structured_grid%ngy
      if (structured_grid%istart > 0) then
        x = structured_grid%local_origin(X_DIRECTION) - &
            0.5d0*structured_grid%dxg_local(1)
      else
        x = structured_grid%local_origin(X_DIRECTION) + &
            0.5d0*structured_grid%dxg_local(1)
      endif
      do i=1, structured_grid%ngx
        ghosted_id = ghosted_id + 1
        grid_x(ghosted_id) = x
        grid_y(ghosted_id) = y
        grid_z(ghosted_id) = z
        if (i < structured_grid%ngx) &
          x = x + &
             0.5d0*(structured_grid%dxg_local(i)+structured_grid%dxg_local(i+1))
      enddo
      if (j < structured_grid%ngy) &
        y = y + &
            0.5d0*(structured_grid%dyg_local(j)+structured_grid%dyg_local(j+1))
    enddo
    if (k < structured_grid%ngz) &
      z = z + &
          0.5d0*(structured_grid%dzg_local(k)+structured_grid%dzg_local(k+1))
  enddo
    
end subroutine StructGridComputeCoord

! ************************************************************************** !

subroutine StructGridGetIJKFromCoordinate(structured_grid,x,y,z,i,j,k)
  ! 
  ! Finds local, non-ghosted i,j,k indices for
  ! grid cell encompassing coordinate
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/27/08
  ! 

  use Option_module
  
  implicit none
    
  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: i, j, k
  PetscInt :: i_local, j_local, k_local
  PetscInt :: i_ghosted, j_ghosted, k_ghosted
  PetscReal :: x, y, z
  
  PetscReal :: x_upper_face, y_upper_face, z_upper_face

  i = -1
  j = -1
  k = -1

  x_upper_face = structured_grid%local_origin(X_DIRECTION)
  i_local = 1
  do i_ghosted=structured_grid%istart,structured_grid%iend
    if (x >= x_upper_face .and. &                   ! since i_ghosted is zero-based
        x <= x_upper_face+structured_grid%dxg_local(i_ghosted+1)) then
      ! test to prevent multiple procs from including a coordinate located on
      ! boundary of local decomposition shared by two procs
      ! if first cell in x-dir on proc
      if (i_ghosted == structured_grid%istart) then
        ! located on upwind boundary and ghosted
          if (Equal(x,x_upper_face) .and. &
            structured_grid%lxs /= structured_grid%gxs) exit
      endif
      i = i_local 
      exit
    endif
    i_local = i_local + 1
    x_upper_face = x_upper_face + structured_grid%dxg_local(i_ghosted+1)
  enddo
  y_upper_face = structured_grid%local_origin(Y_DIRECTION)
  j_local = 1
  do j_ghosted=structured_grid%jstart,structured_grid%jend
    if (y >= y_upper_face .and. &
        y <= y_upper_face+structured_grid%dyg_local(j_ghosted+1)) then
      ! test to prevent multiple procs from including a coordinate located on
      ! boundary of local decomposition shared by two procs
      ! if first cell in y-dir on proc
      if (j_ghosted == structured_grid%jstart) then
        ! located on upwind boundary and ghosted
        if (Equal(y,y_upper_face) .and. &
            structured_grid%lys /= structured_grid%gys) exit
      endif
      j = j_local
      exit
    endif
    j_local = j_local + 1
    y_upper_face = y_upper_face + structured_grid%dyg_local(j_ghosted+1)
  enddo
  z_upper_face = structured_grid%local_origin(Z_DIRECTION)
  k_local = 1
  do k_ghosted=structured_grid%kstart,structured_grid%kend
    if (z >= z_upper_face .and. &
        z <= z_upper_face+structured_grid%dzg_local(k_ghosted+1)) then
      ! test to prevent multiple procs from including a coordinate located on
      ! boundary of local decomposition shared by two procs
      ! if first cell in z-dir on proc
      if (k_ghosted == structured_grid%kstart) then
        ! if located on upwind boundary and ghosted, skip
        if (Equal(z,z_upper_face) .and. &
          structured_grid%lzs /= structured_grid%gzs) exit
      endif          
      k = k_local
      exit
    endif
    k_local = k_local + 1
    z_upper_face = z_upper_face + structured_grid%dzg_local(k_ghosted+1)
  enddo
    
end subroutine StructGridGetIJKFromCoordinate

! ************************************************************************** !

subroutine StructGridGetIJKFromLocalID(structured_grid,local_id,i,j,k)
  ! 
  ! Finds i,j,k indices for grid cell defined by
  ! local_id
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/08
  ! 

  use Option_module
  
  implicit none
  
  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: local_id
  PetscInt :: i, j, k
  
  k= int((local_id-1)/structured_grid%nlxy) + 1
  j= int(mod((local_id-1),structured_grid%nlxy)/structured_grid%nlx) + 1
  i= mod(mod((local_id-1),structured_grid%nlxy),structured_grid%nlx) + 1  
  
end subroutine StructGridGetIJKFromLocalID

! ************************************************************************** !

subroutine StructGridGetIJKFromGhostedID(structured_grid,ghosted_id,i,j,k)
  ! 
  ! Finds i,j,k indices for grid cell defined by
  ! a ghosted id
  ! 
  ! Author: Glenn Hammond
  ! Date: 04/11/08
  ! 

  use Option_module
  
  implicit none
  
  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscInt :: i, j, k
  
  k= int((ghosted_id-1)/structured_grid%ngxy) + 1
  j= int(mod((ghosted_id-1),structured_grid%ngxy)/structured_grid%ngx) + 1
  i= mod(mod((ghosted_id-1),structured_grid%ngxy),structured_grid%ngx) + 1  
  
end subroutine StructGridGetIJKFromGhostedID

! ************************************************************************** !

function StructGridGetLocalIDFromIJK(structured_grid,i,j,k)
  ! 
  ! Finds local_id for grid cell defined by
  ! i,j,k indices
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/11
  ! 

  use Option_module
  
  implicit none
  
  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: i, j, k
  
  PetscInt :: StructGridGetLocalIDFromIJK
  
  StructGridGetLocalIDFromIJK = &
    i+(j-1)*structured_grid%nlx+(k-1)*structured_grid%nlxy
  
end function StructGridGetLocalIDFromIJK

! ************************************************************************** !

function StructGridGetGhostedIDFromIJK(structured_grid,i,j,k)
  ! 
  ! Finds ghosted_id for grid cell defined by
  ! i,j,k indices
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/11
  ! 

  use Option_module
  
  implicit none
  
  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: i, j, k
  
  PetscInt :: StructGridGetGhostedIDFromIJK
  
  StructGridGetGhostedIDFromIJK = &
    i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
  
end function StructGridGetGhostedIDFromIJK

! ************************************************************************** !

function StructGridComputeInternConnect(structured_grid, xc, yc, zc, option)
  ! 
  ! computes internal connectivity of a
  ! structured grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/17/07
  ! 

  use Connection_module
  use Option_module
  
  implicit none

!  PetscReal :: radius(:)
  type(connection_set_type), pointer :: StructGridComputeInternConnect
  type(option_type) :: option
  type(grid_structured_type) :: structured_grid
  PetscReal, pointer :: xc(:),yc(:),zc(:)
  
  PetscReal, parameter :: Pi=3.141592653590d0
  
  PetscInt :: i, j, k, iconn, id_up, id_dn, id_up2, id_dn2
  PetscInt :: nconn
  PetscInt :: lenx, leny, lenz
  PetscInt :: tvd_ghost_offset, ghost_count
  PetscReal :: dist_up, dist_dn
  PetscReal :: r1, r2
  type(connection_set_type), pointer :: connections, connections_2

  PetscErrorCode :: ierr
  PetscReal, pointer :: radius(:)
  PetscInt, pointer :: int_array1(:), int_array2(:),int_array3(:),int_array4(:),int_array5(:),index(:)
  PetscInt :: count

  PetscReal :: tempreal

  radius => xc

  ! the adjustments in the case of AMR are based on the PIMS code adjustments by LC
  nconn = (structured_grid%ngx-1)*structured_grid%nly*structured_grid%nlz+ &
          structured_grid%nlx*(structured_grid%ngy-1)*structured_grid%nlz+ &
          structured_grid%nlx*structured_grid%nly*(structured_grid%ngz-1)
  
  structured_grid%nlmax_faces = 0
  structured_grid%ngmax_faces = 0

  lenx = structured_grid%ngx - 1
  leny = structured_grid%ngy - 1
  lenz = structured_grid%ngz - 1

  connections => ConnectionCreate(nconn,INTERNAL_CONNECTION_TYPE)
  
  ! if using higher order advection, allocate associated arrays
  if (option%itranmode == EXPLICIT_ADVECTION .and. &
      option%transport%tvd_flux_limiter /= 1) then  ! 1 = upwind
    allocate(connections%id_up2(size(connections%id_up)))
    allocate(connections%id_dn2(size(connections%id_dn)))
    connections%id_up2 = 0
    connections%id_dn2 = 0
  endif

  iconn = 0
  tvd_ghost_offset = 0
  ghost_count = 0
  ! x-connections
  if (structured_grid%ngx > 1) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = 1, lenx
              iconn = iconn+1
              id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy
              id_dn = id_up + 1
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              
              if (associated(connections%id_up2)) then
                if (i == 1) then
                  ! id_up indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_up2 = 1 + j + k*structured_grid%nly
                  ghost_count = ghost_count + 1
                  id_up2 = ghost_count
                  connections%id_up2(iconn) = -id_up2
                else
                  connections%id_up2(iconn) = id_up - 1
                endif
                if (i == lenx) then
                  ! id_dn indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_dn2 = 1 + j + k*structured_grid%nly + structured_grid%nlyz
                  id_dn2 = ghost_count + structured_grid%nlyz
                  connections%id_dn2(iconn) = -id_dn2
                else
                  connections%id_dn2(iconn) = id_dn + 1
                endif
              endif
              
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dx(id_up)
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(1,iconn) = 1.d0  ! x component of unit vector
              connections%area(iconn) = structured_grid%dy(id_up)* &
                                        structured_grid%dz(id_up)

#ifdef CLM_PFLOTRAN
              ! need to adjust 'dist' by elevation (z-coordinate) differences, if any, when coupled with CLM
              ! by coupling, grid%z is CLM grid's elevation in meters
              tempreal = abs(zc(id_up) - zc(id_dn))
              if (tempreal>1.d-6) then    !'1.d-6' is the limit of mircometer
                tempreal = tempreal*tempreal+(dist_up+dist_dn)*(dist_up+dist_dn)
                connections%dist(0,iconn) = sqrt(tempreal)
              endif

              if(.not.option%mapping_files) then
                ! dy is the length of height in a isoceles trapezoid grid, so need to adjust the connection (interface) face area
                tempreal = structured_grid%dy(id_dn) - structured_grid%dy(id_up)
                if (tempreal>1.d-6) then    ! 'up' is the short parallel side
                  tempreal = abs(tempreal)*(dist_up/(dist_up+dist_dn))+structured_grid%dy(id_up)
                  connections%area(iconn) = tempreal* &
                                    structured_grid%dz(id_up)

                elseif (tempreal<-1.d-6) then    ! 'dn' is the short parallel side
                  tempreal = abs(tempreal)*(dist_dn/(dist_up+dist_dn))+structured_grid%dy(id_dn)
                  connections%area(iconn) = tempreal* &
                                    structured_grid%dz(id_up)
                endif

              endif
#endif

            enddo
          enddo
        enddo
        tvd_ghost_offset = 2*structured_grid%nlyz ! west & east
        ghost_count = tvd_ghost_offset
      case(CYLINDRICAL_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = 1, lenx
              iconn = iconn+1
              id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy
              id_dn = id_up + 1
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dx(id_up)
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(1,iconn) = 1.d0  ! x component of unit vector
              connections%area(iconn) = 2.d0 * pi * (radius(id_up)+0.5d0*structured_grid%dx(id_up))* &
                                        structured_grid%dz(id_up)
            enddo
          enddo
        enddo
      case(SPHERICAL_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do j = structured_grid%jstart, structured_grid%jend
            do i = 1, lenx
              iconn = iconn+1
              id_up = i + j * structured_grid%ngx + k * structured_grid%ngxy
              id_dn = id_up + 1
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dx(id_up)
              dist_dn = 0.5d0*structured_grid%dx(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(1,iconn) = 1.d0  ! x component of unit vector
              connections%area(iconn) = 4.d0 * pi * (radius(id_up)+0.5d0*structured_grid%dx(id_up))
            enddo
          enddo
        enddo
    end select
  endif

  ! y-connections
  if (structured_grid%ngy > 1) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        do k = structured_grid%kstart, structured_grid%kend
          do i = structured_grid%istart, structured_grid%iend
            do j = 1, leny
              iconn = iconn+1

              id_up = i + 1 + (j-1) * structured_grid%ngx + k * structured_grid%ngxy
              id_dn = id_up + structured_grid%ngx
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              
              if (associated(connections%id_up2)) then
                if (j == 1) then
                  ! id_up indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_up2 = 1 + i + k*structured_grid%nlx + tvd_ghost_offset
                  ghost_count = ghost_count + 1
                  id_up2 = ghost_count
                  connections%id_up2(iconn) = -id_up2
                else
                  connections%id_up2(iconn) = id_up - structured_grid%ngx
                endif
                if (j == leny) then
                  ! id_dn indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_dn2 = 1 + i + k*structured_grid%nlx + &
!                           structured_grid%nlxz + tvd_ghost_offset
                  id_dn2 = ghost_count + structured_grid%nlxz
                  connections%id_dn2(iconn) = -id_dn2
                else
                  connections%id_dn2(iconn) = id_dn + structured_grid%ngx
                endif
              endif
                            
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dy(id_up)
              dist_dn = 0.5d0*structured_grid%dy(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(2,iconn) = 1.d0  ! y component of unit vector
              connections%area(iconn) = structured_grid%dx(id_up)* &
                                    structured_grid%dz(id_up)


#ifdef CLM_PFLOTRAN
              ! need to adjust 'dist' by elevation (z-coordinate) differences, if any, when coupled with CLM
              ! by coupling, grid%z is CLM grid's elevation in meters
              tempreal = abs(zc(id_up) - zc(id_dn))
              if (tempreal>1.d-6) then    !'1.d-6' is mircometer
                tempreal = tempreal*tempreal+(dist_up+dist_dn)*(dist_up+dist_dn)
                connections%dist(0,iconn) = sqrt(tempreal)
              endif

              if(.not.option%mapping_files) then
                ! dx is the mid-length of a isoceles trapezoid grid when from CLM, so need to adjust the connection (interface) face area
                tempreal = structured_grid%dx(id_dn) - structured_grid%dx(id_up)
                if (tempreal>1.d-6) then    ! 'up' is the short parallel side
                  tempreal = abs(tempreal)*(dist_up/(dist_up+dist_dn))+structured_grid%dx(id_up)
                  connections%area(iconn) = tempreal* &
                                    structured_grid%dz(id_up)
                elseif (tempreal<-1.d-6) then    ! 'dn' is the short parallel side
                  tempreal = abs(tempreal)*(dist_dn/(dist_up+dist_dn))+structured_grid%dx(id_dn)
                  connections%area(iconn) = tempreal* &
                                    structured_grid%dz(id_up)
                endif

              endif
#endif

            enddo
          enddo
        enddo
        tvd_ghost_offset = tvd_ghost_offset + &
          2*structured_grid%nlxz ! south & north
        ghost_count = tvd_ghost_offset
      case(CYLINDRICAL_GRID)
        option%io_buffer = 'For cylindrical coordinates, NY must be equal to 1.'
        call printErrMsg(option)
      case(SPHERICAL_GRID)
        option%io_buffer = 'For spherical coordinates, NY must be equal to 1.'
        call printErrMsg(option)
    end select
  endif
      
  ! z-connections
  if (structured_grid%ngz > 1) then
    select case(structured_grid%itype)
      case(CARTESIAN_GRID)
        do j = structured_grid%jstart, structured_grid%jend
          do i = structured_grid%istart, structured_grid%iend
            do k = 1, lenz
              iconn = iconn+1

              id_up = i + 1 + j * structured_grid%ngx + (k-1) * &
                  structured_grid%ngxy 
              id_dn = id_up + structured_grid%ngxy
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              
              if (associated(connections%id_up2)) then
                if (k == 1) then
!                  id_up2 = 1 + i + j*structured_grid%nlx + tvd_ghost_offset
                  ghost_count = ghost_count + 1
                  id_up2 = ghost_count
                  connections%id_up2(iconn) = -id_up2
                else
                  connections%id_up2(iconn) = id_up - structured_grid%ngxy
                endif
                if (k == lenz) then
                  ! id_dn indexes tvd_ghost_vec, see StructGridCreateTVDGhosts() 
!                  id_dn2 = 1 + i + j*structured_grid%nlx + &
!                           structured_grid%nlxy + tvd_ghost_offset
                  id_dn2 = ghost_count + structured_grid%nlxy
                  connections%id_dn2(iconn) = -id_dn2
                else
                  connections%id_dn2(iconn) = id_dn + structured_grid%ngxy
                endif
              endif
                                 
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dz(id_up)
              dist_dn = 0.5d0*structured_grid%dz(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(3,iconn) = 1.d0  ! z component of unit vector
              connections%area(iconn) = structured_grid%dx(id_up) * &
                                        structured_grid%dy(id_up)
            enddo
          enddo
        enddo
      case(CYLINDRICAL_GRID)
        do j = structured_grid%jstart, structured_grid%jend
          do i = structured_grid%istart, structured_grid%iend
            do k = 1, lenz
              iconn = iconn+1
              id_up = i + 1 + j * structured_grid%ngx + (k-1) * &
                  structured_grid%ngxy
              id_dn = id_up + structured_grid%ngxy
              connections%id_up(iconn) = id_up
              connections%id_dn(iconn) = id_dn
              connections%dist(-1:3,iconn) = 0.d0
              dist_up = 0.5d0*structured_grid%dz(id_up)
              dist_dn = 0.5d0*structured_grid%dz(id_dn)
              connections%dist(-1,iconn) = dist_up/(dist_up+dist_dn)
              connections%dist(0,iconn) = dist_up+dist_dn
              connections%dist(3,iconn) = 1.d0  ! z component of unit vector
              ! pi*(r2^2-r1^2)
              r2 = xc(id_up) + 0.5d0*structured_grid%dx(id_up)
              r1 = xc(id_up) - 0.5d0*structured_grid%dx(id_up)
              connections%area(iconn) = pi * dabs(r2*r2 - r1*r1)
            enddo
          enddo
        enddo
      case(SPHERICAL_GRID)
        option%io_buffer = 'For spherical coordinates, NZ must be equal to 1.'
        call printErrMsg(option)
  end select
  endif

  StructGridComputeInternConnect => connections

end function StructGridComputeInternConnect

! ************************************************************************** !

subroutine StructGridPopulateConnection(radius,structured_grid,connection,iface, &
                                        iconn,ghosted_id,option)
  ! 
  ! Computes details of connection (area, dist, etc)
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/09/07
  ! 

  use Connection_module
  use Option_module
  
  implicit none
 
  type(grid_structured_type) :: structured_grid
  type(connection_set_type) :: connection
  PetscInt :: iface
  PetscInt :: iconn
  PetscInt :: ghosted_id
  type(option_type) :: option
  PetscReal :: radius(:)
  
  PetscErrorCode :: ierr
  
  PetscReal, parameter :: Pi=3.141592653590d0
  
  select case(connection%itype)
    case(BOUNDARY_CONNECTION_TYPE)
    
      select case(iface)

        case(WEST_FACE,EAST_FACE)

          select case(structured_grid%itype)
            case(CARTESIAN_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dx(ghosted_id)
              connection%area(iconn) = structured_grid%dy(ghosted_id)* &
                                   structured_grid%dz(ghosted_id)


              if (iface ==  WEST_FACE) then
                connection%dist(1,iconn) = 1.d0
              else
                connection%dist(1,iconn) = -1.d0
              endif
              if (associated(connection%id_dn2)) then
                  connection%id_dn2(iconn) = &
                    StructGetTVDGhostConnection(ghosted_id, &
                                                structured_grid, &
                                                iface,option)
              endif              
            case(CYLINDRICAL_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dx(ghosted_id)
              if (iface ==  WEST_FACE) then
                connection%dist(1,iconn) = 1.d0
                connection%area(iconn) = 2.d0 * pi * (radius(ghosted_id)- &
                                      0.5d0*structured_grid%dx(ghosted_id)) * &
                                      structured_grid%dz(ghosted_id)
              else
                connection%dist(1,iconn) = -1.d0
                connection%area(iconn) = 2.d0 * pi * (radius(ghosted_id)+ &
                                      0.5d0*structured_grid%dx(ghosted_id)) * &
                                      structured_grid%dz(ghosted_id)
              endif
            case(SPHERICAL_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dx(ghosted_id)
              if (iface ==  WEST_FACE) then
                connection%dist(1,iconn) = 1.d0
                connection%area(iconn) = 0.d0
              else
                connection%dist(1,iconn) = -1.d0
                connection%area(iconn) = 4.d0 * pi * (radius(ghosted_id)+ &
                                         0.5d0*structured_grid%dx(ghosted_id))
              endif
          end select

        case(SOUTH_FACE,NORTH_FACE)

          select case(structured_grid%itype)
            case(CARTESIAN_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dy(ghosted_id)
              connection%area(iconn) = structured_grid%dx(ghosted_id)* &
                                   structured_grid%dz(ghosted_id)
              if (iface == SOUTH_FACE) then
                connection%dist(2,iconn) = 1.d0
              else
                connection%dist(2,iconn) = -1.d0
              endif
              if (associated(connection%id_dn2)) then
                  connection%id_dn2(iconn) = &
                    StructGetTVDGhostConnection(ghosted_id, &
                                                structured_grid, &
                                                iface,option)
              endif              
            case(CYLINDRICAL_GRID)
              print *, 'Cylindrical coordinates not applicable.'
              stop
            case(SPHERICAL_GRID)
              print *, 'Spherical coordinates not applicable.'
              stop
          end select

        case(BOTTOM_FACE,TOP_FACE)

          select case(structured_grid%itype)
            case(CARTESIAN_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dz(ghosted_id)
              connection%area(iconn) = structured_grid%dx(ghosted_id)* &
                                   structured_grid%dy(ghosted_id)
              if (structured_grid%invert_z_axis) then
                if (iface == TOP_FACE) then 
                  option%io_buffer = 'Need to ensure that direction of ' // &
                    'inverted z is correct in StructGridPopulateConnection()'
                  call printErrMsg(option)
                  connection%dist(3,iconn) = -1.d0
                else
                  connection%dist(3,iconn) = 1.d0
                endif
              else
                if (iface == BOTTOM_FACE) then 
                  connection%dist(3,iconn) = 1.d0
                else
                  connection%dist(3,iconn) = -1.d0
                endif
              endif
              if (associated(connection%id_dn2)) then
                  connection%id_dn2(iconn) = &
                    StructGetTVDGhostConnection(ghosted_id, &
                                                structured_grid, &
                                                iface,option)
              endif              
            case(CYLINDRICAL_GRID)
              connection%dist(:,iconn) = 0.d0
              connection%dist(0,iconn) = 0.5d0*structured_grid%dz(ghosted_id)
              connection%area(iconn) = 2.d0 * pi * radius(ghosted_id) * &
                                        structured_grid%dx(ghosted_id)
              if (structured_grid%invert_z_axis) then
                if (iface ==  TOP_FACE) then 
                  connection%dist(3,iconn) = 1.d0
                else
                  connection%dist(3,iconn) = -1.d0
                endif
              else
                if (iface ==  TOP_FACE) then 
                  connection%dist(3,iconn) = -1.d0
                else
                  connection%dist(3,iconn) = 1.d0
                endif
              endif
            case(SPHERICAL_GRID)
              option%io_buffer = &
                'Areas for spherical coordinates for z-axis not applicable.'
              call printErrMsg(option)
          end select
      end select
      if (connection%area(iconn) < 1.d-20) then
        write(option%io_buffer,*) connection%id_dn(iconn)
        option%io_buffer = &
          'Zero area in boundary connection at grid cell ' // &
          trim(adjustl(option%io_buffer)) // '.'
        call printErrMsg(option)
      endif
    case(INITIAL_CONNECTION_TYPE)
    case(SRC_SINK_CONNECTION_TYPE)
  end select
  
end subroutine StructGridPopulateConnection

! ************************************************************************** !

subroutine StructGridComputeVolumes(radius,structured_grid,option,nL2G,volume)
  ! 
  ! Computes the volumes of cells in structured grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/25/07
  ! 

#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  
  implicit none

  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: nL2G(:)
  PetscReal :: radius(:)
  Vec :: volume
  
  PetscReal, parameter :: Pi=3.141592653590d0
  
  PetscInt :: local_id, ghosted_id
  PetscReal, pointer :: volume_p(:)
  PetscReal :: r1, r2
  PetscErrorCode :: ierr
  
  call VecGetArrayF90(volume,volume_p, ierr);CHKERRQ(ierr)
  
  select case(structured_grid%itype)
    case(CARTESIAN_GRID)
      do local_id=1, structured_grid%nlmax
        ghosted_id = nL2G(local_id)
        volume_p(local_id) = structured_grid%dx(ghosted_id) * &
                             structured_grid%dy(ghosted_id) * &
                             structured_grid%dz(ghosted_id)
      enddo
    case(CYLINDRICAL_GRID)
      do local_id=1, structured_grid%nlmax
        ghosted_id = nL2G(local_id)
        ! volume = 2 pi r dr dz
        !       = pi * (r2-r1) * (r2+r1) * dz where dr = r2-r1 and 2 r = r2 + r1
        volume_p(local_id) = 2.d0 * pi * radius(ghosted_id) * &
                             structured_grid%dx(ghosted_id) * &
                             structured_grid%dz(ghosted_id)
      enddo
    case(SPHERICAL_GRID)
      do local_id=1, structured_grid%nlmax
        ghosted_id = nL2G(local_id)
        r2 = radius(ghosted_id) + 0.5d0*structured_grid%dx(ghosted_id)
        r1 = radius(ghosted_id) - 0.5d0*structured_grid%dx(ghosted_id)
        volume_p(local_id) = 4.d0/3.d0 * pi * structured_grid%dx(ghosted_id) &
                             * (r2*r2 + r2*r1 + r1*r1)
      enddo
  end select
  
  call VecRestoreArrayF90(volume,volume_p, ierr);CHKERRQ(ierr)
  
  if (option%print_to_screen .and. &
      option%mycommsize > 1 .and. option%mycommsize <= 16) then
    write(*,'(" rank= ",i3,", nlmax= ",i6,", nlx,y,z= ",3i4, &
      & ", nxs,e = ",2i4,", nys,e = ",2i4,", nzs,e = ",2i4)') &
      option%myrank,structured_grid%nlmax,structured_grid%nlx, &
        structured_grid%nly,structured_grid%nlz,structured_grid%lxs, &
        structured_grid%lxe,structured_grid%lys,structured_grid%lye, &
        structured_grid%lzs,structured_grid%lze

    write(*,'(" rank= ",i3,", ngmax= ",i6,", ngx,y,z= ",3i4, &
      & ", ngxs,e= ",2i4,", ngys,e= ",2i4,", ngzs,e= ",2i4)') &
      option%myrank,structured_grid%ngmax,structured_grid%ngx, &
        structured_grid%ngy,structured_grid%ngz,structured_grid%gxs, &
        structured_grid%gxe,structured_grid%gys,structured_grid%gye, &
        structured_grid%gzs,structured_grid%gze
  endif

end subroutine StructGridComputeVolumes

! ************************************************************************** !

subroutine StructGridMapIndices(structured_grid,stencil_type, &
                                nG2L,nL2G,nG2A,option)
  ! 
  ! maps global, local and natural indices of cells
  ! to each other
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Option_module

  implicit none

  type(grid_structured_type) :: structured_grid
  PetscEnum :: stencil_type
  PetscInt, pointer :: nG2L(:), nL2G(:), nG2A(:)
  type(option_type) :: option

  PetscInt :: i, j, k, local_id, ghosted_id, natural_id, count1
  PetscErrorCode :: ierr
  
  allocate(nL2G(structured_grid%nlmax))
  allocate(nG2L(structured_grid%ngmax))
  allocate(nG2A(structured_grid%ngmax))

  structured_grid%istart = structured_grid%lxs-structured_grid%gxs
  structured_grid%jstart = structured_grid%lys-structured_grid%gys
  structured_grid%kstart = structured_grid%lzs-structured_grid%gzs
  structured_grid%iend = structured_grid%istart+structured_grid%nlx-1
  structured_grid%jend = structured_grid%jstart+structured_grid%nly-1
  structured_grid%kend = structured_grid%kstart+structured_grid%nlz-1

  ! Local <-> Ghosted Transformation
  nG2L = 0  ! Must initialize this to zero!
  nL2G = 0

  nG2A = 0

  local_id = 0
  do k=structured_grid%kstart,structured_grid%kend
    do j=structured_grid%jstart,structured_grid%jend
      do i=structured_grid%istart,structured_grid%iend
        local_id = local_id + 1
        ghosted_id = i+j*structured_grid%ngx+k*structured_grid%ngxy+1
        nL2G(local_id) = ghosted_id
        nG2L(ghosted_id) = local_id
      enddo
    enddo
  enddo

  do i=1,structured_grid%ngmax
    j = nG2L(i)
    if (j > 0) then
      k = nL2G(j)
      if (i /= k) then
        print *,'Error in ghost-local node numbering for ghost node =', i
        print *,'node_id_gtol(i) =', j
        print *,'node_id_ltog(node_id_gtol(i)) =', k
        stop
      endif
    endif
  enddo
  ! Local(non ghosted)->Natural(natural order starts from 0)

  ! if STAR stencil, need to set corner ghosted cells to -1
  if (stencil_type == DMDA_STENCIL_STAR) then
    !geh - set corner ghosted nodes to -1
    do k=1,structured_grid%ngz
      do j=1,structured_grid%ngy
        do i=1,structured_grid%ngx
          count1 = 0
          if (i == 1 .and. &
              abs(structured_grid%lxs-structured_grid%gxs) > 0) &
            count1 = count1 + 1
          if (i == structured_grid%ngx .and. &
              abs(structured_grid%gxe-structured_grid%lxe) > 0) &
            count1 = count1 + 1
          if (j == 1 .and. &
              abs(structured_grid%lys-structured_grid%gys) > 0) &
            count1 = count1 + 1
          if (j == structured_grid%ngy .and. &
              abs(structured_grid%gye-structured_grid%lye) > 0) &
            count1 = count1 + 1
          if (k == 1 .and. &
              abs(structured_grid%lzs-structured_grid%gzs) > 0) &
            count1 = count1 + 1
          if (k == structured_grid%ngz .and. &
              abs(structured_grid%gze-structured_grid%lze) > 0) &
            count1 = count1 + 1
          if (count1 > 1) then
            ghosted_id = i+(j-1)*structured_grid%ngx+(k-1)*structured_grid%ngxy
            nG2L(ghosted_id) = -1
          endif
        enddo
      enddo
    enddo
  endif

  ! local ghosted -> natural (1-based)
  ghosted_id = 0
  do k=1,structured_grid%ngz
    do j=1,structured_grid%ngy
      do i=1,structured_grid%ngx
        ghosted_id = ghosted_id + 1
        natural_id = i + structured_grid%gxs + & ! 1-based
                      (j-1+structured_grid%gys)*structured_grid%nx+ &
                      (k-1+structured_grid%gzs)*structured_grid%nxy
        nG2A(ghosted_id) = natural_id
      enddo
    enddo
  enddo

end subroutine StructGridMapIndices

! ************************************************************************** !

subroutine StructGridGetGhostedNeighbors(structured_grid,ghosted_id, &
                                         stencil_type, &
                                         stencil_width_i,stencil_width_j, &
                                         stencil_width_k,x_count,y_count, &
                                         z_count,ghosted_neighbors, &
                                         option)
  ! 
  ! Returns an array of neighboring cells
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/11
  ! 
#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Option_module

  implicit none
  
  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscEnum :: stencil_type
  PetscInt :: stencil_width_i
  PetscInt :: stencil_width_j
  PetscInt :: stencil_width_k
  PetscInt :: x_count
  PetscInt :: y_count
  PetscInt :: z_count
  PetscInt :: ghosted_neighbors(*)

  PetscInt :: i, j, k
  PetscInt :: icount
  PetscInt :: ii, jj, kk

  call StructGridGetIJKFromGhostedID(structured_grid,ghosted_id,i,j,k)
  
  x_count = 0
  y_count = 0
  z_count = 0
  icount = 0
  select case(stencil_type)
    case(DMDA_STENCIL_STAR)
      do ii = max(i-stencil_width_i,1), min(i+stencil_width_i,structured_grid%ngx)
        if (ii /= i) then
          icount = icount + 1
          x_count = x_count + 1
          ghosted_neighbors(icount) = &
            StructGridGetGhostedIDFromIJK(structured_grid,ii,j,k)
        endif
      enddo
      do jj = max(j-stencil_width_j,1), min(j+stencil_width_j,structured_grid%ngy)
        if (jj /= j) then
          icount = icount + 1
          y_count = y_count + 1
          ghosted_neighbors(icount) = &
            StructGridGetGhostedIDFromIJK(structured_grid,i,jj,k)
        endif
      enddo
      do kk = max(k-stencil_width_k,1), min(k+stencil_width_k,structured_grid%ngz)
        if (kk /= k) then
          icount = icount + 1
          z_count = z_count + 1
          ghosted_neighbors(icount) = &
            StructGridGetGhostedIDFromIJK(structured_grid,i,j,kk)
        endif
      enddo
    case(DMDA_STENCIL_BOX)
      option%io_buffer = 'DMDA_STENCIL_BOX not yet supported in ' // &
        'StructGridGetNeighbors.'
      call printErrMsg(option)
  end select

end subroutine StructGridGetGhostedNeighbors

! ************************************************************************** !

subroutine StructGridGetGhostedNeighborsCorners(structured_grid,ghosted_id, &
                                         stencil_type, &
                                         stencil_width_i,stencil_width_j, &
                                         stencil_width_k, icount, &
                                         ghosted_neighbors, &
                                         option)
  ! 
  ! Returns an array of neighboring cells
  ! including the corner nodes
  ! Note that the previous subroutine does not return the corner nodes
  ! 
  ! Author: Satish Karra, LANL
  ! Date: 02/19/12
  ! 

  use Option_module

  implicit none
#include "petsc/finclude/petscdmda.h"
  
  type(grid_structured_type) :: structured_grid
  type(option_type) :: option
  PetscInt :: ghosted_id
  PetscEnum :: stencil_type
  PetscInt :: stencil_width_i
  PetscInt :: stencil_width_j
  PetscInt :: stencil_width_k
  PetscInt :: ghosted_neighbors(*)

  PetscInt :: i, j, k
  PetscInt :: icount
  PetscInt :: ii, jj, kk

  call StructGridGetIJKFromGhostedID(structured_grid,ghosted_id,i,j,k)

  icount = 0
  
  ! gb:08/08/13 Dependence on stencil_type is not necessary.
  !select case(stencil_type)
  !  case(DMDA_STENCIL_STAR)
      do kk = max(k-stencil_width_k,1), &
                min(k+stencil_width_k,structured_grid%ngz)
        do jj = max(j-stencil_width_j,1), &
                  min(j+stencil_width_j,structured_grid%ngy)
          do ii = max(i-stencil_width_i,1), &
                    min(i+stencil_width_i,structured_grid%ngx)
            if (ii == i .and. jj == j .and. kk == k) then
            ! do nothing
            else
              icount = icount + 1
              ghosted_neighbors(icount) = &
              StructGridGetGhostedIDFromIJK(structured_grid,ii,jj,kk)
            endif
          enddo
        enddo          
      enddo
  !  case(DMDA_STENCIL_BOX)
  !    option%io_buffer = 'DMDA_STENCIL_BOX not yet supported in ' // &
  !      'StructGridGetNeighbors.'
  !    call printErrMsg(option)
  !end select

end subroutine StructGridGetGhostedNeighborsCorners

! ************************************************************************** !

subroutine StructGridDestroy(structured_grid)
  ! 
  ! Deallocates a structured grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 
  use Utility_module, only : DeallocateArray
  
  implicit none
  
  type(grid_structured_type), pointer :: structured_grid
  
  PetscErrorCode :: ierr
    
  if (.not.associated(structured_grid)) return
  
  call DeallocateArray(structured_grid%dx_global)
  call DeallocateArray(structured_grid%dy_global)
  call DeallocateArray(structured_grid%dz_global)
  
  call DeallocateArray(structured_grid%dxg_local)
  call DeallocateArray(structured_grid%dyg_local)
  call DeallocateArray(structured_grid%dzg_local)
  
  call DeallocateArray(structured_grid%dx)
  call DeallocateArray(structured_grid%dy)
  call DeallocateArray(structured_grid%dz)
  
  call DeallocateArray(structured_grid%cell_neighbors)
  
  deallocate(structured_grid)
  nullify(structured_grid)

end subroutine StructGridDestroy

! ************************************************************************** !

subroutine StructGridCreateTVDGhosts(structured_grid,ndof,global_vec, &
                                     dm_1dof, &
                                     ghost_vec,scatter_ctx,option)
  ! 
  ! Calculates the TVD ghost vector and the
  ! associated scatter context
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/28/11
  ! 

#include "petsc/finclude/petscdm.h"
  use petscdm
  use Option_module

  implicit none

  type(grid_structured_type) :: structured_grid
  PetscInt :: ndof
  DM :: dm_1dof
   Vec :: global_vec
  Vec :: ghost_vec
  VecScatter :: scatter_ctx
  type(option_type) :: option

  PetscInt :: vector_size
  IS :: is_petsc
  IS :: is_ghost
  PetscInt :: icount, index, offset
  PetscInt :: increment
  PetscInt :: i, j, k
  PetscInt, allocatable :: global_indices_of_local_ghosted(:)
  PetscInt, allocatable :: global_indices_from(:)
  PetscInt, allocatable :: tvd_ghost_indices_to(:)
  ISLocalToGlobalMapping :: mapping_ltog
  PetscViewer :: viewer
  PetscErrorCode :: ierr
  
  ! structured grid has 6 sides to it
  vector_size = 0
  ! east-west
  if (structured_grid%nx > 1) then
    increment = structured_grid%nly*structured_grid%nlz
    vector_size = vector_size + 2*increment
  endif
  ! north-south
  if (structured_grid%ny > 1) then
    increment = structured_grid%nlx*structured_grid%nlz
    vector_size = vector_size + 2*increment
  endif
  ! top-bottom
  if (structured_grid%nz > 1) then
    increment = structured_grid%nlx*structured_grid%nly
    vector_size = vector_size + 2*increment
  endif
  
  if (vector_size == 0) then
    option%io_buffer = 'TVD does not handle a single grid cell.'
    call printErrMsg(option)
  endif
  
  call VecCreateSeq(PETSC_COMM_SELF,vector_size*ndof,ghost_vec, &
                    ierr);CHKERRQ(ierr)
  call VecSetBlockSize(ghost_vec,ndof,ierr);CHKERRQ(ierr)
  
  ! Create an IS composed of the petsc indexing of the ghost cells
  allocate(global_indices_from(vector_size))
  global_indices_from = UNINITIALIZED_INTEGER ! to catch bugs
  allocate(tvd_ghost_indices_to(vector_size))
  do i = 1, vector_size
    tvd_ghost_indices_to(i) = i-1
  enddo
  
  ! GEH: I'm going to play a trick here.  If I know the global index
  ! of my ghost cells, I can calculate the global index of the next
  ! layer of ghost cells in each direction since the block are 
  ! consistent through each dimension of the grid
  allocate(global_indices_of_local_ghosted(structured_grid%ngmax))
  do i = 1, structured_grid%ngmax
    global_indices_of_local_ghosted(i) = i-1
  enddo
  call DMGetLocalToGlobalMapping(dm_1dof,mapping_ltog,ierr);CHKERRQ(ierr)
  ! in and out integer arrays can be the same
  call ISLocalToGlobalMappingApply(mapping_ltog,structured_grid%ngmax, &
                                   global_indices_of_local_ghosted, &
                                   global_indices_of_local_ghosted, &
                                   ierr);CHKERRQ(ierr)
  ! leave global_indices_of_local_ghosted() in zero-based for the below
  
  ! Need to make a list of all indices that will receive updates through
  ! scatter/gather operation. Ghost cells representing physical boundaries
  ! do not need such an update.
  icount = 0

  if (structured_grid%nx > 1) then
    ! west
    offset = 0
    if (structured_grid%lxs /= structured_grid%gxs) offset = -1
    i = structured_grid%istart
    do k = structured_grid%kstart, structured_grid%kend
      do j = structured_grid%jstart, structured_grid%jend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo

    ! east
    offset = 0
    if (structured_grid%lxe /= structured_grid%gxe) offset = 1
    i = structured_grid%iend
    do k = structured_grid%kstart, structured_grid%kend
      do j = structured_grid%jstart, structured_grid%jend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  endif

  if (structured_grid%ny > 1) then
    ! south
    offset = 0
    if (structured_grid%lys /= structured_grid%gys) offset = -structured_grid%ngx
    j = structured_grid%jstart
    do k = structured_grid%kstart, structured_grid%kend
      do i = structured_grid%istart, structured_grid%iend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  
    ! north
    offset = 0
    if (structured_grid%lye /= structured_grid%gye) offset = structured_grid%ngx
    j = structured_grid%jend
    do k = structured_grid%kstart, structured_grid%kend
      do i = structured_grid%istart, structured_grid%iend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  endif
  
  if (structured_grid%nz > 1) then
    ! bottom
    offset = 0
    if (structured_grid%lzs /= structured_grid%gzs) offset = -structured_grid%ngxy
    k = structured_grid%kstart
    do j = structured_grid%jstart, structured_grid%jend
      do i = structured_grid%istart, structured_grid%iend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  
    ! top
    offset = 0
    if (structured_grid%lze /= structured_grid%gze) offset = structured_grid%ngxy
    k = structured_grid%kend
    do j = structured_grid%jstart, structured_grid%jend
      do i = structured_grid%istart, structured_grid%iend
        icount = icount + 1
        index = i + j*structured_grid%ngx + k*structured_grid%ngxy + 1
        global_indices_from(icount) = &
          global_indices_of_local_ghosted(index) + offset
      enddo
    enddo
  endif
  
  deallocate(global_indices_of_local_ghosted)

  if (vector_size /= icount) then
    option%io_buffer = 'Mis-count in TVD ghosting.'
    call printErrMsgByRank(option)
  endif

  ! since global_indices_from was base-zero, global_indices_from is base-zero.
  call ISCreateBlock(option%mycomm,ndof,vector_size, &
                      global_indices_from,PETSC_COPY_VALUES,is_petsc, &
                     ierr);CHKERRQ(ierr)
  deallocate(global_indices_from)

#if TVD_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_petsc_tvd.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_petsc,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  ! already zero-based
  call ISCreateBlock(option%mycomm,ndof,vector_size, &
                      tvd_ghost_indices_to,PETSC_COPY_VALUES,is_ghost, &
                     ierr);CHKERRQ(ierr)
  deallocate(tvd_ghost_indices_to)

#if TVD_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'is_ghost_tvd.out', &
                            viewer,ierr);CHKERRQ(ierr)
  call ISView(is_ghost,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  call VecScatterCreate(global_vec,is_petsc,ghost_vec,is_ghost, &
                        scatter_ctx,ierr);CHKERRQ(ierr)

#if TVD_DEBUG
  call PetscViewerASCIIOpen(option%mycomm,'tvd_ghost_scatter.out',viewer, &
                            ierr);CHKERRQ(ierr)
  call VecScatterView(scatter_ctx,viewer,ierr);CHKERRQ(ierr)
  call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

end subroutine StructGridCreateTVDGhosts  

! ************************************************************************** !

function StructGetTVDGhostConnection(ghosted_id,structured_grid,iface,option)
  ! 
  ! Returns id of tvd ghost cell for connection
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/10/11
  ! 

  use Option_module

  implicit none
  
  PetscInt :: ghosted_id
  type(grid_structured_type) :: structured_grid
  PetscInt :: iface
  type(option_type) :: option
  
  PetscInt :: StructGetTVDGhostConnection
  
  character(len=MAXSTRINGLENGTH) :: string
  PetscInt :: index
  PetscInt :: offset
  PetscInt :: i, j, k
  PetscBool :: error
  
  error = PETSC_FALSE

  select case(iface)
    case(WEST_FACE,EAST_FACE)
      if (structured_grid%ngx > 1) then
        if (iface == WEST_FACE) then
          StructGetTVDGhostConnection = ghosted_id + 1
        else
          StructGetTVDGhostConnection = ghosted_id - 1
        endif
        return
      elseif (structured_grid%nx == 1) then
        option%io_buffer = 'Boundary condition cannot be assigned in X ' // &
          'dimension with nx = 1 with TVD.'
        error = PETSC_TRUE
      endif
    case(SOUTH_FACE,NORTH_FACE)
      if (structured_grid%ngy > 1) then
        if (iface == SOUTH_FACE) then
          StructGetTVDGhostConnection = ghosted_id + structured_grid%ngx
        else
          StructGetTVDGhostConnection = ghosted_id - structured_grid%ngx
        endif
        return
      elseif (structured_grid%ny == 1) then
        option%io_buffer = 'Boundary condition cannot be assigned in Y ' // &
          'dimension with ny = 1 with TVD.'
        error = PETSC_TRUE
      endif
    case(BOTTOM_FACE,TOP_FACE)
      if (structured_grid%ngz > 1) then
        if (iface == BOTTOM_FACE) then
          StructGetTVDGhostConnection = ghosted_id + structured_grid%ngxy
        else
          StructGetTVDGhostConnection = ghosted_id - structured_grid%ngxy
        endif
        return
      elseif (structured_grid%nz == 1) then
        option%io_buffer = 'Boundary condition cannot be assigned in Z ' // &
          'dimension with nz = 1 with TVD.'
        error = PETSC_TRUE
      endif
  end select

  if (error) call printErrMsg(option)
  
  call StructGridGetIJKFromGhostedID(structured_grid,ghosted_id,i,j,k)
  offset = 0
  select case(iface)
    case(WEST_FACE)
      ! ensure that connection is on boundary face
      if (i /= 1) then
        error = PETSC_TRUE
        string = 'WEST'
      endif                       ! this must be in local dimension nly
      index = j + k*structured_grid%nly
    case(EAST_FACE)
      if (i /= structured_grid%ngx) then
        error = PETSC_TRUE
        string = 'EAST'
      endif
      index = j + k*structured_grid%nly + structured_grid%nlyz
    case(SOUTH_FACE)
      if (j /= 1) then
        error = PETSC_TRUE
        string = 'SOUTH'
      endif
      if (structured_grid%nx > 1) then
        offset = 2*structured_grid%nlyz
      endif
      index = i + k*structured_grid%nlx + offset
    case(NORTH_FACE)
      if (j /= structured_grid%ngy) then
        error = PETSC_TRUE
        string = 'NORTH'
      endif
      if (structured_grid%nx > 1) then
        offset = 2*structured_grid%nlyz
      endif
      index = i + k*structured_grid%nlx + structured_grid%nlxz + offset
    case(BOTTOM_FACE)
      if (k /= 1) then
        error = PETSC_TRUE
        string = 'BOTTOM'
      endif
      if (structured_grid%nx > 1) then
        offset = 2*structured_grid%nlyz
      endif
      if (structured_grid%ny > 1) then
        offset = offset + 2*structured_grid%nlxz
      endif
      index = i + j*structured_grid%nlx + offset
    case(TOP_FACE)
      if (k /= structured_grid%ngz) then
        error = PETSC_TRUE
        string = 'TOP'
      endif
      if (structured_grid%nx > 1) then
        offset = 2*structured_grid%nlyz
      endif
      if (structured_grid%ny > 1) then
        offset = offset + 2*structured_grid%nlxz
      endif
      index = i + j*structured_grid%nlx + structured_grid%nlxy + offset
  end select
  
  if (error) then
    write(option%io_buffer, '(''StructGetTVDGhostConnection not on '', a, &
    & ''face for cell:'',3i6)') trim(string), i,j,k
    call printErrMsgByRank(option)
  endif
  
  StructGetTVDGhostConnection = -index

end function StructGetTVDGhostConnection

end module Grid_Structured_module
