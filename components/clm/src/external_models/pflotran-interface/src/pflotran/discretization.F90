module Discretization_module
#include "petsc/finclude/petscdmda.h"
  use petscdmda
  use Grid_module
  use Grid_Structured_module
  use Grid_Unstructured_module
  use Grid_Unstructured_Aux_module
  use Grid_Unstructured_Explicit_module
  use Grid_Unstructured_Polyhedra_module
  use DM_Kludge_module

  use PFLOTRAN_Constants_module

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
#endif

  implicit none

  private

  type, public :: discretization_type
    PetscInt :: itype  ! type of discretization (e.g. structured, unstructured, etc.)
    !geh: note that differentiating between implicit and explicit unstructured 
    !     grids is handled within the grid%itype variable, not discritization%itype
    character(len=MAXWORDLENGTH) :: ctype
    PetscReal :: origin_global(3) ! origin of global domain
    type(grid_type), pointer :: grid  ! pointer to a grid object
    character(len=MAXSTRINGLENGTH) :: filename

    type(dm_ptr_type), pointer :: dmc_nflowdof(:), dmc_ntrandof(:)
      ! Arrays containing hierarchy of coarsened DMs, for use with Galerkin 
      ! multigrid.  Element i of each array is a *finer* DM than element i-1.
    PetscInt :: dm_index_to_ndof(5) ! mapping between a dm_ptr to the number of degrees of freedom
    type(dm_ptr_type), pointer :: dm_1dof
    type(dm_ptr_type), pointer :: dm_nflowdof
    type(dm_ptr_type), pointer :: dm_ntrandof
    type(dm_ptr_type), pointer :: dm_n_stress_strain_dof
    VecScatter :: tvd_ghost_scatter
    
    PetscInt :: stencil_width
    PetscEnum :: stencil_type
    
  end type discretization_type

  public :: DiscretizationCreate, &
            DiscretizationDestroy, &
            DiscretizationReadRequiredCards, &
            DiscretizationRead, &
            DiscretizationCreateVector, &
            DiscretizationDuplicateVector, &         
            DiscretizationCreateJacobian, &
            DiscretizationCreateInterpolation, &
            DiscretizationCreateColoring, &
            DiscretizationGlobalToLocal, &
            DiscretizationLocalToGlobal, &
            DiscretizationLocalToLocal, &
            DiscretizationGlobalToNatural, &
            DiscretizationNaturalToGlobal, &
            DiscretizationGlobalToLocalBegin, &
            DiscretizationGlobalToLocalEnd, &
            DiscretizationLocalToLocalBegin, &
            DiscretizationLocalToLocalEnd, &
            DiscretizGlobalToNaturalBegin, &
            DiscretizGlobalToNaturalEnd, &
            DiscretizNaturalToGlobalBegin, &
            DiscretizNaturalToGlobalEnd, &
            DiscretizationCreateDMs,&
            DiscretizationGetDMPtrFromIndex, &
            DiscretizationUpdateTVDGhosts, &
            DiscretAOApplicationToPetsc, &
            DiscretizationInputRecord, &
            DiscretizationPrintInfo
  
contains

! ************************************************************************** !

function DiscretizationCreate()
  ! 
  ! Creates a structured or unstructured discretization
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  implicit none
  
  type(discretization_type), pointer :: DiscretizationCreate
  
  type(discretization_type), pointer :: discretization
  
  allocate(discretization)
  discretization%ctype = ''
  discretization%itype = 0
  discretization%origin_global = 0.d0
  discretization%filename = ''

  ! nullify DM pointers
  nullify(discretization%dmc_nflowdof)
  nullify(discretization%dmc_ntrandof)
  allocate(discretization%dm_1dof)
  allocate(discretization%dm_nflowdof)
  allocate(discretization%dm_ntrandof)
  allocate(discretization%dm_n_stress_strain_dof)
  discretization%dm_1dof%dm = PETSC_NULL_DM
  discretization%dm_nflowdof%dm = PETSC_NULL_DM
  discretization%dm_ntrandof%dm = PETSC_NULL_DM
  discretization%dm_n_stress_strain_dof%dm = PETSC_NULL_DM
  nullify(discretization%dm_1dof%ugdm)
  nullify(discretization%dm_nflowdof%ugdm)
  nullify(discretization%dm_ntrandof%ugdm)
  nullify(discretization%dm_n_stress_strain_dof%ugdm)
  
  nullify(discretization%grid)
  
  discretization%stencil_width = 1
  discretization%stencil_type = DMDA_STENCIL_STAR

  discretization%tvd_ghost_scatter = PETSC_NULL_VECSCATTER
  
  DiscretizationCreate => discretization

end function DiscretizationCreate

! ************************************************************************** !

subroutine DiscretizationReadRequiredCards(discretization,input,option)
  ! 
  ! Reads a discretization from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(discretization_type),pointer :: discretization
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid, grid2
  type(grid_structured_type), pointer :: str_grid
  type(grid_unstructured_type), pointer :: un_str_grid
  character(len=MAXWORDLENGTH) :: structured_grid_ctype
  character(len=MAXWORDLENGTH) :: unstructured_grid_ctype

  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: structured_grid_itype
  PetscInt :: unstructured_grid_itype
  PetscInt :: nx, ny, nz
  PetscInt :: i
  PetscReal :: tempreal

  nx = 0
  ny = 0
  nz = 0

! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  do
  
    call InputReadPflotranString(input,option)
    if (input%ierr /= 0) exit

    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GRID')
    call StringToUpper(word)
      
    select case(trim(word))
      case('TYPE')
        call InputReadWord(input,option,discretization%ctype,PETSC_TRUE)
        call InputErrorMsg(input,option,'type','GRID')   
        call StringToLower(discretization%ctype)
        select case(trim(discretization%ctype))
          case('structured')
            discretization%itype = STRUCTURED_GRID
            call InputReadWord(input,option,structured_grid_ctype,PETSC_TRUE)
            call InputDefaultMsg(input,option,'grid_structured_type') 
            call StringToLower(structured_grid_ctype)
            select case(trim(structured_grid_ctype))
              case('cartesian')
                structured_grid_itype = CARTESIAN_GRID
              case('cylindrical')
                structured_grid_itype = CYLINDRICAL_GRID
              case('spherical')
                structured_grid_itype = SPHERICAL_GRID
              case default
                structured_grid_itype = CARTESIAN_GRID
                structured_grid_ctype = 'cartesian'
            end select
          case('unstructured','unstructured_explicit','unstructured_polyhedra')
            discretization%itype = UNSTRUCTURED_GRID
            word = discretization%ctype
            discretization%ctype = 'unstructured'
            select case(word)
              case('unstructured')
                unstructured_grid_itype = IMPLICIT_UNSTRUCTURED_GRID
                unstructured_grid_ctype = 'implicit unstructured'
              case('unstructured_explicit')
                unstructured_grid_itype = EXPLICIT_UNSTRUCTURED_GRID
                unstructured_grid_ctype = 'explicit unstructured'
              case('unstructured_polyhedra')
                unstructured_grid_itype = POLYHEDRA_UNSTRUCTURED_GRID
                unstructured_grid_ctype = 'polyhedra unstructured'
            end select
            call InputReadNChars(input,option,discretization%filename, &
                                 MAXSTRINGLENGTH,PETSC_TRUE)
            call InputErrorMsg(input,option,'unstructured filename','GRID')
          case default
            call InputKeywordUnrecognized(word,'discretization type',option)
        end select    
      case('NXYZ')
        call InputReadInt(input,option,nx)
        call InputErrorMsg(input,option,'nx','GRID')
        call InputReadInt(input,option,ny)
        call InputErrorMsg(input,option,'ny','GRID')
        call InputReadInt(input,option,nz)
        call InputErrorMsg(input,option,'nz','GRID')
        if (structured_grid_itype /= CARTESIAN_GRID) then
          ny = 1 ! cylindrical and spherical have 1 cell in Y
          if (structured_grid_itype /= CYLINDRICAL_GRID) nz = 1 ! spherical has 1 cell in Z
        endif
      case('ORIGIN')
        call InputReadDouble(input,option, &
                             discretization%origin_global(X_DIRECTION))
        call InputErrorMsg(input,option,'X direction','Origin')
        call InputReadDouble(input,option, &
                             discretization%origin_global(Y_DIRECTION))
        call InputErrorMsg(input,option,'Y direction','Origin')
        call InputReadDouble(input,option, &
                             discretization%origin_global(Z_DIRECTION))
        call InputErrorMsg(input,option,'Z direction','Origin')        
      case('FILE','GRAVITY','INVERT_Z','MAX_CELLS_SHARING_A_VERTEX',&
           'STENCIL_WIDTH','STENCIL_TYPE','FLUX_METHOD','DOMAIN_FILENAME', &
           'UPWIND_FRACTION_METHOD','PERM_TENSOR_TO_SCALAR_MODEL')
      case('DXYZ','BOUNDS')
        call InputSkipToEND(input,option,word) 
      case default
        call InputKeywordUnrecognized(word,'DISCRETIZATION',option)
    end select 
  enddo  

  if (discretization%itype == NULL_GRID) then
    option%io_buffer = 'Discretization type not defined under ' // &
                       'keyword GRID.' 
    call printErrMsg(option)
  endif
  
  grid => GridCreate()
  select case(discretization%itype)
    case(UNSTRUCTURED_GRID)
      un_str_grid => UGridCreate()
      select case(unstructured_grid_itype)
        case(IMPLICIT_UNSTRUCTURED_GRID)
          if (index(discretization%filename,'.h5') > 0) then
#if !defined(PETSC_HAVE_HDF5)
            option%io_buffer = 'PFLOTRAN must be built with HDF5 ' // &
              'support to read unstructured grid .h5 files'
            call printErrMsg(option)
#else

#ifdef SCORPIO
            call UGridReadHDF5PIOLib(un_str_grid,discretization%filename,option)
#else
            call UGridReadHDF5(un_str_grid,discretization%filename,option)
#endif
! #ifdef SCORPIO

#endif
!#if !defined(PETSC_HAVE_HDF5)

          else
            call UGridRead(un_str_grid,discretization%filename,option)
          endif
          grid%unstructured_grid => un_str_grid
        case(EXPLICIT_UNSTRUCTURED_GRID)
          un_str_grid%explicit_grid => UGridExplicitCreate()
          call UGridExplicitRead(un_str_grid, &
                                 discretization%filename,option)
          grid%unstructured_grid => un_str_grid
        case(POLYHEDRA_UNSTRUCTURED_GRID)
          un_str_grid%polyhedra_grid => UGridPolyhedraCreate()
          if (index(discretization%filename,'.h5') > 0 ) then
            !call UGridPolyhedraReadHDF5(un_str_grid,discretization%filename,option)
            call printErrMsg(option,'Add UGridPolyhedraReadHDF5')
          else
            call UGridPolyhedraRead(un_str_grid,discretization%filename,option)
          endif
          grid%unstructured_grid => un_str_grid
      end select
      grid%itype = unstructured_grid_itype
      grid%ctype = unstructured_grid_ctype
    case(STRUCTURED_GRID)      

#ifdef CLM_PFLOTRAN

      ! override readings from input cards above, if coupled with CLM BUT no mapping files provided
      if (.not.option%mapping_files) then
        !  nodes along X/Y direction
        if (clm_pf_idata%nxclm_mapped <= 0 .or.  &
            clm_pf_idata%nyclm_mapped <= 0) then
          call printErrMsg(option,'nxclm_mapped/nyclm_mapped NOT valid')

        else
          nx = clm_pf_idata%nxclm_mapped
          ny = clm_pf_idata%nyclm_mapped
        endif

        ! X/Y origins
        if (clm_pf_idata%x0clm_global == -9999.d0 .or.  &
            clm_pf_idata%y0clm_global == -9999.d0) then
          call printErrMsg(option,'x0clm_global/y0clm_global NOT valid')

        else
          discretization%origin_global(X_DIRECTION) = clm_pf_idata%x0clm_global
          discretization%origin_global(Y_DIRECTION) = clm_pf_idata%y0clm_global

        endif

      end if

      ! but always over-ride soil (vertical) discretization scheme
      if (clm_pf_idata%nyclm_mapped <= 0.or.  &
          clm_pf_idata%z0clm_global == -9999.d0) then
        call printErrMsg(option,'nx0clm_mapped/z0clm_global NOT valid')

      else
        nz = clm_pf_idata%nzclm_mapped
        discretization%origin_global(Z_DIRECTION) = clm_pf_idata%z0clm_global

      endif

#endif

      if (nx*ny*nz <= 0) &
        call printErrMsg(option,'NXYZ not set correctly for structured grid.')
      str_grid => StructGridCreate()
      str_grid%nx = nx
      str_grid%ny = ny
      str_grid%nz = nz
      str_grid%nxy = str_grid%nx*str_grid%ny
      str_grid%nmax = str_grid%nxy*str_grid%nz
      grid%structured_grid => str_grid
      grid%nmax = str_grid%nmax
      grid%structured_grid%itype = structured_grid_itype
      grid%structured_grid%ctype = structured_grid_ctype
      grid%itype = discretization%itype
      grid%ctype = discretization%ctype
  end select
  discretization%grid => grid
  nullify(grid)

end subroutine DiscretizationReadRequiredCards

! ************************************************************************** !

subroutine DiscretizationRead(discretization,input,option)
  ! 
  ! Reads a discretization from the input file
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Material_Aux_class

  implicit none

  type(option_type), pointer :: option
  type(input_type), pointer :: input
  type(discretization_type),pointer :: discretization
  character(len=MAXWORDLENGTH) :: word
  type(grid_type), pointer :: grid, grid2
  type(grid_structured_type), pointer :: str_grid
  type(grid_unstructured_type), pointer :: un_str_grid
  character(len=MAXWORDLENGTH) :: structured_grid_ctype
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: string

  PetscInt :: structured_grid_itype
  PetscInt :: nx, ny, nz
  PetscInt :: i
  PetscReal :: tempreal

  nx = 0
  ny = 0
  nz = 0

! we initialize the word to blanks to avoid error reported by valgrind
  word = ''

  do
  
    call InputReadPflotranString(input,option)
    if (input%ierr /= 0) exit

    if (InputCheckExit(input,option)) exit

    call InputReadWord(input,option,word,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GRID')
    call StringToUpper(word)
      
    select case(trim(word))
      case('TYPE','NXYZ','ORIGIN','FILE')
      case('DXYZ')
        select case(discretization%itype)
          case(STRUCTURED_GRID)

#ifdef CLM_PFLOTRAN

            !don't read input cards of PF grid, if coupled with CLM but no mapping files provided
            if (.not.option%mapping_files) then

              if (.not.associated(clm_pf_idata%dxclm_global) .or. &
                  .not.associated(clm_pf_idata%dyclm_global) .or. &
                   minval(clm_pf_idata%dxclm_global)<=1.d-6  .or. &
                   minval(clm_pf_idata%dyclm_global)<=1.d-6 )    then
                call printErrMsg(option,'dxclm_global or dyclm_global NOT valid FOR using CLM provided grid')
              endif

              allocate(discretization%grid%structured_grid%dx_global &
                (discretization%grid%structured_grid%nx))

              discretization%grid%structured_grid%dx_global = &
                clm_pf_idata%dxclm_global                               ! unit: longitudal degrees

              allocate(discretization%grid%structured_grid%dy_global &
                (discretization%grid%structured_grid%ny))
              discretization%grid%structured_grid%dy_global = &
                clm_pf_idata%dyclm_global                               ! unit: latitudal degrees

            else  ! the following IS needed, if providing mapping files
              call StructGridReadDXYZ(discretization%grid%structured_grid,input,option)

            endif

            ! always over-riding z-thickness
            if (.not.associated(discretization%grid%structured_grid%dz_global)) then
              allocate(discretization%grid%structured_grid%dz_global &
                 (discretization%grid%structured_grid%nz))
            endif
            discretization%grid%structured_grid%dz_global = &
              clm_pf_idata%dzclm_global                               ! unit: vertical meters

#else
            call StructGridReadDXYZ(discretization%grid%structured_grid,input,option)
#endif

          case default
            call printErrMsg(option,&
                           'Keyword "DXYZ" not supported for unstructured grid')
        end select

#ifdef CLM_PFLOTRAN
        call InputSkipToEND(input,option,'')    ! skip 'GRID/DXYZ ... END' block
#else
        call InputReadPflotranString(input,option) ! read END card
        call InputReadStringErrorMsg(input,option,'DISCRETIZATION,DXYZ,END')
        if (.not.(InputCheckExit(input,option))) then
          option%io_buffer = 'Card DXYZ should include either 3 entires ' // &
                   '(one for each grid direction or NX+NY+NZ entries)'
          call printErrMsg(option)
        endif
#endif

      case('BOUNDS')
        select case(discretization%itype)
          case(STRUCTURED_GRID)
            grid => discretization%grid

            ! read first line and we will split off the legacy approach vs. new
            call InputReadPflotranString(input,option)
            call InputReadStringErrorMsg(input,option, &
                                       'DISCRETIZATION,BOUNDS,Min Coordinates')
            select case(grid%structured_grid%itype)
              case(CARTESIAN_GRID)
                i = 3
              case(CYLINDRICAL_GRID)
                i = 2
              case(SPHERICAL_GRID)
                i = 1
            end select
            call InputReadNDoubles(input,option, &
                                   grid%structured_grid%bounds(:,LOWER), &
                                   i)
            call InputErrorMsg(input,option,'Minimum Coordinate','BOUNDS')
            call InputReadPflotranString(input,option)
            call InputReadStringErrorMsg(input,option, &
                                        'DISCRETIZATION,BOUNDS,Min Coordinates')
            call InputReadNDoubles(input,option, &
                                   grid%structured_grid%bounds(:,UPPER), &
                                   i)
            call InputErrorMsg(input,option,'Maximum Coordinate','BOUNDS')
            if (grid%structured_grid%itype == CYLINDRICAL_GRID) then
              ! 2 values were read in in x and y locations, must move y value
              ! to z as it was really z.
              grid%structured_grid%bounds(Z_DIRECTION,LOWER) = &
                grid%structured_grid%bounds(Y_DIRECTION,LOWER)
              grid%structured_grid%bounds(Z_DIRECTION,UPPER) = &
                grid%structured_grid%bounds(Y_DIRECTION,UPPER)
              ! set y bounds to 0 and 1
              grid%structured_grid%bounds(Y_DIRECTION,LOWER) = 0.d0
              grid%structured_grid%bounds(Y_DIRECTION,UPPER) = 1.d0
            endif
            if (grid%structured_grid%itype == SPHERICAL_GRID) then
              grid%structured_grid%bounds(Y_DIRECTION,LOWER) = 0.d0
              grid%structured_grid%bounds(Y_DIRECTION,UPPER) = 1.d0
              grid%structured_grid%bounds(Z_DIRECTION,LOWER) = 0.d0
              grid%structured_grid%bounds(Z_DIRECTION,UPPER) = 1.d0
            endif
            call InputReadPflotranString(input,option)
            call InputReadStringErrorMsg(input,option, &
                                         'DISCRETIZATION,BOUNDS,END')
            if (.not.(InputCheckExit(input,option))) then
              if (OptionPrintToScreen(option)) then
                if (grid%structured_grid%itype == CARTESIAN_GRID) then
                  print *, 'BOUNDS card for a cartesian structured grid ' // &
                    'must include 4 lines.  I.e.'
                  print *, 'BOUNDS'
                  print *, '  x_min  y_min  z_min'
                  print *, '  x_max  y_max  z_max'
                  print *, 'END'
                else if (grid%structured_grid%itype == CYLINDRICAL_GRID) then
                  print *, 'BOUNDS card for a cylindrical structured grid ' // &
                    'must include 4 lines.  I.e.'
                  print *, 'BOUNDS'
                  print *, '  r_min  z_min'
                  print *, '  r_max  z_max'
                  print *, 'END'
                else if (grid%structured_grid%itype == SPHERICAL_GRID) then
                  print *, 'BOUNDS card for a spherical structured grid ' // &
                    'must include 4 lines.  I.e.'
                  print *, 'BOUNDS'
                  print *, '  r_min'
                  print *, '  r_max'
                  print *, 'END'
                endif
              endif
              stop
            endif            
            discretization%origin_global(X_DIRECTION) = &
              grid%structured_grid%bounds(X_DIRECTION,LOWER)
            discretization%origin_global(Y_DIRECTION) = &
              grid%structured_grid%bounds(Y_DIRECTION,LOWER)
            discretization%origin_global(Z_DIRECTION) = &
              grid%structured_grid%bounds(Z_DIRECTION,LOWER)
        end select
      case ('GRAVITY')
        call InputReadDouble(input,option,option%gravity(X_DIRECTION))
        call InputErrorMsg(input,option,'x-direction','GRAVITY')
        call InputReadDouble(input,option,option%gravity(Y_DIRECTION))
        call InputErrorMsg(input,option,'y-direction','GRAVITY')
        call InputReadDouble(input,option,option%gravity(Z_DIRECTION))
        call InputErrorMsg(input,option,'z-direction','GRAVITY')
        if (option%myrank == option%io_rank .and. &
            option%print_to_screen) &
          write(option%fid_out,'(/," *GRAV",/, &
            & "  gravity    = "," [m/s^2]",3x,1p3e12.4 &
            & )') option%gravity(1:3)
      case ('MAX_CELLS_SHARING_A_VERTEX')
        if (associated(discretization%grid%unstructured_grid)) then
          call InputReadInt(input,option,discretization%grid% &
                            unstructured_grid%max_cells_sharing_a_vertex)
          call InputErrorMsg(input,option,'max_cells_sharing_a_vertex', &
                             'GRID')
        endif          
      case ('INVERT_Z')
      case ('STENCIL_WIDTH')
        call InputReadInt(input,option,discretization%stencil_width)
        call InputErrorMsg(input,option,'stencil_width', &
                           'GRID')
      case ('STENCIL_TYPE')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'stencil type','GRID')
        call StringToUpper(word)
        select case(trim(word))
          case ('BOX')
            discretization%stencil_type = DMDA_STENCIL_BOX
          case ('STAR')
            discretization%stencil_type = DMDA_STENCIL_STAR
          case default
            call InputKeywordUnrecognized(word, &
                   'DISCRETIZATION,stencil type',option)
        end select
      case('DOMAIN_FILENAME')
        select case(discretization%grid%itype)
          case(EXPLICIT_UNSTRUCTURED_GRID)
            call InputReadWord(input,option,word,PETSC_TRUE)
            call InputErrorMsg(input,option,'DOMAIN_FILENAME','GRID')  
            discretization%grid%unstructured_grid% &
              explicit_grid%domain_filename = word
          case default
            option%io_buffer = 'DOMAIN_FILENAME only supported for explicit &
                               &unstructured grids.'
            call printErrMsg(option)
        end select
      case('UPWIND_FRACTION_METHOD')
        if (discretization%itype == STRUCTURED_GRID) then
          option%io_buffer = 'UPWIND_FRACTION_METHOD not supported for &
            &structured grids.'
          call printErrMsg(option)
        endif
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'UPWIND_FRACTION_METHOD','GRID')
        call StringToUpper(word)
        select case(word)
          case('FACE_CENTER_PROJECTION')
            discretization%grid%unstructured_grid%upwind_fraction_method = &
              UGRID_UPWIND_FRACTION_PT_PROJ
          case('CELL_VOLUME')
            discretization%grid%unstructured_grid%upwind_fraction_method = &
              UGRID_UPWIND_FRACTION_CELL_VOL
          case('ABSOLUTE_DISTANCE')
            discretization%grid%unstructured_grid%upwind_fraction_method = &
              UGRID_UPWIND_FRACTION_ABS_DIST
          case default
            call InputKeywordUnrecognized(word,'GRID,UPWIND_FRACTION_METHOD', &
                                          option)
        end select

      case('PERM_TENSOR_TO_SCALAR_MODEL')
        call InputReadWord(input,option,word,PETSC_TRUE)
        call InputErrorMsg(input,option,'PERM_TENSOR_TO_SCALAR_MODEL','GRID')
        call StringToUpper(word)
        select case(word)
          case('LINEAR')
            call MaterialAuxSetPermTensorModel(TENSOR_TO_SCALAR_LINEAR,option)
          case('QUADRATIC')
            call MaterialAuxSetPermTensorModel(TENSOR_TO_SCALAR_QUADRATIC,&
                                              option)
          case default
            call InputKeywordUnrecognized(word,'GRID, PERM_TENSOR_TO_SCALAR_MODEL', &
                                          option)
        end select

      case default
        call InputKeywordUnrecognized(word,'GRID',option)
    end select 
  enddo  

  select case(discretization%itype)
    case(STRUCTURED_GRID)
      if (discretization%grid%structured_grid%invert_z_axis) then
        option%gravity(Z_DIRECTION) = -1.d0*option%gravity(Z_DIRECTION)
      endif
  end select
  
end subroutine DiscretizationRead

! ************************************************************************** !

subroutine DiscretizationCreateDMs(discretization, o_nflowdof, o_ntrandof, &
                                    o_nphase, o_ngeomechdof, o_n_stress_strain_dof, option)

  ! 
  ! creates distributed, parallel meshes/grids
  ! If there are multiple degrees of freedom per grid cell, this will call
  ! DiscretizationCreateDM() multiple times to create the DMs corresponding
  ! to one degree of freedom grid cell and those corresponding to multiple
  ! degrees of freedom per cell.
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/08
  ! 
      
  use Option_module    
      
  implicit none
  
  type(discretization_type) :: discretization
  PetscInt, intent(in) :: o_nflowdof
  PetscInt, intent(in) :: o_ntrandof
  PetscInt, intent(in) :: o_nphase
  PetscInt, intent(in) :: o_ngeomechdof
  PetscInt, intent(in) :: o_n_stress_strain_dof
  type(option_type) :: option
      
  PetscInt :: ndof
  !PetscInt, parameter :: stencil_width = 1
  PetscErrorCode :: ierr
  PetscInt :: i
  type(grid_unstructured_type), pointer :: ugrid

  select case(discretization%itype)
    case(STRUCTURED_GRID)
      discretization%dm_index_to_ndof(ONEDOF) = 1
      discretization%dm_index_to_ndof(NPHASEDOF) = o_nphase
      discretization%dm_index_to_ndof(NFLOWDOF) = o_nflowdof
      discretization%dm_index_to_ndof(NTRANDOF) = o_ntrandof
    case(UNSTRUCTURED_GRID)

    
      select case(discretization%grid%itype)
        case(IMPLICIT_UNSTRUCTURED_GRID)
          ! petsc will call parmetis to calculate the graph/dual
#if !defined(PETSC_HAVE_PARMETIS)
            option%io_buffer = &
             'Must compile with Parmetis in order to use implicit unstructured grids.'
            call printErrMsg(option)
#endif
          call UGridDecompose(discretization%grid%unstructured_grid, &
                              option)
        case(EXPLICIT_UNSTRUCTURED_GRID)
#if !defined(PETSC_HAVE_PARMETIS) && !defined(PETSC_HAVE_PTSCOTCH)
            option%io_buffer = &
             'Must compile with either Parmetis or PTSCOTCH in order to use explicit unstructured grids.'
            call printErrMsg(option)
#endif
          ugrid => discretization%grid%unstructured_grid
          call UGridExplicitDecompose(ugrid,option)
        case(POLYHEDRA_UNSTRUCTURED_GRID)
#if !defined(PETSC_HAVE_PARMETIS)
            option%io_buffer = &
             'Must compile with Parmetis in order to use implicit unstructured grids.'
            call printErrMsg(option)
#endif
          ugrid => discretization%grid%unstructured_grid
          call UGridPolyhedraDecompose(ugrid,option)
      end select
  end select


  !-----------------------------------------------------------------------
  ! Generate the DM objects that will manage communication.
  !-----------------------------------------------------------------------
  ndof = 1
  call DiscretizationCreateDM(discretization,discretization%dm_1dof, &
                              ndof,discretization%stencil_width, &
                              discretization%stencil_type,option)
  
  if (o_nflowdof > 0) then
    ndof = o_nflowdof
    call DiscretizationCreateDM(discretization,discretization%dm_nflowdof, &
                                ndof,discretization%stencil_width, &
                                discretization%stencil_type,option)
  endif
  
  if (o_ntrandof > 0) then
    ndof = o_ntrandof
    call DiscretizationCreateDM(discretization,discretization%dm_ntrandof, &
                                ndof,discretization%stencil_width, &
                                discretization%stencil_type,option)
  endif

  if (o_ngeomechdof > 0) then
    ndof = o_n_stress_strain_dof
    call DiscretizationCreateDM(discretization,discretization%dm_n_stress_strain_dof, &
                                ndof,discretization%stencil_width, &
                                discretization%stencil_type,option)
  endif

  select case(discretization%itype)
    case(STRUCTURED_GRID)
      ! this function must be called to set up str_grid%lxs, etc.
      call StructGridComputeLocalBounds(discretization%grid%structured_grid, &
                                        discretization%dm_1dof%dm,option)    
      discretization%grid%nlmax = discretization%grid%structured_grid%nlmax
      discretization%grid%ngmax = discretization%grid%structured_grid%ngmax
      discretization%grid%global_offset = &
        discretization%grid%structured_grid%global_offset
    case(UNSTRUCTURED_GRID)
      discretization%grid%nmax = discretization%grid%unstructured_grid%nmax
      discretization%grid%nlmax = discretization%grid%unstructured_grid%nlmax
      discretization%grid%ngmax = discretization%grid%unstructured_grid%ngmax
      discretization%grid%global_offset = &
        discretization%grid%unstructured_grid%global_offset
  end select

end subroutine DiscretizationCreateDMs

! ************************************************************************** !

subroutine DiscretizationCreateDM(discretization,dm_ptr,ndof,stencil_width, &
                                  stencil_type,option)
  ! 
  ! creates a distributed, parallel mesh/grid
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/08
  ! 

  use Option_module
  
  implicit none
  
  type(discretization_type) :: discretization
  type(dm_ptr_type), pointer :: dm_ptr
  PetscInt :: ndof
  PetscInt :: stencil_width
  PetscEnum :: stencil_type
  type(option_type) :: option
  PetscErrorCode :: ierr

  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call StructGridCreateDM(discretization%grid%structured_grid, &
                              dm_ptr%dm,ndof,stencil_width,stencil_type,option)
    case(UNSTRUCTURED_GRID)
      call UGridCreateUGDMShell(discretization%grid%unstructured_grid, &
                           dm_ptr%dm,dm_ptr%ugdm,ndof,option)
  end select

end subroutine DiscretizationCreateDM

! ************************************************************************** !

subroutine DiscretizationCreateVector(discretization,dm_index,vector, &
                                      vector_type,option)
  ! 
  ! Creates a global PETSc vector
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 
  use Option_module                                      

  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  Vec :: vector
  PetscInt :: vector_type
  type(option_type) :: option
  PetscInt :: ndof
  PetscErrorCode :: ierr
  
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)

  select case (vector_type)
    case(GLOBAL)
      call DMCreateGlobalVector(dm_ptr%dm,vector,ierr);CHKERRQ(ierr)
    case(LOCAL)
      call DMCreateLocalVector(dm_ptr%dm,vector,ierr);CHKERRQ(ierr)
    case(NATURAL)
      select case(discretization%itype)
        case(STRUCTURED_GRID)
          call DMDACreateNaturalVector(dm_ptr%dm,vector,ierr);CHKERRQ(ierr)
        case(UNSTRUCTURED_GRID)
          call UGridDMCreateVector(discretization%grid%unstructured_grid, &
                                   dm_ptr%ugdm,vector, &
                                   vector_type,option)
        end select
  end select

  call VecSet(vector,0.d0,ierr);CHKERRQ(ierr)
  
end subroutine DiscretizationCreateVector

! ************************************************************************** !

subroutine DiscretizationDuplicateVector(discretization,vector1,vector2)
  ! 
  ! Creates a global PETSc vector
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: vector1
  Vec :: vector2
  
  PetscErrorCode :: ierr
  call VecDuplicate(vector1,vector2,ierr);CHKERRQ(ierr)
  call VecCopy(vector1,vector2,ierr);CHKERRQ(ierr)
  
end subroutine DiscretizationDuplicateVector

! ************************************************************************** !

function DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  ! 
  ! Returns the integer pointer for the DM referenced
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/08/08
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  
  type(dm_ptr_type), pointer :: DiscretizationGetDMPtrFromIndex
  
  select case (dm_index)
    case(ONEDOF)
      DiscretizationGetDMPtrFromIndex => discretization%dm_1dof
    case(NFLOWDOF)
      DiscretizationGetDMPtrFromIndex => discretization%dm_nflowdof
    case(NTRANDOF)
      DiscretizationGetDMPtrFromIndex => discretization%dm_ntrandof
  end select  
  
end function DiscretizationGetDMPtrFromIndex

! ************************************************************************** !

function DiscretizationGetDMCPtrFromIndex(discretization,dm_index)

  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  
  type(dm_ptr_type), pointer :: DiscretizationGetDMCPtrFromIndex(:)
  
  select case (dm_index)
    case(NFLOWDOF)
      DiscretizationGetDMCPtrFromIndex => discretization%dmc_nflowdof
    case(NTRANDOF)
      DiscretizationGetDMCPtrFromIndex => discretization%dmc_ntrandof
  end select  
  
end function DiscretizationGetDMCPtrFromIndex

! ************************************************************************** !

subroutine DiscretizationCreateJacobian(discretization,dm_index,mat_type,Jacobian,option)
  ! 
  ! Creates Jacobian matrix associated with discretization
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  
  implicit none

  type(discretization_type) :: discretization
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  MatType :: mat_type
  Mat :: Jacobian
  type(option_type) :: option
  PetscInt :: ndof, stencilsize
  PetscInt, pointer :: indices(:)
  PetscInt :: ngmax
  PetscInt :: imax, nlevels, ln, npatches, pn, i
  type(dm_ptr_type), pointer :: dm_ptr
  ISLocalToGlobalMapping :: ptmap
  PetscInt :: islocal

  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)


  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DMSetMatType(dm_ptr%dm,mat_type,ierr);CHKERRQ(ierr)
      call DMCreateMatrix(dm_ptr%dm,Jacobian,ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
      call UGridDMCreateJacobian(discretization%grid%unstructured_grid, &
                                 dm_ptr%ugdm,mat_type,Jacobian,option)
  end select
  call MatSetOption(Jacobian,MAT_KEEP_NONZERO_PATTERN,PETSC_FALSE, &
                    ierr);CHKERRQ(ierr)
  call MatSetOption(Jacobian,MAT_ROW_ORIENTED,PETSC_FALSE,ierr);CHKERRQ(ierr)
  call MatSetOption(Jacobian,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE, &
                    ierr);CHKERRQ(ierr)

end subroutine DiscretizationCreateJacobian

! ************************************************************************** !

subroutine DiscretizationCreateInterpolation(discretization,dm_index, &
                                             interpolation,mg_levels_x, &
                                             mg_levels_y, mg_levels_z, &
                                             option)
  ! 
  ! Creates interpolation matrix associated
  ! with discretization for geometric multigrid.
  ! 
  ! Author: Richard Mills
  ! Date: 4/25/08.
  ! 

  use Option_module
  
  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  Mat, pointer :: interpolation(:)
  PetscInt :: mg_levels_x, mg_levels_y, mg_levels_z
  type(option_type) :: option

  PetscInt :: mg_levels
  PetscInt :: refine_x, refine_y, refine_z
  type(dm_ptr_type), pointer :: dm_ptr
  type(dm_ptr_type), pointer :: dmc_ptr(:)
  PetscInt :: i
  type(dm_ptr_type), pointer :: dm_fine_ptr
    ! Used to point to finer-grid DM in the loop that constructst the 
    ! interpolation hierarchy.
  
  mg_levels = max(mg_levels_x, mg_levels_y, mg_levels_z)
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
!  dmc_ptr = DiscretizationGetDMCPtrFromIndex(discretization,dm_index)
  select case (dm_index)
    case(NFLOWDOF)
      allocate(discretization%dmc_nflowdof(mg_levels))
      do i=1, mg_levels
        discretization%dmc_nflowdof(i)%dm = PETSC_NULL_DM
        nullify(discretization%dmc_nflowdof(i)%ugdm)
      enddo
      dmc_ptr => discretization%dmc_nflowdof
    case(NTRANDOF)
      allocate(discretization%dmc_ntrandof(mg_levels))
      do i=1, mg_levels
        discretization%dmc_ntrandof(i)%dm = PETSC_NULL_DM
        nullify(discretization%dmc_ntrandof(i)%ugdm)
      enddo
      dmc_ptr => discretization%dmc_ntrandof
  end select  
   
  allocate(interpolation(mg_levels))

  select case(discretization%itype)
    case(STRUCTURED_GRID)
      dm_fine_ptr => dm_ptr
      refine_x = 2; refine_y = 2; refine_z = 2
      do i=mg_levels-1,1,-1
        ! If number of coarsenings performed so far exceeds mg_levels_x-1, 
        ! set refine_x = 1; likewise for y and z.
        if (i <= mg_levels - mg_levels_x ) refine_x = 1
        if (i <= mg_levels - mg_levels_y ) refine_y = 1
        if (i <= mg_levels - mg_levels_z ) refine_z = 1
        call DMDASetRefinementFactor(dm_fine_ptr%dm, refine_x, refine_y, refine_z, &
                                   ierr);CHKERRQ(ierr)
        call DMDASetInterpolationType(dm_fine_ptr%dm, DMDA_Q0,  &
                                      ierr);CHKERRQ(ierr)
        call DMCoarsen(dm_fine_ptr%dm, option%mycomm, dmc_ptr(i)%dm,  &
                       ierr);CHKERRQ(ierr)
        call DMCreateInterpolation(dmc_ptr(i)%dm, dm_fine_ptr%dm, &
                                   interpolation(i), PETSC_NULL_VEC,  &
                                   ierr);CHKERRQ(ierr)
        dm_fine_ptr => dmc_ptr(i)
      enddo
    case(UNSTRUCTURED_GRID)
  end select

end subroutine DiscretizationCreateInterpolation

! ************************************************************************** !

subroutine DiscretizationCreateColoring(discretization,dm_index,option,coloring)
  ! 
  ! Creates ISColoring for discretization
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  
  implicit none
  
  type(discretization_type) :: discretization
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(option_type) :: option
  ISColoring :: coloring

  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DMSetMatType(dm_ptr%dm,MATBAIJ,ierr);CHKERRQ(ierr)
      call DMCreateColoring(dm_ptr%dm,IS_COLORING_GLOBAL,coloring,&
                            ierr);CHKERRQ(ierr)
      ! I have set the above to use matrix type MATBAIJ, as that is what we 
      ! usually want (note: for DAs with 1 degree of freedom per grid cell, 
      ! the MATAIJ and MATBAIJ colorings should be equivalent).  What we should 
      ! eventually do here is query the type of the Jacobian matrix, but I'm 
      ! not sure of the best way to do this, as this is currently stashed in 
      ! the 'solver' object. --RTM
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizationCreateColoring

! ************************************************************************** !

subroutine DiscretizationGlobalToLocal(discretization,global_vec,local_vec,dm_index)
  ! 
  ! Performs global to local communication with DM
  ! Note that 'dm_index' should correspond to one of the macros defined
  ! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

  implicit none

  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
    
  call DMGlobalToLocalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec, &
                            ierr);CHKERRQ(ierr)
  call DMGlobalToLocalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec, &
                          ierr);CHKERRQ(ierr)
  
end subroutine DiscretizationGlobalToLocal

! ************************************************************************** !

subroutine DiscretizationLocalToGlobal(discretization,local_vec,global_vec,dm_index)
  ! 
  ! Performs local to global communication with DM
  ! Note that 'dm_index' should correspond to one of the macros defined
  ! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/02/08
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMLocalToGlobalBegin(dm_ptr%dm,local_vec,INSERT_VALUES,global_vec, &
                            ierr);CHKERRQ(ierr)
  call DMLocalToGlobalEnd(dm_ptr%dm,local_vec,INSERT_VALUES,global_vec, &
                          ierr);CHKERRQ(ierr)
 
end subroutine DiscretizationLocalToGlobal

! ************************************************************************** !

subroutine DiscretizationLocalToLocal(discretization,local_vec1,local_vec2,dm_index)
  ! 
  ! Performs local to local communication with DM
  ! Some clarification:
  ! A "local to local" operation, in PETSc parlance, refers to communicating
  ! values from a local ghosted vector (in which the ghost points are
  ! irrelevant) and putting those values directly into another ghosted local
  ! vector (in which those ghost points are set correctly).
  ! This uses the same communication pattern as a "global to local" operation,
  ! but a in a "global to local", the originating vector is a PETSc global
  ! vector, not a ghosted local vector.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMLocalToLocalBegin(dm_ptr%dm,local_vec1,INSERT_VALUES,local_vec2, &
                           ierr);CHKERRQ(ierr)
  call DMLocalToLocalEnd(dm_ptr%dm,local_vec1,INSERT_VALUES,local_vec2, &
                         ierr);CHKERRQ(ierr)
  
end subroutine DiscretizationLocalToLocal

! ************************************************************************** !

subroutine DiscretizationGlobalToNatural(discretization,global_vec,natural_vec,dm_index)
  ! 
  ! Performs global to natural communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DMDAGlobalToNaturalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,natural_vec, &
                                    ierr);CHKERRQ(ierr)
      call DMDAGlobalToNaturalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,natural_vec, &
                                  ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
      call VecScatterBegin(dm_ptr%ugdm%scatter_gton,global_vec,natural_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(dm_ptr%ugdm%scatter_gton,global_vec,natural_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  end select
  
end subroutine DiscretizationGlobalToNatural

! ************************************************************************** !

subroutine DiscretizationNaturalToGlobal(discretization,natural_vec,global_vec,dm_index)
  ! 
  ! Performs natural to global communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/08
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DMDANaturalToGlobalBegin(dm_ptr%dm,natural_vec,INSERT_VALUES,global_vec, &
                                    ierr);CHKERRQ(ierr)
      call DMDANaturalToGlobalEnd(dm_ptr%dm,natural_vec,INSERT_VALUES,global_vec, &
                                  ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
      call VecScatterBegin(dm_ptr%ugdm%scatter_gton,natural_vec,global_vec, &
                           INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
      call VecScatterEnd(dm_ptr%ugdm%scatter_gton,natural_vec,global_vec, &
                         INSERT_VALUES,SCATTER_REVERSE,ierr);CHKERRQ(ierr)
  end select
  
end subroutine DiscretizationNaturalToGlobal

! ************************************************************************** !

subroutine DiscretizationGlobalToLocalBegin(discretization,global_vec,local_vec,dm_index)
  ! 
  ! Begins global to local communication with DM
  ! Note that 'dm_index' should correspond to one of the macros defined
  ! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 



  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMGlobalToLocalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec, &
                            ierr);CHKERRQ(ierr)
  
end subroutine DiscretizationGlobalToLocalBegin

! ************************************************************************** !

subroutine DiscretizationGlobalToLocalEnd(discretization,global_vec,local_vec,dm_index)
  ! 
  ! Ends global to local communication with DM
  ! Note that 'dm_index' should correspond to one of the macros defined
  ! in 'definitions.h' such as ONEDOF, NPHASEDOF, etc.  --RTM
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

 

 implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: local_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMGlobalToLocalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,local_vec, &
                          ierr);CHKERRQ(ierr)
 
end subroutine DiscretizationGlobalToLocalEnd

! ************************************************************************** !

subroutine DiscretizationLocalToLocalBegin(discretization,local_vec1,local_vec2,dm_index)
  ! 
  ! Begins local to local communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMLocalToLocalBegin(dm_ptr%dm,local_vec1,INSERT_VALUES,local_vec2, &
                           ierr);CHKERRQ(ierr)

end subroutine DiscretizationLocalToLocalBegin

! ************************************************************************** !

subroutine DiscretizationLocalToLocalEnd(discretization,local_vec1,local_vec2,dm_index)
  ! 
  ! Ends local to local communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/14/07
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: local_vec1
  Vec :: local_vec2
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  call DMLocalToLocalEnd(dm_ptr%dm,local_vec1,INSERT_VALUES,local_vec2, &
                         ierr);CHKERRQ(ierr)

end subroutine DiscretizationLocalToLocalEnd

! ************************************************************************** !

subroutine DiscretizGlobalToNaturalBegin(discretization,global_vec,natural_vec,dm_index)
  ! 
  ! Begins global to natural communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DMDAGlobalToNaturalBegin(dm_ptr%dm,global_vec,INSERT_VALUES,natural_vec, &
                                    ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
      call VecScatterBegin(dm_ptr%ugdm%scatter_gton,global_vec,natural_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  end select
  
end subroutine DiscretizGlobalToNaturalBegin

! ************************************************************************** !

subroutine DiscretizGlobalToNaturalEnd(discretization,global_vec,natural_vec,dm_index)
  ! 
  ! Ends global to natural communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/24/07
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: natural_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DMDAGlobalToNaturalEnd(dm_ptr%dm,global_vec,INSERT_VALUES,natural_vec, &
                                  ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
      call VecScatterEnd(dm_ptr%ugdm%scatter_gton,global_vec,natural_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
  end select
  
end subroutine DiscretizGlobalToNaturalEnd

! ************************************************************************** !

subroutine DiscretizNaturalToGlobalBegin(discretization,natural_vec,global_vec,dm_index)
  ! 
  ! Begins natural to global communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/08
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DMDANaturalToGlobalBegin(dm_ptr%dm,natural_vec,INSERT_VALUES,global_vec, &
                                    ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizNaturalToGlobalBegin

! ************************************************************************** !

subroutine DiscretizNaturalToGlobalEnd(discretization,natural_vec,global_vec,dm_index)
  ! 
  ! Ends natural to global communication with DM
  ! 
  ! Author: Glenn Hammond
  ! Date: 01/12/08
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: natural_vec
  Vec :: global_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  type(dm_ptr_type), pointer :: dm_ptr
  
  dm_ptr => DiscretizationGetDMPtrFromIndex(discretization,dm_index)
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DMDANaturalToGlobalEnd(dm_ptr%dm,natural_vec,INSERT_VALUES,global_vec, &
                                  ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
  end select
  
end subroutine DiscretizNaturalToGlobalEnd

! ************************************************************************** !

subroutine DiscretizationUpdateTVDGhosts(discretization,global_vec, &
                                         tvd_ghost_vec)
  ! 
  ! Updates tvd extended ghost cell values
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/12
  ! 

  implicit none
  
  type(discretization_type) :: discretization
  Vec :: global_vec
  Vec :: tvd_ghost_vec
  PetscInt :: dm_index
  PetscErrorCode :: ierr
  
  call VecScatterBegin(discretization%tvd_ghost_scatter,global_vec, &
                       tvd_ghost_vec,INSERT_VALUES,SCATTER_FORWARD, &
                       ierr);CHKERRQ(ierr)
  call VecScatterEnd(discretization%tvd_ghost_scatter,global_vec, &
                       tvd_ghost_vec,INSERT_VALUES,SCATTER_FORWARD, &
                     ierr);CHKERRQ(ierr)
  
end subroutine DiscretizationUpdateTVDGhosts

! ************************************************************************** !

subroutine DiscretAOApplicationToPetsc(discretization,int_array)
  ! 
  ! Maps application ordering to petsc
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/12/12
  ! 
#include "petsc/finclude/petscdmda.h"
  use petscdmda
  implicit none
  

  
  type(discretization_type) :: discretization
  PetscInt :: int_array(:)
  PetscErrorCode :: ierr
  
  AO :: ao
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      call DMDAGetAO(discretization%dm_1dof%dm,ao,ierr);CHKERRQ(ierr)
    case(UNSTRUCTURED_GRID)
      ao = discretization%grid%unstructured_grid%ao_natural_to_petsc
  end select
  call AOApplicationToPetsc(ao,size(int_array),int_array,ierr);CHKERRQ(ierr)
  
end subroutine DiscretAOApplicationToPetsc

! **************************************************************************** !

subroutine DiscretizationInputRecord(discretization)
  ! 
  ! Prints ingested grid/discretization information
  ! 
  ! Author: Jenn Frederick
  ! Date: 03/30/2016
  ! 

  implicit none

  type(discretization_type), pointer :: discretization

  type(grid_type), pointer :: grid
  character(len=MAXWORDLENGTH) :: word, word1, word2
  PetscInt :: id = INPUT_RECORD_UNIT
  character(len=10) :: Format, iFormat
  
  Format = '(ES14.7)'
  iFormat = '(I10)'

  grid => discretization%grid

  write(id,'(a)') ' '
  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
       &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'GRID'
  write(id,'(a29)',advance='no') 'grid type: '
  select case(grid%itype)
    case(STRUCTURED_GRID)
      write(id,'(a)') trim(grid%ctype)
      write(id,'(a29)',advance='no') ': '
      write(id,'(a)') trim(grid%structured_grid%ctype)
      write(id,'(a29)',advance='no') 'number grid cells X: '
      write(word,iFormat) grid%structured_grid%nx
      write(id,'(a)') adjustl(trim(word)) 
      write(id,'(a29)',advance='no') 'number grid cells Y: '
      write(word,iFormat) grid%structured_grid%ny
      write(id,'(a)') adjustl(trim(word)) 
      write(id,'(a29)',advance='no') 'number grid cells Z: '
      write(word,iFormat) grid%structured_grid%nz
      write(id,'(a)') adjustl(trim(word)) 
      write(id,'(a29)',advance='no') 'delta-X (m): '
      write(id,'(1p10e12.4)') grid%structured_grid%dx_global
      write(id,'(a29)',advance='no') 'delta-Y (m): '
      write(id,'(1p10e12.4)') grid%structured_grid%dy_global
      write(id,'(a29)',advance='no') 'delta-Z (m): '
      write(id,'(1p10e12.4)') grid%structured_grid%dz_global
      write(id,'(a29)',advance='no') 'bounds X: '
      write(word1,Format) grid%structured_grid%bounds(X_DIRECTION,LOWER)
      write(word2,Format) grid%structured_grid%bounds(X_DIRECTION,UPPER)
      write(id,'(a)') adjustl(trim(word1)) // ' ,' // adjustl(trim(word2)) // ' m'
      write(id,'(a29)',advance='no') 'bounds Y: '
      write(word1,Format) grid%structured_grid%bounds(Y_DIRECTION,LOWER)
      write(word2,Format) grid%structured_grid%bounds(Y_DIRECTION,UPPER)
      write(id,'(a)') adjustl(trim(word1)) // ' ,' // adjustl(trim(word2)) // ' m'
      write(id,'(a29)',advance='no') 'bounds Z: '
      write(word1,Format) grid%structured_grid%bounds(Z_DIRECTION,LOWER)
      write(word2,Format) grid%structured_grid%bounds(Z_DIRECTION,UPPER)
      write(id,'(a)') adjustl(trim(word1)) // ' ,' // adjustl(trim(word2)) // ' m'
    case(EXPLICIT_UNSTRUCTURED_GRID,IMPLICIT_UNSTRUCTURED_GRID, &
         POLYHEDRA_UNSTRUCTURED_GRID)
      write(id,'(a)') trim(grid%ctype)
  end select

  write(id,'(a29)',advance='no') 'global origin: '
  write(word,Format) discretization%origin_global(X_DIRECTION)
  write(id,'(a)') '(x) ' // adjustl(trim(word)) // ' m'
  write(id,'(a29)',advance='no') ': '
  write(word,Format) discretization%origin_global(Y_DIRECTION)
  write(id,'(a)') '(y) ' // adjustl(trim(word)) // ' m'
  write(id,'(a29)',advance='no') ': '
  write(word,Format) discretization%origin_global(Z_DIRECTION)
  write(id,'(a)') '(z) ' // adjustl(trim(word)) // ' m'

end subroutine DiscretizationInputRecord

! ************************************************************************** !

subroutine DiscretizationPrintInfo(discretization,grid,option)
  ! 
  ! Deallocates a discretization
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 
  use Option_module
  use Grid_module
  
  implicit none
  
  type(discretization_type) :: discretization
  type(option_type) :: option
  
  type(grid_type), pointer :: grid
  
  grid => discretization%grid
  
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      if (OptionPrintToScreen(option)) then
        write(*,'(/," Requested processors and decomposition = ", &
                 & i5,", npx,y,z= ",3i4)') &
            option%mycommsize,grid%structured_grid%npx, &
            grid%structured_grid%npy,grid%structured_grid%npz
        write(*,'(" Actual decomposition: npx,y,z= ",3i4,/)') &
            grid%structured_grid%npx_final,grid%structured_grid%npy_final, &
            grid%structured_grid%npz_final
      endif
      if (OptionPrintToFile(option)) then
        write(option%fid_out,'(/," Requested processors and decomposition = ", &
                             & i5,", npx,y,z= ",3i4)') &
            option%mycommsize,grid%structured_grid%npx,grid%structured_grid%npy, &
            grid%structured_grid%npz
        write(option%fid_out,'(" Actual decomposition: npx,y,z= ",3i4,/)') &
            grid%structured_grid%npx_final,grid%structured_grid%npy_final, &
            grid%structured_grid%npz_final
      endif
    case default
      if (OptionPrintToScreen(option)) then
        write(*,'(/," Requested processors = ",i5)') option%mycommsize
      endif
      if (OptionPrintToFile(option)) then
        write(option%fid_out,'(/," Requested processors = ",i5)') &
          option%mycommsize
      endif
  end select
  
end subroutine DiscretizationPrintInfo

! ************************************************************************** !

subroutine DiscretizationDestroy(discretization)
  ! 
  ! Deallocates a discretization
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/01/07
  ! 

  implicit none
  
  type(discretization_type), pointer :: discretization
  
  PetscErrorCode :: ierr
  PetscInt :: i
    
  if (.not.associated(discretization)) return
      
  select case(discretization%itype)
    case(STRUCTURED_GRID)
      if (discretization%dm_1dof%dm /= PETSC_NULL_DM) then
        call DMDestroy(discretization%dm_1dof%dm,ierr);CHKERRQ(ierr)
      endif
      discretization%dm_1dof%dm = PETSC_NULL_DM
      if (discretization%dm_nflowdof%dm /= PETSC_NULL_DM) then
        call DMDestroy(discretization%dm_nflowdof%dm,ierr);CHKERRQ(ierr)
      endif
      discretization%dm_nflowdof%dm = PETSC_NULL_DM
      if (discretization%dm_ntrandof%dm /= PETSC_NULL_DM) then
        call DMDestroy(discretization%dm_ntrandof%dm,ierr);CHKERRQ(ierr)
      endif
      discretization%dm_n_stress_strain_dof%dm = PETSC_NULL_DM
      if (discretization%dm_nflowdof%dm /= PETSC_NULL_DM) then
        call DMDestroy(discretization%dm_n_stress_strain_dof%dm, &
                       ierr);CHKERRQ(ierr)
      endif
      discretization%dm_n_stress_strain_dof%dm = PETSC_NULL_DM
      if (associated(discretization%dmc_nflowdof)) then
        do i=1,size(discretization%dmc_nflowdof)
          call DMDestroy(discretization%dmc_nflowdof(i)%dm,ierr);CHKERRQ(ierr)
        enddo
        deallocate(discretization%dmc_nflowdof)
        nullify(discretization%dmc_nflowdof)
      endif
      if (associated(discretization%dmc_ntrandof)) then
        do i=1,size(discretization%dmc_ntrandof)
          call DMDestroy(discretization%dmc_ntrandof(i)%dm,ierr);CHKERRQ(ierr)
        enddo
        deallocate(discretization%dmc_ntrandof)
        nullify(discretization%dmc_ntrandof)
      endif
    case(UNSTRUCTURED_GRID)
      if (associated(discretization%dm_1dof%ugdm)) &
        call UGridDMDestroy(discretization%dm_1dof%ugdm)
      if (associated(discretization%dm_nflowdof%ugdm)) &
        call UGridDMDestroy(discretization%dm_nflowdof%ugdm)
      if (associated(discretization%dm_ntrandof%ugdm)) &
        call UGridDMDestroy(discretization%dm_ntrandof%ugdm)
      if (associated(discretization%dm_n_stress_strain_dof%ugdm)) &
        call UGridDMDestroy(discretization%dm_n_stress_strain_dof%ugdm)

  end select
  if (associated(discretization%dm_1dof)) &
    deallocate(discretization%dm_1dof)
  nullify(discretization%dm_1dof)
  if (associated(discretization%dm_nflowdof)) &
    deallocate(discretization%dm_nflowdof)
  nullify(discretization%dm_nflowdof)
  if (associated(discretization%dm_ntrandof)) &
    deallocate(discretization%dm_ntrandof)
  nullify(discretization%dm_ntrandof)
  if (associated(discretization%dm_n_stress_strain_dof)) &
    deallocate(discretization%dm_n_stress_strain_dof)
  nullify(discretization%dm_n_stress_strain_dof)


  if (discretization%tvd_ghost_scatter /= PETSC_NULL_VECSCATTER) &
    call VecScatterDestroy(discretization%tvd_ghost_scatter)
  
  call GridDestroy(discretization%grid)
  
  deallocate(discretization)
  nullify(discretization)
  
end subroutine DiscretizationDestroy
 
end module Discretization_module
