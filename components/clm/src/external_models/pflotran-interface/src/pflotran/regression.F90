module Regression_module
 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Output_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private
 
  type, public :: regression_type
    type(regression_variable_type), pointer :: variable_list
    PetscInt, pointer :: natural_cell_ids(:)
    PetscInt :: num_cells_per_process
    PetscInt, pointer :: cells_per_process_natural_ids(:)
    PetscBool :: all_cells
    Vec :: natural_cell_id_vec
    Vec :: cells_per_process_vec
    VecScatter :: scatter_natural_cell_id_gtos
    VecScatter :: scatter_cells_per_process_gtos
    type(regression_type), pointer :: next
  end type regression_type

  type, public :: regression_variable_type
    character(len=MAXSTRINGLENGTH) :: name
    type(regression_variable_type), pointer :: next
  end type regression_variable_type
  
  public :: RegressionRead, &
            RegressionCreateMapping, &
            RegressionOutput, &
            RegressionDestroy
  
contains

! ************************************************************************** !

function RegressionCreate()
  ! 
  ! Creates a regression object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/11/12
  ! 
  
  implicit none

  type(regression_type), pointer :: RegressionCreate
  
  type(regression_type), pointer :: regression
  
  allocate(regression)
  nullify(regression%variable_list)
  nullify(regression%natural_cell_ids)
  regression%num_cells_per_process = 0
  regression%all_cells = PETSC_FALSE
  nullify(regression%cells_per_process_natural_ids)
  regression%natural_cell_id_vec = PETSC_NULL_VEC
  regression%cells_per_process_vec =  PETSC_NULL_VEC
  regression%scatter_natural_cell_id_gtos = PETSC_NULL_VECSCATTER
  regression%scatter_cells_per_process_gtos = PETSC_NULL_VECSCATTER
  nullify(regression%next)
  RegressionCreate => regression

end function RegressionCreate

! ************************************************************************** !

function RegressionVariableCreate()
  ! 
  ! Creates a regression variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/11/12
  ! 
  
  implicit none

  type(regression_variable_type), pointer :: RegressionVariableCreate
  
  type(regression_variable_type), pointer :: regression_variable
  
  allocate(regression_variable)
  regression_variable%name = ''
  nullify(regression_variable%next)
  RegressionVariableCreate => regression_variable

end function RegressionVariableCreate

! ************************************************************************** !

subroutine RegressionRead(regression,input,option)
  ! 
  ! Reads in contents of a regression card
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/11/12
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module

  implicit none
  
  type(regression_type), pointer :: regression
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word
  type(regression_variable_type), pointer :: cur_variable, new_variable
  PetscInt :: count, max_cells
  PetscInt, pointer :: int_array(:)
  PetscErrorCode :: ierr

  regression => RegressionCreate()
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','REGRESSION')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('VARIABLES') 
        count = 0
        do 
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'variable','REGRESSION,VARIABLES')
          call StringToUpper(word)
          new_variable => RegressionVariableCreate()
          new_variable%name = word
          if (.not.associated(regression%variable_list)) then
            regression%variable_list => new_variable
          else
            cur_variable%next => new_variable
          endif
          cur_variable => new_variable
        enddo
      case('CELLS')
        max_cells = 100
        allocate(int_array(max_cells))
        count = 0
        do 
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  
          count = count + 1
          if (count > max_cells) then
            call reallocateIntArray(int_array,max_cells)
          endif
          call InputReadInt(input,option,int_array(count))
          call InputErrorMsg(input,option,'natural cell id','REGRESSION,CELLS')
        enddo
        allocate(regression%natural_cell_ids(count))
        regression%natural_cell_ids = int_array(1:count)
        call PetscSortInt(count,regression%natural_cell_ids, &
                          ierr);CHKERRQ(ierr)
        deallocate(int_array)
      case('CELLS_PER_PROCESS')
        call InputReadInt(input,option,regression%num_cells_per_process)
        call InputErrorMsg(input,option,'num cells per process','REGRESSION')
      case('ALL_CELLS')
         regression%all_cells = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(keyword,'REGRESSION',option)
    end select
    
  enddo
  
end subroutine RegressionRead

! ************************************************************************** !

subroutine RegressionCreateMapping(regression,realization)
  ! 
  ! Creates mapping between a natural mpi vec and a
  ! sequential vec on io_rank
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/12/12
  ! 
#include "petsc/finclude/petscvec.h"
  use petscvec
  use Option_module
  use Realization_Subsurface_class
  use Grid_module
  use Discretization_module
  use Utility_module
  
  implicit none

  type(regression_type), pointer :: regression
  class(realization_subsurface_type) :: realization
  
  IS :: is_petsc
  PetscInt, allocatable :: int_array(:)
  PetscInt :: i, upper_bound, lower_bound, count, temp_int
  PetscInt :: local_id
  PetscReal, pointer :: vec_ptr(:)
  character(len=MAXWORDLENGTH) :: word
  Vec :: temp_vec
  VecScatter :: temp_scatter
  IS :: temp_is
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  if (.not.associated(regression)) return
  
  grid => realization%patch%grid
  option => realization%option

  if (regression%all_cells) then
    ! override regression%num_cells_per_process since cells will be duplicated
    regression%num_cells_per_process = 0 
    if (grid%nmax > 100) then
      option%io_buffer = 'Printing regression info for ALL_CELLS not &
        &supported for problem sizes greater than 100 cells.'
      call printErrMsg(option)
    endif
    call DeallocateArray(regression%natural_cell_ids)
    allocate(regression%natural_cell_ids(grid%nmax))
    do i = 1, grid%nmax
      regression%natural_cell_ids(i) = i
    enddo
  endif
  
  ! natural cell ids
  if (associated(regression%natural_cell_ids)) then
    ! ensure that natural ids are within problem domain
    if (maxval(regression%natural_cell_ids) > grid%nmax) then
      option%io_buffer = 'Natural IDs outside problem domain requested ' // &
        'for regression output.  Removing non-existent IDs.'
      call printWrnMsg(option)
      count = 0
      allocate(int_array(size(regression%natural_cell_ids)))
      int_array = 0
      do i = 1, size(regression%natural_cell_ids)
        if (regression%natural_cell_ids(i) <= grid%nmax) then
          count = count + 1
          int_array(count) = regression%natural_cell_ids(i)
        endif
      enddo
      ! reallocate array
      deallocate(regression%natural_cell_ids)
      allocate(regression%natural_cell_ids(count))
      !geh: Since natural_cell_ids and int_array may now be of different sizes,
      !     we need to be explicit about the values to copy.  gfortran has
      !     issues with this while Intel figures it out. Better to be explicit.
      regression%natural_cell_ids = int_array(1:count)
      deallocate(int_array)
    endif
    call VecCreate(PETSC_COMM_SELF,regression%natural_cell_id_vec, &
                   ierr);CHKERRQ(ierr)
    if (option%myrank == option%io_rank) then
      call VecSetSizes(regression%natural_cell_id_vec, &
                       size(regression%natural_cell_ids), &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
    else
      call VecSetSizes(regression%natural_cell_id_vec,0, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
    endif
    call VecSetFromOptions(regression%natural_cell_id_vec,ierr);CHKERRQ(ierr)
  
    if (option%myrank == option%io_rank) then
      count = size(regression%natural_cell_ids)
      ! determine how many of the natural cell ids are local
      allocate(int_array(count))
      int_array = regression%natural_cell_ids
      ! convert to zero based
      int_array = int_array - 1
    else
      count = 0
      allocate(int_array(count))
    endif
    call DiscretAOApplicationToPetsc(realization%discretization,int_array)
  
  ! create IS for global petsc cell ids
    call ISCreateGeneral(option%mycomm,count,int_array,PETSC_COPY_VALUES, &
                         is_petsc,ierr);CHKERRQ(ierr)
    deallocate(int_array)
  
#ifdef REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'is_petsc_natural_cell_id.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call ISView(is_petsc,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

    ! create scatter context
    call VecScatterCreate(realization%field%work,is_petsc, &
                          regression%natural_cell_id_vec,PETSC_NULL_IS, &
                          regression%scatter_natural_cell_id_gtos, &
                          ierr);CHKERRQ(ierr)

    call ISDestroy(is_petsc,ierr);CHKERRQ(ierr)

#ifdef REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'regression_scatter_nat_cell_ids.out',viewer, &
                              ierr);CHKERRQ(ierr)
    call VecScatterView(regression%scatter_natural_cell_id_gtos,viewer, &
                        ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  endif

  if (regression%num_cells_per_process > 0) then
    ! determine minimum number of cells per process
    i = grid%nlmax
    call MPI_Allreduce(i,count,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr)
    if (count < regression%num_cells_per_process) then
      write(word,*) count
      option%io_buffer = 'Number of cells per process for regression file&
        &exceeds minimum number of cells per process.  Truncating to ' // &
        trim(adjustl(word)) // '.'
      call printMsg(option)
      regression%num_cells_per_process = count
    endif
  
    ! cells ids per processor
    call VecCreate(PETSC_COMM_SELF,regression%cells_per_process_vec, &
                   ierr);CHKERRQ(ierr)
    if (option%myrank == option%io_rank) then
      call VecSetSizes(regression%cells_per_process_vec, &
                       regression%num_cells_per_process*option%mycommsize, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
    else
      call VecSetSizes(regression%cells_per_process_vec,ZERO_INTEGER, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
    endif
    call VecSetFromOptions(regression%cells_per_process_vec, &
                           ierr);CHKERRQ(ierr)

    ! create temporary vec to transfer down ids of cells
    call VecCreate(option%mycomm,temp_vec,ierr);CHKERRQ(ierr)
    call VecSetSizes(temp_vec, &
                     regression%num_cells_per_process, &
                     PETSC_DECIDE,ierr);CHKERRQ(ierr)
    call VecSetFromOptions(temp_vec,ierr);CHKERRQ(ierr)
  
    ! calculate interval
    call VecGetArrayF90(temp_vec,vec_ptr,ierr);CHKERRQ(ierr)
    temp_int = grid%nlmax / regression%num_cells_per_process
    do i = 1, regression%num_cells_per_process
      vec_ptr(i) = temp_int*(i-1) + 1 + grid%global_offset
    enddo
    call VecRestoreArrayF90(temp_vec,vec_ptr,ierr);CHKERRQ(ierr)

    ! create temporary scatter to transfer values to io_rank
    if (option%myrank == option%io_rank) then
      count = option%mycommsize*regression%num_cells_per_process
      ! determine how many of the natural cell ids are local
      allocate(int_array(count))
      do i = 1, count
        int_array(i) = i
      enddo
      ! convert to zero based
      int_array = int_array - 1
    else
      count = 0
      allocate(int_array(count))
    endif
    call ISCreateGeneral(option%mycomm,count, &
                         int_array,PETSC_COPY_VALUES,temp_is, &
                         ierr);CHKERRQ(ierr)

    call VecScatterCreate(temp_vec,temp_is, &
                          regression%cells_per_process_vec,PETSC_NULL_IS, &
                          temp_scatter,ierr);CHKERRQ(ierr)
    call ISDestroy(temp_is,ierr);CHKERRQ(ierr)
 
    ! scatter ids to io_rank
    call VecScatterBegin(temp_scatter,temp_vec, &
                         regression%cells_per_process_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(temp_scatter,temp_vec, &
                       regression%cells_per_process_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterDestroy(temp_scatter,ierr);CHKERRQ(ierr)
    call VecDestroy(temp_vec,ierr);CHKERRQ(ierr)
   
    ! transfer cell ids into array for creating new scatter
    if (option%myrank == option%io_rank) then
      count = option%mycommsize*regression%num_cells_per_process
      call VecGetArrayF90(regression%cells_per_process_vec,vec_ptr, &
                          ierr);CHKERRQ(ierr)
      do i = 1, count
        int_array(i) = int(vec_ptr(i)+0.1d0) ! tolerance to ensure int value
      enddo
      call VecRestoreArrayF90(regression%cells_per_process_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
      ! convert to zero based
      int_array = int_array - 1
    endif

    call ISCreateGeneral(option%mycomm,count, &
                         int_array,PETSC_COPY_VALUES,is_petsc, &
                         ierr);CHKERRQ(ierr)
    deallocate(int_array)

#ifdef REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'is_petsc_cells_per_process.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call ISView(is_petsc,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

    call VecScatterCreate(realization%field%work,is_petsc, &
                          regression%cells_per_process_vec, &
                          PETSC_NULL_IS, &
                          regression%scatter_cells_per_process_gtos, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_petsc,ierr);CHKERRQ(ierr)

#ifdef REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'regression_scatter_cells_per_process.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call VecScatterView(regression%scatter_cells_per_process_gtos,viewer, &
                        ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
    ! fill in natural ids of these cells on the io_rank
    if (option%myrank == option%io_rank) then
      allocate(regression%cells_per_process_natural_ids( &
               regression%num_cells_per_process*option%mycommsize))
    endif

    call VecGetArrayF90(realization%field%work,vec_ptr,ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax
      vec_ptr(local_id) = grid%nG2A(grid%nL2G(local_id))
    enddo
    call VecRestoreArrayF90(realization%field%work,vec_ptr,ierr);CHKERRQ(ierr)

    call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                          realization%field%work, &
                          regression%cells_per_process_vec, &
                          INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                       realization%field%work, &
                       regression%cells_per_process_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

    if (option%myrank == option%io_rank) then
      call VecGetArrayF90(regression%cells_per_process_vec,vec_ptr, &
                          ierr);CHKERRQ(ierr)
      regression%cells_per_process_natural_ids(:) = int(vec_ptr(:)+0.1)
      call VecRestoreArrayF90(regression%cells_per_process_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
    endif

  endif
  
end subroutine RegressionCreateMapping

! ************************************************************************** !

subroutine RegressionOutput(regression,realization,flow_timestepper, &
                            tran_timestepper)
  !
  ! Prints regression output through the io_rank
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/12/12
  ! 

  use Realization_Subsurface_class
  use Timestepper_BE_class
  use Option_module
  use Discretization_module
  use Output_module
  use Output_Aux_module
  use Output_Common_module, only : OutputGetCellCenteredVelocities, &
                                   OutputGetVariableArray
  
  implicit none
  
  type(regression_type), pointer :: regression
  class(realization_subsurface_type) :: realization
  ! these must be pointers as they can be null
  class(timestepper_BE_type), pointer :: flow_timestepper
  class(timestepper_BE_type), pointer :: tran_timestepper  
  
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  Vec :: global_vec_vx,global_vec_vy,global_vec_vz
  Vec :: x_vel_natural, y_vel_natural, z_vel_natural
  Vec :: x_vel_process, y_vel_process, z_vel_process
  type(option_type), pointer :: option
  type(output_variable_type), pointer :: cur_variable
  PetscReal, pointer :: vec_ptr(:), y_ptr(:), z_ptr(:)
  PetscInt :: i
  PetscInt :: iphase
  PetscReal :: r_norm, x_norm
  PetscReal :: max, min, mean
  PetscErrorCode :: ierr
  
  if (.not.associated(regression)) return
  
  option => realization%option
  
  if (option%myrank == option%io_rank) then
    string = trim(option%global_prefix) // &
             trim(option%group_prefix) // &  
             '.regression'
    option%io_buffer = '--> write regression output file: ' // trim(string)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=string,action="write")
  endif
  
  call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                  global_vec,GLOBAL,option)  
  call DiscretizationDuplicateVector(realization%discretization,global_vec,global_vec_vx)
  call DiscretizationDuplicateVector(realization%discretization,global_vec,global_vec_vy)
  call DiscretizationDuplicateVector(realization%discretization,global_vec,global_vec_vz)

  cur_variable => realization%output_option%output_snap_variable_list%first
  do 
    if (.not.associated(cur_variable)) exit
    
    call OutputGetVariableArray(realization,global_vec,cur_variable)
    
    call VecMax(global_vec,PETSC_NULL_INTEGER,max,ierr);CHKERRQ(ierr)
    call VecMin(global_vec,PETSC_NULL_INTEGER,min,ierr);CHKERRQ(ierr)
    call VecSum(global_vec,mean,ierr);CHKERRQ(ierr)
    mean = mean / realization%patch%grid%nmax
    
    ! list of natural ids
    if (associated(regression%natural_cell_ids)) then
      call VecScatterBegin(regression%scatter_natural_cell_id_gtos, &
                           global_vec, &
                           regression%natural_cell_id_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(regression%scatter_natural_cell_id_gtos, &
                         global_vec, &
                         regression%natural_cell_id_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    endif
    if (regression%num_cells_per_process > 0) then
      ! cells per process
      call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                           global_vec, &
                           regression%cells_per_process_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                         global_vec, &
                         regression%cells_per_process_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    endif

100 format(i9,': ',es21.13)    
101 format(i9,': ',i9)    
    
    if (option%myrank == option%io_rank) then
      string = OutputVariableToCategoryString(cur_variable%icategory)
      write(OUTPUT_UNIT,'(''-- '',a,'': '',a,'' --'')') &
        trim(string), trim(cur_variable%name)
      
      ! max, min, mean
      if (cur_variable%iformat == 0) then
        write(OUTPUT_UNIT,'(6x,''Max: '',es21.13)') max
        write(OUTPUT_UNIT,'(6x,''Min: '',es21.13)') min
      else
        write(OUTPUT_UNIT,'(6x,''Max: '',i9)') int(max)
        write(OUTPUT_UNIT,'(6x,''Min: '',i9)') int(min)
      endif
      write(OUTPUT_UNIT,'(5x,''Mean: '',es21.13)') mean
      
      ! natural cell ids
      if (associated(regression%natural_cell_ids)) then
        if (size(regression%natural_cell_ids) > 0) then
          call VecGetArrayF90(regression%natural_cell_id_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
          if (cur_variable%iformat == 0) then
            do i = 1, size(regression%natural_cell_ids)
              write(OUTPUT_UNIT,100) &
                regression%natural_cell_ids(i),vec_ptr(i)
            enddo
          else
            do i = 1, size(regression%natural_cell_ids)
              write(OUTPUT_UNIT,101) &
                regression%natural_cell_ids(i),int(vec_ptr(i))
            enddo
          endif
          call VecRestoreArrayF90(regression%natural_cell_id_vec,vec_ptr, &
                                  ierr);CHKERRQ(ierr)
        endif
      endif
      
      ! cell ids per process
      if (regression%num_cells_per_process > 0) then
        call VecGetArrayF90(regression%cells_per_process_vec,vec_ptr, &
                            ierr);CHKERRQ(ierr)
        if (cur_variable%iformat == 0) then
          do i = 1, regression%num_cells_per_process*option%mycommsize
            write(OUTPUT_UNIT,100) &
              regression%cells_per_process_natural_ids(i),vec_ptr(i)
          enddo
        else
          do i = 1, regression%num_cells_per_process*option%mycommsize
            write(OUTPUT_UNIT,101) &
              regression%cells_per_process_natural_ids(i),int(vec_ptr(i))
          enddo
        endif
        call VecRestoreArrayF90(regression%cells_per_process_vec,vec_ptr, &
                                ierr);CHKERRQ(ierr)
      endif
    endif
  
    cur_variable => cur_variable%next
  enddo
  
  ! velocities
  if ((realization%output_option%print_tecplot_vel_cent .or. &
       realization%output_option%print_hdf5_vel_cent) .and. &
      option%nflowdof > 0) then
    if (associated(regression%natural_cell_ids)) then
      call VecDuplicate(regression%natural_cell_id_vec,x_vel_natural, &
                        ierr);CHKERRQ(ierr)
      call VecDuplicate(regression%natural_cell_id_vec,y_vel_natural, &
                        ierr);CHKERRQ(ierr)
      call VecDuplicate(regression%natural_cell_id_vec,z_vel_natural, &
                        ierr);CHKERRQ(ierr)
      call VecZeroEntries(x_vel_natural,ierr);CHKERRQ(ierr)
      call VecZeroEntries(y_vel_natural,ierr);CHKERRQ(ierr)
      call VecZeroEntries(z_vel_natural,ierr);CHKERRQ(ierr)
    endif
    if (regression%num_cells_per_process > 0) then
      call VecDuplicate(regression%cells_per_process_vec,x_vel_process, &
                        ierr);CHKERRQ(ierr)
      call VecDuplicate(regression%cells_per_process_vec,y_vel_process, &
                        ierr);CHKERRQ(ierr)
      call VecDuplicate(regression%cells_per_process_vec,z_vel_process, &
                        ierr);CHKERRQ(ierr)
      call VecZeroEntries(x_vel_process,ierr);CHKERRQ(ierr)
      call VecZeroEntries(y_vel_process,ierr);CHKERRQ(ierr)
      call VecZeroEntries(z_vel_process,ierr);CHKERRQ(ierr)
    endif

    do iphase = 1, option%nphase
      if (associated(regression%natural_cell_ids) .or. &
          regression%num_cells_per_process > 0) then
    
        if (iphase == 1) then
          string = 'LIQUID'
        else
          string = 'GAS'
        endif
        if (option%myrank == option%io_rank) then
          write(OUTPUT_UNIT,'(''-- GENERIC: '',a,'' VELOCITY ['',a, &
                              &''] --'')') &
            trim(string), 'm/' // trim(realization%output_option%tunit)
        endif
    
        ! X
        call OutputGetCellCenteredVelocities(realization,global_vec_vx, &
                                             global_vec_vy,global_vec_vz, &
                                             iphase)
        if (associated(regression%natural_cell_ids)) then
          call VecScatterBegin(regression%scatter_natural_cell_id_gtos, &
                               global_vec_vx,x_vel_natural,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr);CHKERRQ(ierr)
          call VecScatterEnd(regression%scatter_natural_cell_id_gtos, &
                             global_vec_vx,x_vel_natural,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        endif
        if (regression%num_cells_per_process > 0) then
          call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                               global_vec_vx,x_vel_process,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr);CHKERRQ(ierr)
          call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                             global_vec_vx,x_vel_process,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        endif
        ! Y
        if (associated(regression%natural_cell_ids)) then
          call VecScatterBegin(regression%scatter_natural_cell_id_gtos, &
                               global_vec_vy,y_vel_natural,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr);CHKERRQ(ierr)
          call VecScatterEnd(regression%scatter_natural_cell_id_gtos, &
                             global_vec_vy,y_vel_natural,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        endif
        if (regression%num_cells_per_process > 0) then
          call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                               global_vec_vy,y_vel_process,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr);CHKERRQ(ierr)
          call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                             global_vec_vy,y_vel_process,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        endif
        ! Z
        if (associated(regression%natural_cell_ids)) then
          call VecScatterBegin(regression%scatter_natural_cell_id_gtos, &
                               global_vec_vz,z_vel_natural,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr);CHKERRQ(ierr)
          call VecScatterEnd(regression%scatter_natural_cell_id_gtos, &
                             global_vec_vz,z_vel_natural,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        endif
        if (regression%num_cells_per_process > 0) then
          call VecScatterBegin(regression%scatter_cells_per_process_gtos, &
                               global_vec_vz,z_vel_process,INSERT_VALUES, &
                               SCATTER_FORWARD,ierr);CHKERRQ(ierr)
          call VecScatterEnd(regression%scatter_cells_per_process_gtos, &
                             global_vec_vz,z_vel_process,INSERT_VALUES, &
                             SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        endif
      
104 format(i9,': ',3es21.13) 

        ! natural cell ids
        if (option%myrank == option%io_rank) then
          if (associated(regression%natural_cell_ids)) then
            if (size(regression%natural_cell_ids) > 0) then
              call VecGetArrayF90(x_vel_natural,vec_ptr,ierr);CHKERRQ(ierr)
              call VecGetArrayF90(y_vel_natural,y_ptr,ierr);CHKERRQ(ierr)
              call VecGetArrayF90(z_vel_natural,z_ptr,ierr);CHKERRQ(ierr)
              do i = 1, size(regression%natural_cell_ids)
                write(OUTPUT_UNIT,104) &
                  regression%natural_cell_ids(i),vec_ptr(i),y_ptr(i),z_ptr(i)
              enddo
              call VecRestoreArrayF90(x_vel_natural,vec_ptr, &
                                      ierr);CHKERRQ(ierr)
              call VecRestoreArrayF90(y_vel_natural,y_ptr,ierr);CHKERRQ(ierr)
              call VecRestoreArrayF90(z_vel_natural,z_ptr,ierr);CHKERRQ(ierr)
            endif
          endif
      
          ! cell ids per process
          if (regression%num_cells_per_process > 0) then
            call VecGetArrayF90(x_vel_process,vec_ptr,ierr);CHKERRQ(ierr)
            call VecGetArrayF90(y_vel_process,y_ptr,ierr);CHKERRQ(ierr)
            call VecGetArrayF90(z_vel_process,z_ptr,ierr);CHKERRQ(ierr)
            do i = 1, regression%num_cells_per_process*option%mycommsize
              write(OUTPUT_UNIT,104) &
                regression%cells_per_process_natural_ids(i),vec_ptr(i), &
                  y_ptr(i),z_ptr(i)
            enddo
            call VecRestoreArrayF90(x_vel_process,vec_ptr,ierr);CHKERRQ(ierr)
            call VecRestoreArrayF90(y_vel_process,y_ptr,ierr);CHKERRQ(ierr)
            call VecRestoreArrayF90(z_vel_process,z_ptr,ierr);CHKERRQ(ierr)
          endif
        endif
      endif
    enddo

    if (associated(regression%natural_cell_ids)) then
      call VecDestroy(x_vel_natural,ierr);CHKERRQ(ierr)
      call VecDestroy(y_vel_natural,ierr);CHKERRQ(ierr)
      call VecDestroy(z_vel_natural,ierr);CHKERRQ(ierr)
    endif
    if (regression%num_cells_per_process > 0) then
      call VecDestroy(x_vel_process,ierr);CHKERRQ(ierr)
      call VecDestroy(y_vel_process,ierr);CHKERRQ(ierr)
      call VecDestroy(z_vel_process,ierr);CHKERRQ(ierr)
    endif
  endif ! option%nflowdof > 0
  
  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vx,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vy,ierr);CHKERRQ(ierr)
  call VecDestroy(global_vec_vz,ierr);CHKERRQ(ierr)

  ! timestep, newton iteration, solver iteration output
  if (associated(flow_timestepper)) then
    call VecNorm(realization%field%flow_xx,NORM_2,x_norm,ierr);CHKERRQ(ierr)
    call VecNorm(realization%field%flow_r,NORM_2,r_norm,ierr);CHKERRQ(ierr)
    if (option%myrank == option%io_rank) then
      write(OUTPUT_UNIT,'(''-- SOLUTION: Flow --'')')
      write(OUTPUT_UNIT,'(''   Time (seconds): '',es21.13)') &
        flow_timestepper%cumulative_solver_time
      write(OUTPUT_UNIT,'(''   Time Steps: '',i12)') flow_timestepper%steps
      write(OUTPUT_UNIT,'(''   Newton Iterations: '',i12)') &
        flow_timestepper%cumulative_newton_iterations
      write(OUTPUT_UNIT,'(''   Solver Iterations: '',i12)') &
        flow_timestepper%cumulative_linear_iterations
      write(OUTPUT_UNIT,'(''   Time Step Cuts: '',i12)') &
        flow_timestepper%cumulative_time_step_cuts
      write(OUTPUT_UNIT,'(''   Solution 2-Norm: '',es21.13)') x_norm
      write(OUTPUT_UNIT,'(''   Residual 2-Norm: '',es21.13)') r_norm
    endif
  endif
  if (associated(tran_timestepper)) then
    call VecNorm(realization%field%tran_xx,NORM_2,x_norm,ierr);CHKERRQ(ierr)
    if (option%transport%reactive_transport_coupling == GLOBAL_IMPLICIT) then
      call VecNorm(realization%field%tran_r,NORM_2,r_norm,ierr);CHKERRQ(ierr)
    endif
    if (option%myrank == option%io_rank) then
      write(OUTPUT_UNIT,'(''-- SOLUTION: Transport --'')')
      write(OUTPUT_UNIT,'(''   Time (seconds): '',es21.13)') &
        tran_timestepper%cumulative_solver_time
      write(OUTPUT_UNIT,'(''   Time Steps: '',i12)') tran_timestepper%steps
      write(OUTPUT_UNIT,'(''   Newton Iterations: '',i12)') &
        tran_timestepper%cumulative_newton_iterations
      write(OUTPUT_UNIT,'(''   Solver Iterations: '',i12)') &
        tran_timestepper%cumulative_linear_iterations
      write(OUTPUT_UNIT,'(''   Time Step Cuts: '',i12)') &
        tran_timestepper%cumulative_time_step_cuts
      write(OUTPUT_UNIT,'(''   Solution 2-Norm: '',es21.13)') x_norm
      if (option%transport%reactive_transport_coupling == GLOBAL_IMPLICIT) then
        write(OUTPUT_UNIT,'(''   Residual 2-Norm: '',es21.13)') r_norm
      endif
    endif
  endif
  
  close(OUTPUT_UNIT)
  
end subroutine RegressionOutput

! ************************************************************************** !

recursive subroutine RegressionVariableDestroy(regression_variable)
  ! 
  ! Destroys a regression variable object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/11/12
  ! 

  implicit none
  
  type(regression_variable_type), pointer :: regression_variable
  
  if (.not.associated(regression_variable)) return
  
  call RegressionVariableDestroy(regression_variable%next)

  deallocate(regression_variable)
  nullify(regression_variable)
  
end subroutine RegressionVariableDestroy

! ************************************************************************** !

subroutine RegressionDestroy(regression)
  ! 
  ! Destroys a regression object
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/11/12
  ! 

  use Utility_module
  
  implicit none
  
  type(regression_type), pointer :: regression
  
  PetscErrorCode :: ierr
  
  if (.not.associated(regression)) return
  
  call RegressionVariableDestroy(regression%variable_list)
  call DeallocateArray(regression%natural_cell_ids)
  regression%num_cells_per_process = 0
  call DeallocateArray(regression%cells_per_process_natural_ids)
  if (regression%natural_cell_id_vec /= PETSC_NULL_VEC) then
    call VecDestroy(regression%natural_cell_id_vec,ierr);CHKERRQ(ierr)
  endif
  if (regression%cells_per_process_vec /= PETSC_NULL_VEC) then
    call VecDestroy(regression%cells_per_process_vec,ierr);CHKERRQ(ierr)
  endif
  if (regression%scatter_natural_cell_id_gtos /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(regression%scatter_natural_cell_id_gtos, &
                           ierr);CHKERRQ(ierr)
  endif
  if (regression%scatter_cells_per_process_gtos /= PETSC_NULL_VECSCATTER) then
    call VecScatterDestroy(regression%scatter_cells_per_process_gtos, &
                           ierr);CHKERRQ(ierr)
  endif

  deallocate(regression)
  nullify(regression)
  
end subroutine RegressionDestroy

end module Regression_module
