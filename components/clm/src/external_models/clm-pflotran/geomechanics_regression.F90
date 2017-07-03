module Geomechanics_Regression_module
 
  use Output_Aux_module
  
  use PFLOTRAN_Constants_module

  implicit none

  private

#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
 
  type, public :: geomechanics_regression_type
    type(geomechanics_regression_variable_type), pointer :: variable_list
    PetscInt, pointer :: natural_vertex_ids(:)
    PetscInt :: num_vertices_per_process
    PetscInt, pointer :: vertices_per_process_natural_ids(:)
    Vec :: natural_vertex_id_vec
    Vec :: vertices_per_process_vec
    VecScatter :: scatter_natural_vertex_id_gtos
    VecScatter :: scatter_vertices_per_process_gtos
    type(geomechanics_regression_type), pointer :: next
  end type geomechanics_regression_type

  type, public :: geomechanics_regression_variable_type
    character(len=MAXSTRINGLENGTH) :: name
    type(geomechanics_regression_variable_type), pointer :: next
  end type geomechanics_regression_variable_type
  
  public :: GeomechanicsRegressionRead, &
            GeomechanicsRegressionCreateMapping, &
            GeomechanicsRegressionOutput, &
            GeomechanicsRegressionDestroy
  
contains

! ************************************************************************** !

function GeomechanicsRegressionCreate()
  ! 
  ! Creates a geomechanics regression object
  ! 
  ! Author: Satish Karra
  ! Date: 06/01/2016
  ! 
  
  implicit none

  type(geomechanics_regression_type), pointer :: GeomechanicsRegressionCreate
  
  type(geomechanics_regression_type), pointer :: geomechanics_regression
  
  allocate(geomechanics_regression)
  nullify(geomechanics_regression%variable_list)
  nullify(geomechanics_regression%natural_vertex_ids)
  geomechanics_regression%num_vertices_per_process = 0
  nullify(geomechanics_regression%vertices_per_process_natural_ids)
  geomechanics_regression%natural_vertex_id_vec = 0
  geomechanics_regression%vertices_per_process_vec = 0
  geomechanics_regression%scatter_natural_vertex_id_gtos = 0
  geomechanics_regression%scatter_vertices_per_process_gtos = 0
  nullify(geomechanics_regression%next)
  GeomechanicsRegressionCreate => geomechanics_regression

end function GeomechanicsRegressionCreate

! ************************************************************************** !

function GeomechanicsRegressionVariableCreate()
  ! 
  ! Creates a geomechanics_regression variable object
  ! 
  ! Author: Satish Karra
  ! Date: 06/01/2016
  ! 
  
  implicit none

  type(geomechanics_regression_variable_type), pointer :: GeomechanicsRegressionVariableCreate
  
  type(geomechanics_regression_variable_type), pointer :: geomechanics_regression_variable
  
  allocate(geomechanics_regression_variable)
  geomechanics_regression_variable%name = ''
  nullify(geomechanics_regression_variable%next)
  GeomechanicsRegressionVariableCreate => geomechanics_regression_variable

end function GeomechanicsRegressionVariableCreate

! ************************************************************************** !

subroutine GeomechanicsRegressionRead(geomechanics_regression,input,option)
  ! 
  ! Reads in contents of a geomechanics_regression card
  ! 
  ! Author: Satish Karra
  ! Date: 06/22/2016
  ! 

  use Option_module
  use Input_Aux_module
  use String_module
  use Utility_module

  implicit none
  
  type(geomechanics_regression_type), pointer :: geomechanics_regression
  type(input_type), pointer :: input
  type(option_type) :: option
  
  character(len=MAXWORDLENGTH) :: keyword, word
  type(geomechanics_regression_variable_type), pointer :: cur_variable, new_variable
  PetscInt :: count, max_vertices
  PetscInt, pointer :: int_array(:)
  PetscErrorCode :: ierr

  geomechanics_regression => GeomechanicsRegressionCreate()
  
  input%ierr = 0
  do
  
    call InputReadPflotranString(input,option)

    if (InputCheckExit(input,option)) exit  

    call InputReadWord(input,option,keyword,PETSC_TRUE)
    call InputErrorMsg(input,option,'keyword','GEOMECHANICS_REGRESSION')
    call StringToUpper(keyword)   
      
    select case(trim(keyword))
    
      case('VARIABLES') 
        count = 0
        do 
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  

          call InputReadWord(input,option,word,PETSC_TRUE)
          call InputErrorMsg(input,option,'variable','GEOMECHANICS_REGRESSION,VARIABLES')
          call StringToLower(word)
          new_variable => GeomechanicsRegressionVariableCreate()
          new_variable%name = word
          if (.not.associated(geomechanics_regression%variable_list)) then
            geomechanics_regression%variable_list => new_variable
          else
            cur_variable%next => new_variable
          endif
          cur_variable => new_variable
        enddo
      case('VERTICES')
        max_vertices = 100
        allocate(int_array(max_vertices))
        count = 0
        do 
          call InputReadPflotranString(input,option)
          if (InputCheckExit(input,option)) exit  

          count = count + 1
          if (count > max_vertices) then
            call reallocateIntArray(int_array,max_vertices)
          endif
          call InputReadInt(input,option,int_array(count))
          call InputErrorMsg(input,option,'natural vertex id','GEOMECHANICS_REGRESSION,VERTICES')
        enddo
        allocate(geomechanics_regression%natural_vertex_ids(count))
        geomechanics_regression%natural_vertex_ids = int_array(1:count)
        call PetscSortInt(count,geomechanics_regression%natural_vertex_ids, &
                          ierr);CHKERRQ(ierr)
        deallocate(int_array)
      case('VERTICES_PER_PROCESS')
        call InputReadInt(input,option,geomechanics_regression%num_vertices_per_process)
        call InputErrorMsg(input,option,'num vertices per process','GEOMECHANICS_REGRESSION')
      case default
        call InputKeywordUnrecognized(keyword,'GEOMECHANICS_REGRESSION',option)
    end select
    
  enddo

end subroutine GeomechanicsRegressionRead

! ************************************************************************** !

subroutine GeomechanicsRegressionCreateMapping(geomechanics_regression, &
                                               geomechanics_realization)
  ! 
  ! Creates mapping between a natural mpi vec and a
  ! sequential vec on io_rank
  ! 
  ! Author: Satish Karra
  ! Date: 06/22/2016
  ! 

  use Option_module
  use Geomechanics_Realization_class
  use Geomechanics_Grid_Aux_module
  use Geomechanics_Discretization_module
  
  implicit none
  
#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscis.h90"
#include "petsc/finclude/petscviewer.h"

  type(geomechanics_regression_type), pointer :: geomechanics_regression
  class(realization_geomech_type) :: geomechanics_realization
  
  IS :: is_petsc
  PetscInt, allocatable :: int_array(:)
  PetscInt :: i, upper_bound, lower_bound, count, temp_int
  PetscInt :: local_id
  PetscReal, pointer :: vec_ptr(:)
  Vec :: temp_vec
  VecScatter :: temp_scatter
  IS :: temp_is
  PetscViewer :: viewer
  PetscErrorCode :: ierr

  type(geomech_grid_type), pointer :: grid
  type(option_type), pointer :: option
  
  if (.not.associated(geomechanics_regression)) return
  
  grid => geomechanics_realization%geomech_patch%geomech_grid
  option => geomechanics_realization%option
  
  ! natural vertex ids
  if (associated(geomechanics_regression%natural_vertex_ids)) then
    ! ensure that natural ids are within problem domain
    if (maxval(geomechanics_regression%natural_vertex_ids) > &
        grid%nmax_node) then
      option%io_buffer = 'Natural IDs outside geomechanics domain ' // &
        'requested for geomechanics regression output. ' // &
        'Removing non-existent IDs.'
      call printWrnMsg(option)
      count = 0
      allocate(int_array(size(geomechanics_regression%natural_vertex_ids)))
      int_array = 0
      do i = 1, size(geomechanics_regression%natural_vertex_ids)
        if (geomechanics_regression%natural_vertex_ids(i) & 
            <= grid%nmax_node) then
          count = count + 1
          int_array(count) = geomechanics_regression%natural_vertex_ids(i)
        endif
      enddo
      ! reallocate array
      deallocate(geomechanics_regression%natural_vertex_ids)
      allocate(geomechanics_regression%natural_vertex_ids(count))
      !geh: Since natural_vertex_ids and int_array may now be of different sizes,
      !     we need to be explicit about the values to copy.  gfortran has
      !     issues with this while Intel figures it out. Better to be explicit.
      geomechanics_regression%natural_vertex_ids = int_array(1:count)
      deallocate(int_array)
    endif
    ! Create a local vector on io_rank and scatter the nodal data
    ! to that local vector
    call VecCreate(PETSC_COMM_SELF, &
                   geomechanics_regression%natural_vertex_id_vec, &
                   ierr);CHKERRQ(ierr)
    if (option%myrank == option%io_rank) then
      call VecSetSizes(geomechanics_regression%natural_vertex_id_vec, &
                       size(geomechanics_regression%natural_vertex_ids), &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
    else
      call VecSetSizes(geomechanics_regression%natural_vertex_id_vec, &
                       ZERO_INTEGER, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
    endif
    call VecSetFromOptions(geomechanics_regression%natural_vertex_id_vec, &
                           ierr);CHKERRQ(ierr)
  
    if (option%myrank == option%io_rank) then
      count = size(geomechanics_regression%natural_vertex_ids)
      ! determine how many of the natural vertex ids are local
      allocate(int_array(count))
      int_array = geomechanics_regression%natural_vertex_ids
      ! convert to zero based
      int_array = int_array - 1
    else
      count = 0
      allocate(int_array(count))
    endif
    call GeomechDiscretAOApplicationToPetsc(geomechanics_realization% &
                                            geomech_discretization, &
                                            int_array)
  
  ! create IS for global petsc vertex ids
    call ISCreateGeneral(option%mycomm,count,int_array,PETSC_COPY_VALUES, &
                         is_petsc,ierr);CHKERRQ(ierr)
    deallocate(int_array)
  
#ifdef GEOMECHANICS_REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'geomech_is_petsc_natural_vertex_id.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call ISView(is_petsc,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

    ! create scatter context
    call VecScatterCreate(geomechanics_realization%geomech_field% &
                          press,is_petsc, &
                          geomechanics_regression%natural_vertex_id_vec, &
                          PETSC_NULL_OBJECT, &
                          geomechanics_regression% &
                          scatter_natural_vertex_id_gtos, &
                          ierr);CHKERRQ(ierr)

    call ISDestroy(is_petsc,ierr);CHKERRQ(ierr)

#ifdef GEOMECHANICS_REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                      'geomechanics_regression_scatter_nat_vertex_ids.out', &
                      viewer, &
                      ierr);CHKERRQ(ierr)
    call VecScatterView(geomechanics_regression% &
                        scatter_natural_vertex_id_gtos,viewer, &
                        ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

  endif

  if (geomechanics_regression%num_vertices_per_process > 0) then
    ! determine minimum number of vertices per process
    i = grid%nlmax_node
    call MPI_Allreduce(i,count,ONE_INTEGER_MPI,MPIU_INTEGER,MPI_MIN, &
                       option%mycomm,ierr)
    if (count < geomechanics_regression%num_vertices_per_process) then
      option%io_buffer = 'Number of vertices per process for ' // &
        'GeomechanicsRegression exceeds minimum number of vertices ' // & 
        'per process.  Truncating.'
      call printMsg(option)
      geomechanics_regression%num_vertices_per_process = count
    endif
  
    ! vertices ids per processor
    call VecCreate(PETSC_COMM_SELF, &
                   geomechanics_regression%vertices_per_process_vec, &
                   ierr);CHKERRQ(ierr)
    if (option%myrank == option%io_rank) then
      call VecSetSizes(geomechanics_regression%vertices_per_process_vec, &
                       geomechanics_regression%num_vertices_per_process* &
                       option%mycommsize, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
    else
      call VecSetSizes(geomechanics_regression% &
                       vertices_per_process_vec,ZERO_INTEGER, &
                       PETSC_DECIDE,ierr);CHKERRQ(ierr)
    endif
    call VecSetFromOptions(geomechanics_regression%vertices_per_process_vec, &
                           ierr);CHKERRQ(ierr)

    ! create temporary vec to transfer down ids of vertices
    call VecCreate(option%mycomm,temp_vec,ierr);CHKERRQ(ierr)
    call VecSetSizes(temp_vec, &
                     geomechanics_regression%num_vertices_per_process, &
                     PETSC_DECIDE,ierr);CHKERRQ(ierr)
    call VecSetFromOptions(temp_vec,ierr);CHKERRQ(ierr)
  
    ! calculate interval
    call VecGetArrayF90(temp_vec,vec_ptr,ierr);CHKERRQ(ierr)
    temp_int = grid%nlmax_node / &
               geomechanics_regression%num_vertices_per_process
    do i = 1, geomechanics_regression%num_vertices_per_process
      vec_ptr(i) = temp_int*(i-1) + 1 + grid%global_offset
    enddo
    call VecRestoreArrayF90(temp_vec,vec_ptr,ierr);CHKERRQ(ierr)

    ! create temporary scatter to transfer values to io_rank
    if (option%myrank == option%io_rank) then
      count = option%mycommsize* &
              geomechanics_regression%num_vertices_per_process
      ! determine how many of the natural vertex ids are local
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
                          geomechanics_regression%vertices_per_process_vec, &
                          PETSC_NULL_OBJECT, &
                          temp_scatter,ierr);CHKERRQ(ierr)
    call ISDestroy(temp_is,ierr);CHKERRQ(ierr)
 
    ! scatter ids to io_rank
    call VecScatterBegin(temp_scatter,temp_vec, &
                         geomechanics_regression%vertices_per_process_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(temp_scatter,temp_vec, &
                       geomechanics_regression%vertices_per_process_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterDestroy(temp_scatter,ierr);CHKERRQ(ierr)
    call VecDestroy(temp_vec,ierr);CHKERRQ(ierr)
   
    ! transfer vertex ids into array for creating new scatter
    if (option%myrank == option%io_rank) then
      count = option%mycommsize* &
              geomechanics_regression%num_vertices_per_process
      call VecGetArrayF90(geomechanics_regression%vertices_per_process_vec, &
                          vec_ptr,ierr);CHKERRQ(ierr)
      do i = 1, count
        int_array(i) = int(vec_ptr(i)+0.1d0) ! tolerance to ensure int value
      enddo
      call VecRestoreArrayF90(geomechanics_regression% &
                              vertices_per_process_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
      ! convert to zero based
      int_array = int_array - 1
    endif

    call ISCreateGeneral(option%mycomm,count, &
                         int_array,PETSC_COPY_VALUES,is_petsc, &
                         ierr);CHKERRQ(ierr)
    deallocate(int_array)

#ifdef GEOMECHANICS_REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                              'geomech_is_petsc_vertices_per_process.out', &
                              viewer,ierr);CHKERRQ(ierr)
    call ISView(is_petsc,viewer,ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif

    call VecScatterCreate(geomechanics_realization%geomech_field% &
                          press,is_petsc, &
                          geomechanics_regression%vertices_per_process_vec, &
                          PETSC_NULL_OBJECT, &
                          geomechanics_regression% &
                          scatter_vertices_per_process_gtos, &
                          ierr);CHKERRQ(ierr)
    call ISDestroy(is_petsc,ierr);CHKERRQ(ierr)

#ifdef GEOMECHANICS_REGRESSION_DEBUG
    call PetscViewerASCIIOpen(option%mycomm, &
                  'geomechanics_regression_scatter_vertices_per_process.out', &
                  viewer,ierr);CHKERRQ(ierr)
    call VecScatterView(geomechanics_regression% &
                        scatter_vertices_per_process_gtos,viewer, &
                        ierr);CHKERRQ(ierr)
    call PetscViewerDestroy(viewer,ierr);CHKERRQ(ierr)
#endif
  
    ! fill in natural ids of these vertices on the io_rank
    if (option%myrank == option%io_rank) then
      allocate(geomechanics_regression%vertices_per_process_natural_ids( &
               geomechanics_regression%num_vertices_per_process* &
               option%mycommsize))
    endif

    call VecGetArrayF90(geomechanics_realization%geomech_field%press,vec_ptr, &
                        ierr);CHKERRQ(ierr)
    do local_id = 1, grid%nlmax_node
      vec_ptr(local_id) = grid%nG2A(grid%nL2G(local_id))
    enddo
    call VecRestoreArrayF90(geomechanics_realization%geomech_field%press, &
                            vec_ptr,ierr);CHKERRQ(ierr)

    call VecScatterBegin(geomechanics_regression% &
                         scatter_vertices_per_process_gtos, &
                         geomechanics_realization%geomech_field%press, &
                         geomechanics_regression%vertices_per_process_vec, &
                         INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
    call VecScatterEnd(geomechanics_regression% &
                       scatter_vertices_per_process_gtos, &
                       geomechanics_realization%geomech_field%press, &
                       geomechanics_regression%vertices_per_process_vec, &
                       INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)

    if (option%myrank == option%io_rank) then
      call VecGetArrayF90(geomechanics_regression% &
                          vertices_per_process_vec,vec_ptr, &
                          ierr);CHKERRQ(ierr)
      geomechanics_regression%vertices_per_process_natural_ids(:) = &
                                                        int(vec_ptr(:)+0.1)
      call VecRestoreArrayF90(geomechanics_regression% &
                              vertices_per_process_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
    endif

  endif
  
end subroutine GeomechanicsRegressionCreateMapping

! ************************************************************************** !

subroutine GeomechanicsRegressionOutput(geomechanics_regression, &
                                        geomechanics_realization, &
                                        geomechanics_timestepper)
  !
  ! Prints geomechanics regression output through the io_rank
  ! 
  ! Author: Satish Karra
  ! Date: 06/22/2016
  ! 

  use Geomechanics_Realization_class
  use Timestepper_Steady_class
  use Option_module
  use Geomechanics_Discretization_module
  use Output_Geomechanics_module, only : OutputGeomechGetVarFromArray
  
  implicit none
  
  type(geomechanics_regression_type), pointer :: geomechanics_regression
  class(realization_geomech_type) :: geomechanics_realization
  class(timestepper_steady_type), pointer :: geomechanics_timestepper
  ! these must be pointers as they can be null
  character(len=MAXSTRINGLENGTH) :: string
  Vec :: global_vec
  PetscInt :: ivar, isubvar
  type(option_type), pointer :: option
  type(output_variable_type), pointer :: cur_variable
  type(geomechanics_regression_variable_type), pointer :: cur_variable1
  PetscReal, pointer :: vec_ptr(:), y_ptr(:), z_ptr(:)
  PetscInt :: i
  PetscInt :: iphase
  PetscReal :: r_norm, x_norm
  PetscReal :: max, min, mean
  PetscErrorCode :: ierr
  PetscBool :: found

  if (.not.associated(geomechanics_regression)) return
  
  option => geomechanics_realization%option
  
  if (option%myrank == option%io_rank) then
    string = trim(option%global_prefix) // &
             trim(option%group_prefix) // &  
             '.regression'
    option%io_buffer = '--> write geomechanics_regression output file: ' // trim(string)
    call printMsg(option)
    open(unit=OUTPUT_UNIT,file=string,action="write")
  endif
    
  call GeomechDiscretizationCreateVector(geomechanics_realization% &
                                         geomech_discretization,ONEDOF, &
                                         global_vec,GLOBAL,option)
  cur_variable => geomechanics_realization%output_option% &
                  output_snap_variable_list%first

  do
    if (.not.associated(cur_variable)) exit

    found = PETSC_FALSE

    cur_variable1 => geomechanics_regression%variable_list
    do
      if (.not.associated(cur_variable1)) exit
      if (cur_variable%name == cur_variable1%name) then
        found = PETSC_TRUE
        exit
       endif
      cur_variable1 => cur_variable1%next
    enddo

    if (found) then

      ivar = cur_variable%ivar
      isubvar = cur_variable%isubvar

      call OutputGeomechGetVarFromArray(geomechanics_realization,global_vec,ivar,isubvar)
      call VecMax(global_vec,PETSC_NULL_INTEGER,max,ierr);CHKERRQ(ierr)
      call VecMin(global_vec,PETSC_NULL_INTEGER,min,ierr);CHKERRQ(ierr)
      call VecSum(global_vec,mean,ierr);CHKERRQ(ierr)
      mean = mean / geomechanics_realization%geomech_patch%geomech_grid%nmax_node
      
      ! list of natural ids
      if (associated(geomechanics_regression%natural_vertex_ids)) then
        call VecScatterBegin(geomechanics_regression%scatter_natural_vertex_id_gtos, &
                             global_vec, &
                             geomechanics_regression%natural_vertex_id_vec, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        call VecScatterEnd(geomechanics_regression%scatter_natural_vertex_id_gtos, &
                           global_vec, &
                           geomechanics_regression%natural_vertex_id_vec, &
                           INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
      endif
      if (geomechanics_regression%num_vertices_per_process > 0) then
        ! vertices per process
        call VecScatterBegin(geomechanics_regression%scatter_vertices_per_process_gtos, &
                             global_vec, &
                             geomechanics_regression%vertices_per_process_vec, &
                             INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr)
        call VecScatterEnd(geomechanics_regression%scatter_vertices_per_process_gtos, &
                           global_vec, &
                           geomechanics_regression%vertices_per_process_vec, &
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
        
        ! natural vertex ids
        if (associated(geomechanics_regression%natural_vertex_ids)) then
          if (size(geomechanics_regression%natural_vertex_ids) > 0) then
            call VecGetArrayF90(geomechanics_regression%natural_vertex_id_vec,vec_ptr, &
                                ierr);CHKERRQ(ierr)
            if (cur_variable%iformat == 0) then
              do i = 1, size(geomechanics_regression%natural_vertex_ids)
                write(OUTPUT_UNIT,100) &
                  geomechanics_regression%natural_vertex_ids(i),vec_ptr(i)
              enddo
            else
              do i = 1, size(geomechanics_regression%natural_vertex_ids)
                write(OUTPUT_UNIT,101) &
                  geomechanics_regression%natural_vertex_ids(i),int(vec_ptr(i))
              enddo
            endif
            call VecRestoreArrayF90(geomechanics_regression%natural_vertex_id_vec,vec_ptr, &
                                    ierr);CHKERRQ(ierr)
          endif
        endif
        
        ! vertex ids per process
        if (geomechanics_regression%num_vertices_per_process > 0) then
          call VecGetArrayF90(geomechanics_regression%vertices_per_process_vec,vec_ptr, &
                              ierr);CHKERRQ(ierr)
          if (cur_variable%iformat == 0) then
            do i = 1, geomechanics_regression%num_vertices_per_process*option%mycommsize
              write(OUTPUT_UNIT,100) &
                geomechanics_regression%vertices_per_process_natural_ids(i),vec_ptr(i)
            enddo
          else
            do i = 1, geomechanics_regression%num_vertices_per_process*option%mycommsize
              write(OUTPUT_UNIT,101) &
                geomechanics_regression%vertices_per_process_natural_ids(i),int(vec_ptr(i))
            enddo
          endif
          call VecRestoreArrayF90(geomechanics_regression%vertices_per_process_vec,vec_ptr, &
                                  ierr);CHKERRQ(ierr)
        endif
      endif
    endif

    cur_variable => cur_variable%next
  enddo
  
  ! timestep, newton iteration, solver iteration output
  if (associated(geomechanics_timestepper)) then
    call VecNorm(geomechanics_realization%geomech_field%disp_xx, &
                 NORM_2,x_norm,ierr);CHKERRQ(ierr)
    call VecNorm(geomechanics_realization%geomech_field%disp_r, &
                 NORM_2,r_norm,ierr);CHKERRQ(ierr)
    if (option%myrank == option%io_rank) then
      write(OUTPUT_UNIT,'(''-- SOLUTION: Geomechanics --'')')
      write(OUTPUT_UNIT,'(''   Time (seconds): '',es21.13)') &
        geomechanics_timestepper%cumulative_solver_time
      write(OUTPUT_UNIT,'(''   Time Steps: '',i12)') geomechanics_timestepper%steps
      write(OUTPUT_UNIT,'(''   Newton Iterations: '',i12)') &
        geomechanics_timestepper%cumulative_newton_iterations
      write(OUTPUT_UNIT,'(''   Solver Iterations: '',i12)') &
        geomechanics_timestepper%cumulative_linear_iterations
      write(OUTPUT_UNIT,'(''   Time Step Cuts: '',i12)') &
        geomechanics_timestepper%cumulative_time_step_cuts
      write(OUTPUT_UNIT,'(''   Solution 2-Norm: '',es21.13)') x_norm
      write(OUTPUT_UNIT,'(''   Residual 2-Norm: '',es21.13)') r_norm
    endif
  endif

  close(OUTPUT_UNIT)
  
end subroutine GeomechanicsRegressionOutput

! ************************************************************************** !

recursive subroutine GeomechanicsRegressionVariableDestroy( &
                                            geomechanics_regression_variable)
  ! 
  ! Destroys a geomechanics regression variable object
  ! 
  ! Author: Satish Karra
  ! Date: 06/22/2016
  ! 

  implicit none
  
  type(geomechanics_regression_variable_type), pointer :: &
                                              geomechanics_regression_variable
  
  if (.not.associated(geomechanics_regression_variable)) return
  
  call GeomechanicsRegressionVariableDestroy( &
                                        geomechanics_regression_variable%next)

  deallocate(geomechanics_regression_variable)
  nullify(geomechanics_regression_variable)
  
end subroutine GeomechanicsRegressionVariableDestroy

! ************************************************************************** !

subroutine GeomechanicsRegressionDestroy(geomechanics_regression)
  ! 
  ! Destroys a geomechanics regression object
  ! 
  ! Author: Satish Karra
  ! Date: 06/22/2016
  ! 

  use Utility_module
  
  implicit none
  
  type(geomechanics_regression_type), pointer :: geomechanics_regression
  
  PetscErrorCode :: ierr
  
  if (.not.associated(geomechanics_regression)) return
  
  call GeomechanicsRegressionVariableDestroy( &
                                        geomechanics_regression%variable_list)
  call DeallocateArray(geomechanics_regression%natural_vertex_ids)
  geomechanics_regression%num_vertices_per_process = 0
  call DeallocateArray( &
                geomechanics_regression%vertices_per_process_natural_ids)
  if (geomechanics_regression%natural_vertex_id_vec /= 0) then
    call VecDestroy(geomechanics_regression%natural_vertex_id_vec, &
                    ierr);CHKERRQ(ierr)
  endif
  if (geomechanics_regression%vertices_per_process_vec /= 0) then
    call VecDestroy(geomechanics_regression%vertices_per_process_vec, &
                    ierr);CHKERRQ(ierr)
  endif
  if (geomechanics_regression%scatter_natural_vertex_id_gtos /= 0) then
    call VecScatterDestroy( &
                    geomechanics_regression%scatter_natural_vertex_id_gtos, &
                    ierr);CHKERRQ(ierr)
  endif
  if (geomechanics_regression%scatter_vertices_per_process_gtos /= 0) then
    call VecScatterDestroy( &
                  geomechanics_regression%scatter_vertices_per_process_gtos, &
                  ierr);CHKERRQ(ierr)
  endif

  deallocate(geomechanics_regression)
  nullify(geomechanics_regression)
  
end subroutine GeomechanicsRegressionDestroy

end module Geomechanics_Regression_module
