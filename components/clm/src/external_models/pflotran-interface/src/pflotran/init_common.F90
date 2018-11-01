module Init_Common_module
#include "petsc/finclude/petscts.h"
  use petscts
  use PFLOTRAN_Constants_module

  implicit none

  private

  public :: &
!            Init, &
            InitCommonReadRegionFiles, &
#if defined(PETSC_HAVE_HDF5)
            InitCommonReadVelocityField, &
#endif
            InitCommonVerifyAllCouplers, &
            setSurfaceFlowMode, &
            InitCommonAddOutputWaypoints

#if defined(SCORPIO)
  public :: InitCommonCreateIOGroups
#endif  
  
contains

! ************************************************************************** !

subroutine InitReadInputFilenames(option,filenames)
  ! 
  ! Reads filenames for multi-simulation runs
  ! 
  ! Author: Glenn Hammond
  ! Date: 08/11/09
  ! 

  use Option_module
  use Input_Aux_module

  type(option_type) :: option
  character(len=MAXSTRINGLENGTH), pointer :: filenames(:)

  character(len=MAXSTRINGLENGTH) :: string
  character(len=MAXSTRINGLENGTH) :: filename
  PetscInt :: filename_count
  type(input_type), pointer :: input
  PetscBool :: card_found

  input => InputCreate(IN_UNIT,option%input_filename,option)

  string = "FILENAMES"
  call InputFindStringInFile(input,option,string) 

  card_found = PETSC_FALSE
  if (InputError(input)) then
    ! if the FILENAMES card is not included, we will assume that only
    ! filenames exist in the file.
    call InputRewind(input)
  else
    card_found = PETSC_TRUE
  endif
    
  filename_count = 0     
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
  enddo
  
  allocate(filenames(filename_count))
  filenames = ''
  call InputRewind(input)

  if (card_found) then
    string = "FILENAMES"
    call InputFindStringInFile(input,option,string) 
  endif
  
  filename_count = 0     
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit  
    call InputReadNChars(input,option,filename,MAXSTRINGLENGTH,PETSC_FALSE)
    filename_count = filename_count + 1
    filenames(filename_count) = filename
  enddo

  call InputDestroy(input)

end subroutine InitReadInputFilenames

! ************************************************************************** !

subroutine setSurfaceFlowMode(option)
  ! 
  ! Sets the flow mode for surface (richards, th, etc.)
  ! 
  ! Author: Gautam Bisht
  ! Date: 07/30/14
  ! 

  use Option_module
  use String_module

  implicit none 

  type(option_type) :: option
  
  select case(option%iflowmode)
    case(TH_MODE)
      option%nsurfflowdof = TWO_INTEGER
    case default
      write(option%io_buffer,*) option%iflowmode
      option%io_buffer = 'Flow Mode ' // &
        trim(option%io_buffer) // ' not recognized in setSurfaceFlowMode().'
      call printErrMsg(option)
  end select
  
end subroutine setSurfaceFlowMode

! ************************************************************************** !

subroutine InitCommonVerifyAllCouplers(realization)
  ! 
  ! Verifies the connectivity of a coupler
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/8/08
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Coupler_module

  implicit none

  class(realization_subsurface_type) :: realization
  
  type(patch_type), pointer :: cur_patch

  cur_patch => realization%patch_list%first
  do
    if (.not.associated(cur_patch)) exit

      call InitCommonVerifyCoupler(realization,cur_patch, &
                                   cur_patch%initial_condition_list)
      call InitCommonVerifyCoupler(realization,cur_patch, &
                                   cur_patch%boundary_condition_list)
      call InitCommonVerifyCoupler(realization,cur_patch, &
                                   cur_patch%source_sink_list)

    cur_patch => cur_patch%next
  enddo
  
end subroutine InitCommonVerifyAllCouplers

! ************************************************************************** !

subroutine InitCommonVerifyCoupler(realization,patch,coupler_list)
  ! 
  ! Verifies the connectivity of a coupler
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/8/08
  ! 

  use Realization_Subsurface_class
  use Discretization_module
  use Option_module 
  use Coupler_module
  use Condition_module
  use Grid_module
  use Output_module
  use Patch_module

  implicit none

  class(realization_subsurface_type) :: realization
  type(coupler_list_type), pointer :: coupler_list

  type(option_type), pointer :: option
  type(grid_type), pointer :: grid
  type(patch_type), pointer :: patch  
  type(coupler_type), pointer :: coupler
  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: filename
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: iconn, icell, local_id
  Vec :: global_vec
  PetscReal, pointer :: vec_ptr(:)
  PetscErrorCode :: ierr

  patch => realization%patch
  grid => patch%grid
  option => realization%option

  if (.not.associated(coupler_list)) return

  call DiscretizationCreateVector(realization%discretization,ONEDOF, &
                                  global_vec,GLOBAL,option)

  coupler => coupler_list%first

  do
    if (.not.associated(coupler)) exit

    call VecZeroEntries(global_vec,ierr);CHKERRQ(ierr)
    call VecGetArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
    if (associated(coupler%connection_set)) then
      do iconn = 1, coupler%connection_set%num_connections
        local_id = coupler%connection_set%id_dn(iconn)
!        vec_ptr(local_id) = coupler%id
!geh: let's sum the # of connections
         vec_ptr(local_id) = vec_ptr(local_id) + 1
      enddo
    else
      if (associated(coupler%region)) then
        do icell = 1, coupler%region%num_cells
          local_id = coupler%region%cell_ids(icell)
!          vec_ptr(local_id) = coupler%id
         vec_ptr(local_id) = vec_ptr(local_id) + 1
        enddo
      endif
    endif
    call VecRestoreArrayF90(global_vec,vec_ptr,ierr);CHKERRQ(ierr)
    if (len_trim(coupler%flow_condition_name) > 0) then
      dataset_name = coupler%flow_condition_name
    else if (len_trim(coupler%tran_condition_name) > 0) then
      dataset_name = coupler%tran_condition_name
    endif
    !write(word,*) patch%id
    !dataset_name = trim(dataset_name) // '_' // &
    !               trim(coupler%region%name) // '_' // &
    !               trim(adjustl(word))
    !dataset_name = dataset_name(1:28)
    !filename = trim(dataset_name) // '.tec'
    !call OutputVectorTecplot(filename,dataset_name,realization,global_vec)

    coupler => coupler%next
  enddo

  call VecDestroy(global_vec,ierr);CHKERRQ(ierr)

end subroutine InitCommonVerifyCoupler

! ************************************************************************** !

subroutine InitCommonReadRegionFiles(realization)
  ! 
  ! Reads in grid cell ids stored in files
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/03/08
  ! 

  use Realization_Subsurface_class
  use Region_module
  use HDF5_module
  use Option_module

  implicit none

  class(realization_subsurface_type) :: realization

  type(option_type), pointer :: option
  type(region_type), pointer :: region
  PetscBool :: cell_ids_exists
  PetscBool :: face_ids_exists
  PetscBool :: vert_ids_exists
 
  option => realization%option
  region => realization%region_list%first
  do 
    if (.not.associated(region)) exit
    if (len_trim(region%filename) > 1) then
      if (index(region%filename,'.h5') > 0) then
        call HDF5QueryRegionDefinition(region, region%filename, &
                                       realization%option, cell_ids_exists, &
                                       face_ids_exists, vert_ids_exists)
        if ((.not. cell_ids_exists) .and. &
            (.not. face_ids_exists) .and. &
            (.not. vert_ids_exists)) then
          option%io_buffer = '"Regions/' // trim(region%name) // &
                ' is not defined by "Cell Ids" or "Face Ids" or "Vertex Ids".'
          call printErrMsg(option)
        end if
        if (cell_ids_exists .or. face_ids_exists) then
          call HDF5ReadRegionFromFile(realization%patch%grid,region, &
                                      region%filename,option)
        else
          call HDF5ReadRegionDefinedByVertex(realization%option, &
                                             region, region%filename)
        end if
      else if (index(region%filename,'.ss') > 0) then
        region%def_type = DEFINED_BY_SIDESET_UGRID
        region%sideset => RegionCreateSideset()
        call RegionReadFromFile(region%sideset,region%filename, &
                                realization%option)
      else if (index(region%filename,'.ex') > 0) then
        region%def_type = DEFINED_BY_FACE_UGRID_EXP
        call RegionReadFromFile(region%explicit_faceset,region%cell_ids, &
                                region%filename,realization%option)
        region%num_cells = size(region%cell_ids)
      else
        call RegionReadFromFile(region,realization%option, &
                                region%filename)
      endif
    endif
    region => region%next
  enddo

end subroutine InitCommonReadRegionFiles

! ************************************************************************** !

subroutine readVectorFromFile(realization,vector,filename,vector_type)
  ! 
  ! Reads data from a file into an associated vector
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/18/08
  ! 

  use Realization_Subsurface_class
  use Discretization_module
  use Field_module
  use Grid_module
  use Option_module
  use Patch_module
  use Logging_module

  use HDF5_module
  
  implicit none
  
  class(realization_subsurface_type) :: realization
  Vec :: vector
  character(len=MAXWORDLENGTH) :: filename
  PetscInt :: vector_type
  
  type(discretization_type), pointer :: discretization
  type(field_type), pointer :: field
  type(grid_type), pointer :: grid
  type(option_type), pointer :: option
  type(patch_type), pointer :: patch   
  PetscInt :: ghosted_id, natural_id, material_id
  PetscInt :: fid = 86
  PetscInt :: status
  PetscErrorCode :: ierr
  PetscInt :: count, read_count, i
  PetscInt :: flag
  PetscInt, pointer :: indices(:)
  PetscReal, pointer :: values(:)
  PetscInt, parameter :: block_size = 10000
  Vec :: natural_vec, global_vec

  discretization => realization%discretization
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option

  if (index(filename,'.h5') > 0) then
    ! to be taken care of later in readPermeabilitiesFromFile()
  else
    open(unit=fid,file=filename,status="old",iostat=status)
    if (status /= 0) then
      option%io_buffer = 'File: ' // trim(filename) // ' not found.'
      call printErrMsg(option)
    endif
    allocate(values(block_size))
    allocate(indices(block_size))
    call DiscretizationCreateVector(discretization,ONEDOF,natural_vec, &
                                    NATURAL,option)
    count = 0
    do
      if (count >= grid%nmax) exit
      read_count = min(block_size,grid%nmax-count)
      do i=1,read_count
        indices(i) = count+i-1 ! zero-based indexing
      enddo
      ierr = 0
      if (option%myrank == option%io_rank) &
        read(fid,*,iostat=ierr) values(1:read_count)
      flag = ierr
      call MPI_Bcast(flag,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                     option%mycomm,ierr)      
      if (flag /= 0) then
        option%io_buffer = 'Insufficent data in file: ' // filename
        call printErrMsg(option)
      endif
      if (option%myrank == option%io_rank) then
        call VecSetValues(natural_vec,read_count,indices,values,INSERT_VALUES, &
                          ierr);CHKERRQ(ierr)
      endif
      count = count + read_count
    enddo
    call MPI_Bcast(count,ONE_INTEGER_MPI,MPIU_INTEGER,option%io_rank, &
                   option%mycomm,ierr)      
    if (count /= grid%nmax) then
      write(option%io_buffer,'("Number of data in file (",i8, &
      & ") does not match size of vector (",i8,")")') count, grid%nlmax
      call printErrMsg(option)
    endif
    close(fid)
    deallocate(values)
    nullify(values)
    deallocate(indices)
    nullify(indices)
    call VecAssemblyBegin(natural_vec,ierr);CHKERRQ(ierr)
    call VecAssemblyEnd(natural_vec,ierr);CHKERRQ(ierr)
    select case(vector_type)
      case(LOCAL)
        call DiscretizationCreateVector(discretization,ONEDOF,global_vec, &
                                        GLOBAL,option)        
        call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                           global_vec,ONEDOF)  
        call DiscretizationGlobalToLocal(discretization,global_vec, &
                                         vector,ONEDOF)
        call VecDestroy(global_vec,ierr);CHKERRQ(ierr)
      case(GLOBAL)
        call DiscretizationNaturalToGlobal(discretization,natural_vec, &
                                           vector,ONEDOF) 
    end select 
    call VecDestroy(natural_vec,ierr);CHKERRQ(ierr)
  endif
  
end subroutine readVectorFromFile

! ************************************************************************** !

subroutine InitCommonCreateIOGroups(option)
  ! 
  ! Create sub-communicators that are used in initialization
  ! and output HDF5 routines.
  ! 
  ! Author: Vamsi Sripathi
  ! Date: 07/14/09
  ! 

  use Option_module
  use Logging_module

#if defined(SCORPIO)
  use hdf5
#endif

  implicit none

  type(option_type) :: option
  PetscErrorCode :: ierr

#if defined(SCORPIO)

  PetscMPIInt :: numiogroups

  call PetscLogEventBegin(logging%event_create_iogroups,ierr);CHKERRQ(ierr)

  ! Initialize HDF interface to define global constants  
  call h5open_f(ierr)

  if (option%hdf5_read_group_size <= 0) then
    write(option%io_buffer,& 
          '("The keyword HDF5_READ_GROUP_SIZE & 
            & in the input file (pflotran.in) is either not set or &
            & its value is less than or equal to ZERO. &
            & HDF5_READ_GROUP_SIZE =  ",i6)') &
             option%hdf5_read_group_size
    !call printErrMsg(option)
    call printMsg(option)
    ! default is to let one process read and broadcast to everyone
    option%hdf5_read_group_size = option%mycommsize
  endif         
 
  if (option%hdf5_write_group_size <= 0) then
    write(option%io_buffer,& 
          '("The keyword HDF5_WRITE_GROUP_SIZE & 
            &in the input file (pflotran.in) is either not set or &
            &its value is less than or equal to ZERO. &
            &HDF5_WRITE_GROUP_SIZE =  ",i6)') &
             option%hdf5_write_group_size
    !call printErrMsg(option)
    call printMsg(option)
    ! default is to let everyone write separately 
    option%hdf5_write_group_size = 1
  endif                    

  ! create read IO groups
  numiogroups = option%mycommsize/option%hdf5_read_group_size
  call fscorpio_iogroup_init(numiogroups, option%mycomm, &
                             option%ioread_group_id, ierr)

  if ( option%hdf5_read_group_size == option%hdf5_write_group_size ) then
    ! reuse read_group to use for writing too as both groups are same size
    option%iowrite_group_id = option%ioread_group_id
  else   
      ! create write IO groups
      numiogroups = option%mycommsize/option%hdf5_write_group_size
      call fscorpio_iogroup_init(numiogroups, option%mycomm, &
                                 option%iowrite_group_id, ierr)
  end if

    write(option%io_buffer, '(" Read group id :  ", i6)') option%ioread_group_id
    call printMsg(option)      
    write(option%io_buffer, '(" Write group id :  ", i6)') &
      option%iowrite_group_id
    call printMsg(option)      
  call PetscLogEventEnd(logging%event_create_iogroups,ierr);CHKERRQ(ierr)
#endif
! SCORPIO
 
end subroutine InitCommonCreateIOGroups

! ************************************************************************** !

subroutine InitCommonPrintPFLOTRANHeader(option,fid)
  ! 
  ! Initializes pflotran
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/23/07
  ! 

  use Option_module
  
  implicit none
  
  PetscInt :: fid
  
  type(option_type) :: option
  
  write(fid,'(" PFLOTRAN Header")') 
  
end subroutine InitCommonPrintPFLOTRANHeader

! ************************************************************************** !

#if defined(PETSC_HAVE_HDF5)

subroutine InitCommonReadVelocityField(realization)
  ! 
  ! Reads fluxes in for transport with no flow.
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/05/13
  ! 

  use Realization_Subsurface_class
  use Patch_module
  use Field_module
  use Grid_module
  use Option_module
  use Coupler_module
  use Connection_module
  use Discretization_module
  use HDF5_module
  use HDF5_Aux_module

  implicit none
  
  class(realization_subsurface_type) :: realization
  character(len=MAXSTRINGLENGTH) :: filename
  
  type(field_type), pointer :: field
  type(patch_type), pointer :: patch
  type(grid_type), pointer :: grid
  type(discretization_type), pointer :: discretization
  type(option_type), pointer :: option
  character(len=MAXSTRINGLENGTH) :: group_name
  character(len=MAXSTRINGLENGTH) :: dataset_name
  PetscInt :: idir, iconn, sum_connection
  PetscInt :: ghosted_id_up, local_id
  PetscErrorCode :: ierr
  
  PetscReal, pointer :: vec_loc_p(:)
  PetscReal, pointer :: vec_p(:)
  type(coupler_type), pointer :: boundary_condition  
  type(connection_set_list_type), pointer :: connection_set_list
  type(connection_set_type), pointer :: cur_connection_set
  
  field => realization%field
  patch => realization%patch
  grid => patch%grid
  option => realization%option
  discretization => realization%discretization
  
  filename = realization%nonuniform_velocity_filename

  group_name = ''
  do idir = 1, 3
    select case(idir)
      case(1)
        dataset_name = 'Internal Velocity X'
      case(2)
        dataset_name = 'Internal Velocity Y'
      case(3)
        dataset_name = 'Internal Velocity Z'
    end select
    if (.not.HDF5DatasetExists(filename,group_name,dataset_name,option)) then
      option%io_buffer = 'Dataset "' // trim(group_name) // '/' // &
        trim(dataset_name) // &
        '" not found in HDF5 file "' // trim(filename) // '".'
      call printErrMsg(option)
    endif
    call HDF5ReadCellIndexedRealArray(realization,field%work,filename, &
                                      group_name,dataset_name,PETSC_FALSE)
    call DiscretizationGlobalToLocal(discretization,field%work,field%work_loc, &
                                     ONEDOF)
    call VecGetArrayF90(field%work_loc,vec_loc_p,ierr);CHKERRQ(ierr)
    connection_set_list => grid%internal_connection_set_list
    cur_connection_set => connection_set_list%first
    sum_connection = 0  
    do 
      if (.not.associated(cur_connection_set)) exit
      do iconn = 1, cur_connection_set%num_connections
        sum_connection = sum_connection + 1
        ghosted_id_up = cur_connection_set%id_up(iconn)
        if (cur_connection_set%dist(idir,iconn) > 0.9d0) then
          patch%internal_velocities(1,sum_connection) = vec_loc_p(ghosted_id_up)
        endif
      enddo
      cur_connection_set => cur_connection_set%next
    enddo
    call VecRestoreArrayF90(field%work_loc,vec_loc_p,ierr);CHKERRQ(ierr)
  enddo
  
  boundary_condition => patch%boundary_condition_list%first
  sum_connection = 0    
  do 
    if (.not.associated(boundary_condition)) exit
    dataset_name = boundary_condition%name
    if (.not.HDF5DatasetExists(filename,group_name,dataset_name,option)) then
      option%io_buffer = 'Dataset "' // trim(group_name) // '/' // &
        trim(dataset_name) // &
        '" not found in HDF5 file "' // trim(filename) // '".'
      call printErrMsg(option)
    endif
    call HDF5ReadCellIndexedRealArray(realization,field%work,filename, &
                                      group_name,dataset_name,PETSC_FALSE)
    call VecGetArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    cur_connection_set => boundary_condition%connection_set
    do iconn = 1, cur_connection_set%num_connections
      sum_connection = sum_connection + 1
      local_id = cur_connection_set%id_dn(iconn)
      patch%boundary_velocities(1,sum_connection) = vec_p(local_id)
    enddo
    call VecRestoreArrayF90(field%work,vec_p,ierr);CHKERRQ(ierr)
    boundary_condition => boundary_condition%next
  enddo
  
end subroutine InitCommonReadVelocityField

#endif

! ************************************************************************** !

subroutine InitCommonAddOutputWaypoints(option,output_option,waypoint_list)
  ! 
  ! Adds waypoints associated with output options to waypoint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/04/16
  ! 
  use Output_Aux_module
  use Waypoint_module
  use Option_module
  use Utility_module
  
  implicit none
  
  type(option_type) :: option
  type(output_option_type) :: output_option
  type(waypoint_list_type) :: waypoint_list
  
  type(waypoint_type), pointer :: waypoint
  character(len=MAXWORDLENGTH) :: word
  PetscReal :: temp_real
  PetscReal :: final_time
  PetscReal :: num_waypoints, warning_num_waypoints
  PetscInt :: k

  !geh: The repetitive summation of a time increment can result in slight 
  !     error.   The perturbation is designed to allow for a slight shift 
  !     beyond the final time.
  final_time = WaypointListGetFinalTime(waypoint_list)
  temp_real = final_time * 1.d-10
  final_time = final_time + temp_real
  warning_num_waypoints = 15000.0
  
  ! Add waypoints for periodic snapshot output
  if (output_option%periodic_snap_output_time_incr > 0.d0) then
    temp_real = 0.d0
    num_waypoints = final_time / output_option%periodic_snap_output_time_incr
    if ((num_waypoints > warning_num_waypoints) .and. &
        OptionPrintToScreen(option)) then
       write(word,*) floor(num_waypoints)
       write(*,*) 'WARNING: Large number (' // trim(adjustl(word)) // &
                  ') of periodic snapshot output requested.'
      write(*,'(a64)',advance='no') '         Creating periodic output &
                                    &waypoints . . . Progress: 0%-'
    endif
    k = 0
    do
      k = k + 1
      temp_real = temp_real + output_option%periodic_snap_output_time_incr
      if (temp_real > final_time) exit
      waypoint => WaypointCreate()
      waypoint%time = temp_real
      waypoint%print_snap_output = PETSC_TRUE
      call WaypointInsertInList(waypoint,waypoint_list)
      if ((num_waypoints > warning_num_waypoints) .and. &
          OptionPrintToScreen(option)) then
        call PrintProgressBarInt(num_waypoints,TEN_INTEGER,k)
      endif
    enddo
  endif

  ! Add waypoints for periodic observation output
  if (output_option%periodic_obs_output_time_incr > 0.d0) then
    temp_real = 0.d0
    num_waypoints = final_time / output_option%periodic_obs_output_time_incr
    if ((num_waypoints > warning_num_waypoints) .and. &
        OptionPrintToScreen(option)) then
       write(word,*) floor(num_waypoints)
       write(*,*) 'WARNING: Large number (' // trim(adjustl(word)) // &
                  ') of periodic observation output requested.'
      write(*,'(a64)',advance='no') '         Creating periodic output &
                                    &waypoints . . . Progress: 0%-'
    endif
    k = 0
    do
      k = k + 1
      temp_real = temp_real + output_option%periodic_obs_output_time_incr
      if (temp_real > final_time) exit
      waypoint => WaypointCreate()
      waypoint%time = temp_real
      waypoint%print_obs_output = PETSC_TRUE
      call WaypointInsertInList(waypoint,waypoint_list)
      if ((num_waypoints > warning_num_waypoints) .and. &
          OptionPrintToScreen(option)) then
        call PrintProgressBarInt(num_waypoints,TEN_INTEGER,k)
      endif
    enddo
  endif

  ! Add waypoints for periodic mass balance output
  if (output_option%periodic_msbl_output_time_incr > 0.d0) then
    temp_real = 0.d0
    num_waypoints = final_time / output_option%periodic_msbl_output_time_incr
    if ((num_waypoints > warning_num_waypoints) .and. &
        OptionPrintToScreen(option)) then
       write(word,*) floor(num_waypoints)
       write(*,*) 'WARNING: Large number (' // trim(adjustl(word)) // &
                  ') of periodic mass balance output requested.'
      write(*,'(a64)',advance='no') '         Creating periodic output &
                                    &waypoints . . . Progress: 0%-'
    endif
    k = 0
    do
      k = k + 1
      temp_real = temp_real + output_option%periodic_msbl_output_time_incr
      if (temp_real > final_time) exit
      waypoint => WaypointCreate()
      waypoint%time = temp_real
      waypoint%print_msbl_output = PETSC_TRUE
      call WaypointInsertInList(waypoint,waypoint_list)
      if ((num_waypoints > warning_num_waypoints) .and. &
          OptionPrintToScreen(option)) then
        call PrintProgressBarInt(num_waypoints,TEN_INTEGER,k)
      endif
    enddo
  endif 
  
end subroutine InitCommonAddOutputWaypoints

end module Init_Common_module
