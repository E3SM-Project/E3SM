module Waypoint_module
 
  use Option_module
  use PFLOTRAN_Constants_module
  use Utility_module, only : Equal
  
  implicit none
  
  private

#include "petsc/finclude/petscsys.h"

  ! linked-list for waypoints in the simulation
  type, public :: waypoint_type
    PetscReal :: time
    PetscBool :: sync
    PetscBool :: print_snap_output
    PetscBool :: print_obs_output
    PetscBool :: print_msbl_output
    PetscBool :: print_checkpoint
!    type(output_option_type), pointer :: output_option
    PetscBool :: update_conditions
    PetscReal :: dt_max
    PetscBool :: final  ! any waypoint after this will be deleted
    type(waypoint_type), pointer :: prev
    type(waypoint_type), pointer :: next
  end type waypoint_type
  
  type, public :: waypoint_list_type
    PetscInt :: num_waypoints
    type(waypoint_type), pointer :: first
    type(waypoint_type), pointer :: last
    type(waypoint_type), pointer :: array(:)    
  end type waypoint_list_type
  
  interface WaypointCreate
    module procedure WaypointCreate1
    module procedure WaypointCreate2
  end interface  
  
  public :: WaypointCreate, &
            WaypointListCreate, &
            WaypointListDestroy, &
            WaypointInsertInList, &
            WaypointDeleteFromList, &
            WaypointListFillIn, &
            WaypointListCopy, &
            WaypointListMerge, &
            WaypointListCopyAndMerge, &
            WaypointListRemoveExtraWaypnts, &
            WaypointConvertTimes, &
            WaypointReturnAtTime, &
            WaypointSkipToTime, &
            WaypointForceMatchToTime, &
            WaypointListPrint, &
            WaypointListGetFinalTime, &
            WaypointCreateSyncWaypointList, &
            WaypointInputRecord

contains

! ************************************************************************** !

function WaypointCreate1()
  ! 
  ! Creates a simulation waypoint
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none
  
  type(waypoint_type), pointer :: WaypointCreate1
  
  type(waypoint_type), pointer :: waypoint
  
  allocate(waypoint)
  waypoint%time = 0.d0
  waypoint%sync = PETSC_FALSE
  waypoint%print_snap_output = PETSC_FALSE
  waypoint%print_obs_output = PETSC_FALSE
  waypoint%print_msbl_output = PETSC_FALSE
  waypoint%print_checkpoint = PETSC_FALSE
  waypoint%final = PETSC_FALSE
  waypoint%update_conditions = PETSC_FALSE
  waypoint%dt_max = 0.d0
  nullify(waypoint%next)
  nullify(waypoint%prev)
    
  WaypointCreate1 => waypoint
  
end function WaypointCreate1

! ************************************************************************** !

function WaypointCreate2(original_waypoint)
  ! 
  ! Creates a simulation waypoint
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 

  implicit none
  
  type(waypoint_type), pointer :: original_waypoint
  
  type(waypoint_type), pointer :: WaypointCreate2
  
  type(waypoint_type), pointer :: waypoint
  
  waypoint => WaypointCreate()
  waypoint%time = original_waypoint%time
  waypoint%sync = original_waypoint%sync
  waypoint%print_snap_output = original_waypoint%print_snap_output
  waypoint%print_obs_output = original_waypoint%print_obs_output
  waypoint%print_msbl_output = original_waypoint%print_msbl_output
  waypoint%print_checkpoint = original_waypoint%print_checkpoint
  waypoint%final = original_waypoint%final
  waypoint%update_conditions = original_waypoint%update_conditions
  waypoint%dt_max = original_waypoint%dt_max
    
  WaypointCreate2 => waypoint
  
end function WaypointCreate2

! ************************************************************************** !

function WaypointListCreate()
  ! 
  ! Creates a simulation waypoint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 

  implicit none
  
  type(waypoint_list_type), pointer :: WaypointListCreate
  
  type(waypoint_list_type), pointer :: waypoint_list
  
  allocate(waypoint_list)
  nullify(waypoint_list%first)
  nullify(waypoint_list%last)
  nullify(waypoint_list%array)
  waypoint_list%num_waypoints = 0

  WaypointListCreate => waypoint_list
  
end function WaypointListCreate 


! ************************************************************************** !

subroutine WaypointListMerge(waypoint_list1,waypoint_list2,option)
  ! 
  ! Creates a simulation waypoint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/16
  ! 
  use Option_module
  
  implicit none
  
  type(waypoint_list_type), pointer :: waypoint_list1
  type(waypoint_list_type), pointer :: waypoint_list2
  
  type(option_type) :: option
  type(waypoint_type), pointer :: cur_waypoint, next_waypoint
  
  if (.not.associated(waypoint_list1) .and. &
      .not.associated(waypoint_list2)) then
    option%io_buffer = 'Two null waypoints lists.  Send input deck to &
      &pflotran-dev.'
    call printErrMsg(option)
  else if (.not.associated(waypoint_list1)) then
    waypoint_list1 => waypoint_list2
    return
  else if (.not.associated(waypoint_list2)) then
    waypoint_list2 => waypoint_list1
    return
  endif
  
  cur_waypoint => waypoint_list2%first
  do
    if (.not.associated(cur_waypoint)) exit
    next_waypoint => cur_waypoint%next
    nullify(cur_waypoint%next)
    call WaypointInsertInList(cur_waypoint,waypoint_list1)
    cur_waypoint => next_waypoint
    nullify(next_waypoint)
  enddo
  ! must nullify the first waypoint in waypoint_list2 to avoid deleting 
  ! first waypoint which will subsequently delete all waypoints after it 
  ! in waypoint_list1
  nullify(waypoint_list2%first)
  call WaypointListDestroy(waypoint_list2)
  waypoint_list2 => waypoint_list1
  
end subroutine WaypointListMerge 

! ************************************************************************** !

subroutine WaypointListCopyAndMerge(waypoint_list1,waypoint_list2,option)
  ! 
  ! Creates a simulation waypoint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 02/03/16
  ! 
  use Option_module
  
  implicit none
  
  type(waypoint_list_type), pointer :: waypoint_list1
  type(waypoint_list_type), pointer :: waypoint_list2
  
  type(option_type) :: option
 
  type(waypoint_list_type), pointer :: new_waypoint_list

  new_waypoint_list => WaypointListCopy(waypoint_list2)
  call WaypointListMerge(waypoint_list1,new_waypoint_list,option)
  nullify(new_waypoint_list)
  
end subroutine WaypointListCopyAndMerge 

! ************************************************************************** !

subroutine WaypointInsertInList(new_waypoint,waypoint_list)
  ! 
  ! Correctly inserts a waypoing in a list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/09/07
  ! 

  use Utility_module

  type(waypoint_type), pointer :: new_waypoint
  type(waypoint_list_type) :: waypoint_list

  type(waypoint_type), pointer :: waypoint
    
  ! place new waypoint in proper location within list
  waypoint => waypoint_list%first
  if (associated(waypoint)) then ! list exists
    if (new_waypoint%time < waypoint%time) then 
      ! insert at beginning of list
      waypoint_list%first => new_waypoint
      new_waypoint%next => waypoint
      new_waypoint%next%prev => new_waypoint
    else
      ! find its location in the list
      do
        if (Equal(new_waypoint%time,waypoint%time)) then
          call WaypointMerge(waypoint,new_waypoint)
          return ! do not increment num_waypoints at bottom
        else if (associated(waypoint%next)) then 
          if (new_waypoint%time-waypoint%time > 1.d-10 .and. & ! within list
              new_waypoint%time-waypoint%next%time < -1.d-10) then 
            new_waypoint%next => waypoint%next
            new_waypoint%next%prev => new_waypoint
            waypoint%next => new_waypoint
            new_waypoint%prev => waypoint
            exit
          else
            waypoint => waypoint%next
          endif
        else ! at end of list
          waypoint%next => new_waypoint
          new_waypoint%prev => waypoint
          waypoint_list%last => new_waypoint
          exit
        endif
      enddo
    endif
  else
    waypoint_list%first => new_waypoint
    waypoint_list%last => new_waypoint 
  endif
  waypoint_list%num_waypoints = waypoint_list%num_waypoints + 1

end subroutine WaypointInsertInList

! ************************************************************************** !

subroutine WaypointDeleteFromList(obsolete_waypoint,waypoint_list) ! 
  !
  ! Deletes a waypoint in a list
  ! 
  ! Author: Gautam Bisht
  ! Date: 01/20/11
  ! 
  use Utility_module

  implicit none

  type(waypoint_type), pointer :: obsolete_waypoint
  type(waypoint_type), pointer :: waypoint, prev_waypoint
  type(waypoint_list_type) :: waypoint_list

  waypoint => waypoint_list%first

  if (associated(waypoint)) then ! list exists

    ! Is the waypoint to be deleted is the first waypoint?
    if (Equal(waypoint%time,obsolete_waypoint%time)) then
      waypoint_list%first => waypoint%next
      call WaypointDestroy(waypoint)
      waypoint_list%num_waypoints = waypoint_list%num_waypoints - 1
      return
    else

      prev_waypoint => waypoint
      waypoint => waypoint%next
      do
        if (associated(waypoint)) then
          if (dabs(waypoint%time-obsolete_waypoint%time) < 1.d-10) then
            prev_waypoint%next => waypoint%next
            call WaypointDestroy(waypoint)
            waypoint_list%num_waypoints = waypoint_list%num_waypoints - 1
            return
          endif
          prev_waypoint => waypoint
          waypoint => waypoint%next
          cycle
        else
         ! at the end of the list, didn't find obsolete waypoint
          return
        endif
      enddo
    endif
  else
    ! list does not exists
    return
  endif
  
end subroutine WaypointDeleteFromList

! ************************************************************************** !

subroutine WaypointListFillIn(waypoint_list,option)
  ! 
  ! Fills in missing values (e.g. dt_max) in waypoint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/09/07
  ! 
  
  implicit none
  
  type(waypoint_list_type) :: waypoint_list
  type(option_type) :: option
  
  type(waypoint_type), pointer :: waypoint, prev_waypoint
  PetscReal :: dt_max = UNINITIALIZED_DOUBLE
  
  ! find first value of dt_max > 0.d0 in list
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint)) exit
    if (waypoint%dt_max > 1.d-40) then
      dt_max = waypoint%dt_max
      exit
    endif
    waypoint => waypoint%next
  enddo

  if (dt_max <= 1.d-40) then
    option%io_buffer = 'All values of dt_max in input file uninitialized'
    call printErrMsg(option)
  endif
  
  ! assign that value to the first waypoint, if waypoint%dt_max not already > 1.d-40
  waypoint => waypoint_list%first
  if (waypoint%dt_max < 1.d-40) waypoint%dt_max = dt_max
  
  ! fill in missing values
  do
    prev_waypoint => waypoint
    waypoint => waypoint%next
    if (.not.associated(waypoint)) exit 
    if (waypoint%dt_max < 1.d-40) then
      waypoint%dt_max = prev_waypoint%dt_max
    endif
  enddo
  
  ! IMPORTANT NOTE:  The dt_max must be assigned to the "next" waypoint.  The
  ! "current" waypoint in the stepper is always the next waypoint .  Therefore
  ! we must shift all the dt_max entries. 
  waypoint => waypoint_list%last
  ! work backwards
  do
    prev_waypoint => waypoint%prev
    if (.not.associated(prev_waypoint)) exit 
    waypoint%dt_max = prev_waypoint%dt_max
    waypoint => prev_waypoint
  enddo
  
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint)) exit 
    waypoint => waypoint%next
  enddo

end subroutine WaypointListFillIn 

! ************************************************************************** !

subroutine WaypointConvertTimes(waypoint_list,time_conversion)
  ! 
  ! Converts time units to seconds
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/09/07
  ! 

  implicit none
  
  type(waypoint_list_type) :: waypoint_list
  PetscReal :: time_conversion
  
  type(waypoint_type), pointer :: waypoint
  
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint)) exit
    waypoint%time = waypoint%time * time_conversion
    waypoint%dt_max = waypoint%dt_max * time_conversion
    waypoint => waypoint%next
  enddo
  
end subroutine WaypointConvertTimes 

! ************************************************************************** !

subroutine WaypointListRemoveExtraWaypnts(waypoint_list,option)
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/09/07
  ! 

  implicit none
  
  type(waypoint_list_type) :: waypoint_list
  type(option_type) :: option
  
  type(waypoint_type), pointer :: waypoint, prev_waypoint
  
  waypoint => waypoint_list%first
  do
    if (.not.associated(waypoint) .or. waypoint%final) exit
    waypoint => waypoint%next
  enddo
  
  if (associated(waypoint)) then
    prev_waypoint => waypoint
    waypoint => waypoint%next
    nullify(prev_waypoint%next)
  endif
  
  do
    if (.not.associated(waypoint)) exit
    prev_waypoint => waypoint
    waypoint => waypoint%next
    write(option%io_buffer,'("Waypoint at time:", 1pe12.4, &
  &       " is beyond the end of simulation")') &
          prev_waypoint%time
    call printWrnMsg(option)
    call WaypointDestroy(prev_waypoint)   
    waypoint_list%num_waypoints = waypoint_list%num_waypoints - 1
  enddo

end subroutine WaypointListRemoveExtraWaypnts 

! ************************************************************************** !

subroutine WaypointMerge(old_waypoint,new_waypoint)
  ! 
  ! Merges 2 waypoints performing an OR operation on logicals
  ! 
  ! Author: Glenn Hammond
  ! Date: 10/28/03
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none

  type(waypoint_type), pointer :: old_waypoint, new_waypoint

  new_waypoint%time = 0.d0

!    PetscReal :: time
!    PetscBool :: print_output
!    type(output_option_type), pointer :: output_option
!    PetscBool :: update_bcs
!    PetscBool :: update_srcs
!    PetscReal :: dt_max
!    PetscBool :: final  ! any waypoint after this will be deleted
    
  if (old_waypoint%sync .or. new_waypoint%sync) then
    old_waypoint%sync = PETSC_TRUE
  else
    old_waypoint%sync = PETSC_FALSE
  endif

  if (old_waypoint%print_snap_output .or. new_waypoint%print_snap_output) then
    old_waypoint%print_snap_output = PETSC_TRUE
  else
    old_waypoint%print_snap_output = PETSC_FALSE
  endif

  if (old_waypoint%print_obs_output .or. new_waypoint%print_obs_output) then
    old_waypoint%print_obs_output = PETSC_TRUE
  else
    old_waypoint%print_obs_output = PETSC_FALSE
  endif

  if (old_waypoint%print_msbl_output .or. new_waypoint%print_msbl_output) then
    old_waypoint%print_msbl_output = PETSC_TRUE
  else
    old_waypoint%print_msbl_output = PETSC_FALSE
  endif

  if (old_waypoint%update_conditions .or. new_waypoint%update_conditions) then
    old_waypoint%update_conditions = PETSC_TRUE
  else
    old_waypoint%update_conditions = PETSC_FALSE
  endif

  if (new_waypoint%dt_max > 0.d0) then
    old_waypoint%dt_max = new_waypoint%dt_max
  endif
  
  if (old_waypoint%final .or. new_waypoint%final) then
    old_waypoint%final = PETSC_TRUE
  else
    old_waypoint%final = PETSC_FALSE
  endif

  if (old_waypoint%print_checkpoint .or. new_waypoint%print_checkpoint) then
    old_waypoint%print_checkpoint = PETSC_TRUE
  else
    old_waypoint%print_checkpoint = PETSC_FALSE
  endif

  ! deallocate new waypoint
  deallocate(new_waypoint)
  ! point new_waypoint to old
  new_waypoint => old_waypoint

end subroutine WaypointMerge

! ************************************************************************** !

function WaypointReturnAtTime(list,time)
  ! 
  ! Returns a pointer to the first waypoint after time
  ! 
  ! Author: Glenn Hammond
  ! Date: 1/03/08
  ! 

  implicit none

  type(waypoint_list_type), pointer :: list
  PetscReal :: time

  type(waypoint_type), pointer :: WaypointReturnAtTime
  type(waypoint_type), pointer :: waypoint
  
  waypoint => list%first
  do 
    if (.not.associated(waypoint)) exit
    if (waypoint%time > time) exit
    waypoint => waypoint%next
  enddo

  if (associated(waypoint)) then
    WaypointReturnAtTime => waypoint
  else
    nullify(WaypointReturnAtTime)
  endif

end function WaypointReturnAtTime

! ************************************************************************** !

subroutine WaypointSkipToTime(cur_waypoint,time)
  ! 
  ! Skips the waypoint ahead to the correct time.
  ! 
  ! Author: Glenn Hammond
  ! Date: 07/31/13
  ! 

  implicit none

  PetscReal :: time
  type(waypoint_type), pointer :: cur_waypoint
  
  do 
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%time > time) exit
    cur_waypoint => cur_waypoint%next
  enddo

end subroutine WaypointSkipToTime

! ************************************************************************** !

subroutine WaypointListPrint(list,option,output_option)
  ! 
  ! Prints a waypoint
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/20/11
  ! 
  use Output_Aux_module
  use Option_module

  implicit none
  
  type(waypoint_list_type), pointer :: list
  type(option_type) :: option
  type(output_option_type) :: output_option

  type(waypoint_type), pointer :: cur_waypoint
  PetscInt :: icount

  100 format(/)
  110 format(a)
  20 format('  ',a20,':',10i6)

  if (OptionPrintToScreen(option)) then
    write(*,100)
    write(*,110) 'List of Waypoints:'
    write(*,100)
  endif

  if (OptionPrintToFile(option)) then
    write(option%fid_out,100)
    write(option%fid_out,110) 'List of Waypoints:'
    write(option%fid_out,100)
  endif

  icount = 0
  cur_waypoint => list%first
  do 
    if (.not.associated(cur_waypoint)) exit
    call WaypointPrint(cur_waypoint,option,output_option)
    icount = icount + 1
    cur_waypoint => cur_waypoint%next
  enddo

  if (OptionPrintToScreen(option)) then
    write(*,20) 'Total Waypoints:', icount
    write(*,100)
  endif

  if (OptionPrintToFile(option)) then
    write(option%fid_out,20) 'Total Waypoints:', icount
    write(option%fid_out,100)
  endif

end subroutine WaypointListPrint

! ************************************************************************** !

function WaypointListCopy(list)
  ! 
  ! Copies a waypoint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/19/13
  ! 

  use Option_module

  implicit none
  
  type(waypoint_list_type), pointer :: WaypointListCopy
  
  type(waypoint_list_type), pointer :: list
  type(waypoint_type), pointer :: new_waypoint
  type(waypoint_type), pointer :: prev_new_waypoint
  
  type(waypoint_list_type), pointer :: new_list
  type(waypoint_type), pointer :: cur_waypoint

  new_list => WaypointListCreate()
  
  nullify(prev_new_waypoint)
  
  cur_waypoint => list%first
  do 
    if (.not.associated(cur_waypoint)) exit
    new_waypoint => WaypointCreate(cur_waypoint)
    if (associated(prev_new_waypoint)) then
      prev_new_waypoint%next => new_waypoint
    else
      new_list%first => new_waypoint
    endif
    new_list%num_waypoints = new_list%num_waypoints + 1
    prev_new_waypoint => new_waypoint
    nullify(new_waypoint)
    cur_waypoint => cur_waypoint%next
  enddo
  
  WaypointListCopy => new_list

end function WaypointListCopy

! ************************************************************************** !

function WaypointForceMatchToTime(waypoint)
  ! 
  ! Forces a match to waypoint time if condition is
  ! true.
  ! 
  ! Author: Glenn Hammond
  ! Date: 03/19/13
  ! 
#include "petsc/finclude/petscsys.h"
  use petscsys
  implicit none
  
  type(waypoint_type) :: waypoint
  
  PetscBool :: WaypointForceMatchToTime
  
  WaypointForceMatchToTime = PETSC_FALSE

  if (waypoint%sync .or. &
      waypoint%update_conditions .or. &
      waypoint%print_snap_output .or. &
      waypoint%print_obs_output .or. &
      waypoint%print_msbl_output .or. &
      waypoint%print_checkpoint .or. &
      waypoint%final &
      ) then
    WaypointForceMatchToTime = PETSC_TRUE
  endif
  
end function WaypointForceMatchToTime


! ************************************************************************** !

function WaypointCreateSyncWaypointList(waypoint_list)
  !
  ! Creates a list of waypoints for outer synchronization of simulation process
  ! model couplers
  !
  ! Author: Glenn Hammond
  ! Date: 10/08/14
  !

  use Option_module

  implicit none

  type(waypoint_list_type), pointer :: waypoint_list

  type(waypoint_list_type), pointer :: WaypointCreateSyncWaypointList

  type(waypoint_list_type), pointer :: new_waypoint_list
  type(waypoint_type), pointer :: cur_waypoint
  type(waypoint_type), pointer :: new_waypoint

  new_waypoint_list => WaypointListCreate()

  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%sync .or. cur_waypoint%final) then
      new_waypoint => WaypointCreate(cur_waypoint)
      call WaypointInsertInList(new_waypoint,new_waypoint_list)
      if (cur_waypoint%final) exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  WaypointCreateSyncWaypointList => new_waypoint_list

end function WaypointCreateSyncWaypointList

! ************************************************************************** !

subroutine WaypointPrint(waypoint,option,output_option)
  ! 
  ! Prints a waypoint
  ! 
  ! Author: Glenn Hammond
  ! Date: 05/20/11
  ! 
  use Output_Aux_module
  use Option_module

  implicit none
  
  type(waypoint_type), pointer :: waypoint
  type(option_type) :: option
  type(output_option_type) :: output_option

  character(len=MAXSTRINGLENGTH) :: string

  10 format('  ',a20,':',10es13.5)
  30 format('  ',a20,':',10l)
  100 format(/)
  110 format(a)

  if (OptionPrintToScreen(option)) then
    write(*,110) 'Waypoint:'
    write(string,*) 'Time [' // trim(adjustl(output_option%tunit)) // ']'
    write(*,10) trim(string), waypoint%time/output_option%tconv
    write(*,30) 'Sync', waypoint%sync
    write(*,30) 'Print Snapshot Output', waypoint%print_snap_output
    write(*,30) 'Print Observation Output', waypoint%print_obs_output
    write(*,30) 'Print Mass Balance Output', waypoint%print_msbl_output
    write(*,30) 'Print Checkpoint', waypoint%print_checkpoint
    write(*,30) 'Update Conditions', waypoint%update_conditions
    write(string,*) 'Max DT [' // trim(adjustl(output_option%tunit)) // ']'
    write(*,10) trim(string), waypoint%dt_max/output_option%tconv
    write(*,30) 'Final', waypoint%final
    write(*,100)
  endif

  if (OptionPrintToFile(option)) then
    write(option%fid_out,110) 'Waypoint:'
    write(string,*) 'Time [' // trim(adjustl(output_option%tunit)) // ']'
    write(option%fid_out,10) trim(string), waypoint%time/output_option%tconv
    write(option%fid_out,30) 'Sync', waypoint%sync
    write(option%fid_out,30) 'Print Snapshot Output', waypoint%print_snap_output
    write(option%fid_out,30) 'Print Observation Output', &
                                                       waypoint%print_obs_output
    write(option%fid_out,30) 'Print Mass Balance Output', &
                                                      waypoint%print_msbl_output
    write(option%fid_out,30) 'Print Checkpoint', waypoint%print_checkpoint
    write(option%fid_out,30) 'Update Conditions', waypoint%update_conditions
    write(string,*) 'Max DT [' // trim(adjustl(output_option%tunit)) // ']'
    write(option%fid_out,10) trim(string), waypoint%dt_max/output_option%tconv
    write(option%fid_out,30) 'Final', waypoint%final
    write(option%fid_out,100)
  endif
 
end subroutine WaypointPrint

! ************************************************************************** !

subroutine WaypointInputRecord(output_option,waypoint_list)
  !
  ! Prints ingested time information to the input record file.
  !
  ! Author: Jenn Frederick
  ! Date: 05/09/2016
  !
  use Output_Aux_module
  
  implicit none
  
  type(output_option_type), pointer :: output_option
  type(waypoint_list_type), pointer :: waypoint_list
  
  type(waypoint_type), pointer :: cur_waypoint
  character(len=MAXWORDLENGTH) :: word1, word2
  character(len=MAXSTRINGLENGTH) :: string
  PetscReal :: final_time
  PetscReal :: max_dt
  PetscReal :: prev_time
  PetscInt :: id = INPUT_RECORD_UNIT
  character(len=10) :: Format
  
  Format = '(ES14.7)'
  
  write(id,'(a)') ' '
  write(id,'(a)') '---------------------------------------------------------&
                  &-----------------------'
  write(id,'(a29)',advance='no') '---------------------------: '
  write(id,'(a)') 'TIME'
  
  final_time = 0.d0
  prev_time = 0.d0
  max_dt = 0.d0
  
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%final .or. cur_waypoint%time > final_time) then
      final_time = cur_waypoint%time
    endif
    if (.not. Equal(cur_waypoint%dt_max,max_dt)) then
      write(id,'(a29)',advance='no') 'max. timestep: '
      write(word1,Format) cur_waypoint%dt_max/output_option%tconv
      write(word2,Format) prev_time/output_option%tconv
      write(id,'(a)') adjustl(trim(word1)) // ' ' // &
        trim(output_option%tunit) // ' at time ' // adjustl(trim(word2)) &
        // ' ' // trim(output_option%tunit)
    endif
    max_dt = cur_waypoint%dt_max
    prev_time = cur_waypoint%time
    cur_waypoint => cur_waypoint%next
  enddo
  
  write(id,'(a29)',advance='no') 'final time: '
  write(word1,Format) final_time/output_option%tconv
  write(id,'(a)') adjustl(trim(word1)) // ' ' // trim(output_option%tunit)

end subroutine WaypointInputRecord

! ************************************************************************** !

function WaypointListGetFinalTime(waypoint_list)
  ! 
  ! Returns the final time in the waypoint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 06/12/13
  ! 
  implicit none
  
  type(waypoint_list_type) :: waypoint_list
  
  PetscReal :: WaypointListGetFinalTime
  
  type(waypoint_type), pointer :: cur_waypoint

  ! initialize to negative infinity
  WaypointListGetFinalTime = -1.d20
  
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    if (cur_waypoint%final .or. &
        cur_waypoint%time > WaypointListGetFinalTime) then
      WaypointListGetFinalTime = cur_waypoint%time
      if (cur_waypoint%final) exit
    endif
    cur_waypoint => cur_waypoint%next
  enddo
  
end function WaypointListGetFinalTime 

! ************************************************************************** !

subroutine WaypointListDestroy(waypoint_list)
  ! 
  ! Destroys a simulation waypoint list
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/07/07
  ! 

  implicit none
  
  type(waypoint_list_type), pointer :: waypoint_list
  
  type(waypoint_type), pointer :: cur_waypoint, next_waypoint
  
  if (.not.associated(waypoint_list)) return
  
  cur_waypoint => waypoint_list%first
  do
    if (.not.associated(cur_waypoint)) exit
    next_waypoint => cur_waypoint%next
    call WaypointDestroy(cur_waypoint)
    cur_waypoint => next_waypoint
  enddo
  
  nullify(waypoint_list%first)
  nullify(waypoint_list%last)
  if (associated(waypoint_list%array)) deallocate(waypoint_list%array)
  nullify(waypoint_list%array)

  deallocate(waypoint_list)
  nullify(waypoint_list)
  
end subroutine WaypointListDestroy 

! ************************************************************************** !

subroutine WaypointDestroy(waypoint)
  ! 
  ! Deallocates a waypoint
  ! geh: DO NOT make this subroutine recursive as waypoints within lists need to
  ! be destroyed without recursively destroying the remainder of the list.
  ! 
  ! Author: Glenn Hammond
  ! Date: 11/09/07
  ! 

  implicit none
  
  type(waypoint_type), pointer :: waypoint
  
  if (.not.associated(waypoint)) return

  nullify(waypoint%prev)
  nullify(waypoint%next)
  deallocate(waypoint)
  nullify(waypoint)
  
end subroutine WaypointDestroy

end module Waypoint_module
