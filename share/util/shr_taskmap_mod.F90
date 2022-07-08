module shr_taskmap_mod
!-----------------------------------------------------------------------
!
! Purpose:
! Output mapping of MPI tasks to nodes for a specified
! communicator
!
! Methods: 
!  Use mpi_get_processor_name to identify the node that an MPI 
!  task for a given communicator is assigned to. Gather these
!  data to task 0 and then write out the list of MPI 
!  tasks associated with each node using the designated unit
!  number
!
! Author: P. Worley 
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_sys_mod, only: shr_sys_abort

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   include 'mpif.h'
   private
   save                             ! Make the default access private

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public :: &
      shr_taskmap_write             ! write out list of nodes
                                    !  with list of assigned MPI tasks
                                    !  for a given communicator

   CONTAINS

!
!========================================================================
!
   subroutine shr_taskmap_write (unit_num, comm_id, comm_name, &
                                 verbose, no_output, &
                                 save_nnodes, save_task_node_map)

!-----------------------------------------------------------------------
! Purpose: Write out list of nodes used by processes in a given
!          communicator. For each node output the list of MPI tasks
!          assigned to it.
! Author: P. Worley
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
      integer, intent(in)      :: unit_num   ! unit number for output
      integer, intent(in)      :: comm_id    ! MPI communicator
      character(*), intent(in) :: comm_name  ! MPI communicator label

      logical, intent(in), optional :: verbose
                                             ! verbose output flag
                                             !  (Default is .false.)
      logical, intent(in), optional :: no_output
                                             ! no output flag
                                             !  (Default is .false.)
      integer, intent(out), optional :: save_nnodes
                                             ! return number of nodes
      integer, intent(out), optional :: save_task_node_map(:)
                                             ! return task-to-node map

!---------------------------Local Workspace-----------------------------
      integer :: iam                    ! task id in comm_id
      integer :: npes                   ! number of MPI tasks in comm_id
      integer :: ier                    ! return error status    
      integer :: max_len                ! maximum name length
      integer :: length                 ! node name length
      integer :: c, i, j                ! loop indices
      integer :: nnodes                 ! number of nodes
      integer :: start, limit           ! loop bounds
      integer :: head, tail             ! limits of current sequential run
                                        !  of task ids

      ! flag to indicate whether returning number of nodes
      logical :: broadcast_nnodes

      ! flag to indicate whether returning task-to-node mapping
      logical :: broadcast_task_node_map

      ! flag to indicate whether to use verbose or compact output format
      logical :: verbose_output

      ! flag to indicate whether to write out information
      ! (for when want to calculate nnodes and the task_node_map without 
      !  output)
      logical :: output

      ! mapping of tasks to ordered list of nodes
      integer, allocatable :: task_node_map(:)

      ! number of MPI tasks per node
      integer, allocatable :: node_task_cnt(:)    
      integer, allocatable :: node_task_tmpcnt(:)

      ! MPI tasks ordered by nodes to which they are assigned
      integer, allocatable :: node_task_map(:)    

      ! offset into node_task_map for processes assigned to given node
      integer, allocatable :: node_task_offset(:)

      logical :: masterproc             ! masterproc flag
      logical :: done                   ! search completion flag

      ! node names for each mpi task
      character(len=mpi_max_processor_name) :: tmp_name
      character, allocatable :: task_node_name(:)  ! for this task
      character, allocatable :: task_node_names(:) ! for all tasks

      ! node names without duplicates
      character(len=mpi_max_processor_name), allocatable :: node_names(:)  

      ! string versions of numerical values
      character(len=8) :: c_npes       ! number of MPI tasks
      character(len=8) :: c_nnodes     ! number of nodes
      character(len=8) :: c_nodeid     ! node id
      character(len=8) :: c_node_npes  ! number of MPI tasks for a given node
      character(len=8) :: c_taskid     ! MPI task id

      ! routine name, for error reporting
      character(*),parameter :: subname = "(shr_taskmap_write)"

!-----------------------------------------------------------------------
      !
      ! Get my id  
      !
      call mpi_comm_rank (comm_id, iam, ier) 
      if (iam == 0) then 
         masterproc = .true.
      else
         masterproc = .false.
      end if

      !
      ! Get number of MPI tasks
      !
      call mpi_comm_size (comm_id, npes, ier)

      !
      ! Determine whether to use verbose output format
      !
      verbose_output = .false.
      if (present(verbose)) then
         verbose_output = verbose
      endif

      !
      ! Determine whether to output taskmap
      !
      output = .true.
      if (present(no_output)) then
         if (no_output) output = .false.
      endif

      !
      ! Determine whether returning number of nodes
      !
      broadcast_nnodes = .false.
      if (present(save_nnodes)) then
         broadcast_nnodes = .true.
      endif

      !
      ! Determine whether returning task-to-node mapping information
      !
      broadcast_task_node_map = .false.
      if (present(save_task_node_map)) then
         if (size(save_task_node_map) >= npes) then
            broadcast_task_node_map = .true.
         else
            call shr_sys_abort(trim(subname)//': array for task-to-node mapping data too small')
         endif
      endif

      ! 
      ! Allocate arrays for collecting node names
      !
      max_len = mpi_max_processor_name
      allocate ( task_node_name(max_len), stat=ier )
      if (ier /= 0) &
         call shr_sys_abort(trim(subname)//': allocate task_node_name failed')

      allocate ( task_node_names(max_len*npes), stat=ier )
      if (ier /= 0) &
         call shr_sys_abort(trim(subname)//': allocate task_node_names failed')
 
      !
      ! Get node names and send to root. 
      ! (Assume that processor names are node names.)
      !
      call mpi_get_processor_name (tmp_name, length, ier)
      task_node_name(:) = ' '
      do i = 1, length
         task_node_name(i) = tmp_name(i:i)
      end do

      !
      ! Gather node names
      !
      task_node_names(:) = ' '
      call mpi_gather (task_node_name,  max_len, mpi_character, &
                       task_node_names, max_len, mpi_character, &
                       0, comm_id, ier)

      if (masterproc) then
         !
         ! Identify nodes and task/node mapping.
         !
         allocate ( task_node_map(0:npes-1), stat=ier )
         if (ier /= 0) &
            call shr_sys_abort(trim(subname)//': allocate task_node_map failed')
         task_node_map(:) = -1

         allocate ( node_names(0:npes-1), stat=ier )
         if (ier /= 0) &
            call shr_sys_abort(trim(subname)//': allocate node_names failed')
         node_names(:) = ' '

         allocate ( node_task_cnt(0:npes-1), stat=ier )
         if (ier /= 0) &
            call shr_sys_abort(trim(subname)//': allocate node_task_cnt failed')
         node_task_cnt(:) = 0

         do c=1,max_len
            tmp_name(c:c) = task_node_names(c)
         enddo

         node_names(0) = trim(tmp_name)
         task_node_map(0) = 0
         node_task_cnt(0) = 1
         nnodes = 1

         do i=1,npes-1
            do c=1,max_len
               tmp_name(c:c) = task_node_names(i*max_len+c)
            enddo

            j = 0
            done = .false.
            do while ((.not. done) .and. (j < nnodes))
               if (trim(node_names(j)) .eq. trim(tmp_name)) then
                  task_node_map(i) = j
                  node_task_cnt(j) = node_task_cnt(j) + 1
                  done = .true.
               endif
               j = j + 1
            enddo

            if (.not. done) then
               node_names(nnodes) = trim(tmp_name)
               task_node_map(i) = nnodes
               node_task_cnt(nnodes) = 1
               nnodes = nnodes + 1
            endif

         enddo

         !
         ! Identify node/task mapping.
         !
         allocate ( node_task_offset(0:nnodes-1), stat=ier )
         if (ier /= 0) &
            call shr_sys_abort(trim(subname)//': allocate node_task_offset failed')
         node_task_offset(:) = 0

         do j=1,nnodes-1
            node_task_offset(j) = node_task_offset(j-1) + node_task_cnt(j-1)
         enddo

         allocate ( node_task_tmpcnt(0:nnodes-1), stat=ier )
         if (ier /= 0) &
            call shr_sys_abort(trim(subname)//': allocate node_task_tmpcnt failed')
         node_task_tmpcnt(:) = 0

         allocate ( node_task_map(0:npes-1), stat=ier )
         if (ier /= 0) &
            call shr_sys_abort(trim(subname)//': allocate node_task_map failed')
         node_task_map(:) = -1

         do i=0,npes-1
            j = task_node_map(i)
            node_task_map(node_task_offset(j) + node_task_tmpcnt(j)) = i
            node_task_tmpcnt(j) = node_task_tmpcnt(j) + 1
         enddo

         if (output) then
            !
            ! Output node/task mapping
            !
            write(unit_num,100) &
               '--------------------------------------------------------------'
100         format(a)

            write(c_npes,'(i8)') npes
            write(c_nnodes,'(i8)') nnodes
            write(unit_num,101) trim(comm_name), trim(adjustl(c_nnodes)), &
                                trim(adjustl(c_npes))
101         format(a,' communicator : ',a,' nodes, ',a,' MPI tasks')

            write(unit_num,100) &
               'COMMUNICATOR NODE # [NODE NAME] : (# OF MPI TASKS) TASK # LIST'

            do j=0,nnodes-1
               write(c_nodeid,'(i8)') j
               write(c_node_npes,'(i8)') node_task_cnt(j)
               write(unit_num,102,advance='no') &
                  trim(comm_name), trim(adjustl(c_nodeid)), &
                  trim(node_names(j)), trim(adjustl(c_node_npes))
102            format(a,' NODE ',a,' [ ',a,' ] : ( ',a,' MPI TASKS )')

               start = node_task_offset(j)
               limit = start+node_task_cnt(j)-1

               if (verbose_output) then

                  do i=start,limit
                     write(c_taskid,'(i8)') node_task_map(i)
                     write(unit_num,103,advance='no') trim(adjustl(c_taskid))
103                  format(' ',a)
                  enddo

               else

                  head  = node_task_map(start)
                  tail  = head
                  do i=start+1,limit
                     if (node_task_map(i) == tail+1) then
                        tail = tail + 1
                     else
                        write(c_taskid,'(i8)') head
                        write(unit_num,103,advance='no') trim(adjustl(c_taskid))
                        if (head /= tail) then
                           write(c_taskid,'(i8)') tail
                           write(unit_num,104,advance='no') trim(adjustl(c_taskid))
104                        format('-',a)
                        endif
                        head = node_task_map(i)
                        tail = head
                     endif
                  enddo

                  if (node_task_map(limit) == tail) then
                     write(c_taskid,'(i8)') head
                     write(unit_num,103,advance='no') trim(adjustl(c_taskid))
                     if (head /= tail) then
                        write(c_taskid,'(i8)') tail
                        write(unit_num,104,advance='no') trim(adjustl(c_taskid))
                     endif
                  endif

               endif

               write(unit_num,105,advance='no')
105            format(/)
            enddo
            write(unit_num,100) &
               '--------------------------------------------------------------'
         endif

         if (broadcast_nnodes) then
            save_nnodes = nnodes
         endif

         if (broadcast_task_node_map) then
            do i=0,npes-1
               save_task_node_map(i+1) = task_node_map(i)
            enddo
         endif

         deallocate(node_task_map)
         deallocate(node_task_tmpcnt)
         deallocate(node_task_offset)
         deallocate(node_task_cnt)
         deallocate(node_names)
         deallocate(task_node_map)

      endif

      if (broadcast_nnodes) then
         call mpi_bcast(save_nnodes, 1, mpi_integer, 0, comm_id, ier)
      endif

      if (broadcast_task_node_map) then
         call mpi_bcast(save_task_node_map, npes, mpi_integer, 0, comm_id, ier)
      endif

      deallocate(task_node_name)
      deallocate(task_node_names)

   end subroutine shr_taskmap_write

!
!========================================================================
!
end module shr_taskmap_mod
