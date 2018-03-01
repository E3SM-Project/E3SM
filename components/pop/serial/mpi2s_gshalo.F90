
! SERIAL VERSION


!>
!! Routines to initialize and perform the boundary exchange operation.
!<
module mpi2s_gshalo
    use kinds_mod, only: i4,r8, log_kind
    implicit none
    private
    
    include 'mpif.h'
    
    integer(i4), parameter :: NUMexpectMSG  = 101, &
        expectMSG     = 102, &
        UpdateHaloMSG = 103
    
    integer(i4) :: my_task, nprocs
    
    !this is a communication package/schedule
    type, public :: Schedule_t
        integer(i4) :: COMM
        integer(i4) :: nRecv,nSend
        integer(i4) :: lenRecvBuffer,lenSendBuffer
        integer(i4), pointer :: sNeigh(:),ptrSend(:),SendCnt(:)
        integer(i4), pointer :: rNeigh(:),ptrRecv(:),RecvCnt(:)
        integer(i4), pointer :: Rrequest(:), Srequest(:)
        integer(i4), pointer :: Rstatus(:,:), Sstatus(:,:)
        integer(i4), pointer :: halo2send(:) !indir. addressing
        integer(i4), pointer :: recv2halo(:) !indir.  addressing
        real(r8)   , pointer :: bufferRecvDbl(:),bufferSendDbl(:)
        integer(i4), pointer :: bufferRecvInt(:),bufferSendInt(:)

    end type
    
   !this is to keep track of which groups a process belongs to
   !who their neighbors are, etc. 
   type, public :: Groups_t 
      logical (log_kind):: init
      integer(i4) :: numTotal  ! total num overflows, estuaries, etc.
      integer(i4) :: numMyGroups !groups that I am active in
      integer(i4), pointer :: groupIds(:) !this will be size of numMyGroups
      logical (log_kind), pointer :: amMaster(:) ! master of this group
      integer(i4), pointer :: neighbors(:) !these are the neighbors
                                                 ! procs in each of my groups
      integer(i4), pointer:: neighborStarts(:) ! to index into neighbors
                                                   ! size (numMyGroups+1)
      type(Schedule_t), pointer :: commHandle(:) !comm handle schedule (size numMyGroups)
   end type Groups_t



    integer(i4) :: nPrint
    
    public :: mpi2s_gshalo_update_2d_dbl
    public :: mpi2s_gshalo_update_2d_int
    public :: mpi2s_gshalo_init
    public :: mpi2s_gshalo_init_symm_same
    public :: mpi2s_gshalo_global_sum_dbl
    public :: LinearOrderedFind

    integer(i4), public :: timer_mpi2s_halo_init, &
         timer_mpi2s_update_dbl

  contains
    

!***********************************************************************
! !IROUTINE: mpi2s_gshalo_update_2d_dbl

    subroutine mpi2s_gshalo_update_2d_dbl(handle, array, mytag)

! !DESCRIPTION:


! !INPUT PARAMETERS:

      type (Schedule_t), intent(in) :: handle
      real(r8), intent(inout) :: array(:)
      integer(i4), intent(in), optional :: mytag !tag for messages

      !do nothing - array has the correct info
      
      return

    end subroutine mpi2s_gshalo_update_2d_dbl

    
!***********************************************************************
! !IROUTINE: mpi2s_gshalo_update_2d_int

    subroutine mpi2s_gshalo_update_2d_int(handle, array, mytag)

! !DESCRIPTION:


! !INPUT PARAMETERS:
      type (Schedule_t), intent(in) :: handle
      integer(i4), intent(inout) :: array(:)
      integer(i4), intent(in), optional :: mytag !tag for messages
      
      
      !do nothing - array has the correct info

      return


    end subroutine mpi2s_gshalo_update_2d_int


!***********************************************************************
! !IROUTINE: mpi2s_gshalo_global_sum_dbl

    subroutine mpi2s_gshalo_global_sum_dbl (handle, array, a_size, tag)


! !DESCRIPTION:
      !! computes the global sum(s) according to the comm handle
      !! (so using p2p communication) - The handle for the global
      !! sum shoudl be set with the mpi2s_gshalo_init_symm_same

      !! @param handle   Communication handle of type Schedule_t
      !! @array          Array with local sums (will be overwritten
      !!                 with global sums)

      !! @a_size         Number of local sums in array (it may be less
      !!                 then the allocated send buff size - though if
      !!                 it is a lot less, then you should have created 
      !!                 another handle)
      !! @tag            Optional message tag
      
      
! !INPUT PARAMS

      type (Schedule_t), intent(in) :: handle
      real(r8), intent(inout) :: array(:)
      integer(i4), intent (in) :: a_size
      integer(i4), intent(in), optional::  tag


!      print *,'Processor: ',my_task,' in GLOBAL_SUM_DBL'


      !do nothing - array has the correct info

      return

    end subroutine mpi2s_gshalo_global_sum_dbl

!***********************************************************************
! !IROUTINE: mpi2s_gshalo_init_symm_same

    function mpi2s_gshalo_init_symm_same (COMM, nNbors, &
         nTotal, nborProcs, int_msg) result(schedHandle)

! !DESCRIPTION:
      !! Construct the message passing metadata (Schedule_t) necessary to 
      !! perform the MPI based exchange of data in halo region
      !! Here the halo region is symmetric (sends = rcvs), and we 
      !! are sending the same
      !! data to all - useful for the replacement of a global sum when few
      !! procs are involved (So we DO NOT assume that this routine is called
      !! collectively by the COMM procs) 

      !!
      !! @param COMM        MPI communicator
      !! @param nNbors      Number of recv (=send) neighbors
      !!                    
      !! @param nTotal      Number of items (integers) to send to each
      !!                    neighbor
      !! @param nborProcs   A 1D array which contains the MPI rank or task for
      !!                    each neighbor (size nNbors)
      !! @param int_msg     logical, true = sending ints only, false  = real
      !!                    only, both if argument not provided

! !INPUT PARAMETERS:

        integer(i4), intent(in) :: COMM
        integer(i4), intent(in) :: nNbors
        integer(i4), intent(in) :: nTotal
        integer(i4), intent(in) :: nborProcs(nNbors)
        logical(log_kind), intent(in), optional :: int_msg 

! !OUTPUT PARAMETERS:

        type (Schedule_t):: schedHandle

! !local variables



        schedHandle%COMM          = COMM
        schedHandle%nRecv        =  0
        schedHandle%nSend        =  0
        schedHandle%lenSendBuffer = 0
        schedHandle%lenRecvBuffer = 0

           

end function mpi2s_gshalo_init_symm_same


!***********************************************************************
! !IROUTINE: mpi2s_gshalo_init


    function mpi2s_gshalo_init( &
        COMM, maxlinear, nTotal, nActive, LinearGdof, LinearProc, int_msg) &
        result(schedHandle)

! !DESCRIPTION:
      !! Construct the message passing metadata necessary to perform the MPI based
      !! boundary exchanges necessary to update the halo region. 
      !!
      !! @param COMM        MPI communicator
      !! @param maxlinear   The maximum number of gridpoints in the 1D data
      !!                    structure.
      !! @param nTotal      Total number of gridpoints including active and halo
      !!                    points.
      !! @param nActive     The total number of active ocean gridpoints.
      !! @param LinearGdof  A 1D array which contains the global degrees of
      !!                    freedom for this MPI task.
      !! @param LinearProc  A 1D array which contains the MPI rank or task for
      !!                    each gridpoint.
      !! @param int_msg     logical, true = sending ints only, false  = real
      !!                    only, both if argument not provided


! !INPUT PARAMETERS:

        integer(i4), intent(in) :: COMM
        integer(i4), intent(in) :: maxlinear
        integer(i4), intent(in) :: nTotal,nActive
        integer(i4), intent(in) :: LinearGdof(maxlinear)
        integer(i4), intent(in) :: LinearProc(maxlinear)
        logical(log_kind), intent(in), optional :: int_msg 


! !OUTPUT PARAMETERS:

        !! integer(i4) :: schedHandle
        type (Schedule_t):: schedHandle

        
        schedHandle%nRecv        =  0
        schedHandle%nSend        =  0

        schedHandle%COMM          = COMM
        schedHandle%lenSendBuffer = 0
        schedHandle%lenRecvBuffer = 0


    end function mpi2s_gshalo_init

    !***********************************************************************
    !>
    !! This subroutine inserts a unique value into an array in sortered order. 
    !! It also maintains associated satellite data.
    !!
    !! @param array   The array into which a value will be inserted.
    !! @param cnt     Satellite data associated with the array.
    !! @param value   A scalar value for insertion into the array 
    !! @param cvalue  The satellite data for insertion.
    !! @param ierr    Error return code which is not set. [REMOVE]
    !<
    subroutine InsertIntoArray(array,cnt, value,cvalue,ierr)
        integer(i4), intent(inout) :: array(:), cnt(:)
        integer(i4), intent(in)    :: value, cvalue
        integer(i4), intent(out) :: ierr
        integer(i4), allocatable :: tmp(:), ctmp(:)
        logical  :: found
        integer(i4)  :: n,indx
        logical, parameter  :: Debug = .FALSE.

        call LinearOrderedFind(array,value,found,indx)

        if(.not. found) then
            n = SIZE(array)
            allocate(tmp(n),ctmp(n))
            tmp = array
            ctmp  = cnt

            if(indx== 1) then
                array(1)   = value
                cnt(1)     = cvalue
                array(2:n) = tmp(1:n-1)
                cnt(2:n)   = ctmp(1:n-1)
            else
                array(1:indx-1) = tmp(1:indx-1)
                cnt(1:indx-1)   = ctmp(1:indx-1)
                array(indx)     = value
                cnt(indx)       = cvalue
                array(indx+1:n) = tmp(indx:n-1)
                cnt(indx+1:n)   = ctmp(indx:n-1)
            endif

            deallocate(tmp,ctmp)
        else
            cnt(indx) = cnt(indx) + cvalue
        endif
    end subroutine InsertIntoArray

    !***********************************************************************
    !>
    !! This subroutine determines if a scalar value already exists in an array.
    !!
    !! @param array  A 1D integer array to be searched.
    !! @param value  The scalar value for which to search.
    !! @param found  A logical value which indicates the scalar value is found.
    !! @param indx   If found=.true. the location of the scalar in the array.
    !>
    subroutine LinearFind(array,value,found,indx)
        integer(i4), intent(in) :: array(:)
        integer(i4), intent(in) :: value
        logical, intent(out) :: found
        integer(i4), intent(out) :: indx
        integer(i4) :: i, n

        n = SIZE(array)
        found = .FALSE.
        do i=1,n
            if(array(i) == value) then 
                found = .TRUE.
                indx = i
                return
            endif
        enddo
    end subroutine LinearFind

    !***********************************************************************
    !>
    !! This subroutine determines if a scalar value already exists in a sorted
    !! array.
    !!
    !! @param array  A 1D integer array to be searched.
    !! @param value  The scalar value for which to search.
    !! @param found  A logical value which indicates the scalar value is found.
    !! @param indx   If found=.true. the location of the scalar in the array.
    !>
    subroutine LinearOrderedFind(array,value,found,indx)
        integer(i4), intent(in)  :: array(:)
        integer(i4), intent(in)  :: value
        logical, intent(out) :: found
        integer(i4), intent(out) :: indx

        integer(i4) :: n,nz,i
        n = SIZE(array)

        ! =============================================
        ! Array of all zeros... insert at the begining
        ! =============================================
        found = .FALSE.
        if((array(1) < 0 ) .or. (value < array(1)) ) then
            found = .FALSE.
            indx = 1
            return
        endif
        nz = COUNT(array .ge. 0)
        do i=1,nz
            if(array(i) == value) then
                ! ===============
                ! Already in list
                ! ===============
                found = .TRUE.
                indx = i
                return
            else if((array(i) < value) .and. (value < array(i+1)) ) then
                ! =====================================
                ! Insert it into the middle of the list
                ! =====================================
                found = .FALSE.
                indx = i+1
                return
            else if((array(i) < value) .and. (array(i+1) == -1) )then
                ! =====================================
                ! Insert it into the end of the List
                ! =====================================
                found = .FALSE.
                indx = i+1
                return
            endif
        enddo
    end subroutine LinearOrderedFind
end module mpi2s_gshalo
