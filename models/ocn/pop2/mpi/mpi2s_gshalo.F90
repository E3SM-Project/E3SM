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


! !local
      integer(i4) :: lenRecvBuffer,lenSendBuffer
      integer(i4) :: nSend,nRecv
      
      integer(i4) :: src,dest,len,iptr,tag
      
      integer(i4) :: ierr,i
      
      nRecv = handle%nRecv
      lenRecvBuffer=handle%lenRecvBuffer
      
      nSend = handle%nSend
      lenSendBuffer=handle%lenSendBuffer
      
      if (present(mytag)) then
         tag = mytag
      else
         tag = UpdateHaloMSG
      end if

      do i=1,nRecv
         len = handle%RecvCnt(i)
         iptr = handle%ptrRecv(i)
         src  = handle%rNeigh(i)
         call MPI_Irecv(handle%bufferRecvDbl(iptr),len, MPI_REAL8, src, &
              tag, handle%COMM, handle%Rrequest(i), ierr)
      enddo
      
      !=====================================================================
      ! Indirect address into the 1D data structure to create the SendBuffer
      !=====================================================================
      do i=1,lenSendBuffer
         handle%BufferSendDbl(i) = array(handle%halo2send(i))
      enddo
      
      do i=1,nSend
         len  = handle%SendCnt(i)
         iptr = handle%ptrSend(i)
         dest  = handle%sNeigh(i)
         call MPI_Isend(handle%bufferSendDbl(iptr), &
              len, MPI_REAL8, dest, &
              tag, handle%COMM,handle%Srequest(i), ierr)
      enddo
      if(nSend>0) call MPI_Waitall(nSend, handle%Srequest, &
           handle%Sstatus, ierr)
      if(nRecv>0) call MPI_Waitall(nRecv, handle%Rrequest, &
           handle%Rstatus,ierr)
      
      !===============================================================
      ! Indirect address from the RecvBuffer to the 1D data structure
      !===============================================================
      do i=1,lenRecvBuffer
         array(handle%recv2halo(i))=handle%bufferRecvDbl(i)
      enddo
    end subroutine mpi2s_gshalo_update_2d_dbl

    
!***********************************************************************
! !IROUTINE: mpi2s_gshalo_update_2d_int

    subroutine mpi2s_gshalo_update_2d_int(handle, array, mytag)

! !DESCRIPTION:


! !INPUT PARAMETERS:
      type (Schedule_t), intent(in) :: handle
      integer(i4), intent(inout) :: array(:)
      integer(i4), intent(in), optional :: mytag !tag for messages
      
! !local
      integer(i4) :: lenRecvBuffer,lenSendBuffer
      integer(i4) :: nRecv,nSend
      
      integer(i4) :: src,dest,len,iptr,tag
      
      integer(i4) :: ierr,i
      
      logical, parameter :: Debug = .FALSE.
      
      
      nRecv = handle%nRecv
      lenRecvBuffer=handle%lenRecvBuffer
      nSend = handle%nSend
      lenSendBuffer=handle%lenSendBuffer
      
      if(Debug) &
           print *,'IAM: ', my_task, &
           ' GSHalo_update_2d_int: halo2send(1:lenSendBuffer) :=', &
           handle%halo2send(1:lenSendBuffer)
      if(Debug) &
           print *,'IAM: ', my_task, &
           ' GSHalo_update_2d_int: recv2halo(1:lenRecvBuffer) :=', &
           handle%recv2halo(1:lenRecvBuffer)
      
      !=====================================================================
      ! Indirect address into the 1D data structure to create the SendBuffer
      !=====================================================================
      do i=1,lenSendBuffer
         handle%bufferSendInt(i) = array(handle%halo2send(i))
      enddo
      do i=1,nSend
         iptr = handle%ptrSend(i)
         len  = handle%SendCnt(i)
         dest = handle%sNeigh(i)
         if(Debug) &
              print *,'IAM: ',my_task, ' to: ',dest, &
              ' BufferSendInt(:):= ', handle%bufferSendInt(iptr:iptr+len-1)
      enddo
      if (present(mytag)) then
         tag = mytag
      else
         tag = UpdateHaloMSG
      end if
      do i=1,nSend
         len  = handle%SendCnt(i)
         iptr = handle%ptrSend(i)
         dest  = handle%sNeigh(i)
         if(Debug) &
              print *,'IAM: ',my_task,'Posting Send of length ',len, &
              ' to: ',dest
         call MPI_Isend(handle%bufferSendInt(iptr), &
              len, MPI_INTEGER, dest, &
              tag, handle%COMM, handle%Srequest(i), ierr)
      enddo
      
      do i=1,nRecv
         len = handle%RecvCnt(i)
         iptr = handle%ptrRecv(i)
         src  = handle%rNeigh(i)
         if(Debug) &
              print *,'IAM: ',my_task,'Posting Recv of length ',len, &
              ' from: ',src
         call MPI_IRecv(handle%bufferRecvInt(iptr), len, &
              MPI_INTEGER,src, &
              tag, handle%COMM, handle%Rrequest(i), ierr)
      enddo
      if(nSend>0) &
           call MPI_Waitall(nSend,handle%Srequest,handle%Sstatus,ierr)
      if(nRecv>0) &
           call MPI_Waitall(nRecv,handle%Rrequest,handle%Rstatus,ierr)
      
      do i=1,nRecv
         iptr = handle%ptrRecv(i)
         len  = handle%RecvCnt(i)
         src  = handle%rNeigh(i)
         if(Debug) print *,'IAM: ',my_task, ' From: ',src, &
              ' BufferRecvInt(:):= ', handle%bufferRecvInt(iptr:iptr+len-1)
      enddo
      
      !===============================================================
      ! Indirect address from the RecvBuffer to the 1D data structure
      !===============================================================
      do i=1,lenRecvBuffer
         array(handle%recv2halo(i))=handle%bufferRecvInt(i)
      enddo

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

! !local variables      

      real(r8), dimension(:), allocatable :: local_sums
      integer(i4) :: nSends, ls_size, i, j
      integer(i4) :: lenSendBuffer
      real (r8) :: temp_sum
      integer (i4) :: nbor_id
      logical (log_kind) :: mysum_added

!      print *,'Processor: ',my_task,' in GLOBAL_SUM_DBL'

      lenSendBuffer=handle%lenSendBuffer
      nSends =  handle%nSend
      
      if (nSends == 0) then
         return
      end if
      ls_size = (nSends+1)*lenSendBuffer !hold my local sums + 
                                         ! local sums from my neighbors
      
      allocate(local_sums(ls_size))
      !initialize
      local_sums = 0.0
      !now put my local sum(s) in the first position(s)
      local_sums(1:a_size) = array(1:a_size)
      
!      print *,'GS BEGIN Processor: ',my_task,' ARRAY = ', array(:)

!      print *,'GS Processor: ',my_task,' a_size = ', a_size, 'local sums BEFORE: ',local_sums(:)

      !now get my neighbors local sums 
      if (present(tag)) then
         call mpi2s_gshalo_update_2d_dbl(handle, local_sums, tag)
      else
         call mpi2s_gshalo_update_2d_dbl(handle, local_sums)
      end if

!      print *,'GS Processor: ',my_task,' local sums AFTER: ',local_sums(:)


      !now compute the global sums and put back into the array
      !recall that lenSendBuffer is total info size, but we only
      !care about the first a_size items (rest should be zero)
      ! note that the first local_sums (j=0 below) are from this processor
      ! do the sums in order of ascending proc number so that all
      ! procs get the same answer

      do i = 1, a_size
         temp_sum = 0.0
         mysum_added = .false.
         do j = 1, nSends
            nbor_id = handle%sNeigh(j)
            if (nbor_id < my_task .or. mysum_added) then
               temp_sum = temp_sum  + local_sums(j*lenSendBuffer + i) !nbr
            else !need to add local proc's contribution first
               temp_sum = temp_sum + local_sums(i) !local (j=0)
               temp_sum = temp_sum  + local_sums(j*lenSendBuffer + i)!nbr 
               mysum_added = .true.
            endif
         end do !j
         if (mysum_added) then
            array(i) = temp_sum 
         else
            temp_sum = temp_sum + local_sums(i) 
            array(i) = temp_sum
         endif

      end do !i

!      print *,'GS END Processor: ',my_task,' ARRAY = ', array(:)


      !now array contains the a_size global sums

      deallocate(local_sums)

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

        integer(i4) :: i, ierr
        integer(i4) :: nSend,nRecv
        integer(i4), allocatable :: ptrRecv(:)
        integer(i4), allocatable :: RecvCnt(:)
        integer(i4), allocatable :: ptrSend(:)
        integer(i4), allocatable :: rNeigh(:), sNeigh(:)
        integer(i4), allocatable :: SendCnt(:)

        integer(i4) :: lenRecvBuffer,lenSendBuffer


        ! ==============================
        ! Get some MPI information 
        ! ==============================
        call MPI_COMM_SIZE(COMM,nprocs,ierr)
        call MPI_COMM_RANK(COMM,my_task,ierr)
 
!        print *, 'HALO INIT SYMM SAME: IAM: ', my_task
        if (nNbors == 0) then

           schedHandle%COMM          = COMM
           schedHandle%nRecv        =  0
           schedHandle%nSend        =  0
           schedHandle%lenSendBuffer = 0
           schedHandle%lenRecvBuffer = 0
           return
           
        end if

        !allocate
        allocate(sNeigh(nNbors))
        allocate(rNeigh(nNbors))
        allocate(RecvCnt(nNbors))
        allocate(SendCnt(nNbors))
        allocate(ptrRecv(nNbors))
        allocate(ptrSend(nNbors))

        
        !symmetric comm:  #send = #recvs = #neighbors
        nSend=nNbors
        nRecv=nNbors

        rNeigh(1:nNbors) = nborProcs(1:nNbors)
        sNeigh(1:nNbors) = nborProcs(1:nNbors)

        ! send the SAME info (length nTotal) to each nbor
        lenSendBuffer= nTotal
        lenRecvBuffer=nTotal*nNbors

        ptrSend(1:nNbors) = 1
        do i = 1, nNbors
           ptrRecv(i) = nTotal*(i-1) + 1
           RecvCnt(i) = nTotal
           SendCnt(i) = nTotal
        enddo

        ! Pack all the info into the Schedule_t structure
        schedHandle%COMM          = COMM

        schedHandle%nRecv        =  nRecv
        schedHandle%nSend        =  nSend

        schedHandle%lenSendBuffer = lenSendBuffer
        schedHandle%lenRecvBuffer = lenRecvBuffer

        allocate(schedHandle%sNeigh(nSend))
        schedHandle%sNeigh(1:nSend) = sNeigh(1:nSend)
!        print *,'Processor: ',my_task,' has S-neighbors: ',sNeigh(1:nSend)
        deallocate(sNeigh)

        allocate(schedHandle%ptrSend(nSend))
        schedHandle%ptrSend(1:nSend) = ptrSend(1:nSend)
        deallocate(ptrSend)

        allocate(schedHandle%SendCnt(nSend))
        schedHandle%SendCnt(1:nSend)  = SendCnt(1:nSend)
        deallocate(SendCnt)

        allocate(schedHandle%rNeigh(nRecv))
        schedHandle%rNeigh(1:nRecv) = rNeigh(1:nRecv)
!        print *,'Processor: ',my_task,' has R-neighbors: ',rNeigh(1:nRecv)
        deallocate(rNeigh)

        allocate(schedHandle%ptrRecv(nRecv))
        schedHandle%ptrRecv(1:nRecv) = ptrRecv(1:nRecv)
        deallocate(ptrRecv)

        allocate(schedHandle%RecvCnt(nRecv))
        schedHandle%RecvCnt(1:nRecv)=RecvCnt(1:nRecv)
        deallocate(RecvCnt)

        ! allocate status and request arrays
        allocate(schedHandle%Srequest(nSend), schedHandle%Rrequest(nRecv))
        allocate(schedHandle%Sstatus(MPI_STATUS_SIZE, nSend), &
            schedHandle%Rstatus(MPI_STATUS_SIZE, nRecv))

        ! allocate the message buffers
        if (present(int_msg)) then
           if (int_msg) then
              allocate(schedHandle%bufferRecvInt(lenRecvBuffer))
              allocate(schedHandle%bufferSendInt(lenSendBuffer))
           else
              allocate(schedHandle%bufferRecvDbl(lenRecvBuffer))
              allocate(schedHandle%bufferSendDbl(lenSendBuffer))
           end if
        else
           allocate(schedHandle%bufferRecvDbl(lenRecvBuffer))
           allocate(schedHandle%bufferSendDbl(lenSendBuffer))
           
           allocate(schedHandle%bufferRecvInt(lenRecvBuffer))
           allocate(schedHandle%bufferSendInt(lenSendBuffer))
        end if

        ! allocate and set halo2send and recv2halo 
        !(indirect addressing arrays)
        allocate(schedHandle%halo2send(lenSendBuffer))
        allocate(schedHandle%recv2halo(lenRecvBuffer))


        !the array passed in will contain the send items first and
        !then have room for all nTotal recv items from each Nbor
        do i = 1,lenSendBuffer
           schedHandle%halo2send(i) = i
        end do
        do i = 1,lenRecvBuffer
           schedHandle%recv2halo(i) = lenSendBuffer + i
        end do

        ! TO DO:  there should be a destroy function to get 
        ! rid of this storage that has been allocated
        ! that would be called in ovf_destory (DNE yet)


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

! !local

        integer(i4) :: sWords,rWords
        integer(i4) :: sWords_min,sWords_max,sWords_avg,sWords_total

        integer(i4) :: nNeigh
        integer(i4) :: nNeigh_max,nNeigh_min,nNeigh_avg

        integer(i4) :: dest,src,cnt

        integer(i4) :: nSend,nRecv
        integer(i4) :: i,lenHalo,maxNeigh,ierr

        integer(i4), allocatable :: ptrRecv(:)
        integer(i4), allocatable :: RecvCnt(:)

        integer(i4), allocatable :: ptrSend(:)

        integer(i4), allocatable :: rNeigh(:), sNeigh(:)
        integer(i4), allocatable :: SendCnt(:)
        integer(i4), allocatable :: SendTMP(:)

        integer(i4), allocatable :: ptrCnt(:)
        integer(i4) :: lenRecvBuffer,lenSendBuffer
        integer(i4) :: tag
        integer(i4) :: iptr

        logical, parameter :: Debug = .FALSE.
        logical, parameter :: Info = .FALSE.

        integer(i4), allocatable :: tmpBuf(:),rtmpBuf(:)
        integer(i4), allocatable :: sCount(:),sCount2(:)
        integer(i4) :: idx,ig,ip,len
        logical :: found

        integer(i4) :: one
        integer(i4)     :: errorcode,errorlen
        character*(80)  :: errorstring

        ! ==============================
        ! Get some MPI information 
        ! ==============================
        call MPI_COMM_SIZE(COMM,nprocs,ierr)
        call MPI_COMM_RANK(COMM,my_task,ierr)

        lenHalo=nTotal-nActive
        maxNeigh=MIN(lenHalo,nprocs)

        lenSendBuffer=0
        lenRecvBuffer=0

        !=================================
        ! Allocate an array of size nprocs
        !=================================
        allocate(sCount(0:nprocs-1),sCount2(0:nprocs-1))

        allocate(sNeigh(maxNeigh))
        allocate(rNeigh(maxNeigh))
        allocate(RecvCnt(maxNeigh))
        allocate(ptrRecv(maxNeigh))
        rNeigh=-1
        sNeigh=-1
        RecvCnt=0

        nSend=0
        nRecv=0

        ! ==========================
        ! Get a list of my neighbors
        ! ==========================
        one=1
        do i=nActive+1,nTotal
           call InsertIntoArray(rNeigh,RecvCnt,LinearProc(i),one,ierr)
        enddo

        ! ==========================
        ! Count neighbors
        ! ==========================
        nRecv=COUNT(rNeigh .ge. 0)

        sCount=0
        do i=1,nRecv
           sCount(rNeigh(i)) = 1
        enddo
        call MPI_Allreduce(sCount,sCount2,nprocs,MPI_INTEGER,MPI_SUM,COMM,ierr)

        nSend=sCount2(my_task)
        deallocate(sCount,sCount2)
        if(Debug) &
            print *,'Processor: ',my_task,' has R-neighbors: ',rNeigh(1:nRecv)
        rWords = SUM(RecvCnt(1:nRecv))
        if(Info) print *,'IAM: ',my_task,'nNeigh: ',nRecv, &
            ' Bytes received: ', 8*rWords,' : ', 8*RecvCnt(1:nRecv)

        allocate(SendCnt(maxNeigh))
        allocate(SendTMP(maxNeigh))

        SendCnt=0

        if(nRecv>0) then 
            ptrRecv(1)=1
            do i=2,nRecv
                ptrRecv(i) = ptrRecv(i-1)+RecvCnt(i-1)
            enddo
            lenRecvBuffer=ptrRecv(nRecv)+RecvCnt(nRecv)-1
        endif

        if(Debug) &
            print *,'Processor: ',my_task,' Receive pointer: ',ptrRecv(1:nRecv)
        
        ! ===============================================
        ! Pack all the info into the Schedule_t structure
        ! ===============================================
        schedHandle%nRecv        =  nRecv
        schedHandle%nSend        =  nSend

        allocate(schedHandle%ptrRecv(nRecv))
        schedHandle%ptrRecv(1:nRecv) = ptrRecv(1:nRecv)

        allocate(schedHandle%RecvCnt(nRecv))
        schedHandle%RecvCnt(1:nRecv)=RecvCnt(1:nRecv)

        allocate(schedHandle%rNeigh(nRecv))
        schedHandle%rNeigh(1:nRecv) = rNeigh(1:nRecv)

        schedHandle%lenRecvBuffer =  lenRecvBuffer

        allocate(schedHandle%Srequest(nSend),schedHandle%Rrequest(nRecv))
        allocate(schedHandle%Sstatus(MPI_STATUS_SIZE,nSend), &
            schedHandle%Rstatus(MPI_STATUS_SIZE, nRecv))

        !=====================================================
        ! Figure out the number of points to Send to Neighbor
        !=====================================================
        tag = NUMexpectMSG

        do i=1,nSend
            call MPI_Irecv(SendTMP(i), 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                tag, COMM, schedHandle%Srequest(i), ierr)
            if(ierr .ne. MPI_SUCCESS) then 
                errorcode=ierr
                call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
                print *,'GSHalo_init: After call to MPI_Irecv: ',errorstring
            endif
        enddo
        
        do i=1,nRecv
            call MPI_Isend(RecvCnt(i), 1, MPI_INTEGER, rNeigh(i), &
                tag, COMM, schedHandle%Rrequest(i), ierr)
            if(ierr .ne. MPI_SUCCESS) then 
                errorcode=ierr
                call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
                print *,'GSHalo_init: After call to MPI_Isend: ',errorstring
            endif
        enddo
        
        if(nRecv>0) call MPI_Waitall(nRecv, schedHandle%Rrequest, &
             schedHandle%Rstatus, ierr)
        if(nSend>0) call MPI_Waitall(nSend, schedHandle%Srequest, &
             schedHandle%Sstatus, ierr)

        call MPI_Barrier(comm,ierr)

        do i=1,nSend
            cnt = SendTmp(i)
            dest = schedHandle%Sstatus(MPI_SOURCE, i)
            if(Debug) &
                print *,'IAM: ',my_task,' GSHalo_init: Message source is: ', &
                    dest, ' Length: ',cnt
            call InsertIntoArray(sNeigh, SendCnt, dest, cnt, ierr)
        enddo

        if(Debug) &
            print *,'Processor: ',my_task,' has S-neighbors: ',sNeigh(1:nSend)
        if(nSend>0) then 
            sWords = SUM(SendCnt(1:nSend))
        else
            sWords = 0
        endif
        if(Info) &
            print *,'IAM: ',my_task,'nNeigh: ',nSend,' Bytes sent: ', &
                8*sWords,' : ',8*SendCnt(1:nSend)
        nNeigh = nSend

        if(Info) then 
            call MPI_Allreduce( &
                sWords,sWords_total, 1, MPI_INTEGER,MPI_SUM,COMM, ierr)
            call MPI_Allreduce( &
                sWords,sWords_min, 1, MPI_INTEGER,MPI_MIN,COMM, ierr)
            call MPI_Allreduce( &
                sWords,sWords_max, 1, MPI_INTEGER,MPI_MAX,COMM, ierr)
            call MPI_Allreduce( &
                nNeigh,nNeigh_avg, 1, MPI_INTEGER,MPI_SUM,COMM, ierr)
            call MPI_Allreduce( &
                nNeigh,nNeigh_min, 1, MPI_INTEGER,MPI_MIN,COMM, ierr)
            call MPI_Allreduce( &
                nNeigh,nNeigh_max, 1, MPI_INTEGER,MPI_MAX,COMM, ierr)

            if( my_task == 0) then 
                sWords_avg = sWords_total/nprocs
                print *,'gsHalo_init: Bytes {total,min,avg,max}: ', &
                    8*sWords_total,8*sWords_min, &
                    8*sWords_avg,  8*sWords_max 
                nNeigh_avg=nNeigh_avg/nprocs
                print *, 'gsHalo_init: Neighbors {min,avg,max}: ', &
                    nNeigh_min,nNeigh_avg,nNeigh_max
            endif
        endif

        ! ============================================
        ! allocate and fill out the Send pointer array 
        ! ============================================
        allocate(schedHandle%sNeigh(nSend))
        allocate(ptrSend(nSend))
        if(nSend>0) then
            schedHandle%sNeigh(1:nSend) = sNeigh(1:nSend)
            ptrSend(1)=1
            do i=2,nSend
                ptrSend(i) = ptrSend(i-1)+SendCnt(i-1)
            enddo
            lenSendBuffer=ptrSend(nSend)+SendCnt(nSend)-1
        endif

        allocate(ptrCnt(nRecv))
        allocate(tmpBuf(lenRecvBuffer))
        allocate(rtmpBuf(lenSendBuffer))

        ptrCnt=0
        do i=nActive+1,nTotal
            ! find index into Neigh array corresponding to processor
            ! LinearProc(i)
            call LinearOrderedFind(rNeigh(1:nRecv),LinearProc(i),found,ip)
            ! global dof
            ig = LinearGdof(i)
            ! calculate index into message buffer
            idx = ptrRecv(ip)+ptrCnt(ip)
            ! Store the dof in the right spot
            tmpBuf(idx)=ig
            ! Increment the the message buffer offset to processor
            ! LinearProc(i)
            ptrCnt(ip)=ptrCnt(ip)+1
        enddo

        nPrint=lenRecvBuffer
        allocate(schedHandle%halo2send(lenSendBuffer))
        allocate(schedHandle%recv2halo(lenRecvBuffer))

        do i=1,lenRecvBuffer
            ig = tmpBuf(i)
            call LinearFind(LinearGdof(nActive+1:nTotal),ig,found,idx)
            schedHandle%recv2halo(i) = nActive+idx
        enddo

        do i=1,nSend
            len = SendCnt(i)
            src = sNeigh(i)
            if(Debug) print *,'IAM: ',my_task,' Posting Recv from: ',src, &
                ' Length: ',len
        enddo
        call MPI_Barrier(COMM,ierr)

        !============================================
        ! Figure out which points to send to Neighbor 
        !============================================
        tag = expectMSG
        do i=1,nSend
            len = SendCnt(i)
            iptr = ptrSend(i)
            src  = sNeigh(i)
            call MPI_Irecv(rtmpBuf(iptr), len, MPI_INTEGER, src, &
                tag, COMM, schedHandle%Srequest(i), ierr)
        enddo

        do i=1,nRecv
            len  = RecvCnt(i)
            iptr = ptrRecv(i)
            dest = rNeigh(i)
            call MPI_Isend(tmpBuf(iptr), len, MPI_INTEGER, dest, &
                tag, COMM, schedHandle%Rrequest(i),ierr)
        enddo

        if(nSend>0) call MPI_Waitall(nSend, schedHandle%Srequest, &
             schedHandle%Sstatus, ierr)
        if(nRecv>0) call MPI_Waitall(nRecv, schedHandle%Rrequest, &
             schedHandle%Rstatus, ierr)

        do i=1,lenSendBuffer
            ig = rtmpBuf(i)
            call LinearOrderedFind(LinearGdof(1:nActive), ig, found, idx)
            schedHandle%halo2send(i) = idx
        enddo

        ! ========================
        ! Assign the send arrays
        ! ========================
        allocate(schedHandle%SendCnt(nSend))
        schedHandle%SendCnt(1:nSend)  = SendCnt(1:nSend)

        allocate(schedHandle%ptrSend(nSend))
        schedHandle%ptrSend(1:nSend) = ptrSend(1:nSend)

        schedHandle%COMM          = COMM
        schedHandle%lenSendBuffer = lenSendBuffer

        ! ===================================
        ! allocate the message buffers
        if (present(int_msg)) then
           if (int_msg) then
              allocate(schedHandle%bufferRecvInt(lenRecvBuffer))
              allocate(schedHandle%bufferSendInt(lenSendBuffer))
           else
              allocate(schedHandle%bufferRecvDbl(lenRecvBuffer))
              allocate(schedHandle%bufferSendDbl(lenSendBuffer))
           end if
        else
           allocate(schedHandle%bufferRecvDbl(lenRecvBuffer))
           allocate(schedHandle%bufferSendDbl(lenSendBuffer))
           
           allocate(schedHandle%bufferRecvInt(lenRecvBuffer))
           allocate(schedHandle%bufferSendInt(lenSendBuffer))
        end if
        
        if(Debug) print *,'IAM: ',my_task, &
            ' GSHalo_init: finished with subroutine'
        call MPI_Barrier(COMM,ierr)

        !clean up
        deallocate(sNeigh)
        deallocate(rNeigh)
        deallocate(RecvCnt)
        deallocate(ptrRecv)
        deallocate(SendCnt)
        deallocate(SendTMP)
        deallocate(ptrCnt)
        deallocate(tmpBuf)
        deallocate(rtmpBuf)
        deallocate(ptrSend)

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
