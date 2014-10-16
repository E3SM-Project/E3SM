#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#undef DEBUGPART
module schedule_mod 
  use metagraph_mod, only : MetaEdge_t
  use kinds, only : int_kind, iulog
  use schedtype_mod, only : Cycle_t, Schedule_t, schedule 

  implicit none 
  private 

  type, public :: GraphStats_t
     integer(kind=int_kind) :: offnode
     integer(kind=int_kind) :: onnode
     integer(kind=int_kind) :: LB
     integer(kind=int_kind) :: padding 
  end type GraphStats_t

  integer,public,parameter :: BNDRY_EXCHANGE_MESSAGE=10
  integer,private,allocatable,target  :: Global2Local(:)

  integer :: MinNelemd,MaxNelemd

  public :: genEdgeSched              ! Setup the communication schedule for the edge based boundary exchange
  public :: PrintSchedule, PrintCycle
  public :: CheckSchedule
  public :: FindBufferSlot
!  public :: MessageStats
!  public :: PrimMessageStats

contains

  !********************** GENCOMSCHED.F ******************************************
  subroutine genEdgeSched(elem, PartNumber,LSchedule,MetaVertex)
    use element_mod, only : element_t
    use metagraph_mod, only : metavertex_t
    use dimensions_mod, only : nelem, max_neigh_edges
    use gridgraph_mod, only : gridvertex_t, gridedge_t, assignment ( = )
#ifdef _MPI
    use parallel_mod, only : nComPoints, iam, mpi_status_size, rrequest, srequest, &
	status, npackpoints
#else
    use parallel_mod, only : nComPoints, iam
#endif
     use parallel_mod, only : haltmp, iam
    implicit none
    type(element_t), intent(inout)        :: elem(:)
    integer, intent(in)                :: PartNumber
    type (schedule_t), intent(inout)   :: LSchedule
    type (MetaVertex_t),intent(inout)  :: MetaVertex

    integer                       :: lengthP,total_length,lengthp_ghost
    integer                       :: i,j,is,ir,ncycle
    integer                       :: il,ie,ig
    integer                       :: nelemd0
    integer                       :: jmd
    integer                       :: inbr
    integer                       :: nSched
    integer,allocatable           :: tmpP(:,:)
    integer,allocatable           :: tmpP_ghost(:,:)
    integer                       :: nSend,nRecv,nedges 
    integer                       :: icycle
    integer			  :: iSched
    logical, parameter            :: VerbosePrint=.FALSE.
    logical, parameter            :: Debug=.false.
    integer :: ierr
    integer :: l1,l2,l1id,l2id


    nSched=SIZE(schedule)
    ! ================================================
    ! allocate some arrays for the call to MPI_gatherv
    ! ================================================

    MinNelemd = nelem
    MaxNelemd = 0
    ! =====================================================
    ! It looks like this is only used in this routine...
    ! so no need to put it in the schedule data-structure
    ! =====================================================
    allocate(Global2Local(nelem))
    if(Debug) write(iulog,*)'genEdgeSched: point #1'
    iSched = PartNumber

    nelemd0 = MetaVertex%nmembers
    if(VerbosePrint) then
       if(iam .eq. 1)  write(iulog,*)'genEdgeSched: Part # ',i,' has ',nelemd0, ' elements '
    endif
    MaxNelemd = AMAX0(MaxNelemd,nelemd0)
    MinNelemd = AMIN0(MinNelemd,nelemd0)
    if(Debug) write(iulog,*)'genEdgeSched: point #2'

    if(Debug) write(iulog,*)'genEdgeSched: point #3'
    LSchedule%ncycles = MetaVertex%nedges
    LSchedule%nelemd  = nelemd0
    if(Debug) write(iulog,*)'genEdgeSched: point #4'

    !  Note the minus one is for the internal node
    nedges = MetaVertex%nedges
    if(2*(nedges/2) .eq. nedges) then
       nedges = nedges/2
    else
       nedges = (nedges-1)/2
    endif
    LSchedule%nSendCycles = nedges
    LSchedule%nRecvCycles = nedges
    if(Debug) write(iulog,*)'genEdgeSched: point #5'

    ! Temporary array to calculate the Buffer Slot
    allocate(tmpP(2,nedges+1))
    allocate(tmpP_ghost(2,nedges+1))

    tmpP(1,:) = -1
    tmpP(2,:) = 0
    tmpP_ghost(1,:) = -1
    tmpP_ghost(2,:) = 0

    !  Allocate all the cycle structures
    allocate(LSchedule%SendCycle(nedges))
    allocate(LSchedule%RecvCycle(nedges))
    allocate(LSchedule%MoveCycle(1))

    ! Initialize the schedules...
    LSchedule%MoveCycle(1)%ptrP = 0
    LSchedule%MoveCycle(1)%lengthP = 0
    if(Debug) write(iulog,*)'genEdgeSched: point #6'

    !==================================================================
    !  Allocate and initalized the index translation arrays
    Global2Local = -1
    allocate(LSchedule%Local2Global(nelemd0))
    if(Debug) write(iulog,*)'genEdgeSched: point #7'

    do il=1,nelemd0
       ig     = MetaVertex%members(il)%number
       Global2Local(ig)=il
       LSchedule%Local2Global(il)=ig
#ifndef _PREDICT
       elem(il)%desc%putmapP=-1
       elem(il)%desc%getmapP=-1
       elem(il)%desc%putmapP_ghost=-1
       elem(il)%desc%getmapP_ghost=-1
       elem(il)%desc%reverse = .FALSE.
#endif
    enddo
    !==================================================================
    if(Debug) write(iulog,*)'genEdgeSched: point #8'



    total_length = 0
    ncycle = LSchedule%ncycles
    is=1
    ir=1
    do j=1,ncycle
       lengthP     =  MetaVertex%edges(j)%wgtP
       lengthP_ghost     =  MetaVertex%edges(j)%wgtP_ghost

       if((MetaVertex%edges(j)%HeadVertex == PartNumber) .AND. &
            (MetaVertex%edges(j)%TailVertex == PartNumber)) then
          inbr                            = PartNumber
          if(Debug) write(iulog,*)'genEdgeSched: point #9', iam
          LSchedule%MoveCycle%ptrP         = FindBufferSlot(inbr,lengthP,tmpP)
          LSchedule%MoveCycle%ptrP_ghost   = FindBufferSlot(inbr,lengthP_ghost,tmpP_ghost)
          call SetCycle(elem, LSchedule,LSchedule%MoveCycle(1),MetaVertex%edges(j))
          if(Debug) write(iulog,*)'genEdgeSched: point #10',iam
       else if (MetaVertex%edges(j)%TailVertex == PartNumber) then
          inbr                            = MetaVertex%edges(j)%HeadVertex
          if(Debug) write(iulog,*)'genEdgeSched: point #11', iam
          LSchedule%SendCycle(is)%ptrP     = FindBufferSlot(inbr,lengthP,tmpP)
          LSchedule%SendCycle(is)%ptrP_ghost= FindBufferSlot(inbr,lengthP_ghost,tmpP_ghost)
          call SetCycle(elem, LSchedule,LSchedule%SendCycle(is),MetaVertex%edges(j))
          if(Debug) write(iulog,*)'genEdgeSched: point #12',iam
          is = is+1
       else if (MetaVertex%edges(j)%HeadVertex == PartNumber) then
          inbr                            = MetaVertex%edges(j)%TailVertex
          if(Debug) write(iulog,*)'genEdgeSched: point #13',iam
          LSchedule%RecvCycle(ir)%ptrP     = FindBufferSlot(inbr,lengthP,tmpP)
          LSchedule%RecvCycle(ir)%ptrP_ghost= FindBufferSlot(inbr,lengthP_ghost,tmpP_ghost)
          call SetCycle(elem, LSchedule,LSchedule%RecvCycle(ir),MetaVertex%edges(j))
          if(Debug) write(iulog,*)'genEdgeSched: point #14',iam
          ir = ir+1
       endif
    enddo

    deallocate(tmpP)
    deallocate(tmpP_ghost)



    do ie=1,nelemd0
       ! compute number of neighbers for each element
       elem(ie)%desc%actual_neigh_edges=0
       do i=1,max_neigh_edges
          if (elem(ie)%desc%globalID(i)>0) then
             elem(ie)%desc%actual_neigh_edges=elem(ie)%desc%actual_neigh_edges+1
          endif
       enddo

       ! normally, we loop over max_neigh_edges, checking if there is an edge
       ! let's create a mapping so that we can loop over actual_neigh_edges
       ! sort in REVERSE global id order (so the ones with globalID=0 are last)
       do l1 = 1,max_neigh_edges-1
          do l2=l1+1,max_neigh_edges
             l1id=elem(ie)%desc%loc2buf(l1)
             l2id=elem(ie)%desc%loc2buf(l2)
             if (elem(ie)%desc%globalID(l2id) > elem(ie)%desc%globalID(l1id)) then
                ! swap index:
                l1id=elem(ie)%desc%loc2buf(l2)
                elem(ie)%desc%loc2buf(l2)=elem(ie)%desc%loc2buf(l1)
                elem(ie)%desc%loc2buf(l1)=l1id
             endif
          enddo
       enddo




       elem(ie)%vertex     = MetaVertex%members(ie)
       ig                  = MetaVertex%members(ie)%number
       elem(ie)%GlobalId   = ig
       elem(ie)%LocalId    = ie  
#if 0
       call LLInsertEdge(eroot,ig,jmd)
       !DBG write(iulog,*)'After call to LLInsertEdge in schedule: ie,ig ',ie,ig,jmd
#endif
    enddo

    deallocate(Global2Local)
    !S-JMD call CheckSchedule()
#ifdef _MPI
    !================================================================
    !     Allocate a couple of structures for bndry_exchange
    !        done here to remove it from the critical path
    !================================================================

    nComPoints=0

    nSend = nedges
    nRecv = nedges
    allocate(Rrequest(nRecv))
    allocate(Srequest(nSend))
    allocate(status(MPI_STATUS_SIZE,nRecv))

    !===============================================================
    !   Number of communication points ... to be used later to
    !    setup the size of the communication buffer for MPI_Ibsend
    !===============================================================
    do icycle=1,nSend
       nComPoints = nComPoints + LSchedule%SendCycle(icycle)%lengthP
    enddo
    nPackPoints = nComPoints + LSchedule%MoveCycle(1)%lengthP


    !   nbuf = 4*(nv+1)*nelemd*8*4*nlev
    !   write(iulog,*)'before call to allocate(combuffer) ',nbuf
    !   allocate(combuffer(nbuf))
    !   write(iulog,*)'IAM: ',iSched,'Before call to MPI_Buffer_Attach '
    !   call MPI_Buffer_Attach(combuffer,nbuf,ierr)
    !   write(iulog,*)'IAM: ',iSched,'After call to MPI_Buffer_Attach '

#endif
#ifdef DEBUGPART
    call haltmp("genEdgeSched: Just testing the partitioning algorithms")
#endif
  end subroutine genEdgeSched
#ifdef DOTHIS
  subroutine MessageStats(nlyr)
    use kinds, only : real_kind, log_kind
    use dimensions_mod, only : nmpi_per_node, nnodes, npart, nlev, np, nelem
    use perfmodel_mod, only : perf_t, commtime_t , setptopnetworkparams, & ! _EXTERNAL
         setserialparamsexplicit, setserialparamsimplicit
    use control_mod, only : integration
    !-----------------
    implicit none
    integer,intent(in)  ::   nlyr  ! number of 2D layers in the communication

    integer                       :: icycle,ip
    real(kind=real_kind)          :: lb_nelemd,lb_volume
    integer                       :: length,nSend,nRecv
    integer                       :: i

    real(kind=real_kind), allocatable :: Time_total(:),Time_calc(:), Time_comm(:)
    real(kind=real_kind), allocatable :: Time_calc1(:),Time_calc2(:)
    type (commtime_t), allocatable :: bndry1(:),bndry2(:),bndry3(:),bndryf(:)

    real(kind=real_kind) :: Time_comm1,Time_commf, &
         Time_comm2,Time_comm3
    integer(kind=int_kind) ::  bytes_per_point
    real(kind=real_kind) :: time_per_elem
    real(kind=real_kind) :: time_per_iter
    real(kind=real_kind) :: latency_tmp,bandwidth_tmp
    real(kind=real_kind) :: avg_cg_iters

    type (commtime_t)    :: offnode,onnode
    type (perf_t)        :: tcv
    type (GraphStats_t), allocatable  :: count(:)

    real(kind=real_kind) :: Time_serial,Time_parallel,Speedup
    integer(kind=int_kind),allocatable :: offnode_count(:),onnode_count(:)
    integer(kind=int_kind),allocatable :: LB_count(:)
    integer(kind=int_kind) :: edgecut_offnode,edgecut_onnode,edgecut_total
    logical(kind=log_kind) :: FoundNetwork,FoundMachine
    character(len=80) :: networkname,machinename

    integer(kind=int_kind) :: indx_min,indx_comm(1),indx_calc(1)
    integer(kind=8) :: imin_volume, &
         imax_volume
    integer(kind=int_kind),parameter  :: configuration = 5

    integer :: node1,node2,nbr,in
    integer :: nelemd0
    logical, parameter :: Debug = .FALSE.
    logical, parameter :: PredictPerformance = .TRUE.

    if(PredictPerformance) then 
       networkname="protoBGL"
       call SetPtoPNetworkParams(offnode,onnode,networkname,FoundNetwork)
       write(iulog,*)'MessageStats: After SetSerialParamsExplicit FoundNetwork: ',FoundNetwork

       machinename="protoBGL"
       !JMD integration = "explicit"
       write(iulog,*)'MessageStats: integration: ',integration
       if(integration == "explicit") then 
          call SetSerialParamsExplicit(time_per_elem,np,nlev,machinename,FoundMachine)
	  write(iulog,*)'MessageStats: After SetSerialParamsExplicit FoundMachine: ',FoundMachine
       else if(integration == "semi_imp") then 
          call SetSerialParamsImplicit(time_per_elem,time_per_iter,avg_cg_iters, &
               np,nelem,nlev,machinename,FoundMachine)
       else if(integration == "full_imp") then 
	  write(iulog,*)'MessageStats: not set for implicit integration'
       endif
    endif


    allocate(Time_total(npart))
    allocate(Time_calc(npart))
    allocate(Time_comm(npart))
    allocate(Time_calc1(npart))
    allocate(Time_calc2(npart))

    allocate(bndry1(npart))
    allocate(bndry2(npart))
    allocate(bndry3(npart))
    allocate(bndryf(npart))

    write(iulog,*)'MessageStats:  nlev and npart nnodes nmpi_per_node are: ',nlev,npart,nnodes,nmpi_per_node
    if(npart .gt. 1) then
       allocate(count(nnodes))
       call foo(schedule,count)

       !----------------------------------------------------
       ! Call fooCalc for each boundary exchange 
       !----------------------------------------------------
       write(iulog,*)'MessageStats: '
       write(iulog,*)'MessageStats: 		Boundary Exchange #1		'
       bytes_per_point = (3*nlev)*8
       call fooCalc(schedule,bytes_per_point,offnode,onnode,count,bndry1)

       write(iulog,*)'MessageStats: '
       write(iulog,*)'MessageStats: 		Filter Boundary Exchange		'
       bytes_per_point = (2*nlev)*8
       call fooCalc(schedule,bytes_per_point,offnode,onnode,count,bndryf)

       if(integration == "semi_imp") then 
          write(iulog,*)'MessageStats: '
          write(iulog,*)'MessageStats: 		Boundary Exchange #2		'
          bytes_per_point = (2*nlev)*8
          call fooCalc(schedule,bytes_per_point,offnode,onnode,count,bndry2)
          write(iulog,*)'MessageStats: '
          write(iulog,*)'MessageStats: 		Helmholtz Boundary Exchange 	'
          bytes_per_point = (2*nlev)*8
          call fooCalc(schedule,bytes_per_point,offnode,onnode,count,bndry3)
       endif

#ifdef _PREDICT
       if(PredictPerformance) then 
	  !------------------------------------------------------	
	  ! Calculate information about execution time prediction
	  !------------------------------------------------------	
          if(integration == "explicit") then 
             do ip=1,npart
                nelemd0 = Schedule(ip)%nelemd
                Time_calc(ip) = time_per_elem*nelemd0
                if(filter_freq > 0) then 
                   Time_commf = (bndryf(ip)%latency + bndryf(ip)%bandwidth)/real(filter_freq,kind=real_kind)
                else
                   Time_commf = 0
                endif
                Time_comm(ip)   = Time_commf + bndry1(ip)%latency + bndry1(ip)%bandwidth
             enddo
	  else if(integration == "semi_imp") then 
             do ip=1,npart
                nelemd0 = Schedule(ip)%nelemd
                Time_calc1(ip) = time_per_elem*nelemd0
                Time_calc2(ip) = time_per_iter*nelemd0*avg_cg_iters
                Time_calc(ip)  = Time_calc1(ip) + Time_calc2(ip)
                if(filter_freq > 0) then 
                   Time_commf = (bndryf(ip)%latency + bndryf(ip)%bandwidth)/real(filter_freq,kind=real_kind)
                else
                   Time_commf = 0
                endif
                Time_comm(ip)  = Time_commf &
                     + bndry1(ip)%latency + bndry1(ip)%bandwidth   &
                     +	bndry2(ip)%latency + bndry2(ip)%bandwidth   &
                     + avg_cg_iters*(bndry3(ip)%latency + bndry3(ip)%bandwidth)
             enddo
	  endif

	  !-----------------------------------------------------------
	  ! Print out information about the execution time prediction	
	  !-----------------------------------------------------------
          indx_calc = maxloc(Time_Calc)
          indx_comm = maxloc(Time_Comm)
          Time_parallel = Time_Calc(indx_calc(1)) + Time_comm(indx_comm(1))



          write(iulog,*)'MessageStats: '
          write(iulog,*)'MessageStats: 		Execution Time Prediction		'
          if(FoundMachine) write(iulog,*)'MessageStats: For Machine: ',TRIM(machinename)
          if(FoundNetwork) write(iulog,*)'MessageStats: With Network: ',TRIM(networkname)

          if(integration == "explicit" ) then 
             Time_serial = real(nelem,kind=real_kind)*time_per_elem
             Speedup=Time_serial/Time_parallel
	     write(iulog,201) Time_serial,1.0D-6*Time_serial*(secpday/tstep)
	     !---------------------------------------
	     ! First Boundary exchange in advance
	     !---------------------------------------
             if(FoundNetwork) write(iulog,65) bndry1(indx_comm(1))%latency + bndry1(indx_comm(1))%bandwidth, &
                  bndry1(indx_comm(1))%latency,bndry1(indx_comm(1))%bandwidth
	     !---------------------------------------
	     ! Filter Boundary exchange in advance
	     !---------------------------------------
             if(filter_freq > 0) then 
                Time_commf   = (bndryf(indx_comm(1))%latency  &
                     + bndryf(indx_comm(1))%bandwidth)/real(filter_freq,kind=real_kind)
                if(FoundNetwork) write(iulog,68) Time_commf, bndryf(indx_comm(1))%latency,bndryf(indx_comm(1))%bandwidth
             endif
	     !---------------------------------------
             !  Total Communication Cost
	     !---------------------------------------
             if(FoundMachine .and. FoundNetwork) &
                  write(iulog,55) Time_Parallel,Time_calc(indx_calc(1)), Time_comm(indx_comm(1))
          else if(integration ==  "semi_imp") then 
             Time_serial = real(nelem,kind=real_kind)*(time_per_elem + avg_cg_iters*time_per_iter)
             Speedup     = Time_serial/Time_parallel
	     write(iulog,202) avg_cg_iters
	     write(iulog,201) Time_serial,1.0D-6*Time_serial*(secpday/tstep)
             write(iulog,35) MAXVAL(Time_calc1)
             write(iulog,45) MAXVAL(Time_calc2)
	     !---------------------------------------
	     ! Filter Boundary exchange in advance_si
	     !---------------------------------------
             if(filter_freq > 0) then 
                Time_commf   = (bndryf(indx_comm(1))%latency  &
                     + bndryf(indx_comm(1))%bandwidth)/real(filter_freq,kind=real_kind)
                if(FoundNetwork) write(iulog,68) Time_commf, bndryf(indx_comm(1))%latency,bndryf(indx_comm(1))%bandwidth
             endif
	     !---------------------------------------
	     ! First Boundary exchange in advance_si
	     !---------------------------------------
             Time_comm1   = bndry1(indx_comm(1))%latency + bndry1(indx_comm(1))%bandwidth
             if(FoundNetwork) write(iulog,65) Time_comm1, bndry1(indx_comm(1))%latency,bndry1(indx_comm(1))%bandwidth

	     !---------------------------------------
	     ! Second Boundary exchange in advance_si
	     !---------------------------------------
             Time_comm2 = bndry2(indx_comm(1))%latency + bndry2(indx_comm(1))%bandwidth
             if(FoundNetwork) write(iulog,66) Time_comm2, bndry2(indx_comm(1))%latency,bndry2(indx_comm(1))%bandwidth

	     !---------------------------------------
	     ! Boundary exchange in solver_mod
	     !---------------------------------------
             Time_comm3 = bndry3(indx_comm(1))%latency + bndry3(indx_comm(1))%bandwidth
             if(FoundNetwork)  then
		write(iulog,67) avg_cg_iters*Time_comm3, avg_cg_iters*bndry3(indx_comm(1))%latency, &
                     avg_cg_iters*bndry3(indx_comm(1))%bandwidth
             endif

	     !---------------------------------------
             !  Total Communication Cost
	     !---------------------------------------
             if(FoundMachine .and. FoundNetwork) &
                  write(iulog,55) Time_Parallel,Time_calc(indx_calc(1)), &
                  Time_comm1 + Time_commf + Time_comm2 + avg_cg_iters*Time_comm3
          endif

          if(FoundMachine .and. FoundNetwork) write(iulog,95) Speedup
       endif
       call haltmp("MessageStats: Just testing the partitioning algorithms")
#endif
    endif
110 format(1x,A,I10,I10,f12.2,f12.2)
100 format(1x,A,I8,I8,f10.2,f8.5)
35  format(" MessageStats: Predicted SI STEP       Time (usec/step) ",f12.1)
45  format(" MessageStats: Predicted PCG SOLVER    Time (usec/step) ",f12.1)
55  format(" MessageStats: Predicted TOTAL      Time (usec/step) total:= ",f12.1," calc:= ",f12.1," comm:= ",f12.1)
65  format(" MessageStats: Predicted BNDRY #1   Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
165 format(" MessageStats: Predicted BNDRY #1 ",i5,"  Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
66  format(" MessageStats: Predicted BNDRY #2   Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
67  format(" MessageStats: Predicted HELM BNDRY Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
68  format(" MessageStats: Predicted FLTR BNDRY Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
75  format(" MessageStats: Edgecut   total: ",i10," onnode: ",i10," offnode: ",i10)
95  format(" MessageStats: Speedup:  ",f10.1)
201 format(" MessageStats: Serial Execution Time (usec/step)",f12.1," (sec/day) ",f12.1)
202 format(" MessageStats: Average CG iterations: ",f6.3)

120 format(A,I8,I8,f8.5)

  end subroutine MessageStats

  subroutine PrimMessageStats(nlyr)
    use kinds, only : real_kind, log_kind
    use perfmodel_mod, only : perf_t, commtime_t , setptopnetworkparams, & ! _EXTERNAL
         setserialparamsprimexplicit
    use dimensions_mod, only : nmpi_per_node, nnodes, npart, nlev, np, nelem, nelemd
    use control_mod, only : integration
    !-----------------
    implicit none
    integer,intent(in)  ::   nlyr  ! number of 2D layers in the communication

    integer                       :: icycle,ip
    real(kind=real_kind)          :: lb_nelemd,lb_volume
    integer                       :: length,nSend,nRecv
    integer                       :: i

    real(kind=real_kind), allocatable :: Time_total(:),Time_calc(:), Time_comm(:)
    real(kind=real_kind), allocatable :: Time_calc1(:),Time_calc2(:)
    type (commtime_t), allocatable :: bndry1(:),bndry2(:),bndry3(:),bndryf(:)

    real(kind=real_kind) :: Time_comm1,Time_commf, &
         Time_comm2,Time_comm3
    integer(kind=int_kind) ::  bytes_per_point
    real(kind=real_kind) :: time_per_elem
    real(kind=real_kind) :: time_per_iter
    real(kind=real_kind) :: latency_tmp,bandwidth_tmp
    real(kind=real_kind) :: avg_cg_iters

    type (commtime_t)    :: offnode,onnode
    type (perf_t)        :: tcv
    type (GraphStats_t), allocatable  :: count(:)

    real(kind=real_kind) :: Time_serial,Time_parallel,Speedup
    integer(kind=int_kind),allocatable :: offnode_count(:),onnode_count(:)
    integer(kind=int_kind),allocatable :: LB_count(:)
    integer(kind=int_kind) :: edgecut_offnode,edgecut_onnode,edgecut_total
    logical(kind=log_kind) :: FoundNetwork,FoundMachine
    character(len=80) :: networkname,machinename

    integer(kind=int_kind) :: indx_min,indx_comm(1),indx_calc(1)
    integer(kind=8) :: imin_volume, &
         imax_volume
    integer(kind=int_kind),parameter  :: configuration = 5

    integer :: node1,node2,nbr,in
    integer :: nelemd0
    logical, parameter :: Debug = .FALSE.
    logical, parameter :: PredictPerformance = .TRUE.

    if(PredictPerformance) then 
       networkname="protoBGL"
       call SetPtoPNetworkParams(offnode,onnode,networkname,FoundNetwork)
       write(iulog,*)'PrimMessageStats: After SetSerialParamsExplicit FoundNetwork: ',FoundNetwork

       machinename="protoBGL"
       integration = "explicit"
       write(iulog,*)'PrimMessageStats: integration: ',integration
       if(integration == "explicit") then 
          call SetSerialParamsPrimExplicit(time_per_elem,np,nlev,machinename,FoundMachine)
	  write(iulog,*)'PrimMessageStats: After SetSerialParamsPrimExplicit FoundMachine: ',FoundMachine
       endif
    endif


    allocate(Time_total(npart))
    allocate(Time_calc(npart))
    allocate(Time_comm(npart))
    allocate(Time_calc1(npart))
    allocate(Time_calc2(npart))

    allocate(bndry1(npart))
    allocate(bndry2(npart))
    allocate(bndry3(npart))
    allocate(bndryf(npart))

    write(iulog,*)'PrimMessageStats:  nlev and npart nnodes nmpi_per_node are: ',nlev,npart,nnodes,nmpi_per_node
    if(npart .gt. 1) then
       allocate(count(nnodes))
       call foo(schedule,count)

       !----------------------------------------------------
       ! Call fooCalc for each boundary exchange 
       !----------------------------------------------------
       write(iulog,*)'PrimMessageStats: '
       write(iulog,*)'PrimMessageStats: 		Boundary Exchange #1		'
       bytes_per_point = 3*8
       call fooCalc(schedule,bytes_per_point,offnode,onnode,count,bndry1)

       write(iulog,*)'PrimMessageStats: '
       write(iulog,*)'PrimMessageStats: 		Boundary Exchange #2		'
       bytes_per_point = (3*nlev+1)*8
       call fooCalc(schedule,bytes_per_point,offnode,onnode,count,bndry2)

       write(iulog,*)'PrimMessageStats: '
       write(iulog,*)'PrimMessageStats: 		Filter Boundary Exchange		'
       bytes_per_point = (3*nlev+1)*8
       call fooCalc(schedule,bytes_per_point,offnode,onnode,count,bndryf)

#ifdef _PREDICT
       if(PredictPerformance) then 
	  !------------------------------------------------------	
	  ! Calculate information about execution time prediction
	  !------------------------------------------------------	
          if(integration == "explicit") then 
             do ip=1,npart
                nelemd0 = Schedule(ip)%nelemd
                Time_calc(ip) = time_per_elem*nelemd0
                if(filter_freq > 0) then 
                   Time_commf = (bndryf(ip)%latency + bndryf(ip)%bandwidth)/real(filter_freq,kind=real_kind)
                else
                   Time_commf = 0
                endif
                Time_comm(ip)   = Time_commf +     bndry1(ip)%latency + bndry1(ip)%bandwidth &
                     + 4.*(bndry2(ip)%latency + bndry2(ip)%bandwidth)
             enddo
	  endif

	  !-----------------------------------------------------------
	  ! Print out information about the execution time prediction	
	  !-----------------------------------------------------------
          indx_calc = maxloc(Time_Calc)
          indx_comm = maxloc(Time_Comm)
          Time_parallel = Time_Calc(indx_calc(1)) + Time_comm(indx_comm(1))



          write(iulog,*)'PrimMessageStats: '
          write(iulog,*)'PrimMessageStats: 		Execution Time Prediction		'
          if(FoundMachine) write(iulog,*)'PrimMessageStats: For Machine: ',TRIM(machinename)
          if(FoundNetwork) write(iulog,*)'PrimMessageStats: With Network: ',TRIM(networkname)

          if(integration == "explicit" ) then 
             Time_serial = real(nelem,kind=real_kind)*time_per_elem
             Speedup=Time_serial/Time_parallel
	     write(iulog,201) Time_serial,1.0D-6*Time_serial*(secpday/tstep)
	     !---------------------------------------
	     ! Small Boundary exchange in prim_advance
	     !---------------------------------------
             if(FoundNetwork) write(iulog,65) bndry1(indx_comm(1))%latency + bndry1(indx_comm(1))%bandwidth, &
                  bndry1(indx_comm(1))%latency,bndry1(indx_comm(1))%bandwidth

	     !-------------------------------------------
	     ! 4 Large Boundary exchanges in prim_advance
	     !-------------------------------------------
             if(FoundNetwork) write(iulog,65) 4.*(bndry2(indx_comm(1))%latency + bndry2(indx_comm(1))%bandwidth), &
                  4.*bndry2(indx_comm(1))%latency,4.*bndry1(indx_comm(1))%bandwidth
	     !---------------------------------------
	     ! Filter Boundary exchange in advance
	     !---------------------------------------
             if(filter_freq > 0) then 
                Time_commf   = (bndryf(indx_comm(1))%latency  &
                     + bndryf(indx_comm(1))%bandwidth)/real(filter_freq,kind=real_kind)
                if(FoundNetwork) write(iulog,68) Time_commf, bndryf(indx_comm(1))%latency,bndryf(indx_comm(1))%bandwidth
             endif
	     !---------------------------------------
             !  Total Communication Cost
	     !---------------------------------------
             if(FoundMachine .and. FoundNetwork) &
                  write(iulog,55) Time_Parallel,Time_calc(indx_calc(1)), Time_comm(indx_comm(1))
          else if(integration ==  "semi_imp") then 
             write(iulog,*)'PrimMessageStats: not yet support for Semi-Implicit time integration'
          endif

          if(FoundMachine .and. FoundNetwork) write(iulog,95) Speedup
       endif
       call haltmp("PrimMessageStats: Just testing the partitioning algorithms")
#endif
    endif
110 format(1x,A,I10,I10,f12.2,f12.2)
100 format(1x,A,I8,I8,f10.2,f8.5)
35  format(" PrimMessageStats: Predicted SI STEP       Time (usec/step) ",f12.1)
45  format(" PrimMessageStats: Predicted PCG SOLVER    Time (usec/step) ",f12.1)
55  format(" PrimMessageStats: Predicted TOTAL      Time (usec/step) total:= ",f12.1," calc:= ",f12.1," comm:= ",f12.1)
65  format(" PrimMessageStats: Predicted BNDRY #1   Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
165 format(" PrimMessageStats: Predicted BNDRY #1 ",i5,"  Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
66  format(" PrimMessageStats: Predicted BNDRY #2   Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
67  format(" PrimMessageStats: Predicted HELM BNDRY Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
68  format(" PrimMessageStats: Predicted FLTR BNDRY Time (usec/step) Total:=",f12.1," latency:= ",f12.1," band:= ",f12.1)
75  format(" PrimMessageStats: Edgecut   total: ",i10," onnode: ",i10," offnode: ",i10)
95  format(" PrimMessageStats: Speedup:  ",f10.1)
201 format(" PrimMessageStats: Serial Execution Time (usec/step)",f12.1," (sec/day) ",f12.1)
202 format(" PrimMessageStats: Average CG iterations: ",f6.3)

120 format(A,I8,I8,f8.5)

  end subroutine PrimMessageStats
#endif

  subroutine CheckSchedule()
    implicit none 

    integer                             :: i,nSched,nbufferwords_1,nbufferwords_2
    type (Schedule_t), pointer          :: pSchedule

    nSched = SIZE(Schedule)

    do i=1,nSched
       pSchedule => Schedule(i)
       nbufferwords_1 = SUM(pSchedule%SendCycle(:)%lengthP)
       nbufferwords_2 = SUM(pSchedule%RecvCycle(:)%lengthP)
       if(nbufferwords_1 .ne. nbufferwords_2) then 
          write (*,100) i,nbufferwords_1, nbufferwords_2
       endif
    enddo
100 format('CheckSchedule: ERR IAM:',I3,' SIZEOF(SendBuffer):',I10,' != SIZEOF(RecvBuffer) :',I10)
110 format('CheckSchedule: ERR IAM:',I3,' LENGTH(SendBuffer):',I10,' != LENGTH(RecvBuffer) :',I10)

  end subroutine CheckSchedule
  subroutine PrintSchedule(Schedule)
    use gridgraph_mod, only : printgridedge

    implicit none
    type (Schedule_t),intent(in),target   :: Schedule(:)
    type (Schedule_t), pointer            :: pSchedule
    type (Cycle_t),pointer                :: pCycle

    integer               :: i,j,nSched

    nSched = SIZE(Schedule)

    write(*,*) '------NEW SCHEDULE FORMAT---------------------'
    do i=1,nSched
       pSchedule => Schedule(i)
       write(*,*)
       write(*,*) '----------------------------------------------'
       write(*,90) i,pSchedule%ncycles
       write(*,*) '----------------------------------------------'
       write(*,*) '-----------SEND-------------------------------'
       do j=1,pSchedule%nSendCycles
          pCycle => pSchedule%SendCycle(j)
          call PrintCycle(pCycle)
          call PrintGridEdge(pCycle%edge%members)
       enddo
       write(*,*) '-----------RECV-------------------------------'
       do j=1,pSchedule%nRecvCycles
          pCycle => pSchedule%RecvCycle(j)
          call PrintCycle(pCycle)
          call PrintGridEdge(pCycle%edge%members)
       enddo
       write(*,*) '-----------MOVE-------------------------------'
       pCycle => pSchedule%MoveCycle(1)
       call PrintCycle(pCycle)
       call PrintGridEdge(pCycle%edge%members)
    enddo
90  format('NODE # ',I2,2x,'NCYCLES ',I2)
97  format(10x,'EDGE #',I2,2x,'TYPE ',I1,2x,'G.EDGES',I4,2x,'WORDS ',I5,2x, &
         'SRC ',I3,2x,'DEST ',I3,2x,'PTR ',I4)
100 format(15x,I4,5x,I3,1x,'(',I1,') --',I1,'--> ',I3,1x,'(',I1,')')


  end subroutine PrintSchedule

#ifdef DOTHIS
  subroutine fooCalc(Schedule,bytes_per_point,offnode,onnode,count,Time)
    use kinds, only : real_kind
    use perfmodel_mod, only : commtime_t, perf_t ! _EXTERNAL
    use flops_mod, only : min_number,max_number, avg_number, & ! _EXTERNAL
         min_volume,max_volume, avg_volume, &
         min_message,max_message, avg_message, &
         LocalComVolume,TotalComVolume, tnsend
    use parallel_mod, only : iam, ncompoints
    use dimensions_mod, only : nmpi_per_node, npart
    ! results:: TnSend,{min,max}_{volume,number},Time_{total,calc,latency,bandwidth}
    ! inputs:: results + Schedule + bytes_per_point 

    type (Schedule_t),intent(in) :: Schedule(:)
    integer(kind=int_kind),intent(in) :: bytes_per_point
    type (commtime_t),intent(in)    :: offnode,onnode
    type (GraphStats_t), intent(in)   :: count(:)
    type (commtime_t),intent(inout)   :: Time(:)

    real(kind=real_kind)   :: bandwidth_tmp,latency_tmp
    integer(kind=int_kind) :: nbr,node1,node2,nsend,ip,icycle,nrecv
    integer(kind=int_kind) :: length
    real(kind=real_kind) :: lb_volume
    integer(kind=int_kind) :: imax_volume,imin_volume
    type (perf_t)          :: tcv
    integer                :: nSched
    integer                :: nelemd0
    logical, parameter     :: Debug = .FALSE.

    nSched=SIZE(Schedule)
    min_message = 1000000
    min_volume  = 1000000
    max_message = 0
    max_volume  = 0
    TotalComVolume = 0
    TnSend = 0

    tcv%onnode  = 0.0
    tcv%offnode = 0.0
    do ip=1,nSched
       nSend = Schedule(ip)%nSendCycles
       nRecv = Schedule(ip)%nRecvCycles
       nelemd0 = Schedule(ip)%nelemd
       nComPoints = 0
       Time(ip)%latency = 0
       Time(ip)%bandwidth = 0
       do icycle=1,nSend
          length = Schedule(ip)%SendCycle(icycle)%lengthV
          nComPoints = nComPoints + length
          min_message = MIN(min_message, bytes_per_point*length)
          max_message = MAX(max_message, bytes_per_point*length)
       enddo
       do icycle=1,nRecv
          length = Schedule(ip)%RecvCycle(icycle)%lengthV
          nbr = Schedule(ip)%RecvCycle(icycle)%source
          node1 = ((ip-1)/nmpi_per_node) + 1
          node2 = ((nbr-1)/nmpi_per_node) + 1
          if(node1 .eq. node2 )  then
             ! Message is set to an on-node processor
             tcv%onnode = tcv%onnode + length*bytes_per_point
             bandwidth_tmp = onnode%bandwidth
             latency_tmp = onnode%latency
          else
             ! Message is set to an off-node processor
             tcv%offnode = tcv%offnode + length*bytes_per_point
             bandwidth_tmp = real(count(node1)%offnode,kind=real_kind)*offnode%bandwidth
             latency_tmp = offnode%latency
             if(Debug) write(iulog,*)'fooCalc: node # ',node1,' offnode count is ',real(count(node1)%offnode,kind=real_kind)
             if(Debug) write(iulog,*)'fooCalc: node # ',node1,' offnode ={total,effective} ',offnode%bandwidth,bandwidth_tmp

          endif
          Time(ip)%latency = Time(ip)%latency + latency_tmp
          if(Debug) then 
             write(iulog,*)'fooCalc: node # ',node1,' Contribution in time ',length*bytes_per_point,  &
                  length*bytes_per_point*bandwidth_tmp
          endif
          Time(ip)%bandwidth = Time(ip)%bandwidth + length*bytes_per_point*bandwidth_tmp
       enddo
       !JMD  Note: I am removing the size of the communication buffer out of the computation
       LocalComVolume = bytes_per_point*nComPoints
       min_volume = MIN(min_volume,LocalComVolume)
       max_volume = MAX(max_volume,LocalComVolume)
       TnSend     = TnSend + nSend

       TotalComVolume = TotalComVolume + LocalComVolume
    enddo

    !------------------------------------------------------
    ! This Outputs information about each boundary exchange 
    !------------------------------------------------------
    avg_volume  = dble(TotalComVolume)/dble(npart)
    lb_volume  = (max_volume-avg_volume)/max_volume
    avg_message = dble(TotalComVolume)/dble(TnSend)
    imin_volume = NINT(min_volume)
    imax_volume = NINT(max_volume)
    if(iam .eq. 1) then 

       write(iulog,*)'MessageStats: Total Message volume is (Mbytes): ',1.0e-6*TotalComVolume
#if 1
       write(iulog,*) 'MessageStats: imin_volume',imin_volume
       write(iulog,*) 'MessageStats: imax_volume',imax_volume
       write(iulog,*) 'MessageStats: avg_volume',avg_volume
       write(iulog,*) 'MessageStats: lb_volume',lb_volume
#else
       !JMD For some reason on AIX this line breaks when the numbers are high
       write(iulog,110) 'MessageStats: Single process volume (bytes) {MIN,MAX,AVG,LB} =: ', &
            imin_volume, imax_volume, avg_volume,lb_volume
#endif
       write(iulog,*) 'MessageStats: Single message (bytes) {MIN,MAX,AVG} =: ',  &
            min_message,max_message,avg_message
       write(iulog,85) 1.0E-3*TotalComVolume,1.0E-3*tcv%onnode,1.0E-3*tcv%offnode
    endif
85  format(" MessageStats: TCV    total: ",f10.1," onnode: ",f10.1," offnode: ",f10.1)
110 format(1x,A,I10,I10,f12.2,f12.2)

  end subroutine fooCalc

  subroutine foo(Schedule,count)
    use kinds, only : real_kind, log_kind
    use flops_mod, only : min_number,max_number, avg_number ! _EXTERNAL
    use dimensions_mod, only : nmpi_per_node, npart, nelem
    use parallel_mod, only : iam
    type (Schedule_t),intent(in) :: Schedule(:)
    type (GraphStats_t),intent(inout) :: count(:)

    integer :: node2,nbr,ip,node1,icycle,nrecv,nsend
    logical :: offnode,found_offnode,found_onnode,found_lbnode

    integer(kind=int_kind) :: edgecut_total
    real(kind=real_kind)   :: lb_nelemd,avg_nelemd
    integer(kind=int_kind) :: edgecut_offnode,edgecut_onnode
    integer(kind=int_kind) :: TnSend
    integer :: nSched
    integer :: nelemd0

    logical(kind=log_kind),parameter :: Debug=.FALSE.

    nSched=SIZE(Schedule)
    count(:)%offnode=0
    count(:)%onnode=0
    count(:)%LB=0
    edgecut_offnode=0
    edgecut_onnode=0
    min_number  = 1000000
    max_number  = 0
    TnSend=0
    do ip=1,nSched
       node1=((ip-1)/nmpi_per_node)+1
       found_offnode=.FALSE.
       found_onnode=.FALSE.
       found_lbnode=.FALSE.
       nRecv = Schedule(ip)%nRecvCycles
       nSend = Schedule(ip)%nSendCycles
       nelemd0 = Schedule(ip)%nelemd
       do icycle=1,nRecv
          nbr = Schedule(ip)%RecvCycle(icycle)%source
          node2=((nbr-1)/nmpi_per_node)+1
          if(node2 .ne. node1 ) then
             ! Message is set off-node
             if(nelemd0 .eq. MaxNelemd) found_lbnode = .TRUE.
             found_offnode=.TRUE.
             edgecut_offnode = edgecut_offnode + 1
          else
             ! Message is set to an on-node processor 
             edgecut_onnode = edgecut_onnode + 1
          endif
          if(node2 .eq. node1 ) found_onnode = .TRUE.
       enddo
       if(found_offnode) count(node1)%offnode = count(node1)%offnode + 1
       if(found_onnode) count(node1)%onnode = count(node1)%onnode + 1
       if(found_lbnode) count(node1)%LB = count(node1)%LB + 1
       min_number = MIN(min_number,nSend)
       max_number = MAX(max_number,nSend)
       TnSend=TnSend+nSend
    enddo
    if(Debug) write(iulog,*)'count(:)%offnode : ',count(:)%offnode
    if(Debug) write(iulog,*)'count(:)%LB : ',count(:)%LB

    !-------------------------------------------------
    ! This Outputs general information about the Graph 
    !-------------------------------------------------
    edgecut_offnode = edgecut_offnode/2
    edgecut_onnode  = edgecut_onnode/2
    edgecut_total   = edgecut_offnode + edgecut_onnode
    avg_nelemd=dble(nelem)/dble(npart)
    lb_nelemd=(MaxNelemd-avg_nelemd)/MaxNelemd
    avg_number  = dble(TnSend)/dble(npart)
    if(iam .eq. 1) then
       write(iulog,*)'MessageStats: Number of MPI processes are:  ',npart
       write(iulog,100) 'MessageStats: nelemd is {MIN,MAX,AVG,LB}: ',MinNelemd,MaxNelemd,avg_nelemd,lb_nelemd
       write(iulog,*) 'MessageStats: number of neighbors  {MIN,MAX,AVG} ',min_number,max_number,avg_number
       write(iulog,75) edgecut_total,edgecut_onnode,edgecut_offnode
    endif

100 format(1x,A,I8,I8,f10.2,f8.5)
75  format(" MessageStats: Edgecut   total: ",i10," onnode: ",i10," offnode: ",i10)

  end subroutine foo
#endif
  subroutine PrintCycle(Cycle)

    implicit none 
    type (Cycle_t),intent(in),target  ::  Cycle

    write(*,97) Cycle%edge%number,Cycle%type,Cycle%edge%nmembers, &
         Cycle%lengthP,Cycle%source, Cycle%dest,Cycle%ptrP

97  format(5x,'METAEDGE #',I2,2x,'TYPE ',I1,2x,'G.EDGES',I4,2x,'WORDS ',I5,2x, &
         'SRC ',I3,2x,'DEST ',I3,2x,'PTR ',I5)

  end subroutine PrintCycle

  subroutine SetCycle(elem, schedule,Cycle,Edge)
    use element_mod, only : element_t
    use dimensions_mod, only : max_corner_elem, max_neigh_edges
    use parallel_mod, only : abortmp, iam   
    implicit none 

    type(element_t), intent(inout)        :: elem(:)
    type (Schedule_t),intent(inout)          :: Schedule
    type (Cycle_t),intent(inout)             :: Cycle
    type (MetaEdge_t),intent(in),target      :: Edge
    integer                                  :: i,il,face, loc, dir


#ifndef _PREDICT
    do i=1,Edge%nmembers
       !   Setup send index
       il                     = Global2Local(Edge%members(i)%tail%number)
       face                   = Edge%members(i)%tail_face
       !need to convert the location of corner elements for getmap and putmap
       if (face.ge.5) then ! if a corner element
          dir = Edge%members(i)%tail_dir
          loc = MOD(dir,max_corner_elem) !this is the location within that direction
          dir = (dir - loc)/max_corner_elem !this is the direction (1-8)
          loc = dir + (dir-5)*(max_corner_elem-1)+loc
       else
          loc = face
       end if

       if(il .gt. 0) then 
          elem(il)%desc%putmapP(loc) = Edge%edgeptrP(i) + Cycle%ptrP - 1  ! offset, so start at 0
          elem(il)%desc%putmapP_ghost(loc) = Edge%edgeptrP_ghost(i) + Cycle%ptrP_ghost  ! index, start at 1
          elem(il)%desc%reverse(loc) = Edge%members(i)%reverse
       endif



       !   Setup receive index
       il                     = Global2Local(Edge%members(i)%head%number)
       face                   = Edge%members(i)%head_face
       !need to convert the location of corner elements for getmap and putmap
       if (face.ge.5) then !its a corner
          dir = Edge%members(i)%head_dir
          loc = MOD(dir,max_corner_elem) !this is the location within that direction
          dir = (dir - loc)/max_corner_elem !this is the direction (1-8)
          loc = dir + (dir-5)*(max_corner_elem-1)+loc
          if(loc>max_neigh_edges) then
             print *,__FILE__,__LINE__,iam,face,i,max_corner_elem,max_neigh_edges,edge%members(i)%head_face
             call abortmp('max_neigh_edges set too low.')
          end if
       else
          loc = face
       end if

       if(il .gt. 0) then 
          elem(il)%desc%getmapP(loc) = Edge%edgeptrP(i) + Cycle%ptrP - 1
          elem(il)%desc%getmapP_ghost(loc) = Edge%edgeptrP_ghost(i) + Cycle%ptrP_ghost 
          elem(il)%desc%globalID(loc) = Edge%members(i)%tail%number
       endif


    enddo
#endif
    Cycle%edge   => Edge
    Cycle%type   = Edge%type
    Cycle%dest   = Edge%HeadVertex
    Cycle%source = Edge%TailVertex
    Cycle%tag    = BNDRY_EXCHANGE_MESSAGE
    Cycle%lengthP = Edge%wgtP
    Cycle%lengthP_ghost = Edge%wgtP_ghost

  end subroutine SetCycle

  function FindBufferSlot(inbr,length,tmp) result(ptr)

    integer                          :: ptr
    integer,intent(in)               :: inbr,length
    integer,intent(inout)    :: tmp(:,:)

    integer                          :: i,n

    n = SIZE(tmp,2)

    ptr = 0
    do i=1,n
       if( tmp(1,i) == inbr) then 
          ptr = tmp(2,i)
          return	
       endif
       if( tmp(1,i) == -1 ) then  
          tmp(1,i) = inbr
          if(i .eq. 1) tmp(2,i) = 1
          ptr = tmp(2,i)
          if(i .ne. n) tmp(2,i+1) = ptr +length
          return
       endif
    enddo

  end function FindBufferSlot

end module schedule_mod
