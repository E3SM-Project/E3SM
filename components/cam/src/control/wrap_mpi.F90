!----------------------------------------------------------------------- 
! 
! Purpose:
!
! 	Wrapper routines for the MPI (Message Passing) library for the
!	distributed memory (SPMD) version of the code. Also data with
!	"shorthand" names for the MPI data types.
!
! Entry points:
!      mpibarrier             Calls mpi_barrier
!      mpifinalize            Calls mpi_finalize
!      mpipack_size           Calls mpi_pack
!      mpipack                Calls mpi_pack
!      mpiunpack              Calls mpi_unpack
!      mpisendrecv            Calls mpi_sendrecv
!      mpiisend               Calls mpi_isend
!      mpiirsend              Calls mpi_irsend
!      mpiissend              Calls mpi_issend
!      mpiirecv               Calls mpi_irecv
!      mpiwait                Calls mpi_wait
!      mpiwaitall             Calls mpi_waitall
!      mpisend                Calls mpi_send
!      mpirsend               Calls mpi_rsend
!      mpissend               Calls mpi_ssend
!      mpirecv                Calls mpi_recv
!      mpigather              Calls mpi_gather
!      mpigatherv             Calls mpi_gatherv
!      mpigathervr4           Calls mpi_gatherv for real*4 data
!      mpigathervint          Calls mpi_gatherv for integer data
!      mpisum                 Calls mpi_sum
!      mpiscatter             Calls mpi_scatter
!      mpiscatterv            Calls mpi_scatterv
!      mpibcast               Calls mpi_bcast
!      mpiallmaxint           Calls mpi_allreduce on integer vector with mpi_max operator
!      mpialltoallv           Calls mpi_alltoallv
!      mpialltoallint         Calls mpi_alltoall for integer data
!      mpiallgatherv          Calls mpi_allgatherv
!      mpiallgatherint        Calls mpi_allgatherv for integer data
!      mpiwincreate           Calls mpi_win_create and mpi_win_fence
!
! Author: Many
! 
!-----------------------------------------------------------------------
!

!
! Disable the use of the MPI ready send protocol by default, to
! address recurrent issues with poor performance or incorrect 
! functionality in MPI libraries. When support is known to be robust,
! or for experimentation, can be re-enabled by defining the CPP token
! _USE_MPI_RSEND during the build process.
!
#ifndef _USE_MPI_RSEND
#define mpi_rsend mpi_send
#define mpi_irsend mpi_isend
#endif

!
! Compile these routines only when SPMD is defined
!
#if (defined SPMD)

!****************************************************************

   subroutine mpibarrier (comm)
!
! MPI barrier, have threads wait until all threads have reached this point
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog

   implicit none

   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call mpi_barrier (comm, ier)
   if (ier.ne.mpi_success) then
      write(iulog,*)'mpi_barrier failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpibarrier
 
!****************************************************************
 
   subroutine mpifinalize
!
! End of all MPI communication
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog

   implicit none

   integer ier   !MP error code
 
   call mpi_finalize (ier)
   if (ier.ne.mpi_success) then
      write(iulog,*)'mpi_finalize failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpifinalize
 
!****************************************************************
 
   subroutine mpipack_size (incount, datatype, comm, size)
!
! Returns the size of the packed data
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog

   implicit none

   integer, intent(in):: incount
   integer, intent(in):: datatype
   integer, intent(in):: comm
   integer, intent(out):: size
 
   integer ier   !MP error code
 
   call mpi_pack_size (incount, datatype, comm, size, ier)
   if (ier.ne.mpi_success) then
      write(iulog,*)'mpi_pack_size failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpipack_size
 
!****************************************************************
 
   subroutine mpipack (inbuf, incount, datatype, outbuf, outsize,    &
                       position, comm)
!
! Pack the data and send it.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog

   implicit none

   real(r8), intent(in):: inbuf(*)
   real(r8), intent(out):: outbuf(*)
   integer, intent(in):: incount
   integer, intent(in):: datatype
   integer, intent(out):: outsize
   integer, intent(inout):: position
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call mpi_pack (inbuf, incount, datatype, outbuf, outsize,         &
                  position, comm, ier)
   if (ier.ne.mpi_success) then
      write(iulog,*)'mpi_pack failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpipack
 
!****************************************************************
 
   subroutine mpiunpack (inbuf, insize, position, outbuf, outcount,  &
                         datatype, comm)
!
! Un-packs the data from the packed receive buffer
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog

   implicit none

   real(r8), intent(in):: inbuf(*)
   real(r8), intent(out):: outbuf(*)
   integer, intent(in):: insize
   integer, intent(inout):: position
   integer, intent(in):: outcount
   integer, intent(in):: datatype
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
   call mpi_unpack (inbuf, insize, position, outbuf, outcount,       &
                    datatype, comm, ier)
   if (ier.ne.mpi_success) then
      write(iulog,*)'mpi_unpack failed ier=',ier
      call endrun
   end if
 
   return
   end subroutine mpiunpack
 
!****************************************************************
 
   subroutine mpisendrecv (sendbuf, sendcount, sendtype, dest, sendtag,  &
                           recvbuf, recvcount, recvtype, source,recvtag, &
                           comm)
!
! Blocking send and receive.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real(r8), intent(in):: sendbuf(*)
   real(r8), intent(out):: recvbuf(*)
   integer, intent(in):: sendcount
   integer, intent(in):: sendtype
   integer, intent(in):: dest
   integer, intent(in):: sendtag
   integer, intent(in):: recvcount
   integer, intent(in):: recvtype
   integer, intent(in):: source
   integer, intent(in):: recvtag
   integer, intent(in):: comm
 
   integer :: status(MPI_STATUS_SIZE)
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_sendrecv')
#endif
   call mpi_sendrecv (sendbuf, sendcount, sendtype, dest, sendtag,   &
                      recvbuf, recvcount, recvtype, source, recvtag, &
                      comm, status, ier)
   if (ier.ne.mpi_success) then
      write(iulog,*)'mpi_sendrecv failed ier=',ier
      call endrun
   end if
!
! ASSUME nrecv = nsend for stats gathering purposes.  This is not actually
! correct, but its the best we can do since recvcount is a Max number
!
   nsend = nsend + 1
   nrecv = nrecv + 1
   nwsend = nwsend + sendcount
   nwrecv = nwrecv + sendcount

#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_sendrecv')
#endif
 
   return
   end subroutine mpisendrecv
 
!****************************************************************
 
   subroutine mpiisend (buf, count, datatype, dest, tag, comm, request)
!
! Does a non-blocking send.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_isend')
#endif
   call mpi_isend (buf, count, datatype, dest, tag, comm, request, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_isend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_isend')
#endif
 
   return
   end subroutine mpiisend
 
!****************************************************************
 
   subroutine mpiirsend (buf, count, datatype, dest, tag, comm, request)
!
! Does a non-blocking ready send.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_irsend')
#endif
   call mpi_irsend (buf, count, datatype, dest, tag, comm, request, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_irsend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_irsend')
#endif
 
   return
   end subroutine mpiirsend
 
!****************************************************************
 
   subroutine mpiissend (buf, count, datatype, dest, tag, comm, request)
!
! Does a non-blocking synchronous send.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_issend')
#endif
   call mpi_issend (buf, count, datatype, dest, tag, comm, request, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_issend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_issend')
#endif
 
   return
   end subroutine mpiissend
 
!****************************************************************
 
   subroutine mpiirecv (buf, count, datatype, source, tag, comm, request)
!
! Does a non-blocking receive.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(out):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: source
   integer, intent(in):: tag
   integer, intent(in):: comm
   integer, intent(out):: request
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_irecv')
#endif
   call mpi_irecv (buf, count, datatype, source, tag, comm, request, ier )
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_irecv failed ier=',ier
      call endrun
   end if
   nrecv = nrecv + 1
   nwrecv = nwrecv + count
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_irecv')
#endif
 
   return
   end subroutine mpiirecv
 
!****************************************************************
 
   subroutine mpiwait (request, status)
!
! Waits for a nonblocking operation to complete.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   integer, intent(inout):: request
   integer, intent(out):: status
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_wait')
#endif
   call mpi_wait (request, status, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_wait failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_wait')
#endif
 
   return
   end subroutine mpiwait
 
!****************************************************************
 
   subroutine mpiwaitall (count, array_of_requests, array_of_statuses)
!
! Waits for a collection of nonblocking operations to complete.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   integer, intent(in):: count
   integer, intent(inout):: array_of_requests(*)
   integer, intent(out):: array_of_statuses(*)
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_waitall')
#endif
   call mpi_waitall (count, array_of_requests, array_of_statuses, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_waitall failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_waitall')
#endif
 
   return
   end subroutine mpiwaitall
 
!****************************************************************
 
   subroutine mpisend (buf, count, datatype, dest, tag, comm)
!
! Does a blocking send
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_send')
#endif
   call mpi_send (buf, count, datatype, dest, tag, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_send failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_send')
#endif
 
   return
   end subroutine mpisend
 
!****************************************************************
 
   subroutine mpirsend (buf, count, datatype, dest, tag, comm)
!
! Does a blocking ready send
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_rsend')
#endif
   call mpi_rsend (buf, count, datatype, dest, tag, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_rsend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_rsend')
#endif
 
   return
   end subroutine mpirsend
 
!****************************************************************
 
   subroutine mpissend (buf, count, datatype, dest, tag, comm)
!
! Does a blocking synchronous send
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: dest
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_ssend')
#endif
   call mpi_ssend (buf, count, datatype, dest, tag, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_ssend failed ier=',ier
      call endrun
   end if
   nsend = nsend + 1
   nwsend = nwsend + count
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_ssend')
#endif
 
   return
   end subroutine mpissend
 
!****************************************************************
 
   subroutine mpirecv (buf, count, datatype, source, tag, comm)
!
! Does a blocking receive
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(out):: buf(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: source
   integer, intent(in):: tag
   integer, intent(in):: comm
 
   integer status (MPI_STATUS_SIZE) ! Status of message
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_recv')
#endif
   call mpi_recv (buf, count, datatype, source, tag, comm, status, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_recv failed ier=',ier
      call endrun
   end if
   nrecv = nrecv + 1
   nwrecv = nwrecv + count
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_recv')
#endif
 
   return
   end subroutine mpirecv
 
!****************************************************************
 
   subroutine mpigather (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
                         recvtype, root, comm)
!
! Collects different messages from each thread on masterproc
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in):: sendbuf(*)
   real (r8), intent(out):: recvbuf(*)
   integer, intent(in):: sendcnt
   integer, intent(in):: sendtype
   integer, intent(in):: recvcnt
   integer, intent(in):: recvtype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_gather')
#endif
   call mpi_gather (sendbuf, sendcnt, sendtype,                      &
                    recvbuf, recvcnt, recvtype, root, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_gather failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_gather')
#endif
 
   return
   end subroutine mpigather
 
!****************************************************************
 
   subroutine mpigatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, &
                          displs, recvtype, root, comm)
!
! Collects different messages from each thread on masterproc
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils, only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
 
   integer ier   ! MPI error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_gatherv')
#endif
   call mpi_gatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, &
                     root, comm, ier)
   if (ier /= mpi_success) then
      write(iulog,*)'mpi_gatherv failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_gatherv')
#endif

   return
   end subroutine mpigatherv
 
!****************************************************************
 
   subroutine mpigathervr4 (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, &
                          displs, recvtype, root, comm)
!
! Collects different messages from each thread on masterproc
!
   use shr_kind_mod, only: r4 => shr_kind_r4, r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils, only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r4), intent(in)  :: sendbuf(*)
   real (r4), intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
 
   integer ier   ! MPI error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_gatherv')
#endif
   call mpi_gatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, &
                     root, comm, ier)
   if (ier /= mpi_success) then
      write(iulog,*)'mpi_gatherv failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_gatherv')
#endif

   return
   end subroutine mpigathervr4
 
!****************************************************************
 
   subroutine mpigathervint (sendbuf, sendcnt, sendtype, recvbuf, &
                             recvcnts, displs, recvtype, root, comm)
!
! Collects different messages from each thread on masterproc
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils, only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   integer, intent(in)  :: sendbuf(*)
   integer, intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
 
   integer ier   ! MPI error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_gatherv')
#endif
   call mpi_gatherv (sendbuf, sendcnt, sendtype, recvbuf, recvcnts, displs, recvtype, &
                     root, comm, ier)
   if (ier /= mpi_success) then
      write(iulog,*)'mpi_gatherv failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_gatherv')
#endif

   return
   end subroutine mpigathervint
 
!****************************************************************
 
   subroutine mpisum (sendbuf, recvbuf, cnt, datatype, root, comm)
!
! Sums sendbuf across all processors on communicator, returning 
! result to root.
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in):: sendbuf(*)
   real (r8), intent(out):: recvbuf(*)
   integer, intent(in):: cnt
   integer, intent(in):: datatype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_reduce')
#endif
   call mpi_reduce (sendbuf, recvbuf, cnt, datatype, mpi_sum, &
                    root, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_reduce failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_reduce')
#endif
 
   return
   end subroutine mpisum
 
!****************************************************************
 
   subroutine mpiscatter (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
                          recvtype, root, comm)
!
! Sends different messages from masterproc to each thread
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8),intent(in):: sendbuf(*)
   real (r8), intent(out):: recvbuf(*)
   integer,intent(in):: sendcnt
   integer,intent(in):: sendtype
   integer,intent(in):: recvcnt
   integer,intent(in):: recvtype
   integer,intent(in):: root
   integer,intent(in):: comm
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_scatter')
#endif
   call mpi_scatter (sendbuf, sendcnt, sendtype, recvbuf, recvcnt, &
                     recvtype, root, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_scatter failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_scatter')
#endif
 
   return
   end subroutine mpiscatter
 
!****************************************************************
 
   subroutine mpiscatterv (sendbuf, sendcnts, displs, sendtype, recvbuf, &
                           recvcnt, recvtype, root, comm)
!
! Sends different messages from masterproc to each thread
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: displs(*)
   integer, intent(in) :: sendcnts(*)
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnt
   integer, intent(in) :: recvtype
   integer, intent(in) :: root
   integer, intent(in) :: comm
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_scatter')
#endif
   call mpi_scatterv (sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, &
                      recvtype, root, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_scatter failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_scatter')
#endif
 
   return
   end subroutine mpiscatterv
 
!****************************************************************
 
   subroutine mpibcast (buffer, count, datatype, root, comm )
!
! Broadcasts a message from masterproc to all threads
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(inout):: buffer(*)
   integer, intent(in):: count
   integer, intent(in):: datatype
   integer, intent(in):: root
   integer, intent(in):: comm
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_bcast')
#endif
   call mpi_bcast (buffer, count, datatype, root, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_bcast failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_bcast')
#endif
 
   return
   end subroutine mpibcast
!****************************************************************
 
   subroutine mpiallmaxint (sendbuf, recvbuf, count, comm)
!
! Allreduce integer vector maximum
! 
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   integer, intent(in)  :: sendbuf(*)
   integer, intent(out) :: recvbuf(*)
   integer, intent(in)  :: count
   integer, intent(in)  :: comm
 
   integer :: ier              ! MPI error code

#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_allreduce')
#endif
   call mpi_allreduce (sendbuf, recvbuf, count, mpiint, &
                       mpimax, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_allreduce failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_allreduce')
#endif

   return
   end subroutine mpiallmaxint

!****************************************************************
 
   subroutine mpialltoallv (sendbuf, sendcnts, sdispls, sendtype, &
                            recvbuf, recvcnts, rdispls, recvtype, &
                            comm)
!
! All-to-all scatter/gather
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: sdispls(*)
   integer, intent(in) :: sendcnts(*)
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: rdispls(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: comm
 
   integer :: ier              ! MPI error code

#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_alltoallv')
#endif
   call mpi_alltoallv (sendbuf, sendcnts, sdispls, sendtype, &
                       recvbuf, recvcnts, rdispls, recvtype, &
                       comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_alltoallv failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_alltoallv')
#endif

   return
   end subroutine mpialltoallv
!****************************************************************
 
   subroutine mpialltoallint (sendbuf, sendcnt, recvbuf, recvcnt, &
                              comm)
!
! All-to-all scatter/gather
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   integer, intent(in)  :: sendbuf(*)
   integer, intent(in)  :: sendcnt
   integer, intent(out) :: recvbuf(*)
   integer, intent(in)  :: recvcnt
   integer, intent(in)  :: comm
 
   integer :: ier              ! MPI error code

#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_alltoallint')
#endif
   call mpi_alltoall (sendbuf, sendcnt, mpiint, &
                      recvbuf, recvcnt, mpiint, &
                      comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_alltoallint failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_alltoallint')
#endif

   return
   end subroutine mpialltoallint

!****************************************************************
 
   subroutine mpiallgatherv (sendbuf, sendcnt, sendtype, &
                             recvbuf, recvcnts, rdispls, recvtype, &
                             comm)
!
! Collect data from each task and broadcast resulting
! vector to all tasks
! 
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real (r8), intent(in)  :: sendbuf(*)
   real (r8), intent(out) :: recvbuf(*)
   integer, intent(in) :: sendcnt
   integer, intent(in) :: sendtype
   integer, intent(in) :: recvcnts(*)
   integer, intent(in) :: rdispls(*)
   integer, intent(in) :: recvtype
   integer, intent(in) :: comm
 
   integer ier   !MP error code
 
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_allgatherv')
#endif
   call mpi_allgatherv (sendbuf, sendcnt, sendtype, &
                        recvbuf, recvcnts, rdispls, recvtype, &
                        comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_allgatherv failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_allgatherv')
#endif
 
   return
   end subroutine mpiallgatherv
!****************************************************************
 
   subroutine mpiallgatherint (sendbuf, scount, recvbuf, rcount, comm)
!
! Collects integer data from each task and broadcasts resulting
! vector to all tasks
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   integer, intent(in)  :: sendbuf(*)
   integer, intent(out) :: recvbuf(*)
   integer, intent(in)  :: scount
   integer, intent(in)  :: rcount
   integer, intent(in)  :: comm
 
   integer ier   !MP error code

#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_allgather')
#endif
   call mpi_allgather (sendbuf, scount, mpiint, recvbuf, rcount, &
                       mpiint, comm, ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_allgather failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_allgather')
#endif
 
   return
   end subroutine mpiallgatherint

!****************************************************************

   subroutine mpiwincreate(base,size,comm,win)
!
! Creates window for MPI2 one-sided commands
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use mpishorthand
   use cam_abortutils,   only: endrun
   use cam_logfile,  only: iulog
#if defined( WRAP_MPI_TIMING )
   use perf_mod
#endif

   implicit none

   real(r8), intent(in)  :: base(*)
   integer,  intent(in)  :: size
   integer,  intent(in)  :: comm
   integer,  intent(out) :: win
!
#ifdef MPI2
   integer(kind=MPI_ADDRESS_KIND) :: size8
   integer :: ier, info
!
#if defined( WRAP_MPI_TIMING )
   call t_startf ('mpi_win_create')
#endif
   info = MPI_INFO_NULL
   size8 = size
   call mpi_win_create(base,size8,8,info,comm,win,ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_win_create failed ier=',ier
      call endrun
   end if
   call mpi_win_fence(0,win,ier)
   if (ier/=mpi_success) then
      write(iulog,*)'mpi_win_fence failed ier=',ier
      call endrun
   end if
#if defined( WRAP_MPI_TIMING )
   call t_stopf ('mpi_win_create')
#endif
#endif

   return
   end subroutine mpiwincreate
!****************************************************************
!
! If SPMD is not turned on
!
#else
   subroutine wrap_mpi
   use cam_abortutils, only: endrun
   implicit none
!
! A unused stub routine to make the compiler happy when SPMD is
! turned off (which means you don't need anything in this file).
!
   call endrun ('(WRAP_MPI): This should not be called at all')
   end subroutine wrap_mpi
#endif

