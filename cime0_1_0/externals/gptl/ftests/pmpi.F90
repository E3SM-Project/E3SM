#undef DEBUG
module myvars
  integer :: iam
  integer :: commsize
end module myvars

program pmpi
  use myvars
  use gptl

  implicit none

#include <mpif.h>

  integer, parameter :: tag = 98
  integer, parameter :: count = 100000

  integer :: i, j, ret
  integer :: val
  integer :: comm = MPI_COMM_WORLD
  integer :: sendbuf(0:count-1)
  integer :: recvbuf(0:count-1)
  integer :: sum
  integer :: status(MPI_STATUS_SIZE)
  integer :: sendreq, recvreq
  integer :: dest
  integer :: source
  integer :: rdispls(0:count-1)
  integer :: sdispls(0:count-1)

  integer :: kount, of           ! for gptlquery
  real(8) :: wc, usr, sys ! for gptlquery
  integer(8) :: pc               ! for gptlquery

  integer, allocatable :: atoabufsend(:)
  integer, allocatable :: atoabufrecv(:)
  integer, allocatable :: gsbufsend(:,:)      ! gather/scatter buffer send
  integer, allocatable :: gsbufrecv(:,:)      ! gather/scatter buffer recv
  integer, allocatable :: recvcounts(:)
  integer, allocatable :: sendcounts(:)
  integer, allocatable :: atoacounts(:)
  integer, allocatable :: atoadispls(:)

  logical :: flag
  integer :: debugflag = 1

  ret = gptlsetoption (gptloverhead, 0)
  ret = gptlsetoption (gptlpercent, 0)
  ret = gptlsetoption (gptlabort_on_error, 1)
  ret = gptlsetoption (gptlsync_mpi, 1)

#if ( ! defined HAVE_IARGCGETARG )
  ret = gptlinitialize ()
  ret = gptlstart ("total")
#endif

  call mpi_init (ret)

! For debugging, go into infinite loop so debugger can attach and reset
#ifdef DEBUG
  do while (debugflag == 1)
  end do
#endif

  call mpi_comm_rank (comm, iam, ret)
  call mpi_comm_size (comm, commsize, ret)
  if (iam == 0) write(6,*)'commsize is ', commsize

  do i=0,count-1
    sendbuf(i) = iam
  end do

  dest = mod ((iam + 1), commsize)
  source = iam - 1
  if (source < 0) then
    source = commsize - 1
  end if
!
! mpi_send
! mpi_recv
! mpi_probe
!
  recvbuf(:) = -1
  if (mod (commsize, 2) == 0) then
    if (iam == 0) then
      write(6,*)'Testing send, recv, probe...'
    end if

    if (mod (iam, 2) == 0) then
      call mpi_send (sendbuf, count, MPI_INTEGER, dest, tag, comm, ret)
      call mpi_recv (recvbuf, count, MPI_INTEGER, source, tag, comm, status, ret)
    else
      call mpi_probe (source, tag, comm, status, ret)
      if (ret /= MPI_SUCCESS) then
        write(6,*) "iam=", iam, " mpi_probe: bad return"
        call mpi_abort (MPI_COMM_WORLD, -1, ret)
      end if
      call mpi_recv (recvbuf, count, MPI_INTEGER, source, tag, comm, status, ret)
      call mpi_send (sendbuf, count, MPI_INTEGER, dest, tag, comm, ret)
    end if
    call chkbuf ('mpi_send + mpi_recv', recvbuf(:), count, source)

    if (iam == 0) then
      write(6,*)'Success'
      write(6,*)'Testing ssend...'
    end if
!
! mpi_ssend
!
    recvbuf(:) = -1
    if (mod (iam, 2) == 0) then
      call mpi_ssend (sendbuf, count, MPI_INTEGER, dest, tag, comm, ret)
      call mpi_recv (recvbuf, count, MPI_INTEGER, source, tag, comm, status, ret)
    else
      call mpi_recv (recvbuf, count, MPI_INTEGER, source, tag, comm, status, ret)
      call mpi_ssend (sendbuf, count, MPI_INTEGER, dest, tag, comm, ret)
    end if
    call chkbuf ('mpi_send + mpi_recv', recvbuf(:), count, source)
    if (iam == 0) then
      write(6,*)'Success'
      write(6,*)'Testing sendrecv...'
    end if
  else
    if (iam == 0) write(6,*)'NOTE: commsize=',commsize,' is odd so wont test ', &
                            'send, recv, probe, ssend'
  end if
!
! mpi_sendrecv
! 
  recvbuf(:) = -1
  call mpi_sendrecv (sendbuf, count, MPI_INTEGER, dest, tag, &
                     recvbuf, count, MPI_INTEGER, source, tag, &
                     comm, status, ret)
  call chkbuf ('mpi_sendrecv', recvbuf(:), count, source)
  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing irecv, isend, issend, iprobe, itest, wait, waitall...'
  end if
!
! mpi_irecv
! mpi_isend
! mpi_issend
! mpi_iprobe
! mpi_test
! mpi_wait
! mpi_waitall
!
  recvbuf(:) = -1
  call mpi_irecv (recvbuf, count, MPI_INTEGER, source, tag, &
                  comm, recvreq, ret)
  call mpi_iprobe (source, tag, comm, flag, status, ret)
  call mpi_test (recvreq, flag, status, ret)
  call mpi_isend (sendbuf, count, MPI_INTEGER, dest, tag, &
                  comm, sendreq, ret)
  call mpi_wait (recvreq, status, ret)
  call mpi_wait (sendreq, status, ret)
  call chkbuf ("mpi_wait", recvbuf(:), count, source)

  recvbuf(:) = -1
  call mpi_irecv (recvbuf, count, MPI_INTEGER, source, tag, &
                  comm, recvreq, ret)
  call mpi_iprobe (source, tag, comm, flag, status, ret)
  call mpi_test (recvreq, flag, status, ret)
  call mpi_issend (sendbuf, count, MPI_INTEGER, dest, tag, &
                   comm, sendreq, ret)
  call mpi_wait (recvreq, status, ret)
  call mpi_wait (sendreq, status, ret)
  call chkbuf ("mpi_wait", recvbuf(:), count, source)

  recvbuf(:) = -1
  call mpi_irecv (recvbuf, count, MPI_INTEGER, source, tag, &
                  comm, recvreq, ret)
  call mpi_isend (sendbuf, count, MPI_INTEGER, dest, tag, &
                  comm, sendreq, ret)
  call mpi_waitall (1, recvreq, status, ret)
  call mpi_waitall (1, sendreq, status, ret)
  call chkbuf ("mpi_waitall", recvbuf(:), count, source)

  call mpi_barrier (comm, ret)
  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing bcast...'
  end if
!
! mpi_bcast
!
  call mpi_bcast (sendbuf, count, MPI_INTEGER, 0, comm, ret)
  call chkbuf ("mpi_bcast", sendbuf(:), count, 0)
  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing allreduce...'
  end if
!
! mpi_allreduce: need to reset sendbuf due to bcast just done
!
  do i=0,count-1
    sendbuf(i) = iam
  end do

  recvbuf(:) = -1
  call mpi_allreduce (sendbuf, recvbuf, count, MPI_INTEGER, MPI_SUM, comm, ret)
  sum = 0.
  do i=0,commsize-1
    sum = sum + i
  end do
  call chkbuf ("mpi_allreduce", recvbuf(:), count, sum)
  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing gather...'
  end if

  allocate (gsbufsend(0:count-1,0:commsize-1))
  allocate (gsbufrecv(0:count-1,0:commsize-1))
  allocate (recvcounts(0:commsize-1))
  allocate (sendcounts(0:commsize-1))
!
! mpi_gather
!
  gsbufrecv(:,:) = -1
  call mpi_gather (sendbuf, count, MPI_INTEGER, &
                   gsbufrecv, count, MPI_INTEGER, 0, comm, ret)
  if (iam == 0) then
    do j=1,commsize-1
      call chkbuf ("mpi_gather", gsbufrecv(:,j), count, j)
    end do
    write(6,*)'Success'
    write(6,*)'Testing gatherv...'
  end if
!
! mpi_gatherv: make just like mpi_gather for simplicity
!
  gsbufrecv(:,:) = -1
  recvcounts(:) = count
  rdispls(0) = 0
  do j=1,commsize-1
    rdispls(j) = rdispls(j-1) + recvcounts(j-1)
  end do
  call mpi_gatherv (sendbuf, count, MPI_INTEGER, &
                    gsbufrecv, recvcounts, rdispls, &
                    MPI_INTEGER, 0, comm, ret)
  if (iam == 0) then
    do j=1,commsize-1
      call chkbuf ("mpi_gatherv", gsbufrecv(:,j), count, j)
    end do
    write(6,*)'Success'
    write(6,*)'Testing scatter...'
  end if
!
! mpi_scatter
!
  if (iam == 0) then
    do j=0,commsize-1
      gsbufsend(:,j) = j
    end do
  else
    do j=0,commsize-1
      gsbufsend(:,j) = -1
    end do
  end if
  recvbuf(:) = -1
  call mpi_scatter (gsbufsend, count, MPI_INTEGER, recvbuf, count, MPI_INTEGER, &
                    0, comm, ret)
  call chkbuf ("mpi_scatter", recvbuf(:), count, iam)
  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing scatterv...'
  end if
!
! mpi_scatterv: make just like mpi_scatter for simplicity.
!
  if (iam == 0) then
    do j=0,commsize-1
      gsbufsend(:,j) = j
    end do
  else
    gsbufsend(:,:) = -1
  end if
  sendcounts(:) = count
  sdispls(0) = 0
  do j=1,commsize-1
    sdispls(j) = sdispls(j-1) + sendcounts(j-1)
  end do
  recvbuf(:) = -1
  call mpi_scatterv (gsbufsend, sendcounts, sdispls, &
                     MPI_INTEGER, recvbuf, count, &
                     MPI_INTEGER, 0, comm, ret)
  call chkbuf ("mpi_scatterv", recvbuf(:), count, iam)
  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing alltoall...'
  end if
!
! mpi_alltoall
!
  allocate (atoabufsend(0:commsize-1))
  allocate (atoabufrecv(0:commsize-1))
  allocate (atoacounts(0:commsize-1))
  allocate (atoadispls(0:commsize-1))
  do j=0,commsize-1
    atoabufsend(j) = j
  end do
  atoabufrecv(:) = -1
  call mpi_alltoall (atoabufsend, 1, MPI_INTEGER, atoabufrecv, 1, MPI_INTEGER, comm, ret)
  call chkbuf ("mpi_alltoall", atoabufrecv(:), 1, iam)
  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing alltoallv...'
  end if
!
! mpi_alltoallv
!
  atoabufrecv(:) = -1
  atoacounts(:) = 1
  atoadispls(0) = 0
  do j=1,commsize-1
    atoadispls(j) = atoadispls(j-1) + atoacounts(j-1)
  end do
  
  call mpi_alltoallv (atoabufsend, atoacounts, atoadispls, MPI_INTEGER, &
                      atoabufrecv, atoacounts, atoadispls, MPI_INTEGER, comm, ret)
  call chkbuf ("mpi_alltoall", atoabufrecv(:), 1, iam)

  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing reduce...'
  end if
!
! mpi_reduce
!
  call mpi_reduce (sendbuf, recvbuf, count, MPI_INTEGER, MPI_SUM, 0, comm, ret)
  if (iam == 0) then
    sum = 0.
    do i=0,commsize-1
      sum = sum + i
    end do
    call chkbuf ("mpi_reduce", recvbuf(:), count, sum)
  end if

  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing allgather...'
  end if
!
! mpi_allgather
!
  gsbufrecv(:,:) = -1
  call mpi_allgather (sendbuf, count, MPI_INTEGER, &
                      gsbufrecv, count, MPI_INTEGER, comm, ret)
  do j=0,commsize-1
    call chkbuf ("mpi_allgather", gsbufrecv(:,j), count, j)
  end do

  if (iam == 0) then
    write(6,*)'Success'
    write(6,*)'Testing allgatherv...'
  end if
!
! mpi_allgatherv: Make just like mpi_allgather for simplicity
!
  gsbufrecv(:,:) = -1
  recvcounts(:) = count
  call mpi_allgatherv (sendbuf, count, MPI_INTEGER, &
                       gsbufrecv, recvcounts, rdispls, &
                       MPI_INTEGER, comm, ret)
  do j=0,commsize-1
    call chkbuf ("mpi_allgatherv", gsbufrecv(:,j), count, j)
  end do

  if (iam == 0) then
    write(6,*)'Success. Calling finalize'
  end if
!
! mpi_finalize
!
  call mpi_finalize (ret)

#if ( defined HAVE_IARGCGETARG )
  if (iam == 0) then
    write(6,*)'Testing for auto-generated MPI_Init_thru_Finalize region...'
    ret = gptlquery ('MPI_Init_thru_Finalize', 0, kount, of, wc, usr, sys, pc, 0)
    if (ret == 0) then
      write(6,*)'Success'
    else
      write(6,*)'Failure'
    end if
  end if
#else
  ret = gptlstop ("total")
  ret = gptlpr (iam)
#endif

  stop 0
end program pmpi

subroutine chkbuf (msg, recvbuf, count, val)
  use myvars
  implicit none

#include <mpif.h>

  character(len=*), intent(in) :: msg

  integer, intent(in) :: count
  integer, intent(in) :: recvbuf(0:count-1)
  integer, intent(in) :: val

  integer :: i
  integer :: ret
  do i=0,count-1
    if (recvbuf(i) /= val) then
      write(6,*) "iam=", iam, msg, " bad recvbuf(", i,")=",recvbuf(i), "/= ", val
      call mpi_abort (MPI_COMM_WORLD, -1, ret)
    end if
  end do
end subroutine chkbuf
