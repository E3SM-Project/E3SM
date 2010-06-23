program pioasync
  use pio   !_EXTERNAL
  use mpi   !_EXTERNAL
  use pio_support, only : CheckMPIReturn  ! _EXTERNAL

  implicit none

  integer :: ierr, my_task(3), nprocs(3)
  integer :: mpigrp_world, mpigrp_compute, mpigrp_io
  integer :: pelist(3,1), val, ioroot

  integer :: mpi_comm_compute, mpi_comm_io, mpi_icomm_cio
  type(iosystem_desc_t) :: iosystem


  call mpi_init(ierr)
  call CheckMPIReturn('Call to MPI_INIT()',ierr,__FILE__,__LINE__)

  nprocs = -1
  my_task = -1

  
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_task(1),ierr)
  call CheckMPIReturn('Call to MPI_COMM_RANK()',ierr,__FILE__,__LINE__)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs(1),ierr)
  call CheckMPIReturn('Call to MPI_COMM_SIZE()',ierr,__FILE__,__LINE__)
  

  ! Split into compute and io comms

  mpi_comm_compute = MPI_COMM_NULL
  mpi_comm_io = MPI_COMM_NULL

#ifdef DOTHIS 
  if(my_task(1) < nprocs(1)/2) then
     call mpi_comm_split(mpi_comm_world, 0, 0, mpi_comm_compute, ierr)
  else
     call mpi_comm_split(mpi_comm_world, 1, 0, mpi_comm_io, ierr)
  end if
 
#else


  pelist(1,1) = 0
  pelist(2,1) = nprocs(1)/2-1
  pelist(3,1) = 1

  call mpi_comm_group(mpi_comm_world, mpigrp_world, ierr)
     call CheckMPIReturn('Call to MPI_COMM_group()',ierr,__FILE__,__LINE__)

  call mpi_group_range_incl(mpigrp_world, 1, pelist, mpigrp_compute, ierr)
     call CheckMPIReturn('Call to MPI_group_range_incl()',ierr,__FILE__,__LINE__)

  call mpi_comm_create(mpi_comm_world, mpigrp_compute, mpi_comm_compute, ierr)
     call CheckMPIReturn('Call to MPI_COMM_create()',ierr,__FILE__,__LINE__)


  pelist(1,1) = nprocs(1)/2
  pelist(2,1) = nprocs(1)-1
  pelist(3,1) = 1

  print *, pelist


  call mpi_group_range_incl(mpigrp_world, 1, pelist, mpigrp_io, ierr)

 call CheckMPIReturn('Call to MPI_group_range_incl()',ierr,__FILE__,__LINE__)
  call mpi_comm_create(mpi_comm_world, mpigrp_io, mpi_comm_io, ierr)

 call CheckMPIReturn('Call to MPI_COMM_create()',ierr,__FILE__,__LINE__)
#endif

 print *, mpi_comm_world, mpi_comm_compute, mpi_comm_io

 call mpi_barrier(mpi_comm_world, ierr)



  if(my_task(1) < nprocs(1)/2) then
     call MPI_COMM_RANK(MPI_COMM_compute,my_task(2),ierr)
     call CheckMPIReturn('Call to MPI_COMM_RANK()',ierr,__FILE__,__LINE__)
     call MPI_COMM_SIZE(MPI_COMM_compute,nprocs(2),ierr)
     call CheckMPIReturn('Call to MPI_COMM_SIZE()',ierr,__FILE__,__LINE__)

     call mpi_intercomm_create(mpi_comm_compute, 0, mpi_comm_world, nprocs(1)/2, 1, mpi_icomm_cio, ierr)
     ioroot = 0
     val = -1
  else
     call MPI_COMM_RANK(MPI_COMM_IO,my_task(3),ierr)
     call CheckMPIReturn('Call to MPI_COMM_RANK()',ierr,__FILE__,__LINE__)
     call MPI_COMM_SIZE(MPI_COMM_IO,nprocs(3),ierr)
     call CheckMPIReturn('Call to MPI_COMM_SIZE()',ierr,__FILE__,__LINE__)

     call mpi_intercomm_create(mpi_comm_io, 0, mpi_comm_world, 0, 1, mpi_icomm_cio, ierr)

     if(my_task(3)==0) then
        val = 18
        ioroot = MPI_ROOT
     else
        val = -1
        ioroot = MPI_PROC_NULL
     end if


  end if

  print *, __LINE__,my_task, ierr
  print *, __LINE__,nprocs

! A simple check of the intercomm - val should be 18 on mpi_comm_compute and -1 on mpi_comm_io except on the 
! root task of that comm

  call mpi_bcast(val, 1, mpi_integer, ioroot, mpi_icomm_cio, ierr)

  print *, __LINE__,my_task(1), val


  call pio_setdebuglevel(3)

!  This flavor of pio_init will take the compute comm and the io comm
!  and create an internal intercomm

  call pio_init(mpi_comm_compute, mpi_comm_io, mpi_icomm_cio, iosystem)

! only the tasks in mpi_comm_compute should return from pio_init - really?


  call piotest(iosystem)



  call mpi_finalize(ierr)


end program pioasync

subroutine piotest(iosystem)
  use mpi !_EXTERNAL
  use pio !_EXTERNAL
  implicit none
  type(iosystem_desc_t) :: iosystem
  type(file_desc_t) :: file
  character(len=20) :: fname='afilename.nc'
  integer :: ierr, msg=999

  ierr= pio_createfile(iosystem, file, pio_iotype_netcdf, fname)

  print *,__FILE__,__LINE__, iosystem%compmaster, iosystem%intercomm


  call mpi_bcast(msg, 1, mpi_integer, iosystem%compmaster, iosystem%intercomm, ierr)

  print *,__FILE__,__LINE__


end subroutine piotest
