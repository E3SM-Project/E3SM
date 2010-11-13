module utils_mod

  use pio   ! _EXTERNAL
  use kinds_mod
  implicit none
  private

  public :: WriteHeader, split_comm

contains

!>
!! @private
!! @brief Writes netcdf header information for testpio. 
!! @param File @copydoc file_desc_t
!! @param nx
!! @param ny
!! @param nz
!! @param dimid_x
!! @param dimid_y
!! @param dimid_z
!<
subroutine WriteHeader(File,nx,ny,nz,dimid_x,dimid_y,dimid_z)

       type (File_desc_t), intent(inout) :: File
       integer(i4), intent(in) :: nx,ny,nz
       integer(i4), intent(out) :: dimid_x,dimid_y,dimid_z

       integer(i4) :: itmp,iostat

       iostat = PIO_put_att(File,pio_global,'title','Test NetCDF file')
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error writing TITLE to netCDF file'
       endif

       iostat = PIO_put_att(File,pio_global,'ivalue', 4)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error writing iVALUE to netCDF file'
       endif

       iostat = PIO_def_dim(File,'X',nx,dimid_x)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension X for netCDF file'
       endif

       iostat = PIO_def_dim(File,'Y',ny,dimid_y)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension Y for netCDF file'
       endif

       iostat = PIO_def_dim(File,'Z',nz,dimid_z)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension Z for netCDF file'
       endif

end subroutine WriteHeader

subroutine ReadHeader(File,nx,ny,nz,dimid_x,dimid_y,dimid_z)
!! NOTE: This subroutine is screwed up, it should be get, not put (tcraig 2/18/09)

       type (File_desc_t), intent(inout) :: File
       integer(i4), intent(in) :: nx,ny,nz
       integer(i4), intent(out) :: dimid_x,dimid_y,dimid_z

       integer(i4) :: itmp,iostat

       iostat = PIO_put_att(File,pio_global,'title','Test NetCDF file')
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error writing TITLE to netCDF file'
       endif

       itmp = 4
       iostat = PIO_put_att(File,pio_global,'ivalue',(/itmp/))
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error writing iVALUE to netCDF file'
       endif

       iostat = PIO_def_dim(File,'X',nx,dimid_x)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension X for netCDF file'
       endif

       iostat = PIO_def_dim(File,'Y',ny,dimid_y)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension Y for netCDF file'
       endif

       iostat = PIO_def_dim(File,'Z',nz,dimid_z)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension Z for netCDF file'
       endif

end subroutine ReadHeader

subroutine split_comm(initial_comm, nprocs, num_iotasks, stride, base, mpi_comm_compute, mpi_comm_io, intercomm)
  use pio_support !_EXTERNAL
  implicit none

  integer, intent(in) :: initial_comm, nprocs, num_iotasks, stride, base
  integer, intent(out) :: mpi_comm_compute, mpi_comm_io, intercomm

  integer :: ierr
  integer :: pelist(3,1), mpigrp_init, mpigrp_io, mpigrp_compute
  include 'mpif.h' !_EXTERNAL
#ifndef _MPISERIAL
  mpi_comm_compute = MPI_COMM_NULL
  mpi_comm_io = MPI_COMM_NULL

  pelist(1,1) = base
  pelist(2,1) = min(nprocs-1,num_iotasks*stride)
  pelist(3,1) = stride

  call mpi_comm_group(initial_comm, mpigrp_init, ierr)

  call mpi_group_range_incl(mpigrp_init, 1, pelist, mpigrp_io, ierr)

  call mpi_group_range_excl(mpigrp_init, 1, pelist, mpigrp_compute, ierr)

  call mpi_comm_create(initial_comm, mpigrp_compute, mpi_comm_compute, ierr)

  call mpi_comm_create(initial_comm, mpigrp_io, mpi_comm_io, ierr)

  if(mpi_comm_compute/=MPI_COMM_NULL) then
     call mpi_intercomm_create(mpi_comm_compute, 0, initial_comm, base, 1, intercomm, ierr)
  else if(mpi_comm_io/=MPI_COMM_NULL) then
     if(base==0) then
        if(stride>1) then
           call mpi_intercomm_create(mpi_comm_io, 0, initial_comm, 1, 1, intercomm, ierr)
        else
           call mpi_intercomm_create(mpi_comm_io, 0, initial_comm, num_iotasks, 1, intercomm, ierr)
        end if
     else
        call mpi_intercomm_create(mpi_comm_io, 0, initial_comm, 0, 1, intercomm, ierr)
     end if
  else
    call piodie(__FILE__,__LINE__)
  end if
#endif
end subroutine split_comm




end module utils_mod
