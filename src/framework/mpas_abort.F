! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_abort

   contains

   !-----------------------------------------------------------------------
   !  routine mpas_dmpar_global_abort
   !
   !> \brief Forces the exit of all processes in MPI_COMM_WORLD
   !> \author Michael Duda
   !> \date   26 March 2013
   !> \details
   !>  This routine writes the specified message to standard error and to 
   !>  a per-process log file named log.????.abort. However, if the optional
   !>  argument deferredAbort is set to .true., messages will be written to
   !>  standard error and to the log.????.abort files, but MPI tasks will not
   !>  be termintated. This allows code to write several lines of messages
   !>  before exiting.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_dmpar_global_abort(mesg, deferredAbort)!{{{
   
      use mpas_kind_types, only : StrKIND
      use mpas_io_units, only : stderrUnit, stdoutUnit, mpas_new_unit
      use mpas_threading, only : mpas_threading_get_thread_num
   
#ifdef _MPI
#ifndef NOMPIMOD
      use mpi
#endif
#endif
   
      implicit none
   
#ifdef _MPI
#ifdef NOMPIMOD
      include 'mpif.h'
#endif
#endif
   
      character(len=*), intent(in) :: mesg !< Input: Abort message
      logical, intent(in), optional :: deferredAbort !< Input: Defer call to abort until later
   
      integer :: threadNum
   
#ifdef _MPI
      integer :: mpi_ierr, mpi_errcode, my_proc_id, nprocs
#endif
   
      character(len=StrKIND) :: errorFile
      integer :: errorUnit
      logical :: local_deferredAbort
   
      if (present(deferredAbort)) then
         local_deferredAbort = deferredAbort
      else
         local_deferredAbort = .false.
      end if
   
      threadNum = mpas_threading_get_thread_num()
   
      if ( threadNum == 0 ) then
         write(stderrUnit,*) trim(mesg)
      end if
   
   
      if (.not. local_deferredAbort) then
         close(stderrUnit)
         close(stdoutUnit)
      end if
   
#ifdef _MPI
      call MPI_Comm_rank(MPI_COMM_WORLD, my_proc_id, mpi_ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nprocs, mpi_ierr)
      if (nprocs < 1E5) then
         write(errorFile,fmt='(a,i4.4,a)') 'log.', my_proc_id, '.abort'
      else if (nprocs < 1E6) then
         write(errorFile,fmt='(a,i5.5,a)') 'log.', my_proc_id, '.abort'
      else if (nprocs < 1E7) then
         write(errorFile,fmt='(a,i6.6,a)') 'log.', my_proc_id, '.abort'
      else if (nprocs < 1E8) then
         write(errorFile,fmt='(a,i7.7,a)') 'log.', my_proc_id, '.abort'
      else if (nprocs < 1E9) then
         write(errorFile,fmt='(a,i8.8,a)') 'log.', my_proc_id, '.abort'
      else
         write(errorFile,fmt='(a,i9.9,a)') 'log.', my_proc_id, '.abort'
      end if
#else
      errorFile = 'log.abort'
#endif
   
      if ( threadNum == 0 ) then
         call mpas_new_unit(errorUnit)
         open(unit=errorUnit, file=trim(errorFile), form='formatted', position='append')
         write(errorUnit,*) trim(mesg)
         close(errorUnit)
      end if
   
      if (.not. local_deferredAbort) then
#ifdef _MPI
         call MPI_Abort(MPI_COMM_WORLD, mpi_errcode, mpi_ierr)
#else
         stop
#endif
      end if
   
   end subroutine mpas_dmpar_global_abort!}}}

end module mpas_abort
