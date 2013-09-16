! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!-----------------------------------------------------------------------
!  mpas_framework
!
!> \brief MPAS Framework routines
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This module contains all routines related to the general MPAS framework interface.
!
!-----------------------------------------------------------------------
module mpas_framework

   use mpas_dmpar
   use mpas_grid_types
   use mpas_io_input
   use mpas_io_output
   use mpas_configure
   use mpas_timer
   use mpas_timekeeping
   use mpas_io
   use mpas_io_units


   contains

!-----------------------------------------------------------------------
!  routine mpas_framework_init
!
!> \brief MPAS framework initialization routine.
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine initializes the MPAS framework. It calls routines related to initializing different parts of MPAS, that are housed within the framework.
!
!-----------------------------------------------------------------------
   subroutine mpas_framework_init(dminfo, domain, mpi_comm, nml_filename, io_system, calendar, stdoutUnit, stderrUnit)!{{{

      implicit none

      type (dm_info), pointer :: dminfo
      type (domain_type), pointer :: domain
      integer, intent(in), optional :: mpi_comm

      character (len=*), optional :: nml_filename
      type (iosystem_desc_t), optional, pointer :: io_system
      character(len=*), intent(in), optional :: calendar

      integer, intent(in), optional :: stdoutUnit, stderrUnit

      integer :: pio_num_iotasks
      integer :: pio_stride

      allocate(dminfo)
      call mpas_io_units_init(stdoutUnit, stderrUnit)
      call mpas_dmpar_init(dminfo, mpi_comm)

      call mpas_read_namelist(dminfo, nml_filename)

      call mpas_allocate_domain(domain, dminfo)
      
      if(present(calendar)) then
        call mpas_timekeeping_init(calendar)
      else
        call mpas_timekeeping_init(config_calendar_type)
      end if

      pio_num_iotasks = config_pio_num_iotasks
      pio_stride = config_pio_stride
      if (pio_num_iotasks == 0) then
         pio_num_iotasks = domain % dminfo % nprocs
      end if
      call MPAS_io_init(dminfo, pio_num_iotasks, pio_stride, io_system)

   end subroutine mpas_framework_init!}}}

!-----------------------------------------------------------------------
!  routine mpas_framework_finalize
!
!> \brief MPAS framework finalization routine.
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This routine finalizes the MPAS framework. It calls routines related to finalizing different parts of MPAS, that are housed within the framework.
!
!-----------------------------------------------------------------------  
   subroutine mpas_framework_finalize(dminfo, domain, io_system)!{{{
  
      implicit none

      type (dm_info), pointer :: dminfo
      type (domain_type), pointer :: domain
      type (iosystem_desc_t), optional, pointer :: io_system

      call MPAS_io_finalize(io_system)

      call mpas_deallocate_domain(domain)

      call mpas_dmpar_finalize(dminfo)

      call mpas_timekeeping_finalize()

   end subroutine mpas_framework_finalize!}}}

end module mpas_framework
