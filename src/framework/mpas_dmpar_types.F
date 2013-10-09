! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!-----------------------------------------------------------------------
!  mpas_dmpar_types
!
!> \brief MPAS Communication Type Definitions
!> \author Michael Duda, Doug Jacobsen
!> \date   03/26/13
!> \details
!>  This module defines all communication related derived data types
!
!-----------------------------------------------------------------------
module mpas_dmpar_types

   use mpas_kind_types

   type dm_info
     integer :: nprocs, my_proc_id, comm, info
     logical :: using_external_comm
   end type dm_info


   type mpas_exchange_list
     integer :: endPointID
     integer :: nlist
     integer, dimension(:), pointer :: srcList
     integer, dimension(:), pointer :: destList
     type (mpas_exchange_list), pointer :: next
   end type mpas_exchange_list


   type mpas_exchange_list_pointer
     type (mpas_exchange_list), pointer :: exchList
   end type mpas_exchange_list_pointer


   type mpas_multihalo_exchange_list
     type (mpas_exchange_list_pointer), dimension(:), pointer :: halos
   end type mpas_multihalo_exchange_list


   type mpas_communication_list
     integer :: procID
     integer :: nlist
     real (kind=RKIND), dimension(:), pointer :: rbuffer
     integer, dimension(:), pointer :: ibuffer
     integer :: reqID
     type (mpas_communication_list), pointer :: next
   end type mpas_communication_list


end module mpas_dmpar_types
