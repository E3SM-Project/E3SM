! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_framework_forcing
!
!> \brief Declare the forcing group linked list head
!> \author Adrian K. Turner
!> \date   19th May 2017
!> \details
!>  Declare the forcing group linked list head
!
!-----------------------------------------------------------------------

module ocn_framework_forcing

  use mpas_derived_types
  use mpas_forcing

  implicit none
  private

  type(MPAS_forcing_group_type), pointer, public :: &
       forcingGroupHead

end module ocn_framework_forcing
