!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_restart_gcm.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glimmer_restart_gcm

 !TODO - Should this be moved to the Glint directory?  Used only by glint_initialise.

!BOP
! !MODULE: glimmer_restart_gcm

! !DESCRIPTION:
!  Contains routines for specialized glimmer restarts called by GCMs
!
! !REVISION HISTORY:
!
! !USES:

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: glimmer_read_restart_gcm

!----------------------------------------------------------------------
!
!   module variables
!
!----------------------------------------------------------------------

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glimmer_read_restart_gcm
! !INTERFACE:

   subroutine glimmer_read_restart_gcm(model, restart_filename)

    use glide_types
    implicit none
    type(glide_global_type), intent(inout) :: model
    character(*),            intent(in   ) :: restart_filename

    ! local variables
    type(glimmer_nc_input),  pointer :: ic => null()

    ! create the input unit
    allocate(ic)
    ic%get_time_slice = 1
    ic%nc%filename    = trim(restart_filename)
    ic%nc%vars        = ' restart '
    ic%nc%restartfile = .true.
    ic%nc%vars_copy   = ic%nc%vars

    ! add the input unit to the model
    ! note that the model will do the actual reading of data
    model%funits%in_first => ic

  end subroutine glimmer_read_restart_gcm

!-----------------------------------------------------------------------

end module glimmer_restart_gcm

!-----------------------------------------------------------------------
