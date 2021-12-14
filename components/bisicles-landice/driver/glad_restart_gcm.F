!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_restart_gcm.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glad_restart_gcm

!BOP
! !MODULE: glad_restart_gcm

! !DESCRIPTION:
!  Contains routines for specialized restarts called by GCMs
!
! !REVISION HISTORY:
!
! !USES:

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: glad_read_restart_gcm

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
! !IROUTINE: glad_read_restart_gcm
! !INTERFACE:

   subroutine glad_read_restart_gcm(model, restart_filename)

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

  end subroutine glad_read_restart_gcm

!-----------------------------------------------------------------------

end module glad_restart_gcm

!-----------------------------------------------------------------------
