! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
! Routines for permuting arrays: here one (x,y,z) -> (z,x,y) and (x,y,z) -> (z,y,x)
!   Routines are just the front end to kernels
!
! -------------------------------------------------------------------------------------------------
module mo_util_reorder
  use mo_rte_kind, only: wp
  use mo_reorder_kernels, only: reorder_123x312_kernel, reorder_123x321_kernel
  implicit none
  private
  public :: reorder123x312, reorder123x321
contains
  ! -------------------------------------------------------------------------------------------------
  !
  ! (x,y,z) -> (z,x,y)
  !
  function reorder123x312(array)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=1),size(array,dim=2)) :: reorder123x312

    call reorder_123x312_kernel(size(array,dim=1), size(array,dim=2), size(array,dim=3), array, reorder123x312)
  end function reorder123x312
  ! -------------------------------------------------------------------------------------------------
  !
  ! (x,y,z) -> (z,y,x) 
  !
  function reorder123x321(array)
    real(wp), dimension(:,:,:), intent(in) :: array
    real(wp), dimension(size(array,dim=3),size(array,dim=2),size(array,dim=1)) :: reorder123x321

    call reorder_123x321_kernel(size(array,dim=1), size(array,dim=2), size(array,dim=3), array, reorder123x321)
  end function reorder123x321
  ! -------------------------------------------------------------------------------------------------
end module
