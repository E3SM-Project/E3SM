! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Kernels to permute arrays

module mo_reorder_kernels
  use mo_rte_kind,      only: wp
  implicit none
contains
  ! ----------------------------------------------------------------------------
  subroutine reorder_123x312_kernel(d1, d2, d3, array_in, array_out) &
      bind(C, name = "reorder_123x312_kernel")
    integer,                         intent( in) :: d1, d2, d3
    real(wp), dimension(d1, d2, d3), intent( in) :: array_in
    real(wp), dimension(d3, d1, d2), intent(out) :: array_out

    integer :: i1, i2, i3

    do i2 = 1, d2
      do i1 = 1, d1
        do i3 = 1, d3
          array_out(i3,i1,i2) = array_in(i1,i2,i3)
        end do
      end do
    end do
  end subroutine reorder_123x312_kernel
  ! ----------------------------------------------------------------------------
  subroutine reorder_123x321_kernel(d1, d2, d3, array_in, array_out) & 
      bind(C, name="reorder_123x321_kernel")
    integer,                         intent( in) :: d1, d2, d3
    real(wp), dimension(d1, d2, d3), intent( in) :: array_in
    real(wp), dimension(d3, d2, d1), intent(out) :: array_out

    integer :: i1, i2, i3

    do i1 = 1, d1
      do i2 = 1, d2
        do i3 = 1, d3
          array_out(i3,i2,i1) = array_in(i1,i2,i3)
        end do
      end do
    end do
  end subroutine reorder_123x321_kernel
  ! ----------------------------------------------------------------------------
end module mo_reorder_kernels
