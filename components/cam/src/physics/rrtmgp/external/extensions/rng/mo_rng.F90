! Module: mo_rng

! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!

module mo_rng
  implicit none

  ! --- abstract base class: random number generator
  ! used as an abstract class
  type, abstract, public :: ty_rng
  contains
    procedure(abstract_get_random_vec),      deferred :: get_random_vec
    procedure(abstract_get_random_vec_mask), deferred :: get_random_vec_mask
    procedure(abstract_init_rng),            deferred :: init_rng
    procedure(abstract_end_rng),             deferred :: end_rng
    generic, public :: get_random => get_random_vec, get_random_vec_mask
    generic, public :: init       => init_rng
    generic, public :: end        => end_rng
  end type

  ! --- Deferred procedures
  !   These define the interface that can be expected and need to be implmented in any derived class
  !
  abstract interface
    ! -----------------------------------------------------------------------------------
    ! Provide num random numbers following a uniform distribution between 0 and 1
    !
    function abstract_get_random_vec(this, num)
      use mo_rte_kind, only : dp
      import ty_rng
      class(ty_rng)        :: this
      integer,  intent(in) :: num
      real(DP), dimension(num) :: abstract_get_random_vec
    end function abstract_get_random_vec
    ! -----------------------------------------------------------------------------------
    ! Provide random numbers for the TRUE elements of MASK
    !
    function abstract_get_random_vec_mask(this, mask)
      use mo_rte_kind, only : dp
      import ty_rng
      class(ty_rng)                     :: this
      logical, dimension(:), intent(in) :: mask
      real(DP), dimension(size(mask)) :: abstract_get_random_vec_mask
    end function abstract_get_random_vec_mask
    ! -----------------------------------------------------------------------------------
    ! Initialize the random number state.
    !
    subroutine abstract_init_rng(this, seeds)
      import ty_rng
      class(ty_rng) :: this
      integer,  dimension(:), intent(in) :: seeds
    end subroutine abstract_init_rng
    ! -----------------------------------------------------------------------------------
    ! Release any resources associated with the RNG. Could do nothing.
    !
    subroutine abstract_end_rng(this)
      import ty_rng
      class(ty_rng) :: this
    end subroutine abstract_end_rng
    ! -----------------------------------------------------------------------------------
  end interface

end
