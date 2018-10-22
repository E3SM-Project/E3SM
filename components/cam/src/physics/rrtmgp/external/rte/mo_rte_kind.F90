! This code is part of Radiative Transfer for Energetics (RTE)
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

! This module provides the Fortran KIND parameters for REAL and INTEGER variables.
!   By default we use constant from the ISO C binding and use double precision for working.
!   If the host model in which RRTGMP is embedded has defined these constants elsewhere
!   the model definitions can be used instead by renaming. For example,
! use  mo_model_kind, only wp => dp, ...
!   where the syntax is local_name => original_name
!   and all the local names need to be defined

module mo_rte_kind
  use, intrinsic :: iso_c_binding, only: c_float, c_double, c_long, c_int, c_bool
  implicit none
  integer, parameter :: dp = c_double, sp = c_float, i8 = c_long, i4 = c_int
  !
  ! Floating point working precision
  !
  integer, parameter :: wp = dp

  !
  ! Logical - for use with kernels
  !
  integer, parameter :: wl = c_bool

end module mo_rte_kind
