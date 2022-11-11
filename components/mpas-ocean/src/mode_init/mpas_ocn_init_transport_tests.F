! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_transport_tests
!
!> \brief MPAS ocean initialize case -- tracer transport tests
!> \author Andrew Bradley and Peter Bosler (Sandia National Laboratories)
!> \date   12/20/2020
!> \details
!>  This module contains the routines for initializing the
!>  the 2D and 3D spherical advection test cases from DCMIP.
!
!-----------------------------------------------------------------------
module ocn_init_transport_tests
  use mpas_kind_types
  use mpas_io_units
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_stream_manager
  use mpas_dmpar

  use ocn_config
  use ocn_transport_tests

  implicit none
  private
  public :: ocn_init_setup_transport_tests, &
            ocn_init_validate_transport_tests

  contains

!***********************************************************************
!
!  routine ocn_init_setup_transport_tests
!
!> \brief   Setup for advection test cases
!> \author  Andrew Bradley
!> \date    12/20/2020
!> \details
!>  This routine sets up the initial conditions for the suite of time-reversible advection
!>  tests.
!
!-----------------------------------------------------------------------
  subroutine ocn_init_setup_transport_tests(domain, iErr)
    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    character (len=StrKIND), pointer :: config_init_configuration

#ifdef MPAS_DEBUG
    call mpas_log_write("entered subroutine ocn_init_setup_transport_tests")
#endif

    call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)

    iErr = 0
    if (config_init_configuration /= trim('transport_tests')) return
    call ocn_transport_test_init_setup(domain, iErr)
#ifdef MPAS_DEBUG
    call mpas_log_write("exiting subroutine ocn_init_setup_transport_tests")
#endif
  end subroutine ocn_init_setup_transport_tests

!***********************************************************************
!
!  routine ocn_init_validate_transport_tests
!
!> \brief   Validation for advection test cases
!> \author  Andrew Bradley
!> \date    12/20/2020
!> \details
!>  This routine validates the configuration options for the time-reversible
!>  advection test cases.
!
!-----------------------------------------------------------------------
  subroutine ocn_init_validate_transport_tests(configPool, packagePool, iocontext, iErr)
    type (mpas_pool_type), intent(inout) :: configPool
    type (mpas_pool_type), intent(inout) :: packagePool
    type (mpas_io_context_type), intent(inout) :: iocontext
    integer, intent(out) :: iErr

    character (len=StrKIND), pointer :: config_init_configuration

#ifdef MPAS_DEBUG
    call mpas_log_write("entered subroutine ocn_init_validate_transport_tests")
#endif

    call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)
    iErr = 0
    if (trim(config_init_configuration) /= trim('transport_tests')) return


    call ocn_transport_test_init_validate(configPool, packagePool, iocontext, iErr)
  end subroutine ocn_init_validate_transport_tests

end module ocn_init_transport_tests
