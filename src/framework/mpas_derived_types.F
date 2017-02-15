! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!***********************************************************************
!
!  mpas_derived_types
!
!> \brief   MPAS Derived data types
!> \author  Doug Jacobsen, Michael Duda, Matt Hoffman
!> \date    03/10/2015
!> \details
!> This module defines derived data types related to fields, and variable structures.
!> Routines for allocating and deallocating these types are located in their
!> corresponding modules.  All of the framework derived types are collected in this
!> module so that they can all be "use"d from this module to avoid circular
!> dependencies.
!> This module also includes public parameters related to the corresponding modules.
!> Including them here provides easy access to them throughout the code without
!> needing to use the module they "belong" to.
!
!-----------------------------------------------------------------------
module mpas_derived_types

   use mpas_kind_types

   use pio
   use pio_types

   use ESMF
   use ESMF_BaseMod
   use ESMF_Stubs
   use ESMF_CalendarMod
   use ESMF_ClockMod
   use ESMF_TimeMod
   use ESMF_TimeIntervalMod

#include "mpas_attlist_types.inc"

#include "mpas_hash_types.inc"

#include "mpas_dmpar_types.inc"

#include "mpas_field_types.inc"

#include "mpas_pool_types.inc"

#include "mpas_particle_list_types.inc"

#include "mpas_io_types.inc"

#include "mpas_io_streams_types.inc"

#include "mpas_stream_list_types.inc"

#include "mpas_stream_manager_types.inc"

#include "mpas_timekeeping_types.inc"

#include "mpas_timer_types.inc"

#include "mpas_block_types.inc"

#include "mpas_decomp_types.inc"

#include "mpas_domain_types.inc"

#include "mpas_core_types.inc"

#include "mpas_forcing_types.inc"

   contains

end module mpas_derived_types
