! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_config
!
!> \brief MALI specific config
!> \details
!>  This module contains config specific to MALI.
!
!-----------------------------------------------------------------------

module li_config

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_kind_types

   implicit none
   public
   save

#include "../inc/config_declare.inc"

!***********************************************************************

contains

!***********************************************************************
!
!  routine li_config_init
!
!> \brief   Initializes the MALI config
!> \details
!>  This routine sets up config for use in MALI.
!
!-----------------------------------------------------------------------
   subroutine li_config_init(configPool)!{{{
       type (mpas_pool_type), pointer :: configPool

#include "../inc/config_get.inc"

   end subroutine li_config_init!}}}

!***********************************************************************

end module li_config

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
