! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_config
!
!> \brief MPAS ocean specific config
!> \details
!>  This module contains config specific to the ocean model.
!
!-----------------------------------------------------------------------

module ocn_config

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
!  routine ocn_config_init
!
!> \brief   Initializes the ocean config
!> \details
!>  This routine sets up config for use in the ocean model.
!
!-----------------------------------------------------------------------
   subroutine ocn_config_init(configPool)!{{{
       type (mpas_pool_type), pointer :: configPool

#include "../inc/config_get.inc"

   end subroutine ocn_config_init!}}}

!***********************************************************************

end module ocn_config

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
