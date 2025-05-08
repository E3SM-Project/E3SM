! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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

   function ocn_build_log_filename(prefix, identifier) result(filename)!{{{
      character (len=*), intent(in) :: prefix
      integer, intent(in) :: identifier

      character (len=StrKIND) :: filename

      character (len=StrKIND) :: identifierString

      if ( identifier .lt. 10 ) then
         write(identifierString, '(I1)') identifier
      else if ( identifier .lt. 100 ) then
         write(identifierString, '(I2)') identifier
      else if ( identifier .lt. 1000 ) then
         write(identifierString, '(I3)') identifier
      else if ( identifier .lt. 10000 ) then
         write(identifierString, '(I4)') identifier
      else if ( identifier .lt. 100000 ) then
         write(identifierString, '(I5)') identifier
      else if ( identifier .lt. 1000000 ) then
         write(identifierString, '(I6)') identifier
      else if ( identifier .lt. 10000000 ) then
         write(identifierString, '(I7)') identifier
      else
         write(identifierString, '(I99)') identifier
      end if

      filename = trim(prefix) // trim(identifierString)

   end function ocn_build_log_filename!}}}

!***********************************************************************

end module ocn_config

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
