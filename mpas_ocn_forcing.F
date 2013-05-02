! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_forcing
!
!> \brief MPAS ocean forcing
!> \author Doug Jacobsen
!> \date   04/25/12
!> \version SVN:$Id:$
!> \details
!>  This module contains driver routines for building the forcing arrays.
!
!-----------------------------------------------------------------------

module ocn_forcing

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_timekeeping
   use ocn_forcing_coupled

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_forcing_build_arrays, &
             ocn_forcing_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND) :: attenuationCoefficient

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_forcing_build_arrays
!
!> \brief   Determines the forcing arrays.
!> \author  Doug Jacobsen
!> \date    12/13/12
!> \version SVN:$Id$
!> \details 
!>  This routine computes the forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_build_arrays(timeStamp, mesh, forcing, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type(MPAS_Time_type), intent(in) :: timeStamp

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(inout) :: &
         mesh          !< Input: mesh information

      type (forcing_type), intent(inout) :: forcing !< Input: State information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      call ocn_forcing_coupled_build_arrays(timeStamp, mesh, forcing, err)

   !--------------------------------------------------------------------

   end subroutine ocn_forcing_build_arrays!}}}

!***********************************************************************
!
!  routine ocn_forcing_init
!
!> \brief   Initializes forcing module
!> \author  Doug Jacobsen
!> \date    12/13/12
!> \version SVN:$Id$
!> \details 
!>  This routine initializes the forcing modules.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      integer :: err1

      err = 0

      attenuationCoefficient = config_flux_attenuation_coefficient

      call ocn_forcing_coupled_init(err1)

      err = ior(err,err1)

   end subroutine ocn_forcing_init!}}}

!***********************************************************************

end module ocn_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
