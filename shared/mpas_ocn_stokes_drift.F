! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_stokes_drift
!
!> \brief  Routines related to calulationing Stokes drift profiles
!> \author Steven Brus 
!> \date   June 2021
!> \details
!>  This module contains the routines for calculating wave-induced 
!>  Stokes drift profiles from an active wave model.
!
!-----------------------------------------------------------------------

module ocn_stokes_drift

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_constants
   use ocn_config
   use mpas_log

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

   public :: ocn_stokes_drift_reconstruct, &
             ocn_stokes_drift_average, &
             ocn_stokes_drift_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_stokes_drift_reconstruct
!
!> \brief   Reconstructs depth-dependent Stokes drift profile 
!> \author  Steven Brus
!> \date    June 2021 
!> \details
!>   This routine reconstructs the depth-dependent Stokes drift profile 
!>   using partitioned surface Stokes drift in 6 wavenumber bins from 
!>   an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_stokes_drift_reconstruct()!{{{
      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------


      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_reconstruct!}}}

!***********************************************************************
!
!  routine ocn_stokes_drift_average
!
!> \brief   Averages Stokes drift profile over a specified depth
!> \author  Steven Brus
!> \date    June 2021 
!> \details
!>   This routine averages the Stokes drift profile over a specified 
!>   depth using partitioned surface Stokes drift in 6 wavenumber bins
!>   from an active wave model. 
!
!-----------------------------------------------------------------------

   subroutine ocn_stokes_drift_reconstruct()!{{{
      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------


      !-----------------------------------------------------------------

   end subroutine ocn_stokes_drift_reconstruct!}}}

!***********************************************************************
!
!  routine ocn_stokes_drift_init
!
!> \brief   Initializes Stokes drift computations 
!> \author  Steven Brus
!> \date    June 2021
!> \details
!>  This routine initializes quantities related to Stokes drift-related
!>  computations
!
!----------------------------------------------------------------------

   subroutine ocn_stokes_drift_init(domain, err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: &
         domain            ! domain information primarily for vert grid

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_stokes_drift_init!}}}



!***********************************************************************

end module ocn_stokes_drift

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
