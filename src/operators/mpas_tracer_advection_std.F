! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_tracer_advection_std
!
!> \brief MPAS tracer advection driver (non-monotonic/fct)
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains driver routine for tracer advection tendencies
!>  as well as the routines for setting up advection coefficients and 
!>  initialization of the advection routines.
!
!-----------------------------------------------------------------------
module mpas_tracer_advection_std

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar
   use mpas_timer

   use mpas_tracer_advection_std_hadv
   use mpas_tracer_advection_std_vadv

   implicit none
   private
   save

   public :: mpas_tracer_advection_std_tend, &
             mpas_tracer_advection_std_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_tracer_advection_std_tend
!
!> \brief MPAS standard tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine is the driver routine for the standard computation of 
!>  tracer advection tendencies.
!
!-----------------------------------------------------------------------
   subroutine mpas_tracer_advection_std_tend(tracers, normalThicknessFlux, w, verticalCellSize, grid, tend)!{{{

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: Tracer tendency
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: Tracer values
      real (kind=RKIND), dimension(:,:), intent(in) :: normalThicknessFlux !< Input: Thickness weighted horizontal velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: w !< Input: Vertical Velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of a cell
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      call mpas_timer_start("tracer-hadv", .false.)
      call mpas_tracer_advection_std_hadv_tend(tracers, normalThicknessFlux, grid, tend)
      call mpas_timer_stop("tracer-hadv")
      call mpas_timer_start("tracer-vadv", .false.)
      call mpas_tracer_advection_std_vadv_tend(tracers, w, verticalCellSize, grid, tend)
      call mpas_timer_stop("tracer-vadv")

   end subroutine mpas_tracer_advection_std_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_tracer_advection_std_init
!
!> \brief MPAS standard tracer advection initialization
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine is the driver routine for the initializtion of the standard 
!>  tracer advection routines.
!
!-----------------------------------------------------------------------
   subroutine mpas_tracer_advection_std_init(err)!{{{
      integer, intent(inout) :: err !< Input: Error Flags

      integer :: err_tmp

      err = 0

      call mpas_tracer_advection_std_hadv_init(err_tmp)
      err = ior(err, err_tmp)
      call mpas_tracer_advection_std_vadv_init(err_tmp)
      err = ior(err, err_tmp)

   end subroutine mpas_tracer_advection_std_init!}}}

end module mpas_tracer_advection_std
