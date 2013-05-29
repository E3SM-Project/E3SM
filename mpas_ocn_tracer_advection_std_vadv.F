! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_std_vadv
!
!> \brief MPAS ocean vertical tracer advection driver (non-monotonic/fct)
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains driver routines for vertical tracer advection tendencies
!>  and initialization of the advection routines.
!
!-----------------------------------------------------------------------
module mpas_ocn_tracer_advection_std_vadv

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar

   use mpas_ocn_tracer_advection_std_vadv2
   use mpas_ocn_tracer_advection_std_vadv3
   use mpas_ocn_tracer_advection_std_vadv4

   implicit none
   private
   save

   public :: mpas_ocn_tracer_advection_std_vadv_tend, &
             mpas_ocn_tracer_advection_std_vadv_init

   logical :: order2, order3, order4

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_std_vadv_tend
!
!> \brief MPAS ocean standard vertical tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine is the driver routine for the standard computation of 
!>  vertical tracer advection tendencies.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_std_vadv_tend(tracers, w, verticalCellSize, grid, tend)!{{{

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: tend !< Input/Output: Tracer Tendency
      real (kind=RKIND), dimension(:,:,:), intent(in) :: tracers !< Input: Tracer Values
      real (kind=RKIND), dimension(:,:), intent(in) :: w !< Input: Vertical Velocity
      real (kind=RKIND), dimension(:,:), intent(in) :: verticalCellSize !< Input: Distance between vertical interfaces of cell
      type (mesh_type), intent(in) :: grid !< Input: Grid information

      if(order2) then
        call mpas_ocn_tracer_advection_std_vadv2_tend(tracers, w, verticalCellSize, grid, tend) 
      else if(order3) then
        call mpas_ocn_tracer_advection_std_vadv3_tend(tracers, w, verticalCellSize, grid, tend) 
      else if(order4) then
        call mpas_ocn_tracer_advection_std_vadv4_tend(tracers, w, verticalCellSize, grid, tend) 
      endif

   end subroutine mpas_ocn_tracer_advection_std_vadv_tend!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine mpas_ocn_tracer_advection_std_vadv_init
!
!> \brief MPAS ocean standard vertical tracer advection tendency
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This routine initializes the vertical tracer advection tendencies.
!
!-----------------------------------------------------------------------
   subroutine mpas_ocn_tracer_advection_std_vadv_init(err)!{{{
     integer, intent(inout) :: err !< Input/Output: Error flag

     err = 0

     order2 = .false.
     order3 = .false.
     order4 = .false.

     if (config_vert_tracer_adv_order == 2) then
         order2 = .true.
     else if (config_vert_tracer_adv_order == 3) then
         order3 = .true.
     else if (config_vert_tracer_adv_order == 4) then
         order4 = .true.
     else
         print *, 'invalid value for config_tracer_vadv_order'
         print *, 'options are 2, 3, or 4'
         err = 1
     endif
        
   end subroutine mpas_ocn_tracer_advection_std_vadv_init!}}}

end module mpas_ocn_tracer_advection_std_vadv

