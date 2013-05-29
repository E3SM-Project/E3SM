! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_advection_helpers
!
!> \brief MPAS ocean tracer advection helper functions
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This module contains helper functions tracer advection.
!
!-----------------------------------------------------------------------
module mpas_ocn_tracer_advection_helpers

   use mpas_kind_types

   implicit none
   save

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  function mpas_ocn_tracer_advection_vflux4
!
!> \brief MPAS ocean 4th order vertical tracer advection stencil
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This function provides the stencil for 4th order vertical advection of tracers.
!
!-----------------------------------------------------------------------
   real function mpas_ocn_tracer_advection_vflux4(q_im2, q_im1, q_i, q_ip1, w)!{{{
        real (kind=RKIND), intent(in) :: q_im2 !< Input: Tracer value at index i-2
        real (kind=RKIND), intent(in) :: q_im1 !< Input: Tracer value at index i-1
        real (kind=RKIND), intent(in) :: q_i !< Input: Tracer value at index i
        real (kind=RKIND), intent(in) :: q_ip1 !< Input: Tracer value at index i+1
        real (kind=RKIND), intent(in) :: w !< Input: vertical veloicity
        mpas_ocn_tracer_advection_vflux4 = w*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0
   end function!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  function mpas_ocn_tracer_advection_vflux3
!
!> \brief MPAS ocean 3rd order vertical tracer advection stencil
!> \author Doug Jacobsen
!> \date   03/09/12
!> \details
!>  This function provides the stencil for 3rd order vertical advection of tracers.
!
!-----------------------------------------------------------------------
   real function mpas_ocn_tracer_advection_vflux3( q_im2, q_im1, q_i, q_ip1, w, coef)!{{{
        real (kind=RKIND), intent(in) :: q_im2 !< Input: Tracer value at index i-2
        real (kind=RKIND), intent(in) :: q_im1 !< Input: Tracer value at index i-1
        real (kind=RKIND), intent(in) :: q_i !< Input: Tracer value at index i
        real (kind=RKIND), intent(in) :: q_ip1 !< Input: Tracer value at index i+1
        real (kind=RKIND), intent(in) :: w !< Input: vertical veloicity
        real (kind=RKIND), intent(in) :: coef !< Input: Advection coefficient

        !dwj 02/21/12 flux3 is different in ocean and atmosphere
        !flux3 = (u * (7.0 * (q_i + q_im1) - (q_ip1 + q_im2)) + coef * abs(u) * ((q_ip1 - q_im2) - 3.0*(q_i-q_im1)))/12.0
        mpas_ocn_tracer_advection_vflux3 = (w * (7.0 * (q_i + q_im1) - (q_ip1 + q_im2)) - coef * abs(w) * ((q_ip1 - q_im2) - 3.0*(q_i-q_im1)))/12.0
   end function!}}}

end module mpas_ocn_tracer_advection_helpers
