! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_two_equation_structure_functions
!
!> \brief MPAS ocean structure functions calculation 
!> \author Luke Van Roekel
!> \date   July 2021
!> \details
!>  This module contains the computation for structure functions  
!>  for the two equation turbulence model viscosity/diffusivity calculation
!>
!
!-----------------------------------------------------------------------

module ocn_two_equation_structure_functions

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading
   use mpas_constants

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

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

   public :: ocn_canuto_structure_function_compute

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_structure_function_compute
!
!> \brief   Computes structure functions for diffusivity and viscosity  
!> \author  Luke Van Roekel 
!> \date    July 2021
!> \details
!>  This routine computes the structure functions for the viscosity and
!>  diffusivity in the two equation model, follows Canuto et al (2001)
!>  Canuto, V. M., Howard, A., Cheng, Y., & Dubovikov, M. S. (2001). 
!>     Ocean Turbulence. Part I: One-Point Closure Model—Momentum and 
!>     Heat Vertical Diffusivities, Journal of Physical Oceanography, 
!>     31(6), 1413-1426
!
!-----------------------------------------------------------------------

   subroutine ocn_canuto_structure_function_compute(Sm, Sh, tke, glsPsi, c_p, &
                  c_m, c_n, c_mu)!{{{
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
                        Sm, Sh ! structure functions

      real (kind=RKIND), dimension(:,:), intent(in) :: &
                        tke, glsPsi ! input variables for calculation

      real (kind=RKIND), intent(in) :: c_p, c_m, c_n, c_mu

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      real(kind=RKIND), parameter :: lambda1 = 0.107_RKIND,  &
                                     lambda2 = 0.0032_RKIND, &
                                     lambda3 = 0.0864_RKIND, &
                                     lambda4 = 0.12_RKIND,   &
                                     lambda5 = 11.9_RKIND,   &
                                     lambda6 = 0.4_RKIND,    &
                                     lambda7 = 0.0_RKIND,    &
                                     lambda8 = 0.48_RKIND

      real(kind=RKIND) :: d0, d1, d2, d3, d4, d5, riMid, bvfMid

      real(kind=RKIND) :: tau, D, dissipation, s0, s1, s2, s3, s4, s5, s6

      integer :: k, iCell

      d0 = 3.0_RKIND*lambda5**2.0
      d1 = lambda5*(7.0_RKIND*lambda4 + 3.0_RKIND*lambda8)
      d2 = lambda5**2.0*(3.0_RKIND*lambda3**2.0 - lambda2**2.0) -  &
           0.75_RKIND*(lambda6**2.0 - lambda7**2.0)
      d3 = lambda4*(4.0_RKIND*lambda4 + 3.0_RKIND*lambda8)
      d4 = lambda4*(lambda2*lambda6 - 3.0_RKIND*lambda3*lambda7 - &
           lambda5*(lambda2**2.0 - lambda3**2.0)) + lambda5*      &
           lambda8*(3.0_RKIND*lambda3**2.0 - lambda2**2.0)
      d5 = 0.25_RKIND*(lambda2**2.0 - 3.0_RKIND*lambda3**2.0)*    &
           (lambda6**2.0 - lambda7**2.0)

      s0 = 1.5_RKIND*lambda1*lambda5**2.0
      s1 = -lambda4*(lambda6 + lambda7) + 2.0_RKIND*lambda4*lambda5* &
              (lambda1 - 1.0_RKIND/3.0_RKIND*lambda2 - lambda3) +    &
              1.5_RKIND*lambda1*lambda5*lambda8
      s2 = -3.0_RKIND/8.0_RKIND*lambda1*(lambda6**2.0 - lambda7**2.0)
      s4 = 2.0_RKIND*lambda5
      s5 = 2.0_RKIND*lambda4
      s6 = 2.0_RKIND/3.0_RKIND*lambda5*(3.0_RKIND*lambda3**2.0 - lambda2**2.0) - &
            0.5_RKIND*lambda5*lambda1*(3.0_RKIND*lambda3 - lambda2) + &
            3.0_RKIND/4.0_RKIND*lambda1*(lambda6 - lambda7)

      do iCell=1,nCellsOwned
         do k=1,maxLevelCell(iCell)-1
            dissipation = c_mu**(3.0+c_p/c_n)*tke(k,iCell)**(1.5+c_m/c_n)* &
                                        glsPsi(k,iCell)**(-1.0/c_n)
            tau = 2.0_RKIND*tke(k,iCell)/(dissipation + 1.0E-15_RKIND)
            bvfMid = 0.5_RKIND*(bruntVaisalaFreqTop(k,iCell) +              &
                                bruntVaisalaFreqTop(k+1,iCell))
            riMid = 0.5_RKIND*(RiTopOfCell(k,iCell) + RiTopOfCell(k+1,iCell))
            D = d0 + d1*tau**2.0*bvfMid + d2*tau**2.0* &
                bvfMid / (1.0E-15_RKIND + riMid) &
                + d3*tau**4.0*bvfMid**2.0 + d4*tau**4.0* &
                bvfMid**2.0 / (1.0E-15_RKIND +           &
                riMid) + d5*tau**4.0*bvfMid**2.0 / &
                (1.0E-15_RKIND + riMid**2.0)
            Sm(k,iCell) = (s0 + s1*tau**2.0*bvfMid + &
                          s2*tau**2.0*bvfMid /       &
                          (1.0E-15_RKIND + riMid)) / (D + 1.0E-15_RKIND)
            Sh(k,iCell) = (s4 + s5*tau**2.0*bvfMid + &
                          s6*tau**2.0*bvfMid /       &
                          (1.0E-15_RKIND + riMid)) / (D + 1.0E-15_RKIND)
        end do
        Sm(maxLevelCell(iCell):nVertLevels,iCell) = 0.0_RKIND
        Sh(maxLevelCell(iCell):nVertLevels,iCell) = 0.0_RKIND
     end do

   end subroutine ocn_canuto_structure_function_compute

end module ocn_two_equation_structure_functions

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
