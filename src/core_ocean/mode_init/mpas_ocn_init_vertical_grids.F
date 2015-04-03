! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_vertical_grids
!
!> \brief MPAS ocean vertical grid generator
!> \author Doug Jacobsen
!> \date   03/20/2015
!> \details
!>  This module contains the routines for generating
!>  vertical grids.
!
!-----------------------------------------------------------------------
module ocn_init_vertical_grids

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timer

   use ocn_constants

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

   public :: ocn_generate_vertical_grid

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

   !***********************************************************************
   !
   !  routine ocn_generate_vertical_grid
   !
   !> \brief   Vertical grid generator driver
   !> \author  Doug Jacobsen
   !> \date    03/20/2015
   !> \details 
   !>  This routine is a driver for generating vertical grids. It calls a private
   !>  module routine based on the value of the input argument gridType.
   !>  The output array layerInterfaces will contain values between 1 and 0
   !>  representing the relative locations of layer interfaces.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_generate_vertical_grid(gridType, layerInterfaces)!{{{
      implicit none

      character (len=*), intent(in) :: gridType
      real (kind=RKIND), dimension(:), intent(out) :: layerInterfaces

      if ( trim(gridType) == 'uniform' ) then
         call ocn_generate_uniform_vertical_grid(layerInterfaces)
      else if ( trim(gridType) == '60layerPHC' ) then
         call ocn_generate_60layerPHC_vertical_grid(layerInterfaces)
      else if ( trim(gridType) == '42layerWOCE' ) then
         call ocn_generate_42layerWOCE_vertical_grid(layerInterfaces)
      else
         write(stderrUnit, *) ' WARNING: '//trim(gridType)//' is an invalid vertical grid choice. No vertical grid will be generated.'
      end if

   end subroutine ocn_generate_vertical_grid!}}}

   !***********************************************************************
   !
   !  routine ocn_generate_uniform_vertical_grid
   !
   !> \brief   Uniform Vertical grid generator
   !> \author  Doug Jacobsen
   !> \date    03/20/2015
   !> \details 
   !>  This routine generates a uniform vertical grid.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_generate_uniform_vertical_grid(layerInterfaces)!{{{
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: layerInterfaces

      real (kind=RKIND) :: layerSpacing
      integer :: nInterfaces, iInterface

      write(stderrUnit,* ) ' ---- Generating uniform vertical grid ---- '

      nInterfaces = size(layerInterfaces, dim=1)
      layerSpacing = 1.0_RKIND / (nInterfaces - 1)

      layerInterfaces(1) = 0.0_RKIND

      do iInterface = 2, nInterfaces
         layerInterfaces(iInterface) = layerInterfaces(iInterface-1) + layerSpacing
      end do

   end subroutine ocn_generate_uniform_vertical_grid!}}}

   !***********************************************************************
   !
   !  routine ocn_generate_60layerPHC_vertical_grid
   !
   !> \brief   60 layer PHC vertical grid generator
   !> \author  Doug Jacobsen
   !> \date    03/20/2015
   !> \details 
   !>  This routine generates a 60 layer vertical grid based on the PHC data set.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_generate_60layerPHC_vertical_grid(layerInterfaces)!{{{
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: layerInterfaces

      real (kind=RKIND) :: maxInterfaceLocation
      integer :: nInterfaces, iInterface

      nInterfaces = size(layerInterfaces, dim=1)

      if ( nInterfaces /= 61 ) then
         call mpas_dmpar_global_abort("ERROR: Vertical grid must have 60 layers to apply 60 Layer PHC grid. Exiting...")
      end if

      layerInterfaces(1) = 0.0_RKIND
      layerInterfaces(2) = 500_RKIND
      layerInterfaces(3) = 1500_RKIND
      layerInterfaces(4) = 2500_RKIND
      layerInterfaces(5) = 3500_RKIND
      layerInterfaces(6) = 4500_RKIND
      layerInterfaces(7) = 5500_RKIND
      layerInterfaces(8) = 6500_RKIND
      layerInterfaces(9) = 7500_RKIND
      layerInterfaces(10) = 8500_RKIND
      layerInterfaces(11) = 9500_RKIND
      layerInterfaces(12) = 10500_RKIND
      layerInterfaces(13) = 11500_RKIND
      layerInterfaces(14) = 12500_RKIND
      layerInterfaces(15) = 13500_RKIND
      layerInterfaces(16) = 14500_RKIND
      layerInterfaces(17) = 15500_RKIND
      layerInterfaces(18) = 16509.83984375_RKIND
      layerInterfaces(19) = 17547.904296875_RKIND
      layerInterfaces(20) = 18629.125_RKIND
      layerInterfaces(21) = 19766.025390625_RKIND
      layerInterfaces(22) = 20971.134765625_RKIND
      layerInterfaces(23) = 22257.826171875_RKIND
      layerInterfaces(24) = 23640.880859375_RKIND
      layerInterfaces(25) = 25137.013671875_RKIND
      layerInterfaces(26) = 26765.416015625_RKIND
      layerInterfaces(27) = 28548.361328125_RKIND
      layerInterfaces(28) = 30511.91796875_RKIND
      layerInterfaces(29) = 32686.794921875_RKIND
      layerInterfaces(30) = 35109.34375_RKIND
      layerInterfaces(31) = 37822.75390625_RKIND
      layerInterfaces(32) = 40878.4609375_RKIND
      layerInterfaces(33) = 44337.765625_RKIND
      layerInterfaces(34) = 48273.66796875_RKIND
      layerInterfaces(35) = 52772.796875_RKIND
      layerInterfaces(36) = 57937.28515625_RKIND
      layerInterfaces(37) = 63886.2578125_RKIND
      layerInterfaces(38) = 70756.328125_RKIND
      layerInterfaces(39) = 78700.25_RKIND
      layerInterfaces(40) = 87882.5234375_RKIND
      layerInterfaces(41) = 98470.5859375_RKIND
      layerInterfaces(42) = 110620.421875_RKIND
      layerInterfaces(43) = 124456.6953125_RKIND
      layerInterfaces(44) = 140049.71875_RKIND
      layerInterfaces(45) = 157394.640625_RKIND
      layerInterfaces(46) = 176400.328125_RKIND
      layerInterfaces(47) = 196894.421875_RKIND
      layerInterfaces(48) = 218645.65625_RKIND
      layerInterfaces(49) = 241397.15625_RKIND
      layerInterfaces(50) = 264900.125_RKIND
      layerInterfaces(51) = 288938.46875_RKIND
      layerInterfaces(52) = 313340.46875_RKIND
      layerInterfaces(53) = 337979.375_RKIND
      layerInterfaces(54) = 362767.0625_RKIND
      layerInterfaces(55) = 387645.21875_RKIND
      layerInterfaces(56) = 412576.84375_RKIND
      layerInterfaces(57) = 437539.28125_RKIND
      layerInterfaces(58) = 462519.0625_RKIND
      layerInterfaces(59) = 487508.375_RKIND
      layerInterfaces(60) = 512502.84375_RKIND
      layerInterfaces(61) = 537500_RKIND

      maxInterfaceLocation = maxval(layerInterfaces)

      layerInterfaces(:) = layerInterfaces(:) / maxInterfaceLocation

   end subroutine ocn_generate_60layerPHC_vertical_grid!}}}

   !***********************************************************************
   !
   !  routine ocn_generate_42layerWOCE_vertical_grid
   !
   !> \brief   42 layer WOCE vertical grid generator
   !> \author  Doug Jacobsen
   !> \date    03/20/2015
   !> \details 
   !>  This routine generates a 42 layer vertical grid based on the WOCE data set.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_generate_42layerWOCE_vertical_grid(layerInterfaces)!{{{
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: layerInterfaces

      real (kind=RKIND) :: maxInterfaceLocation
      integer :: nInterfaces, iInterface

      nInterfaces = size(layerInterfaces, dim=1)

      if ( nInterfaces /= 43 ) then
         call mpas_dmpar_global_abort("ERROR: Vertical grid must have 60 layers to apply 60 Layer PHC grid. Exiting...")
      end if

      layerInterfaces(1) = 0.0_RKIND
      layerInterfaces(2) = 5.00622_RKIND
      layerInterfaces(3) = 15.06873_RKIND
      layerInterfaces(4) = 25.28343_RKIND
      layerInterfaces(5) = 35.75849_RKIND
      layerInterfaces(6) = 46.61269_RKIND
      layerInterfaces(7) = 57.98099_RKIND
      layerInterfaces(8) = 70.02139_RKIND
      layerInterfaces(9) = 82.92409_RKIND
      layerInterfaces(10) = 96.92413_RKIND
      layerInterfaces(11) = 112.3189_RKIND
      layerInterfaces(12) = 129.4936_RKIND
      layerInterfaces(13) = 148.9582_RKIND
      layerInterfaces(14) = 171.4044_RKIND
      layerInterfaces(15) = 197.7919_RKIND
      layerInterfaces(16) = 229.4842_RKIND
      layerInterfaces(17) = 268.4617_RKIND
      layerInterfaces(18) = 317.6501_RKIND
      layerInterfaces(19) = 381.3864_RKIND
      layerInterfaces(20) = 465.9132_RKIND
      layerInterfaces(21) = 579.3073_RKIND
      layerInterfaces(22) = 729.3513_RKIND
      layerInterfaces(23) = 918.3723_RKIND
      layerInterfaces(24) = 1139.153_RKIND
      layerInterfaces(25) = 1378.574_RKIND
      layerInterfaces(26) = 1625.7_RKIND
      layerInterfaces(27) = 1875.106_RKIND
      layerInterfaces(28) = 2125.011_RKIND
      layerInterfaces(29) = 2375_RKIND
      layerInterfaces(30) = 2624.999_RKIND
      layerInterfaces(31) = 2874.999_RKIND
      layerInterfaces(32) = 3124.999_RKIND
      layerInterfaces(33) = 3374.999_RKIND
      layerInterfaces(34) = 3624.999_RKIND
      layerInterfaces(35) = 3874.999_RKIND
      layerInterfaces(36) = 4124.999_RKIND
      layerInterfaces(37) = 4374.999_RKIND
      layerInterfaces(38) = 4624.999_RKIND
      layerInterfaces(39) = 4874.999_RKIND
      layerInterfaces(40) = 5124.999_RKIND
      layerInterfaces(41) = 5374.999_RKIND
      layerInterfaces(42) = 5624.999_RKIND
      layerInterfaces(43) = 5874.999_RKIND

      maxInterfaceLocation = maxval(layerInterfaces)

      layerInterfaces(:) = layerInterfaces(:) / maxInterfaceLocation

   end subroutine ocn_generate_42layerWOCE_vertical_grid!}}}


!***********************************************************************

end module ocn_init_vertical_grids

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
