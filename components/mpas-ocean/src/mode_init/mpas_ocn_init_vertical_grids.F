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
   use mpas_io

   use ocn_constants
   use ocn_config
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

   public :: ocn_generate_vertical_grid, &
             ocn_compute_z_level_layerThickness, &
             ocn_compute_zMid_from_layerThickness, &
             ocn_compute_z_star_layerThickness, &
             ocn_alter_bottomDepth_for_pbcs, &
             ocn_alter_ssh_for_ptcs, &
             ocn_compute_Haney_number, &
             ocn_init_vertical_grid, &
             ocn_init_vertical_grid_with_max_rx1

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
   !>  The output array interfaceLocations will contain values between
   !>  0 being the top of top layer and 1 being the bottom of bottom layer
   !
   !-----------------------------------------------------------------------
   subroutine ocn_generate_vertical_grid(gridType, interfaceLocations, configPool)!{{{
      implicit none

      character (len=*), intent(in) :: gridType
      real (kind=RKIND), dimension(:), intent(out) :: interfaceLocations
      type (mpas_pool_type), optional, intent(in) :: configPool !< Input: Pool with namelist options

      if ( trim(gridType) == 'uniform' ) then
         call ocn_generate_uniform_vertical_grid(interfaceLocations)
      else if ( trim(gridType) == '60layerPHC' ) then
         call ocn_generate_60layerPHC_vertical_grid(interfaceLocations)
      else if ( trim(gridType) == '42layerWOCE' ) then
         call ocn_generate_42layerWOCE_vertical_grid(interfaceLocations)
      else if ( trim(gridType) == '100layerE3SMv1' ) then
         call ocn_generate_100layerE3SMv1_vertical_grid(interfaceLocations)
      else if ( trim(gridType) == '1dCVTgenerator' ) then
         call ocn_generate_1dCVT_vertical_grid(interfaceLocations)
      else
         call mpas_log_write( ' WARNING: '//trim(gridType)//' is an invalid vertical grid choice. No vertical ' &
                           // 'grid will be generated.')
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
   subroutine ocn_generate_uniform_vertical_grid(interfaceLocations)!{{{
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: interfaceLocations

      real (kind=RKIND) :: layerSpacing
      integer :: nInterfaces, iInterface

      call mpas_log_write( ' ---- Generating uniform vertical grid ---- ')

      nInterfaces = size(interfaceLocations, dim=1)
      layerSpacing = 1.0_RKIND / (nInterfaces - 1)

      interfaceLocations(1) = 0.0_RKIND

      do iInterface = 2, nInterfaces
         interfaceLocations(iInterface) = interfaceLocations(iInterface-1) + layerSpacing
      end do

   end subroutine ocn_generate_uniform_vertical_grid!}}}

   !***********************************************************************
   !
   !  routine ocn_generate_60layerPHC_vertical_grid
   !
   !> \brief    60 layer PHC vertical grid generator
   !> \author   Doug Jacobsen, Xylar Asay-Davis
   !> \date     03/20/2015
   !> \modified 01/11/2022
   !> \details
   !>  This routine generates a 60 layer vertical grid based on the PHC data set.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_generate_60layerPHC_vertical_grid(interfaceLocations)!{{{
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: interfaceLocations

      real (kind=RKIND) :: maxInterfaceLocation
      integer :: nInterfaces, iInterface

      nInterfaces = size(interfaceLocations, dim=1)

      if ( nInterfaces /= 61 ) then
         call mpas_log_write("MPAS-ocean: Vertical grid must have 60 layers to "// &
                                      "apply 60 Layer PHC grid. Exiting...", MPAS_LOG_CRIT)
      end if

      interfaceLocations(1) = 0.0_RKIND
      interfaceLocations(2) = 10.0_RKIND
      interfaceLocations(3) = 20.0_RKIND
      interfaceLocations(4) = 30.0_RKIND
      interfaceLocations(5) = 40.0_RKIND
      interfaceLocations(6) = 50.0_RKIND
      interfaceLocations(7) = 60.0_RKIND
      interfaceLocations(8) = 70.0_RKIND
      interfaceLocations(9) = 80.0_RKIND
      interfaceLocations(10) = 90.0_RKIND
      interfaceLocations(11) = 100.0_RKIND
      interfaceLocations(12) = 110.0_RKIND
      interfaceLocations(13) = 120.0_RKIND
      interfaceLocations(14) = 130.0_RKIND
      interfaceLocations(15) = 140.0_RKIND
      interfaceLocations(16) = 150.0_RKIND
      interfaceLocations(17) = 160.0_RKIND
      interfaceLocations(18) = 170.19677734375_RKIND
      interfaceLocations(19) = 180.76129150390625_RKIND
      interfaceLocations(20) = 191.82119750976562_RKIND
      interfaceLocations(21) = 203.49929809570312_RKIND
      interfaceLocations(22) = 215.92340087890625_RKIND
      interfaceLocations(23) = 229.23312377929688_RKIND
      interfaceLocations(24) = 243.58447265625_RKIND
      interfaceLocations(25) = 259.1557922363281_RKIND
      interfaceLocations(26) = 276.1524963378906_RKIND
      interfaceLocations(27) = 294.8147277832031_RKIND
      interfaceLocations(28) = 315.4236145019531_RKIND
      interfaceLocations(29) = 338.3122863769531_RKIND
      interfaceLocations(30) = 363.8746032714844_RKIND
      interfaceLocations(31) = 392.5804748535156_RKIND
      interfaceLocations(32) = 424.9887390136719_RKIND
      interfaceLocations(33) = 461.7665710449219_RKIND
      interfaceLocations(34) = 503.7067565917969_RKIND
      interfaceLocations(35) = 551.7491760253906_RKIND
      interfaceLocations(36) = 606.9965515136719_RKIND
      interfaceLocations(37) = 670.7285461425781_RKIND
      interfaceLocations(38) = 744.3980407714844_RKIND
      interfaceLocations(39) = 829.6069641113281_RKIND
      interfaceLocations(40) = 928.0434265136719_RKIND
      interfaceLocations(41) = 1041.3681945800781_RKIND
      interfaceLocations(42) = 1171.0402526855469_RKIND
      interfaceLocations(43) = 1318.0935363769531_RKIND
      interfaceLocations(44) = 1482.9008483886719_RKIND
      interfaceLocations(45) = 1664.9919738769531_RKIND
      interfaceLocations(46) = 1863.0146179199219_RKIND
      interfaceLocations(47) = 2074.873809814453_RKIND
      interfaceLocations(48) = 2298.039276123047_RKIND
      interfaceLocations(49) = 2529.903594970703_RKIND
      interfaceLocations(50) = 2768.098846435547_RKIND
      interfaceLocations(51) = 3010.670196533203_RKIND
      interfaceLocations(52) = 3256.138885498047_RKIND
      interfaceLocations(53) = 3503.448516845703_RKIND
      interfaceLocations(54) = 3751.892791748047_RKIND
      interfaceLocations(55) = 4001.011505126953_RKIND
      interfaceLocations(56) = 4250.525604248047_RKIND
      interfaceLocations(57) = 4500.259552001953_RKIND
      interfaceLocations(58) = 4750.121307373047_RKIND
      interfaceLocations(59) = 5000.045684814453_RKIND
      interfaceLocations(60) = 5250.010955810547_RKIND
      interfaceLocations(61) = 5499.989044189453_RKIND

      maxInterfaceLocation = maxval(interfaceLocations)

      interfaceLocations(:) = interfaceLocations(:) / maxInterfaceLocation

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
   subroutine ocn_generate_42layerWOCE_vertical_grid(interfaceLocations)!{{{
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: interfaceLocations

      real (kind=RKIND) :: maxInterfaceLocation
      integer :: nInterfaces, iInterface

      nInterfaces = size(interfaceLocations, dim=1)

      if ( nInterfaces /= 43 ) then
         call mpas_log_write("MPAS-ocean: Vertical grid must have 42 layers "// &
                                      "to apply 42 Layer WOCE grid. Exiting...", MPAS_LOG_CRIT)
      end if

      interfaceLocations(1) = 0.0_RKIND
      interfaceLocations(2) = 5.00622_RKIND
      interfaceLocations(3) = 15.06873_RKIND
      interfaceLocations(4) = 25.28343_RKIND
      interfaceLocations(5) = 35.75849_RKIND
      interfaceLocations(6) = 46.61269_RKIND
      interfaceLocations(7) = 57.98099_RKIND
      interfaceLocations(8) = 70.02139_RKIND
      interfaceLocations(9) = 82.92409_RKIND
      interfaceLocations(10) = 96.92413_RKIND
      interfaceLocations(11) = 112.3189_RKIND
      interfaceLocations(12) = 129.4936_RKIND
      interfaceLocations(13) = 148.9582_RKIND
      interfaceLocations(14) = 171.4044_RKIND
      interfaceLocations(15) = 197.7919_RKIND
      interfaceLocations(16) = 229.4842_RKIND
      interfaceLocations(17) = 268.4617_RKIND
      interfaceLocations(18) = 317.6501_RKIND
      interfaceLocations(19) = 381.3864_RKIND
      interfaceLocations(20) = 465.9132_RKIND
      interfaceLocations(21) = 579.3073_RKIND
      interfaceLocations(22) = 729.3513_RKIND
      interfaceLocations(23) = 918.3723_RKIND
      interfaceLocations(24) = 1139.153_RKIND
      interfaceLocations(25) = 1378.574_RKIND
      interfaceLocations(26) = 1625.7_RKIND
      interfaceLocations(27) = 1875.106_RKIND
      interfaceLocations(28) = 2125.011_RKIND
      interfaceLocations(29) = 2375_RKIND
      interfaceLocations(30) = 2624.999_RKIND
      interfaceLocations(31) = 2874.999_RKIND
      interfaceLocations(32) = 3124.999_RKIND
      interfaceLocations(33) = 3374.999_RKIND
      interfaceLocations(34) = 3624.999_RKIND
      interfaceLocations(35) = 3874.999_RKIND
      interfaceLocations(36) = 4124.999_RKIND
      interfaceLocations(37) = 4374.999_RKIND
      interfaceLocations(38) = 4624.999_RKIND
      interfaceLocations(39) = 4874.999_RKIND
      interfaceLocations(40) = 5124.999_RKIND
      interfaceLocations(41) = 5374.999_RKIND
      interfaceLocations(42) = 5624.999_RKIND
      interfaceLocations(43) = 5874.999_RKIND

      maxInterfaceLocation = maxval(interfaceLocations)

      interfaceLocations(:) = interfaceLocations(:) / maxInterfaceLocation

   end subroutine ocn_generate_42layerWOCE_vertical_grid!}}}


   !***********************************************************************
   !
   !  routine ocn_generate_100layerE3SMv1_vertical_grid
   !
   !> \brief   100 vertical layer vertical grid generator for E3SM v1
   !> \author  Todd Ringler
   !> \date    04/23/2015
   !> \details
   !>  This routine generates a 100 layer grid
   !
   !-----------------------------------------------------------------------
   subroutine ocn_generate_100layerE3SMv1_vertical_grid(interfaceLocations)!{{{
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: interfaceLocations

      real (kind=RKIND) :: maxInterfaceLocation
      integer :: nInterfaces, iInterface

      nInterfaces = size(interfaceLocations, dim=1)

      if ( nInterfaces /= 101 ) then
         call mpas_log_write("MPAS-ocean: Vertical grid must have 100 layers to "// &
                                      "apply 100 Layer PHC grid. Exiting...", MPAS_LOG_CRIT)
      end if

      interfaceLocations(  1) =  0.0000E+00_RKIND
      interfaceLocations(  2) =  0.1510E+01_RKIND
      interfaceLocations(  3) =  0.3135E+01_RKIND
      interfaceLocations(  4) =  0.4882E+01_RKIND
      interfaceLocations(  5) =  0.6761E+01_RKIND
      interfaceLocations(  6) =  0.8779E+01_RKIND
      interfaceLocations(  7) =  0.1095E+02_RKIND
      interfaceLocations(  8) =  0.1327E+02_RKIND
      interfaceLocations(  9) =  0.1577E+02_RKIND
      interfaceLocations( 10) =  0.1845E+02_RKIND
      interfaceLocations( 11) =  0.2132E+02_RKIND
      interfaceLocations( 12) =  0.2440E+02_RKIND
      interfaceLocations( 13) =  0.2769E+02_RKIND
      interfaceLocations( 14) =  0.3122E+02_RKIND
      interfaceLocations( 15) =  0.3500E+02_RKIND
      interfaceLocations( 16) =  0.3904E+02_RKIND
      interfaceLocations( 17) =  0.4335E+02_RKIND
      interfaceLocations( 18) =  0.4797E+02_RKIND
      interfaceLocations( 19) =  0.5289E+02_RKIND
      interfaceLocations( 20) =  0.5815E+02_RKIND
      interfaceLocations( 21) =  0.6377E+02_RKIND
      interfaceLocations( 22) =  0.6975E+02_RKIND
      interfaceLocations( 23) =  0.7614E+02_RKIND
      interfaceLocations( 24) =  0.8294E+02_RKIND
      interfaceLocations( 25) =  0.9018E+02_RKIND
      interfaceLocations( 26) =  0.9790E+02_RKIND
      interfaceLocations( 27) =  0.1061E+03_RKIND
      interfaceLocations( 28) =  0.1148E+03_RKIND
      interfaceLocations( 29) =  0.1241E+03_RKIND
      interfaceLocations( 30) =  0.1340E+03_RKIND
      interfaceLocations( 31) =  0.1445E+03_RKIND
      interfaceLocations( 32) =  0.1556E+03_RKIND
      interfaceLocations( 33) =  0.1674E+03_RKIND
      interfaceLocations( 34) =  0.1799E+03_RKIND
      interfaceLocations( 35) =  0.1932E+03_RKIND
      interfaceLocations( 36) =  0.2072E+03_RKIND
      interfaceLocations( 37) =  0.2221E+03_RKIND
      interfaceLocations( 38) =  0.2379E+03_RKIND
      interfaceLocations( 39) =  0.2546E+03_RKIND
      interfaceLocations( 40) =  0.2722E+03_RKIND
      interfaceLocations( 41) =  0.2909E+03_RKIND
      interfaceLocations( 42) =  0.3106E+03_RKIND
      interfaceLocations( 43) =  0.3314E+03_RKIND
      interfaceLocations( 44) =  0.3534E+03_RKIND
      interfaceLocations( 45) =  0.3766E+03_RKIND
      interfaceLocations( 46) =  0.4011E+03_RKIND
      interfaceLocations( 47) =  0.4269E+03_RKIND
      interfaceLocations( 48) =  0.4541E+03_RKIND
      interfaceLocations( 49) =  0.4827E+03_RKIND
      interfaceLocations( 50) =  0.5128E+03_RKIND
      interfaceLocations( 51) =  0.5445E+03_RKIND
      interfaceLocations( 52) =  0.5779E+03_RKIND
      interfaceLocations( 53) =  0.6130E+03_RKIND
      interfaceLocations( 54) =  0.6498E+03_RKIND
      interfaceLocations( 55) =  0.6885E+03_RKIND
      interfaceLocations( 56) =  0.7291E+03_RKIND
      interfaceLocations( 57) =  0.7717E+03_RKIND
      interfaceLocations( 58) =  0.8164E+03_RKIND
      interfaceLocations( 59) =  0.8633E+03_RKIND
      interfaceLocations( 60) =  0.9124E+03_RKIND
      interfaceLocations( 61) =  0.9638E+03_RKIND
      interfaceLocations( 62) =  0.1018E+04_RKIND
      interfaceLocations( 63) =  0.1074E+04_RKIND
      interfaceLocations( 64) =  0.1133E+04_RKIND
      interfaceLocations( 65) =  0.1194E+04_RKIND
      interfaceLocations( 66) =  0.1259E+04_RKIND
      interfaceLocations( 67) =  0.1326E+04_RKIND
      interfaceLocations( 68) =  0.1396E+04_RKIND
      interfaceLocations( 69) =  0.1469E+04_RKIND
      interfaceLocations( 70) =  0.1546E+04_RKIND
      interfaceLocations( 71) =  0.1625E+04_RKIND
      interfaceLocations( 72) =  0.1708E+04_RKIND
      interfaceLocations( 73) =  0.1794E+04_RKIND
      interfaceLocations( 74) =  0.1884E+04_RKIND
      interfaceLocations( 75) =  0.1978E+04_RKIND
      interfaceLocations( 76) =  0.2075E+04_RKIND
      interfaceLocations( 77) =  0.2176E+04_RKIND
      interfaceLocations( 78) =  0.2281E+04_RKIND
      interfaceLocations( 79) =  0.2390E+04_RKIND
      interfaceLocations( 80) =  0.2503E+04_RKIND
      interfaceLocations( 81) =  0.2620E+04_RKIND
      interfaceLocations( 82) =  0.2742E+04_RKIND
      interfaceLocations( 83) =  0.2868E+04_RKIND
      interfaceLocations( 84) =  0.2998E+04_RKIND
      interfaceLocations( 85) =  0.3134E+04_RKIND
      interfaceLocations( 86) =  0.3274E+04_RKIND
      interfaceLocations( 87) =  0.3418E+04_RKIND
      interfaceLocations( 88) =  0.3568E+04_RKIND
      interfaceLocations( 89) =  0.3723E+04_RKIND
      interfaceLocations( 90) =  0.3882E+04_RKIND
      interfaceLocations( 91) =  0.4047E+04_RKIND
      interfaceLocations( 92) =  0.4218E+04_RKIND
      interfaceLocations( 93) =  0.4393E+04_RKIND
      interfaceLocations( 94) =  0.4574E+04_RKIND
      interfaceLocations( 95) =  0.4761E+04_RKIND
      interfaceLocations( 96) =  0.4953E+04_RKIND
      interfaceLocations( 97) =  0.5151E+04_RKIND
      interfaceLocations( 98) =  0.5354E+04_RKIND
      interfaceLocations( 99) =  0.5564E+04_RKIND
      interfaceLocations(100) =  0.5779E+04_RKIND
      interfaceLocations(101) =  0.6000E+04_RKIND

      maxInterfaceLocation = maxval(interfaceLocations)

      interfaceLocations(:) = interfaceLocations(:) / maxInterfaceLocation

   end subroutine ocn_generate_100layerE3SMv1_vertical_grid!}}}


!***********************************************************************
!
!  routine ocn_generate_1dCVT_vertical_grid
!
!> \brief   1D CVT vertical grid generator
!> \author  Juan A. Saenz
!> \date    09/10/2015
!> \details
!>  This routine generates a vertical grid with total depth = 1.
!>  This code is adapted from Todd's cvt_1d code.
!
!-----------------------------------------------------------------------

   subroutine ocn_generate_1dCVT_vertical_grid(interfaceLocations)!{{{

      real (kind=RKIND), dimension(:), intent(out) :: interfaceLocations

      integer :: k
      integer :: nInterfaces, nVertLevels

      real (kind=RKIND) :: stretch
      real (kind=RKIND) :: dz
      real (kind=RKIND) :: maxInterfaceLocation

      nInterfaces = size(interfaceLocations, dim=1)
      nVertLevels = nInterfaces - 1

      ! compute profile starting at top and stretch dz as we move down
      dz = config_1dCVTgenerator_dzSeed
      interfaceLocations(1) = 0.0_RKIND
      interfaceLocations(2) = dz
      do k=2,nVertLevels
         stretch = config_1dCVTgenerator_stretch1 + (config_1dCVTgenerator_stretch2-config_1dCVTgenerator_stretch1)*k/nVertLevels
         dz = stretch*dz
         interfaceLocations(k+1) = interfaceLocations(k) + dz
      enddo

      ! normalize so that positions span 0 to 1
      maxInterfaceLocation = maxval(interfaceLocations)
      interfaceLocations(:) = interfaceLocations(:) / maxInterfaceLocation

   end subroutine ocn_generate_1dCVT_vertical_grid!}}}


!***********************************************************************
!
!  routine ocn_compute_z_level_layerThickness
!
!> \brief   Compute auxiliary z-variables from bottomDepth
!> \author  Mark Petersen, Xylar Asay-Davis
!> \date    10/17/2015
!> \details
!>  This routine computes auxiliary z-variables from bottomDepth
!
!-----------------------------------------------------------------------

    subroutine ocn_compute_z_level_layerThickness(layerThickness, refBottomDepth, ssh, bottomDepth, &
                                                  minLevelCell, maxLevelCell, nVertLevels, iErr)!{{{
      real (kind=RKIND), dimension(nVertLevels), intent(out) :: layerThickness
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: refBottomDepth
      real (kind=RKIND), intent(in) :: ssh, bottomDepth
      integer, intent(in) :: minLevelCell, maxLevelCell, nVertLevels
      integer, intent(out) :: iErr

      integer :: k

      iErr = 0

      layerThickness(:) = 0.0_RKIND  ! = MPAS_REAL_FILLVAL

      if (maxLevelCell<=0) return

      ! compute the z-level layer thickness
      if (maxLevelCell==minLevelCell) then
         layerThickness(minLevelCell) = bottomDepth + ssh
      else
         k = minLevelCell
         layerThickness(k) = refBottomDepth(k) + ssh

         do k = minLevelCell+1, maxLevelCell-1
            layerThickness(k) = refBottomDepth(k) - refBottomDepth(k-1)
         end do

         k = maxLevelCell
         layerThickness(k) = bottomDepth - refBottomDepth(k-1)

      endif

    end subroutine ocn_compute_z_level_layerThickness  !}}}

!***********************************************************************
!
!  routine ocn_compute_zMid_from_layerThickness
!
!> \brief   Compute zMid from layerThickness
!> \author  Xylar Asay-Davis
!> \date    02/14/2021
!> \details
!>  This routine computes zMid from ssh and layerThickness for any vertical coordinate
!
!-----------------------------------------------------------------------

    subroutine ocn_compute_zMid_from_layerThickness(zMid, layerThickness, ssh, minLevelCell, &
                                                    maxLevelCell, nVertLevels, iErr)!{{{
      real (kind=RKIND), dimension(nVertLevels), intent(out) :: zMid
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: layerThickness
      real (kind=RKIND), intent(in) :: ssh
      integer, intent(in) :: minLevelCell, maxLevelCell, nVertLevels
      integer, intent(out) :: iErr

      integer :: k
      real (kind=RKIND) :: zTop

      iErr = 0

      zMid(:) = 0.0_RKIND  ! MPAS_REAL_FILLVAL

      if (maxLevelCell<=0) return

      zTop = ssh
      ! compute zMid based on the layer thickness
      do k = minLevelCell, maxLevelCell
         zMid(k) = zTop - 0.5_RKIND * layerThickness(k)
         zTop = zTop - layerThickness(k)
      end do

    end subroutine ocn_compute_zMid_from_layerThickness  !}}}


!***********************************************************************
!
!  routine ocn_compute_z_star_layerThickness
!
!> \brief   Compute z-star layerThickness
!> \author  Xylar Asay-Davis
!> \date    02/14/2021
!> \details
!>  This routine computes z-star layerThickness from restingThickness, bottomDepth, and ssh
!
!-----------------------------------------------------------------------

    subroutine ocn_compute_z_star_layerThickness(layerThickness, restingThickness, ssh, bottomDepth, &
                                                 minLevelCell, maxLevelCell, nVertLevels, iErr)!{{{
      real (kind=RKIND), dimension(nVertLevels), intent(out) :: layerThickness
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: restingThickness
      real (kind=RKIND), intent(in) :: ssh, bottomDepth
      integer, intent(in) :: minLevelCell, maxLevelCell, nVertLevels
      integer, intent(out) :: iErr

      integer :: k
      real (kind=RKIND) :: layerStretch

      iErr = 0

      layerThickness(:) = 0.0_RKIND  ! = MPAS_REAL_FILLVAL

      if (maxLevelCell<=0) return

      layerStretch = (ssh + bottomDepth)/bottomDepth
      do k = minLevelCell, maxLevelCell
         layerThickness(k) = layerStretch * restingThickness(k)
      end do

    end subroutine ocn_compute_z_star_layerThickness  !}}}


!***********************************************************************
!
!  routine ocn_alter_bottomDepth_for_pbcs
!
!> \brief   Alter bottom depth for partial bottom cells
!> \author  Mark Petersen
!> \date    10/19/2015
!> \details
!>  This routine alters the bottom depth in a single column based on partial cell settings
!
!-----------------------------------------------------------------------
    subroutine ocn_alter_bottomDepth_for_pbcs(bottomDepth, refBottomDepth, maxLevelCell, iErr)

      real (kind=RKIND), intent(inout) :: bottomDepth
      integer, intent(inout) :: maxLevelCell
      real (kind=RKIND), dimension(:), intent(in) :: refBottomDepth
      integer, intent(out) :: iErr
      integer :: k
      real (kind=RKIND) :: minBottomDepth, minBottomDepthMid

      iErr = 0

      if (maxLevelCell > 1) then
         if (config_alter_ICs_for_pcs) then

            if (config_pc_alteration_type .eq. 'partial_cell') then
               ! Change value of maxLevelCell for partial bottom cells
               k = maxLevelCell
               minBottomDepth = refBottomDepth(k) - (1.0-config_min_pc_fraction)*(refBottomDepth(k) - refBottomDepth(k-1))
               minBottomDepthMid = 0.5_RKIND*(minBottomDepth + refBottomDepth(k-1))
               if (bottomDepth .lt. minBottomDepthMid) then
                  ! Round up to cell above
                  maxLevelCell = maxLevelCell - 1
                  bottomDepth = refBottomDepth(maxLevelCell)
               else if (bottomDepth .lt. minBottomDepth) then
                  ! Round down cell to the min_pbc_fraction.
                  bottomDepth = minBottomDepth
               end if
            elseif (config_pc_alteration_type .eq. 'full_cell') then
               bottomDepth = refBottomDepth(maxLevelCell)
            else
               call mpas_log_write( ' Error: Incorrect choice of config_pc_alteration_type: '// config_pc_alteration_type)
               iErr = 1
            endif
         endif
      endif

    end subroutine ocn_alter_bottomDepth_for_pbcs

!***********************************************************************
!
!  routine ocn_alter_ssh_for_ptcs
!
!> \brief   Alter sea surface height for partial top cells
!> \author  Xylar Asay-Davis
!> \date    02/14/2021
!> \details
!>  This routine alters the sea surface height in a single column based on partial cell settings
!
!-----------------------------------------------------------------------
    subroutine ocn_alter_ssh_for_ptcs(ssh, refBottomDepth, minLevelCell, iErr)

      real (kind=RKIND), intent(inout) :: ssh
      integer, intent(inout) :: minLevelCell
      real (kind=RKIND), dimension(:), intent(in) :: refBottomDepth
      integer, intent(out) :: iErr
      real (kind=RKIND) :: minSSH, minSSHMid, zTop, zBot

      iErr = 0

      if (.not. config_alter_ICs_for_pcs) return

      if (minLevelCell < 1) return

      if (minLevelCell == 1) then
         zTop = 0.0_RKIND
      else
         zTop = -refBottomDepth(minLevelCell-1)
      end if

      zBot = -refBottomDepth(minLevelCell)

      if (config_pc_alteration_type == 'partial_cell') then
         ! Change value of minLevelCell for partial top cells
         minSSH = zBot + config_min_pc_fraction*(zTop - zBot)
         minSSHMid = 0.5_RKIND*(minSSH + zBot)
         if (ssh < minSSHMid) then
            ! Move down to the cell below
            minLevelCell = minLevelCell + 1
            ssh = zBot
         else if (ssh < minSSH) then
            ! Move up to min_pc_fraction.
            ssh = minSSH
         end if
      else if (config_pc_alteration_type == 'full_cell') then
         ssh = zTop
      else
         call mpas_log_write(' Error: Incorrect choice of config_pc_alteration_type: '// config_pc_alteration_type)
         iErr = 1
      endif

    end subroutine ocn_alter_ssh_for_ptcs


!***********************************************************************
!
!  routine ocn_compute_Haney_number
!
!> \brief   computes the Haney number (rx1)
!> \author  Xylar Asay-Davis
!> \date    11/20/2015
!> \details
!>  This routine computes the Haney number (rx1), which is a measure of
!>  hydrostatic consistency
!
!-----------------------------------------------------------------------
    subroutine ocn_compute_Haney_number(domain, iErr)

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, statePool
      real (kind=RKIND), dimension(:), pointer :: ssh, rx1MaxLevel

      integer, pointer :: nCells, nVertLevels, nEdges
      integer, dimension(:), pointer :: minLevelCell
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer :: iEdge, c1, c2, k, minLevelEdge, maxLevelEdge

      real (kind=RKIND) :: dzVert1, dzVert2, dzEdgeK, dzEdgeKp1, rx1, localMaxRx1Edge

      iErr = 0

      localMaxRx1Edge = 0.0_RKIND

      block_ptr => domain % blocklist
      call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      allocate(rx1MaxLevel(nVertLevels))
      rx1MaxLevel(:) = 0.0_RKIND


      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

        call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

        call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

        rx1Edge(:,:) = 0.0_RKIND
        rx1Cell(:,:) = 0.0_RKIND
        rx1MaxEdge(:) = 0.0_RKIND
        rx1MaxCell(:) = 0.0_RKIND
        do iEdge = 1,nEdges
          c1 = cellsOnEdge(1,iEdge)
          c2 = cellsOnEdge(2,iEdge)
          ! not a valid edge
          if((c1 > nCells) .or. (c2 > nCells)) cycle
          minLevelEdge = max(minLevelCell(c1), minLevelCell(c2))
          maxLevelEdge = min(maxLevelCell(c1), maxLevelCell(c2))
          do k = minLevelEdge,maxLevelEdge
            if(k == minLevelEdge) then
              dzVert1 = 2.0_RKIND*(ssh(c1)-zMid(k,c1))
              dzVert2 = 2.0_RKIND*(ssh(c2)-zMid(k,c2))
              dzEdgeK = ssh(c2)-ssh(c1)
            else
              dzVert1 = zMid(k-1,c1)-zMid(k,c1)
              dzVert2 = zMid(k-1,c2)-zMid(k,c2)
              dzEdgeK = zMid(k-1,c2)-zMid(k-1,c1)
            end if
            dzEdgeKp1 = zMid(k,c2)-zMid(k,c1)

            rx1 = abs(dzEdgeK+dzEdgeKp1)/(dzVert1+dzVert2)

            rx1Edge(k,iEdge) = rx1
            rx1MaxLevel(k) = max(rx1MaxLevel(k),rx1)
            rx1Cell(k,c1) = max(rx1Cell(k,c1),rx1)
            rx1Cell(k,c2) = max(rx1Cell(k,c2),rx1)

            rx1MaxEdge(iEdge) = max(rx1MaxEdge(iEdge),rx1)
            rx1MaxCell(c2) = max(rx1MaxCell(c2),rx1)
            rx1MaxCell(c1) = max(rx1MaxCell(c1),rx1)
          end do
        end do

        localMaxRx1Edge = max(localMaxRx1Edge,maxval(rx1MaxEdge))

        block_ptr => block_ptr % next
      end do
      do k = 1,nVertLevels
         call mpas_dmpar_max_real(domain % dminfo, rx1MaxLevel(k), globalRx1Max)
         call mpas_log_write ('  max of rx1 in level $i :  $r', intArgs=(/ k /),  realArgs=(/ globalRx1Max /))
      end do

      call mpas_dmpar_max_real(domain % dminfo, localMaxRx1Edge, globalRx1Max)
      call mpas_log_write ('global max of rx1: $r', realArgs=(/ globalRx1Max /))

      deallocate(rx1MaxLevel)

    end subroutine ocn_compute_Haney_number

!***********************************************************************
!
!  routine ocn_init_vertical_grid
!
!> \brief   Initialize z-star vertical grid based on SSH
!> \author  Xylar Asay-Davis
!> \date    8/8/2016
!> \details
!>  This routine sets up the vertical grid (layerThickness, zMid and
!>  restingThickness) needed for computing land-ice pressure from SSH.
!>  bottomDepth, refBottomDepth, must have been computed by the test case
!>  before calling this routine.  For coordinates other than z-level,
!>  minLevelCell and maxLevelCell must also be computed before calling this
!>  routine.  For z-level, they are computed here.  For tests with ice-shelf
!>  cavities, ssh must be computed before calling this routine. If the test
!>  case is using a z-star or z-level coordinate, this routine will take care
!>  of setting  up partial bottom (and, for z-level, top) cells by calling
!>  ocn_alter_bottomDepth_for_pbcs() and ocn_alter_ssh_for_ptcs() . If using
!>  the Haney-number-constrained coordinate, another approach is used to handle
!>  maintaining a minimum layer thickness.
!-----------------------------------------------------------------------

   subroutine ocn_init_vertical_grid(domain, iErr)!{{{

   !--------------------------------------------------------------------

     type (domain_type), intent(inout) :: domain
     integer, intent(out) :: iErr

     type (block_type), pointer :: block_ptr

     type (mpas_pool_type), pointer :: meshPool, statePool, verticalMeshPool

     ! Define dimension pointers
     integer, pointer :: nCells, nVertLevels

     ! Define variable pointers
     integer, dimension(:), pointer :: minLevelCell, maxLevelCell
     real (kind=RKIND), dimension(:), pointer :: refBottomDepth, bottomDepth, ssh
     real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness

     integer :: iCell, k

     logical :: useZStar, useZLevel

     real (kind=RKIND) :: depth, layerStretch

   !--------------------------------------------------------------------

     iErr = 0

     if (config_init_vertical_grid_type /= 'z-level' .and. &
         config_init_vertical_grid_type /= 'z-star' .and. &
         config_init_vertical_grid_type /= 'haney-number') then

         call mpas_log_write('Unexpected value for config_init_vertical_grid_type ' // &
                             trim(config_init_vertical_grid_type), MPAS_LOG_CRIT)

     end if

     if(config_init_vertical_grid_type == 'haney-number') then
       ! We already know the ssh and landIcePressure we want to use.
       ! Compute the layer thicknesses and zMid based on topography and ssh.
       ! Use rx1 constraint to recompute the vertical grid.
       call ocn_init_vertical_grid_with_max_rx1(domain, iErr)

       if(iErr .ne. 0) then
         call mpas_log_write( 'ocn_init_vertical_grid_with_max_rx1 failed.', MPAS_LOG_CRIT)
         return
       end if

     else if(config_init_vertical_grid_type == 'z-star') then

       block_ptr => domain % blocklist
       do while(associated(block_ptr))

         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

         call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

         do iCell = 1, nCells
           ! alter pbcs
           call ocn_alter_bottomDepth_for_pbcs(bottomDepth(iCell), refBottomDepth, maxLevelCell(iCell), iErr)
           if(iErr .ne. 0) then
             call mpas_log_write( 'ocn_alter_bottomDepth_for_pbcs failed.', MPAS_LOG_CRIT)
             return
           end if

           ! restingThickness is z-level, with ssh = 0
           call ocn_compute_z_level_layerThickness(restingThickness(:,iCell), refBottomDepth, 0.0_RKIND, &
                                                   bottomDepth(iCell), minLevelCell(iCell),              &
                                                   maxLevelCell(iCell), nVertLevels, iErr)

           if(iErr .ne. 0) then
             call mpas_log_write( 'ocn_compute_z_level_layerThickness failed.', MPAS_LOG_CRIT)
             return
           end if

           ! stretch restingThickness to get the z-star layerThickness
           call ocn_compute_z_star_layerThickness(layerThickness(:,iCell), restingThickness(:,iCell),  &
                                                  ssh(iCell), bottomDepth(iCell), minLevelCell(iCell), &
                                                  maxLevelCell(iCell), nVertLevels, iErr)

           if(iErr .ne. 0) then
             call mpas_log_write( 'ocn_compute_z_star_layerThickness failed.', MPAS_LOG_CRIT)
             return
           end if

           ! compute zMid
           call ocn_compute_zMid_from_layerThickness(zMid(:,iCell), layerThickness(:,iCell), ssh(iCell), &
                                                     minLevelCell(iCell), maxLevelCell(iCell), nVertLevels, iErr)

           if(iErr .ne. 0) then
             call mpas_log_write( 'ocn_compute_zMid_from_layerThickness failed.', MPAS_LOG_CRIT)
             return
           end if
         end do !iCell

         block_ptr => block_ptr % next
       end do !block_ptr
     else  ! config_init_vertical_grid_type == 'z-level'

       block_ptr => domain % blocklist
       do while(associated(block_ptr))

         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
         call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

         call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
         call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
         call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

         call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
         call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

         call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

         ! compute minLevelCell and maxLevelCell based on refBottomDepth, ssh and bottomDepth
         do iCell = 1, nCells
           maxLevelCell(iCell) = -1
           minLevelCell(iCell) = -1

           if ((bottomDepth(iCell) == 0.0_RKIND) .or. (bottomDepth(iCell) <= -ssh(iCell))) then
              cycle
           end if

           do k = 1, nVertLevels
             depth = refBottomDepth(k)
             if (depth >= bottomDepth(iCell)) then
               maxLevelCell(iCell) = k
               exit
             end if
           end do

           do k = 1, nVertLevels
             depth = -refBottomDepth(k)
             if (depth < ssh(iCell)) then
               minLevelCell(iCell) = k
               exit
             end if
           end do
         end do

         do iCell = 1, nCells
           ! alter pbcs
           call ocn_alter_bottomDepth_for_pbcs(bottomDepth(iCell), refBottomDepth, maxLevelCell(iCell), iErr)
           if(iErr .ne. 0) then
             call mpas_log_write( 'ocn_alter_bottomDepth_for_pbcs failed.', MPAS_LOG_CRIT)
             return
           end if

           ! alter ptcs
           call ocn_alter_ssh_for_ptcs(ssh(iCell), refBottomDepth, minLevelCell(iCell), iErr)
           if(iErr .ne. 0) then
             call mpas_log_write( 'ocn_alter_ssh_for_ptcs failed.', MPAS_LOG_CRIT)
             return
           end if

           ! layerThickness is z-level
           call ocn_compute_z_level_layerThickness(layerThickness(:,iCell), refBottomDepth, ssh(iCell), &
                                                   bottomDepth(iCell), minLevelCell(iCell),             &
                                                   maxLevelCell(iCell), nVertLevels, iErr)

           if(iErr .ne. 0) then
             call mpas_log_write( 'ocn_compute_z_level_layerThickness failed.', MPAS_LOG_CRIT)
             return
           end if

           ! cases with SSH > 0 require a resting thickness that is stretched to compensate for the SSH
           restingThickness(:, iCell) = 0.0_RKIND

           if (maxLevelCell(iCell)  > 0) then
             layerStretch = (ssh(iCell) + bottomDepth(iCell))/bottomDepth(iCell)
             do k = minLevelCell(iCell), maxLevelCell(iCell)
               restingThickness(k, iCell) = layerThickness(k, iCell)/layerStretch
             end do
           end if

           ! compute zMid
           call ocn_compute_zMid_from_layerThickness(zMid(:,iCell), layerThickness(:,iCell), ssh(iCell),    &
                                                     minLevelCell(iCell), maxLevelCell(iCell), nVertLevels, &
                                                     iErr)

           if(iErr .ne. 0) then
             call mpas_log_write( 'ocn_compute_zMid_from_layerThickness failed.', MPAS_LOG_CRIT)
             return
           end if
         end do !iCell

         block_ptr => block_ptr % next
       end do !block_ptr
     end if

   end subroutine ocn_init_vertical_grid

!***********************************************************************
!
!  routine ocn_init_vertical_grid_with_max_rx1
!
!> \brief   re-initializes the vertical grid so rx1 < rx1Max
!> \author  Xylar Asay-Davis
!> \date    11/23/2015
!> \details
!>  This routine re-initializes the vertical grid (layerThickness,
!>  restingThickness maxLevelCell, zMid) so that the Haney number is
!>  less than a maximum value (rx1 < rx1Max). ssh and bottomDepth should
!>  have been initialized before calling this routine.  bottomDepth will
!>  be modified for full or partial bottom cells in this routine, so
!>  this step should not be performed before calling this routine
!
!-----------------------------------------------------------------------
    subroutine ocn_init_vertical_grid_with_max_rx1(domain, iErr)

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, statePool, verticalMeshPool, scratchPool, forcingPool

      integer, pointer :: config_rx1_outer_iter_count, config_rx1_inner_iter_count, &
                          config_rx1_horiz_smooth_open_ocean_cells, config_rx1_min_levels
      real (kind=RKIND), pointer :: config_rx1_max, config_rx1_horiz_smooth_weight, &
                           config_rx1_vert_smooth_weight, config_rx1_slope_weight, &
                           config_rx1_zstar_weight, config_rx1_init_inner_weight, &
                           config_rx1_min_layer_thickness

      type (field2DReal), pointer :: zInterfaceField, goalStretchField, goalWeightField
      type (field1DReal), pointer :: zTopField, zBotField, zBotNewField
      type (field1DInteger), pointer :: maxLevelCellField
      type (field1DInteger), pointer :: smoothingMaskNewField

      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness, zInterface, &
                                                    goalStretch, goalWeight
      integer, dimension(:), pointer :: landIceMask
      real (kind=RKIND), dimension(:), pointer :: ssh, bottomDepth, refBottomDepth, zTop, zBot, zBotNew, &
                                        refLayerThickness

      integer, pointer :: nCells, nVertLevels, nEdges, nCellsSolve
      integer, dimension(:), pointer :: maxLevelCell, cullCell, nEdgesOnCell, smoothingMaskNew
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnCell, edgesOnCell

      integer :: iCell, iEdge, coc, c1, c2, k, maxLevelEdge, iSmooth, iterIndex, maxLevelNeighbors

      real (kind=RKIND) :: dzEdgeK, dzEdgeKp1, dzEdgeMean, dzVertGoal, &
                           zMean, weight, rx1Goal, dzVertMean, &
                           zMidNext, frac, stretch, localMaxRx1Edge, &
                           localStretchMax, localStretchMin

      logical :: moveInterface

      iErr = 0

      call mpas_pool_get_config(domain % configs, 'config_rx1_outer_iter_count', config_rx1_outer_iter_count)
      call mpas_pool_get_config(domain % configs, 'config_rx1_inner_iter_count', config_rx1_inner_iter_count)
      call mpas_pool_get_config(domain % configs, 'config_rx1_init_inner_weight', config_rx1_init_inner_weight)
      call mpas_pool_get_config(domain % configs, 'config_rx1_max', config_rx1_max)
      call mpas_pool_get_config(domain % configs, 'config_rx1_horiz_smooth_weight', config_rx1_horiz_smooth_weight)
      call mpas_pool_get_config(domain % configs, 'config_rx1_vert_smooth_weight', config_rx1_vert_smooth_weight)
      call mpas_pool_get_config(domain % configs, 'config_rx1_slope_weight', config_rx1_slope_weight)
      call mpas_pool_get_config(domain % configs, 'config_rx1_zstar_weight', config_rx1_zstar_weight)
      call mpas_pool_get_config(domain % configs, 'config_rx1_horiz_smooth_open_ocean_cells', &
                                config_rx1_horiz_smooth_open_ocean_cells)
      call mpas_pool_get_config(domain % configs, 'config_rx1_min_levels', config_rx1_min_levels)
      call mpas_pool_get_config(domain % configs, 'config_rx1_min_layer_thickness', config_rx1_min_layer_thickness)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      ! allocate scratch variables that persist across blocks
      call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
      call mpas_pool_get_field(scratchPool, 'zInterfaceScratch', zInterfaceField)
      call mpas_pool_get_field(scratchPool, 'zTopScratch', zTopField)
      call mpas_pool_get_field(scratchPool, 'zBotScratch', zBotField)
      call mpas_allocate_scratch_field(zInterfaceField, .false.)
      call mpas_allocate_scratch_field(zTopField, .false.)
      call mpas_allocate_scratch_field(zBotField, .false.)

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(meshPool, 'cullCell', cullCell)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)

        maxLevelCell(nCells+1) = -1
        do iCell = 1, nCells
          if(associated(cullCell)) then
            if(cullCell(iCell) == 1) then
              ! we need to know to ignore this cell later
              maxLevelCell(iCell) = -1
            end if
          end if
        end do

        ! initialize the smoothing mask to valid cells under land ice
        smoothingMask(:) = 0
        where((maxLevelCell(:) > 0) .and. (landIceMask(:) == 1))
          smoothingMask(:) = 1
        end where

        block_ptr => block_ptr % next
      end do !block_ptr

      ! expand the smoothing mask to neighbors of land-ice cells
      do iSmooth = 1, config_rx1_horiz_smooth_open_ocean_cells
        block_ptr => domain % blocklist
        do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

          call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
          call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

          call mpas_pool_get_field(scratchPool, 'smoothingMaskNewScratch', smoothingMaskNewField)
          call mpas_allocate_scratch_field(smoothingMaskNewField, .true.)
          call mpas_pool_get_array(scratchPool, 'smoothingMaskNewScratch', smoothingMaskNew)

          smoothingMaskNew(:) = smoothingMask(:)

          ! expand the mask to neighbors
          do iCell = 1, nCells
            if(smoothingMask(iCell) == 0) cycle

            do iEdge = 1, nEdgesOnCell(iCell)
              coc = cellsOnCell(iEdge,iCell)
              if(coc > 0) then
                if((maxLevelCell(coc) > 0) .and. (smoothingMaskNew(coc) == 0)) then
                  ! we have a neighbor of a cell being smoothed so this one should also be smoothed
                  smoothingMaskNew(coc) = 1
                end if
              end if
            end do !iEdge
          end do !iCell

          smoothingMask(:) = smoothingMaskNew(:)
          call mpas_deallocate_scratch_field(smoothingMaskNewField, .true.)

          block_ptr => block_ptr % next
        end do !block_ptr

        ! do halo update on smoothingMask
        call mpas_dmpar_exch_halo_field(smoothingMaskField)
      end do !iSmooth

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
        call mpas_pool_get_array(scratchPool, 'zInterfaceScratch', zInterface)

        call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

        do iCell = 1, nCells
          if(maxLevelCell(iCell) <= 0) then
            ! this is land
            zInterface(:,iCell) = 0.0_RKIND
            verticalStretch(:,iCell) = 1.0_RKIND
            cycle
          end if

          ! initialize zInterface to z* without PBCs and extended to nVertLevels+1 (beyond bottomDepth)

          ! don't let bottomDepth go below the z-level grid
          bottomDepth(iCell) = min(bottomDepth(iCell), refBottomDepth(nVertLevels))

          ! lower bottomDepth if the whole column is thinner than the minimum
          !bottomDepth(iCell) = max(bottomDepth(iCell), -ssh(iCell) + config_rx1_min_layer_thickness*config_rx1_min_levels)

          verticalStretch(:,iCell) = (ssh(iCell) + bottomDepth(iCell))/bottomDepth(iCell)

          zInterface(1,iCell) = ssh(iCell)
          do k = 1, nVertLevels
            zInterface(k+1,iCell) = zInterface(k,iCell) - verticalStretch(k,iCell)*refLayerThickness(k)
          end do

        end do !iCell

        block_ptr => block_ptr % next
      end do !block_ptr


      do iterIndex = 1, config_rx1_outer_iter_count

        ! smooth/nudge twice so changes propagate further per halo update
        do iSmooth = 1, 2
          block_ptr => domain % blocklist
          do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
            call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

            call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
            call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
            call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

            call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
            call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
            call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
            call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

            call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

            call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

            call mpas_pool_get_field(scratchPool, 'goalStretchScratch', goalStretchField)
            call mpas_allocate_scratch_field(goalStretchField, .true.)
            call mpas_pool_get_array(scratchPool, 'goalStretchScratch', goalStretch)
            call mpas_pool_get_field(scratchPool, 'goalWeightScratch', goalWeightField)
            call mpas_allocate_scratch_field(goalWeightField, .true.)
            call mpas_pool_get_array(scratchPool, 'goalWeightScratch', goalWeight)

            goalStretch(:,:) = 0.0_RKIND
            goalWeight(:,:) = 0.0_RKIND

            do iCell = 1, nCells
              if(smoothingMask(iCell) == 0) cycle

              do iEdge = 1, nEdgesOnCell(iCell)
                coc = cellsOnCell(iEdge,iCell)
                if(maxLevelCell(coc) <= 0) cycle

                ! average horiz neighbors
                goalStretch(:, iCell) = goalStretch(:, iCell) + config_rx1_horiz_smooth_weight*verticalStretch(:, coc)
                goalWeight(:, iCell) = goalWeight(:, iCell) + config_rx1_horiz_smooth_weight

                ! change stretch toward flatter level interfaces
                do k = 1, nVertLevels
                  ! we want to try to move the bottom interface to zMean so it's more level
                  zMean = 0.5_RKIND*(zInterface(k+1,iCell) + zInterface(k+1,coc))
                  ! frac is the factor by which we want to modify the stretch above (and including) this cell
                  frac = (ssh(iCell)-zMean)/(ssh(iCell)-zInterface(k+1,iCell))
                  goalStretch(1:k, iCell) = goalStretch(1:k, iCell) &
                                          + config_rx1_slope_weight*frac*verticalStretch(1:k, iCell)
                  goalWeight(1:k, iCell) = goalWeight(1:k, iCell) + config_rx1_slope_weight
                end do
              end do !iEdge

              ! include this cell in average
              goalStretch(:, iCell) = goalStretch(:, iCell) + verticalStretch(:, iCell)
              goalWeight(:, iCell) = goalWeight(:, iCell) + 1.0_RKIND

              ! average vert neighbors
              goalStretch(1:nVertLevels-1, iCell) = goalStretch(1:nVertLevels-1, iCell) &
                + config_rx1_vert_smooth_weight*verticalStretch(2:nVertLevels, iCell)
              goalWeight(1:nVertLevels-1, iCell) = goalWeight(1:nVertLevels-1, iCell) &
                + config_rx1_vert_smooth_weight
              goalStretch(2:nVertLevels, iCell) = goalStretch(2:nVertLevels, iCell) &
                + config_rx1_vert_smooth_weight*verticalStretch(1:nVertLevels-1, iCell)
              goalWeight(2:nVertLevels, iCell) = goalWeight(2:nVertLevels, iCell) &
                + config_rx1_vert_smooth_weight

              ! nudge toward z-star
              stretch = (ssh(iCell) + bottomDepth(iCell))/bottomDepth(iCell)
              goalStretch(:, iCell) = goalStretch(:, iCell) + config_rx1_zstar_weight*stretch
              goalWeight(:, iCell) = goalWeight(:, iCell) + config_rx1_zstar_weight

            end do !iCell

            do iCell = 1, nCells
              if(smoothingMask(iCell) == 1) then
                do k = 1, nVertLevels
                  ! minimum allowed stretch
                  stretch = config_rx1_min_layer_thickness/refLayerThickness(k)
                  verticalStretch(k,iCell) = max(stretch, goalStretch(k,iCell)/goalWeight(k,iCell))
                end do
              end if
            end do

            call mpas_deallocate_scratch_field(goalStretchField, .true.)
            call mpas_deallocate_scratch_field(goalWeightField, .true.)

            block_ptr => block_ptr % next
          end do !block_ptr
        end do !iSmooth

        ! do a halo exchange on verticalStretch
        call mpas_dmpar_exch_halo_field(verticalStretchField)

        block_ptr => domain % blocklist
        do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_array(scratchPool, 'zTopScratch', zTop)
          call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

          zTop(:) = ssh(:)
          where(smoothingMask == 1)
            ! start with maxLevelCell == nVertLevels+1; we will update it when we encounter the bottom
            maxLevelCell(:) = nVertLevels+1
          end where

          block_ptr => block_ptr % next
        end do

        do k = 1, nVertLevels
          block_ptr => domain % blocklist
          do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

            call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

            call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
            call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

            call mpas_pool_get_array(scratchPool, 'zInterfaceScratch', zInterface)
            call mpas_pool_get_array(scratchPool, 'zBotScratch', zBot)
            call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

            zInterface(k+1,:) = zInterface(k,:) - verticalStretch(k,:)*refLayerThickness(k)

            ! match interfaces to bottomDepth
            if(k >= config_rx1_min_levels) then
              moveInterface = .false.
              if(config_rx1_outer_iter_count == 1) then
                weight = 1.0_RKIND
              else
                weight = (iterIndex - 1.0_RKIND)/(config_rx1_outer_iter_count - 1.0_RKIND)
              end if

              do iCell = 1, nCells
                if((smoothingMask(iCell) == 0) .or. (maxLevelCell(iCell) .ne. nVertLevels+1)) cycle

                if(-bottomDepth(iCell) > zBot(iCell)) then
                  ! we missed maxLevelCell, either because k < config_rx1_min_levels or
                  ! because our guess at the layer depth wasn't accurate.  We can't move
                  ! the layer interface to match bottomDepth but we can at least set
                  ! maxLevelCell to the appropriate value
                  maxLevelCell(iCell) = k
                  cycle
                end if

                if(k == nVertLevels) then
                  ! match the bottom of the layer to bottomDepth
                  moveInterface = .true.
                else
                  zInterface(k+1,iCell) = zInterface(k,iCell) - verticalStretch(k,iCell)*refLayerThickness(k)
                  ! our current best guess at the mid depth of the next layer
                  zMidNext = zInterface(k+1,iCell) - 0.5_RKIND*verticalStretch(k+1,iCell)*refLayerThickness(k+1)
                  moveInterface = -bottomDepth(iCell) >= zMidNext
                end if
                ! relax toward bottomDepth with increasing strength with each iteration
                if(moveInterface) then
                  zInterface(k+1,iCell) = min((1.0_RKIND - weight)*zInterface(k+1,iCell) + weight*(-bottomDepth(iCell)), &
                                              zInterface(k,iCell) - config_rx1_min_layer_thickness)
                  maxLevelCell(iCell) = k
                end if
              end do
            end if
            zBot(:) = 0.5_RKIND*(zInterface(k,:) + zInterface(k+1,:))

            block_ptr => block_ptr % next
          end do !block_ptr

          if(k == 1) then
            ! rx1 is allowed to get twice as big in the top layer because we're only looking at half a layer
            rx1Goal = 2.0_RKIND*config_rx1_max
          else
            rx1Goal = config_rx1_max
          end if
          call constrain_rx1_layer(domain, config_rx1_inner_iter_count, config_rx1_init_inner_weight, &
                                   rx1Goal, k > config_rx1_min_levels, iErr)

          ! update zInterface, zTop
          block_ptr => domain % blocklist
          do while(associated(block_ptr))
            call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)
            call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

            call mpas_pool_get_array(scratchPool, 'zInterfaceScratch', zInterface)
            call mpas_pool_get_array(scratchPool, 'zTopScratch', zTop)
            call mpas_pool_get_array(scratchPool, 'zBotScratch', zBot)
            call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

            where (smoothingMask(:) == 1)
              verticalStretch(k,:) = 2.0_RKIND*(zInterface(k,:) - zBot(:))/refLayerThickness(k)
              zInterface(k+1,:) = 2.0_RKIND*zBot(:) - zInterface(k,:)
            end where

            zTop(:) = zBot(:)
            zMid(k,:) = zBot(:)

            block_ptr => block_ptr % next
          end do !block_ptr

        end do !k

        block_ptr => domain % blocklist
        do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)

          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
          call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
          call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)

          ! Fill bathymetry holes, i.e. single cells that are deeper than all neighbors.
          ! These cells can collect tracer extrema and are not able to use
          ! horizontal diffusion or advection to clear them out. Reduce pits to
          ! make them level with the next deepest neighbor cell.

          maxLevelCell(nCells+1) = -1

          do iCell = 1, nCellsSolve
             maxLevelNeighbors = 0
             do iEdge = 1, nEdgesOnCell(iCell)
                coc = cellsOnCell(iEdge, iCell)
                maxLevelNeighbors = max(maxLevelNeighbors, maxLevelCell(coc))
             end do

             if (maxLevelCell(iCell) > maxLevelNeighbors) then
                maxLevelCell(iCell) = maxLevelNeighbors
             end if
          end do

          where((smoothingMask(:) == 1) .and. (maxLevelCell(:) == nVertLevels+1))
            ! we never found the bottom, so it must be the last level
            maxLevelCell(:) = nVertLevels
          end where

          block_ptr => block_ptr % next
        end do

       call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
       call mpas_pool_get_field(meshPool, 'maxLevelCell', maxLevelCellField)
       call mpas_dmpar_exch_halo_field(maxLevelCellField)

        localMaxRx1Edge = -1e30_RKIND
        localStretchMin = 1e30_RKIND
        localStretchMax = -1e30_RKIND

        block_ptr => domain % blocklist
        do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
          call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

          call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

          call mpas_pool_get_array(statePool, 'ssh', ssh, 1)

          ! compute rx1Edge so we can determine which cells need to be thickened
          rx1Edge(:,:) = 0.0_RKIND
          do iEdge = 1,nEdges
            c1 = cellsOnEdge(1,iEdge)
            c2 = cellsOnEdge(2,iEdge)
            if((c1 <= 0) .or. (c2 <= 0) .or. (c1 > nCells) .or. (c2 > nCells)) cycle
            if((maxLevelCell(c1) <= 0) .or. (maxLevelCell(c2) <= 0)) cycle

            ! at the top level, use ssh instead of zMid
            dzVertMean = 0.5_RKIND*(ssh(c1)-zMid(1,c1)+ssh(c2)-zMid(1,c2))
            dzEdgeMean = 0.5_RKIND*abs(ssh(c2)-ssh(c1)+zMid(1,c2)-zMid(1,c1))
            ! a factor of 0.5 because ssh is at top interface, not middle of the previous layer
            rx1Edge(1,iEdge) = 0.5_RKIND*dzEdgeMean/dzVertMean
            maxLevelEdge = min(maxLevelCell(c1), maxLevelCell(c2))
            do k = 1, maxLevelEdge-1
              dzVertMean = 0.5_RKIND*(zMid(k,c1)-zMid(k+1,c1)+zMid(k,c2)-zMid(k+1,c2))
              dzEdgeMean = 0.5_RKIND*abs(zMid(k,c2)-zMid(k,c1)+zMid(k+1,c2)-zMid(k+1,c1))
              rx1Edge(k+1,iEdge) = dzEdgeMean/dzVertMean
            end do
          end do

          localMaxRx1Edge = max(localMaxRx1Edge,maxval(rx1Edge))

          do iCell = 1, nCells
            if(maxLevelCell(iCell) <= 0) cycle
            do k = 1, maxLevelCell(iCell)
              localStretchMax = max(localStretchMax, verticalStretch(k,iCell))
              localStretchMin = min(localStretchMin, verticalStretch(k,iCell))
            end do
          end do

          block_ptr => block_ptr % next
        end do
        call mpas_dmpar_max_real(domain % dminfo, localMaxRx1Edge, globalRx1Max)
        call mpas_log_write (' iter:  $i global max of rx1 $r', intArgs=(/ iterIndex /),  realArgs=(/ globalRx1Max /))
        call mpas_dmpar_min_real(domain % dminfo, localStretchMin, globalVerticalStretchMin)
        call mpas_log_write ('            global min of verticalStretch: $r', realArgs=(/ globalVerticalStretchMin /))
        call mpas_dmpar_max_real(domain % dminfo, localStretchMax, globalVerticalStretchMax)
        call mpas_log_write ('            global max of verticalStretch: $r', realArgs=(/ globalVerticalStretchMax /))

      end do !iterIndex

      ! compute maxLevelCell, zMid and restingThickness; update bottomDepth and layerThickness
      ! for full or partial bottom cells (if requested)
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
        call mpas_pool_get_array(statePool, 'ssh', ssh, 1)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
        call mpas_pool_get_array(scratchPool, 'zInterfaceScratch', zInterface)
        call mpas_pool_get_array(verticalMeshPool, 'refLayerThickness', refLayerThickness)

        ! compute zMid, layerThickness and restingThickness
        do iCell = 1, nCells
          if(maxLevelCell(iCell) == -1) then
            bottomDepth(iCell) = 0.0_RKIND
            zMid(:,iCell) = 0.0_RKIND
            layerThickness(:,iCell) = 0.0_RKIND
            restingThickness(:,iCell) = 0.0_RKIND
            cycle
          end if


          if(smoothingMask(iCell) == 0) then
            ! initialize with PBCs
            call ocn_alter_bottomDepth_for_pbcs(bottomDepth(iCell), refBottomDepth, maxLevelCell(iCell), iErr)
            if(iErr .ne. 0) then
              call mpas_log_write( 'ocn_alter_bottomDepth_for_pbcs failed.', MPAS_LOG_CRIT)
              return
            end if
          else
            ! we may not have been able to put the layer interface exactly at bottomDepth,
            ! either because bottomDepth was too shallow for the minimum number of layers
            ! or because contracting the layer would have led to rx1 > rx1Max
            bottomDepth(iCell) = -zInterface(maxLevelCell(iCell)+1,iCell)
          end if

          ! don't allow layers to go below -bottomDepth
          zInterface(:,iCell) = max(zInterface(:,iCell),-bottomDepth(iCell))

          zMid(:,iCell) = 0.5_RKIND*(zInterface(1:nVertLevels,iCell) + zInterface(2:nVertLevels+1,iCell))
          layerThickness(:,iCell) = zInterface(1:nVertLevels,iCell) - zInterface(2:nVertLevels+1,iCell)
          verticalStretch(:,iCell) = layerThickness(:,iCell)/refLayerThickness(:)

          !restingThickness can be computed by "undoing" the z* stretch
          stretch = (ssh(iCell) + bottomDepth(iCell))/bottomDepth(iCell)
          zInterface(:,iCell) = (zInterface(:,iCell) + bottomDepth(iCell))/stretch - bottomDepth(iCell)
          restingThickness(:,iCell) = zInterface(1:nVertLevels,iCell) - zInterface(2:nVertLevels+1,iCell)
        end do

        block_ptr => block_ptr % next
      end do !block_ptr

      call mpas_deallocate_scratch_field(zInterfaceField, .false.)
      call mpas_deallocate_scratch_field(zTopField, .false.)
      call mpas_deallocate_scratch_field(zBotField, .false.)

    end subroutine ocn_init_vertical_grid_with_max_rx1

!***********************************************************************
!
!  routine constrain_rx1_layer
!
!> \brief   modify zBot in a layer so rx1 <= rx1Goal
!> \author  Xylar Asay-Davis
!> \date    05/09/2015
!> \details
!>  This routine is used to iteratively constrains zBot such that
!>  rx1 <= rx1Goal.  zBot is nudged toward a goal field initially
!>  with weight initIterWeight and finally with a weight of 1.0.
!>
!
!-----------------------------------------------------------------------

    subroutine constrain_rx1_layer(domain, iterCount, initIterWeight, rx1Goal, checkBelowBottom, iErr)

      type (domain_type), intent(inout) :: domain
      integer, intent(in) :: iterCount
      real (kind=RKIND), intent(in) :: initIterWeight, rx1Goal
      logical :: checkBelowBottom
      integer, intent(out) :: iErr

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, scratchPool

      type (field1DReal), pointer :: zBotField, zBotNewField

      real (kind=RKIND), dimension(:), pointer :: zTop, zBot, zBotNew, bottomDepth

      integer, pointer :: nCells, nVertLevels, nEdges
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer :: iCell, iEdge, c1, c2, iterIndex

      real (kind=RKIND) :: dzEdgeK, dzEdgeKp1, dzEdgeMean, dzVertGoal, weight, zBotEdge, &
                           deltaZBot

      iErr = 0

      call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
      call mpas_pool_get_field(scratchPool, 'zBotScratch', zBotField)

      do iterIndex = 1, iterCount
        ! next, adjust zBot toward zBotNew (with rx1 < rx1Max)
        if(iterCount == 1) then
          weight = 1.0_RKIND
        else
          weight = (iterIndex - 1.0_RKIND)/(iterCount - 1.0_RKIND)
        end if
        ! the weight goes from initIterWeight for the first iteration to 1.0 for the last
        weight = (1.0_RKIND - weight)*initIterWeight + weight

        block_ptr => domain % blocklist
        do while(associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call mpas_pool_get_subpool(block_ptr % structs, 'scratch', scratchPool)

          call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
          call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
          call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

          call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
          call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
          call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

          call mpas_pool_get_array(scratchPool, 'zTopScratch', zTop)
          call mpas_pool_get_array(scratchPool, 'zBotScratch', zBot)

          call mpas_pool_get_field(scratchPool, 'zBotNewScratch', zBotNewField)
          call mpas_allocate_scratch_field(zBotNewField, .true.)
          call mpas_pool_get_array(scratchPool, 'zBotNewScratch', zBotNew)

          zBotNew(:) = zBot(:)
          do iEdge = 1,nEdges
            c1 = cellsOnEdge(1,iEdge)
            c2 = cellsOnEdge(2,iEdge)
            if((c1 <= 0) .or. (c1 > nCells) .or. (c2 <= 0) .or. (c2 > nCells)) cycle
            if((maxLevelCell(c1) <= 0) .or. (maxLevelCell(c2) <= 0)) cycle
            if((smoothingMask(c1) == 0) .and. (smoothingMask(c2) == 0)) cycle

            if(checkBelowBottom) then
              ! if both cells are definitely below the bathymetry, no need to constrain rx1
              if((zTop(c1) < -bottomDepth(c1)) .and. (zTop(c2) < -bottomDepth(c2))) cycle
            end if

            dzEdgeK = zTop(c2)-zTop(c1)
            dzEdgeKp1 = zBot(c2)-zBot(c1)
            dzEdgeMean = 0.5_RKIND*abs(dzEdgeK+dzEdgeKp1)
            dzVertGoal = dzEdgeMean/rx1Goal
            zBotEdge = 0.5_RKIND*(zTop(c1)+zTop(c2)) - dzVertGoal

            ! Once iteration has converged, we want 0.5_RKIND*(zBot(c1) + zBot(c2)) <= zBotEdge
            deltaZBot = (2.0_RKIND*zBotEdge - zBot(c2) - zBotNew(c1))
            if(deltaZBot < 0.0_RKIND) then
              zBotNew(c1) = zBotNew(c1) + deltaZBot
            end if
            deltaZBot = (2.0_RKIND*zBotEdge - zBot(c1) - zBotNew(c2))
            if(deltaZBot < 0.0_RKIND) then
              zBotNew(c2) = zBotNew(c2) + deltaZBot
            end if
          end do !iEdge

          where(smoothingMask(:) == 1)
            zBot(:) = (1.0_RKIND - weight)*zBot(:) + weight*zBotNew(:)
          end where

          call mpas_deallocate_scratch_field(zBotNewField, .true.)

          block_ptr => block_ptr % next
        end do !block_ptr

        ! do halo update on zBot
        call mpas_dmpar_exch_halo_field(zBotField)

      end do !iterIndex

    end subroutine constrain_rx1_layer

!***********************************************************************

end module ocn_init_vertical_grids


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker et ts=3 tw=132
