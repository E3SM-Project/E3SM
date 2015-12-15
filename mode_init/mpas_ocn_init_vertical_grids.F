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

   public :: ocn_generate_vertical_grid, &
             ocn_compute_layerThickness_zMid_from_bottomDepth, &
             ocn_alter_bottomDepth_for_pbcs, &
             ocn_compute_Haney_number

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
      else if ( trim(gridType) == '100layerACMEv1' ) then
         call ocn_generate_100layerACMEv1_vertical_grid(interfaceLocations)
      else if ( trim(gridType) == '1dCVTgenerator' ) then
         if (.not. present(configPool)) then
            call mpas_dmpar_global_abort("ERROR: requesting a 1d CVT vertical grid generation without passing the corresponding parameters. Exiting...")
         else
            call ocn_generate_1dCVT_vertical_grid(configPool, interfaceLocations)
         end if
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
   subroutine ocn_generate_uniform_vertical_grid(interfaceLocations)!{{{
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: interfaceLocations

      real (kind=RKIND) :: layerSpacing
      integer :: nInterfaces, iInterface

      write(stderrUnit,* ) ' ---- Generating uniform vertical grid ---- '

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
   !> \brief   60 layer PHC vertical grid generator
   !> \author  Doug Jacobsen
   !> \date    03/20/2015
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
         call mpas_dmpar_global_abort("ERROR: Vertical grid must have 60 layers to apply 60 Layer PHC grid. Exiting...")
      end if

      interfaceLocations(1) = 0.0_RKIND
      interfaceLocations(2) = 500_RKIND
      interfaceLocations(3) = 1500_RKIND
      interfaceLocations(4) = 2500_RKIND
      interfaceLocations(5) = 3500_RKIND
      interfaceLocations(6) = 4500_RKIND
      interfaceLocations(7) = 5500_RKIND
      interfaceLocations(8) = 6500_RKIND
      interfaceLocations(9) = 7500_RKIND
      interfaceLocations(10) = 8500_RKIND
      interfaceLocations(11) = 9500_RKIND
      interfaceLocations(12) = 10500_RKIND
      interfaceLocations(13) = 11500_RKIND
      interfaceLocations(14) = 12500_RKIND
      interfaceLocations(15) = 13500_RKIND
      interfaceLocations(16) = 14500_RKIND
      interfaceLocations(17) = 15500_RKIND
      interfaceLocations(18) = 16509.83984375_RKIND
      interfaceLocations(19) = 17547.904296875_RKIND
      interfaceLocations(20) = 18629.125_RKIND
      interfaceLocations(21) = 19766.025390625_RKIND
      interfaceLocations(22) = 20971.134765625_RKIND
      interfaceLocations(23) = 22257.826171875_RKIND
      interfaceLocations(24) = 23640.880859375_RKIND
      interfaceLocations(25) = 25137.013671875_RKIND
      interfaceLocations(26) = 26765.416015625_RKIND
      interfaceLocations(27) = 28548.361328125_RKIND
      interfaceLocations(28) = 30511.91796875_RKIND
      interfaceLocations(29) = 32686.794921875_RKIND
      interfaceLocations(30) = 35109.34375_RKIND
      interfaceLocations(31) = 37822.75390625_RKIND
      interfaceLocations(32) = 40878.4609375_RKIND
      interfaceLocations(33) = 44337.765625_RKIND
      interfaceLocations(34) = 48273.66796875_RKIND
      interfaceLocations(35) = 52772.796875_RKIND
      interfaceLocations(36) = 57937.28515625_RKIND
      interfaceLocations(37) = 63886.2578125_RKIND
      interfaceLocations(38) = 70756.328125_RKIND
      interfaceLocations(39) = 78700.25_RKIND
      interfaceLocations(40) = 87882.5234375_RKIND
      interfaceLocations(41) = 98470.5859375_RKIND
      interfaceLocations(42) = 110620.421875_RKIND
      interfaceLocations(43) = 124456.6953125_RKIND
      interfaceLocations(44) = 140049.71875_RKIND
      interfaceLocations(45) = 157394.640625_RKIND
      interfaceLocations(46) = 176400.328125_RKIND
      interfaceLocations(47) = 196894.421875_RKIND
      interfaceLocations(48) = 218645.65625_RKIND
      interfaceLocations(49) = 241397.15625_RKIND
      interfaceLocations(50) = 264900.125_RKIND
      interfaceLocations(51) = 288938.46875_RKIND
      interfaceLocations(52) = 313340.46875_RKIND
      interfaceLocations(53) = 337979.375_RKIND
      interfaceLocations(54) = 362767.0625_RKIND
      interfaceLocations(55) = 387645.21875_RKIND
      interfaceLocations(56) = 412576.84375_RKIND
      interfaceLocations(57) = 437539.28125_RKIND
      interfaceLocations(58) = 462519.0625_RKIND
      interfaceLocations(59) = 487508.375_RKIND
      interfaceLocations(60) = 512502.84375_RKIND
      interfaceLocations(61) = 537500_RKIND

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
         call mpas_dmpar_global_abort("ERROR: Vertical grid must have 60 layers to apply 60 Layer PHC grid. Exiting...")
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
   !  routine ocn_generate_100layerACMEv1_vertical_grid
   !
   !> \brief   100 vertical layer vertical grid generator for ACME v1
   !> \author  Todd Ringler
   !> \date    04/23/2015
   !> \details
   !>  This routine generates a 100 layer grid 
   !
   !-----------------------------------------------------------------------
   subroutine ocn_generate_100layerACMEv1_vertical_grid(interfaceLocations)!{{{
      implicit none

      real (kind=RKIND), dimension(:), intent(out) :: interfaceLocations

      real (kind=RKIND) :: maxInterfaceLocation
      integer :: nInterfaces, iInterface

      nInterfaces = size(interfaceLocations, dim=1)

      if ( nInterfaces /= 101 ) then
         call mpas_dmpar_global_abort("ERROR: Vertical grid must have 100 layers to apply 100 Layer PHC grid. Exiting...")
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

   end subroutine ocn_generate_100layerACMEv1_vertical_grid!}}}


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

   subroutine ocn_generate_1dCVT_vertical_grid(configPool, interfaceLocations)!{{{

      type (mpas_pool_type), intent(in) :: configPool
      real (kind=RKIND), dimension(:), intent(out) :: interfaceLocations

      integer :: k
      integer :: nInterfaces, nVertLevels
      real (kind=RKIND) :: stretch1
      real (kind=RKIND) :: stretch2
      real (kind=RKIND) :: dzSeed

      real (kind=RKIND) :: stretch
      real (kind=RKIND) :: dz
      real (kind=RKIND) :: maxInterfaceLocation

      real (kind=RKIND), pointer :: config_1dCVTgenerator_stretch1
      real (kind=RKIND), pointer :: config_1dCVTgenerator_stretch2
      real (kind=RKIND), pointer :: config_1dCVTgenerator_dzSeed

      call mpas_pool_get_config(configPool, 'config_1dCVTgenerator_stretch1', config_1dCVTgenerator_stretch1)
      call mpas_pool_get_config(configPool, 'config_1dCVTgenerator_stretch2', config_1dCVTgenerator_stretch2)
      call mpas_pool_get_config(configPool, 'config_1dCVTgenerator_dzSeed', config_1dCVTgenerator_dzSeed)

      stretch1 = config_1dCVTgenerator_stretch1
      stretch2 = config_1dCVTgenerator_stretch2
      dzSeed = config_1dCVTgenerator_dzSeed

      nInterfaces = size(interfaceLocations, dim=1)
      nVertLevels = nInterfaces - 1

      ! compute profile starting at top and stretch dz as we move down
      dz = dzSeed
      interfaceLocations(1) = 0.0_RKIND
      interfaceLocations(2) = dz
      do k=2,nVertLevels
         stretch = stretch1 + (stretch2-stretch1)*k/nVertLevels
         dz = stretch*dz
         interfaceLocations(k+1) = interfaceLocations(k) + dz
      enddo
      
      ! normalize so that positions span 0 to 1
      maxInterfaceLocation = maxval(interfaceLocations)
      interfaceLocations(:) = interfaceLocations(:) / maxInterfaceLocation
      
   end subroutine ocn_generate_1dCVT_vertical_grid!}}}


!***********************************************************************
!
!  routine ocn_compute_layerThickness_zMid_from_bottomDepth
!
!> \brief   Compute auxiliary z-variables from bottomDepth
!> \author  Mark Petersen
!> \date    10/17/2015
!> \details 
!>  This routine computes auxiliary z-variables from bottomDepth
!
!-----------------------------------------------------------------------

    subroutine ocn_compute_layerThickness_zMid_from_bottomDepth(layerThickness,zMid,refBottomDepth,bottomDepth, &
                                                                maxLevelCell,nVertLevels,iErr,restingThickness,ssh)!{{{
      real (kind=RKIND), dimension(nVertLevels), intent(out) :: layerThickness, zMid
      real (kind=RKIND), dimension(nVertLevels), intent(in) :: refBottomDepth
      real (kind=RKIND), intent(in) :: bottomDepth
      integer, intent(in) :: maxLevelCell, nVertLevels
      integer, intent(out) :: iErr
      real (kind=RKIND), dimension(nVertLevels), intent(out), optional :: restingThickness
      real (kind=RKIND), intent(in), optional :: ssh

      integer :: k
      real (kind=RKIND) :: layerStretch, zTop

      iErr = 0

      layerThickness(:) = 0.0_RKIND
      zMid(:) = 0.0_RKIND

      if(present(ssh) .and. .not. present(restingThickness)) then
         write (stderrUnit,*) ' Error: ssh present but restingThickness not present in ocn_compute_layerThickness_zMid_from_bottomDepth'
         iErr = 1
         return
      end if

      if (maxLevelCell<=0) return

      ! first, compute the resting layer thickness (same as layer thickness if ssh not present)
      if (maxLevelCell==1) then
         layerThickness(1) = bottomDepth
      else
         layerThickness(1) = refBottomDepth(1)

         do k = 2, maxLevelCell-1
            layerThickness(k) = refBottomDepth(k) - refBottomDepth(k-1)
         end do

         k = maxLevelCell
         layerThickness(k) = bottomDepth - refBottomDepth(k-1)

      endif

      zTop = 0.0_RKIND
      ! copy to layerThickness to restingThickness
      if (present(restingThickness)) then
         restingThickness(:) = layerThickness(:)
         ! stretch layers if ssh is present
         if(present(ssh)) then
            layerStretch = (ssh + bottomDepth)/bottomDepth
            zTop = ssh
            do k=1,maxLevelCell
               layerThickness(k) = layerStretch*restingThickness(k)
            end do
         end if
      end if

      ! compute zMid based on the layer thickness
      do k = 1, maxLevelCell
         zMid(k) = zTop - 0.5_RKIND*layerThickness(k)
         zTop = zTop - layerThickness(k)
      end do

    end subroutine ocn_compute_layerThickness_zMid_from_bottomDepth  !}}}


!***********************************************************************
!
!  routine ocn_alter_bottomDepth_for_pbcs
!
!> \brief   Alter bottom depth for partial bottom cells
!> \author  Mark Petersen
!> \date    10/19/2015
!> \details 
!>  This routine alters the bottom depth in a single column based on pbc settings
!
!-----------------------------------------------------------------------
    subroutine ocn_alter_bottomDepth_for_pbcs(bottomDepth, refBottomDepth, maxLevelCell, iErr)

      real (kind=RKIND), intent(inout) :: bottomDepth
      integer, intent(inout) :: maxLevelCell
      real (kind=RKIND), dimension(maxLevelCell), intent(in) :: refBottomDepth
      integer, intent(out) :: iErr
      integer :: k

      logical, pointer :: config_alter_ICs_for_pbcs
      real (kind=RKIND) :: minBottomDepth, minBottomDepthMid
      real (kind=RKIND), pointer :: config_min_pbc_fraction
      character (len=StrKIND), pointer :: config_pbc_alteration_type
      call mpas_pool_get_config(ocnConfigs, 'config_alter_ICs_for_pbcs', config_alter_ICs_for_pbcs)
      call mpas_pool_get_config(ocnConfigs, 'config_pbc_alteration_type', config_pbc_alteration_type)
      call mpas_pool_get_config(ocnConfigs, 'config_min_pbc_fraction', config_min_pbc_fraction)

      iErr = 0

      if (maxLevelCell > 1) then
         if (config_alter_ICs_for_pbcs) then

            if (config_pbc_alteration_type .eq. 'partial_cell') then
               ! Change value of maxLevelCell for partial bottom cells
               k = maxLevelCell
               minBottomDepth = refBottomDepth(k) - (1.0-config_min_pbc_fraction)*(refBottomDepth(k) - refBottomDepth(k-1))
               minBottomDepthMid = 0.5*(minBottomDepth + refBottomDepth(k-1))
               if (bottomDepth .lt. minBottomDepthMid) then
                  ! Round up to cell above
                  maxLevelCell = maxLevelCell - 1
                  bottomDepth = refBottomDepth(maxLevelCell)
               else if (bottomDepth .lt. minBottomDepth) then
                  ! Round down cell to the min_pbc_fraction.
                  bottomDepth = minBottomDepth
               end if
            elseif (config_pbc_alteration_type .eq. 'full_cell') then
               bottomDepth = refBottomDepth(maxLevelCell)
            else
               write (stderrUnit,*) ' Error: Incorrect choice of config_pbc_alteration_type: ', config_pbc_alteration_type
               iErr = 1
            endif
         endif
      endif

    end subroutine ocn_alter_bottomDepth_for_pbcs

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
    subroutine ocn_compute_Haney_number(dminfo, meshPool, diagnosticsPool, iErr)

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(inout) :: diagnosticsPool
      integer, intent(out) :: iErr
      type (dm_info), intent(in) :: dminfo

      real (kind=RKIND), dimension(:,:), pointer :: zMid
      real (kind=RKIND), dimension(:,:), pointer :: rx1Edge, rx1Cell
      real (kind=RKIND), dimension(:), pointer :: rx1MaxEdge, rx1MaxCell
      real (kind=RKIND), pointer :: globalMaxRx1Edge

      integer, pointer :: nCells, nVertLevels, nEdges
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer :: iEdge, c1, c2, k, maxLevelEdge

      real (kind=RKIND) :: dzVert1, dzVert2, dzEdgeK, dzEdgeKp1, rx1, localMaxRx1Edge

      iErr = 0

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)

      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

      call mpas_pool_get_array(diagnosticsPool, 'rx1Edge', rx1Edge)
      call mpas_pool_get_array(diagnosticsPool, 'rx1Cell', rx1Cell)
      call mpas_pool_get_array(diagnosticsPool, 'rx1MaxEdge', rx1MaxEdge)
      call mpas_pool_get_array(diagnosticsPool, 'rx1MaxCell', rx1MaxCell)
      call mpas_pool_get_array(diagnosticsPool, 'globalMaxRx1Edge', globalMaxRx1Edge)

      rx1Edge(:,:) = 0.0_RKIND
      rx1Cell(:,:) = 0.0_RKIND
      rx1MaxEdge(:) = 0.0_RKIND
      rx1MaxCell(:) = 0.0_RKIND
      do iEdge = 1,nEdges
        c1 = cellsOnEdge(1,iEdge)
        c2 = cellsOnEdge(2,iEdge)
        maxLevelEdge = min(maxLevelCell(c1), maxLevelCell(c2))
        do k = 1,maxLevelEdge-1
          dzVert1 = zMid(k,c1)-zMid(k+1,c1)
          dzVert2 = zMid(k,c2)-zMid(k+1,c2)
          dzEdgeK = zMid(k,c2)-zMid(k,c1)
          dzEdgeKp1 = zMid(k+1,c2)-zMid(k+1,c1)

          rx1 = abs(dzEdgeK+dzEdgeKp1)/(dzVert1+dzVert2)

          rx1Edge(k,iEdge) = rx1
          rx1Cell(k,c1) = max(rx1Cell(k,c1),rx1)
          rx1Cell(k,c2) = max(rx1Cell(k,c2),rx1)

          rx1MaxEdge(iEdge) = max(rx1MaxEdge(iEdge),rx1)
          rx1MaxCell(c2) = max(rx1MaxCell(c2),rx1)
          rx1MaxCell(c1) = max(rx1MaxCell(c1),rx1)
        end do
      end do

      localMaxRx1Edge = maxval(rx1MaxEdge)
      call mpas_dmpar_max_real(dminfo, localMaxRx1Edge, globalMaxRx1Edge)
      write (stdoutUnit,'(a, es10.2)') ' global max of rx1Edge:', globalMaxRx1Edge

    end subroutine ocn_compute_Haney_number

!***********************************************************************

end module ocn_init_vertical_grids

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
