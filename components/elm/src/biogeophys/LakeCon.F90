module LakeCon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing constants and parameters for the Lake code 
  ! (CLM4-LISSS, documented in Subin et al. 2011, JAMES)
  ! Also contains time constant variables for Lake code
  ! Created by Zack Subin, 2011
  !
  ! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use decompMod    , only : bounds_type
  use clm_varctl   , only : iulog
  use spmdMod      , only : masterproc
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: LakeConInit
  !-----------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Lake Model non-tuneable constants 
  !------------------------------------------------------------------

  ! temperature of maximum water density (K)
  ! This is from Hostetler and Bartlein (1990); more updated sources suggest 277.13 K.
  real(r8), parameter :: tdmax = 277._r8   

  !------------------------------------------------------------------
  ! Lake Model tuneable constants 
  !------------------------------------------------------------------

  ! lake emissivity. This is used for both frozen and unfrozen lakes. 
  ! This is pulled in from CLM4 and the reference is unclear.
  real(r8), parameter :: emg_lake = 0.97_r8

  ! The fraction of the visible (e.g. vis not nir from atm) sunlight
  ! absorbed in ~1 m of water (the surface layer za_lake).
  ! This is roughly the fraction over 700 nm but may depend on the details
  ! of atmospheric radiative transfer. As long as NIR = 700 nm and up, this can be zero.
  real(r8) :: betavis = 0.0_r8            

  ! Momentum Roughness length over frozen lakes without snow  (m)
  ! Typical value found in the literature, and consistent with Mironov expressions.
  ! See e.g. Morris EM 1989, Andreas EL 1987, Guest & Davidson 1991 (as cited in Vavrus 1996)
  real(r8), parameter :: z0frzlake = 0.001_r8  

  ! Base of surface light absorption layer for lakes (m)
  real(r8), parameter :: za_lake = 0.6_r8           

  ! For calculating prognostic roughness length
  real(r8), parameter :: cur0    = 0.01_r8  ! min. Charnock parameter
  real(r8), parameter :: cus     = 0.1_r8   ! empirical constant for roughness under smooth flow
  real(r8), parameter :: curm    = 0.1_r8   ! maximum Charnock parameter

  ! The following will be set in initLake based on namelists. !TODO - fix this commend
  real(r8)            :: fcrit              ! critical dimensionless fetch for Charnock parameter.
  real(r8)            :: minz0lake          ! (m) Minimum allowed roughness length for unfrozen lakes.

  ! For calculating enhanced diffusivity
  real(r8), parameter :: n2min = 7.5e-5_r8 ! (s^-2) (yields diffusivity about 6 times km) ! Fang & Stefan 1996

  ! Note, this will be adjusted in initLake if the timestep is not 1800 s.
  ! Lake top numerics can oscillate with 0.01m top layer and 1800 s timestep.
  ! The problem is that the surface flux is fixed during the calculation of the top
  ! layer temperature in the diffusion and not corrected for the tendency of the top layer.
  ! This thickness will be added to all minimum and maximum snow layer thicknesses compared to that used over non-lakes.
  ! Analysis of the CFL condition suggests that the minimum snow layer thickness for 1800 s needs
  ! to be at least ~1.2 cm for the bulk snow values of conductivity and heat capacity
  ! and as much as 2.3 cm for pure ice.
  ! Alternatively, a check could be done in LakeTemperature in case
  ! t_grnd(c) - t_soisno(c,snl(c)+1) changed sign after the Crank-Nicholson step.
  ! Such an approach, while perhaps allowing additional snow layer resolution, has not been tested.
  ! The approach used over non-lakes is to have a first-order surface flux correction.
  ! We choose not to do that here because t_grnd can vary independently of the top model
  ! layer temperature, while it is fixed to the top layer temperature if tbot > tfrz and
  ! the lake is frozen, or if there is an unstable density gradient in the top unfrozen lake layer.
  real(r8)            :: lsadz = 0.03_r8 ! m

  !! The following will be set in initLake based on namelists.
  real(r8)            :: pudz          ! (m) Optional minimum total ice thickness required to allow lake puddling.
  ! Currently used for sensitivity tests only.
  real(r8)            :: depthcrit     ! (m) Depth beneath which to increase mixing. See discussion in Subin et al. 2011
  real(r8)            :: mixfact       ! Mixing increase factor.

  !!!!!!!!!!!
  ! Namelists (some of these have not been extensively tested and are hardwired to default values currently).
  !!!!!!!!!!!

  ! used in LakeFluxes
  ! true => use old fcrit & minz0 as per Subin et al 2011 form
  ! See initLakeMod for details. Difference is very small for
  ! small lakes and negligible for large lakes. Currently hardwired off.
  logical,  public :: lake_use_old_fcrit_minz0 = .false. 

  ! used in LakeTemperature
  ! Increase mixing by a large factor for deep lakes
  ! Crude but enhanced performance at all 4 deep lakes tested.
  ! See Subin et al 2011 (JAMES) for details

  ! (m) minimum lake depth to invoke deepmixing
  real(r8), public :: deepmixing_depthcrit = 25._r8     

  ! factor to increase mixing by
  real(r8), public :: deepmixing_mixfact   = 10._r8     

  ! true => Suppress enhanced diffusion. Small differences.
  ! Currently hardwired .false.
  ! See Subin et al 2011 for details.
  ! Enhanced diffusion is intended for under ice and at large depths.
  ! It is a much smaller change on its own than the "deepmixing"
  ! above, but it increases the effect of deepmixing under ice and for large depths.
  logical,  public :: lake_no_ed = .false.              

  ! puddling (not extensively tested and currently hardwired off)
  ! used in LakeTemperature and SurfaceAlbedo

  ! true => suppress convection when greater than minimum amount
  ! of ice is present. This also effectively sets lake_no_melt_icealb.
  logical,  public :: lakepuddling = .false.            

  ! (m) minimum amount of total ice nominal thickness before
  ! convection is suppressed
  real(r8), public :: lake_puddle_thick = 0.2_r8        
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine LakeConInit()
    !
    ! !DESCRIPTION:
    ! Initialize time invariant variables for S Lake code 
    !------------------------------------------------------------------------

    if (masterproc) write (iulog,*) 'Attempting to initialize time invariant variables for lakes'

    ! Set LakeCon constants according to namelist fields
    if (lake_use_old_fcrit_minz0) then
       ! critical dimensionless fetch for Charnock parameter. From Vickers & Mahrt 1997
       ! but converted to use u instead of u* (Form used in Subin et al. 2011)
       fcrit   = 22._r8 

       ! (m) Minimum allowed roughness length for unfrozen lakes.
       ! (Used in Subin et al. 2011)
       minz0lake = 1.e-5_r8        
    else
       ! Vickers & Mahrt 1997
       fcrit   = 100._r8  

       ! (m) Minimum allowed roughness length for unfrozen lakes.
       ! Now set low so it is only to avoid floating point exceptions.
       minz0lake = 1.e-10_r8       
    end if

    if (lakepuddling) then
       ! (m) Minimum total ice thickness required to allow lake puddling. Default is 0.2m.
       ! This option has not been extensively tested.
       ! This option turns on lake_no_melt_icealb, as the decrease in albedo will be based
       ! on whether there is water over nice, not purely a function of ice top temperature.
       pudz = lake_puddle_thick    
    end if

    ! (m) Depth beneath which to increase mixing. See discussion in Subin et al. 2011
    depthcrit = deepmixing_depthcrit 

    ! Mixing increase factor. ! Defaults are 25 m, increase by 10.
    ! Note some other namelists will be used directly in lake physics during model integration.
    mixfact = deepmixing_mixfact 

    if (masterproc) write (iulog,*) 'Successfully initialized time invariant variables for lakes'

  end subroutine LakeConInit

end module LakeCon
