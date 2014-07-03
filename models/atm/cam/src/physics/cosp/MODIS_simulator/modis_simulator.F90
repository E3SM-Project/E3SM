! (c) 2009-2010, Regents of the Unversity of Colorado
!   Author: Robert Pincus, Cooperative Institute for Research in the Environmental Sciences
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!

!
! History:
!   May 2009 - Robert Pincus - Initial version
!   June 2009 - Steve Platnick and Robert Pincus - Simple radiative transfer for size retrievals
!   August 2009 - Robert Pincus - Consistency and bug fixes suggested by Rick Hemler (GFDL) 
!   November 2009 - Robert Pincus - Bux fixes and speed-ups after experience with Rick Hemler using AM2 (GFDL) 
!   January 2010 - Robert Pincus - Added high, middle, low cloud fractions 
!

!
! Notes on using the MODIS simulator: 
!  *) You may provide either layer-by-layer values of optical thickness at 0.67 and 2.1 microns, or 
!     optical thickness at 0.67 microns and ice- and liquid-water contents (in consistent units of 
!     your choosing)
!  *) Required input also includes the optical thickness and cloud top pressure 
!     derived from the ISCCP simulator run with parameter top_height = 1. 
!  *) Cloud particle sizes are specified as radii, measured in meters, though within the module we 
!     use units of microns. Where particle sizes are outside the bounds used in the MODIS retrieval
!     libraries (parameters re_water_min, re_ice_min, etc.) the simulator returns missing values (re_fill)

!
! When error conditions are encountered this code calls the function complain_and_die, supplied at the 
!   bottom of this module. Users probably want to replace this with something more graceful. 
!
module mod_modis_sim
  USE MOD_COSP_CONSTANTS, only: R_UNDEF
  USE abortutils,         only: complain_and_die => endrun
  implicit none
  ! ------------------------------
  ! Algorithmic parameters
  !
 
  real, parameter :: ice_density          = 0.93               ! liquid density is 1.  
  !
  ! Retrieval parameters
  !
  real, parameter :: min_OpticalThickness = 0.3,             & ! Minimum detectable optical thickness
                     CO2Slicing_PressureLimit = 700. * 100., & ! Cloud with higher pressures use thermal methods, units Pa
                     CO2Slicing_TauLimit = 1.,               & ! How deep into the cloud does CO2 slicing see? 
                     phase_TauLimit      = 1.,               & ! How deep into the cloud does the phase detection see?
                     size_TauLimit       = 2.,               & ! Depth of the re retreivals
                     phaseDiscrimination_Threshold = 0.7       ! What fraction of total extincton needs to be 
                                                               !  in a single category to make phase discrim. work? 
  real,    parameter :: re_fill= -999.
  integer, parameter :: phaseIsNone = 0, phaseIsLiquid = 1, phaseIsIce = 2, phaseIsUndetermined = 3
  
  logical, parameter :: useSimpleReScheme = .false. 
  !
  ! These are the limits of the libraries for the MODIS collection 5 algorithms 
  !   They are also the limits used in the fits for g and w0
  !
  real,    parameter :: re_water_min= 4., re_water_max= 30., re_ice_min= 5., re_ice_max= 90.
  integer, parameter :: num_trial_res = 15             ! increase to make the linear pseudo-retrieval of size more accurate
  logical, parameter :: use_two_re_iterations = .false. ! do two retrieval iterations? 
  
  !
  ! Precompute near-IR optical params vs size for retrieval scheme
  !
  integer, private :: i

! When Fortran 2003 is supported (e.g. Lahey support is dropped), "real" can once again be used to convert num_trial_res-1.
  real, dimension(num_trial_res), parameter :: & 
        trial_re_w = re_water_min + (re_water_max - re_water_min)/(num_trial_res-1.) * (/ (i - 1, i = 1, num_trial_res) /), &
        trial_re_i = re_ice_min   + (re_ice_max -   re_ice_min)/(num_trial_res-1.) * (/ (i - 1, i = 1, num_trial_res) /)
  
  ! Can't initialze these during compilation, but do in before looping columns in retrievals
  real, dimension(num_trial_res) ::  g_w, g_i, w0_w, w0_i
  ! ------------------------------
  ! Bin boundaries for the joint optical thickness/cloud top pressure histogram
  !
  integer, parameter :: numTauHistogramBins = 6, numPressureHistogramBins = 7

  real, private :: dummy_real 
  real, dimension(numTauHistogramBins + 1),      parameter :: &
    tauHistogramBoundaries = (/ min_OpticalThickness, 1.3, 3.6, 9.4, 23., 60., 100000. /) 
  real, dimension(numPressureHistogramBins + 1), parameter :: & ! Units Pa 
    pressureHistogramBoundaries = (/ 0., 180., 310., 440., 560., 680., 800., 1000. /) * 100. 
  real, parameter :: highCloudPressureLimit = 440. * 100., lowCloudPressureLimit = 680.  * 100.
  !
  ! For output - nominal bin centers and  bin boundaries. On output pressure bins are highest to lowest. 
  !
  integer, private :: k, l

  !! +jek the next line of code is causing the lehey compiler to issue a warning but it compiles anyways.
  !! "Element types in array constructor must be the same."
  real, parameter, dimension(2, numTauHistogramBins) ::   &
    nominalTauHistogramBoundaries =                       &
        reshape(source = (/ tauHistogramBoundaries(1),    &
                            ((tauHistogramBoundaries(k), l = 1, 2), k = 2, numTauHistogramBins), &
                            100000. /),                   &
                shape = (/2,  numTauHistogramBins /) )
  real, parameter, dimension(numTauHistogramBins) ::                    &
    nominalTauHistogramCenters = (nominalTauHistogramBoundaries(1, :) + &
                                  nominalTauHistogramBoundaries(2, :) ) / 2.
  
  real, parameter, dimension(2, numPressureHistogramBins) :: &
    nominalPressureHistogramBoundaries =                     &
        reshape(source = (/ 100000.,                         &
                            ((pressureHistogramBoundaries(k), l = 1, 2), k = numPressureHistogramBins, 2, -1), &
                            0.  /), &
                shape = (/2,  numPressureHistogramBins /) )
  real, parameter, dimension(numPressureHistogramBins) ::                         &
    nominalPressureHistogramCenters = (nominalPressureHistogramBoundaries(1, :) + &
                                       nominalPressureHistogramBoundaries(2, :) ) / 2.
  ! ------------------------------
  ! There are two ways to call the MODIS simulator: 
  !  1) Provide total optical thickness and liquid/ice water content and we'll partition tau in 
  !     subroutine modis_L2_simulator_oneTau, or 
  !  2) Provide ice and liquid optical depths in each layer
  !
  interface modis_L2_simulator
    module procedure modis_L2_simulator_oneTau, modis_L2_simulator_twoTaus
  end interface 
contains
  !------------------------------------------------------------------------------------------------
  ! MODIS simulator using specified liquid and ice optical thickness in each layer 
  !
  !   Note: this simulator operates on all points; to match MODIS itself night-time 
  !     points should be excluded
  !
  !   Note: the simulator requires as input the optical thickness and cloud top pressure 
  !     derived from the ISCCP simulator run with parameter top_height = 1. 
  !     If cloud top pressure is higher than about 700 mb, MODIS can't use CO2 slicing 
  !     and reverts to a thermal algorithm much like ISCCP's. Rather than replicate that 
  !     alogrithm in this simulator we simply report the values from the ISCCP simulator. 
  !
  subroutine modis_L2_simulator_twoTaus(                                       &
                                temp, pressureLayers, pressureLevels,          &
                                liquid_opticalThickness, ice_opticalThickness, snow_opticalThickness, &
                                waterSize, iceSize, snowSize,                  & 
                                isccpTau, isccpCloudTopPressure,               &
                                retrievedPhase, retrievedCloudTopPressure, retrievedTau, retrievedSize)

    ! Grid-mean quantities at layer centers, starting at the model top
    !   dimension nLayers
    real, dimension(:),    intent(in ) :: temp,           & ! Temperature, K
                                          pressureLayers, & ! Pressure, Pa
                                          pressureLevels    ! Pressure at layer edges, Pa (dimension nLayers + 1) 
    ! Sub-column quantities
    !   dimension  nSubcols, nLayers
    real, dimension(:, :), intent(in ) :: liquid_opticalThickness, & ! Layer optical thickness @ 0.67 microns due to liquid
                                          ice_opticalThickness,    & ! ditto, due to ice
                                          snow_opticalThickness      ! ditto for snow
    real, dimension(:, :), intent(in ) :: waterSize,        & ! Cloud drop effective radius, microns
                                          iceSize,          & ! Cloud ice effective radius, microns
                                          snowSize
                                          
    ! Cloud properties retrieved from ISCCP using top_height = 1
    !    dimension nSubcols
    real, dimension(:),    intent(in ) :: isccpTau, &           ! Column-integrated optical thickness 
                                          isccpCloudTopPressure ! ISCCP-retrieved cloud top pressure (Pa) 

    ! Properties retrieved by MODIS
    !   dimension nSubcols
    integer, dimension(:), intent(out) :: retrievedPhase               ! liquid/ice/other - integer, defined in module header
    real,    dimension(:), intent(out) :: retrievedCloudTopPressure, & ! units of pressureLayers
                                          retrievedTau,              & ! unitless
                                          retrievedSize                ! microns 
    ! ---------------------------------------------------
    ! Local variables
    logical, dimension(size(retrievedTau))                     :: cloudMask
    real,    dimension(size(waterSize, 1), size(waterSize, 2)) :: tauLiquidFraction, tauTotal
    real    :: integratedLiquidFraction
    integer :: i, nSubcols, nLevels

    ! ---------------------------------------------------
    nSubcols = size(liquid_opticalThickness, 1)
    nLevels  = size(liquid_opticalThickness, 2) 
 
    !
    ! Initial error checks 
    !   
    if(any((/ size(ice_opticalThickness, 1), size(waterSize, 1), size(iceSize, 1), &
              size(snow_opticalThickness, 1), size(snowSize, 1),        &
              size(isccpTau), size(isccpCloudTopPressure),              &
              size(retrievedPhase), size(retrievedCloudTopPressure),    &
              size(retrievedTau), size(retrievedSize) /) /= nSubcols )) &
       call complain_and_die("Failure in MODIS simulator: Differing number of subcolumns in one or more arrays") 
    
    if(any((/ size(ice_opticalThickness, 2), size(snow_opticalThickness, 2), & 
              size(waterSize, 2), size(iceSize, 2), size(snowSize, 2),       &
              size(temp), size(pressureLayers), size(pressureLevels)-1 /) /= nLevels )) &
       call complain_and_die("Failure in MODIS simulator: Differing number of levels in one or more arrays") 
       
    if(any( (/ any(temp <= 0.), any(pressureLayers <= 0.),  &
               any(liquid_opticalThickness < 0.),           &
               any(ice_opticalThickness < 0.),              &
               any(snow_opticalThickness < 0.),             &
               any(waterSize < 0.), any(iceSize < 0.), any(snowSize < 0.) /) )) &
       call complain_and_die("Failure in MODIS simulator: Input values out of bounds") 
             
    ! ---------------------------------------------------
    !
    ! Compute the total optical thickness and the proportion due to liquid in each cell
    !
    tauTotal(:, :) = liquid_opticalThickness(:, :) + ice_opticalThickness(:, :) + snow_opticalThickness(:, :)
    where(tauTotal > 0.) 
      tauLiquidFraction(:, :) = liquid_opticalThickness(:, :) / tauTotal(:, :) 
    elsewhere
      tauLiquidFraction(:, :) = 0. 
    end  where 
    
    !
    ! Optical depth retrieval 
    !   This is simply a sum over the optical thickness in each layer 
    !   It should agree with the ISCCP values after min values have been excluded 
    !
    retrievedTau(:) = sum(tauTotal(:, :), dim = 2)

    !
    ! Cloud detection - does optical thickness exceed detection threshold? 
    !
    cloudMask = retrievedTau(:) >= min_OpticalThickness
    
    !
    ! Initialize initial estimates for size retrievals
    !
    if(any(cloudMask) .and. .not. useSimpleReScheme) then 
      g_w(:)  = get_g_nir(  phaseIsLiquid, trial_re_w(:))
      w0_w(:) = get_ssa_nir(phaseIsLiquid, trial_re_w(:))
      g_i(:)  = get_g_nir(  phaseIsIce,    trial_re_i(:))
      w0_i(:) = get_ssa_nir(phaseIsIce,    trial_re_i(:))
    end if 
    
    do i = 1, nSubCols
      if(cloudMask(i)) then 
        !
        ! Cloud top pressure determination 
        !   MODIS uses CO2 slicing for clouds with tops above about 700 mb and thermal methods for clouds
        !   lower than that. 
        !  For CO2 slicing we report the optical-depth weighted pressure, integrating to a specified 
        !    optical depth
        ! This assumes linear variation in p between levels. Linear in ln(p) is probably better, 
        !   though we'd need to deal with the lowest pressure gracefully. 
        !
        retrievedCloudTopPressure(i) = cloud_top_pressure((/ 0., tauTotal(i, :) /), &
                                                          pressureLevels,           &
                                                          CO2Slicing_TauLimit)  
        
        
        !
        ! Phase determination - determine fraction of total tau that's liquid 
        ! When ice and water contribute about equally to the extinction we can't tell 
        !   what the phase is 
        !
        integratedLiquidFraction = weight_by_extinction(tauTotal(i, :),          &
                                                        tauLiquidFraction(i, :), &
                                                        phase_TauLimit)
        if(integratedLiquidFraction >= phaseDiscrimination_Threshold) then 
          retrievedPhase(i) = phaseIsLiquid
        else if (integratedLiquidFraction <= 1.- phaseDiscrimination_Threshold) then 
          retrievedPhase(i) = phaseIsIce
        else 
          retrievedPhase(i) = phaseIsUndetermined
        end if 
        
        !
        ! Size determination 
        !
        if(useSimpleReScheme) then 
          !   This is the extinction-weighted size considering only the phase we've chosen 
          !
          if(retrievedPhase(i) == phaseIsIce) then 
            retrievedSize(i) = weight_by_extinction(ice_opticalThickness(i, :),  &
                                                    iceSize(i, :), &
                                                    (1. - integratedLiquidFraction) * size_TauLimit)
  
          else if(retrievedPhase(i) == phaseIsLiquid) then 
            retrievedSize(i) = weight_by_extinction(liquid_opticalThickness(i, :), &
                                                    waterSize(i, :), &
                                                    integratedLiquidFraction * size_TauLimit)
  
          else
            retrievedSize(i) = 0. 
          end if 
        else
          retrievedSize(i) = 1.0e-06*retrieve_re(retrievedPhase(i), retrievedTau(i), &
               obs_Refl_nir = compute_nir_reflectance(liquid_opticalThickness(i, :), waterSize(i, :)*1.0e6, & 
               ice_opticalThickness(i, :),      iceSize(i, :)*1.0e6, &
               snow_opticalThickness(i, :),    snowSize(i, :)*1.0e6))
        end if 
      else 
        !
        ! Values when we don't think there's a cloud. 
        !
        retrievedCloudTopPressure(i) = R_UNDEF 
        retrievedPhase(i) = phaseIsNone
        retrievedSize(i) = R_UNDEF 
        retrievedTau(i) = R_UNDEF 
      end if
    end do
    where((retrievedSize(:) < 0.).and.(retrievedSize(:) /= R_UNDEF)) retrievedSize(:) = 1.0e-06*re_fill

    ! We use the ISCCP-derived CTP for low clouds, since the ISCCP simulator ICARUS 
    !   mimics what MODIS does to first order. 
    !   Of course, ISCCP cloud top pressures are in mb. 
    !   
    where(cloudMask(:) .and. retrievedCloudTopPressure(:) > CO2Slicing_PressureLimit) &
      retrievedCloudTopPressure(:) = isccpCloudTopPressure * 100. 
    
  end subroutine modis_L2_simulator_twoTaus
  !------------------------------------------------------------------------------------------------
  !
  ! MODIS simulator: provide a single optical thickness and the cloud ice and liquid contents; 
  !   we'll partition this into ice and liquid optical thickness and call the full MODIS simulator 
  ! 
  subroutine modis_L2_simulator_oneTau(                                         &
                                temp, pressureLayers, pressureLevels,           &
                                opticalThickness, cloudWater, cloudIce, cloudSnow, &
                                waterSize, iceSize, snowSize,                      & 
                                isccpTau, isccpCloudTopPressure,                &
                                retrievedPhase, retrievedCloudTopPressure, retrievedTau, retrievedSize)
    ! Grid-mean quantities at layer centers, 
    !   dimension nLayers
    real, dimension(:),    intent(in ) :: temp,           & ! Temperature, K
                                          pressureLayers, & ! Pressure, Pa
                                          pressureLevels    ! Pressure at layer edges, Pa (dimension nLayers + 1) 
    ! Sub-column quantities
    !   dimension nLayers, nSubcols
    real, dimension(:, :), intent(in ) :: opticalThickness, & ! Layer optical thickness @ 0.67 microns
                                          cloudWater,       & ! Cloud water content, arbitrary units
                                          cloudIce,         & ! Cloud water content, same units as cloudWater
                                          cloudSnow           ! Snow content, same units as cloudWater
    real, dimension(:, :), intent(in ) :: waterSize,        & ! Cloud drop effective radius, microns
                                          iceSize,          & ! Cloud ice effective radius, microns
                                          snowSize            ! Snow effective radius, microns 

    ! Cloud properties retrieved from ISCCP using top_height = 1
    !    dimension nSubcols
    
    real, dimension(:),    intent(in ) :: isccpTau, &           ! Column-integrated optical thickness 
                                          isccpCloudTopPressure ! ISCCP-retrieved cloud top pressure (Pa) 

    ! Properties retrieved by MODIS
    !   dimension nSubcols
    integer, dimension(:), intent(out) :: retrievedPhase               ! liquid/ice/other - integer
    real,    dimension(:), intent(out) :: retrievedCloudTopPressure, & ! units of pressureLayers
                                          retrievedTau,              & ! unitless
                                          retrievedSize                ! microns (or whatever units 
                                                                       !   waterSize and iceSize are supplied in)
    ! ---------------------------------------------------
    ! Local variables
    real, dimension(size(opticalThickness, 1), size(opticalThickness, 2)) :: & 
           liquid_opticalThickness, ice_opticalThickness, snow_opticalThickness, totalExtinction
    logical, dimension(size(opticalThickness, 1), size(opticalThickness, 2)) :: inconsist 
    
    ! ---------------------------------------------------
    
    ! 
    ! Geometic optics limit - tau as LWP/re  (proportional to LWC/re) 
    !
    !!+jek floating point exception here.... likely because dividing by zero.

!    totalExtinction = merge(cloudWater(:, :)/waterSize(:, :),                0., cloudWater(:, :) > 0.) + &
!                      merge(cloudIce  (:, :)/(ice_density * iceSize (:, :)), 0., cloudIce  (:, :) > 0.) + &
!                      merge(cloudSnow (:, :)/(ice_density * snowSize(:, :)), 0., cloudSnow (:, :) > 0.) 
    
   ! 
   ! RP - compute total exinction in steps to avoid FPEs in merge statements
   !   We assume that (waterSize > 0. .eqv. cloudWater > 0.) at all times 
   !
   where(waterSize(:, :) > 0.) 
      totalExtinction(:, :) = cloudWater(:, :)/waterSize(:, :)
    elsewhere
      totalExtinction(:, :) = 0. 
    end where 
    where(  iceSize(:, :) > 0.) &
      totalExtinction(:, :) = totalExtinction(:, :) + cloudIce  (:, :)/(ice_density * iceSize (:, :))
    where(  snowSize(:, :) > 0.) &
      totalExtinction(:, :) = totalExtinction(:, :) + cloudSnow (:, :)/(ice_density * snowSize(:, :))
    

!    liquid_opticalThickness(:, :) = merge(opticalThickness(:, :) * cloudWater(:, :) / &
!                                          (             waterSize(:, :) * totalExtinction(:, :)), &
!                                          0.,                                         &
!                                          cloudWater(:, :) > 0.)
!       ice_opticalThickness(:, :) = merge(opticalThickness(:, :) * cloudIce  (:, :) / &
!                                          (ice_density * iceSize (:, :) * totalExtinction(:, :)), &
!                                          0.,                                         &
!                                          cloudIce(:, :) > 0.) 
!      snow_opticalThickness(:, :) = merge(opticalThickness(:, :) * cloudSnow  (:, :) / &
!                                          (ice_density * snowSize (:, :) * totalExtinction(:, :)), &
!                                          0.,                                          &
!                                          cloudsnow(:, :) > 0.) 
    
   where((waterSize(:, :) > 0.) .and. (totalExtinction(:, :) > 0.))   !!+jek
     liquid_opticalThickness(:, :) = opticalThickness(:, :) * cloudWater(:, :) / &
                                          (              waterSize(:, :) * totalExtinction(:, :))
    elsewhere
      liquid_opticalThickness(:, :) = 0. 
    end where 

    where(  (iceSize(:, :) > 0.) .and. (totalExtinction(:, :) > 0.))   !!+jek
        ice_opticalThickness(:, :) = opticalThickness(:, :) * cloudIce  (:, :) / &
                                          (ice_density *   iceSize(:, :) * totalExtinction(:, :))
    elsewhere
       ice_opticalThickness(:, :) = 0. 
    end where 

    where( (snowSize(:, :) > 0.) .and. (totalExtinction(:, :) > 0.)) !!+jek
      snow_opticalThickness(:, :) = opticalThickness(:, :) * cloudSnow (:, :) / &
                                          (ice_density *  snowSize(:, :) * totalExtinction(:, :))
    elsewhere
      snow_opticalThickness(:, :) = 0. 
    end where 
 
    call modis_L2_simulator_twoTaus(temp, pressureLayers, pressureLevels,          &
                                    liquid_opticalThickness, ice_opticalThickness, snow_opticalThickness, &
                                    waterSize, iceSize, snowSize,                  & 
                                    isccpTau, isccpCloudTopPressure,               &
                                    retrievedPhase, retrievedCloudTopPressure, retrievedTau, retrievedSize)
                                
  end subroutine modis_L2_simulator_oneTau
  !------------------------------------------------------------------------------------------------
  subroutine modis_L3_simulator(phase, cloud_top_pressure, optical_thickness, particle_size,            &
       Cloud_Fraction_Total_Mean,       Cloud_Fraction_Water_Mean,       Cloud_Fraction_Ice_Mean,       &
       Cloud_Fraction_High_Mean,        Cloud_Fraction_Mid_Mean,         Cloud_Fraction_Low_Mean,       &
       Optical_Thickness_Total_Mean,    Optical_Thickness_Water_Mean,    Optical_Thickness_Ice_Mean,    &
       Optical_Thickness_Total_MeanLog10, Optical_Thickness_Water_MeanLog10, Optical_Thickness_Ice_MeanLog10, &
                                        Cloud_Particle_Size_Water_Mean,  Cloud_Particle_Size_Ice_Mean,  &
       Cloud_Top_Pressure_Total_Mean,                                                                   &
                                        Liquid_Water_Path_Mean,          Ice_Water_Path_Mean,           &    
       Optical_Thickness_vs_Cloud_Top_Pressure)
    !
    ! Inputs; dimension nPoints, nSubcols
    !
    integer, dimension(:, :),   intent(in)  :: phase
    real,    dimension(:, :),   intent(in)  :: cloud_top_pressure, optical_thickness, particle_size
    !
    ! Outputs; dimension nPoints
    !
    real,    dimension(:),      intent(out) :: &
       Cloud_Fraction_Total_Mean,       Cloud_Fraction_Water_Mean,       Cloud_Fraction_Ice_Mean,       &
       Cloud_Fraction_High_Mean,        Cloud_Fraction_Mid_Mean,         Cloud_Fraction_Low_Mean,       &
       Optical_Thickness_Total_Mean,    Optical_Thickness_Water_Mean,    Optical_Thickness_Ice_Mean,    &
       Optical_Thickness_Total_MeanLog10, Optical_Thickness_Water_MeanLog10, Optical_Thickness_Ice_MeanLog10, &
                                        Cloud_Particle_Size_Water_Mean,  Cloud_Particle_Size_Ice_Mean,  &
       Cloud_Top_Pressure_Total_Mean,                                                                   &
                                        Liquid_Water_Path_Mean,          Ice_Water_Path_Mean
    ! tau/ctp histogram; dimensions nPoints, numTauHistogramBins , numPressureHistogramBins 
    real,    dimension(:, :, :), intent(out) :: Optical_Thickness_vs_Cloud_Top_Pressure
    ! ---------------------------
    ! Local variables
    !
    real, parameter :: LWP_conversion = 2./3. * 1000. ! MKS units  
    integer :: i, j
    integer :: nPoints, nSubcols 
    logical, dimension(size(phase, 1), size(phase, 2)) :: &
      cloudMask, waterCloudMask, iceCloudMask, validRetrievalMask
    logical, dimension(size(phase, 1), size(phase, 2), numTauHistogramBins     ) :: tauMask
    logical, dimension(size(phase, 1), size(phase, 2), numPressureHistogramBins) :: pressureMask
!+jek fix from brian eaton
    real, dimension(size(phase, 1), size(phase, 2)) :: opt_thick_log

    ! ---------------------------
    
    nPoints  = size(phase, 1) 
    nSubcols = size(phase, 2) 
    !
    ! Array conformance checks
    !
    ! Ampersands below denote 132 character line limit. Do not write a longer line.
    if(any( (/ size(cloud_top_pressure, 1),     size(optical_thickness, 1),            size(particle_size, 1),                     &
         size(Cloud_Fraction_Total_Mean),       size(Cloud_Fraction_Water_Mean),       size(Cloud_Fraction_Ice_Mean),              &
         size(Cloud_Fraction_High_Mean),        size(Cloud_Fraction_Mid_Mean),         size(Cloud_Fraction_Low_Mean),              &
         size(Optical_Thickness_Total_Mean),    size(Optical_Thickness_Water_Mean),    size(Optical_Thickness_Ice_Mean),           &
         size(Optical_Thickness_Total_MeanLog10), size(Optical_Thickness_Water_MeanLog10), size(Optical_Thickness_Ice_MeanLog10),  &
                                                size(Cloud_Particle_Size_Water_Mean),  size(Cloud_Particle_Size_Ice_Mean),         &
         size(Cloud_Top_Pressure_Total_Mean),                                                                                      &
                                                size(Liquid_Water_Path_Mean),          size(Ice_Water_Path_Mean) /) /= nPoints))   &
      call complain_and_die("Failure in MODIS simulator: Some L3 arrays have wrong number of grid points") 
    if(any( (/ size(cloud_top_pressure, 2), size(optical_thickness, 2), size(particle_size, 2) /)  /= nSubcols)) &
      call complain_and_die("Failure in MODIS simulator: Some L3 arrays have wrong number of subcolumns") 
    
    
    !
    ! Include only those pixels with successful retrievals in the statistics 
    !
    validRetrievalMask(:, :) = particle_size(:, :) > 0.
    cloudMask      = phase(:, :) /= phaseIsNone   .and. validRetrievalMask(:, :)
    waterCloudMask = phase(:, :) == phaseIsLiquid .and. validRetrievalMask(:, :)
    iceCloudMask   = phase(:, :) == phaseIsIce    .and. validRetrievalMask(:, :)
    !
    ! Use these as pixel counts at first 
    !
    Cloud_Fraction_Total_Mean(:) = real(count(cloudMask,      dim = 2))
    Cloud_Fraction_Water_Mean(:) = real(count(waterCloudMask, dim = 2))
    Cloud_Fraction_Ice_Mean(:)   = real(count(iceCloudMask,   dim = 2))
    
    Cloud_Fraction_High_Mean(:) = real(count(cloudMask .and. cloud_top_pressure <= highCloudPressureLimit, dim = 2)) 
    Cloud_Fraction_Low_Mean(:)  = real(count(cloudMask .and. cloud_top_pressure >  lowCloudPressureLimit,  dim = 2)) 
    Cloud_Fraction_Mid_Mean(:)  = Cloud_Fraction_Total_Mean(:) - Cloud_Fraction_High_Mean(:) - Cloud_Fraction_Low_Mean(:)
    
    !
    ! Don't want to divide by 0, even though the sums will be 0 where the pixel counts are 0. 
    !
    where (Cloud_Fraction_Total_Mean == 0) Cloud_Fraction_Total_Mean = -1. 
    where (Cloud_Fraction_Water_Mean == 0) Cloud_Fraction_Water_Mean = -1.
    where (Cloud_Fraction_Ice_Mean   == 0) Cloud_Fraction_Ice_Mean   = -1.
    
    Optical_Thickness_Total_Mean = sum(optical_thickness, mask = cloudMask,      dim = 2) / Cloud_Fraction_Total_Mean(:) 
    Optical_Thickness_Water_Mean = sum(optical_thickness, mask = waterCloudMask, dim = 2) / Cloud_Fraction_Water_Mean(:)
    Optical_Thickness_Ice_Mean   = sum(optical_thickness, mask = iceCloudMask,   dim = 2) / Cloud_Fraction_Ice_Mean(:)

    opt_thick_log(:,:) = -1.e30

    where (cloudMask) opt_thick_log = log10(optical_thickness)
    Optical_Thickness_Total_MeanLog10 = sum(opt_thick_log, mask = cloudMask,      dim = 2) / Cloud_Fraction_Total_Mean(:)
    Optical_Thickness_Water_MeanLog10 = sum(opt_thick_log, mask = waterCloudMask, dim = 2) / Cloud_Fraction_Water_Mean(:)
    Optical_Thickness_Ice_MeanLog10   = sum(opt_thick_log, mask = iceCloudMask,   dim = 2) / Cloud_Fraction_Ice_Mean(:)

    Cloud_Particle_Size_Water_Mean = sum(particle_size, mask = waterCloudMask, dim = 2) / Cloud_Fraction_Water_Mean(:)
    Cloud_Particle_Size_Ice_Mean   = sum(particle_size, mask = iceCloudMask,   dim = 2) / Cloud_Fraction_Ice_Mean(:)
    
    Cloud_Top_Pressure_Total_Mean = sum(cloud_top_pressure, mask = cloudMask, dim = 2) / max(1, count(cloudMask, dim = 2))
    
    Liquid_Water_Path_Mean = LWP_conversion &
                             * sum(particle_size * optical_thickness, mask = waterCloudMask, dim = 2) &
                             / Cloud_Fraction_Water_Mean(:)
    Ice_Water_Path_Mean    = LWP_conversion * ice_density &
                             * sum(particle_size * optical_thickness, mask = iceCloudMask,   dim = 2) &
                             / Cloud_Fraction_Ice_Mean(:)

    !
    ! Normalize pixel counts to fraction
    !   The first three cloud fractions have been set to -1 in cloud-free areas, so set those places to 0.
    ! 
    Cloud_Fraction_Total_Mean(:) = max(0., Cloud_Fraction_Total_Mean(:)/nSubcols)
    Cloud_Fraction_Water_Mean(:) = max(0., Cloud_Fraction_Water_Mean(:)/nSubcols)
    Cloud_Fraction_Ice_Mean(:)   = max(0., Cloud_Fraction_Ice_Mean(:)  /nSubcols)
    
    Cloud_Fraction_High_Mean(:)  = Cloud_Fraction_High_Mean(:) /nSubcols
    Cloud_Fraction_Mid_Mean(:)   = Cloud_Fraction_Mid_Mean(:)  /nSubcols
    Cloud_Fraction_Low_Mean(:)   = Cloud_Fraction_Low_Mean(:)  /nSubcols
    
    ! ----
    ! Joint histogram 
    ! 
    do i = 1, numTauHistogramBins 
      where(cloudMask(:, :)) 
        tauMask(:, :, i) = optical_thickness(:, :) >= tauHistogramBoundaries(i) .and. &
                           optical_thickness(:, :) <  tauHistogramBoundaries(i+1)
      elsewhere
        tauMask(:, :, i) = .false.
      end where
    end do 

    do i = 1, numPressureHistogramBins 
      where(cloudMask(:, :)) 
        pressureMask(:, :, i) = cloud_top_pressure(:, :) >= pressureHistogramBoundaries(i) .and. &
                                cloud_top_pressure(:, :) <  pressureHistogramBoundaries(i+1)
      elsewhere
        pressureMask(:, :, i) = .false.
      end where
    end do 
    
    do i = 1, numPressureHistogramBins
      do j = 1, numTauHistogramBins
        Optical_Thickness_vs_Cloud_Top_Pressure(:, j, i) = & 
          real(count(tauMask(:, :, j) .and. pressureMask(:, :, i), dim = 2)) / real(nSubcols)
      end do 
    end do 
    
  end subroutine modis_L3_simulator
  !------------------------------------------------------------------------------------------------
  function cloud_top_pressure(tauIncrement, pressure, tauLimit) 
    real, dimension(:), intent(in) :: tauIncrement, pressure
    real,               intent(in) :: tauLimit
    real                           :: cloud_top_pressure
    !
    ! Find the extinction-weighted pressure. Assume that pressure varies linearly between 
    !   layers and use the trapezoidal rule.
    !
    
    real :: deltaX, totalTau, totalProduct
    integer :: i 
    
    totalTau = 0.; totalProduct = 0. 
    do i = 2, size(tauIncrement)
      if(totalTau + tauIncrement(i) > tauLimit) then 
        deltaX = tauLimit - totalTau
        totalTau = totalTau + deltaX
        !
        ! Result for trapezoidal rule when you take less than a full step
        !   tauIncrement is a layer-integrated value
        !
        totalProduct = totalProduct           &
                     + pressure(i-1) * deltaX &
                     + (pressure(i) - pressure(i-1)) * deltaX**2/(2. * tauIncrement(i)) 
      else
        totalTau =     totalTau     + tauIncrement(i) 
        totalProduct = totalProduct + tauIncrement(i) * (pressure(i) + pressure(i-1)) / 2.
      end if 
      if(totalTau >= tauLimit) exit
    end do 
    cloud_top_pressure = totalProduct/totalTau
  end function cloud_top_pressure
  !------------------------------------------------------------------------------------------------
  function weight_by_extinction(tauIncrement, f, tauLimit) 
    real, dimension(:), intent(in) :: tauIncrement, f
    real,               intent(in) :: tauLimit
    real                           :: weight_by_extinction
    !
    ! Find the extinction-weighted value of f(tau), assuming constant f within each layer
    !
    
    real    :: deltaX, totalTau, totalProduct
    integer :: i 
    
    totalTau = 0.; totalProduct = 0. 
    do i = 1, size(tauIncrement)
      if(totalTau + tauIncrement(i) > tauLimit) then 
        deltaX       = tauLimit - totalTau
        totalTau     = totalTau     + deltaX
        totalProduct = totalProduct + deltaX * f(i) 
      else
        totalTau     = totalTau     + tauIncrement(i) 
        totalProduct = totalProduct + tauIncrement(i) * f(i) 
      end if 
      if(totalTau >= tauLimit) exit
    end do 
    weight_by_extinction = totalProduct/totalTau
  end function weight_by_extinction
  !------------------------------------------------------------------------------------------------
  pure function compute_nir_reflectance(water_tau, water_size, ice_tau, ice_size, snow_tau, snow_size) 
    real, dimension(:), intent(in) :: water_tau, water_size, ice_tau, ice_size, snow_tau, snow_size
    real                           :: compute_nir_reflectance
    
    real, dimension(size(water_tau)) :: water_g, water_w0, ice_g, ice_w0, snow_g, snow_w0, &
                                        tau, g, w0
    !----------------------------------------
    water_g(:)  = get_g_nir(  phaseIsLiquid, water_size) 
    water_w0(:) = get_ssa_nir(phaseIsLiquid, water_size) 
    ice_g(:)    = get_g_nir(  phaseIsIce,    ice_size) 
    ice_w0(:)   = get_ssa_nir(phaseIsIce,    ice_size) 
    snow_g(:)   = get_g_nir(  phaseIsIce,    snow_size) 
    snow_w0(:)  = get_ssa_nir(phaseIsIce,    snow_size) 
    !
    ! Combine ice and water optical properties
    !
    g(:) = 0; w0(:) = 0. 
    tau(:) = ice_tau(:) + water_tau(:) + snow_tau(:)
    where (tau(:) > 0) 
      g(:)  = (water_tau(:) * water_g(:) + &
                 ice_tau(:) *   ice_g(:) + &
                snow_tau(:) *  snow_g(:) ) / & 
              tau(:) 
      w0(:) = (water_tau(:) * water_g(:) * water_w0(:) + &
                 ice_tau(:) *   ice_g(:) *   ice_w0(:) + &
                snow_tau(:) *  snow_g(:) *  snow_w0(:) ) / &
              (g(:) * tau(:))
    end where
    
    compute_nir_reflectance = compute_toa_reflectace(tau, g, w0)
  end function compute_nir_reflectance
  !------------------------------------------------------------------------------------------------
  ! Retreivals
  !------------------------------------------------------------------------------------------------
  elemental function retrieve_re (phase, tau, obs_Refl_nir)
      integer, intent(in) :: phase
      real,    intent(in) :: tau, obs_Refl_nir
      real                :: retrieve_re
      !
      ! Finds the re that produces the minimum mis-match between predicted and observed reflectance in 
      !   MODIS band 7 (near IR)
      ! Uses 
      !  fits for asymmetry parameter g(re) and single scattering albedo w0(re) based on MODIS tables 
      !  two-stream for layer reflectance and transmittance as a function of optical thickness tau, g, and w0
      !  adding-doubling for total reflectance 
      !  
      !
      !
      ! Local variables
      !
      real, parameter :: min_distance_to_boundary = 0.01
      real    :: re_min, re_max, delta_re
      integer :: i 
      
      real, dimension(num_trial_res) :: trial_re, g, w0, predicted_Refl_nir
      ! --------------------------
    
    if(any(phase == (/ phaseIsLiquid, phaseIsUndetermined, phaseIsIce /))) then 
      if (phase == phaseIsLiquid .OR. phase == phaseIsUndetermined) then
        re_min = re_water_min
        re_max = re_water_max
        trial_re(:) = trial_re_w
        g(:)   = g_w(:) 
        w0(:)  = w0_w(:)
      else
        re_min = re_ice_min
        re_max = re_ice_max
        trial_re(:) = trial_re_i
        g(:)   = g_i(:) 
        w0(:)  = w0_i(:)
      end if
      !
      ! 1st attempt at index: w/coarse re resolution
      !
      predicted_Refl_nir(:) = two_stream_reflectance(tau, g(:), w0(:))
      retrieve_re = interpolate_to_min(trial_re(:), predicted_Refl_nir(:), obs_Refl_nir) 
      !
      ! If first retrieval works, can try 2nd iteration using greater re resolution 
      !
      if(use_two_re_iterations .and. retrieve_re > 0.) then
        re_min = retrieve_re - delta_re
        re_max = retrieve_re + delta_re
        delta_re = (re_max - re_min)/real(num_trial_res-1)
  
        trial_re(:) = re_min + delta_re * (/ (i - 1, i = 1, num_trial_res) /) 
        g(:)  = get_g_nir(  phase, trial_re(:))
        w0(:) = get_ssa_nir(phase, trial_re(:))
        predicted_Refl_nir(:) = two_stream_reflectance(tau, g(:), w0(:))
        retrieve_re = interpolate_to_min(trial_re(:), predicted_Refl_nir(:), obs_Refl_nir) 
      end if
    else 
      retrieve_re = re_fill
    end if 
    
  end function retrieve_re
  ! --------------------------------------------
  pure function interpolate_to_min(x, y, yobs)
    real, dimension(:), intent(in) :: x, y 
    real,               intent(in) :: yobs
    real                           :: interpolate_to_min
    ! 
    ! Given a set of values of y as y(x), find the value of x that minimizes abs(y - yobs)
    !   y must be monotonic in x
    !
    real, dimension(size(x)) :: diff
    integer                  :: nPoints, minDiffLoc, lowerBound, upperBound
    ! ---------------------------------
    nPoints = size(y)
    diff(:) = y(:) - yobs
    minDiffLoc = minloc(abs(diff), dim = 1) 
    
    if(minDiffLoc == 1) then 
      lowerBound = minDiffLoc
      upperBound = minDiffLoc + 1
    else if(minDiffLoc == nPoints) then
      lowerBound = minDiffLoc - 1
      upperBound = minDiffLoc
    else
      if(diff(minDiffLoc-1) * diff(minDiffLoc) < 0) then
        lowerBound = minDiffLoc-1
        upperBound = minDiffLoc
      else 
        lowerBound = minDiffLoc
        upperBound = minDiffLoc + 1
      end if 
    end if 
    
    if(diff(lowerBound) * diff(upperBound) < 0) then     
      !
      ! Interpolate the root position linearly if we bracket the root
      !
      interpolate_to_min = x(upperBound) - & 
                           diff(upperBound) * (x(upperBound) - x(lowerBound)) / (diff(upperBound) - diff(lowerBound))
    else 
      interpolate_to_min = re_fill
    end if 
    

  end function interpolate_to_min
  ! --------------------------------------------
  ! Optical properties
  ! --------------------------------------------
  elemental function get_g_nir (phase, re)
    !
    ! Polynomial fit for asummetry parameter g in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick
    !

    integer, intent(in) :: phase
    real,    intent(in) :: re
    real :: get_g_nir 
    
    real, dimension(3), parameter :: ice_coefficients   = (/ 0.7432,  4.5563e-3, -2.8697e-5 /), & 
                               small_water_coefficients = (/ 0.8027, -1.0496e-2,  1.7071e-3 /), & 
                                 big_water_coefficients = (/ 0.7931,  5.3087e-3, -7.4995e-5 /) 
    
    ! approx. fits from MODIS Collection 5 LUT scattering calculations
    if(phase == phaseIsLiquid) then
      if(re < 8.) then 
        get_g_nir = fit_to_quadratic(re, small_water_coefficients)
        if(re < re_water_min) get_g_nir = fit_to_quadratic(re_water_min, small_water_coefficients)
      else
        get_g_nir = fit_to_quadratic(re,   big_water_coefficients)
        if(re > re_water_max) get_g_nir = fit_to_quadratic(re_water_max, big_water_coefficients)
      end if 
    else
      get_g_nir = fit_to_quadratic(re, ice_coefficients)
      if(re < re_ice_min) get_g_nir = fit_to_quadratic(re_ice_min, ice_coefficients)
      if(re > re_ice_max) get_g_nir = fit_to_quadratic(re_ice_max, ice_coefficients)
    end if 
    
  end function get_g_nir

  ! --------------------------------------------
    elemental function get_ssa_nir (phase, re)
        integer, intent(in) :: phase
        real,    intent(in) :: re
        real                :: get_ssa_nir
        !
        ! Polynomial fit for single scattering albedo in MODIS band 7 (near IR) as a function 
        !   of size for ice and water
        ! Fits from Steve Platnick
        !
        
        real, dimension(4), parameter :: ice_coefficients   = (/ 0.9994, -4.5199e-3, 3.9370e-5, -1.5235e-7 /)
        real, dimension(3), parameter :: water_coefficients = (/ 1.0008, -2.5626e-3, 1.6024e-5 /) 
        
        ! approx. fits from MODIS Collection 5 LUT scattering calculations
        if(phase == phaseIsLiquid) then
          get_ssa_nir = fit_to_quadratic(re, water_coefficients)
          if(re < re_water_min) get_ssa_nir = fit_to_quadratic(re_water_min, water_coefficients)
          if(re > re_water_max) get_ssa_nir = fit_to_quadratic(re_water_max, water_coefficients)
        else
          get_ssa_nir = fit_to_cubic(re, ice_coefficients)
          if(re < re_ice_min) get_ssa_nir = fit_to_cubic(re_ice_min, ice_coefficients)
          if(re > re_ice_max) get_ssa_nir = fit_to_cubic(re_ice_max, ice_coefficients)
        end if 

    end function get_ssa_nir
   ! --------------------------------------------
  pure function fit_to_cubic(x, coefficients) 
    real,               intent(in) :: x
    real, dimension(:), intent(in) :: coefficients
    real                           :: fit_to_cubic
    
    
    fit_to_cubic = coefficients(1) + x * (coefficients(2) + x * (coefficients(3) + x * coefficients(4)))
 end function fit_to_cubic
   ! --------------------------------------------
  pure function fit_to_quadratic(x, coefficients) 
    real,               intent(in) :: x
    real, dimension(:), intent(in) :: coefficients
    real                           :: fit_to_quadratic
    
    
    fit_to_quadratic = coefficients(1) + x * (coefficients(2) + x * (coefficients(3)))
 end function fit_to_quadratic
  ! --------------------------------------------
  ! Radiative transfer
  ! --------------------------------------------
  pure function compute_toa_reflectace(tau, g, w0)
    real, dimension(:), intent(in) :: tau, g, w0
    real                           :: compute_toa_reflectace
    
    logical, dimension(size(tau))         :: cloudMask
    integer, dimension(count(tau(:) > 0)) :: cloudIndicies
    real,    dimension(count(tau(:) > 0)) :: Refl,     Trans
    real                                  :: Refl_tot, Trans_tot
    integer                               :: i
    ! ---------------------------------------
    !
    ! This wrapper reports reflectance only and strips out non-cloudy elements from the calculation
    !
    cloudMask = tau(:) > 0. 
    cloudIndicies = pack((/ (i, i = 1, size(tau)) /), mask = cloudMask) 
    do i = 1, size(cloudIndicies)
      call two_stream(tau(cloudIndicies(i)), g(cloudIndicies(i)), w0(cloudIndicies(i)), Refl(i), Trans(i))
    end do 
                    
    call adding_doubling(Refl(:), Trans(:), Refl_tot, Trans_tot)  
    
    compute_toa_reflectace = Refl_tot
    
  end function compute_toa_reflectace
  ! --------------------------------------------
  pure subroutine two_stream(tauint, gint, w0int, ref, tra) 
    real, intent(in)  :: tauint, gint, w0int
    real, intent(out) :: ref, tra
    !
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    !
    ! ------------------------
    ! Local variables 
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    integer, parameter :: beam = 2
    real,    parameter :: xmu = 0.866, minConservativeW0 = 0.9999999
    real :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
            rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den, th
    !
    ! Compute reflectance and transmittance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    !
    f   = gint**2
    tau = (1 - w0int * f) * tauint
    w0  = (1 - f) * w0int / (1 - w0int * f)
    g   = (gint - f) / (1 - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7 - w0* (4 + 3 * g)) / 4.0
    gamma2 = -(1 - w0* (4 - 3 * g)) / 4.0
    gamma3 =  (2 - 3*g*xmu) / 4.0
    gamma4 =   1 - gamma3

    if (w0int > minConservativeW0) then
      ! Conservative scattering
      if (beam == 1) then
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-exp(-tau/xmu)))
          ref = rh / (1 + gamma1 * tau)
          tra = 1 - ref       
      else if(beam == 2) then
          ref = gamma1*tau/(1 + gamma1*tau)
          tra = 1 - ref
      endif
    else
      ! Non-conservative scattering
      a1 = gamma1 * gamma4 + gamma2 * gamma3
      a2 = gamma1 * gamma3 + gamma2 * gamma4

      rk = sqrt(gamma1**2 - gamma2**2)
      
      r1 = (1 - rk * xmu) * (a2 + rk * gamma3)
      r2 = (1 + rk * xmu) * (a2 - rk * gamma3)
      r3 = 2 * rk *(gamma3 - a2 * xmu)
      r4 = (1 - (rk * xmu)**2) * (rk + gamma1)
      r5 = (1 - (rk * xmu)**2) * (rk - gamma1)
      
      t1 = (1 + rk * xmu) * (a1 + rk * gamma4)
      t2 = (1 - rk * xmu) * (a1 - rk * gamma4)
      t3 = 2 * rk * (gamma4 + a1 * xmu)
      t4 = r4
      t5 = r5

      beta = -r5 / r4         
      
      e1 = min(rk * tau, 500.) 
      e2 = min(tau / xmu, 500.) 
      
      if (beam == 1) then
         den = r4 * exp(e1) + r5 * exp(-e1)
         ref  = w0*(r1*exp(e1)-r2*exp(-e1)-r3*exp(-e2))/den
         den = t4 * exp(e1) + t5 * exp(-e1)
         th  = exp(-e2)
         tra = th-th*w0*(t1*exp(e1)-t2*exp(-e1)-t3*exp(e2))/den
      elseif (beam == 2) then
         ef1 = exp(-e1)
         ef2 = exp(-2*e1)
         ref = (gamma2*(1-ef2))/((rk+gamma1)*(1-beta*ef2))
         tra = (2*rk*ef1)/((rk+gamma1)*(1-beta*ef2))
      endif
    end if
  end subroutine two_stream
  ! --------------------------------------------------
  elemental function two_stream_reflectance(tauint, gint, w0int) 
    real, intent(in) :: tauint, gint, w0int
    real             :: two_stream_reflectance
    !
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    !
    ! ------------------------
    ! Local variables 
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    integer, parameter :: beam = 2
    real,    parameter :: xmu = 0.866, minConservativeW0 = 0.9999999
    real :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
            rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den
    ! ------------------------


    f   = gint**2
    tau = (1 - w0int * f) * tauint
    w0  = (1 - f) * w0int / (1 - w0int * f)
    g   = (gint - f) / (1 - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7 - w0* (4 + 3 * g)) / 4.0
    gamma2 = -(1 - w0* (4 - 3 * g)) / 4.0
    gamma3 =  (2 - 3*g*xmu) / 4.0
    gamma4 =   1 - gamma3

    if (w0int > minConservativeW0) then
      ! Conservative scattering
      if (beam == 1) then
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-exp(-tau/xmu)))
          two_stream_reflectance = rh / (1 + gamma1 * tau)
      elseif (beam == 2) then
          two_stream_reflectance = gamma1*tau/(1 + gamma1*tau)
      endif
        
    else    !

        ! Non-conservative scattering
         a1 = gamma1 * gamma4 + gamma2 * gamma3
         a2 = gamma1 * gamma3 + gamma2 * gamma4

         rk = sqrt(gamma1**2 - gamma2**2)
         
         r1 = (1 - rk * xmu) * (a2 + rk * gamma3)
         r2 = (1 + rk * xmu) * (a2 - rk * gamma3)
         r3 = 2 * rk *(gamma3 - a2 * xmu)
         r4 = (1 - (rk * xmu)**2) * (rk + gamma1)
         r5 = (1 - (rk * xmu)**2) * (rk - gamma1)
         
         t1 = (1 + rk * xmu) * (a1 + rk * gamma4)
         t2 = (1 - rk * xmu) * (a1 - rk * gamma4)
         t3 = 2 * rk * (gamma4 + a1 * xmu)
         t4 = r4
         t5 = r5

         beta = -r5 / r4         
         
         e1 = min(rk * tau, 500.) 
         e2 = min(tau / xmu, 500.) 
         
         if (beam == 1) then
           den = r4 * exp(e1) + r5 * exp(-e1)
           two_stream_reflectance  = w0*(r1*exp(e1)-r2*exp(-e1)-r3*exp(-e2))/den
         elseif (beam == 2) then
           ef1 = exp(-e1)
           ef2 = exp(-2*e1)
           two_stream_reflectance = (gamma2*(1-ef2))/((rk+gamma1)*(1-beta*ef2))
         endif
           
      end if
  end function two_stream_reflectance 
  ! --------------------------------------------
    pure subroutine adding_doubling (Refl, Tran, Refl_tot, Tran_tot)      
      real,    dimension(:), intent(in)  :: Refl,     Tran
      real,                  intent(out) :: Refl_tot, Tran_tot
      !
      ! Use adding/doubling formulas to compute total reflectance and transmittance from layer values
      !
      
      integer :: i
      real, dimension(size(Refl)) :: Refl_cumulative, Tran_cumulative
      
      Refl_cumulative(1) = Refl(1); Tran_cumulative(1) = Tran(1)    
      
      do i=2, size(Refl)
          ! place (add) previous combined layer(s) reflectance on top of layer i, w/black surface (or ignoring surface):
          Refl_cumulative(i) = Refl_cumulative(i-1) + Refl(i)*(Tran_cumulative(i-1)**2)/(1 - Refl_cumulative(i-1) * Refl(i))
          Tran_cumulative(i) = (Tran_cumulative(i-1)*Tran(i)) / (1 - Refl_cumulative(i-1) * Refl(i))
      end do
      
      Refl_tot = Refl_cumulative(size(Refl))
      Tran_tot = Tran_cumulative(size(Refl))

    end subroutine adding_doubling
  ! --------------------------------------------------
!  subroutine complain_and_die(message) 
!    character(len = *), intent(in) :: message
    
!    write(6, *) "Failure in MODIS simulator" 
!    write(6, *)  trim(message)
     !! Flush is a Fortran 2003 feature.
!    flush(6)
!    stop
!  end subroutine complain_and_die
  !------------------------------------------------------------------------------------------------
end module mod_modis_sim
