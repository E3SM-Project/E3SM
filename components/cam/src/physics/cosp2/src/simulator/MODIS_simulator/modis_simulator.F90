! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History
! May 2009:      Robert Pincus - Initial version
! June 2009:     Steve Platnick and Robert Pincus - Simple radiative transfer for size 
!                retrievals
! August 2009:   Robert Pincus - Consistency and bug fixes suggested by Rick Hemler (GFDL) 
! November 2009: Robert Pincus - Bux fixes and speed-ups after experience with Rick Hemler 
!                using AM2 (GFDL) 
! January 2010:  Robert Pincus - Added high, middle, low cloud fractions
! May 2015:      Dustin Swales - Modified for COSPv2.0
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Notes on using the MODIS simulator: 
!  *) You may provide either layer-by-layer values of optical thickness at 0.67 and 2.1
!     microns, or optical thickness at 0.67 microns and ice- and liquid-water contents 
!     (in consistent units of your choosing)
!  *) Required input also includes the optical thickness and cloud top pressure 
!     derived from the ISCCP simulator run with parameter top_height = 1. 
!  *) Cloud particle sizes are specified as radii, measured in meters, though within the 
!     module we use units of microns. Where particle sizes are outside the bounds used in 
!     the MODIS retrieval libraries (parameters re_water_min, re_ice_min, etc.) the 
!     simulator returns missing values (re_fill)
!
! When error conditions are encountered this code calls the function complain_and_die, 
! supplied at the bottom of this module. Users probably want to replace this with 
! something more graceful. 
!
module mod_modis_sim
  USE MOD_COSP_CONFIG, only: R_UNDEF,modis_histTau,modis_histPres,numMODISTauBins,       &
                             numMODISPresBins,numMODISReffIceBins,numMODISReffLiqBins,   &
                             modis_histReffIce,modis_histReffLiq
  USE COSP_KINDS,      ONLY: wp
  use MOD_COSP_STATS,  ONLY: hist2D

  implicit none
  ! ##########################################################################
  ! Retrieval parameters
   integer, parameter :: &
        num_trial_res = 15              ! Increase to make the linear pseudo-retrieval of size more accurate
   
   real(wp) :: &
       min_OpticalThickness,          & ! Minimum detectable optical thickness
       CO2Slicing_PressureLimit,      & ! Cloud with higher pressures use thermal methods, units Pa
       CO2Slicing_TauLimit,           & ! How deep into the cloud does CO2 slicing see? 
       phase_TauLimit,                & ! How deep into the cloud does the phase detection see?
       size_TauLimit,                 & ! Depth of the re retreivals
       phaseDiscrimination_Threshold, & ! What fraction of total extincton needs to be in a single
                                        ! category to make phase discrim. work? 
       re_fill,                       & !
       re_water_min,                  & ! Minimum effective radius (liquid)
       re_water_max,                  & ! Maximum effective radius (liquid)
       re_ice_min,                    & ! Minimum effective radius (ice)
       re_ice_max,                    & ! Minimum effective radius (ice)
       highCloudPressureLimit,        & ! High cloud pressure limit (Pa)
       lowCloudPressureLimit            ! Low cloud pressure limit (Pa)
  integer :: &
       phaseIsNone,                   & !
       phaseIsLiquid,                 & !
       phaseIsIce,                    & !
       phaseIsUndetermined              !
 
  real(wp),dimension(num_trial_res) :: &
       trial_re_w, & ! Near-IR optical params vs size for retrieval scheme (liquid)
       trial_re_i    ! Near-IR optical params vs size for retrieval scheme (ice)
  real(wp),dimension(num_trial_res) :: &
       g_w,        & ! Assymettry parameter for size retrieval (liquid)
       g_i,        & ! Assymettry parameter for size retrieval (ice)
       w0_w,       & ! Single-scattering albedo for size retrieval (liquid)
       w0_i          ! Single-scattering albedo for size retrieval (ice)
  ! Algorithmic parameters
  real(wp),parameter :: &
     ice_density = 0.93_wp ! Liquid density is 1. 
      
contains
  ! ########################################################################################
  ! MODIS simulator using specified liquid and ice optical thickness in each layer 
  !
  ! Note: this simulator operates on all points; to match MODIS itself night-time 
  !       points should be excluded
  !
  ! Note: the simulator requires as input the optical thickness and cloud top pressure 
  !       derived from the ISCCP simulator run with parameter top_height = 1. 
  !       If cloud top pressure is higher than about 700 mb, MODIS can't use CO2 slicing 
  !       and reverts to a thermal algorithm much like ISCCP's. Rather than replicate that 
  !       alogrithm in this simulator we simply report the values from the ISCCP simulator. 
  ! ########################################################################################
  subroutine modis_subcolumn(nSubCols, nLevels, pressureLevels, optical_thickness,       & 
                         tauLiquidFraction, g, w0,isccpCloudTopPressure,                 &
                         retrievedPhase, retrievedCloudTopPressure,                      &
                         retrievedTau,   retrievedSize)

    ! INPUTS
    integer,intent(in) :: &
         nSubCols,                  & ! Number of subcolumns
         nLevels                      ! Number of levels         
    real(wp),dimension(nLevels+1),intent(in) :: &
         pressureLevels               ! Gridmean pressure at layer edges (Pa)                  
    real(wp),dimension(nSubCols,nLevels),intent(in) :: &
         optical_thickness,         & ! Subcolumn optical thickness @ 0.67 microns.
         tauLiquidFraction,         & ! Liquid water fraction
         g,                         & ! Subcolumn assymetry parameter  
         w0                           ! Subcolumn single-scattering albedo 
    real(wp),dimension(nSubCols),intent(in) :: &
         isccpCloudTopPressure        ! ISCCP retrieved cloud top pressure (Pa)

    ! OUTPUTS
    integer, dimension(nSubCols), intent(inout) :: &
         retrievedPhase               ! MODIS retrieved phase (liquid/ice/other)              
    real(wp),dimension(nSubCols), intent(inout) :: &
         retrievedCloudTopPressure, & ! MODIS retrieved CTP (Pa)
         retrievedTau,              & ! MODIS retrieved optical depth (unitless)              
         retrievedSize                ! MODIS retrieved particle size (microns)              

    ! LOCAL VARIABLES
    logical, dimension(nSubCols)      :: &
         cloudMask
    real(wp)                          :: &
         integratedLiquidFraction,       &
         obs_Refl_nir
    real(wp),dimension(num_trial_res) :: &
         predicted_Refl_nir
    integer                           :: &
         i

    ! ########################################################################################
    !                           Optical depth retrieval 
    ! This is simply a sum over the optical thickness in each layer. 
    ! It should agree with the ISCCP values after min values have been excluded.
    ! ########################################################################################
    retrievedTau(1:nSubCols) = sum(optical_thickness(1:nSubCols,1:nLevels), dim = 2)

    ! ########################################################################################
    !                                 Cloud detection
    ! does optical thickness exceed detection threshold? 
    ! ########################################################################################
    cloudMask = retrievedTau(1:nSubCols) >= min_OpticalThickness
    
    do i = 1, nSubCols
       if(cloudMask(i)) then 
          ! ##################################################################################
          !                       Cloud top pressure determination 
          ! MODIS uses CO2 slicing for clouds with tops above about 700 mb and thermal 
          ! methods for clouds lower than that. For CO2 slicing we report the optical-depth
          ! weighted pressure, integrating to a specified optical depth.
          ! This assumes linear variation in p between levels. Linear in ln(p) is probably 
          ! better, though we'd need to deal with the lowest pressure gracefully. 
          ! ##################################################################################
          retrievedCloudTopPressure(i) = cloud_top_pressure(nLevels,(/ 0._wp, optical_thickness(i,1:nLevels) /), &
                                                            pressureLevels(1:nLevels),CO2Slicing_TauLimit)  
        
          ! ##################################################################################
          !                               Phase determination 
          ! Determine fraction of total tau that's liquid when ice and water contribute about 
          ! equally to the extinction we can't tell what the phase is.
          ! ##################################################################################
          integratedLiquidFraction = weight_by_extinction(nLevels,optical_thickness(i,1:nLevels),       &
                                                          tauLiquidFraction(i, 1:nLevels), &
                                                          phase_TauLimit)
          if(integratedLiquidFraction >= phaseDiscrimination_Threshold) then 
             retrievedPhase(i) = phaseIsLiquid
          else if (integratedLiquidFraction <= 1._wp- phaseDiscrimination_Threshold) then 
             retrievedPhase(i) = phaseIsIce
          else 
             retrievedPhase(i) = phaseIsUndetermined
          end if
        
          ! ##################################################################################
          !                                 Size determination 
          ! ##################################################################################
          
          ! Compute observed reflectance
          obs_Refl_nir = compute_toa_reflectace(nLevels,optical_thickness(i,1:nLevels), g(i,1:nLevels), w0(i,1:nLevels))

          ! Compute predicted reflectance
          if(any(retrievedPhase(i) == (/ phaseIsLiquid, phaseIsUndetermined, phaseIsIce /))) then 
             if (retrievedPhase(i) == phaseIsLiquid .OR. retrievedPhase(i) == phaseIsUndetermined) then
                predicted_Refl_nir(1:num_trial_res) = two_stream_reflectance(retrievedTau(i), &
                     g_w(1:num_trial_res), w0_w(1:num_trial_res))
                retrievedSize(i) = 1.0e-06_wp*interpolate_to_min(trial_re_w(1:num_trial_res), &
                     predicted_Refl_nir(1:num_trial_res), obs_Refl_nir)
             else
                predicted_Refl_nir(1:num_trial_res) = two_stream_reflectance(retrievedTau(i), &
                     g_i(1:num_trial_res), w0_i(1:num_trial_res))
                retrievedSize(i) = 1.0e-06_wp*interpolate_to_min(trial_re_i(1:num_trial_res), &
                     predicted_Refl_nir(1:num_trial_res), obs_Refl_nir)
             endif
          else 
             retrievedSize(i) = re_fill
          endif
       else   
          ! Values when we don't think there's a cloud. 
          retrievedCloudTopPressure(i) = R_UNDEF 
          retrievedPhase(i)            = phaseIsNone
          retrievedSize(i)             = R_UNDEF 
          retrievedTau(i)              = R_UNDEF 
       end if
    end do
    where((retrievedSize(1:nSubCols) < 0.).and.(retrievedSize(1:nSubCols) /= R_UNDEF)) &
         retrievedSize(1:nSubCols) = 1.0e-06_wp*re_fill

    ! We use the ISCCP-derived CTP for low clouds, since the ISCCP simulator ICARUS 
    ! mimics what MODIS does to first order. 
    ! Of course, ISCCP cloud top pressures are in mb.   
    where(cloudMask(1:nSubCols) .and. retrievedCloudTopPressure(1:nSubCols) > CO2Slicing_PressureLimit) &
         retrievedCloudTopPressure(1:nSubCols) = isccpCloudTopPressure! * 100._wp
    
  end subroutine modis_subcolumn

  ! ########################################################################################
  subroutine modis_column(nPoints,nSubCols,phase, cloud_top_pressure, optical_thickness, particle_size,     &
       Cloud_Fraction_Total_Mean,         Cloud_Fraction_Water_Mean,         Cloud_Fraction_Ice_Mean,        &
       Cloud_Fraction_High_Mean,          Cloud_Fraction_Mid_Mean,           Cloud_Fraction_Low_Mean,        &
       Optical_Thickness_Total_Mean,      Optical_Thickness_Water_Mean,      Optical_Thickness_Ice_Mean,     &
       Optical_Thickness_Total_MeanLog10, Optical_Thickness_Water_MeanLog10, Optical_Thickness_Ice_MeanLog10,&
       Cloud_Particle_Size_Water_Mean,    Cloud_Particle_Size_Ice_Mean,      Cloud_Top_Pressure_Total_Mean,  &
       Liquid_Water_Path_Mean,            Ice_Water_Path_Mean,                                               &    
       Optical_Thickness_vs_Cloud_Top_Pressure,Optical_Thickness_vs_ReffIce,Optical_Thickness_vs_ReffLiq)
    
    ! INPUTS
    integer,intent(in) :: &
         nPoints,                           & ! Number of horizontal gridpoints
         nSubCols                             ! Number of subcolumns
    integer,intent(in), dimension(nPoints, nSubCols) ::  &
         phase                             
    real(wp),intent(in),dimension(nPoints, nSubCols) ::  &
         cloud_top_pressure,                &
         optical_thickness,                 &
         particle_size
 
    ! OUTPUTS 
    real(wp),intent(inout),dimension(nPoints)  ::   & !
         Cloud_Fraction_Total_Mean,         & !
         Cloud_Fraction_Water_Mean,         & !
         Cloud_Fraction_Ice_Mean,           & !
         Cloud_Fraction_High_Mean,          & !
         Cloud_Fraction_Mid_Mean,           & !
         Cloud_Fraction_Low_Mean,           & !
         Optical_Thickness_Total_Mean,      & !
         Optical_Thickness_Water_Mean,      & !
         Optical_Thickness_Ice_Mean,        & !
         Optical_Thickness_Total_MeanLog10, & !
         Optical_Thickness_Water_MeanLog10, & !
         Optical_Thickness_Ice_MeanLog10,   & !
         Cloud_Particle_Size_Water_Mean,    & !
         Cloud_Particle_Size_Ice_Mean,      & !
         Cloud_Top_Pressure_Total_Mean,     & !
         Liquid_Water_Path_Mean,            & !
         Ice_Water_Path_Mean                  !
    real(wp),intent(inout),dimension(nPoints,numMODISTauBins,numMODISPresBins) :: &
         Optical_Thickness_vs_Cloud_Top_Pressure
    real(wp),intent(inout),dimension(nPoints,numMODISTauBins,numMODISReffIceBins) :: &    
         Optical_Thickness_vs_ReffIce
    real(wp),intent(inout),dimension(nPoints,numMODISTauBins,numMODISReffLiqBins) :: &    
         Optical_Thickness_vs_ReffLiq         

    ! LOCAL VARIABLES
    real(wp), parameter :: &
         LWP_conversion = 2._wp/3._wp * 1000._wp ! MKS units  
    integer :: j
    logical, dimension(nPoints,nSubCols) :: &
         cloudMask,      &
         waterCloudMask, &
         iceCloudMask,   &
         validRetrievalMask
    real(wp),dimension(nPoints,nSubCols) :: &
         tauWRK,ctpWRK,reffIceWRK,reffLiqWRK

    ! ########################################################################################
    ! Include only those pixels with successful retrievals in the statistics 
    ! ########################################################################################
    validRetrievalMask(1:nPoints,1:nSubCols) = particle_size(1:nPoints,1:nSubCols) > 0.
    cloudMask(1:nPoints,1:nSubCols) = phase(1:nPoints,1:nSubCols) /= phaseIsNone .and.       &
         validRetrievalMask(1:nPoints,1:nSubCols)
    waterCloudMask(1:nPoints,1:nSubCols) = phase(1:nPoints,1:nSubCols) == phaseIsLiquid .and. &
         validRetrievalMask(1:nPoints,1:nSubCols)
    iceCloudMask(1:nPoints,1:nSubCols)   = phase(1:nPoints,1:nSubCols) == phaseIsIce .and.    &
         validRetrievalMask(1:nPoints,1:nSubCols)
    
    ! ########################################################################################
    ! Use these as pixel counts at first 
    ! ########################################################################################
    Cloud_Fraction_Total_Mean(1:nPoints) = real(count(cloudMask,      dim = 2))
    Cloud_Fraction_Water_Mean(1:nPoints) = real(count(waterCloudMask, dim = 2))
    Cloud_Fraction_Ice_Mean(1:nPoints)   = real(count(iceCloudMask,   dim = 2))
    Cloud_Fraction_High_Mean(1:nPoints)  = real(count(cloudMask .and. cloud_top_pressure <=          &
                                           highCloudPressureLimit, dim = 2)) 
    Cloud_Fraction_Low_Mean(1:nPoints)   = real(count(cloudMask .and. cloud_top_pressure >           &
                                           lowCloudPressureLimit,  dim = 2)) 
    Cloud_Fraction_Mid_Mean(1:nPoints)   = Cloud_Fraction_Total_Mean(1:nPoints) - Cloud_Fraction_High_Mean(1:nPoints)&
                                           - Cloud_Fraction_Low_Mean(1:nPoints)

    ! ########################################################################################
    ! Compute column amounts.
    ! ########################################################################################
    where(Cloud_Fraction_Total_Mean(1:nPoints) > 0)
       Optical_Thickness_Total_Mean(1:nPoints) = sum(optical_thickness, mask = cloudMask,      dim = 2) / &
            Cloud_Fraction_Total_Mean(1:nPoints)
       Optical_Thickness_Total_MeanLog10(1:nPoints) = sum(log10(abs(optical_thickness)), mask = cloudMask, &
            dim = 2) / Cloud_Fraction_Total_Mean(1:nPoints)
    elsewhere
       Optical_Thickness_Total_Mean      = R_UNDEF
       Optical_Thickness_Total_MeanLog10 = R_UNDEF
    endwhere
    where(Cloud_Fraction_Water_Mean(1:nPoints) > 0)
       Optical_Thickness_Water_Mean(1:nPoints) = sum(optical_thickness, mask = waterCloudMask, dim = 2) / &
            Cloud_Fraction_Water_Mean(1:nPoints)
       Liquid_Water_Path_Mean(1:nPoints) = LWP_conversion*sum(particle_size*optical_thickness, &
            mask=waterCloudMask,dim=2)/Cloud_Fraction_Water_Mean(1:nPoints)
       Optical_Thickness_Water_MeanLog10(1:nPoints) = sum(log10(abs(optical_thickness)), mask = waterCloudMask,&
            dim = 2) / Cloud_Fraction_Water_Mean(1:nPoints)
       Cloud_Particle_Size_Water_Mean(1:nPoints) = sum(particle_size, mask = waterCloudMask, dim = 2) / &
            Cloud_Fraction_Water_Mean(1:nPoints)
    elsewhere
       Optical_Thickness_Water_Mean      = R_UNDEF
       Optical_Thickness_Water_MeanLog10 = R_UNDEF
       Cloud_Particle_Size_Water_Mean    = R_UNDEF
       Liquid_Water_Path_Mean            = R_UNDEF
    endwhere
    where(Cloud_Fraction_Ice_Mean(1:nPoints) > 0)
       Optical_Thickness_Ice_Mean(1:nPoints)   = sum(optical_thickness, mask = iceCloudMask,   dim = 2) / &
            Cloud_Fraction_Ice_Mean(1:nPoints)
       Ice_Water_Path_Mean(1:nPoints) = LWP_conversion * ice_density*sum(particle_size*optical_thickness,&
            mask=iceCloudMask,dim = 2) /Cloud_Fraction_Ice_Mean(1:nPoints) 
       Optical_Thickness_Ice_MeanLog10(1:nPoints) = sum(log10(abs(optical_thickness)), mask = iceCloudMask,&
            dim = 2) / Cloud_Fraction_Ice_Mean(1:nPoints)
       Cloud_Particle_Size_Ice_Mean(1:nPoints) = sum(particle_size, mask = iceCloudMask,   dim = 2) / &
            Cloud_Fraction_Ice_Mean(1:nPoints)    
    elsewhere
       Optical_Thickness_Ice_Mean        = R_UNDEF
       Optical_Thickness_Ice_MeanLog10   = R_UNDEF
       Cloud_Particle_Size_Ice_Mean      = R_UNDEF
       Ice_Water_Path_Mean               = R_UNDEF
    endwhere
    Cloud_Top_Pressure_Total_Mean  = sum(cloud_top_pressure, mask = cloudMask, dim = 2) / &
                                     max(1, count(cloudMask, dim = 2))

    ! ########################################################################################
    ! Normalize pixel counts to fraction. 
    ! ########################################################################################
    Cloud_Fraction_High_Mean(1:nPoints)  = Cloud_Fraction_High_Mean(1:nPoints)  /nSubcols
    Cloud_Fraction_Mid_Mean(1:nPoints)   = Cloud_Fraction_Mid_Mean(1:nPoints)   /nSubcols
    Cloud_Fraction_Low_Mean(1:nPoints)   = Cloud_Fraction_Low_Mean(1:nPoints)   /nSubcols
    Cloud_Fraction_Total_Mean(1:nPoints) = Cloud_Fraction_Total_Mean(1:nPoints) /nSubcols
    Cloud_Fraction_Ice_Mean(1:nPoints)   = Cloud_Fraction_Ice_Mean(1:nPoints)   /nSubcols
    Cloud_Fraction_Water_Mean(1:nPoints) = Cloud_Fraction_Water_Mean(1:nPoints) /nSubcols
    
    ! ########################################################################################
    ! Joint histograms
    ! ########################################################################################
    ! Loop over all points
    tauWRK(1:nPoints,1:nSubCols)     = optical_thickness(1:nPoints,1:nSubCols)
    ctpWRK(1:nPoints,1:nSubCols)     = cloud_top_pressure(1:nPoints,1:nSubCols)
    reffIceWRK(1:nPoints,1:nSubCols) = merge(particle_size,R_UNDEF,iceCloudMask)
    reffLiqWRK(1:nPoints,1:nSubCols) = merge(particle_size,R_UNDEF,waterCloudMask)
    do j=1,nPoints

       ! Fill clear and optically thin subcolumns with fill
       where(.not. cloudMask(j,1:nSubCols)) 
          tauWRK(j,1:nSubCols) = -999._wp
          ctpWRK(j,1:nSubCols) = -999._wp
       endwhere
       ! Joint histogram of tau/CTP
       call hist2D(tauWRK(j,1:nSubCols),ctpWRK(j,1:nSubCols),nSubCols,&
                   modis_histTau,numMODISTauBins,&
                   modis_histPres,numMODISPresBins,&
                   Optical_Thickness_vs_Cloud_Top_Pressure(j,1:numMODISTauBins,1:numMODISPresBins))
       ! Joint histogram of tau/ReffICE
       call hist2D(tauWRK(j,1:nSubCols),reffIceWrk(j,1:nSubCols),nSubCols,               &
                   modis_histTau,numMODISTauBins,modis_histReffIce,         &
                   numMODISReffIceBins, Optical_Thickness_vs_ReffIce(j,1:numMODISTauBins,1:numMODISReffIceBins))
       ! Joint histogram of tau/ReffLIQ
       call hist2D(tauWRK(j,1:nSubCols),reffLiqWrk(j,1:nSubCols),nSubCols,               &
                   modis_histTau,numMODISTauBins,modis_histReffLiq,         &
                   numMODISReffLiqBins, Optical_Thickness_vs_ReffLiq(j,1:numMODISTauBins,1:numMODISReffLiqBins))                   

    enddo   
    Optical_Thickness_vs_Cloud_Top_Pressure(1:nPoints,1:numMODISTauBins,1:numMODISPresBins) = &
         Optical_Thickness_vs_Cloud_Top_Pressure(1:nPoints,1:numMODISTauBins,1:numMODISPresBins)/nSubCols
    Optical_Thickness_vs_ReffIce(1:nPoints,1:numMODISTauBins,1:numMODISReffIceBins) = &
         Optical_Thickness_vs_ReffIce(1:nPoints,1:numMODISTauBins,1:numMODISReffIceBins)/nSubCols
    Optical_Thickness_vs_ReffLiq(1:nPoints,1:numMODISTauBins,1:numMODISReffLiqBins) = &
         Optical_Thickness_vs_ReffLiq(1:nPoints,1:numMODISTauBins,1:numMODISReffLiqBins)/nSubCols 
                 

    ! Unit conversion
    where(Optical_Thickness_vs_Cloud_Top_Pressure /= R_UNDEF) &
      Optical_Thickness_vs_Cloud_Top_Pressure = Optical_Thickness_vs_Cloud_Top_Pressure*100._wp
    where(Optical_Thickness_vs_ReffIce /= R_UNDEF) Optical_Thickness_vs_ReffIce = Optical_Thickness_vs_ReffIce*100._wp
    where(Optical_Thickness_vs_ReffLiq /= R_UNDEF) Optical_Thickness_vs_ReffLiq = Optical_Thickness_vs_ReffLiq*100._wp
    where(Cloud_Fraction_Total_Mean /= R_UNDEF) Cloud_Fraction_Total_Mean = Cloud_Fraction_Total_Mean*100._wp
    where(Cloud_Fraction_Water_Mean /= R_UNDEF) Cloud_Fraction_Water_Mean = Cloud_Fraction_Water_Mean*100._wp
    where(Cloud_Fraction_Ice_Mean /= R_UNDEF)   Cloud_Fraction_Ice_Mean = Cloud_Fraction_Ice_Mean*100._wp
    where(Cloud_Fraction_High_Mean /= R_UNDEF)  Cloud_Fraction_High_Mean = Cloud_Fraction_High_Mean*100._wp
    where(Cloud_Fraction_Mid_Mean /= R_UNDEF)   Cloud_Fraction_Mid_Mean = Cloud_Fraction_Mid_Mean*100._wp
    where(Cloud_Fraction_Low_Mean /= R_UNDEF)   Cloud_Fraction_Low_Mean = Cloud_Fraction_Low_Mean*100._wp

  end subroutine modis_column

  ! ########################################################################################
  function cloud_top_pressure(nLevels,tauIncrement, pressure, tauLimit) 
    ! INPUTS
    integer, intent(in)                    :: nLevels
    real(wp),intent(in),dimension(nLevels) :: tauIncrement, pressure
    real(wp),intent(in)                    :: tauLimit
    ! OUTPUTS
    real(wp)                               :: cloud_top_pressure
    ! LOCAL VARIABLES
    real(wp)                               :: deltaX, totalTau, totalProduct
    integer                                :: i 
    
    ! Find the extinction-weighted pressure. Assume that pressure varies linearly between 
    !   layers and use the trapezoidal rule.
    totalTau = 0._wp; totalProduct = 0._wp
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
                     + (pressure(i) - pressure(i-1)) * deltaX**2/(2._wp * tauIncrement(i)) 
      else
        totalTau =     totalTau     + tauIncrement(i) 
        totalProduct = totalProduct + tauIncrement(i) * (pressure(i) + pressure(i-1)) / 2._wp
      end if 
      if(totalTau >= tauLimit) exit
    end do 

    if (totalTau > 0._wp) then
       cloud_top_pressure = totalProduct/totalTau
    else
       cloud_top_pressure = 0._wp
    endif
    
  end function cloud_top_pressure

  ! ########################################################################################
  function weight_by_extinction(nLevels,tauIncrement, f, tauLimit) 
    ! INPUTS
    integer, intent(in)                    :: nLevels
    real(wp),intent(in),dimension(nLevels) :: tauIncrement, f
    real(wp),intent(in)                    :: tauLimit
    ! OUTPUTS
    real(wp)                               :: weight_by_extinction
    ! LOCAL VARIABLES
    real(wp)                               :: deltaX, totalTau, totalProduct
    integer                                :: i 
    
    ! Find the extinction-weighted value of f(tau), assuming constant f within each layer
    totalTau = 0._wp; totalProduct = 0._wp
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

    if (totalTau > 0._wp) then
       weight_by_extinction = totalProduct/totalTau
    else
       weight_by_extinction = 0._wp
    endif
    
  end function weight_by_extinction

  ! ########################################################################################
  pure function interpolate_to_min(x, y, yobs)
    ! INPUTS
    real(wp),intent(in),dimension(num_trial_res) :: x, y 
    real(wp),intent(in)                          :: yobs
    ! OUTPUTS
    real(wp)                                     :: interpolate_to_min
    ! LOCAL VARIABLES
    real(wp), dimension(num_trial_res)           :: diff
    integer                                      :: nPoints, minDiffLoc, lowerBound, upperBound
    
    ! Given a set of values of y as y(x), find the value of x that minimizes abs(y - yobs)
    !   y must be monotonic in x
 
    nPoints = size(y)
    diff(1:num_trial_res) = y(1:num_trial_res) - yobs
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

  ! ########################################################################################
  ! Optical properties
  ! ########################################################################################
  elemental function get_g_nir_old (phase, re)
    ! Polynomial fit for asummetry parameter g in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick

    ! INPUTS
    integer, intent(in) :: phase
    real(wp),intent(in) :: re
    ! OUTPUTS
    real(wp)            :: get_g_nir_old 
    ! LOCAL VARIABLES(parameters)
    real(wp), dimension(3), parameter :: &
         ice_coefficients         = (/ 0.7432,  4.5563e-3, -2.8697e-5 /), & 
         small_water_coefficients = (/ 0.8027, -1.0496e-2,  1.7071e-3 /), & 
         big_water_coefficients   = (/ 0.7931,  5.3087e-3, -7.4995e-5 /) 
   
    ! approx. fits from MODIS Collection 5 LUT scattering calculations
    if(phase == phaseIsLiquid) then
      if(re < 8.) then 
        get_g_nir_old = fit_to_quadratic(re, small_water_coefficients)
        if(re < re_water_min) get_g_nir_old = fit_to_quadratic(re_water_min, small_water_coefficients)
      else
        get_g_nir_old = fit_to_quadratic(re,   big_water_coefficients)
        if(re > re_water_max) get_g_nir_old = fit_to_quadratic(re_water_max, big_water_coefficients)
      end if 
    else
      get_g_nir_old = fit_to_quadratic(re, ice_coefficients)
      if(re < re_ice_min) get_g_nir_old = fit_to_quadratic(re_ice_min, ice_coefficients)
      if(re > re_ice_max) get_g_nir_old = fit_to_quadratic(re_ice_max, ice_coefficients)
    end if 
    
  end function get_g_nir_old

  ! ########################################################################################
  elemental function get_ssa_nir_old (phase, re)
    ! Polynomial fit for single scattering albedo in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick
    
    ! INPUTS
    integer, intent(in) :: phase
    real(wp),intent(in) :: re
    ! OUTPUTS
    real(wp)            :: get_ssa_nir_old
    ! LOCAL VARIABLES (parameters)
    real(wp), dimension(4), parameter :: ice_coefficients   = (/ 0.9994, -4.5199e-3, 3.9370e-5, -1.5235e-7 /)
    real(wp), dimension(3), parameter :: water_coefficients = (/ 1.0008, -2.5626e-3, 1.6024e-5 /) 
    
    ! approx. fits from MODIS Collection 5 LUT scattering calculations
    if(phase == phaseIsLiquid) then
       get_ssa_nir_old = fit_to_quadratic(re, water_coefficients)
       if(re < re_water_min) get_ssa_nir_old = fit_to_quadratic(re_water_min, water_coefficients)
       if(re > re_water_max) get_ssa_nir_old = fit_to_quadratic(re_water_max, water_coefficients)
    else
       get_ssa_nir_old = fit_to_cubic(re, ice_coefficients)
       if(re < re_ice_min) get_ssa_nir_old = fit_to_cubic(re_ice_min, ice_coefficients)
       if(re > re_ice_max) get_ssa_nir_old = fit_to_cubic(re_ice_max, ice_coefficients)
    end if
    
  end function get_ssa_nir_old
  
  elemental function get_g_nir (phase, re)
    !
    ! Polynomial fit for asummetry parameter g in MODIS band 7 (near IR) as a function 
    !   of size for ice and water
    ! Fits from Steve Platnick
    !

    integer, intent(in) :: phase
    real(wp),    intent(in) :: re
    real(wp) :: get_g_nir 

    real(wp), dimension(3), parameter :: ice_coefficients         = (/ 0.7490, 6.5153e-3, -5.4136e-5 /), &
                                         small_water_coefficients = (/ 1.0364, -8.8800e-2, 7.0000e-3 /)
    real(wp), dimension(4), parameter :: big_water_coefficients   = (/ 0.6035, 2.8993e-2, -1.1051e-3, 1.5134e-5 /)

    ! approx. fits from MODIS Collection 6 LUT scattering calculations for 3.7 µm channel size retrievals
    if(phase == phaseIsLiquid) then 
       if(re < 7.) then
          get_g_nir = fit_to_quadratic(re, small_water_coefficients)
          if(re < re_water_min) get_g_nir = fit_to_quadratic(re_water_min, small_water_coefficients)
       else
          get_g_nir = fit_to_cubic(re, big_water_coefficients)
          if(re > re_water_max) get_g_nir = fit_to_cubic(re_water_max, big_water_coefficients)
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
        real(wp),    intent(in) :: re
        real(wp)                :: get_ssa_nir
        !
        ! Polynomial fit for single scattering albedo in MODIS band 7 (near IR) as a function 
        !   of size for ice and water
        ! Fits from Steve Platnick
        !
        real(wp), dimension(4), parameter :: ice_coefficients   = (/ 0.9625, -1.8069e-2, 3.3281e-4,-2.2865e-6/)
        real(wp), dimension(3), parameter :: water_coefficients = (/ 1.0044, -1.1397e-2, 1.3300e-4 /)
        
        ! approx. fits from MODIS Collection 6 LUT scattering calculations
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

  

  ! ########################################################################################
  pure function fit_to_cubic(x, coefficients) 
    ! INPUTS
    real(wp),               intent(in) :: x
    real(wp), dimension(4), intent(in) :: coefficients
    ! OUTPUTS
    real(wp)                           :: fit_to_cubic  
    
    fit_to_cubic = coefficients(1) + x * (coefficients(2) + x * (coefficients(3) + x * coefficients(4)))
  end function fit_to_cubic
    
  ! ########################################################################################
  pure function fit_to_quadratic(x, coefficients) 
    ! INPUTS
    real(wp),               intent(in) :: x
    real(wp), dimension(3), intent(in) :: coefficients
    ! OUTPUTS
    real(wp)                           :: fit_to_quadratic
    
    fit_to_quadratic = coefficients(1) + x * (coefficients(2) + x * (coefficients(3)))
  end function fit_to_quadratic

  ! ########################################################################################
  ! Radiative transfer
  ! ########################################################################################
  pure function compute_toa_reflectace(nLevels,tau, g, w0)
    ! This wrapper reports reflectance only and strips out non-cloudy elements from the 
    ! calculation
    
    ! INPUTS
    integer,intent(in)                     :: nLevels
    real(wp),intent(in),dimension(nLevels) :: tau, g, w0
    ! OUTPUTS
    real(wp)                               :: compute_toa_reflectace
    ! LOCAL VARIABLES
    logical, dimension(nLevels)                   :: cloudMask
    integer, dimension(count(tau(1:nLevels) > 0)) :: cloudIndicies
    real(wp),dimension(count(tau(1:nLevels) > 0)) :: Refl,Trans
    real(wp)                                      :: Refl_tot, Trans_tot
    integer                                       :: i

    cloudMask(1:nLevels) = tau(1:nLevels) > 0. 
    cloudIndicies = pack((/ (i, i = 1, nLevels) /), mask = cloudMask) 
    do i = 1, size(cloudIndicies)
       call two_stream(tau(cloudIndicies(i)), g(cloudIndicies(i)), w0(cloudIndicies(i)), Refl(i), Trans(i))
    end do
    
    call adding_doubling(count(tau(1:nLevels) > 0),Refl(:), Trans(:), Refl_tot, Trans_tot)  
    
    compute_toa_reflectace = Refl_tot
    
  end function compute_toa_reflectace
 
  ! ########################################################################################
  pure subroutine two_stream(tauint, gint, w0int, ref, tra) 
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    ! INPUTS
    real(wp), intent(in)  :: tauint, gint, w0int
    ! OUTPUTS
    real(wp), intent(out) :: ref, tra
    ! LOCAL VARIABLES
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    integer, parameter :: beam = 2
    real(wp),parameter :: xmu = 0.866, minConservativeW0 = 0.9999999
    real(wp) :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
         rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den, th
    
    ! Compute reflectance and transmittance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    f   = gint**2
    tau = (1._wp - w0int * f) * tauint
    w0  = (1._wp - f) * w0int / (1._wp - w0int * f)
    g   = (gint - f) / (1._wp - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7._wp - w0* (4._wp + 3._wp * g)) / 4._wp
    gamma2 = -(1._wp - w0* (4._wp - 3._wp * g)) / 4._wp
    gamma3 =  (2._wp - 3._wp*g*xmu) / 4._wp
    gamma4 =   1._wp - gamma3

    if (w0int > minConservativeW0) then
      ! Conservative scattering
      if (beam == 1) then
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-exp(-tau/xmu)))
  
          ref = rh / (1._wp + gamma1 * tau)
          tra = 1._wp - ref       
      else if(beam == 2) then
          ref = gamma1*tau/(1._wp + gamma1*tau)
          tra = 1._wp - ref
      endif
    else
      ! Non-conservative scattering
      a1 = gamma1 * gamma4 + gamma2 * gamma3
      a2 = gamma1 * gamma3 + gamma2 * gamma4

      rk = sqrt(gamma1**2 - gamma2**2)
      
      r1 = (1._wp - rk * xmu) * (a2 + rk * gamma3)
      r2 = (1._wp + rk * xmu) * (a2 - rk * gamma3)
      r3 = 2._wp * rk *(gamma3 - a2 * xmu)
      r4 = (1._wp - (rk * xmu)**2) * (rk + gamma1)
      r5 = (1._wp - (rk * xmu)**2) * (rk - gamma1)
      
      t1 = (1._wp + rk * xmu) * (a1 + rk * gamma4)
      t2 = (1._wp - rk * xmu) * (a1 - rk * gamma4)
      t3 = 2._wp * rk * (gamma4 + a1 * xmu)
      t4 = r4
      t5 = r5

      beta = -r5 / r4         
  
      e1 = min(rk * tau, 500._wp) 
      e2 = min(tau / xmu, 500._wp) 
      
      if (beam == 1) then
         den = r4 * exp(e1) + r5 * exp(-e1)
         ref  = w0*(r1*exp(e1)-r2*exp(-e1)-r3*exp(-e2))/den
         den = t4 * exp(e1) + t5 * exp(-e1)
         th  = exp(-e2)
         tra = th-th*w0*(t1*exp(e1)-t2*exp(-e1)-t3*exp(e2))/den
      elseif (beam == 2) then
         ef1 = exp(-e1)
         ef2 = exp(-2*e1)
         ref = (gamma2*(1._wp-ef2))/((rk+gamma1)*(1._wp-beta*ef2))
         tra = (2._wp*rk*ef1)/((rk+gamma1)*(1._wp-beta*ef2))
      endif
    end if
  end subroutine two_stream

  ! ########################################################################################
  elemental function two_stream_reflectance(tauint, gint, w0int)
    ! Compute reflectance in a single layer using the two stream approximation 
    !   The code itself is from Lazaros Oreopoulos via Steve Platnick 
    
    ! INPUTS
    real(wp), intent(in) :: tauint, gint, w0int
    ! OUTPUTS
    real(wp)             :: two_stream_reflectance
    ! LOCAL VARIABLES
    !   for delta Eddington code
    !   xmu, gamma3, and gamma4 only used for collimated beam approximation (i.e., beam=1)
    integer, parameter :: beam = 2
    real(wp),parameter :: xmu = 0.866, minConservativeW0 = 0.9999999
    real(wp) :: tau, w0, g, f, gamma1, gamma2, gamma3, gamma4, &
         rh, a1, a2, rk, r1, r2, r3, r4, r5, t1, t2, t3, t4, t5, beta, e1, e2, ef1, ef2, den

    f   = gint**2
    tau = (1._wp - w0int * f) * tauint
    w0  = (1._wp - f) * w0int / (1._wp - w0int * f)
    g   = (gint - f) / (1._wp - f)

    ! delta-Eddington (Joseph et al. 1976)
    gamma1 =  (7._wp - w0* (4._wp + 3._wp * g)) / 4._wp
    gamma2 = -(1._wp - w0* (4._wp - 3._wp * g)) / 4._wp
    gamma3 =  (2._wp - 3._wp*g*xmu) / 4._wp
    gamma4 =   1._wp - gamma3

    if (w0int > minConservativeW0) then
      ! Conservative scattering
      if (beam == 1) then
          rh = (gamma1*tau+(gamma3-gamma1*xmu)*(1-exp(-tau/xmu)))
          two_stream_reflectance = rh / (1._wp + gamma1 * tau)
      elseif (beam == 2) then
          two_stream_reflectance = gamma1*tau/(1._wp + gamma1*tau)
      endif
        
    else    !

        ! Non-conservative scattering
         a1 = gamma1 * gamma4 + gamma2 * gamma3
         a2 = gamma1 * gamma3 + gamma2 * gamma4

         rk = sqrt(gamma1**2 - gamma2**2)
         
         r1 = (1._wp - rk * xmu) * (a2 + rk * gamma3)
         r2 = (1._wp + rk * xmu) * (a2 - rk * gamma3)
         r3 = 2._wp * rk *(gamma3 - a2 * xmu)
         r4 = (1._wp - (rk * xmu)**2) * (rk + gamma1)
         r5 = (1._wp - (rk * xmu)**2) * (rk - gamma1)
         
         t1 = (1._wp + rk * xmu) * (a1 + rk * gamma4)
         t2 = (1._wp - rk * xmu) * (a1 - rk * gamma4)
         t3 = 2._wp * rk * (gamma4 + a1 * xmu)
         t4 = r4
         t5 = r5

         beta = -r5 / r4         
         
         e1 = min(rk * tau, 500._wp) 
         e2 = min(tau / xmu, 500._wp) 
         
         if (beam == 1) then
           den = r4 * exp(e1) + r5 * exp(-e1)
           two_stream_reflectance  = w0*(r1*exp(e1)-r2*exp(-e1)-r3*exp(-e2))/den
         elseif (beam == 2) then
           ef1 = exp(-e1)
           ef2 = exp(-2*e1)
           two_stream_reflectance = (gamma2*(1._wp-ef2))/((rk+gamma1)*(1._wp-beta*ef2))
         endif
           
      end if
  end function two_stream_reflectance 

  ! ########################################################################################
  pure subroutine adding_doubling (npts,Refl, Tran, Refl_tot, Tran_tot)      
    ! Use adding/doubling formulas to compute total reflectance and transmittance from 
    ! layer values
    
    ! INPUTS
    integer,intent(in)                  :: npts
    real(wp),intent(in),dimension(npts) :: Refl,Tran
    ! OUTPUTS
    real(wp),intent(out)                :: Refl_tot, Tran_tot
    ! LOCAL VARIABLES
    integer :: i
    real(wp), dimension(npts) :: Refl_cumulative, Tran_cumulative
    
    Refl_cumulative(1) = Refl(1)
    Tran_cumulative(1) = Tran(1)    
    
    do i=2, npts
       ! place (add) previous combined layer(s) reflectance on top of layer i, w/black surface (or ignoring surface):
       Refl_cumulative(i) = Refl_cumulative(i-1) + Refl(i)*(Tran_cumulative(i-1)**2)/(1._wp - Refl_cumulative(i-1) * Refl(i))
       Tran_cumulative(i) = (Tran_cumulative(i-1)*Tran(i)) / (1._wp - Refl_cumulative(i-1) * Refl(i))
    end do
    
    Refl_tot = Refl_cumulative(size(Refl))
    Tran_tot = Tran_cumulative(size(Refl))
    
  end subroutine adding_doubling

end module mod_modis_sim
