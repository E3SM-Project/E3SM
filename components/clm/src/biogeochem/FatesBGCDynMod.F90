module FatesBGCDynMod
   
   ! ==============================================================================================
   ! This module creates a pathway to call the belowground biogeochemistry 
   ! code as driven by the FATES vegetation model but bypassing the aboveground 
   ! CN vegetation code.  It is modeled after the CNDriverMod in its call sequence and 
   ! functionality.
   ! ==============================================================================================
   
   use shr_kind_mod, only : r8 => shr_kind_r8
   use perf_mod    , only : t_startf, t_stopf
   use shr_log_mod , only : errMsg => shr_log_errMsg
   use abortutils  , only : endrun
   
   implicit none

   public :: FatesBGCDyn
   
   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   
contains
   
   
   !-----------------------------------------------------------------------
   subroutine FatesBGCDyn(bounds,        &
         num_soilc, filter_soilc, num_soilp, filter_soilp, &
         carbonflux_vars, carbonstate_vars, cnstate_vars, &
         c13_carbonflux_vars, c13_carbonstate_vars,  &
         c14_carbonflux_vars, c14_carbonstate_vars,  &
         canopystate_vars, soilstate_vars, temperature_vars, &
         ch4_vars, nitrogenflux_vars, nitrogenstate_vars, &
         phosphorusstate_vars, phosphorusflux_vars, &
         alm_fates, crop_vars)
      
      use clm_varctl             , only : use_c13, use_c14, use_fates
      use decompMod              , only : bounds_type
      use clm_varpar             , only : nlevgrnd, nlevdecomp_full 
      use clm_varpar             , only : nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools 
      use clm_varctl             , only : use_century_decomp
      use DecompCascadeBGCMod  , only : decomp_rate_constants_bgc
      use DecompCascadeCNMod   , only : decomp_rate_constants_cn
      use CarbonStateUpdate1Mod     , only : CarbonStateUpdate1
      use SoilLittVertTranspMod, only : SoilLittVertTransp
      use PrecisionControlMod    , only : PrecisionControl
      use CNCarbonFluxType       , only : carbonflux_type
      use CNCarbonStateType      , only : carbonstate_type
      use CNStateType            , only : cnstate_type
      use CanopyStateType        , only : canopystate_type
      use SoilStateType          , only : soilstate_type
      use TemperatureType        , only : temperature_type
      use CNNitrogenFluxType     , only : nitrogenflux_type
      use CNNitrogenStateType    , only : nitrogenstate_type
      use CH4Mod                 , only : ch4_type
      use PhosphorusStateType    , only : phosphorusstate_type
      use PhosphorusFluxType     , only : phosphorusflux_type
      use CNDecompCascadeConType , only : decomp_cascade_con
      use CLMFatesInterfaceMod   , only : hlm_fates_interface_type
      use CropType               , only : crop_type

    !
    ! !ARGUMENTS:
    type(bounds_type)         , intent(in)    :: bounds  
    integer                   , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                   , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                   , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    type(carbonflux_type)     , intent(inout) :: carbonflux_vars
    type(carbonstate_type)    , intent(inout) :: carbonstate_vars
    type(cnstate_type)        , intent(inout) :: cnstate_vars
    type(carbonflux_type)     , intent(inout) :: c13_carbonflux_vars
    type(carbonstate_type)    , intent(inout) :: c13_carbonstate_vars
    type(carbonflux_type)     , intent(inout) :: c14_carbonflux_vars
    type(carbonstate_type)    , intent(inout) :: c14_carbonstate_vars
    type(canopystate_type)    , intent(in)    :: canopystate_vars
    type(soilstate_type)      , intent(in)    :: soilstate_vars
    type(temperature_type)    , intent(inout) :: temperature_vars 
    type(ch4_type)            , intent(in)    :: ch4_vars
    type(nitrogenflux_type)   , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)  , intent(inout) :: nitrogenstate_vars
    type(phosphorusstate_type), intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type) , intent(inout) :: phosphorusflux_vars
    type(hlm_fates_interface_type),intent(inout) :: alm_fates
    type(crop_type)          , intent(inout)  :: crop_vars
    
    !
    ! !LOCAL VARIABLES:
    integer :: k,j,fc,c
    real(r8):: p_decomp_cpool_loss(bounds%begc:bounds%endc, &
          1:nlevdecomp,1:ndecomp_cascade_transitions)       !potential C loss from one pool to another
    ! For methane code
    real(r8):: hrsum(bounds%begc:bounds%endc,1:nlevdecomp)  !sum of HR (gC/m2/s)
    !-----------------------------------------------------------------------

    associate( &
      ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)       
      phr_vr                      => carbonflux_vars%phr_vr_col, &       
      ! Input:  [real(r8) (:,:,:) ]  respired fraction in decomposition step (frac)
      rf_decomp_cascade           => cnstate_vars%rf_decomp_cascade_col, & 
      ! Input:  [real(r8) (:,:,:) ]
      ! vertically-resolved decomposing (litter, cwd, soil) c pools (gC/m3)
      decomp_cpools_vr            => carbonstate_vars%decomp_cpools_vr_col, & 
      ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec) 
      decomp_k                    => carbonflux_vars%decomp_k_col, &
      ! Input:  [real(r8) (:,:,:) ]  what fraction of 
      ! C leaving a given pool passes through a 
      ! given transition (frac)  
      pathfrac_decomp_cascade     => cnstate_vars%pathfrac_decomp_cascade_col, &  
      ! Output:  [real(r8) (:,:)   ]  fraction by which decomposition is limited by moisture availability
      w_scalar                    => carbonflux_vars%w_scalar_col, &   
      ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
      decomp_cascade_hr_vr        => carbonflux_vars%decomp_cascade_hr_vr_col, & 
      ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
      decomp_cascade_ctransfer_vr => carbonflux_vars%decomp_cascade_ctransfer_vr_col, & 
      ! Input:  [integer  (:)     ]  which pool is C taken from for a given 
      ! decomposition  step
      cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool, & 
      ! Input:  [integer  (:)     ]  which pool is C addedto for a given 
      ! decomposition step
      cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool )


    ! --------------------------------------------------
    ! zero the column-level C and N fluxes
    ! --------------------------------------------------
    
    call t_startf('BGCZero')

    call carbonflux_vars%SetValues(&
          num_soilp, filter_soilp, 0._r8, num_soilc, filter_soilc, 0._r8)
    if ( use_c13 ) then
       call c13_carbonflux_vars%SetValues(&
             num_soilp, filter_soilp, 0._r8, num_soilc, filter_soilc, 0._r8)
    end if
    if ( use_c14 ) then
       call c14_carbonflux_vars%SetValues(&
             num_soilp, filter_soilp, 0._r8, num_soilc, filter_soilc, 0._r8)
    end if

    call t_stopf('BGCZero')

    ! --------------------------------------------------
    ! Nitrogen Deposition, Fixation and Respiration
    ! --------------------------------------------------

    ! call t_startf('CNDeposition')
    ! call NitrogenDeposition(bounds, &
    !      atm2lnd_inst, soilbiogeochem_nitrogenflux_inst)
    ! call t_stopf('CNDeposition')
    ! if (crop_prog) then
    !    call NitrogenFert(bounds, num_soilc,filter_soilc, &
    !         cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenflux_inst)

    !    call  CNSoyfix (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
    !         waterstate_inst, crop_inst, cnveg_state_inst, cnveg_nitrogenflux_inst , &
    !         soilbiogeochem_state_inst, soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst)
    ! end if

    !--------------------------------------------
    ! Soil Biogeochemistry
    !--------------------------------------------
    if (use_century_decomp) then
       call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
               canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars, cnstate_vars)
    else
       call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc, &
             canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars, cnstate_vars)
    end if

    ! SoilBiogeochemPotential() in CLM
    ! Add up potential hr for methane calculations 

    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) &
                   * decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k)
          end do
       end do
      end do
    do j = 1,nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          phr_vr(c,j) = 0._r8
       end do
    end do
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             phr_vr(c,j) = phr_vr(c,j) + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
          end do
       end do
    end do


    !--------------------------------------------
    ! Resolve the competition between plants and soil heterotrophs 
    ! for available soil mineral N resource 
    !--------------------------------------------
    ! will add this back in when integrtating hte nutirent cycles


    !--------------------------------------------
    ! Calculate litter and soil decomposition rate
    !--------------------------------------------

    ! Calculation of actual immobilization and decomp rates, following
    ! resolution of plant/heterotroph  competition for mineral N (previously inlined in SoilLittDecompAllocation in SoilLittDecompMod)
    !  call SoilBiogeochemDecomp() in CLM
    call t_startf('SoilBiogeochemDecomp')

    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             !
             decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
             !
             decomp_cascade_ctransfer_vr(c,j,k) = (1._r8 - rf_decomp_cascade(c,j,k)) * p_decomp_cpool_loss(c,j,k)
             !
          end do
       end do
    end do
      
    call t_stopf('SoilBiogeochemDecomp')


    !--------------------------------------------
    ! Update1
    !--------------------------------------------

    call t_startf('BNGCUpdate1')


    ! -------------------------------------------------------
    ! Pass in FATES boundary conditions, ie litter fluxes
    ! -------------------------------------------------------
    
    call alm_fates%UpdateLitterFluxes(bounds,carbonflux_vars)

    ! Update all prognostic carbon state variables (except for gap-phase mortality and fire fluxes)

    call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            crop_vars, carbonflux_vars, carbonstate_vars)

    call t_stopf('BNGCUpdate1')

    !--------------------------------------------
    ! Calculate vertical mixing of soil and litter pools
    !--------------------------------------------

    call t_startf('SoilBiogeochemLittVertTransp')

    call SoilLittVertTransp(bounds, &
          num_soilc, filter_soilc, &
          canopystate_vars, cnstate_vars,                               &
          carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, &
          carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,    &
          nitrogenstate_vars, nitrogenflux_vars,&
          phosphorusstate_vars,phosphorusflux_vars)

    call t_stopf('SoilBiogeochemLittVertTransp')

    call t_startf('BGCsum')
    
    ! Set controls on very low values in critical state variables 
    ! Added some new logical filters to prevent
    ! above ground precision control calculations with use_fates, as well
    ! bypass on nitrogen calculations
    call PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp, &
          carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars,       &
          nitrogenstate_vars,phosphorusstate_vars)

    
    call FatesBGCSummary(bounds, num_soilc, filter_soilc,carbonflux_vars, carbonstate_vars)

    ! ----------------------------------------------
    ! calculate balance checks on entire carbon cycle (FATES + BGC)
    ! ----------------------------------------------

    call alm_fates%wrap_bgc_summary( bounds, carbonflux_vars, carbonstate_vars)

    call t_stopf('BGCsum')


    end associate

  end subroutine FatesBGCDyn

  ! ======================================================================================

  subroutine FatesBGCSummary(bounds, num_soilc, filter_soilc, carbonflux_vars, carbonstate_vars )

    use decompMod              , only : bounds_type
    use CNCarbonFluxType       , only : carbonflux_type
    use CNCarbonStateType      , only : carbonstate_type
    use clm_varcon             , only : dzsoi_decomp, zisoi
    use clm_varpar             , only : nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools 
    use CNDecompCascadeConType , only : decomp_cascade_con
    
    type(bounds_type)         , intent(in)    :: bounds  
    integer                   , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                   , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    type(carbonflux_type)     , intent(inout) :: carbonflux_vars
    type(carbonstate_type)    , intent(inout) :: carbonstate_vars

    real(r8) :: maxdepth
    integer  :: k,j,fc,c,l


    ! ------------------------------------------------------------------------------------
    !
    ! Summarize Carbon Fluxes
    !
    ! ------------------------------------------------------------------------------------

    do fc = 1,num_soilc
       c = filter_soilc(fc)
       carbonflux_vars%som_c_leached_col(c) = 0._r8
    end do
    
    ! vertically integrate HR and decomposition cascade fluxes
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             carbonflux_vars%decomp_cascade_hr_col(c,k) = &
                  carbonflux_vars%decomp_cascade_hr_col(c,k) + &
                  carbonflux_vars%decomp_cascade_hr_vr_col(c,j,k) * dzsoi_decomp(j) 
             
             carbonflux_vars%decomp_cascade_ctransfer_col(c,k) = &
                  carbonflux_vars%decomp_cascade_ctransfer_col(c,k) + &
                  carbonflux_vars%decomp_cascade_ctransfer_vr_col(c,j,k) * dzsoi_decomp(j) 
          end do
       end do
    end do
    
    ! total heterotrophic respiration, vertically resolved (HR)
    do j = 1,nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          carbonflux_vars%hr_vr_col(c,j) = 0._r8
       end do
    end do
    do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             carbonflux_vars%hr_vr_col(c,j) = &
                  carbonflux_vars%hr_vr_col(c,j) + &
                  carbonflux_vars%decomp_cascade_hr_vr_col(c,j,k)
          end do
       end do
    end do

    ! add up all vertical transport tendency terms and calculate total som leaching loss as the sum of these
    do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          carbonflux_vars%decomp_cpools_leached_col(c,l) = 0._r8
       end do
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             carbonflux_vars%decomp_cpools_leached_col(c,l) = carbonflux_vars%decomp_cpools_leached_col(c,l) + &
                  carbonflux_vars%decomp_cpools_transport_tendency_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          carbonflux_vars%som_c_leached_col(c) = carbonflux_vars%som_c_leached_col(c) + &
               carbonflux_vars%decomp_cpools_leached_col(c,l)
       end do
    end do

    ! soil organic matter heterotrophic respiration 
    associate(is_soil => decomp_cascade_con%is_soil) ! TRUE => pool is a soil pool  
      do k = 1, ndecomp_cascade_transitions
         if ( is_soil(decomp_cascade_con%cascade_donor_pool(k)) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               carbonflux_vars%somhr_col(c) = carbonflux_vars%somhr_col(c) + &
                    carbonflux_vars%decomp_cascade_hr_col(c,k)
            end do
         end if
      end do
    end associate

    ! litter heterotrophic respiration (LITHR)
    associate(is_litter => decomp_cascade_con%is_litter) ! TRUE => pool is a litter pool
      do k = 1, ndecomp_cascade_transitions
         if ( is_litter(decomp_cascade_con%cascade_donor_pool(k)) ) then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               carbonflux_vars%lithr_col(c) = carbonflux_vars%lithr_col(c) + &
                    carbonflux_vars%decomp_cascade_hr_col(c,k)
            end do
         end if
      end do
    end associate

    ! total heterotrophic respiration (HR)
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       
          carbonflux_vars%hr_col(c) = &
               carbonflux_vars%lithr_col(c) + &
               carbonflux_vars%somhr_col(c)
       
    end do


    ! ------------------------------------------------------------------------------------
    !
    ! Summarize Carbon States
    !
    ! ------------------------------------------------------------------------------------

    ! vertically integrate each of the decomposing C pools
    do l = 1, ndecomp_pools
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          carbonstate_vars%decomp_cpools_col(c,l) = 0._r8
       end do
    end do
    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             carbonstate_vars%decomp_cpools_col(c,l) = &
                  carbonstate_vars%decomp_cpools_col(c,l) + &
                  carbonstate_vars%decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
          end do
       end do
    end do

    if ( nlevdecomp > 1) then

       ! vertically integrate each of the decomposing C pools to 1 meter
       maxdepth = 1._r8
       do l = 1, ndecomp_pools
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             carbonstate_vars%decomp_cpools_1m_col(c,l) = 0._r8
          end do
       end do
       do l = 1, ndecomp_pools
          do j = 1, nlevdecomp
             if ( zisoi(j) <= maxdepth ) then
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   carbonstate_vars%decomp_cpools_1m_col(c,l) = &
                        carbonstate_vars%decomp_cpools_1m_col(c,l) + &
                        carbonstate_vars%decomp_cpools_vr_col(c,j,l) * dzsoi_decomp(j)
                end do
             elseif ( zisoi(j-1) < maxdepth ) then
                do fc = 1,num_soilc
                   c = filter_soilc(fc)
                   carbonstate_vars%decomp_cpools_1m_col(c,l) = &
                        carbonstate_vars%decomp_cpools_1m_col(c,l) + &
                        carbonstate_vars%decomp_cpools_vr_col(c,j,l) * (maxdepth - zisoi(j-1))
                end do
             endif
          end do
       end do

    endif

    ! truncation carbon
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       carbonstate_vars%ctrunc_col(c) = 0._r8
    end do
    do j = 1, nlevdecomp
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          carbonstate_vars%ctrunc_col(c) = &
               carbonstate_vars%ctrunc_col(c) + &
               carbonstate_vars%ctrunc_vr_col(c,j) * dzsoi_decomp(j)
       end do
    end do

    ! total litter carbon in the top meter (TOTLITC_1m)
    if ( nlevdecomp > 1) then
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          carbonstate_vars%totlitc_1m_col(c) = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_litter(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                carbonstate_vars%totlitc_1m_col(c) = carbonstate_vars%totlitc_1m_col(c) + &
                     carbonstate_vars%decomp_cpools_1m_col(c,l)
             end do
          endif
       end do
    end if

    ! total soil organic matter carbon in the top meter (TOTSOMC_1m)
    if ( nlevdecomp > 1) then
       do fc = 1,num_soilc
          c = filter_soilc(fc)
          carbonstate_vars%totsomc_1m_col(c) = 0._r8
       end do
       do l = 1, ndecomp_pools
          if ( decomp_cascade_con%is_soil(l) ) then
             do fc = 1,num_soilc
                c = filter_soilc(fc)
                carbonstate_vars%totsomc_1m_col(c) = carbonstate_vars%totsomc_1m_col(c) + &
                     carbonstate_vars%decomp_cpools_1m_col(c,l)
             end do
          end if
       end do
    end if

    ! total litter carbon (TOTLITC)
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       carbonstate_vars%totlitc_col(c) = 0._r8
    end do
    do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             carbonstate_vars%totlitc_col(c) = carbonstate_vars%totlitc_col(c) + &
                  carbonstate_vars%decomp_cpools_col(c,l)
          end do
       endif
    end do

    ! total soil organic matter carbon (TOTSOMC)
    do fc = 1,num_soilc
       c = filter_soilc(fc)
       carbonstate_vars%totsomc_col(c) = 0._r8
    end do
    do l = 1, ndecomp_pools
       if ( decomp_cascade_con%is_soil(l) ) then
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             carbonstate_vars%totsomc_col(c) = carbonstate_vars%totsomc_col(c) + &
                  carbonstate_vars%decomp_cpools_col(c,l)
          end do
       end if
    end do
    
    return
  end subroutine FatesBGCSummary


end  module FatesBGCDynMod
