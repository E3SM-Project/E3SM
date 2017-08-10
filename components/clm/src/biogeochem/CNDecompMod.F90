module CNDecompMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in litter and soil decomposition model
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_const_mod          , only : SHR_CONST_TKFRZ
  use decompMod              , only : bounds_type
  use perf_mod               , only : t_startf, t_stopf
  use clm_varctl             , only : iulog, use_nitrif_denitrif, use_lch4, use_century_decomp
  use clm_varcon             , only : dzsoi_decomp
  use clm_varpar             , only : nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
  use CNDecompCascadeCNMod   , only : decomp_rate_constants_cn
  use CNDecompCascadeBGCMod  , only : decomp_rate_constants_bgc
  use CNNitrifDenitrifMod    , only : nitrif_denitrif
  use CNVerticalProfileMod   , only : decomp_vertprofiles
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  use CNNitrogenStateType    , only : nitrogenstate_type
  use CNNitrogenFluxType     , only : nitrogenflux_type
  !!  add phosphorus  -X. YANG
  use PhosphorusStateType    , only : phosphorusstate_type
  use PhosphorusFluxType     , only : phosphorusflux_type
  use clm_varctl             , only : nu_com

  use CNCarbonStateType      , only : carbonstate_type
  use CNCarbonFluxType       , only : carbonflux_type
  use PhotosynthesisType     , only : photosyns_type
  use CanopyStateType        , only : canopystate_type
  use SoilStateType          , only : soilstate_type
  use TemperatureType        , only : temperature_type
  use WaterStateType         , only : waterstate_type
  use ch4Mod                 , only : ch4_type
  use cropType               , only : crop_type
  !! bgc interface & pflotran:
  use clm_varctl             , only : use_bgc_interface, use_pflotran, pf_cmode
  !
  implicit none
  save
  private :: CNvariables_nan4pf  !pflotran
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  public :: readCNDecompParams
  public :: CNDecompAlloc
  ! pflotran
  public :: CNDecompAlloc2
  !
  type, private :: CNDecompParamsType
     real(r8) :: dnp         !denitrification proportion
  end type CNDecompParamsType

  type(CNDecompParamsType),private ::  CNDecompParamsInst
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readCNDecompParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read parameters
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io
    use abortutils   , only: endrun
    use shr_log_mod  , only: errMsg => shr_log_errMsg

    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNDecompParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    tString='dnp'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNDecompParamsInst%dnp=tempr 

  end subroutine readCNDecompParams

!!-------------------------------------------------------------------------------------------------
  subroutine CNDecompAlloc (bounds, num_soilc, filter_soilc,    &
                num_soilp, filter_soilp,                        &
                canopystate_vars, soilstate_vars,               &
                temperature_vars, waterstate_vars,              &
                cnstate_vars, ch4_vars,                         &
                carbonstate_vars, carbonflux_vars,              &
                nitrogenstate_vars, nitrogenflux_vars,          &
                phosphorusstate_vars,phosphorusflux_vars)

    !!-----------------------------------------------------------------------------
    !! DESCRIPTION:
    !! Modified for bgc-interface (wgs): 9/12/2015
    !! clm-bgc soil Module, can be called through clm_bgc_interface
    !! ONLY includes SOM decomposition & nitrification/denitrification (if use_nitrif_denitrif)
    !! CNAllocaiton is divided into 3 subroutines:
    !! (1) CNAllocation1_PlantNPDemand  is called in CNEcosystemDynNoLeaching1
    !! (2) CNAllocation2_ResolveNPLimit is called in CNDecompAlloc (this subroutine)
    !! (3) CNAllocation3_PlantCNPAlloc  is called in CNDecompAlloc2
    !!-----------------------------------------------------------------------------

    ! !USES:
!    use CNAllocationMod , only: CNAllocation
    use CNAllocationMod , only: CNAllocation2_ResolveNPLimit !! Phase-2 of CNAllocation
    !
    ! !ARGUMENT:
    type(bounds_type)        , intent(in)    :: bounds   
    integer                  , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    integer                  , intent(in)    :: num_soilp          ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)    ! filter for soil patches
!    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
!    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
!    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    !! add phosphorus --
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
!    type(crop_type)          , intent(in)    :: crop_vars
    !
    ! !LOCAL VARIABLES:
    integer :: c,j,k,l,m                                                                               !indices
    integer :: fc                                                                                      !lake filter column index
    real(r8):: p_decomp_cpool_loss(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential C loss from one pool to another
    !real(r8):: pmnf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential mineral N flux, from one pool to another
    !real(r8):: pmpf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential mineral P flux, from one pool to another
    real(r8):: immob(bounds%begc:bounds%endc,1:nlevdecomp)                                             !potential N immobilization
    real(r8):: immob_p(bounds%begc:bounds%endc,1:nlevdecomp)                                             !potential P immobilization
    real(r8):: ratio                                                                                   !temporary variable
    real(r8):: dnp                                                                                     !denitrification proportion
    real(r8):: cn_decomp_pools(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8):: cp_decomp_pools(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    real(r8):: cp_decomp_pools_new(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)               !C:P ratio of new SOM
    integer, parameter :: i_atm = 0
    ! For methane code
    real(r8):: phr_vr(bounds%begc:bounds%endc,1:nlevdecomp)                                            !potential HR (gC/m3/s)
    real(r8):: hrsum(bounds%begc:bounds%endc,1:nlevdecomp)                                             !sum of HR (gC/m2/s)
    !-----------------------------------------------------------------------
   
    associate(                                                                                           &
         cascade_donor_pool               =>    decomp_cascade_con%cascade_donor_pool                  , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool            =>    decomp_cascade_con%cascade_receiver_pool               , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step
         floating_cn_ratio_decomp_pools   =>    decomp_cascade_con%floating_cn_ratio_decomp_pools      , & ! Input:  [logical  (:)     ]  TRUE => pool has fixed C:N ratio                   
         floating_cp_ratio_decomp_pools   =>    decomp_cascade_con%floating_cp_ratio_decomp_pools      , & ! Input:  [logical  (:)     ]  TRUE => pool has fixed C:P ratio                   
         initial_cn_ratio                 =>    decomp_cascade_con%initial_cn_ratio                    , & ! Input:  [real(r8) (:)     ]  c:n ratio for initialization of pools             
         initial_cp_ratio                 =>    decomp_cascade_con%initial_cp_ratio                    , & ! Input:  [real(r8) (:)     ]  c:p ratio for initialization of pools             
  
         is_cwd                           =>    decomp_cascade_con%is_cwd                              , &
!
         fpi_vr                           =>    cnstate_vars%fpi_vr_col                                , & ! Input:  [real(r8) (:,:)   ]  fraction of potential immobilization for N(no units) 
         fpi_p_vr                         =>    cnstate_vars%fpi_p_vr_col                              , & ! Input:  [real(r8) (:,:)   ]  fraction of potential immobilization for P (no units) 
         rf_decomp_cascade                =>    cnstate_vars%rf_decomp_cascade_col                     , & ! Input:  [real(r8) (:,:,:) ]  respired fraction in decomposition step (frac)
         pathfrac_decomp_cascade          =>    cnstate_vars%pathfrac_decomp_cascade_col               , & ! Input:  [real(r8) (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)

         decomp_npools_vr                 =>    nitrogenstate_vars%decomp_npools_vr_col                , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
         decomp_ppools_vr                 =>    phosphorusstate_vars%decomp_ppools_vr_col              , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) P pools

         decomp_cpools_vr                 =>    carbonstate_vars%decomp_cpools_vr_col                  , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools

         w_scalar                         =>    carbonflux_vars%w_scalar_col                           , & ! Input:  [real(r8) (:,:)   ]  fraction by which decomposition is limited by moisture availability
         
         decomp_cascade_ntransfer_vr      =>    nitrogenflux_vars%decomp_cascade_ntransfer_vr_col      , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
         decomp_cascade_sminn_flux_vr     =>    nitrogenflux_vars%decomp_cascade_sminn_flux_vr_col     , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
         potential_immob_vr               =>    nitrogenflux_vars%potential_immob_vr_col               , & ! Output: [real(r8) (:,:)   ]                                                  
         sminn_to_denit_decomp_cascade_vr =>    nitrogenflux_vars%sminn_to_denit_decomp_cascade_vr_col , & ! Output: [real(r8) (:,:,:) ] 
         gross_nmin_vr                    =>    nitrogenflux_vars%gross_nmin_vr_col                    , & ! Output: [real(r8) (:,:)   ]                                                  
         net_nmin_vr                      =>    nitrogenflux_vars%net_nmin_vr_col                      , & ! Output: [real(r8) (:,:)   ]                                                  
         gross_nmin                       =>    nitrogenflux_vars%gross_nmin_col                       , & ! Output: [real(r8) (:)     ]  gross rate of N mineralization (gN/m2/s)          
         net_nmin                         =>    nitrogenflux_vars%net_nmin_col                         , & ! Output: [real(r8) (:)     ]  net rate of N mineralization (gN/m2/s)            
         !!! add phosphorus  
         decomp_cascade_ptransfer_vr      =>    phosphorusflux_vars%decomp_cascade_ptransfer_vr_col    , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
         decomp_cascade_sminp_flux_vr     =>    phosphorusflux_vars%decomp_cascade_sminp_flux_vr_col   , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
         potential_immob_p_vr             =>    phosphorusflux_vars%potential_immob_p_vr_col           , & ! Output: [real(r8) (:,:)   ]
         gross_pmin_vr                    =>    phosphorusflux_vars%gross_pmin_vr_col                  , & ! Output: [real(r8) (:,:)   ]
         net_pmin_vr                      =>    phosphorusflux_vars%net_pmin_vr_col                    , & ! Output: [real(r8) (:,:)   ]
         gross_pmin                       =>    phosphorusflux_vars%gross_pmin_col                     , & ! Output: [real(r8) (:)     ]  gross rate of P mineralization (gP/m2/s)
         net_pmin                         =>    phosphorusflux_vars%net_pmin_col                       , & ! Output: [real(r8) (:)     ]  net rate of P mineralization (gP/m2/s)

         decomp_cascade_hr_vr             =>    carbonflux_vars%decomp_cascade_hr_vr_col               , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      =>    carbonflux_vars%decomp_cascade_ctransfer_vr_col        , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_k                         =>    carbonflux_vars%decomp_k_col                           , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)      
         phr_vr                           =>    carbonflux_vars%phr_vr_col                             , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)                           
         fphr                             =>    carbonflux_vars%fphr_col                               , & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic
         pmnf_decomp_cascade              =>    nitrogenflux_vars%pmnf_decomp_cascade                  , &
         pmpf_decomp_cascade              =>    phosphorusflux_vars%pmpf_decomp_cascade                , & 
         soil_n_immob_flux                =>    nitrogenflux_vars%soil_n_immob_flux                    , &
         soil_n_immob_flux_vr             =>    nitrogenflux_vars%soil_n_immob_flux_vr                 , &
         soil_n_grossmin_flux             =>    nitrogenflux_vars%soil_n_grossmin_flux                 , &
         soil_p_immob_flux                =>    phosphorusflux_vars%soil_p_immob_flux                  , &
         soil_p_immob_flux_vr             =>    phosphorusflux_vars%soil_p_immob_flux_vr               , &
         soil_p_grossmin_flux             =>    phosphorusflux_vars%soil_p_grossmin_flux               , &
         actual_immob_vr                  =>    nitrogenflux_vars%actual_immob_vr_col                  , &
         actual_immob_p_vr                =>    phosphorusflux_vars%actual_immob_p_vr_col                &
         )

!!-------------------------------------------------------------------------------------------------
      ! calculate decomposition rates (originally called in CNEcosystemDynNoLeaching1)
      if (use_century_decomp) then
          call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
               canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars)
      else
          call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc, &
               canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars, cnstate_vars)
      end if
!!-------------------------------------------------------------------------------------------------

      ! set initial values for potential C and N fluxes
      p_decomp_cpool_loss(bounds%begc : bounds%endc, :, :) = 0._r8
      pmnf_decomp_cascade(bounds%begc : bounds%endc, :, :) = 0._r8
      pmpf_decomp_cascade(bounds%begc : bounds%endc, :, :) = 0._r8    !! initial values for potential P fluxes

      ! column loop to calculate potential decomp rates and total immobilization
      ! demand.

      !! calculate c:n ratios of applicable pools
      do l = 1, ndecomp_pools
         if ( floating_cn_ratio_decomp_pools(l) ) then
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if ( decomp_npools_vr(c,j,l) > 0._r8 ) then
                     cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_npools_vr(c,j,l)
                  end if
               end do
            end do
         else
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
               end do
            end do
         end if
      end do

      !! calculate c:p ratios of applicable pools
      do l = 1, ndecomp_pools
         if ( floating_cp_ratio_decomp_pools(l) ) then
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if ( decomp_ppools_vr(c,j,l) > 0._r8 ) then
                     cp_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_ppools_vr(c,j,l)
                  end if
               end do
            end do
         else
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  cp_decomp_pools(c,j,l) = initial_cp_ratio(l)
               end do
            end do
         end if
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               cp_decomp_pools_new(c,j,l) = initial_cp_ratio(l)
            end do
         end do
      end do

      ! calculate the non-nitrogen-limited fluxes
      ! these fluxes include the  "/ dt" term to put them on a
      ! per second basis, since the rate constants have been
      ! calculated on a per timestep basis.

      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. &
                    decomp_k(c,j,cascade_donor_pool(k)) > 0._r8 ) then
                  p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) &
                       * decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k)
                  if ( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)) ) then  !! not transition of cwd to litter

                     if (cascade_receiver_pool(k) /= i_atm ) then  ! not 100% respiration
                        ratio = 0._r8

                        if (decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._r8) then
                           ratio = cn_decomp_pools(c,j,cascade_receiver_pool(k))/cn_decomp_pools(c,j,cascade_donor_pool(k))
                        endif

                        pmnf_decomp_cascade(c,j,k) = (p_decomp_cpool_loss(c,j,k) * (1.0_r8 - rf_decomp_cascade(c,j,k) - ratio) &
                             / cn_decomp_pools(c,j,cascade_receiver_pool(k)) )

                     else   ! 100% respiration
                        pmnf_decomp_cascade(c,j,k) = - p_decomp_cpool_loss(c,j,k) / cn_decomp_pools(c,j,cascade_donor_pool(k))
                     endif

                  else   ! CWD -> litter
                     pmnf_decomp_cascade(c,j,k) = 0._r8
                  end if
                  !!! add phosphorus fluxes
                  if ( .not. is_cwd( (cascade_receiver_pool(k)) ) ) then  !! not transition of cwd to litter

                     if (cascade_receiver_pool(k) /= i_atm ) then  ! not 100% respiration
                        ratio = 0._r8

                        if (decomp_ppools_vr(c,j,cascade_donor_pool(k)) > 0._r8) then
!                           ratio = cp_decomp_pools(c,j,cascade_receiver_pool(k))/cp_decomp_pools(c,j,cascade_donor_pool(k))
                           ratio = cp_decomp_pools_new(c,j,cascade_receiver_pool(k))/cp_decomp_pools(c,j,cascade_donor_pool(k))
                        endif

                        pmpf_decomp_cascade(c,j,k) = (p_decomp_cpool_loss(c,j,k) * (1.0_r8 - rf_decomp_cascade(c,j,k) - ratio) &
                             / cp_decomp_pools_new(c,j,cascade_receiver_pool(k)) )

                     else   ! 100% respiration
                        pmpf_decomp_cascade(c,j,k) = - p_decomp_cpool_loss(c,j,k) / cp_decomp_pools(c,j,cascade_donor_pool(k))
                     endif

                  else   ! CWD -> litter
                     pmpf_decomp_cascade(c,j,k) = 0._r8
                  end if

               end if
            end do
         end do
      end do

      ! Sum up all the potential immobilization fluxes (positive pmnf flux)
      ! and all the mineralization fluxes (negative pmnf flux)
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            immob(c,j) = 0._r8
            immob_p(c,j) = 0._r8
         end do
      end do
      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if (pmnf_decomp_cascade(c,j,k) > 0._r8) then
                  immob(c,j) = immob(c,j) + pmnf_decomp_cascade(c,j,k)
               else
                  gross_nmin_vr(c,j) = gross_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
               end if
               if (pmpf_decomp_cascade(c,j,k) > 0._r8) then
                  immob_p(c,j) = immob_p(c,j) + pmpf_decomp_cascade(c,j,k)
               else
                  gross_pmin_vr(c,j) = gross_pmin_vr(c,j) - pmpf_decomp_cascade(c,j,k)
               end if
            end do
         end do
      end do

      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            potential_immob_vr(c,j) = immob(c,j)
            potential_immob_p_vr(c,j) = immob_p(c,j)
         end do
      end do

      ! Add up potential hr for methane calculations
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


!!-------------------------------------------------------------------------------------------------
!! 'decomp_vertprofiles' (calc nfixation_prof) is moved to CNEcosystemDynNoLeaching1
!! 'nfixation_prof' is used to 'calc_nuptake_prof' & 'calc_puptake_prof', which are called in CNAllocation1,2,3
!      call decomp_vertprofiles(bounds, &
!           num_soilc, filter_soilc, num_soilp, filter_soilp, &
!           soilstate_vars, canopystate_vars, cnstate_vars)
!!-------------------------------------------------------------------------------------------------
      
      if (use_nitrif_denitrif) then ! calculate nitrification and denitrification rates
         call nitrif_denitrif(bounds, &
              num_soilc, filter_soilc, &
              soilstate_vars, waterstate_vars, temperature_vars, ch4_vars, &
              carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars)
      end if

      ! now that potential N immobilization is known, call allocation
      ! to resolve the competition between plants and soil heterotrophs
      ! for available soil mineral N resource.
      ! in addition, calculate fpi_vr, fpi_p_vr, & fgp
      call t_startf('CNAllocation - phase-2')
      call CNAllocation2_ResolveNPLimit(bounds,                     &
               num_soilc, filter_soilc, num_soilp, filter_soilp,    &
               cnstate_vars,                                        &
               carbonstate_vars, carbonflux_vars,                   &
               nitrogenstate_vars, nitrogenflux_vars,               &
               phosphorusstate_vars,phosphorusflux_vars,            &
               soilstate_vars,waterstate_vars)
      call t_stopf('CNAllocation - phase-2')

      
      ! column loop to calculate actual immobilization and decomp rates, following
      ! resolution of plant/heterotroph  competition for mineral N

!!-------------------------------------------------------------------------------------------------
!! comment out c:n,c:p ratios calculation, they have been calculated at the beginning of this subroutine
!!-------------------------------------------------------------------------------------------------
      ! calculate c:n ratios of applicable pools
!      do l = 1, ndecomp_pools
!         if ( floating_cn_ratio_decomp_pools(l) ) then
!            do j = 1,nlevdecomp
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  if ( decomp_npools_vr(c,j,l) > 0._r8 ) then
!                     cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_npools_vr(c,j,l)
!                  end if
!               end do
!            end do
!         else
!            do j = 1,nlevdecomp
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
!               end do
!            end do
!         end if
!      end do
!
!
!      !! calculate c:p ratios of applicable pools
!      do l = 1, ndecomp_pools
!         if ( floating_cp_ratio_decomp_pools(l) ) then
!            do j = 1,nlevdecomp
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  if ( decomp_ppools_vr(c,j,l) > 0._r8 ) then
!                     cp_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_ppools_vr(c,j,l)
!                  end if
!               end do
!            end do
!         else
!            do j = 1,nlevdecomp
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!                  cp_decomp_pools(c,j,l) = initial_cp_ratio(l)
!               end do
!            end do
!         end if
!      end do
!!-------------------------------------------------------------------------------------------------

      ! upon return from CNAllocation, the fraction of potential immobilization
      ! has been set (cnstate_vars%fpi_vr_col). now finish the decomp calculations.
      ! Only the immobilization steps are limited by fpi_vr (pmnf > 0)
      ! Also calculate denitrification losses as a simple proportion
      ! of mineralization flux.
       
      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._r8) then
                  if ( pmnf_decomp_cascade(c,j,k) > 0._r8 .and. pmpf_decomp_cascade(c,j,k) > 0._r8 ) then    ! N and P co-limitation
                     p_decomp_cpool_loss(c,j,k) = p_decomp_cpool_loss(c,j,k) * min( fpi_vr(c,j),fpi_p_vr(c,j) )
                     pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * min( fpi_vr(c,j),fpi_p_vr(c,j) ) 
                     pmpf_decomp_cascade(c,j,k) = pmpf_decomp_cascade(c,j,k) * min( fpi_vr(c,j),fpi_p_vr(c,j) )   !!! immobilization step
                     if (.not. use_nitrif_denitrif) then
                        sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
                     end if
                  elseif ( pmnf_decomp_cascade(c,j,k) > 0._r8 .and. pmpf_decomp_cascade(c,j,k) <= 0._r8 ) then  ! N limitation 
                     p_decomp_cpool_loss(c,j,k) = p_decomp_cpool_loss(c,j,k) * fpi_vr(c,j)
                     pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * fpi_vr(c,j)
                     pmpf_decomp_cascade(c,j,k) = pmpf_decomp_cascade(c,j,k) * fpi_vr(c,j) !!! immobilization step
                     if (.not. use_nitrif_denitrif) then
                        sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
                     end if
                  elseif ( pmnf_decomp_cascade(c,j,k) < 0._r8 .and. pmpf_decomp_cascade(c,j,k) >  0._r8 ) then  ! P limitation 
                     p_decomp_cpool_loss(c,j,k) = p_decomp_cpool_loss(c,j,k) * fpi_p_vr(c,j)
                     pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * fpi_p_vr(c,j)
                     pmpf_decomp_cascade(c,j,k) = pmpf_decomp_cascade(c,j,k) * fpi_p_vr(c,j) !!! immobilization step

                     if (.not. use_nitrif_denitrif) then
                        sminn_to_denit_decomp_cascade_vr(c,j,k) = -CNDecompParamsInst%dnp * pmnf_decomp_cascade(c,j,k)
                     end if
                  elseif ( pmnf_decomp_cascade(c,j,k) < 0._r8 .and. pmpf_decomp_cascade(c,j,k) <=  0._r8 ) then  ! No limitation 
                     if (.not. use_nitrif_denitrif) then
                        sminn_to_denit_decomp_cascade_vr(c,j,k) = -CNDecompParamsInst%dnp * pmnf_decomp_cascade(c,j,k)
                     end if

                  end if
                  decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
                  decomp_cascade_ctransfer_vr(c,j,k) = (1._r8 - rf_decomp_cascade(c,j,k)) * p_decomp_cpool_loss(c,j,k)
                  if (decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. cascade_receiver_pool(k) /= i_atm) then
                     decomp_cascade_ntransfer_vr(c,j,k) = p_decomp_cpool_loss(c,j,k) / cn_decomp_pools(c,j,cascade_donor_pool(k))
                  else
                     decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
                  endif
                  !!! phosphorus fluxes
                  if (decomp_ppools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. cascade_receiver_pool(k) /= i_atm) then
                     decomp_cascade_ptransfer_vr(c,j,k) = p_decomp_cpool_loss(c,j,k) / cp_decomp_pools(c,j,cascade_donor_pool(k))
                  else
                     decomp_cascade_ptransfer_vr(c,j,k) = 0._r8
                  endif
                  if ( cascade_receiver_pool(k) /= 0 ) then
                     decomp_cascade_sminn_flux_vr(c,j,k) = pmnf_decomp_cascade(c,j,k)
                     decomp_cascade_sminp_flux_vr(c,j,k) = pmpf_decomp_cascade(c,j,k)
                  else  ! keep sign convention negative for terminal pools
                     decomp_cascade_sminn_flux_vr(c,j,k) = - pmnf_decomp_cascade(c,j,k)
                     decomp_cascade_sminp_flux_vr(c,j,k) = - pmpf_decomp_cascade(c,j,k)
                  endif
                  net_nmin_vr(c,j) = net_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
                  net_pmin_vr(c,j) = net_pmin_vr(c,j) - pmpf_decomp_cascade(c,j,k)
               else
                  decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
                  decomp_cascade_ptransfer_vr(c,j,k) = 0._r8
                  if (.not. use_nitrif_denitrif) then
                     sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
                  end if
                  decomp_cascade_sminn_flux_vr(c,j,k) = 0._r8
                  decomp_cascade_sminp_flux_vr(c,j,k) = 0._r8
               end if

            end do
         end do
      end do
     
      if (nu_com .ne. 'RD') then
      do fc = 1,num_soilc
          c = filter_soilc(fc)
          soil_n_immob_flux(c) =0.0_r8
          soil_p_immob_flux(c) = 0.0_r8
          soil_n_grossmin_flux(c) = 0.0_r8
          soil_p_grossmin_flux(c) = 0.0_r8
          do j = 1,nlevdecomp
              soil_n_immob_flux_vr(c,j) = 0.0_r8
              soil_p_immob_flux_vr(c,j) = 0.0_r8
              gross_nmin_vr(c,j) = 0.0_r8
              gross_pmin_vr(c,j) = 0.0_r8
          end do
      end do
      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
          	   if (pmnf_decomp_cascade(c,j,k) > 0._r8) then 
                   soil_n_immob_flux(c) = soil_n_immob_flux(c) + pmnf_decomp_cascade(c,j,k)*dzsoi_decomp(j)
                   soil_n_immob_flux_vr(c,j) = soil_n_immob_flux_vr(c,j) + pmnf_decomp_cascade(c,j,k)
               else
                   soil_n_grossmin_flux(c) = soil_n_grossmin_flux(c) -1.0_r8*pmnf_decomp_cascade(c,j,k)*dzsoi_decomp(j)
                   gross_nmin_vr(c,j) = gross_nmin_vr(c,j) - 1.0_r8*pmnf_decomp_cascade(c,j,k)
               end if
               if (pmpf_decomp_cascade(c,j,k) > 0._r8) then 
                   soil_p_immob_flux(c) = soil_p_immob_flux(c) + pmpf_decomp_cascade(c,j,k)*dzsoi_decomp(j)
                   soil_p_immob_flux_vr(c,j) = soil_p_immob_flux_vr(c,j) + pmpf_decomp_cascade(c,j,k)
               else
                   soil_p_grossmin_flux(c) = soil_p_grossmin_flux(c) -1.0_r8*pmpf_decomp_cascade(c,j,k)*dzsoi_decomp(j)
                   gross_pmin_vr(c,j) = gross_pmin_vr(c,j) - 1.0_r8*pmpf_decomp_cascade(c,j,k)
               end if
             end do
          end do
       end do
       end if

      if (use_lch4) then
         ! Calculate total fraction of potential HR, for methane code
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               hrsum(c,j) = 0._r8
            end do
         end do
         do k = 1, ndecomp_cascade_transitions
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  hrsum(c,j) = hrsum(c,j) + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
               end do
            end do
         end do

         ! Nitrogen limitation / (low)-moisture limitation
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if (phr_vr(c,j) > 0._r8) then
                  fphr(c,j) = hrsum(c,j) / phr_vr(c,j) * w_scalar(c,j)
                  fphr(c,j) = max(fphr(c,j), 0.01_r8) ! Prevent overflow errors for 0 respiration
               else
                  fphr(c,j) = 1._r8
               end if
            end do
         end do
      end if

      ! vertically integrate net and gross mineralization fluxes for diagnostic output
      ! moved to CNDecompAlloc2
!      do j = 1,nlevdecomp
!         do fc = 1,num_soilc
!            c = filter_soilc(fc)
!            net_nmin(c) = net_nmin(c) + net_nmin_vr(c,j) * dzsoi_decomp(j)
!            gross_nmin(c) = gross_nmin(c) + gross_nmin_vr(c,j) * dzsoi_decomp(j)
!            net_pmin(c) = net_pmin(c) + net_pmin_vr(c,j) * dzsoi_decomp(j)
!            gross_pmin(c) = gross_pmin(c) + gross_pmin_vr(c,j) * dzsoi_decomp(j)
!         end do
!      end do

    end associate

  end subroutine CNDecompAlloc

!!-------------------------------------------------------------------------------------------------

  subroutine CNDecompAlloc2 (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,   &
       photosyns_vars, canopystate_vars, soilstate_vars, temperature_vars,               &
       waterstate_vars, cnstate_vars, ch4_vars,                                          &
       carbonstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,      &
       nitrogenstate_vars, nitrogenflux_vars, crop_vars, atm2lnd_vars,                   &
       phosphorusstate_vars,phosphorusflux_vars)
    !!-----------------------------------------------------------------------------
    !! DESCRIPTION:
    !! bgc interface & pflotran:
    !! (1) Simplified codes of CNDecompAlloc subroutine for coupling with pflotran
    !! (2) call CNAllocation3_PlantCNPAlloc
    !! (3) calculate net_nmin(c), gross_nmin(c), net_pmin(c), gross_pmin(c)
    !!-----------------------------------------------------------------------------

    ! !USES:
    use CNAllocationMod , only: CNAllocation3_PlantCNPAlloc !! Phase-3 of CNAllocation
    use atm2lndType     , only: atm2lnd_type
    use clm_time_manager, only: get_step_size
!    use clm_varpar      , only: nlevdecomp, ndecomp_pools

    !
    ! !ARGUMENT:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    integer                  , intent(in)    :: num_soilp          ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)    ! filter for soil patches
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars
    !! add phosphorus --
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, j                                     ! indices
!    real(r8):: col_plant_ndemand(bounds%begc:bounds%endc)   ! column-level vertically-integrated plant N demand (gN/m2/s)
    real(r8):: dt                                           ! time step (seconds)

    ! For methane code
    real(r8):: hrsum(bounds%begc:bounds%endc,1:nlevdecomp)                                             !sum of HR (gC/m2/s)
    !-----------------------------------------------------------------------

    associate(                                                                                      &
         gross_nmin_vr                    =>    nitrogenflux_vars%gross_nmin_vr_col                    , & ! Output: [real(r8) (:,:)   ]
         net_nmin_vr                      =>    nitrogenflux_vars%net_nmin_vr_col                      , & ! Output: [real(r8) (:,:)   ]
         gross_nmin                       =>    nitrogenflux_vars%gross_nmin_col                       , & ! Output: [real(r8) (:)     ]  gross rate of N mineralization (gN/m2/s)
         net_nmin                         =>    nitrogenflux_vars%net_nmin_col                         , & ! Output: [real(r8) (:)     ]  net rate of N mineralization (gN/m2/s)
         !phosphorus
         gross_pmin_vr                    =>    phosphorusflux_vars%gross_pmin_vr_col                  , & ! Output: [real(r8) (:,:)   ]
         net_pmin_vr                      =>    phosphorusflux_vars%net_pmin_vr_col                    , & ! Output: [real(r8) (:,:)   ]
         gross_pmin                       =>    phosphorusflux_vars%gross_pmin_col                     , & ! Output: [real(r8) (:)     ]  gross rate of P mineralization (gP/m2/s)
         net_pmin                         =>    phosphorusflux_vars%net_pmin_col                       , & ! Output: [real(r8) (:)     ]  net rate of P mineralization (gP/m2/s)


         fpi_vr                           =>    cnstate_vars%fpi_vr_col                                , & ! Output:  [real(r8) (:,:)   ]  fraction of potential immobilization (no units)
         potential_immob_vr               =>    nitrogenflux_vars%potential_immob_vr_col               , & ! Input:
         actual_immob_vr                  =>    nitrogenflux_vars%actual_immob_vr_col                  , & ! Input:

         fpg                              =>    cnstate_vars%fpg_col                                   , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         sminn_to_plant                   =>    nitrogenflux_vars%sminn_to_plant_col                   , & ! Output: [real(r8) (:)     ]  col N uptake (gN/m2/s)
         sminn_to_plant_vr                =>    nitrogenflux_vars%sminn_to_plant_vr_col                , & ! Input:  [real(r8) (:,:)    ]  vertically-resolved N uptake (gN/m3/s)
         col_plant_ndemand_vr             =>    nitrogenflux_vars%plant_ndemand_vr_col                 , & ! Input:  [real(r8) (:)     ]  col N uptake (gN/m2/s)

         plant_ndemand_col                =>    nitrogenflux_vars%plant_ndemand_col                    , & ! Output:  [real(r8) (:,:) ]
         plant_pdemand_col                =>    phosphorusflux_vars%plant_pdemand_col                  , & ! Output:  [real(r8) (:,:) ]

         w_scalar                         =>    carbonflux_vars%w_scalar_col                           , & ! Input:  [real(r8) (:,:)   ]  fraction by which decomposition is limited by moisture availability
         decomp_cascade_hr_vr             =>    carbonflux_vars%decomp_cascade_hr_vr_col               , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      =>    carbonflux_vars%decomp_cascade_ctransfer_vr_col        , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_k                         =>    carbonflux_vars%decomp_k_col                           , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         hr_vr                            =>    carbonflux_vars%hr_vr_col                              , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
         phr_vr                           =>    carbonflux_vars%phr_vr_col                             , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
         fphr                             =>    carbonflux_vars%fphr_col                                 & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic
         )

      ! set time steps
      dt = real( get_step_size(), r8 )

      ! MUST have already updated needed bgc variables from PFLOTRAN by this point
!!------------------------------------------------------------------
! moved to EcosystemDynNoLeaching1
!      call decomp_vertprofiles(bounds, &
!           num_soilc, filter_soilc, num_soilp, filter_soilp, &
!           soilstate_vars, canopystate_vars, cnstate_vars)

!!------------------------------------------------------------------
    if(use_bgc_interface.and.use_pflotran.and.pf_cmode) then
      ! fpg calculation
      do fc=1,num_soilc
         c = filter_soilc(fc)
         sminn_to_plant(c)       = 0._r8
         !! local var 'col_plant_ndemand' is replaced by fluxtype%plant_ndemand_col
!         col_plant_ndemand(c)    = 0._r8
         do j = 1, nlevdecomp           ! sum up actual and potential column-level N fluxes to plant
            sminn_to_plant(c)    = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
!            col_plant_ndemand(c) = col_plant_ndemand(c) + col_plant_ndemand_vr(c,j) * dzsoi_decomp(j)
         end do
      end do
      do fc=1,num_soilc
         c = filter_soilc(fc)
         ! calculate the fraction of potential growth that can be
         ! acheived with the N available to plants
         if (plant_ndemand_col(c) > 0.0_r8) then
            fpg(c) = max(0._r8,sminn_to_plant(c)) / plant_ndemand_col(c)
            fpg(c) = min(1._r8, fpg(c))
         else
            fpg(c) = 1.0_r8
         end if
      end do

      ! fpi calculation
      do fc=1,num_soilc
         c = filter_soilc(fc)
         do j = 1, nlevdecomp
            if (potential_immob_vr(c,j) > 0.0_r8) then
                fpi_vr(c,j) = actual_immob_vr(c,j) / potential_immob_vr(c,j)
            else
                fpi_vr(c,j) = 0.0_r8
            end if
         end do

      end do

      if (use_lch4) then
         ! Add up potential hr for methane calculations
         ! potential hr is not available from PFLOTRAN, so here temporarily as actual hr
         ! in the end, this methane module will be moving into PFLOTRAN as well

         do j = 1,nlevdecomp
            do fc = 1, num_soilc
               c = filter_soilc(fc)
               phr_vr(c,j) = hr_vr(c,j)
            end do
         end do

         ! Calculate total fraction of potential HR, for methane code
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               hrsum(c,j) = hr_vr(c,j)
            end do
         end do

         ! Nitrogen limitation / (low)-moisture limitation
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if (phr_vr(c,j) > 0._r8) then
                  fphr(c,j) = hrsum(c,j) / phr_vr(c,j) * w_scalar(c,j)
                  fphr(c,j) = max(fphr(c,j), 0.01_r8) ! Prevent overflow errors for 0 respiration
               else
                  fphr(c,j) = 1._r8
               end if
            end do
         end do
      end if

      ! needs to zero CLM-CN variables NOT available from pflotran bgc coupling
      call CNvariables_nan4pf(bounds, num_soilc, filter_soilc,   &
                        carbonflux_vars, nitrogenflux_vars)
    end if !!if(use_bgc_interface.and.use_pflotran.and.pf_cmode)

!!------------------------------------------------------------------
      ! phase-3 Allocation for plants
      call t_startf('CNAllocation - phase-3')
      call CNAllocation3_PlantCNPAlloc (bounds                      , &
                num_soilc, filter_soilc, num_soilp, filter_soilp    , &
                canopystate_vars                                    , &
                cnstate_vars, carbonstate_vars, carbonflux_vars     , &
                c13_carbonflux_vars, c14_carbonflux_vars            , &
                nitrogenstate_vars, nitrogenflux_vars               , &
                phosphorusstate_vars, phosphorusflux_vars)
      call t_stopf('CNAllocation - phase-3')
!!------------------------------------------------------------------

      ! vertically integrate net and gross mineralization fluxes for diagnostic output
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            net_nmin(c) = net_nmin(c) + net_nmin_vr(c,j) * dzsoi_decomp(j)
            gross_nmin(c) = gross_nmin(c) + gross_nmin_vr(c,j) * dzsoi_decomp(j)
            !phosphorus
            net_pmin(c) = net_pmin(c) + net_pmin_vr(c,j) * dzsoi_decomp(j)
            gross_pmin(c) = gross_pmin(c) + gross_pmin_vr(c,j) * dzsoi_decomp(j)
         end do
      end do
    end associate
  end subroutine CNDecompAlloc2
 
!!-------------------------------------------------------------------------------------------------
  !
  subroutine CNvariables_nan4pf (bounds, num_soilc, filter_soilc, &
            carbonflux_vars, nitrogenflux_vars)
  !
  !DESCRIPTION:
  !  CN variables not available from PFLOTRAN, some of which may be output and may cause issues,
  !  if not properly set.
  !
  !USES:

    use clm_varpar   , only: nlevdecomp, ndecomp_cascade_transitions
   !
   !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
   !
   !CALLED FROM:
   !
   !LOCAL VARIABLES:
   integer :: c,j,k,fc          !indices
   !
   !-----------------------------------------------------------------------
   associate (&
         decomp_cascade_hr_vr             =>    carbonflux_vars%decomp_cascade_hr_vr_col               , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      =>    carbonflux_vars%decomp_cascade_ctransfer_vr_col        , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ntransfer_vr      =>    nitrogenflux_vars%decomp_cascade_ntransfer_vr_col        & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
   )
   ! set zeros for those variables NOT available from PFLOTRAN
   ! (TODO) will check 'zero' or '-9999' is better.

   do fc = 1,num_soilc
      c = filter_soilc(fc)
      do j = 1,nlevdecomp
         do k = 1, ndecomp_cascade_transitions
            decomp_cascade_ctransfer_vr(c,j,k) = 0._r8
            decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
            decomp_cascade_hr_vr(c,j,k)        = 0._r8
         end do
      end do

   end do

 end associate
 end subroutine CNvariables_nan4pf

!!-------------------------------------------------------------------------------------------------

end module CNDecompMod
