module CNDecompMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in litter and soil decomposition model
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod         , only: r8 => shr_kind_r8
  use shr_const_mod        , only: SHR_CONST_TKFRZ
  use clm_varctl           , only: iulog, use_nitrif_denitrif, use_lch4, use_century_decomp
  use clm_varcon           , only: dzsoi_decomp
  use CNDecompCascadeCNMod , only: decomp_rate_constants_cn
  use CNDecompCascadeBGCMod, only: decomp_rate_constants_bgc
  use CNNitrifDenitrifMod  , only: nitrif_denitrif
  use CNVerticalProfileMod , only: decomp_vertprofiles
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNDecompAlloc
  public :: readCNDecompParams
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

  !-----------------------------------------------------------------------
  subroutine CNDecompAlloc (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !USES:
    use clmtype
    use CNAllocationMod , only: CNAllocation
    use clm_varpar, only: nlevdecomp,ndecomp_cascade_transitions,ndecomp_pools
    use decompMod , only: bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds   ! bounds
    integer, intent(in) :: num_soilc          ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:)    ! filter for soil columns
    integer, intent(in) :: num_soilp          ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:)    ! filter for soil pfts
    !
    ! !LOCAL VARIABLES:
    integer :: c,j,k,l,m          !indices
    integer :: fc           !lake filter column index
    real(r8):: p_decomp_cpool_loss(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential C loss from one pool to another
    real(r8):: pmnf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential mineral N flux, from one pool to another
    real(r8):: immob(bounds%begc:bounds%endc,1:nlevdecomp)        !potential N immobilization
    real(r8):: ratio        !temporary variable
    real(r8):: dnp          !denitrification proportion
    real(r8):: cn_decomp_pools(bounds%begc:bounds%endc,1:nlevdecomp,1:ndecomp_pools)
    integer, parameter :: i_atm = 0
    ! For methane code
    real(r8):: phr_vr(bounds%begc:bounds%endc,1:nlevdecomp)       !potential HR (gC/m3/s)
    real(r8):: hrsum(bounds%begc:bounds%endc,1:nlevdecomp)        !sum of HR (gC/m2/s)
    !-----------------------------------------------------------------------
   
   associate(& 
   decomp_cpools_vr                    =>    ccs%decomp_cpools_vr                        , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools
   decomp_cascade_hr_vr                =>    ccf%decomp_cascade_hr_vr                    , & ! InOut:  [real(r8) (:,:,:)]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   decomp_cascade_ctransfer_vr         =>    ccf%decomp_cascade_ctransfer_vr             , & ! InOut:  [real(r8) (:,:,:)]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
   decomp_npools_vr                    =>    cns%decomp_npools_vr                        , & ! Input:  [real(r8) (:,:,:)]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   decomp_cascade_ntransfer_vr         =>    cnf%decomp_cascade_ntransfer_vr             , & ! InOut:  [real(r8) (:,:,:)]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
   decomp_cascade_sminn_flux_vr        =>    cnf%decomp_cascade_sminn_flux_vr            , & ! InOut:  [real(r8) (:,:,:)]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
   fpi_vr                              =>    cps%fpi_vr                                  , & ! InOut:  [real(r8) (:,:)]  fraction of potential immobilization (no units) 
   potential_immob_vr                  =>    cnf%potential_immob_vr                      , & ! InOut:  [real(r8) (:,:)]                                                  
   decomp_k                            =>    ccf%decomp_k                                , & ! InOut:  [real(r8) (:,:,:)]  rate constant for decomposition (1./sec)      
   rf_decomp_cascade                   =>    cps%rf_decomp_cascade                       , & ! InOut:  [real(r8) (:,:,:)]  respired fraction in decomposition step (frac)
   cascade_donor_pool                  =>    decomp_cascade_con%cascade_donor_pool       , & ! InOut:  [integer (:)]  which pool is C taken from for a given decomposition step
   cascade_receiver_pool               =>    decomp_cascade_con%cascade_receiver_pool    , & ! InOut:  [integer (:)]  which pool is C added to for a given decomposition step
   pathfrac_decomp_cascade             =>    cps%pathfrac_decomp_cascade                 , & ! InOut:  [real(r8) (:,:,:)]  what fraction of C leaving a given pool passes through a given transition (frac)
   floating_cn_ratio_decomp_pools      =>    decomp_cascade_con%floating_cn_ratio_decomp_pools  , & ! InOut:  [logical (:)]  TRUE => pool has fixed C:N ratio                   
   initial_cn_ratio                    =>    decomp_cascade_con%initial_cn_ratio         , & ! Output: [real(r8) (:)]  c:n ratio for initialization of pools             
   altmax_indx                         =>    cps%altmax_indx                             , & ! Output: [integer (:)]  maximum annual depth of thaw                       
   altmax_lastyear_indx                =>    cps%altmax_lastyear_indx                    , & ! Output: [integer (:)]  prior year maximum annual depth of thaw            
   sminn_to_denit_decomp_cascade_vr    =>    cnf%sminn_to_denit_decomp_cascade_vr        , & ! InOut:  [real(r8) (:,:,:)] 
   phr_vr                              =>    ccf%phr_vr                                  , & ! Output: [real(r8) (:,:)] potential HR (gC/m3/s)                           
   gross_nmin_vr                       =>    cnf%gross_nmin_vr                           , & ! InOut:  [real(r8) (:,:)]                                                  
   net_nmin_vr                         =>    cnf%net_nmin_vr                             , & ! InOut:  [real(r8) (:,:)]                                                  
   gross_nmin                          =>    cnf%gross_nmin                              , & ! InOut:  [real(r8) (:)]  gross rate of N mineralization (gN/m2/s)          
   net_nmin                            =>    cnf%net_nmin                                , & ! InOut:  [real(r8) (:)]  net rate of N mineralization (gN/m2/s)            
   fphr                                =>    cch4%fphr                                   , & ! InOut:  [real(r8) (:,:)]  fraction of potential SOM + LITTER heterotrophic respiration
   w_scalar                            =>    ccf%w_scalar                                , & ! InOut:  [real(r8) (:,:)]  fraction by which decomposition is limited by moisture availability
   rootfr                              =>    pps%rootfr                                  , & ! Input:  [real(r8) (:,:)]  fraction of roots in each soil layer  (nlevgrnd)
   clandunit                           =>   col%landunit                                 , & ! Input:  [integer (:)]  index into landunit level quantities               
   itypelun                            =>    lun%itype                                     & ! Input:  [integer (:)]  landunit type                                      
   )
   
   if (use_century_decomp) then
      call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc)
   else
      call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc)
   end if
   
   ! set initial values for potential C and N fluxes
   p_decomp_cpool_loss(bounds%begc : bounds%endc, :, :) = 0._r8
   pmnf_decomp_cascade(bounds%begc : bounds%endc, :, :) = 0._r8
   
   ! column loop to calculate potential decomp rates and total immobilization
   ! demand.
   
   
   !! calculate c:n ratios of applicable pools
   do l = 1, ndecomp_pools
      if ( floating_cn_ratio_decomp_pools(l) ) then
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if ( decomp_npools_vr(c,j,l) .gt. 0._r8 ) then
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
   
   ! calculate the non-nitrogen-limited fluxes
   ! these fluxes include the  "/ dt" term to put them on a
   ! per second basis, since the rate constants have been
   ! calculated on a per timestep basis.
   
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. &
               decomp_k(c,j,cascade_donor_pool(k)) .gt. 0._r8 ) then
                  p_decomp_cpool_loss(c,j,k) = decomp_cpools_vr(c,j,cascade_donor_pool(k)) &
                     * decomp_k(c,j,cascade_donor_pool(k))  * pathfrac_decomp_cascade(c,j,k)
               
               if ( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)) ) then  !! not transition of cwd to litter
                  
                  if (cascade_receiver_pool(k) .ne. i_atm ) then  ! not 100% respiration
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
         end do
      end do
   end do
   
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         potential_immob_vr(c,j) = immob(c,j)
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
   
   call decomp_vertprofiles(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
   
   if (use_nitrif_denitrif) then
      ! calculate nitrification and denitrification rates
      call nitrif_denitrif(bounds, num_soilc, filter_soilc)
   end if
   
   ! now that potential N immobilization is known, call allocation
   ! to resolve the competition between plants and soil heterotrophs
   ! for available soil mineral N resource.
   
   call CNAllocation(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
   
   ! column loop to calculate actual immobilization and decomp rates, following
   ! resolution of plant/heterotroph  competition for mineral N
   
   
   ! calculate c:n ratios of applicable pools
   do l = 1, ndecomp_pools
      if ( floating_cn_ratio_decomp_pools(l) ) then
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if ( decomp_npools_vr(c,j,l) .gt. 0._r8 ) then
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
   
   ! upon return from CNAllocation, the fraction of potential immobilization
   ! has been set (cps%fpi_vr). now finish the decomp calculations.
   ! Only the immobilization steps are limited by fpi_vr (pmnf > 0)
   ! Also calculate denitrification losses as a simple proportion
   ! of mineralization flux.
   
   do k = 1, ndecomp_cascade_transitions
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) .gt. 0._r8) then
               if ( pmnf_decomp_cascade(c,j,k) .gt. 0._r8 ) then
                  p_decomp_cpool_loss(c,j,k) = p_decomp_cpool_loss(c,j,k) * fpi_vr(c,j)
                  pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * fpi_vr(c,j)
                  if (.not. use_nitrif_denitrif) then
                     sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
                  end if
               else
                  if (.not. use_nitrif_denitrif) then
                     sminn_to_denit_decomp_cascade_vr(c,j,k) = -CNDecompParamsInst%dnp * pmnf_decomp_cascade(c,j,k)
                  end if
               end if
               decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
               decomp_cascade_ctransfer_vr(c,j,k) = (1._r8 - rf_decomp_cascade(c,j,k)) * p_decomp_cpool_loss(c,j,k)
               if (decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. cascade_receiver_pool(k) .ne. i_atm) then
                  decomp_cascade_ntransfer_vr(c,j,k) = p_decomp_cpool_loss(c,j,k) / cn_decomp_pools(c,j,cascade_donor_pool(k))
               else
                  decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
               endif
               if ( cascade_receiver_pool(k) .ne. 0 ) then
                  decomp_cascade_sminn_flux_vr(c,j,k) = pmnf_decomp_cascade(c,j,k)
               else  ! keep sign convention negative for terminal pools
                  decomp_cascade_sminn_flux_vr(c,j,k) = - pmnf_decomp_cascade(c,j,k)
               endif
               net_nmin_vr(c,j) = net_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
            else
               decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
               if (.not. use_nitrif_denitrif) then
                  sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
               end if
               decomp_cascade_sminn_flux_vr(c,j,k) = 0._r8
            end if
            
         end do
      end do
   end do
   
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
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         net_nmin(c) = net_nmin(c) + net_nmin_vr(c,j) * dzsoi_decomp(j)
         gross_nmin(c) = gross_nmin(c) + gross_nmin_vr(c,j) * dzsoi_decomp(j)   
      end do
   end do
   
 end associate
 end subroutine CNDecompAlloc
 
end module CNDecompMod
