module SoilLittDecompMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in litter and soil decomposition model
  ! for coupled carbon-nitrogen code.
  !
  ! !USES:
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_const_mod          , only : SHR_CONST_TKFRZ
  use decompMod              , only : bounds_type
  use perfMod_GPU
  use elm_varctl             , only : iulog, use_lch4, use_century_decomp
  use elm_varcon             , only : dzsoi_decomp
  use elm_varpar             , only : nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
  use DecompCascadeCNMod     , only : decomp_rate_constants_cn
  use DecompCascadeBGCMod    , only : decomp_rate_constants_bgc
  use NitrifDenitrifMod      , only : nitrif_denitrif
  use VerticalProfileMod     , only : decomp_vertprofiles
  use CNDecompCascadeConType , only : decomp_cascade_con
  use CNStateType            , only : cnstate_type
  !!  add phosphorus  -X. YANG
  use elm_varctl             , only : nu_com

  use PhotosynthesisType     , only : photosyns_type
  use CanopyStateType        , only : canopystate_type
  use SoilStateType          , only : soilstate_type
  use CH4Mod                 , only : ch4_type
  use cropType               , only : crop_type
  use ColumnDataType         , only : col_cs, col_cf
  use ColumnDataType         , only : col_ns, col_nf
  use ColumnDataType         , only : col_ps, col_pf
  use VegetationDataType     , only : veg_ps, veg_pf
  ! clm interface & pflotran:
  use elm_varctl             , only : use_elm_interface, use_pflotran, pf_cmode
  use elm_varctl             , only : use_cn, use_fates
  !
  implicit none
  save
  private :: CNvariables_nan4pf  !pflotran
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  public :: readSoilLittDecompParams
  public :: SoilLittDecompAlloc
  ! pflotran
  public :: SoilLittDecompAlloc2
  !
  type, public :: CNDecompParamsType
     real(r8), pointer :: dnp => null()         !denitrification proportion
  end type CNDecompParamsType

  type(CNDecompParamsType)  , public ::  CNDecompParamsInst
  !$acc declare create(CNDecompParamsInst)
  !-----------------------------------------------------------------------

contains
   !-----------------------------------------------------------------------
    subroutine readSoilLittDecompParams ( ncid)
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
      allocate(CNDecompParamsInst%dnp)
      tString='dnp'
      call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
      CNDecompParamsInst%dnp=tempr

    end subroutine readSoilLittDecompParams

!-----------------------------------------------
  subroutine SoilLittDecompAlloc (bounds, num_soilc, filter_soilc,    &
                num_soilp, filter_soilp,                        &
                canopystate_vars, soilstate_vars,               &
                cnstate_vars, ch4_vars, dtime)

    !-----------------------------------------------------------------------------
    ! DESCRIPTION:
    ! Modified for clm_interface: 9/12/2015
    ! clm-bgc soil Module, can be called through clm_bgc_interface
    ! ONLY includes SOM decomposition & nitrification/denitrification
    ! CNAllocaiton is divided into 3 subroutines:
    ! (1) Allocation1_PlantNPDemand  is called in EcosystemDynNoLeaching1
    ! (2) Allocation2_ResolveNPLimit is called in SoilLittDecompAlloc (this subroutine)
    ! (3) Allocation3_PlantCNPAlloc  is called in SoilLittDecompAlloc2
    !-----------------------------------------------------------------------------

    ! !USES:
    use AllocationMod , only: Allocation2_ResolveNPLimit ! Phase-2 of CNAllocation
    ! !ARGUMENT:
    type(bounds_type)        , intent(in)    :: bounds 
    integer                  , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    integer                  , intent(in)    :: num_soilp          ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)    ! filter for soil patches
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    real(r8),   intent(in)    :: dtime
    !
    ! !LOCAL VARIABLES:
    integer :: c,j,k,l,m                                                                 !indices
    integer :: fc                                                                        !lake filter column index
    real(r8):: p_decomp_cpool_loss(num_soilc,1:nlevdecomp,1:ndecomp_cascade_transitions) !potential C loss from one pool to another
    real(r8):: immob(num_soilc,1:nlevdecomp)                                             !potential N immobilization
    real(r8):: immob_p(num_soilc,1:nlevdecomp)                                           !potential P immobilization
    real(r8):: ratio                                                                     !temporary variable
    real(r8):: dnp                                                                       !denitrification proportion
    real(r8):: cn_decomp_pools(num_soilc,1:nlevdecomp,1:ndecomp_pools)
    real(r8):: cp_decomp_pools(num_soilc,1:nlevdecomp,1:ndecomp_pools)
    real(r8):: cp_decomp_pools_new(num_soilc,1:nlevdecomp,1:ndecomp_pools)              !C:P ratio of new SOM
    integer, parameter :: i_atm = 0
    integer :: k_donor_pool
    ! For methane code
    !real(r8):: hrsum(num_soilc,1:nlevdecomp)          !sum of HR (gC/m2/s)
    real :: startt, stopt
    real(r8) :: sum_1, sum_2, sum_3, sum_4

    character(len=64) :: event
    !-----------------------------------------------------------------------

    associate(                                                                                           &
         cascade_donor_pool               =>    decomp_cascade_con%cascade_donor_pool                  , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool            =>    decomp_cascade_con%cascade_receiver_pool               , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step
         floating_cn_ratio_decomp_pools   =>    decomp_cascade_con%floating_cn_ratio_decomp_pools      , & ! Input:  [logical  (:)     ]  TRUE => pool has fixed C:N ratio
         floating_cp_ratio_decomp_pools   =>    decomp_cascade_con%floating_cp_ratio_decomp_pools      , & ! Input:  [logical  (:)     ]  TRUE => pool has fixed C:P ratio
         initial_cn_ratio                 =>    decomp_cascade_con%initial_cn_ratio                    , & ! Input:  [real(r8) (:)     ]  c:n ratio for initialization of pools
         initial_cp_ratio                 =>    decomp_cascade_con%initial_cp_ratio                    , & ! Input:  [real(r8) (:)     ]  c:p ratio for initialization of pools

         is_cwd                           =>    decomp_cascade_con%is_cwd                 , &
         fpi_vr                           =>    cnstate_vars%fpi_vr_col                   , & ! Input:  [real(r8) (:,:)   ]  fraction of potential immobilization for N(no units)
         fpi_p_vr                         =>    cnstate_vars%fpi_p_vr_col                 , & ! Input:  [real(r8) (:,:)   ]  fraction of potential immobilization for P (no units)
         rf_decomp_cascade                =>    cnstate_vars%rf_decomp_cascade_col        , & ! Input:  [real(r8) (:,:,:) ]  respired fraction in decomposition step (frac)
         pathfrac_decomp_cascade          =>    cnstate_vars%pathfrac_decomp_cascade_col  , & ! Input:  [real(r8) (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)

         decomp_npools_vr                 =>    col_ns%decomp_npools_vr                , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
         decomp_ppools_vr                 =>    col_ps%decomp_ppools_vr              , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) P pools

         decomp_cpools_vr                 =>    col_cs%decomp_cpools_vr                  , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools

         w_scalar                         =>    col_cf%w_scalar                           , & ! Input:  [real(r8) (:,:)   ]  fraction by which decomposition is limited by moisture availability

         decomp_cascade_ntransfer_vr      =>    col_nf%decomp_cascade_ntransfer_vr      , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
         decomp_cascade_sminn_flux_vr     =>    col_nf%decomp_cascade_sminn_flux_vr     , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
         potential_immob_vr               =>    col_nf%potential_immob_vr               , & ! Output: [real(r8) (:,:)   ]
         sminn_to_denit_decomp_cascade_vr =>    col_nf%sminn_to_denit_decomp_cascade_vr , & ! Output: [real(r8) (:,:,:) ]
         gross_nmin_vr                    =>    col_nf%gross_nmin_vr                    , & ! Output: [real(r8) (:,:)   ]
         net_nmin_vr                      =>    col_nf%net_nmin_vr                      , & ! Output: [real(r8) (:,:)   ]
         gross_nmin                       =>    col_nf%gross_nmin                       , & ! Output: [real(r8) (:)     ]  gross rate of N mineralization (gN/m2/s)
         net_nmin                         =>    col_nf%net_nmin                         , & ! Output: [real(r8) (:)     ]  net rate of N mineralization (gN/m2/s)
         ! add phosphorus
         decomp_cascade_ptransfer_vr      =>    col_pf%decomp_cascade_ptransfer_vr    , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
         decomp_cascade_sminp_flux_vr     =>    col_pf%decomp_cascade_sminp_flux_vr   , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
         potential_immob_p_vr             =>    col_pf%potential_immob_p_vr           , & ! Output: [real(r8) (:,:)   ]
         gross_pmin_vr                    =>    col_pf%gross_pmin_vr                  , & ! Output: [real(r8) (:,:)   ]
         net_pmin_vr                      =>    col_pf%net_pmin_vr                    , & ! Output: [real(r8) (:,:)   ]
         gross_pmin                       =>    col_pf%gross_pmin                     , & ! Output: [real(r8) (:)     ]  gross rate of P mineralization (gP/m2/s)
         net_pmin                         =>    col_pf%net_pmin                       , & ! Output: [real(r8) (:)     ]  net rate of P mineralization (gP/m2/s)

         decomp_cascade_hr_vr             =>    col_cf%decomp_cascade_hr_vr           , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      =>    col_cf%decomp_cascade_ctransfer_vr    , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_k                         =>    col_cf%decomp_k                       , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         phr_vr                           =>    col_cf%phr_vr                         , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
         fphr                             =>    col_cf%fphr                           , & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic
         pmnf_decomp_cascade              =>    col_nf%pmnf_decomp_cascade            , &
         pmpf_decomp_cascade              =>    col_pf%pmpf_decomp_cascade            , &
         soil_n_immob_flux                =>    col_nf%soil_n_immob_flux              , &
         soil_n_immob_flux_vr             =>    col_nf%soil_n_immob_flux_vr           , &
         soil_n_grossmin_flux             =>    col_nf%soil_n_grossmin_flux           , &
         soil_p_immob_flux                =>    col_pf%soil_p_immob_flux              , &
         soil_p_immob_flux_vr             =>    col_pf%soil_p_immob_flux_vr           , &
         soil_p_grossmin_flux             =>    col_pf%soil_p_grossmin_flux           , &
         actual_immob_vr                  =>    col_nf%actual_immob_vr                , &
         actual_immob_p_vr                =>    col_pf%actual_immob_p_vr                &
         )

      !-------------------------------------------------------------------------------------------------
      ! call decomp_rate_constants_bgc() or decomp_rate_constants_cn(): now called in EcosystemDynNoLeaching1
      !-------------------------------------------------------------------------------------------------

      ! column loop to calculate potential decomp rates and total immobilization
      ! demand.
      call cpu_time(startt) 
      !$acc enter data create(cn_decomp_pools(:num_soilc,1:nlevdecomp,1:ndecomp_pools), &
      !$acc                                     p_decomp_cpool_loss(:num_soilc,1:nlevdecomp,1:ndecomp_cascade_transitions))
      !$acc enter data create(cp_decomp_pools(:num_soilc,1:nlevdecomp,1:ndecomp_pools),&
      !$acc                   immob(:num_soilc,1:nlevdecomp) ,immob_p(:num_soilc,1:nlevdecomp)  ) 
      !$acc enter data create(cp_decomp_pools_new(:num_soilc,1:nlevdecomp,1:ndecomp_pools)) 
      call cpu_time(stopt) 
      write(iulog,*) "SoilLittDecompAlloc1 :: create",(stopt-startt)*1.E+3,"ms" 
      !! calculate c:n ratios of applicable pools

      !$acc parallel loop gang independent collapse(2) default(present)
      do l = 1, ndecomp_pools
         do j = 1,nlevdecomp
            !$acc loop vector independent private(c)
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if ( floating_cn_ratio_decomp_pools(l) ) then
                  if ( decomp_npools_vr(c,j,l) > 0._r8 ) then
                      cn_decomp_pools(fc,j,l) = decomp_cpools_vr(c,j,l) / decomp_npools_vr(c,j,l)
                  end if
               else
                  cn_decomp_pools(fc,j,l) = initial_cn_ratio(l)
               end if
               !! calculate c:p ratios of applicable pools
               if ( floating_cp_ratio_decomp_pools(l) ) then
                   if ( decomp_ppools_vr(c,j,l) > 0._r8 ) then
                      cp_decomp_pools(fc,j,l) = decomp_cpools_vr(c,j,l) / decomp_ppools_vr(c,j,l)
                   end if
                else
                   cp_decomp_pools(fc,j,l) = initial_cp_ratio(l)
                end if
              cp_decomp_pools_new(fc,j,l) = initial_cp_ratio(l)
            end do
         end do
      end do

      ! calculate the non-nitrogen-limited fluxes
      ! these fluxes include the  "/ dt" term to put them on a
      ! per second basis, since the rate constants have been
      ! calculated on a per timestep basis.
       !$acc parallel loop  gang independent collapse(2) default(present)
       do k = 1, ndecomp_cascade_transitions
          do j = 1,nlevdecomp
             !$acc loop vector private(c,k_donor_pool,ratio)
             do fc = 1,num_soilc
                k_donor_pool = cascade_donor_pool(k)
                c = filter_soilc(fc)
                ! set initial values for potential C and N fluxes
                p_decomp_cpool_loss(fc,j,k) = 0._r8
                pmnf_decomp_cascade(c ,j,k) = 0._r8
                pmpf_decomp_cascade(c ,j,k) = 0._r8    !! initial values for potential P fluxes

                if (decomp_cpools_vr(c,j,k_donor_pool) > 0._r8 .and. &
                     decomp_k(c,j,k_donor_pool) > 0._r8 ) then

                   p_decomp_cpool_loss(fc,j,k) = decomp_cpools_vr(c,j,k_donor_pool) &
                            * decomp_k(c,j,k_donor_pool)  * pathfrac_decomp_cascade(c,j,k)
                   if ( .not. floating_cn_ratio_decomp_pools(cascade_receiver_pool(k)) ) then  !! not transition of cwd to litter

                      if (cascade_receiver_pool(k) /= i_atm ) then  ! not 100% respiration
                         ratio = 0._r8
                         if (decomp_npools_vr(c,j,k_donor_pool) > 0._r8) then
                            ratio = cn_decomp_pools(fc,j,cascade_receiver_pool(k))/cn_decomp_pools(fc,j,k_donor_pool)
                         endif
                         pmnf_decomp_cascade(c,j,k) = (p_decomp_cpool_loss(fc,j,k) * (1.0_r8 - rf_decomp_cascade(c,j,k) - ratio) &
                              / cn_decomp_pools(fc,j,cascade_receiver_pool(k)) )

                      else   ! 100% respiration
                         pmnf_decomp_cascade(c,j,k) = - p_decomp_cpool_loss(fc,j,k) / cn_decomp_pools(fc,j,k_donor_pool)
                      endif

                   else   ! CWD -> litter
                      pmnf_decomp_cascade(c,j,k) = 0._r8
                   end if
                   !!! add phosphorus fluxes
                   if ( .not. is_cwd( (cascade_receiver_pool(k)) ) ) then  !! not transition of cwd to litter
                      if (cascade_receiver_pool(k) /= i_atm ) then  ! not 100% respiration
                         ratio = 0._r8
                         if (decomp_ppools_vr(c,j,k_donor_pool) > 0._r8) then
                            ratio = cp_decomp_pools_new(fc,j,cascade_receiver_pool(k))/cp_decomp_pools(fc,j,k_donor_pool)
                         endif
                         pmpf_decomp_cascade(c,j,k) = (p_decomp_cpool_loss(fc,j,k) * (1.0_r8 - rf_decomp_cascade(c,j,k) - ratio) &
                              / cp_decomp_pools_new(fc,j,cascade_receiver_pool(k)) )
                      else   ! 100% respiration
                         pmpf_decomp_cascade(c,j,k) = - p_decomp_cpool_loss(fc,j,k) / cp_decomp_pools(fc,j,k_donor_pool)
                      endif
                   else   ! CWD -> litter
                      pmpf_decomp_cascade(c,j,k) = 0._r8
                   end if
                end if
             end do
          end do
       end do
       
      !$acc enter data create(sum_1,sum_2,sum_3,sum_4)
      !$acc  parallel loop independent gang collapse(2) default(present) private(c,sum_1,sum_2,sum_3,sum_4)
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            sum_1 = 0._r8
            sum_2 = 0._r8
            sum_3 = 0._r8
            sum_4 = 0._r8
            !$acc loop vector reduction(+:sum_1,sum_2,sum_3,sum_4)
            do k = 1, ndecomp_cascade_transitions
               if (pmnf_decomp_cascade(c,j,k) > 0._r8) then
                  sum_1 = sum_1 + pmnf_decomp_cascade(c,j,k)
               else
                  sum_2 = sum_2 - pmnf_decomp_cascade(c,j,k)
               end if
               if (pmpf_decomp_cascade(c,j,k) > 0._r8) then
                  sum_3 = sum_3 + pmpf_decomp_cascade(c,j,k)
               else
                  sum_4 = sum_4 - pmpf_decomp_cascade(c,j,k)
               end if
            end do
            immob(fc,j) = sum_1
            gross_nmin_vr(c,j) = sum_2
            immob_p(fc,j) = sum_3
            gross_pmin_vr(c,j) = sum_4
         end do
      end do

      !$acc parallel loop independent gang worker collapse(2) private(c,sum_1) default(present)
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            sum_1 = 0._r8
            !$acc loop vector reduction(+:sum_1)
            do k = 1, ndecomp_cascade_transitions
               sum_1 = sum_1 + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(fc,j,k)
            end do
            phr_vr(c,j) = sum_1
         end do
      end do
      
      !$acc parallel loop independent gang  vector collapse(2) default(present) 
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            potential_immob_vr(c,j) = immob(fc,j)
            potential_immob_p_vr(c,j) = immob_p(fc,j)
         end do
      end do

      !-------------------------------------------------------------------------------------------------
      ! 'call decomp_vertprofiles()' (calc nfixation_prof) is moved to EcosystemDynNoLeaching1
      ! 'nfixation_prof' is used in 'calc_nuptake_prof' & 'calc_puptake_prof', which are called in Allocation1,2,3
      !-------------------------------------------------------------------------------------------------
      call nitrif_denitrif( num_soilc,filter_soilc, &
                     soilstate_vars,  ch4_vars )

      ! now that potential N immobilization is known, call allocation
      ! to resolve the competition between plants and soil heterotrophs
      ! for available soil mineral N resource.
      ! in addition, calculate fpi_vr, fpi_p_vr, & fgp
      call Allocation2_ResolveNPLimit(bounds,num_soilc, filter_soilc, &
               cnstate_vars, soilstate_vars, dtime)
      ! column loop to calculate actual immobilization and decomp rates, following
      ! resolution of plant/heterotroph  competition for mineral N
      !-------------------------------------------------------------------------------------------------
      ! delete c:n,c:p ratios calculation, they have been calculated at the beginning of this subroutine
      !-------------------------------------------------------------------------------------------------
      ! upon return from CNAllocation, the fraction of potential immobilization
      ! has been set (cnstate_vars%fpi_vr_col). now finish the decomp calculations.
      ! Only the immobilization steps are limited by fpi_vr (pmnf > 0)
      ! Also calculate denitrification losses as a simple proportion
      ! of mineralization flux.

      !$acc parallel loop independent gang worker collapse(2) default(present)
      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            !$acc loop vector private(c,k_donor_pool)
            do fc = 1,num_soilc
               k_donor_pool = cascade_donor_pool(k)
               c = filter_soilc(fc)

               if (decomp_cpools_vr(c,j,k_donor_pool) > 0._r8) then
                  if ( pmnf_decomp_cascade(c,j,k) > 0._r8 .and. pmpf_decomp_cascade(c,j,k) > 0._r8 ) then    ! N and P co-limitation
                     p_decomp_cpool_loss(fc,j,k) = p_decomp_cpool_loss(fc,j,k) * min( fpi_vr(c,j),fpi_p_vr(c,j) )
                     pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * min( fpi_vr(c,j),fpi_p_vr(c,j) )
                     pmpf_decomp_cascade(c,j,k) = pmpf_decomp_cascade(c,j,k) * min( fpi_vr(c,j),fpi_p_vr(c,j) )   !!! immobilization step
                  elseif ( pmnf_decomp_cascade(c,j,k) > 0._r8 .and. pmpf_decomp_cascade(c,j,k) <= 0._r8 ) then  ! N limitation
                     p_decomp_cpool_loss(fc,j,k) = p_decomp_cpool_loss(fc,j,k) * fpi_vr(c,j)
                     pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * fpi_vr(c,j)
                     pmpf_decomp_cascade(c,j,k) = pmpf_decomp_cascade(c,j,k) * fpi_vr(c,j) !!! immobilization step
                  elseif ( pmnf_decomp_cascade(c,j,k) <= 0._r8 .and. pmpf_decomp_cascade(c,j,k) >  0._r8 ) then  ! P limitation
                     p_decomp_cpool_loss(fc,j,k) = p_decomp_cpool_loss(fc,j,k) * fpi_p_vr(c,j)
                     pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * fpi_p_vr(c,j)
                     pmpf_decomp_cascade(c,j,k) = pmpf_decomp_cascade(c,j,k) * fpi_p_vr(c,j) !!! immobilization step
                  end if
                  decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(fc,j,k)
                  decomp_cascade_ctransfer_vr(c,j,k) = (1._r8 - rf_decomp_cascade(c,j,k)) * p_decomp_cpool_loss(fc,j,k)
                  if (decomp_npools_vr(c,j,k_donor_pool) > 0._r8 .and. cascade_receiver_pool(k) /= i_atm) then
                     decomp_cascade_ntransfer_vr(c,j,k) = p_decomp_cpool_loss(fc,j,k) / cn_decomp_pools(fc,j,k_donor_pool)
                  else
                     decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
                  endif
                  !!! phosphorus fluxes
                  if (decomp_ppools_vr(c,j,k_donor_pool) > 0._r8 .and. cascade_receiver_pool(k) /= i_atm) then
                     decomp_cascade_ptransfer_vr(c,j,k) = p_decomp_cpool_loss(fc,j,k) / cp_decomp_pools(fc,j,k_donor_pool)
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

               else
                  decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
                  decomp_cascade_ptransfer_vr(c,j,k) = 0._r8
                  decomp_cascade_sminn_flux_vr(c,j,k) = 0._r8
                  decomp_cascade_sminp_flux_vr(c,j,k) = 0._r8
               end if

            end do
         end do
      end do


      !$acc parallel loop independent gang worker collapse(2) private(c,sum_1,sum_2) default(present)
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            sum_1 = 0._r8
            sum_2 = 0._r8
            !$acc loop vector reduction(+:sum_1,sum_2)
            do k = 1, ndecomp_cascade_transitions
               sum_1 = sum_1 - pmnf_decomp_cascade(c,j,k)
               sum_2 = sum_2 - pmpf_decomp_cascade(c,j,k)
            end do
            net_nmin_vr(c,j) = sum_1
            net_pmin_vr(c,j) = sum_2
         end do
      end do

      if (nu_com .eq. 'RD') then
         !$acc parallel loop independent gang worker collapse(2) private(c,sum_1,sum_2)
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               sum_1 = 0.0_r8
               sum_2 = 0.0_r8
               !$acc loop vector reduction(+:sum_1,sum_2)
               do k = 1, ndecomp_cascade_transitions
                  if (pmnf_decomp_cascade(c,j,k) <= 0._r8) then
                      sum_1 = sum_1 - 1.0_r8*pmnf_decomp_cascade(c,j,k)
                  end if
                  if (pmpf_decomp_cascade(c,j,k) <= 0._r8) then
                      sum_2 = sum_2 - 1.0_r8*pmpf_decomp_cascade(c,j,k)
                  end if
               end do
               gross_nmin_vr(c,j) = sum_1
               gross_pmin_vr(c,j) = sum_2
             end do
          end do
      end if

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
         !$acc parallel loop independent gang worker collapse(2) private(c,sum_1) default(present)
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               ! Calculate total fraction of potential HR, for methane code
               c = filter_soilc(fc)
               sum_1 = 0.0_r8
               if(phr_vr(c,j) > 0._r8) then
                  !$acc loop vector reduction(+:sum_1)
                  do k = 1, ndecomp_cascade_transitions
                     sum_1 = sum_1 + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(fc,j,k)
                  end do
                  ! Nitrogen limitation / (low)-moisture limitation
                  fphr(c,j) = sum_1/phr_vr(c,j)* w_scalar(c,j)
                  fphr(c,j) = max(fphr(c,j),0.01_r8) ! Prevent overflow errors for 0 respiration
               else
                  fphr(c,j) = 1._r8
               end if
            end do
         end do

      end if

      !$acc exit data delete(sum_1,sum_2,sum_3,sum_4)
      !$acc exit data delete(cn_decomp_pools(:,:,:), p_decomp_cpool_loss(:,:,:), &
      !$acc   cp_decomp_pools(:,:,:),immob(:,:) ,immob_p(:,:) , &
      !$acc   cp_decomp_pools_new(:,:,:) )

    end associate

  end subroutine SoilLittDecompAlloc

!-------------------------------------------------------------------------------------------------

  subroutine SoilLittDecompAlloc2 ( num_soilc, filter_soilc, num_soilp, filter_soilp,   &
        canopystate_vars, soilstate_vars,          &
       cnstate_vars, crop_vars,  dt)
    !-----------------------------------------------------------------------------
    ! DESCRIPTION:
    ! bgc interface & pflotran:
    ! (1) Simplified codes of SoilLittDecompAlloc subroutine for coupling with pflotran
    ! (2) call Allocation3_PlantCNPAlloc
    ! (3) calculate net_nmin(c), gross_nmin(c), net_pmin(c), gross_pmin(c)
    !-----------------------------------------------------------------------------

    ! !USES:
    use AllocationMod , only: Allocation3_PlantCNPAlloc ! Phase-3 of CNAllocation
    !
    ! !ARGUMENT:
    integer                  , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    integer                  , intent(in)    :: num_soilp          ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)    ! filter for soil patches
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(crop_type)          , intent(inout) :: crop_vars
    real(r8), intent(in) :: dt                                           ! time step (seconds)

    !
    ! !LOCAL VARIABLES:
    integer :: fc, c, j                                     ! indices
    !real(r8):: col_plant_ndemand(bounds%begc:bounds%endc)   ! column-level vertically-integrated plant N demand (gN/m2/s)
    real(r8) :: smin_nh4_to_plant_vr_loc(1:num_soilc,1:nlevdecomp)
    real(r8) :: smin_no3_to_plant_vr_loc(1:num_soilc,1:nlevdecomp)

    ! For methane code
    real(r8):: hrsum(1:num_soilc,1:nlevdecomp) !sum of HR (gC/m2/s)
    real(r8):: sum1,sum2,sum3,sum4

    character(len=64) :: event
    !-----------------------------------------------------------------------

    associate(                                                                                      &
         gross_nmin_vr                    =>    col_nf%gross_nmin_vr                    , & ! Output: [real(r8) (:,:)   ]
         net_nmin_vr                      =>    col_nf%net_nmin_vr                      , & ! Output: [real(r8) (:,:)   ]
         gross_nmin                       =>    col_nf%gross_nmin                       , & ! Output: [real(r8) (:)     ]  gross rate of N mineralization (gN/m2/s)
         net_nmin                         =>    col_nf%net_nmin                         , & ! Output: [real(r8) (:)     ]  net rate of N mineralization (gN/m2/s)
         !phosphorus
         gross_pmin_vr                    =>    col_pf%gross_pmin_vr                  , & ! Output: [real(r8) (:,:)   ]
         net_pmin_vr                      =>    col_pf%net_pmin_vr                    , & ! Output: [real(r8) (:,:)   ]
         gross_pmin                       =>    col_pf%gross_pmin                     , & ! Output: [real(r8) (:)     ]  gross rate of P mineralization (gP/m2/s)
         net_pmin                         =>    col_pf%net_pmin                       , & ! Output: [real(r8) (:)     ]  net rate of P mineralization (gP/m2/s)

         fpi_vr                           =>    cnstate_vars%fpi_vr_col                                , & ! Output:  [real(r8) (:,:)   ]  fraction of potential immobilization (no units)
         fpi                              =>    cnstate_vars%fpi_col                                   , & ! Output: [real(r8) (:)   ]  fraction of potential immobilization (no units)
         potential_immob_vr               =>    col_nf%potential_immob_vr               , & ! Input:
         actual_immob_vr                  =>    col_nf%actual_immob_vr                  , & ! Input:
         potential_immob                  =>    col_nf%potential_immob                  , & ! Output: [real(r8) (:)   ]
         actual_immob                     =>    col_nf%actual_immob                     , & ! Output: [real(r8) (:)   ]

         fpg                              =>    cnstate_vars%fpg_col                                   , & ! Output: [real(r8) (:)   ]  fraction of potential gpp (no units)
         sminn_to_plant                   =>    col_nf%sminn_to_plant                   , & ! Output: [real(r8) (:)     ]  col N uptake (gN/m2/s)
         sminn_to_plant_vr                =>    col_nf%sminn_to_plant_vr                , & ! Input:  [real(r8) (:,:)    ]  vertically-resolved N uptake (gN/m3/s)

         smin_no3_to_plant_vr             =>    col_nf%smin_no3_to_plant_vr             , & ! Output: [real(r8) (:,:) ]
         smin_nh4_to_plant_vr             =>    col_nf%smin_nh4_to_plant_vr             , & ! Output: [real(r8) (:,:) ]

         col_plant_ndemand_vr             =>    col_nf%plant_ndemand_vr                 , & ! Input:  [real(r8) (:)     ]  col N uptake (gN/m2/s)

         plant_ndemand_col                =>    col_nf%plant_ndemand                    , & ! Output:  [real(r8) (:,:) ]

         w_scalar                         =>    col_cf%w_scalar                           , & ! Input:  [real(r8) (:,:)   ]  fraction by which decomposition is limited by moisture availability
         decomp_cascade_hr_vr             =>    col_cf%decomp_cascade_hr_vr               , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      =>    col_cf%decomp_cascade_ctransfer_vr        , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_k                         =>    col_cf%decomp_k                           , & ! Output: [real(r8) (:,:,:) ]  rate constant for decomposition (1./sec)
         hr_vr                            =>    col_cf%hr_vr                              , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
         phr_vr                           =>    col_cf%phr_vr                             , & ! Output: [real(r8) (:,:)   ]  potential HR (gC/m3/s)
         fphr                             =>    col_cf%fphr                               , & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic

         smin_no3_vr                      =>    col_ns%smin_no3_vr                     , &
         smin_nh4_vr                      =>    col_ns%smin_nh4_vr                       &
         )

      ! set time steps
      !------------------------------------------------------------------
      ! 'call decomp_vertprofiles()' moved to EcosystemDynNoLeaching1
      !------------------------------------------------------------------

      ! MUST have already updated needed bgc variables from PFLOTRAN by this point
      if(use_elm_interface.and.use_pflotran.and.pf_cmode) then
         smin_nh4_to_plant_vr_loc(1:num_soilc,:) = 0._r8
         smin_no3_to_plant_vr_loc(1:num_soilc,:) = 0._r8

         ! fpg calculation
         do fc=1,num_soilc
            c = filter_soilc(fc)
            sminn_to_plant(c)       = 0._r8
            do j = 1, nlevdecomp           ! sum up actual and potential column-level N fluxes to plant
               sminn_to_plant(c)    = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
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
            potential_immob(c) = 0._r8
            actual_immob(c)    = 0._r8
            do j = 1, nlevdecomp
               if (potential_immob_vr(c,j) > 0.0_r8) then
                  fpi_vr(c,j) = actual_immob_vr(c,j) / potential_immob_vr(c,j)
                  potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j)*dzsoi_decomp(j)
                  actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j)*dzsoi_decomp(j)
               else
                  fpi_vr(c,j) = 0.0_r8
               end if
            end do
         end do
         do fc=1,num_soilc
            c = filter_soilc(fc)
            if (potential_immob(c) > 0.0_r8) then
               fpi(c) = max(0._r8,actual_immob(c)) / potential_immob(c)
               fpi(c) = min(1._r8, fpi(c))
            else
               fpi(c) = 1.0_r8
            end if
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
                  hrsum(fc,j) = hr_vr(c,j)
               end do
            end do

            ! Nitrogen limitation / (low)-moisture limitation
            do j = 1,nlevdecomp
               do fc = 1,num_soilc
                  c = filter_soilc(fc)
                  if (phr_vr(c,j) > 0._r8) then
                     fphr(c,j) = hrsum(fc,j) / phr_vr(c,j) * w_scalar(c,j)
                     fphr(c,j) = max(fphr(c,j), 0.01_r8) ! Prevent overflow errors for 0 respiration
                  else
                     fphr(c,j) = 1._r8
                  end if
               end do
            end do
         end if

         ! ! needs to zero CLM-CNP variables NOT available from pflotran bgc coupling
         ! call CNvariables_nan4pf(num_soilc, filter_soilc, &
         !                num_soilp, filter_soilp)

         ! save variables before updating
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            do j = 1,nlevdecomp
                smin_no3_to_plant_vr_loc(fc,j) = smin_no3_to_plant_vr(c,j)
                smin_nh4_to_plant_vr_loc(fc,j) = smin_nh4_to_plant_vr(c,j)
            end do
         end do

      end if !if(use_elm_interface.and.use_pflotran.and.pf_cmode)
      !------------------------------------------------------------------
      ! phase-3 Allocation for plants
      if(.not.use_fates)then
        ! event = 'CNAllocation - phase-3'
        !call t_start_lnd(event)
        call Allocation3_PlantCNPAlloc( &
                  num_soilc, filter_soilc, num_soilp, filter_soilp    , &
                  canopystate_vars                                    , &
                  cnstate_vars, crop_vars, dt)
        !call t_stop_lnd(event)
      end if
      !------------------------------------------------------------------
    if(use_pflotran.and.pf_cmode) then
      ! in Allocation3_PlantCNPAlloc():
      ! smin_nh4_to_plant_vr(c,j), smin_no3_to_plant_vr(c,j), sminn_to_plant_vr(c,j) may be adjusted
      ! therefore, we need to update smin_no3_vr(c,j) & smin_nh4_vr(c,j)
      do fc = 1,num_soilc
           c = filter_soilc(fc)
           do j = 1,nlevdecomp
               smin_no3_vr(c,j) = smin_no3_vr(c,j) - (smin_no3_to_plant_vr(c,j) - smin_no3_to_plant_vr_loc(fc,j))*dt
               smin_nh4_vr(c,j) = smin_nh4_vr(c,j) - (smin_nh4_to_plant_vr(c,j) - smin_nh4_to_plant_vr_loc(fc,j))*dt
               smin_no3_vr(c,j) = max(0._r8, smin_no3_vr(c,j))
               smin_nh4_vr(c,j) = max(0._r8, smin_nh4_vr(c,j))
            end do
      end do
    end if !(use_pflotran.and.pf_cmode)
    !------------------------------------------------------------------
      ! vertically integrate net and gross mineralization fluxes for diagnostic output
      !$acc enter data create(sum1,sum2,sum3,sum4)

      !$acc parallel loop independent gang worker private(c,sum1,sum2,sum3,sum4)
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         sum1 = 0.0_r8; sum2 =0.0_r8;
         sum3 = 0.0_r8; sum4 = 0.0_r8;
         !$acc loop vector reduction(+:sum1,sum2,sum3,sum4)
         do j = 1,nlevdecomp
            sum1 = sum1 + net_nmin_vr(c,j) * dzsoi_decomp(j)
            sum2 = sum2 + gross_nmin_vr(c,j) * dzsoi_decomp(j)
            !phosphorus
            sum3 = sum3 + net_pmin_vr(c,j) * dzsoi_decomp(j)
            sum4 = sum4 + gross_pmin_vr(c,j) * dzsoi_decomp(j)
         end do
         net_nmin(c)   = sum1
         gross_nmin(c) = sum2
         net_pmin(c)   = sum3
         gross_pmin(c) = sum4
      end do
      !$acc exit data delete(sum1,sum2,sum3,sum4) 
    end associate

  end subroutine SoilLittDecompAlloc2

  !-------------------------------------------------------------------------------------------------
  !
  subroutine CNvariables_nan4pf (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
  !
  !DESCRIPTION:
  !  CN variables not available from PFLOTRAN, some of which may be output and may cause issues,
  !  if not properly set.
  !
  !USES:
    use elm_varctl   , only: carbon_only, carbonnitrogen_only
    use elm_varpar   , only: nlevdecomp, ndecomp_cascade_transitions
    use ColumnDataType     , only : col_ps_setvalues, col_pf_setvalues
    use VegetationDataType , only : veg_ps_setvalues, veg_pf_setvalues
   !
   !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc          ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)    ! filter for soil columns
    integer                  , intent(in)    :: num_soilp          ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)    ! filter for soil patches
   !
   !CALLED FROM:
   !
   !LOCAL VARIABLES:
   integer :: c,j,k,fc          !indices
   !
   !-----------------------------------------------------------------------
   associate (&
         decomp_cascade_hr_vr             =>    col_cf%decomp_cascade_hr_vr               , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      =>    col_cf%decomp_cascade_ctransfer_vr        , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ntransfer_vr      =>    col_nf%decomp_cascade_ntransfer_vr        & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
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

   ! pflotran not yet support phosphous cycle
   if ( carbon_only .or.  carbonnitrogen_only  ) then
      call veg_ps_SetValues(veg_ps,num_patch=num_soilp,  filter_patch=filter_soilp,  value_patch=0._r8)
      call col_ps_SetValues(col_ps,num_column=num_soilc, filter_column=filter_soilc, value_column=0._r8)

      call veg_pf_setvalues(veg_pf,num_patch=num_soilp,  filter_patch=filter_soilp,  value_patch=0._r8)
      call col_pf_setvalues(col_pf,num_column=num_soilc, filter_column=filter_soilc, value_column=0._r8)
   end if


  end associate
  end subroutine CNvariables_nan4pf

end module SoilLittDecompMod
