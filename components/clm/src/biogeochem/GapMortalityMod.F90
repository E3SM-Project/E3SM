module GapMortalityMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in gap mortality for coupled carbon
  ! nitrogen code.
  ! add phosphorus fluxes - X.YANG
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use abortutils          , only : endrun
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use VegetationPropertiesType      , only : veg_vp
  use CNDVType            , only : dgvs_type
  use CNStateType         , only : cnstate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use ColumnType          , only : col_pp                
  use VegetationType           , only : veg_pp                

  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type

  use clm_varctl          , only : nu_com

  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: GapMortality
  public :: readGapMortParams

  type, private :: CNGapMortParamsType
      real(r8):: am     ! mortality rate based on annual rate, fractional mortality (1/yr)
      real(r8):: k_mort ! coeff. of growth efficiency in mortality equation
  end type CNGapMortParamsType

  type(CNGapMortParamsType),private ::  CNGapMortParamsInst
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readGapMortParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read in parameters
    !
    ! !USES:
    use ncdio_pio  , only : file_desc_t,ncd_io
    !
    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNGapMortParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    tString='r_mort'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNGapMortParamsInst%am=tempr

    tString='k_mort'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNGapMortParamsInst%k_mort=tempr   
    
  end subroutine readGapMortParams

  !-----------------------------------------------------------------------
  subroutine GapMortality (&
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       dgvs_vars, cnstate_vars, &
       carbonstate_vars, nitrogenstate_vars, carbonflux_vars,nitrogenflux_vars,&
       phosphorusstate_vars,phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! Gap-phase mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use clm_time_manager , only: get_days_per_year
    use clm_varcon       , only: secspday
    use pftvarcon        , only: npcropmin
    use clm_varctl       , only: use_cndv, spinup_state, spinup_mortality_factor
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(dgvs_type)          , intent(inout) :: dgvs_vars
    type(cnstate_type)       , intent(inout)    :: cnstate_vars
    type(carbonstate_type)   , intent(inout)    :: carbonstate_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars

    type(phosphorusstate_type) , intent(in) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars

    !
    ! !LOCAL VARIABLES:
    integer :: p             ! patch index
    integer :: fp            ! patch filter index
    real(r8):: am            ! rate for fractional mortality (1/yr)
    real(r8):: m             ! rate for fractional mortality (1/s)
    real(r8):: mort_max      ! asymptotic max mortality rate (/yr)
    real(r8):: k_mort = 0.3  ! coeff of growth efficiency in mortality equation
    !-----------------------------------------------------------------------

    associate(                                                                                              & 
         ivt                                 =>    veg_pp%itype                                              , & ! Input:  [integer  (:) ]  pft vegetation type                                

         woody                               =>    veg_vp%woody                                       , & ! Input:  [real(r8) (:) ]  binary flag for woody lifeform                    
         
         greffic                             =>    dgvs_vars%greffic_patch                                , & ! Input:  [real(r8) (:) ]                                                    
         heatstress                          =>    dgvs_vars%heatstress_patch                             , & ! Input:  [real(r8) (:) ]                                                    
         nind                                =>    dgvs_vars%nind_patch                                     & ! Output: [real(r8) (:) ]  number of individuals (#/m2) added by F. Li and S. Levis
         )

      ! set the mortality rate based on annual rate
      am = CNGapMortParamsInst%am
      ! set coeff of growth efficiency in mortality equation 
      k_mort = CNGapMortParamsInst%k_mort

      if (nu_com .eq. 'RD') then
          call mortality_rate_soilorder(num_soilp,filter_soilp,cnstate_vars)
      end if


      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         if (use_cndv) then
            ! Stress mortality from lpj's subr Mortality.

            if (woody(ivt(p)) == 1._r8) then

               if (ivt(p) == 8) then
                  mort_max = 0.03_r8 ! BDT boreal
               else
                  mort_max = 0.01_r8 ! original value for all patches
               end if

               ! heatstress and greffic calculated in Establishment once/yr

               ! Mortality rate inversely related to growth efficiency
               ! (Prentice et al 1993)
               am = mort_max / (1._r8 + k_mort * greffic(p))

               ! Mortality rate inversely related to growth efficiency
               ! (Prentice et al 1993)
               am = mort_max / (1._r8 + k_mort * greffic(p))

               am = min(1._r8, am + heatstress(p))
            else ! lpj didn't set this for grasses; cn does
               ! set the mortality rate based on annual rate
               am = CNGapMortParamsInst%am
            end if

         end if

         if (nu_com .eq. 'RD') then
             am = cnstate_vars%r_mort_cal_patch(p)
         end if


         m  = am/(get_days_per_year() * secspday)

         !------------------------------------------------------
         ! patch-level gap mortality carbon fluxes
         !------------------------------------------------------

         ! displayed pools
         carbonflux_vars%m_leafc_to_litter_patch(p)               = carbonstate_vars%leafc_patch(p)               * m
         carbonflux_vars%m_frootc_to_litter_patch(p)              = carbonstate_vars%frootc_patch(p)              * m
         carbonflux_vars%m_livestemc_to_litter_patch(p)           = carbonstate_vars%livestemc_patch(p)           * m
         carbonflux_vars%m_deadstemc_to_litter_patch(p)           = carbonstate_vars%deadstemc_patch(p)           * m
         carbonflux_vars%m_livecrootc_to_litter_patch(p)          = carbonstate_vars%livecrootc_patch(p)          * m
         carbonflux_vars%m_deadcrootc_to_litter_patch(p)          = carbonstate_vars%deadcrootc_patch(p)          * m
         if (spinup_state >= 1) then 
           carbonflux_vars%m_deadstemc_to_litter_patch(p)         = carbonstate_vars%deadstemc_patch(p)*m * spinup_mortality_factor
           carbonflux_vars%m_deadcrootc_to_litter_patch(p)        = carbonstate_vars%deadcrootc_patch(p)*m * spinup_mortality_factor
         end if

         ! storage pools
         carbonflux_vars%m_leafc_storage_to_litter_patch(p)       = carbonstate_vars%leafc_storage_patch(p)       * m
         carbonflux_vars%m_frootc_storage_to_litter_patch(p)      = carbonstate_vars%frootc_storage_patch(p)      * m
         carbonflux_vars%m_livestemc_storage_to_litter_patch(p)   = carbonstate_vars%livestemc_storage_patch(p)   * m
         carbonflux_vars%m_deadstemc_storage_to_litter_patch(p)   = carbonstate_vars%deadstemc_storage_patch(p)   * m
         carbonflux_vars%m_livecrootc_storage_to_litter_patch(p)  = carbonstate_vars%livecrootc_storage_patch(p)  * m
         carbonflux_vars%m_deadcrootc_storage_to_litter_patch(p)  = carbonstate_vars%deadcrootc_storage_patch(p)  * m
         carbonflux_vars%m_gresp_storage_to_litter_patch(p)       = carbonstate_vars%gresp_storage_patch(p)       * m
         carbonflux_vars%m_cpool_to_litter_patch(p)               = carbonstate_vars%cpool_patch(p)               * m
 
         ! transfer pools
         carbonflux_vars%m_leafc_xfer_to_litter_patch(p)          = carbonstate_vars%leafc_xfer_patch(p)          * m
         carbonflux_vars%m_frootc_xfer_to_litter_patch(p)         = carbonstate_vars%frootc_xfer_patch(p)         * m
         carbonflux_vars%m_livestemc_xfer_to_litter_patch(p)      = carbonstate_vars%livestemc_xfer_patch(p)      * m
         carbonflux_vars%m_deadstemc_xfer_to_litter_patch(p)      = carbonstate_vars%deadstemc_xfer_patch(p)      * m
         carbonflux_vars%m_livecrootc_xfer_to_litter_patch(p)     = carbonstate_vars%livecrootc_xfer_patch(p)     * m
         carbonflux_vars%m_deadcrootc_xfer_to_litter_patch(p)     = carbonstate_vars%deadcrootc_xfer_patch(p)     * m
         carbonflux_vars%m_gresp_xfer_to_litter_patch(p)          = carbonstate_vars%gresp_xfer_patch(p)          * m

         !------------------------------------------------------
         ! patch-level gap mortality nitrogen fluxes
         !------------------------------------------------------

         ! displayed pools
         nitrogenflux_vars%m_leafn_to_litter_patch(p)               = nitrogenstate_vars%leafn_patch(p)               * m
         nitrogenflux_vars%m_frootn_to_litter_patch(p)              = nitrogenstate_vars%frootn_patch(p)              * m
         nitrogenflux_vars%m_livestemn_to_litter_patch(p)           = nitrogenstate_vars%livestemn_patch(p)           * m
         nitrogenflux_vars%m_deadstemn_to_litter_patch(p)           = nitrogenstate_vars%deadstemn_patch(p)           * m
         nitrogenflux_vars%m_livecrootn_to_litter_patch(p)          = nitrogenstate_vars%livecrootn_patch(p)          * m
         nitrogenflux_vars%m_deadcrootn_to_litter_patch(p)          = nitrogenstate_vars%deadcrootn_patch(p)          * m
         if (ivt(p) < npcropmin) then
            nitrogenflux_vars%m_retransn_to_litter_patch(p) = nitrogenstate_vars%retransn_patch(p) * m
         end if
         nitrogenflux_vars%m_npool_to_litter_patch(p)               = nitrogenstate_vars%npool_patch(p)               * m

         if (spinup_state >= 1) then
           nitrogenflux_vars%m_deadstemn_to_litter_patch(p)         = nitrogenstate_vars%deadstemn_patch(p)  * m &
                * spinup_mortality_factor
           nitrogenflux_vars%m_deadcrootn_to_litter_patch(p)        = nitrogenstate_vars%deadcrootn_patch(p) * m &
                * spinup_mortality_factor
         end if
   
         ! storage pools
         nitrogenflux_vars%m_leafn_storage_to_litter_patch(p)       = nitrogenstate_vars%leafn_storage_patch(p)       * m
         nitrogenflux_vars%m_frootn_storage_to_litter_patch(p)      = nitrogenstate_vars%frootn_storage_patch(p)      * m
         nitrogenflux_vars%m_livestemn_storage_to_litter_patch(p)   = nitrogenstate_vars%livestemn_storage_patch(p)   * m
         nitrogenflux_vars%m_deadstemn_storage_to_litter_patch(p)   = nitrogenstate_vars%deadstemn_storage_patch(p)   * m
         nitrogenflux_vars%m_livecrootn_storage_to_litter_patch(p)  = nitrogenstate_vars%livecrootn_storage_patch(p)  * m
         nitrogenflux_vars%m_deadcrootn_storage_to_litter_patch(p)  = nitrogenstate_vars%deadcrootn_storage_patch(p)  * m

         ! transfer pools
         nitrogenflux_vars%m_leafn_xfer_to_litter_patch(p)          = nitrogenstate_vars%leafn_xfer_patch(p)          * m
         nitrogenflux_vars%m_frootn_xfer_to_litter_patch(p)         = nitrogenstate_vars%frootn_xfer_patch(p)         * m
         nitrogenflux_vars%m_livestemn_xfer_to_litter_patch(p)      = nitrogenstate_vars%livestemn_xfer_patch(p)      * m
         nitrogenflux_vars%m_deadstemn_xfer_to_litter_patch(p)      = nitrogenstate_vars%deadstemn_xfer_patch(p)      * m
         nitrogenflux_vars%m_livecrootn_xfer_to_litter_patch(p)     = nitrogenstate_vars%livecrootn_xfer_patch(p)     * m
         nitrogenflux_vars%m_deadcrootn_xfer_to_litter_patch(p)     = nitrogenstate_vars%deadcrootn_xfer_patch(p)     * m

         !------------------------------------------------------
         ! patch-level gap mortality phosphorus fluxes
         !------------------------------------------------------

         ! displayed pools
         phosphorusflux_vars%m_leafp_to_litter_patch(p)               = phosphorusstate_vars%leafp_patch(p)               * m
         phosphorusflux_vars%m_frootp_to_litter_patch(p)              = phosphorusstate_vars%frootp_patch(p)              * m
         phosphorusflux_vars%m_livestemp_to_litter_patch(p)           = phosphorusstate_vars%livestemp_patch(p)           * m
         phosphorusflux_vars%m_deadstemp_to_litter_patch(p)           = phosphorusstate_vars%deadstemp_patch(p)           * m
         phosphorusflux_vars%m_livecrootp_to_litter_patch(p)          = phosphorusstate_vars%livecrootp_patch(p)          * m
         phosphorusflux_vars%m_deadcrootp_to_litter_patch(p)          = phosphorusstate_vars%deadcrootp_patch(p)          * m
         if (ivt(p) < npcropmin) then
            phosphorusflux_vars%m_retransp_to_litter_patch(p) = phosphorusstate_vars%retransp_patch(p) * m
         end if
         phosphorusflux_vars%m_ppool_to_litter_patch(p)               = phosphorusstate_vars%ppool_patch(p)               * m

         if (spinup_state >= 1) then
           phosphorusflux_vars%m_deadstemp_to_litter_patch(p)         = phosphorusstate_vars%deadstemp_patch(p)  * m &
                * spinup_mortality_factor
           phosphorusflux_vars%m_deadcrootp_to_litter_patch(p)        = phosphorusstate_vars%deadcrootp_patch(p) * m &
                * spinup_mortality_factor
         end if
           
         ! storage pools
         phosphorusflux_vars%m_leafp_storage_to_litter_patch(p)       = phosphorusstate_vars%leafp_storage_patch(p)       * m
         phosphorusflux_vars%m_frootp_storage_to_litter_patch(p)      = phosphorusstate_vars%frootp_storage_patch(p)      * m
         phosphorusflux_vars%m_livestemp_storage_to_litter_patch(p)   = phosphorusstate_vars%livestemp_storage_patch(p)   * m
         phosphorusflux_vars%m_deadstemp_storage_to_litter_patch(p)   = phosphorusstate_vars%deadstemp_storage_patch(p)   * m
         phosphorusflux_vars%m_livecrootp_storage_to_litter_patch(p)  = phosphorusstate_vars%livecrootp_storage_patch(p)  * m
         phosphorusflux_vars%m_deadcrootp_storage_to_litter_patch(p)  = phosphorusstate_vars%deadcrootp_storage_patch(p)  * m

         ! transfer pools
         phosphorusflux_vars%m_leafp_xfer_to_litter_patch(p)          = phosphorusstate_vars%leafp_xfer_patch(p)          * m
         phosphorusflux_vars%m_frootp_xfer_to_litter_patch(p)         = phosphorusstate_vars%frootp_xfer_patch(p)         * m
         phosphorusflux_vars%m_livestemp_xfer_to_litter_patch(p)      = phosphorusstate_vars%livestemp_xfer_patch(p)      * m
         phosphorusflux_vars%m_deadstemp_xfer_to_litter_patch(p)      = phosphorusstate_vars%deadstemp_xfer_patch(p)      * m
         phosphorusflux_vars%m_livecrootp_xfer_to_litter_patch(p)     = phosphorusstate_vars%livecrootp_xfer_patch(p)     * m
         phosphorusflux_vars%m_deadcrootp_xfer_to_litter_patch(p)     = phosphorusstate_vars%deadcrootp_xfer_patch(p)     * m

         ! added by F. Li and S. Levis
         if (use_cndv) then
            if (woody(ivt(p)) == 1._r8)then
               if (carbonstate_vars%livestemc_patch(p) + carbonstate_vars%deadstemc_patch(p)> 0._r8)then
                  nind(p)=nind(p)*(1._r8-m)
               else
                  nind(p) = 0._r8 
               end if
            end if
         end if

      end do ! end of pft loop

      ! gather all pft-level litterfall fluxes to the column
      ! for litter C and N inputs

      call CNGapPftToColumn(num_soilc, filter_soilc, &
           cnstate_vars, carbonflux_vars, nitrogenflux_vars,phosphorusflux_vars)

    end associate
  end subroutine GapMortality

  !-----------------------------------------------------------------------
  subroutine CNGapPftToColumn ( &
       num_soilc, filter_soilc, &
       cnstate_vars, carbonflux_vars, nitrogenflux_vars,phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! called in the middle of CNGapMoratlity to gather all pft-level gap mortality fluxes
    ! to the column level and assign them to the three litter pools
    !
    ! !USES:
    use clm_varpar , only : maxpatch_pft, nlevdecomp
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! soil column filter
    type(cnstate_type)      , intent(in)    :: cnstate_vars
    type(carbonflux_type)   , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,pi,p,j               ! indices
    !-----------------------------------------------------------------------

    associate(                                                                                              & 
         ivt                                 =>    veg_pp%itype                                              , & ! Input:  [integer  (:)   ]  pft vegetation type                                
         wtcol                               =>    veg_pp%wtcol                                              , & ! Input:  [real(r8) (:)   ]  pft weight relative to column (0-1)               
         
         lf_flab                             =>    veg_vp%lf_flab                                     , & ! Input:  [real(r8) (:)   ]  leaf litter labile fraction                       
         lf_fcel                             =>    veg_vp%lf_fcel                                     , & ! Input:  [real(r8) (:)   ]  leaf litter cellulose fraction                    
         lf_flig                             =>    veg_vp%lf_flig                                     , & ! Input:  [real(r8) (:)   ]  leaf litter lignin fraction                       
         fr_flab                             =>    veg_vp%fr_flab                                     , & ! Input:  [real(r8) (:)   ]  fine root litter labile fraction                  
         fr_fcel                             =>    veg_vp%fr_fcel                                     , & ! Input:  [real(r8) (:)   ]  fine root litter cellulose fraction               
         fr_flig                             =>    veg_vp%fr_flig                                     , & ! Input:  [real(r8) (:)   ]  fine root litter lignin fraction                  
         
         leaf_prof                           =>    cnstate_vars%leaf_prof_patch                           , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof                          =>    cnstate_vars%froot_prof_patch                          , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
         croot_prof                          =>    cnstate_vars%croot_prof_patch                          , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
         stem_prof                           =>    cnstate_vars%stem_prof_patch                           , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          
         
         m_leafc_to_litter                   =>    carbonflux_vars%m_leafc_to_litter_patch                , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_to_litter                  =>    carbonflux_vars%m_frootc_to_litter_patch               , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_to_litter               =>    carbonflux_vars%m_livestemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_to_litter               =>    carbonflux_vars%m_deadstemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_to_litter              =>    carbonflux_vars%m_livecrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_to_litter              =>    carbonflux_vars%m_deadcrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafc_storage_to_litter           =>    carbonflux_vars%m_leafc_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_storage_to_litter          =>    carbonflux_vars%m_frootc_storage_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_storage_to_litter       =>    carbonflux_vars%m_livestemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_storage_to_litter       =>    carbonflux_vars%m_deadstemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_storage_to_litter      =>    carbonflux_vars%m_livecrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_storage_to_litter      =>    carbonflux_vars%m_deadcrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
         m_gresp_storage_to_litter           =>    carbonflux_vars%m_gresp_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafc_xfer_to_litter              =>    carbonflux_vars%m_leafc_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_xfer_to_litter             =>    carbonflux_vars%m_frootc_xfer_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_xfer_to_litter          =>    carbonflux_vars%m_livestemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_xfer_to_litter          =>    carbonflux_vars%m_deadstemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_xfer_to_litter         =>    carbonflux_vars%m_livecrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_xfer_to_litter         =>    carbonflux_vars%m_deadcrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_gresp_xfer_to_litter              =>    carbonflux_vars%m_gresp_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_cpool_to_litter                   =>    carbonflux_vars%m_cpool_to_litter_patch                , & ! Input:  [real(r8) (:)   ]               

         m_leafn_to_litter                   =>    nitrogenflux_vars%m_leafn_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_to_litter                  =>    nitrogenflux_vars%m_frootn_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_to_litter               =>    nitrogenflux_vars%m_livestemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_to_litter               =>    nitrogenflux_vars%m_deadstemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_to_litter              =>    nitrogenflux_vars%m_livecrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_to_litter              =>    nitrogenflux_vars%m_deadcrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_retransn_to_litter                =>    nitrogenflux_vars%m_retransn_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_npool_to_litter                   =>    nitrogenflux_vars%m_npool_to_litter_patch              , & ! Input:  [real(r8) (:)   ]
         m_leafn_storage_to_litter           =>    nitrogenflux_vars%m_leafn_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_storage_to_litter          =>    nitrogenflux_vars%m_frootn_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_storage_to_litter       =>    nitrogenflux_vars%m_livestemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_storage_to_litter       =>    nitrogenflux_vars%m_deadstemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_storage_to_litter      =>    nitrogenflux_vars%m_livecrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_storage_to_litter      =>    nitrogenflux_vars%m_deadcrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafn_xfer_to_litter              =>    nitrogenflux_vars%m_leafn_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_xfer_to_litter             =>    nitrogenflux_vars%m_frootn_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_xfer_to_litter          =>    nitrogenflux_vars%m_livestemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_xfer_to_litter          =>    nitrogenflux_vars%m_deadstemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_xfer_to_litter         =>    nitrogenflux_vars%m_livecrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_xfer_to_litter         =>    nitrogenflux_vars%m_deadcrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         
         !! add phosphorus  -X.YANG
         m_leafp_to_litter                   =>    phosphorusflux_vars%m_leafp_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootp_to_litter                  =>    phosphorusflux_vars%m_frootp_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemp_to_litter               =>    phosphorusflux_vars%m_livestemp_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemp_to_litter               =>    phosphorusflux_vars%m_deadstemp_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootp_to_litter              =>    phosphorusflux_vars%m_livecrootp_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootp_to_litter              =>    phosphorusflux_vars%m_deadcrootp_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_retransp_to_litter                =>    phosphorusflux_vars%m_retransp_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_ppool_to_litter                   =>    phosphorusflux_vars%m_ppool_to_litter_patch              , & ! Input:  [real(r8) (:)   ]         
         m_leafp_storage_to_litter           =>    phosphorusflux_vars%m_leafp_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootp_storage_to_litter          =>    phosphorusflux_vars%m_frootp_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemp_storage_to_litter       =>    phosphorusflux_vars%m_livestemp_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemp_storage_to_litter       =>    phosphorusflux_vars%m_deadstemp_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootp_storage_to_litter      =>    phosphorusflux_vars%m_livecrootp_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootp_storage_to_litter      =>    phosphorusflux_vars%m_deadcrootp_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafp_xfer_to_litter              =>    phosphorusflux_vars%m_leafp_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootp_xfer_to_litter             =>    phosphorusflux_vars%m_frootp_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemp_xfer_to_litter          =>    phosphorusflux_vars%m_livestemp_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemp_xfer_to_litter          =>    phosphorusflux_vars%m_deadstemp_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootp_xfer_to_litter         =>    phosphorusflux_vars%m_livecrootp_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootp_xfer_to_litter         =>    phosphorusflux_vars%m_deadcrootp_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    

         gap_mortality_c_to_litr_met_c       =>    carbonflux_vars%gap_mortality_c_to_litr_met_c_col      , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
         gap_mortality_c_to_litr_cel_c       =>    carbonflux_vars%gap_mortality_c_to_litr_cel_c_col      , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
         gap_mortality_c_to_litr_lig_c       =>    carbonflux_vars%gap_mortality_c_to_litr_lig_c_col      , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
         gap_mortality_c_to_cwdc             =>    carbonflux_vars%gap_mortality_c_to_cwdc_col            , & ! InOut:  [real(r8) (:,:) ]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
         
         gap_mortality_n_to_litr_met_n       =>    nitrogenflux_vars%gap_mortality_n_to_litr_met_n_col    , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
         gap_mortality_n_to_litr_cel_n       =>    nitrogenflux_vars%gap_mortality_n_to_litr_cel_n_col    , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
         gap_mortality_n_to_litr_lig_n       =>    nitrogenflux_vars%gap_mortality_n_to_litr_lig_n_col    , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
         gap_mortality_n_to_cwdn             =>    nitrogenflux_vars%gap_mortality_n_to_cwdn_col          ,  & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with gap mortality to CWD pool (gN/m3/s)

         gap_mortality_p_to_litr_met_p       =>    phosphorusflux_vars%gap_mortality_p_to_litr_met_p_col    , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
         gap_mortality_p_to_litr_cel_p       =>    phosphorusflux_vars%gap_mortality_p_to_litr_cel_p_col    , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
         gap_mortality_p_to_litr_lig_p       =>    phosphorusflux_vars%gap_mortality_p_to_litr_lig_p_col    , & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
         gap_mortality_p_to_cwdp             =>    phosphorusflux_vars%gap_mortality_p_to_cwdp_col            & ! InOut:  [real(r8) (:,:) ]  N fluxes associated with gap mortality to CWD pool (gN/m3/s)

         )

      do j = 1,nlevdecomp
         do pi = 1,maxpatch_pft
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (pi <=  col_pp%npfts(c)) then
                  p = col_pp%pfti(c) + pi - 1

                  if (veg_pp%active(p)) then

                     ! leaf gap mortality carbon fluxes
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          m_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_c_to_litr_cel_c(c,j) = gap_mortality_c_to_litr_cel_c(c,j) + &
                          m_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_c_to_litr_lig_c(c,j) = gap_mortality_c_to_litr_lig_c(c,j) + &
                          m_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! fine root gap mortality carbon fluxes
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          m_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_c_to_litr_cel_c(c,j) = gap_mortality_c_to_litr_cel_c(c,j) + &
                          m_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_c_to_litr_lig_c(c,j) = gap_mortality_c_to_litr_lig_c(c,j) + &
                          m_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! wood gap mortality carbon fluxes
                     gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                          (m_livestemc_to_litter(p) + m_deadstemc_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                          (m_livecrootc_to_litter(p) + m_deadcrootc_to_litter(p)) * wtcol(p) * croot_prof(p,j)
                     ! storage gap mortality carbon fluxes
                     gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                          (m_cpool_to_litter(p) + m_leafc_storage_to_litter(p) + m_gresp_storage_to_litter(p)) * wtcol(p)&
                          * leaf_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j)     = gap_mortality_c_to_litr_met_c(c,j)     + &
                          m_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                          (m_livestemc_storage_to_litter(p) + m_deadstemc_storage_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          (m_livecrootc_storage_to_litter(p) + m_deadcrootc_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! transfer gap mortality carbon fluxes
                     gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                          (m_leafc_xfer_to_litter(p) + m_gresp_xfer_to_litter(p))     * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j)     = gap_mortality_c_to_litr_met_c(c,j)     + &
                          m_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                          (m_livestemc_xfer_to_litter(p) + m_deadstemc_xfer_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          (m_livecrootc_xfer_to_litter(p) + m_deadcrootc_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! leaf gap mortality nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_cel_n(c,j) = gap_mortality_n_to_litr_cel_n(c,j) + &
                          m_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_lig_n(c,j) = gap_mortality_n_to_litr_lig_n(c,j) + &
                          m_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! fine root litter nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_cel_n(c,j) = gap_mortality_n_to_litr_cel_n(c,j) + &
                          m_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_lig_n(c,j) = gap_mortality_n_to_litr_lig_n(c,j) + &
                          m_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! wood gap mortality nitrogen fluxes
                     gap_mortality_n_to_cwdn(c,j)  = gap_mortality_n_to_cwdn(c,j)  + &
                          (m_livestemn_to_litter(p) + m_deadstemn_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j) + &
                          (m_livecrootn_to_litter(p) + m_deadcrootn_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! retranslocated N pool gap mortality fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                     ! storage N pool gap mortality fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_npool_to_litter(p) * wtcol(p) * leaf_prof(p,j)

                     ! storage gap mortality nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j)      = gap_mortality_n_to_litr_met_n(c,j)      + &
                          m_leafn_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j)     = gap_mortality_n_to_litr_met_n(c,j)     + &
                          m_frootn_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j)  = gap_mortality_n_to_litr_met_n(c,j)  + &
                          (m_livestemn_storage_to_litter(p) + m_deadstemn_storage_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livecrootn_storage_to_litter(p) + m_deadcrootn_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! transfer gap mortality nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j)      = gap_mortality_n_to_litr_met_n(c,j)      + &
                          m_leafn_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j)     = gap_mortality_n_to_litr_met_n(c,j)     + &
                          m_frootn_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j)  = gap_mortality_n_to_litr_met_n(c,j)  + &
                          (m_livestemn_xfer_to_litter(p) + m_deadstemn_xfer_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livecrootn_xfer_to_litter(p) + m_deadcrootn_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! leaf gap mortality phosphorus fluxes
                     gap_mortality_p_to_litr_met_p(c,j) = gap_mortality_p_to_litr_met_p(c,j) + &
                          m_leafp_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_p_to_litr_cel_p(c,j) = gap_mortality_p_to_litr_cel_p(c,j) + &
                          m_leafp_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_p_to_litr_lig_p(c,j) = gap_mortality_p_to_litr_lig_p(c,j) + &
                          m_leafp_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! fine root litter phosphorus fluxes
                     gap_mortality_p_to_litr_met_p(c,j) = gap_mortality_p_to_litr_met_p(c,j) + &
                          m_frootp_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_p_to_litr_cel_p(c,j) = gap_mortality_p_to_litr_cel_p(c,j) + &
                          m_frootp_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_p_to_litr_lig_p(c,j) = gap_mortality_p_to_litr_lig_p(c,j) + &
                          m_frootp_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! wood gap mortality phosphorus fluxes
                     gap_mortality_p_to_cwdp(c,j)  = gap_mortality_p_to_cwdp(c,j)  + &
                          (m_livestemp_to_litter(p) + m_deadstemp_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_p_to_cwdp(c,j) = gap_mortality_p_to_cwdp(c,j) + &
                          (m_livecrootp_to_litter(p) + m_deadcrootp_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! retranslocated N pool gap mortality fluxes
                     gap_mortality_p_to_litr_met_p(c,j) = gap_mortality_p_to_litr_met_p(c,j) + &
                          m_retransp_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_p_to_litr_met_p(c,j) = gap_mortality_p_to_litr_met_p(c,j) + &
                          m_ppool_to_litter(p) * wtcol(p) * leaf_prof(p,j)


                     ! storage gap mortality phosphorus fluxes
                     gap_mortality_p_to_litr_met_p(c,j)      = gap_mortality_p_to_litr_met_p(c,j)      + &
                          m_leafp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_p_to_litr_met_p(c,j)     = gap_mortality_p_to_litr_met_p(c,j)     + &
                          m_frootp_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                     gap_mortality_p_to_litr_met_p(c,j)  = gap_mortality_p_to_litr_met_p(c,j)  + &
                          (m_livestemp_storage_to_litter(p) + m_deadstemp_storage_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_p_to_litr_met_p(c,j) = gap_mortality_p_to_litr_met_p(c,j) + &
                          (m_livecrootp_storage_to_litter(p) + m_deadcrootp_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! transfer gap mortality phosphorus fluxes
                     gap_mortality_p_to_litr_met_p(c,j)      = gap_mortality_p_to_litr_met_p(c,j)      + &
                          m_leafp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_p_to_litr_met_p(c,j)     = gap_mortality_p_to_litr_met_p(c,j)     + &
                          m_frootp_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                     gap_mortality_p_to_litr_met_p(c,j)  = gap_mortality_p_to_litr_met_p(c,j)  + &
                          (m_livestemp_xfer_to_litter(p) + m_deadstemp_xfer_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_p_to_litr_met_p(c,j) = gap_mortality_p_to_litr_met_p(c,j) + &
                          (m_livecrootp_xfer_to_litter(p) + m_deadcrootp_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                  end if
               end if

            end do
         end do
      end do

    end associate

  end subroutine CNGapPftToColumn

  subroutine mortality_rate_soilorder(&
       num_soilp, filter_soilp, &
       cnstate_vars)
    !
    ! !DESCRIPTION:
    ! !this surroutine is to calculate mortality rate based on soil order

    ! USES
    use pftvarcon       , only: nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree
    use soilorder_varcon, only: r_mort_soilorder

    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(cnstate_type)       , intent(inout)    :: cnstate_vars

    ! local variables
    integer :: p,c,fp


    associate(                                                      &
       ivt            =>    veg_pp%itype                             , & ! Input:[integer  (:)   ]  patch vegetation type                                
       isoilorder     =>    cnstate_vars%isoilorder               , &
       r_mort_cal     =>    cnstate_vars%r_mort_cal_patch )

       ! loop over the patches
       do fp = 1,num_soilp
          p = filter_soilp(fp)
          c = veg_pp%column(p)
               if( veg_pp%itype(p) == nbrdlf_evr_trp_tree .or. veg_pp%itype(p) == nbrdlf_dcd_trp_tree )then
                   r_mort_cal(p) = r_mort_soilorder( isoilorder(c) )
               else
                   r_mort_cal(p) = 0.02_r8                 ! Default mortality rate 
               endif
       end do

     end associate


  end subroutine mortality_rate_soilorder


end module GapMortalityMod
