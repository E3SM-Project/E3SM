module CNGapMortalityMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in gap mortality for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod            , only : r8 => shr_kind_r8
  use decompMod               , only : bounds_type
  use abortutils              , only : endrun
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use pftconMod               , only : pftcon
  use CNDVType                , only : dgvs_type
  use CNVegCarbonStateType    , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType     , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType  , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType   , only : cnveg_nitrogenflux_type
  use ColumnType              , only : col                
  use PatchType               , only : patch                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams
  public :: CNGapMortality

  type, private :: params_type
     real(r8):: am     ! mortality rate based on annual rate, fractional mortality (1/yr)
     real(r8):: k_mort ! coeff. of growth efficiency in mortality equation
  end type params_type
  !
  type(params_type), private :: params_inst
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CNGap_PatchToColumn
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
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
    params_inst%am=tempr

    tString='k_mort'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    params_inst%k_mort=tempr   
    
  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine CNGapMortality (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       dgvs_inst, cnveg_carbonstate_inst, cnveg_nitrogenstate_inst, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
       leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch)
    !
    ! !DESCRIPTION:
    ! Gap-phase mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use clm_time_manager , only: get_days_per_year
    use clm_varpar       , only: nlevdecomp_full
    use clm_varcon       , only: secspday
    use clm_varctl       , only: use_cndv
    use pftconMod        , only: npcropmin
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! column filter for soil points
    integer                         , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    type(dgvs_type)                 , intent(inout) :: dgvs_inst
    type(cnveg_carbonstate_type)    , intent(in)    :: cnveg_carbonstate_inst
    type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    real(r8)                        , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
    !
    ! !LOCAL VARIABLES:
    integer :: p             ! patch index
    integer :: fp            ! patch filter index
    real(r8):: am            ! rate for fractional mortality (1/yr)
    real(r8):: m             ! rate for fractional mortality (1/s)
    real(r8):: mort_max      ! asymptotic max mortality rate (/yr)
    real(r8):: k_mort = 0.3  ! coeff of growth efficiency in mortality equation
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(leaf_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(froot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(croot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(stem_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), errMsg(__FILE__, __LINE__))

    associate(                                         & 
         ivt        =>    patch%itype                  , & ! Input:  [integer  (:) ]  patch vegetation type                                

         woody      =>    pftcon%woody               , & ! Input:  binary flag for woody lifeform                    
         
         greffic    =>    dgvs_inst%greffic_patch    , & ! Input:  [real(r8) (:) ]                                                    
         heatstress =>    dgvs_inst%heatstress_patch , & ! Input:  [real(r8) (:) ]                                                    
         nind       =>    dgvs_inst%nind_patch         & ! Output: [real(r8) (:) ]  number of individuals (#/m2) added by F. Li and S. Levis
         )

      ! set the mortality rate based on annual rate
      am = params_inst%am
      ! set coeff of growth efficiency in mortality equation 
      k_mort = params_inst%k_mort

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
               am = params_inst%am
            end if

         end if

         m  = am/(get_days_per_year() * secspday)

         !------------------------------------------------------
         ! patch-level gap mortality carbon fluxes
         !------------------------------------------------------

         ! displayed pools
         cnveg_carbonflux_inst%m_leafc_to_litter_patch(p)               = cnveg_carbonstate_inst%leafc_patch(p)               * m
         cnveg_carbonflux_inst%m_frootc_to_litter_patch(p)              = cnveg_carbonstate_inst%frootc_patch(p)              * m
         cnveg_carbonflux_inst%m_livestemc_to_litter_patch(p)           = cnveg_carbonstate_inst%livestemc_patch(p)           * m
         cnveg_carbonflux_inst%m_deadstemc_to_litter_patch(p)           = cnveg_carbonstate_inst%deadstemc_patch(p)           * m
         cnveg_carbonflux_inst%m_livecrootc_to_litter_patch(p)          = cnveg_carbonstate_inst%livecrootc_patch(p)          * m
         cnveg_carbonflux_inst%m_deadcrootc_to_litter_patch(p)          = cnveg_carbonstate_inst%deadcrootc_patch(p)          * m

         ! storage pools
         cnveg_carbonflux_inst%m_leafc_storage_to_litter_patch(p)       = cnveg_carbonstate_inst%leafc_storage_patch(p)       * m
         cnveg_carbonflux_inst%m_frootc_storage_to_litter_patch(p)      = cnveg_carbonstate_inst%frootc_storage_patch(p)      * m
         cnveg_carbonflux_inst%m_livestemc_storage_to_litter_patch(p)   = cnveg_carbonstate_inst%livestemc_storage_patch(p)   * m
         cnveg_carbonflux_inst%m_deadstemc_storage_to_litter_patch(p)   = cnveg_carbonstate_inst%deadstemc_storage_patch(p)   * m
         cnveg_carbonflux_inst%m_livecrootc_storage_to_litter_patch(p)  = cnveg_carbonstate_inst%livecrootc_storage_patch(p)  * m
         cnveg_carbonflux_inst%m_deadcrootc_storage_to_litter_patch(p)  = cnveg_carbonstate_inst%deadcrootc_storage_patch(p)  * m
         cnveg_carbonflux_inst%m_gresp_storage_to_litter_patch(p)       = cnveg_carbonstate_inst%gresp_storage_patch(p)       * m

         ! transfer pools
         cnveg_carbonflux_inst%m_leafc_xfer_to_litter_patch(p)          = cnveg_carbonstate_inst%leafc_xfer_patch(p)          * m
         cnveg_carbonflux_inst%m_frootc_xfer_to_litter_patch(p)         = cnveg_carbonstate_inst%frootc_xfer_patch(p)         * m
         cnveg_carbonflux_inst%m_livestemc_xfer_to_litter_patch(p)      = cnveg_carbonstate_inst%livestemc_xfer_patch(p)      * m
         cnveg_carbonflux_inst%m_deadstemc_xfer_to_litter_patch(p)      = cnveg_carbonstate_inst%deadstemc_xfer_patch(p)      * m
         cnveg_carbonflux_inst%m_livecrootc_xfer_to_litter_patch(p)     = cnveg_carbonstate_inst%livecrootc_xfer_patch(p)     * m
         cnveg_carbonflux_inst%m_deadcrootc_xfer_to_litter_patch(p)     = cnveg_carbonstate_inst%deadcrootc_xfer_patch(p)     * m
         cnveg_carbonflux_inst%m_gresp_xfer_to_litter_patch(p)          = cnveg_carbonstate_inst%gresp_xfer_patch(p)          * m

         !------------------------------------------------------
         ! patch-level gap mortality nitrogen fluxes
         !------------------------------------------------------

         ! displayed pools
         cnveg_nitrogenflux_inst%m_leafn_to_litter_patch(p)               = cnveg_nitrogenstate_inst%leafn_patch(p)      * m
         cnveg_nitrogenflux_inst%m_frootn_to_litter_patch(p)              = cnveg_nitrogenstate_inst%frootn_patch(p)     * m
         cnveg_nitrogenflux_inst%m_livestemn_to_litter_patch(p)           = cnveg_nitrogenstate_inst%livestemn_patch(p)  * m
         cnveg_nitrogenflux_inst%m_deadstemn_to_litter_patch(p)           = cnveg_nitrogenstate_inst%deadstemn_patch(p)  * m
         cnveg_nitrogenflux_inst%m_livecrootn_to_litter_patch(p)          = cnveg_nitrogenstate_inst%livecrootn_patch(p) * m
         cnveg_nitrogenflux_inst%m_deadcrootn_to_litter_patch(p)          = cnveg_nitrogenstate_inst%deadcrootn_patch(p) * m
         if (ivt(p) < npcropmin) then
            cnveg_nitrogenflux_inst%m_retransn_to_litter_patch(p) = cnveg_nitrogenstate_inst%retransn_patch(p) * m
         end if
            
         ! storage pools
         cnveg_nitrogenflux_inst%m_leafn_storage_to_litter_patch(p)       = cnveg_nitrogenstate_inst%leafn_storage_patch(p)      * m
         cnveg_nitrogenflux_inst%m_frootn_storage_to_litter_patch(p)      = cnveg_nitrogenstate_inst%frootn_storage_patch(p)     * m
         cnveg_nitrogenflux_inst%m_livestemn_storage_to_litter_patch(p)   = cnveg_nitrogenstate_inst%livestemn_storage_patch(p)  * m
         cnveg_nitrogenflux_inst%m_deadstemn_storage_to_litter_patch(p)   = cnveg_nitrogenstate_inst%deadstemn_storage_patch(p)  * m
         cnveg_nitrogenflux_inst%m_livecrootn_storage_to_litter_patch(p)  = cnveg_nitrogenstate_inst%livecrootn_storage_patch(p) * m
         cnveg_nitrogenflux_inst%m_deadcrootn_storage_to_litter_patch(p)  = cnveg_nitrogenstate_inst%deadcrootn_storage_patch(p) * m

         ! transfer pools
         cnveg_nitrogenflux_inst%m_leafn_xfer_to_litter_patch(p)          = cnveg_nitrogenstate_inst%leafn_xfer_patch(p)         * m
         cnveg_nitrogenflux_inst%m_frootn_xfer_to_litter_patch(p)         = cnveg_nitrogenstate_inst%frootn_xfer_patch(p)        * m
         cnveg_nitrogenflux_inst%m_livestemn_xfer_to_litter_patch(p)      = cnveg_nitrogenstate_inst%livestemn_xfer_patch(p)     * m
         cnveg_nitrogenflux_inst%m_deadstemn_xfer_to_litter_patch(p)      = cnveg_nitrogenstate_inst%deadstemn_xfer_patch(p)     * m
         cnveg_nitrogenflux_inst%m_livecrootn_xfer_to_litter_patch(p)     = cnveg_nitrogenstate_inst%livecrootn_xfer_patch(p)    * m
         cnveg_nitrogenflux_inst%m_deadcrootn_xfer_to_litter_patch(p)     = cnveg_nitrogenstate_inst%deadcrootn_xfer_patch(p)    * m

         ! added by F. Li and S. Levis
         if (use_cndv) then
            if (woody(ivt(p)) == 1._r8)then
               if (cnveg_carbonstate_inst%livestemc_patch(p) + cnveg_carbonstate_inst%deadstemc_patch(p)> 0._r8)then
                  nind(p)=nind(p)*(1._r8-m)
               else
                  nind(p) = 0._r8 
               end if
            end if
         end if

      end do ! end of patch loop

      ! gather all patch-level litterfall fluxes to the column
      ! for litter C and N inputs

      call CNGap_PatchToColumn(bounds, num_soilc, filter_soilc, &
           cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
           leaf_prof_patch(bounds%begp:bounds%endp, 1:nlevdecomp_full), &
           froot_prof_patch(bounds%begp:bounds%endp, 1:nlevdecomp_full), & 
           croot_prof_patch(bounds%begp:bounds%endp, 1:nlevdecomp_full), &
           stem_prof_patch(bounds%begp:bounds%endp, 1:nlevdecomp_full))

    end associate

  end subroutine CNGapMortality

  !-----------------------------------------------------------------------
  subroutine CNGap_PatchToColumn (bounds, num_soilc, filter_soilc, &
       cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
       leaf_prof_patch, froot_prof_patch, croot_prof_patch, stem_prof_patch)
    !
    ! !DESCRIPTION:
    ! gathers all patch-level gap mortality fluxes to the column level and
    ! assigns them to the three litter pools
    !
    ! !USES:
    use clm_varpar , only : maxpatch_pft, nlevdecomp, nlevdecomp_full
    !
    ! !ARGUMENTS:
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                         , intent(in)    :: filter_soilc(:) ! soil column filter
    type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
    real(r8)                        , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: croot_prof_patch(bounds%begp:,1:)
    real(r8)                        , intent(in)    :: stem_prof_patch(bounds%begp:,1:)
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,pi,p,j               ! indices
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(leaf_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(froot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(croot_prof_patch)  == (/bounds%endp,nlevdecomp_full/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(stem_prof_patch)   == (/bounds%endp,nlevdecomp_full/)), errMsg(__FILE__, __LINE__))

    associate(                                                                                                 & 
         leaf_prof                           => leaf_prof_patch                                              , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves                         
         froot_prof                          => froot_prof_patch                                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
         croot_prof                          => croot_prof_patch                                             , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of coarse roots                   
         stem_prof                           => stem_prof_patch                                              , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of stems                          

         ivt                                 => patch%itype                                                    , & ! Input:  [integer  (:)   ]  patch vegetation type                                
         wtcol                               => patch%wtcol                                                    , & ! Input:  [real(r8) (:)   ]  patch weight relative to column (0-1)               
         
         lf_flab                             => pftcon%lf_flab                                               , & ! Input:  [real(r8) (:)   ]  leaf litter labile fraction                       
         lf_fcel                             => pftcon%lf_fcel                                               , & ! Input:  [real(r8) (:)   ]  leaf litter cellulose fraction                    
         lf_flig                             => pftcon%lf_flig                                               , & ! Input:  [real(r8) (:)   ]  leaf litter lignin fraction                       
         fr_flab                             => pftcon%fr_flab                                               , & ! Input:  [real(r8) (:)   ]  fine root litter labile fraction                  
         fr_fcel                             => pftcon%fr_fcel                                               , & ! Input:  [real(r8) (:)   ]  fine root litter cellulose fraction               
         fr_flig                             => pftcon%fr_flig                                               , & ! Input:  [real(r8) (:)   ]  fine root litter lignin fraction                  
         
         m_leafc_to_litter                   => cnveg_carbonflux_inst%m_leafc_to_litter_patch                , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_to_litter                  => cnveg_carbonflux_inst%m_frootc_to_litter_patch               , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_to_litter               => cnveg_carbonflux_inst%m_livestemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_to_litter               => cnveg_carbonflux_inst%m_deadstemc_to_litter_patch            , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_to_litter              => cnveg_carbonflux_inst%m_livecrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_to_litter              => cnveg_carbonflux_inst%m_deadcrootc_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafc_storage_to_litter           => cnveg_carbonflux_inst%m_leafc_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_storage_to_litter          => cnveg_carbonflux_inst%m_frootc_storage_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_storage_to_litter       => cnveg_carbonflux_inst%m_livestemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_storage_to_litter       => cnveg_carbonflux_inst%m_deadstemc_storage_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_storage_to_litter      => cnveg_carbonflux_inst%m_livecrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_storage_to_litter      => cnveg_carbonflux_inst%m_deadcrootc_storage_to_litter_patch   , & ! Input:  [real(r8) (:)   ]                                                    
         m_gresp_storage_to_litter           => cnveg_carbonflux_inst%m_gresp_storage_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafc_xfer_to_litter              => cnveg_carbonflux_inst%m_leafc_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootc_xfer_to_litter             => cnveg_carbonflux_inst%m_frootc_xfer_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemc_xfer_to_litter          => cnveg_carbonflux_inst%m_livestemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemc_xfer_to_litter          => cnveg_carbonflux_inst%m_deadstemc_xfer_to_litter_patch       , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootc_xfer_to_litter         => cnveg_carbonflux_inst%m_livecrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootc_xfer_to_litter         => cnveg_carbonflux_inst%m_deadcrootc_xfer_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_gresp_xfer_to_litter              => cnveg_carbonflux_inst%m_gresp_xfer_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         gap_mortality_c_to_litr_met_c       => cnveg_carbonflux_inst%gap_mortality_c_to_litr_met_c_col      , & ! Output: [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
         gap_mortality_c_to_litr_cel_c       => cnveg_carbonflux_inst%gap_mortality_c_to_litr_cel_c_col      , & ! Output: [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
         gap_mortality_c_to_litr_lig_c       => cnveg_carbonflux_inst%gap_mortality_c_to_litr_lig_c_col      , & ! Output: [real(r8) (:,:) ]  C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
         gap_mortality_c_to_cwdc             => cnveg_carbonflux_inst%gap_mortality_c_to_cwdc_col            , & ! Output: [real(r8) (:,:) ]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
         
         m_leafn_to_litter                   => cnveg_nitrogenflux_inst%m_leafn_to_litter_patch              , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_to_litter                  => cnveg_nitrogenflux_inst%m_frootn_to_litter_patch             , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_to_litter               => cnveg_nitrogenflux_inst%m_livestemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_to_litter               => cnveg_nitrogenflux_inst%m_deadstemn_to_litter_patch          , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_to_litter              => cnveg_nitrogenflux_inst%m_livecrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_to_litter              => cnveg_nitrogenflux_inst%m_deadcrootn_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_retransn_to_litter                => cnveg_nitrogenflux_inst%m_retransn_to_litter_patch           , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafn_storage_to_litter           => cnveg_nitrogenflux_inst%m_leafn_storage_to_litter_patch      , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_storage_to_litter          => cnveg_nitrogenflux_inst%m_frootn_storage_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_storage_to_litter       => cnveg_nitrogenflux_inst%m_livestemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_storage_to_litter       => cnveg_nitrogenflux_inst%m_deadstemn_storage_to_litter_patch  , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_storage_to_litter      => cnveg_nitrogenflux_inst%m_livecrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_storage_to_litter      => cnveg_nitrogenflux_inst%m_deadcrootn_storage_to_litter_patch , & ! Input:  [real(r8) (:)   ]                                                    
         m_leafn_xfer_to_litter              => cnveg_nitrogenflux_inst%m_leafn_xfer_to_litter_patch         , & ! Input:  [real(r8) (:)   ]                                                    
         m_frootn_xfer_to_litter             => cnveg_nitrogenflux_inst%m_frootn_xfer_to_litter_patch        , & ! Input:  [real(r8) (:)   ]                                                    
         m_livestemn_xfer_to_litter          => cnveg_nitrogenflux_inst%m_livestemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadstemn_xfer_to_litter          => cnveg_nitrogenflux_inst%m_deadstemn_xfer_to_litter_patch     , & ! Input:  [real(r8) (:)   ]                                                    
         m_livecrootn_xfer_to_litter         => cnveg_nitrogenflux_inst%m_livecrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         m_deadcrootn_xfer_to_litter         => cnveg_nitrogenflux_inst%m_deadcrootn_xfer_to_litter_patch    , & ! Input:  [real(r8) (:)   ]                                                    
         gap_mortality_n_to_litr_met_n       => cnveg_nitrogenflux_inst%gap_mortality_n_to_litr_met_n_col    , & ! Output: [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
         gap_mortality_n_to_litr_cel_n       => cnveg_nitrogenflux_inst%gap_mortality_n_to_litr_cel_n_col    , & ! Output: [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
         gap_mortality_n_to_litr_lig_n       => cnveg_nitrogenflux_inst%gap_mortality_n_to_litr_lig_n_col    , & ! Output: [real(r8) (:,:) ]  N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
         gap_mortality_n_to_cwdn             => cnveg_nitrogenflux_inst%gap_mortality_n_to_cwdn_col            & ! Output: [real(r8) (:,:) ]  N fluxes associated with gap mortality to CWD pool (gN/m3/s)
         )

      do j = 1,nlevdecomp
         do pi = 1,maxpatch_pft
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (pi <=  col%npatches(c)) then
                  p = col%patchi(c) + pi - 1

                  if (patch%active(p)) then

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
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          (m_leafc_storage_to_litter(p) + m_gresp_storage_to_litter(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          m_frootc_storage_to_litter(p) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j)  + &
                          (m_livestemc_storage_to_litter(p) + m_deadstemc_storage_to_litter(p)) * wtcol(p) * stem_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          (m_livecrootc_storage_to_litter(p) + m_deadcrootc_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! transfer gap mortality carbon fluxes
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          (m_leafc_xfer_to_litter(p) + m_gresp_xfer_to_litter(p)) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                          m_frootc_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j)
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
                     gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j)  + &
                          (m_livestemn_to_litter(p) + m_deadstemn_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j) + &
                          (m_livecrootn_to_litter(p) + m_deadcrootn_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! retranslocated N pool gap mortality fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)

                     ! storage gap mortality nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_leafn_storage_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_frootn_storage_to_litter(p) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j)  = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livestemn_storage_to_litter(p) + m_deadstemn_storage_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livecrootn_storage_to_litter(p) + m_deadcrootn_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j)

                     ! transfer gap mortality nitrogen fluxes
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_leafn_xfer_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          m_frootn_xfer_to_litter(p) * wtcol(p) * froot_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livestemn_xfer_to_litter(p) + m_deadstemn_xfer_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                     gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                          (m_livecrootn_xfer_to_litter(p) + m_deadcrootn_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)


                  end if
               end if

            end do
         end do
      end do

    end associate

  end subroutine CNGap_PatchToColumn

end module CNGapMortalityMod
