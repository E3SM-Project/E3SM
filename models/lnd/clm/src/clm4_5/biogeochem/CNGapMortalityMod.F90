module CNGapMortalityMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in gap mortality for coupled carbon
  ! nitrogen code.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils  , only: endrun
  use shr_log_mod , only: errMsg => shr_log_errMsg
  implicit none
  save
  private

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNGapMortality
  public :: readCNGapMortParams

  type, private :: CNGapMortParamsType
      real(r8):: am     ! mortality rate based on annual rate, fractional mortality (1/yr)
      real(r8):: k_mort ! coeff. of growth efficiency in mortality equation
  end type CNGapMortParamsType

  type(CNGapMortParamsType),private ::  CNGapMortParamsInst
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readCNGapMortParams ( ncid )
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
    
  end subroutine readCNGapMortParams

  !-----------------------------------------------------------------------
  subroutine CNGapMortality (num_soilc, filter_soilc, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! Gap-phase mortality routine for coupled carbon-nitrogen code (CN)
    !
    ! !USES:
    use clmtype
    use clm_time_manager, only: get_days_per_year
    use clm_varcon      , only: secspday
    use pftvarcon       , only: npcropmin
    use clm_varctl      , only: use_cndv
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:) ! column filter for soil points
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:) ! pft filter for soil points
    !
    ! !LOCAL VARIABLES:
    integer :: p             ! pft index
    integer :: fp            ! pft filter index
    real(r8):: am            ! rate for fractional mortality (1/yr)
    real(r8):: m             ! rate for fractional mortality (1/s)
    real(r8):: mort_max      ! asymptotic max mortality rate (/yr)
    real(r8):: k_mort = 0.3  ! coeff of growth efficiency in mortality equation
    !-----------------------------------------------------------------------

   associate(& 
   woody                               =>    pftcon%woody                                , & ! Input:  [real(r8) (:)]  binary flag for woody lifeform                    
   ivt                                 =>    pft%itype                                   , & ! Input:  [integer (:)]  pft vegetation type                                
   leafc                               =>    pcs%leafc                                   , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C                                    
   frootc                              =>    pcs%frootc                                  , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C                               
   livestemc                           =>    pcs%livestemc                               , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C                               
   deadstemc                           =>    pcs%deadstemc                               , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C                               
   livecrootc                          =>    pcs%livecrootc                              , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C                        
   deadcrootc                          =>    pcs%deadcrootc                              , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C                        
   leafc_storage                       =>    pcs%leafc_storage                           , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C storage                            
   frootc_storage                      =>    pcs%frootc_storage                          , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C storage                       
   livestemc_storage                   =>    pcs%livestemc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C storage                       
   deadstemc_storage                   =>    pcs%deadstemc_storage                       , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C storage                       
   livecrootc_storage                  =>    pcs%livecrootc_storage                      , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C storage                
   deadcrootc_storage                  =>    pcs%deadcrootc_storage                      , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C storage                
   gresp_storage                       =>    pcs%gresp_storage                           , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration storage                
   leafc_xfer                          =>    pcs%leafc_xfer                              , & ! Input:  [real(r8) (:)]  (gC/m2) leaf C transfer                           
   frootc_xfer                         =>    pcs%frootc_xfer                             , & ! Input:  [real(r8) (:)]  (gC/m2) fine root C transfer                      
   livestemc_xfer                      =>    pcs%livestemc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) live stem C transfer                      
   deadstemc_xfer                      =>    pcs%deadstemc_xfer                          , & ! Input:  [real(r8) (:)]  (gC/m2) dead stem C transfer                      
   livecrootc_xfer                     =>    pcs%livecrootc_xfer                         , & ! Input:  [real(r8) (:)]  (gC/m2) live coarse root C transfer               
   deadcrootc_xfer                     =>    pcs%deadcrootc_xfer                         , & ! Input:  [real(r8) (:)]  (gC/m2) dead coarse root C transfer               
   gresp_xfer                          =>    pcs%gresp_xfer                              , & ! Input:  [real(r8) (:)]  (gC/m2) growth respiration transfer               
   leafn                               =>    pns%leafn                                   , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N                                    
   frootn                              =>    pns%frootn                                  , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N                               
   livestemn                           =>    pns%livestemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N                               
   deadstemn                           =>    pns%deadstemn                               , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N                               
   livecrootn                          =>    pns%livecrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N                        
   deadcrootn                          =>    pns%deadcrootn                              , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N                        
   retransn                            =>    pns%retransn                                , & ! Input:  [real(r8) (:)]  (gN/m2) plant pool of retranslocated N            
   leafn_storage                       =>    pns%leafn_storage                           , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N storage                            
   frootn_storage                      =>    pns%frootn_storage                          , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N storage                       
   livestemn_storage                   =>    pns%livestemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N storage                       
   deadstemn_storage                   =>    pns%deadstemn_storage                       , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N storage                       
   livecrootn_storage                  =>    pns%livecrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N storage                
   deadcrootn_storage                  =>    pns%deadcrootn_storage                      , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N storage                
   leafn_xfer                          =>    pns%leafn_xfer                              , & ! Input:  [real(r8) (:)]  (gN/m2) leaf N transfer                           
   frootn_xfer                         =>    pns%frootn_xfer                             , & ! Input:  [real(r8) (:)]  (gN/m2) fine root N transfer                      
   livestemn_xfer                      =>    pns%livestemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) live stem N transfer                      
   deadstemn_xfer                      =>    pns%deadstemn_xfer                          , & ! Input:  [real(r8) (:)]  (gN/m2) dead stem N transfer                      
   livecrootn_xfer                     =>    pns%livecrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) live coarse root N transfer               
   deadcrootn_xfer                     =>    pns%deadcrootn_xfer                         , & ! Input:  [real(r8) (:)]  (gN/m2) dead coarse root N transfer               
   m_leafc_to_litter                   =>    pcf%m_leafc_to_litter                       , & ! Output: [real(r8) (:)]                                                    
   m_frootc_to_litter                  =>    pcf%m_frootc_to_litter                      , & ! Output: [real(r8) (:)]                                                    
   m_livestemc_to_litter               =>    pcf%m_livestemc_to_litter                   , & ! Output: [real(r8) (:)]                                                    
   m_deadstemc_to_litter               =>    pcf%m_deadstemc_to_litter                   , & ! Output: [real(r8) (:)]                                                    
   m_livecrootc_to_litter              =>    pcf%m_livecrootc_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   m_deadcrootc_to_litter              =>    pcf%m_deadcrootc_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   m_leafc_storage_to_litter           =>    pcf%m_leafc_storage_to_litter               , & ! Output: [real(r8) (:)]                                                    
   m_frootc_storage_to_litter          =>    pcf%m_frootc_storage_to_litter              , & ! Output: [real(r8) (:)]                                                    
   m_livestemc_storage_to_litter       =>    pcf%m_livestemc_storage_to_litter           , & ! Output: [real(r8) (:)]                                                    
   m_deadstemc_storage_to_litter       =>    pcf%m_deadstemc_storage_to_litter           , & ! Output: [real(r8) (:)]                                                    
   m_livecrootc_storage_to_litter      =>    pcf%m_livecrootc_storage_to_litter          , & ! Output: [real(r8) (:)]                                                    
   m_deadcrootc_storage_to_litter      =>    pcf%m_deadcrootc_storage_to_litter          , & ! Output: [real(r8) (:)]                                                    
   m_gresp_storage_to_litter           =>    pcf%m_gresp_storage_to_litter               , & ! Output: [real(r8) (:)]                                                    
   m_leafc_xfer_to_litter              =>    pcf%m_leafc_xfer_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   m_frootc_xfer_to_litter             =>    pcf%m_frootc_xfer_to_litter                 , & ! Output: [real(r8) (:)]                                                    
   m_livestemc_xfer_to_litter          =>    pcf%m_livestemc_xfer_to_litter              , & ! Output: [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_litter          =>    pcf%m_deadstemc_xfer_to_litter              , & ! Output: [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_litter         =>    pcf%m_livecrootc_xfer_to_litter             , & ! Output: [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_litter         =>    pcf%m_deadcrootc_xfer_to_litter             , & ! Output: [real(r8) (:)]                                                    
   m_gresp_xfer_to_litter              =>    pcf%m_gresp_xfer_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   m_leafn_to_litter                   =>    pnf%m_leafn_to_litter                       , & ! Output: [real(r8) (:)]                                                    
   m_frootn_to_litter                  =>    pnf%m_frootn_to_litter                      , & ! Output: [real(r8) (:)]                                                    
   m_livestemn_to_litter               =>    pnf%m_livestemn_to_litter                   , & ! Output: [real(r8) (:)]                                                    
   m_deadstemn_to_litter               =>    pnf%m_deadstemn_to_litter                   , & ! Output: [real(r8) (:)]                                                    
   m_livecrootn_to_litter              =>    pnf%m_livecrootn_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   m_deadcrootn_to_litter              =>    pnf%m_deadcrootn_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   m_retransn_to_litter                =>    pnf%m_retransn_to_litter                    , & ! Output: [real(r8) (:)]                                                    
   m_leafn_storage_to_litter           =>    pnf%m_leafn_storage_to_litter               , & ! Output: [real(r8) (:)]                                                    
   m_frootn_storage_to_litter          =>    pnf%m_frootn_storage_to_litter              , & ! Output: [real(r8) (:)]                                                    
   m_livestemn_storage_to_litter       =>    pnf%m_livestemn_storage_to_litter           , & ! Output: [real(r8) (:)]                                                    
   m_deadstemn_storage_to_litter       =>    pnf%m_deadstemn_storage_to_litter           , & ! Output: [real(r8) (:)]                                                    
   m_livecrootn_storage_to_litter      =>    pnf%m_livecrootn_storage_to_litter          , & ! Output: [real(r8) (:)]                                                    
   m_deadcrootn_storage_to_litter      =>    pnf%m_deadcrootn_storage_to_litter          , & ! Output: [real(r8) (:)]                                                    
   m_leafn_xfer_to_litter              =>    pnf%m_leafn_xfer_to_litter                  , & ! Output: [real(r8) (:)]                                                    
   m_frootn_xfer_to_litter             =>    pnf%m_frootn_xfer_to_litter                 , & ! Output: [real(r8) (:)]                                                    
   m_livestemn_xfer_to_litter          =>    pnf%m_livestemn_xfer_to_litter              , & ! Output: [real(r8) (:)]                                                    
   m_deadstemn_xfer_to_litter          =>    pnf%m_deadstemn_xfer_to_litter              , & ! Output: [real(r8) (:)]                                                    
   m_livecrootn_xfer_to_litter         =>    pnf%m_livecrootn_xfer_to_litter             , & ! Output: [real(r8) (:)]                                                    
   greffic                             =>    pdgvs%greffic                               , & ! Input:  [real(r8) (:)]                                                    
   heatstress                          =>    pdgvs%heatstress                            , & ! Input:  [real(r8) (:)]                                                    
   nind                                =>    pdgvs%nind                                  , & ! Input:  [real(r8) (:)]  number of individuals (#/m2) added by F. Li and S. Levis
   m_deadcrootn_xfer_to_litter         =>    pnf%m_deadcrootn_xfer_to_litter               & ! Output: [real(r8) (:)]                                                    
   )

   ! set the mortality rate based on annual rate
   am = CNGapMortParamsInst%am
   ! set coeff of growth efficiency in mortality equation 
   k_mort = CNGapMortParamsInst%k_mort

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      if (use_cndv) then
         ! Stress mortality from lpj's subr Mortality.
         
         if (woody(ivt(p)) == 1._r8) then

            if (ivt(p) == 8) then
               mort_max = 0.03_r8 ! BDT boreal
            else
               mort_max = 0.01_r8 ! original value for all pfts
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

      m  = am/(get_days_per_year() * secspday)

      ! pft-level gap mortality carbon fluxes
      ! displayed pools
      m_leafc_to_litter(p)               = leafc(p)               * m
      m_frootc_to_litter(p)              = frootc(p)              * m
      m_livestemc_to_litter(p)           = livestemc(p)           * m
      m_deadstemc_to_litter(p)           = deadstemc(p)           * m
      m_livecrootc_to_litter(p)          = livecrootc(p)          * m
      m_deadcrootc_to_litter(p)          = deadcrootc(p)          * m

      ! storage pools
      m_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
      m_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
      m_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
      m_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
      m_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
      m_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
      m_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

      ! transfer pools
      m_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
      m_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
      m_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
      m_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
      m_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
      m_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
      m_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

      ! pft-level gap mortality nitrogen fluxes
      ! displayed pools
      m_leafn_to_litter(p)               = leafn(p)               * m
      m_frootn_to_litter(p)              = frootn(p)              * m
      m_livestemn_to_litter(p)           = livestemn(p)           * m
      m_deadstemn_to_litter(p)           = deadstemn(p)           * m
      m_livecrootn_to_litter(p)          = livecrootn(p)          * m
      m_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
      if (ivt(p) < npcropmin) m_retransn_to_litter(p) = retransn(p) * m

      ! storage pools
      m_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
      m_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
      m_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
      m_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
      m_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
      m_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

      ! transfer pools
      m_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
      m_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
      m_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
      m_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
      m_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
      m_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m

      ! added by F. Li and S. Levis
      if (use_cndv) then
         if (woody(ivt(p)) == 1._r8)then
            if (livestemc(p)+deadstemc(p)> 0._r8)then
               nind(p)=nind(p)*(1._r8-m)
            else
               nind(p) = 0._r8 
            end if
         end if
      end if

   end do ! end of pft loop

   ! gather all pft-level litterfall fluxes to the column
   ! for litter C and N inputs

   call CNGapPftToColumn(num_soilc, filter_soilc)

   end associate
end subroutine CNGapMortality

!-----------------------------------------------------------------------
subroutine CNGapPftToColumn (num_soilc, filter_soilc)
  !
  ! !DESCRIPTION:
  ! called in the middle of CNGapMoratlity to gather all pft-level gap mortality fluxes
  ! to the column level and assign them to the three litter pools
  !
  ! !USES:
  use clmtype
  use clm_varpar, only : maxpatch_pft, nlevdecomp
  !
  ! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! soil column filter
  !
  ! !LOCAL VARIABLES:
  integer :: fc,c,pi,p,j               ! indices
  !-----------------------------------------------------------------------

   associate(& 
   lf_flab                             =>    pftcon%lf_flab                              , & ! Input:  [real(r8) (:)]  leaf litter labile fraction                       
   lf_fcel                             =>    pftcon%lf_fcel                              , & ! Input:  [real(r8) (:)]  leaf litter cellulose fraction                    
   lf_flig                             =>    pftcon%lf_flig                              , & ! Input:  [real(r8) (:)]  leaf litter lignin fraction                       
   fr_flab                             =>    pftcon%fr_flab                              , & ! Input:  [real(r8) (:)]  fine root litter labile fraction                  
   fr_fcel                             =>    pftcon%fr_fcel                              , & ! Input:  [real(r8) (:)]  fine root litter cellulose fraction               
   fr_flig                             =>    pftcon%fr_flig                              , & ! Input:  [real(r8) (:)]  fine root litter lignin fraction                  
   npfts                               =>    col%npfts                                   , & ! Input:  [integer (:)]  number of pfts for each column                     
   pfti                                =>    col%pfti                                    , & ! Input:  [integer (:)]  beginning pft index for each column                
   pactive                             =>    pft%active                                  , & ! Input:  [logical (:)]  true=>do computations on this pft 
   ivt                                 =>    pft%itype                                   , & ! Input:  [integer (:)]  pft vegetation type                                
   wtcol                               =>    pft%wtcol                                   , & ! Input:  [real(r8) (:)]  pft weight relative to column (0-1)               
   m_leafc_to_litter                   =>    pcf%m_leafc_to_litter                       , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_to_litter                  =>    pcf%m_frootc_to_litter                      , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_litter               =>    pcf%m_livestemc_to_litter                   , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_to_litter               =>    pcf%m_deadstemc_to_litter                   , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_litter              =>    pcf%m_livecrootc_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_to_litter              =>    pcf%m_deadcrootc_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_storage_to_litter           =>    pcf%m_leafc_storage_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_storage_to_litter          =>    pcf%m_frootc_storage_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_storage_to_litter       =>    pcf%m_livestemc_storage_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_storage_to_litter       =>    pcf%m_deadstemc_storage_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_storage_to_litter      =>    pcf%m_livecrootc_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_storage_to_litter      =>    pcf%m_deadcrootc_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_storage_to_litter           =>    pcf%m_gresp_storage_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_xfer_to_litter              =>    pcf%m_leafc_xfer_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_xfer_to_litter             =>    pcf%m_frootc_xfer_to_litter                 , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_xfer_to_litter          =>    pcf%m_livestemc_xfer_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_litter          =>    pcf%m_deadstemc_xfer_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_litter         =>    pcf%m_livecrootc_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_litter         =>    pcf%m_deadcrootc_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_xfer_to_litter              =>    pcf%m_gresp_xfer_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_to_litter                   =>    pnf%m_leafn_to_litter                       , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_to_litter                  =>    pnf%m_frootn_to_litter                      , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_to_litter               =>    pnf%m_livestemn_to_litter                   , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_to_litter               =>    pnf%m_deadstemn_to_litter                   , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_to_litter              =>    pnf%m_livecrootn_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_to_litter              =>    pnf%m_deadcrootn_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_retransn_to_litter                =>    pnf%m_retransn_to_litter                    , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_storage_to_litter           =>    pnf%m_leafn_storage_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_storage_to_litter          =>    pnf%m_frootn_storage_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_storage_to_litter       =>    pnf%m_livestemn_storage_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_storage_to_litter       =>    pnf%m_deadstemn_storage_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_storage_to_litter      =>    pnf%m_livecrootn_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_storage_to_litter      =>    pnf%m_deadcrootn_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_leafn_xfer_to_litter              =>    pnf%m_leafn_xfer_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   m_frootn_xfer_to_litter             =>    pnf%m_frootn_xfer_to_litter                 , & ! Input:  [real(r8) (:)]                                                    
   m_livestemn_xfer_to_litter          =>    pnf%m_livestemn_xfer_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemn_xfer_to_litter          =>    pnf%m_deadstemn_xfer_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootn_xfer_to_litter         =>    pnf%m_livecrootn_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootn_xfer_to_litter         =>    pnf%m_deadcrootn_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   gap_mortality_c_to_litr_met_c       =>    ccf%gap_mortality_c_to_litr_met_c           , & ! Input:  [real(r8) (:,:)]  C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
   gap_mortality_c_to_litr_cel_c       =>    ccf%gap_mortality_c_to_litr_cel_c           , & ! Input:  [real(r8) (:,:)]  C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
   gap_mortality_c_to_litr_lig_c       =>    ccf%gap_mortality_c_to_litr_lig_c           , & ! Input:  [real(r8) (:,:)]  C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
   gap_mortality_c_to_cwdc             =>    ccf%gap_mortality_c_to_cwdc                 , & ! Input:  [real(r8) (:,:)]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
   gap_mortality_n_to_litr_met_n       =>    cnf%gap_mortality_n_to_litr_met_n           , & ! Input:  [real(r8) (:,:)]  N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
   gap_mortality_n_to_litr_cel_n       =>    cnf%gap_mortality_n_to_litr_cel_n           , & ! Input:  [real(r8) (:,:)]  N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
   gap_mortality_n_to_litr_lig_n       =>    cnf%gap_mortality_n_to_litr_lig_n           , & ! Input:  [real(r8) (:,:)]  N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
   gap_mortality_n_to_cwdn             =>    cnf%gap_mortality_n_to_cwdn                 , & ! Input:  [real(r8) (:,:)]  N fluxes associated with gap mortality to CWD pool (gN/m3/s)
   leaf_prof                           =>    pps%leaf_prof                               , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of leaves                         
   froot_prof                          =>    pps%froot_prof                              , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of fine roots                     
   croot_prof                          =>    pps%croot_prof                              , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of coarse roots                   
   stem_prof                           =>    pps%stem_prof                                 & ! InOut:  [real(r8) (:,:)]  (1/m) profile of stems                          
   )


   do j = 1,nlevdecomp
      do pi = 1,maxpatch_pft
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (pi <=  npfts(c)) then
               p = pfti(c) + pi - 1
               
               if (pactive(p)) then
                  
                  
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
                       (m_leafc_storage_to_litter(p) + m_gresp_storage_to_litter(p))      * wtcol(p) * leaf_prof(p,j)
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
                  
                  
               end if
            end if
            
         end do
      end do
   end do

 end associate
end subroutine CNGapPftToColumn

end module CNGapMortalityMod
