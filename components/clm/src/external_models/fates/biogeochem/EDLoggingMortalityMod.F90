
module EDLoggingMortalityMod

   ! ====================================================================================
   !  Purpose: 1. create logging mortalities: 
   !           (a)logging mortality (cohort level)
   !           (b)collateral mortality (cohort level)
   !           (c)infrastructure mortality (cohort level)
   !           2. move the logged trunk fluxes from live into product pool 
   !           3. move logging-associated mortality fluxes from live to CWD
   !           4. keep carbon balance (in ed_total_balance_check)
   !
   !  Yi Xu
   !  Date: 2017
   ! ====================================================================================

   use FatesConstantsMod , only : r8 => fates_r8
   use EDTypesMod        , only : ed_cohort_type
   use EDTypesMod        , only : ed_patch_type
   use EDTypesMod        , only : ncwd
   use EDTypesMod        , only : ed_site_type
   use EDTypesMod        , only : ed_resources_management_type
   use EDTypesMod        , only : dtype_ilog
   use EDTypesMod        , only : dtype_ifall
   use EDTypesMod        , only : dtype_ifire
   use EDPftvarcon       , only : EDPftvarcon_inst
   use EDParamsMod       , only : logging_event_code
   use EDParamsMod       , only : logging_dbhmin
   use EDParamsMod       , only : logging_collateral_frac 
   use EDParamsMod       , only : logging_direct_frac
   use EDParamsMod       , only : logging_mechanical_frac 
   use EDParamsMod       , only : ED_val_understorey_death
   use FatesInterfaceMod , only : hlm_current_year
   use FatesInterfaceMod , only : hlm_current_month
   use FatesInterfaceMod , only : hlm_current_day
   use FatesInterfaceMod , only : hlm_model_day
   use FatesInterfaceMod , only : hlm_day_of_year 
   use FatesInterfaceMod , only : hlm_days_per_year
   use FatesInterfaceMod , only : hlm_use_logging
   use FatesConstantsMod , only : itrue,ifalse
   use FatesGlobals      , only : endrun => fates_endrun 
   use FatesGlobals      , only : fates_log
   use shr_log_mod       , only : errMsg => shr_log_errMsg

   implicit none
   private

   logical, protected :: logging_time   ! If true, logging should be 
                                        ! performed during the current time-step

   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   
   public :: LoggingMortality_frac
   public :: logging_litter_fluxes
   public :: logging_time
   public :: IsItLoggingTime

contains

   subroutine IsItLoggingTime(is_master,currentSite)

      ! -------------------------------------------------------------------------------
      ! This subroutine determines if the current dynamics step should enact
      ! the logging module.
      ! This is done by comparing the current model time to the logging event
      ! ids.  If there is a match, it is logging time.
      ! -------------------------------------------------------------------------------
     
      integer, intent(in) :: is_master
      type(ed_site_type), intent(inout), target :: currentSite     ! site structure

      integer :: icode   ! Integer equivalent of the event code (parameter file only allows reals)
      integer :: log_date  ! Day of month for logging exctracted from event code
      integer :: log_month ! Month of year for logging extraced from event code
      integer :: log_year  ! Year for logging extracted from event code
      character(len=64) :: fmt = '(a,i2.2,a,i2.2,a,i4.4)'

      logging_time = .false.
      icode = int(logging_event_code)

      if(hlm_use_logging.eq.ifalse) return

      if(icode .eq. 1) then
         ! Logging is turned off
         logging_time = .false.

      else if(icode .eq. 2) then
         ! Logging event on the first step
         if( hlm_model_day.eq.1 ) then
            logging_time = .true.
         end if

      else if(icode .eq. 3) then
         ! Logging event every day
         logging_time = .true.

      else if(icode .eq. 4) then
         ! logging event once a month
         if(hlm_current_day.eq.1  ) then
            logging_time = .true.
         end if

      else if(icode < 0 .and. icode > -366) then
         ! Logging event every year on specific day of year
         if(hlm_day_of_year .eq. icode  ) then
            logging_time = .true.
         end if

      else if(icode > 10000 ) then
         ! Specific Event: YYYYMMDD
         log_date  = icode - int(100* floor(real(icode)/100))
         log_year  = floor(real(icode)/10000)
         log_month = floor(real(icode)/100) - log_year*100

         if( hlm_current_day.eq.log_date    .and. &
               hlm_current_month.eq.log_month .and. &
               hlm_current_year.eq.log_year ) then
            logging_time = .true.
         end if
      else 
         ! Bad logging event flag
         write(fates_log(),*) 'An invalid logging code was specified in fates_params'
         write(fates_log(),*) 'Check EDLoggingMortalityMod.F90:IsItLoggingTime()'
         write(fates_log(),*) 'for a breakdown of the valid codes and change'
         write(fates_log(),*) 'fates_logging_event_code in the file accordingly.'
         write(fates_log(),*) 'exiting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      ! Initialize some site level diagnostics that are calculated for each event
      currentSite%resources_management%delta_litter_stock  = 0.0_r8
      currentSite%resources_management%delta_biomass_stock = 0.0_r8
      currentSite%resources_management%delta_individual    = 0.0_r8

      if(logging_time .and. (is_master.eq.itrue) ) then
         write(fates_log(),fmt) 'Logging Event Enacted on date: ', &
               hlm_current_month,'-',hlm_current_day,'-',hlm_current_year
      end if
      return
   end subroutine IsItLoggingTime

   ! ======================================================================================

   subroutine LoggingMortality_frac( pft_i, dbh, lmort_logging,lmort_collateral,lmort_infra )

      ! Arguments
      integer,  intent(in)  :: pft_i            ! pft index 
      real(r8), intent(in)  :: dbh              ! diameter at breast height (cm)
      real(r8), intent(out) :: lmort_logging    ! direct (harvestable) mortality fraction
      real(r8), intent(out) :: lmort_collateral ! collateral damage mortality fraction
      real(r8), intent(out) :: lmort_infra      ! infrastructure mortality fraction

      ! Parameters
      real(r8), parameter   :: adjustment = 1.0 ! adjustment for mortality rates

      if (logging_time) then 
         if(EDPftvarcon_inst%woody(pft_i) == 1)then ! only set logging rates for trees

            ! Pass logging rates to cohort level 

            if (dbh >= logging_dbhmin ) then
               lmort_logging = logging_direct_frac * adjustment
               lmort_collateral = logging_collateral_frac * adjustment
            else
               lmort_logging = 0.0_r8 
               lmort_collateral = 0.0_r8
            end if
           
            lmort_infra      = logging_mechanical_frac * adjustment
            !damage rates for size class < & > threshold_size need to be specified seperately

            ! Collateral damage to smaller plants below the direct logging size threshold
            ! will be applied via "understory_death" via the disturbance algorithm

         else
            lmort_logging    = 0.0_r8
            lmort_collateral = 0.0_r8
            lmort_infra      = 0.0_r8
         end if
      else 
         lmort_logging    = 0.0_r8
         lmort_collateral = 0.0_r8
         lmort_infra      = 0.0_r8
      end if

   end subroutine LoggingMortality_frac

   ! ============================================================================

   subroutine logging_litter_fluxes(currentSite, currentPatch, newPatch, patch_site_areadis)

      ! -------------------------------------------------------------------------------------------
      !
      !  DESCRIPTION:
      !  Carbon going from ongoing mortality into CWD pools. 
      !  This module includes only those fluxes associated with a disturbance generated by logging.
      !  Purpose: 
      !	  1) move logging-associated carbon to CWD and litter pool
      !   2) move the logging trunk from live into product pool 
      !   3) generate fluxes used in carbon balance checking
      !  E.g,:
      !  Remove trunk of logged trees from litter/CWD
      !  Add other parts of logged trees and all parts of collaterally and mechanically 
      !  damaged trees into CWD/litter  
      !
      !  This routine is only called if logging disturbance is the dominant disturbance.
      !
      !
      !  Note: The litter losses due to disturbance in the logging case is almost
      !        exactly like the natural tree-fall case.  The big differences are that
      !        the mortality rates governing the fluxes, follow a different rule set.
      !        We also compute an export flux (product) that does not go to litter.  
      !
      !  Trunk Product Flux: Only usable wood is exported from a site.  This is the above-ground
      !                      portion of the bole, and only boles associated with direct-logging,
      !                      not inftrastructure or collateral damage mortality.
      !        
      ! -------------------------------------------------------------------------------------------


      !USES:
      use SFParamsMod,  only : SF_val_cwd_frac
      use EDtypesMod,   only : area
      use EDtypesMod,   only : ed_site_type
      use EDtypesMod,   only : ed_patch_type
      use EDtypesMod,   only : ed_cohort_type
      use EDGrowthFunctionsMod, only : c_area

      ! !ARGUMENTS:
      type(ed_site_type)  , intent(inout), target  :: currentSite 
      type(ed_patch_type) , intent(inout), target  :: currentPatch
      type(ed_patch_type) , intent(inout), target  :: newPatch
      real(r8)            , intent(in)             :: patch_site_areadis

      !LOCAL VARIABLES:
      type(ed_cohort_type), pointer :: currentCohort
      real(r8) :: litter_area         ! area over which to distribute this litter (m2/site). 
      real(r8) :: np_mult             ! Fraction of the new patch which came from the current patch
      real(r8) :: direct_dead         ! Mortality count through direct logging
      real(r8) :: indirect_dead       ! Mortality count through: impacts, infrastructure and collateral damage
      real(r8) :: trunk_product_site  ! flux of carbon in trunk products exported off site      [ kgC/site ] 
                                      ! (note we are accumulating over the patch, but scale is site level)
      real(r8) :: delta_litter_stock  ! flux of carbon in total litter flux                     [ kgC/site ]
      real(r8) :: delta_biomass_stock ! total flux of carbon through mortality (litter+product) [ kgC/site ]
      real(r8) :: delta_individual    ! change in plant number through mortality [ plants/site ]
      real(r8) :: cwd_litter_density  ! Component woody biomass transferred through mortality [kgC/m2]
                                      ! (works with canopy_mortality_woody_litter, breaks into CWD partition
                                      !  and converts units to /m2)
      real(r8) :: woody_litter        ! Woody biomass transferred through mortality [kgC/site]
      real(r8) :: leaf_litter         ! Leafy biomass transferred through mortality [kgC/site]
      real(r8) :: root_litter         ! Rooty + storage biomass transferred through mort [kgC/site]
      real(r8) :: agb_frac            ! local copy of the above ground biomass fraction [fraction]
      integer  :: p                   ! pft index
      integer  :: c                   ! cwd index


      ! Zero some site level accumulator diagnsotics
      trunk_product_site  = 0.0_r8
      delta_litter_stock  = 0.0_r8
      delta_biomass_stock = 0.0_r8
      delta_individual    = 0.0_r8
      

      currentCohort => currentPatch%shortest
      do while(associated(currentCohort))       
         p = currentCohort%pft
         
        
         if(currentCohort%canopy_layer == 1)then         
            direct_dead   = currentCohort%n * currentCohort%lmort_logging
            indirect_dead = currentCohort%n * &
                  (currentCohort%lmort_collateral + currentCohort%lmort_infra)

         else
            if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then
               direct_dead   = 0.0_r8
               indirect_dead = ED_val_understorey_death * currentCohort%n * &
                     (patch_site_areadis/currentPatch%area)  !kgC/site/day
            else
               ! If the cohort of interest is grass, it will not experience
               ! any mortality associated with the logging disturbance
               direct_dead   = 0.0_r8
               indirect_dead = 0.0_r8
            end if
         end if

         agb_frac    = EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)
         litter_area = currentPatch%area 
         np_mult     = patch_site_areadis/newPatch%area
         
         ! ----------------------------------------------------------------------------------------
         ! Handle woody litter flux for non-bole components of biomass
         ! This litter is distributed between the current and new patches, &
         ! not to any other patches. This is really the eventually area of the current patch &
         ! (currentPatch%area-patch_site_areadis) +patch_site_areadis...
         ! For the new patch, only some fraction of its land area (patch_areadis/np%area) is 
         ! derived from the current patch, so we need to multiply by patch_areadis/np%area
         ! ----------------------------------------------------------------------------------------
         
         do c = 1,ncwd-1
            woody_litter = (direct_dead+indirect_dead) * &
                  (currentCohort%bdead+currentCohort%bsw)
            cwd_litter_density = SF_val_CWD_frac(c) * woody_litter / litter_area
            
            newPatch%cwd_ag(c)     = newPatch%cwd_ag(c)     + agb_frac * cwd_litter_density * np_mult
            currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + agb_frac * cwd_litter_density
            newPatch%cwd_bg(c)     = newPatch%cwd_bg(c)     + (1._r8-agb_frac) * cwd_litter_density * np_mult
            currentPatch%cwd_bg(c) = currentPatch%cwd_bg(c) + (1._r8-agb_frac) * cwd_litter_density 
            
            ! Diagnostics on fluxes into the AG and BG CWD pools
            currentSite%CWD_AG_diagnostic_input_carbonflux(c) =       &
                  currentSite%CWD_AG_diagnostic_input_carbonflux(c) + &
                  SF_val_CWD_frac(c) * woody_litter * hlm_days_per_year * agb_frac/ AREA 

            currentSite%CWD_BG_diagnostic_input_carbonflux(c) =       &
                  currentSite%CWD_BG_diagnostic_input_carbonflux(c) + &
                  SF_val_CWD_frac(c) * woody_litter * hlm_days_per_year * (1.0_r8 - agb_frac) / AREA

            ! Diagnostic specific to resource management code
            delta_litter_stock  = delta_litter_stock  + woody_litter * SF_val_CWD_frac(c)

         enddo
         
         ! ----------------------------------------------------------------------------------------
         ! Handle litter flux for the boles of infrastucture and collateral damage mort
         ! In this case the boles from direct logging are exported off-site and are not added 
         ! to the litter pools.  That is why we handle this outside the loop above. Only the 
         ! collateral damange and infrastructure logging is applied to bole litter
         ! ----------------------------------------------------------------------------------------

         woody_litter =  indirect_dead * (currentCohort%bdead+currentCohort%bsw)
         
         cwd_litter_density = SF_val_CWD_frac(ncwd) * woody_litter / litter_area
         
         newPatch%cwd_ag(ncwd)     = newPatch%cwd_ag(ncwd)     + agb_frac * cwd_litter_density * np_mult
         currentPatch%cwd_ag(ncwd) = currentPatch%cwd_ag(ncwd) + agb_frac * cwd_litter_density

         newPatch%cwd_bg(ncwd)     = newPatch%cwd_bg(ncwd)     + (1._r8-agb_frac) * cwd_litter_density * np_mult 
         currentPatch%cwd_bg(ncwd) = currentPatch%cwd_bg(ncwd) + (1._r8-agb_frac) * cwd_litter_density 
         
         currentSite%CWD_AG_diagnostic_input_carbonflux(ncwd) =       &
               currentSite%CWD_AG_diagnostic_input_carbonflux(ncwd) + &
               SF_val_CWD_frac(ncwd) * woody_litter * hlm_days_per_year * agb_frac/ AREA 

         currentSite%CWD_BG_diagnostic_input_carbonflux(ncwd) =       &
               currentSite%CWD_BG_diagnostic_input_carbonflux(ncwd) + &
               SF_val_CWD_frac(ncwd) * woody_litter * hlm_days_per_year * (1.0_r8 - agb_frac) / AREA

         delta_litter_stock  = delta_litter_stock  + woody_litter * SF_val_CWD_frac(ncwd)

         ! ----------------------------------------------------------------------------------------
         ! Handle litter flux for the belowground portion of directly logged boles
         ! ----------------------------------------------------------------------------------------

         woody_litter =  direct_dead * (currentCohort%bdead+currentCohort%bsw)
         cwd_litter_density = SF_val_CWD_frac(ncwd) * woody_litter / litter_area
         
         newPatch%cwd_bg(ncwd)     = newPatch%cwd_bg(ncwd)     + &
               (1._r8-agb_frac) * cwd_litter_density * np_mult 

         currentPatch%cwd_bg(ncwd) = currentPatch%cwd_bg(ncwd) + &
               (1._r8-agb_frac) * cwd_litter_density 
         
         currentSite%CWD_BG_diagnostic_input_carbonflux(ncwd) =       &
               currentSite%CWD_BG_diagnostic_input_carbonflux(ncwd) + &
               SF_val_CWD_frac(ncwd) * woody_litter * hlm_days_per_year * (1.0_r8 - agb_frac) / AREA
         
         
         ! ----------------------------------------------------------------------------------------
         ! Handle harvest (export, flux-out) flux for the above ground boles 
         ! In this case the boles from direct logging are exported off-site and are not added 
         ! to the litter pools.  That is why we handle this outside the loop above. Only the 
         ! collateral damange and infrastructure logging is applied to litter
         ! 
         ! Losses to the system as a whole, for C-balancing (kGC/site/day)
         ! Site level product, (kgC/site, accumulated over simulation)
         ! ----------------------------------------------------------------------------------------
         
         trunk_product_site = trunk_product_site + &
               SF_val_CWD_frac(ncwd) * agb_frac * direct_dead * (currentCohort%bdead+currentCohort%bsw)


         ! ----------------------------------------------------------------------------------------
         ! Handle fluxes of leaf, root and storage carbon into litter pools. 
         !  (none of these are exported)
         ! ----------------------------------------------------------------------------------------

         leaf_litter = (direct_dead+indirect_dead)*currentCohort%bl
         root_litter = (direct_dead+indirect_dead)*(currentCohort%br+currentCohort%bstore)

         newPatch%leaf_litter(p) = newPatch%leaf_litter(p) + leaf_litter / litter_area * np_mult
         newPatch%root_litter(p) = newPatch%root_litter(p) + root_litter / litter_area * np_mult 
         
         currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + leaf_litter / litter_area
         currentPatch%root_litter(p) = currentPatch%root_litter(p) + root_litter / litter_area
         
         ! track as diagnostic fluxes
         currentSite%leaf_litter_diagnostic_input_carbonflux(p) =       &
               currentSite%leaf_litter_diagnostic_input_carbonflux(p) + &
               leaf_litter * hlm_days_per_year / AREA
         
         currentSite%root_litter_diagnostic_input_carbonflux(p) =       &
               currentSite%root_litter_diagnostic_input_carbonflux(p) + &
               root_litter * hlm_days_per_year / AREA


         ! Logging specific diagnostics
         ! ----------------------------------------------------------------------------------------

         ! Note that litter stock also has terms above in the CWD loop
         delta_litter_stock  = delta_litter_stock  + &
                               leaf_litter         + &
                               root_litter

         delta_biomass_stock = delta_biomass_stock + &
                               leaf_litter         + &
                               root_litter         + &
                               (direct_dead+indirect_dead) * (currentCohort%bdead+currentCohort%bsw)

         delta_individual    = delta_individual    + &
                               direct_dead         + &
                               indirect_dead

         currentCohort => currentCohort%taller
      end do

      ! Update the amount of carbon exported from the site through logging
      ! operations.  Currently we assume only above-ground portion
      ! of the tree bole that experienced "direct" logging is exported
      ! This portion is known as "trunk_product_site


      currentSite%flux_out = currentSite%flux_out + trunk_product_site

      currentSite%resources_management%trunk_product_site  = &
           currentSite%resources_management%trunk_product_site + &
           trunk_product_site

      currentSite%resources_management%delta_litter_stock  = &
           currentSite%resources_management%delta_litter_stock + &
           delta_litter_stock

      currentSite%resources_management%delta_biomass_stock = &
           currentSite%resources_management%delta_biomass_stock + &
           delta_biomass_stock

      currentSite%resources_management%delta_individual    = &
           currentSite%resources_management%delta_individual + &
           delta_individual

      currentCohort => newPatch%shortest
      do while(associated(currentCohort))
         currentCohort%c_area = c_area(currentCohort)
         currentCohort => currentCohort%taller
      enddo

   end subroutine logging_litter_fluxes

end module EDLoggingMortalityMod
