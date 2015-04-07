module EDCohortDynamicsMod
  !
  ! !DESCRIPTION:
  ! Cohort stuctures in ED. 
  !
  ! !USES: 
  use shr_kind_mod          , only : r8 => shr_kind_r8;
  use clm_varctl            , only : iulog 
  use pftconMod             , only : pftcon
  use EDEcophysContype      , only : EDecophyscon
  use EDGrowthFunctionsMod  , only : c_area, tree_lai
  use EDTypesMod            , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod            , only : fusetol, nclmax
  use EDtypesMod            , only : ncwd, numcohortsperpatch, udata
  !
  implicit none
  private
  !
  public :: create_cohort
  public :: zero_cohort
  public :: nan_cohort
  public :: terminate_cohorts
  public :: fuse_cohorts
  public :: insert_cohort
  public :: sort_cohorts
  public :: copy_cohort
  public :: count_cohorts
  public :: countCohorts
  public :: allocate_live_biomass

  ! 10/30/09: Created by Rosie Fisher
  !-------------------------------------------------------------------------------------!

contains

  !-------------------------------------------------------------------------------------!
  subroutine create_cohort(patchptr, pft, nn, hite, dbh, &
       balive, bdead, bstore, laimemory, status, ctrim, clayer)
    !
    ! !DESCRIPTION:
    ! create new cohort
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout), pointer :: patchptr
    integer,  intent(in)   :: pft       ! Cohort Plant Functional Type
    integer,  intent(in)   :: clayer    ! canopy status of cohort (1 = canopy, 2 = understorey, etc.)
    integer,  intent(in)   :: status    ! growth status of plant  (2 = leaves on , 1 = leaves off)
    real(r8), intent(in)   :: nn        ! number of individuals in cohort per 'area' (10000m2 default)
    real(r8), intent(in)   :: hite      ! height: meters
    real(r8), intent(in)   :: dbh       ! dbh: cm
    real(r8), intent(in)   :: balive    ! total living biomass: kGC per indiv
    real(r8), intent(in)   :: bdead     ! total dead biomass: kGC per indiv
    real(r8), intent(in)   :: bstore    ! stored carbon: kGC per indiv
    real(r8), intent(in)   :: laimemory ! target leaf biomass- set from previous year: kGC per indiv
    real(r8), intent(in)   :: ctrim     ! What is the fraction of the maximum leaf biomass that we are targeting? :-
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: new_cohort         ! Pointer to New Cohort structure.
    type(ed_cohort_type), pointer :: storesmallcohort 
    type(ed_cohort_type), pointer :: storebigcohort   
    integer :: tnull,snull                      ! are the tallest and shortest cohorts allocate
    !----------------------------------------------------------------------
 
    allocate(new_cohort)
    udata%cohort_number = udata%cohort_number + 1 !give each cohort a unique number for checking cohort fusing routine.
    
    call nan_cohort(new_cohort)  ! Make everything in the cohort not-a-number
    call zero_cohort(new_cohort) ! Zero things that need to be zeroed. 

    !**********************/
    ! Define cohort state variable
    !**********************/
 
    new_cohort%indexnumber  = udata%cohort_number
    new_cohort%siteptr      => patchptr%siteptr
    new_cohort%patchptr     => patchptr
    new_cohort%pft          = pft     
    new_cohort%status_coh   = status
    new_cohort%n            = nn
    new_cohort%hite         = hite
    new_cohort%dbh          = dbh
    new_cohort%canopy_trim  = ctrim
    new_cohort%canopy_layer = clayer
    new_cohort%laimemory    = laimemory
    new_cohort%bdead        = bdead
    new_cohort%balive       = balive
    new_cohort%bstore       = bstore

    if (new_cohort%dbh <= 0_r8 .or. new_cohort%n == 0._r8 .or. new_cohort%pft == 0 &
          .or. new_cohort%canopy_trim <= 0_r8 .or. new_cohort%balive <= 0._r8) then
       write(iulog,*) 'ED: something is zero in create_cohort',new_cohort%indexnumber,new_cohort%dbh,new_cohort%n, &
       new_cohort%pft,new_cohort%canopy_trim,new_cohort%balive
    endif
    if (new_cohort%siteptr%status==2.and.pftcon%season_decid(pft) == 1) then
      new_cohort%laimemory = 0.0_r8
    endif
    if (new_cohort%siteptr%dstatus==2.and.pftcon%stress_decid(pft) == 1) then
      new_cohort%laimemory = 0.0_r8
    endif
    
    ! Calculate live biomass allocation
    call allocate_live_biomass(new_cohort)

    ! Assign canopy extent and depth
    new_cohort%c_area  = c_area(new_cohort)
    new_cohort%treelai = tree_lai(new_cohort)
    new_cohort%lai     = new_cohort%treelai * new_cohort%c_area/patchptr%area
    new_cohort%treesai = 0.0_r8 !FIX(RF,032414)   

    ! Put cohort at the right place in the linked list
    storebigcohort   => patchptr%tallest
    storesmallcohort => patchptr%shortest 

    if (associated(patchptr%tallest)) then
       tnull = 0
    else
       tnull = 1
       patchptr%tallest => new_cohort
    endif

    if (associated(patchptr%shortest)) then
       snull = 0
    else
       snull = 1
       patchptr%shortest => new_cohort 
    endif

    call insert_cohort(new_cohort, patchptr%tallest, patchptr%shortest, tnull, snull, &
         storebigcohort, storesmallcohort)

    patchptr%tallest  => storebigcohort 
    patchptr%shortest => storesmallcohort

  end subroutine create_cohort

  !-------------------------------------------------------------------------------------!
  subroutine allocate_live_biomass(cc_p)
    !
    ! !DESCRIPTION:
    ! Divide alive biomass between leaf, root and sapwood parts. 
    ! Needs to be called whenver balive changes. 
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type (ed_cohort_type), intent(inout), target  :: cc_p ! current cohort pointer
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type), pointer :: currentCohort
    real(r8)  :: leaf_frac    ! fraction of live biomass in leaves
    real(r8)  :: ideal_balive ! theoretical ideal (root and stem) biomass for deciduous trees with leaves off. 
                              ! accounts for the fact that live biomass may decline in the off-season, 
                              ! making leaf_memory unrealistic.    
    real(r8)  :: ratio_balive ! ratio between root+shoot biomass now and root+shoot biomass when leaves fell off. 
                                            
    integer   :: ft           ! functional type
    integer   :: leaves_off_switch
    !----------------------------------------------------------------------

    currentCohort => cc_p
    ft = currentcohort%pft
    leaf_frac = 1.0_r8/(1.0_r8 + EDecophyscon%sapwood_ratio(ft) * currentcohort%hite + pftcon%froot_leaf(ft))     

    currentcohort%bl = currentcohort%balive*leaf_frac    
    ratio_balive = 1.0_r8
    !for deciduous trees, there are no leaves  
    
    if (pftcon%evergreen(ft) == 1) then
       currentcohort%laimemory = 0._r8
       currentcohort%status_coh    = 2      
    endif    

    !diagnore the root and stem biomass from the functional balance hypothesis. This is used when the leaves are 
    !fully on. 
    currentcohort%br  = pftcon%froot_leaf(ft) * (currentcohort%balive + currentcohort%laimemory) * leaf_frac
    currentcohort%bsw = EDecophyscon%sapwood_ratio(ft) * currentcohort%hite *(currentcohort%balive + &
         currentcohort%laimemory)*leaf_frac 
    
    leaves_off_switch = 0
    if (currentcohort%status_coh == 1.and.pftcon%stress_decid(ft) == 1.and.currentcohort%siteptr%dstatus==1) then !no leaves 
      leaves_off_switch = 1 !drought decid
    endif
    if (currentcohort%status_coh == 1.and.pftcon%season_decid(ft) == 1.and.currentcohort%siteptr%status==1) then !no leaves
      leaves_off_switch = 1 !cold decid
    endif
  
    if (leaves_off_switch==1) then

    !the purpose of this section is to figure out the root and stem biomass when the leaves are off
    !at this point, we know the former leaf mass (laimemory) and the current alive mass
    !because balive may decline in the off-season, we need to adjust the root and stem biomass that are predicted
    !from the laimemory, for the fact that we now might not have enough live biomass to support the hypothesized root mass
    !thus, we use 'ratio_balive' to adjust br and bsw. Apologies that this is so complicated! RF
       currentcohort%bl  = 0.0_r8
       ideal_balive      = currentcohort%laimemory * pftcon%froot_leaf(ft) +  &
            currentcohort%laimemory*  EDecophyscon%sapwood_ratio(ft) * currentcohort%hite
       currentcohort%br  = pftcon%froot_leaf(ft) * (ideal_balive + currentcohort%laimemory) * leaf_frac
       currentcohort%bsw = EDecophyscon%sapwood_ratio(ft) * currentcohort%hite *(ideal_balive + &
            currentcohort%laimemory)*leaf_frac 
     
       ratio_balive           = currentcohort%balive / ideal_balive
       currentcohort%br       = currentcohort%br  * ratio_balive
       currentcohort%bsw      = currentcohort%bsw * ratio_balive               
    endif

    
    if (abs(currentcohort%balive -currentcohort%bl- currentcohort%br - currentcohort%bsw)>1e-12) then
       write(iulog,*) 'issue with carbon allocation in create_cohort',&
       currentcohort%balive -currentcohort%bl- currentcohort%br - currentcohort%bsw, currentcohort%status_coh,currentcohort%balive
       write(iulog,*) 'actual vs predicted balive',ideal_balive,currentcohort%balive ,ratio_balive,leaf_frac
       write(iulog,*) 'leaf,root,stem',currentcohort%bl,currentcohort%br,currentcohort%bsw
    endif
    currentCohort%b   = currentCohort%bdead + currentCohort%balive

  end subroutine allocate_live_biomass

  !-------------------------------------------------------------------------------------!
  subroutine nan_cohort(cc_p)
    !
    ! !DESCRIPTION:
    !  Make all the cohort variables NaN so they aren't used before defined.   
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)  
    !
    ! !ARGUMENTS    
    type (ed_cohort_type), intent(inout), target  :: cc_p
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type)   , pointer             :: currentCohort
    !----------------------------------------------------------------------

    currentCohort => cc_p

    currentCohort%taller      => null()       ! pointer to next tallest cohort     
    currentCohort%shorter     => null()       ! pointer to next shorter cohort     
    currentCohort%patchptr    => null()       ! pointer to patch that cohort is in
    currentCohort%siteptr     => null()       ! pointer to site that cohort is in

    nullify(currentCohort%taller) 
    nullify(currentCohort%shorter) 
    nullify(currentCohort%patchptr) 
    nullify(currentCohort%siteptr) 

    ! VEGETATION STRUCTURE
    currentCohort%pft                = 999 ! pft number                           
    currentCohort%indexnumber        = 999 ! unique number for each cohort. (within clump?)
    currentCohort%canopy_layer       = 999 ! canopy status of cohort (1 = canopy, 2 = understorey, etc.)   
    currentCohort%NV                 = 999 ! Number of leaf layers: -
    currentCohort%status_coh         = 999 ! growth status of plant  (2 = leaves on , 1 = leaves off)

    currentCohort%n                  = nan ! number of individuals in cohort per 'area' (10000m2 default)     
    currentCohort%dbh                = nan ! 'diameter at breast height' in cm                            
    currentCohort%hite               = nan ! height: meters                   
    currentCohort%balive             = nan ! total living biomass: kGC per indiv  
    currentCohort%bdead              = nan ! dead biomass:  kGC per indiv          
    currentCohort%bstore             = nan ! stored carbon: kGC per indiv
    currentCohort%laimemory          = nan ! target leaf biomass- set from previous year: kGC per indiv
    currentCohort%b                  = nan ! total biomass: kGC per indiv        
    currentCohort%bsw                = nan ! sapwood in stem and roots: kGC per indiv
    currentCohort%bl                 = nan ! leaf biomass: kGC per indiv       
    currentCohort%br                 = nan ! fine root biomass: kGC per indiv
    currentCohort%lai                = nan ! leaf area index of cohort   m2/m2      
    currentCohort%sai                = nan ! stem area index of cohort   m2/m2
    currentCohort%gscan              = nan ! Stomatal resistance of cohort. 
    currentCohort%canopy_trim        = nan ! What is the fraction of the maximum leaf biomass that we are targeting? :-
    currentCohort%leaf_cost          = nan ! How much does it cost to maintain leaves: kgC/m2/year-1
    currentCohort%excl_weight        = nan ! How much of this cohort is demoted each year, as a proportion of all cohorts:-
    currentCohort%prom_weight        = nan ! How much of this cohort is promoted each year, as a proportion of all cohorts:-
    currentCohort%c_area             = nan ! areal extent of canopy (m2)
    currentCohort%treelai            = nan ! lai of tree (total leaf area (m2) / canopy area (m2)
    currentCohort%treesai            = nan ! stem area index of tree (total stem area (m2) / canopy area (m2)

    ! CARBON FLUXES 
    currentCohort%gpp                = nan ! GPP:  kgC/indiv/year
    currentCohort%gpp_clm            = nan ! GPP:  kgC/indiv/timestep
    currentCohort%gpp_acc            = nan ! GPP:  kgC/indiv/day         
    currentCohort%npp                = nan ! NPP:  kgC/indiv/year
    currentCohort%npp_clm            = nan ! NPP:  kGC/indiv/timestep
    currentCohort%npp_acc            = nan ! NPP:  kgC/indiv/day  
    currentCohort%year_net_uptake(:) = nan ! Net uptake of individual leaf layers kgC/m2/year
    currentCohort%ts_net_uptake(:)   = nan ! Net uptake of individual leaf layers kgC/m2/s
    currentCohort%resp               = nan ! RESP: kgC/indiv/year
    currentCohort%resp_clm           = nan ! RESP: kgC/indiv/timestep
    currentCohort%resp_acc           = nan ! RESP: kGC/cohort/day

    !RESPIRATION
    currentCohort%rd                 = nan
    currentCohort%resp_m             = nan ! Maintenance respiration.  kGC/cohort/year
    currentCohort%resp_g             = nan ! Growth respiration.       kGC/cohort/year
    currentCohort%livestem_mr        = nan ! Live stem maintenance respiration. kgC/indiv/s-1 
    currentCohort%livecroot_mr       = nan ! Coarse root maintenance respiration. kgC/indiv/s-1 
    currentCohort%froot_mr           = nan ! Fine root maintenance respiration. kgC/indiv/s-1 

    ! ALLOCATION
    currentCohort%md                 = nan ! plant maintenance demand: kgC/indiv/year
    currentCohort%leaf_md            = nan ! leaf  maintenance demand: kgC/indiv/year
    currentCohort%root_md            = nan ! root  maintenance demand: kgC/indiv/year
    currentCohort%carbon_balance     = nan ! carbon remaining for growth and storage: kg/indiv/year
    currentCohort%dmort              = nan ! proportional mortality rate. (year-1)
    currentCohort%seed_prod          = nan ! reproduction seed and clonal: KgC/indiv/year
    currentCohort%c_area             = nan ! areal extent of canopy (m2)
    currentCohort%treelai            = nan ! lai of tree (total leaf area (m2) / canopy area (m2)
    currentCohort%treesai            = nan ! stem area index of tree (total stem area (m2) / canopy area (m2)
    currentCohort%leaf_litter        = nan ! leaf litter from phenology: KgC/m2
    currentCohort%woody_turnover     = nan ! amount of wood lost each day: kgC/indiv/year. Currently set to zero.

    ! NITROGEN POOLS      
    currentCohort%livestemn          = nan ! live stem nitrogen       : KgN/invid
    currentCohort%livecrootn         = nan ! live coarse root nitrogen: KgN/invid
    currentCohort%frootn             = nan ! fine root  nitrogen      : KgN/invid

    ! VARIABLES NEEDED FOR INTEGRATION 
    currentCohort%dndt               = nan ! time derivative of cohort size 
    currentCohort%dhdt               = nan ! time derivative of height 
    currentCohort%ddbhdt             = nan ! time derivative of dbh 
    currentCohort%dbalivedt          = nan ! time derivative of total living biomass 
    currentCohort%dbdeaddt           = nan ! time derivative of dead biomass 
    currentCohort%dbstoredt          = nan ! time derivative of stored biomass 
    currentCohort%storage_flux       = nan ! flux from npp into bstore

    ! FIRE
    currentCohort%cfa                = nan ! proportion of crown affected by fire
    currentCohort%cambial_mort       = nan ! probability that trees dies due to cambial char P&R (1986)
    currentCohort%crownfire_mort     = nan ! probability of tree post-fire mortality due to crown scorch
    currentCohort%fire_mort          = nan ! post-fire mortality from cambial and crown damage assuming two are independent

  end subroutine nan_cohort

  !-------------------------------------------------------------------------------------!
  subroutine zero_cohort(cc_p)
    !
    ! !DESCRIPTION:
    ! Zero variables that need to be accounted for if 
    ! this cohort is altered before they are defined.       
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type (ed_cohort_type), intent(inout), target  :: cc_p
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type)   , pointer             :: currentCohort
    !----------------------------------------------------------------------

    currentCohort => cc_p

    currentCohort%NV                 = 0    
    currentCohort%status_coh         = 0    
    currentCohort%rd                 = 0._r8
    currentCohort%resp_m             = 0._r8 
    currentCohort%resp_g             = 0._r8
    currentCohort%livestem_mr        = 0._r8
    currentCohort%livecroot_mr       = 0._r8
    currentCohort%froot_mr           = 0._r8
    currentCohort%fire_mort          = 0._r8 
    currentcohort%npp_acc            = 0._r8
    currentcohort%gpp_acc            = 0._r8
    currentcohort%resp_acc           = 0._r8
    currentcohort%npp_clm            = 0._r8
    currentcohort%gpp_clm            = 0._r8
    currentcohort%resp_clm           = 0._r8
    currentcohort%resp               = 0._r8
    currentcohort%carbon_balance     = 0._r8
    currentcohort%leaf_litter        = 0._r8
    currentcohort%year_net_uptake(:) = 999 ! this needs to be 999, or trimming of new cohorts will break. 
    currentcohort%ts_net_uptake(:)   = 0._r8
    currentcohort%seed_prod          = 0._r8
    currentcohort%cfa                = 0._r8 
    currentcohort%md                 = 0._r8
    currentcohort%root_md            = 0._r8
    currentcohort%leaf_md            = 0._r8
    currentcohort%npp                = 0._r8 
    currentcohort%gpp                = 0._r8  
    currentcohort%storage_flux       = 0._r8   
    currentcohort%dmort              = 0._r8 
    currentcohort%gscan              = 0._r8 
    currentcohort%treesai            = 0._r8  

  end subroutine zero_cohort

  !-------------------------------------------------------------------------------------!
  subroutine terminate_cohorts( patchptr )
    !
    ! !DESCRIPTION:
    ! terminates cohorts when they get too small      
    !
    ! !USES:
    use EDParamsMod, only : ED_val_ag_biomass
    use SFParamsMod, only : SF_val_CWD_frac
    !
    ! !ARGUMENTS    
    type (ed_patch_type), intent(inout), target :: patchptr
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type)  , pointer :: currentPatch
    type (ed_cohort_type) , pointer :: currentCohort
    type (ed_cohort_type) , pointer :: nextc
    integer :: terminate   ! do we terminate (1) or not (0) 
    integer :: c           ! counter for litter size class. 
    !----------------------------------------------------------------------

    currentPatch  => patchptr
    currentCohort => currentPatch%tallest  

    do while (associated(currentCohort))
       nextc      => currentCohort%shorter    
       terminate = 0 

       ! Not enough n or dbh
       if (currentCohort%n/currentPatch%area <= 0.00001_r8 .or. currentCohort%dbh <  &
            0.00001_r8.and.currentCohort%bstore < 0._r8) then
          terminate = 1
          !  write(iulog,*) 'terminating cohorts 1',currentCohort%n/currentPatch%area,currentCohort%dbh
       endif

       ! In the third canopy layer
       if (currentCohort%canopy_layer > NCLMAX) then 
          terminate = 1
          ! write(iulog,*) 'terminating cohorts 2', currentCohort%canopy_layer
       endif

       ! live biomass pools are terminally depleted
       if (currentCohort%balive < 1e-10_r8 .or. currentCohort%bstore < 1e-10_r8) then 
          terminate = 1  
          ! write(iulog,*) 'terminating cohorts 3', currentCohort%balive,currentCohort%bstore
       endif

       ! Total cohort biomass is negative
       if (currentCohort%balive+currentCohort%bdead+currentCohort%bstore < 0._r8) then
          terminate = 1
          ! write(iulog,*) 'terminating cohorts 4', currentCohort%balive, currentCohort%bstore, currentCohort%bdead, &
          ! currentCohort%balive+currentCohort%bdead+&
          ! currentCohort%bstore, currentCohort%n
       endif


       if (terminate == 1) then 
          if (.not. associated(currentCohort%taller)) then
             currentPatch%tallest => currentCohort%shorter
          else 
             currentCohort%taller%shorter => currentCohort%shorter
          endif
          if (.not. associated(currentCohort%shorter)) then
             currentPatch%shortest => currentCohort%taller
          else 
             currentCohort%shorter%taller => currentCohort%taller
          endif

          !put the litter from the terminated cohorts straight into the fragmenting pools
          if (currentCohort%n.gt.0.0_r8) then
             do c=1,ncwd

                currentPatch%CWD_AG(c)  = currentPatch%CWD_AG(c) + currentCohort%n*(currentCohort%bdead+currentCohort%bsw) / &
                     currentPatch%area &
                     * SF_val_CWD_frac(c) * ED_val_ag_biomass 
                currentPatch%CWD_BG(c)  = currentPatch%CWD_BG(c) + currentCohort%n*(currentCohort%bdead+currentCohort%bsw) / &
                     currentPatch%area &
                     * SF_val_CWD_frac(c) * (1.0_r8 -  ED_val_ag_biomass) 
             enddo

             currentPatch%leaf_litter(currentCohort%pft) = currentPatch%leaf_litter(currentCohort%pft) + currentCohort%n* &
                  (currentCohort%bl)/currentPatch%area
             currentPatch%root_litter(currentCohort%pft) = currentPatch%root_litter(currentCohort%pft) + currentCohort%n* &
                  (currentCohort%br+currentCohort%bstore)/currentPatch%area 

             deallocate(currentCohort)     
          endif
       endif
       currentCohort => nextc
    enddo

  end subroutine terminate_cohorts

  !-------------------------------------------------------------------------------------!
  subroutine fuse_cohorts(patchptr)  
    !
    ! !DESCRIPTION:
    ! Join similar cohorts to reduce total number            
    !
    ! !USES:
    use clm_varpar  , only :  nlevcan_ed
    !
    ! !ARGUMENTS    
    type (ed_patch_type), intent(inout), target :: patchptr
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type)  , pointer :: currentPatch
    type (ed_cohort_type) , pointer :: currentCohort, nextc, nextnextc
    integer  :: i  
    integer  :: fusion_took_place
    integer  :: maxcohorts !maximum total no of cohorts. Needs to be >numpft_edx2 
    integer  :: iterate    !do we need to keep fusing to get below maxcohorts?
    integer  :: nocohorts
    real(r8) :: newn
    real(r8) :: diff
    real(r8) :: dynamic_fusion_tolerance
    !----------------------------------------------------------------------
    
    !set initial fusion tolerance
    dynamic_fusion_tolerance = fusetol
   
    !This needs to be a function of the canopy layer, because otherwise, at canopy closure
    !the number of cohorts doubles and very dissimilar cohorts are fused together
    !because c_area and biomass are non-linear with dbh, this causes several mass inconsistancies
    !in theory, all of this routine therefore causes minor losses of C and area, but these are below 
    !detection limit normally. 
    iterate = 1
    fusion_took_place = 0   
    currentPatch => patchptr
    maxcohorts =  currentPatch%NCL_p * numCohortsPerPatch  
    !---------------------------------------------------------------------!
    !  Keep doing this until nocohorts <= maxcohorts                         !
    !---------------------------------------------------------------------!
    if (associated(currentPatch%shortest)) then  
       do while(iterate == 1)

         currentCohort => currentPatch%tallest

         !CHANGED FROM C VERSION  loop from tallest to smallest, fusing if they are similar
         do while (currentCohort%indexnumber /= currentPatch%shortest%indexnumber)  
          nextc => currentPatch%tallest

          do while (associated(nextc))
             nextnextc => nextc%shorter                      
             diff = abs((currentCohort%dbh - nextc%dbh)/(0.5*(currentCohort%dbh + nextc%dbh)))  

             !Criteria used to divide up the height continuum into different cohorts.

             if (diff < dynamic_fusion_tolerance) then

                if (currentCohort%indexnumber /= nextc%indexnumber) then

                   if (currentCohort%pft == nextc%pft) then              

                      ! check cohorts in same c. layer. before fusing
                      if (currentCohort%canopy_layer == nextc%canopy_layer) then 
                         fusion_took_place = 1         
                         newn = currentCohort%n + nextc%n    ! sum individuals in both cohorts.     

                         currentCohort%balive    = (currentCohort%n*currentCohort%balive    + nextc%n*nextc%balive)/newn
                         currentCohort%bdead     = (currentCohort%n*currentCohort%bdead     + nextc%n*nextc%bdead)/newn
                         currentCohort%bstore    = (currentCohort%n*currentCohort%bstore    + nextc%n*nextc%bstore)/newn   
                         currentCohort%seed_prod = (currentCohort%n*currentCohort%seed_prod + nextc%n*nextc%seed_prod)/newn  
                         currentCohort%root_md   = (currentCohort%n*currentCohort%root_md   + nextc%n*nextc%root_md)/newn   
                         currentCohort%leaf_md   = (currentCohort%n*currentCohort%leaf_md   + nextc%n*nextc%leaf_md)/newn  
                         currentCohort%laimemory = (currentCohort%n*currentCohort%laimemory + nextc%n*nextc%laimemory)/newn
                         currentCohort%md        = (currentCohort%n*currentCohort%md        + nextc%n*nextc%md)/newn

                         currentCohort%carbon_balance = (currentCohort%n*currentCohort%carbon_balance + &
                              nextc%n*nextc%carbon_balance)/newn
                         currentCohort%storage_flux = (currentCohort%n*currentCohort%storage_flux + &
                              nextc%n*nextc%storage_flux)/newn               

                         currentCohort%b           = (currentCohort%n*currentCohort%b           + nextc%n*nextc%b)/newn
                         currentCohort%bsw         = (currentCohort%n*currentCohort%bsw         + nextc%n*nextc%bsw)/newn
                         currentCohort%bl          = (currentCohort%n*currentCohort%bl          + nextc%n*nextc%bl)/newn
                         currentCohort%br          = (currentCohort%n*currentCohort%br          + nextc%n*nextc%br)/newn
                         currentCohort%hite        = (currentCohort%n*currentCohort%hite        + nextc%n*nextc%hite)/newn         
                         currentCohort%dbh         = (currentCohort%n*currentCohort%dbh         + nextc%n*nextc%dbh)/newn
                         currentCohort%gpp_acc     = (currentCohort%n*currentCohort%gpp_acc     + nextc%n*nextc%gpp_acc)/newn
                         currentCohort%npp_acc     = (currentCohort%n*currentCohort%npp_acc     + nextc%n*nextc%npp_acc)/newn
                         currentCohort%resp_acc    = (currentCohort%n*currentCohort%resp_acc    + nextc%n*nextc%resp_acc)/newn
                         currentCohort%resp        = (currentCohort%n*currentCohort%resp        + nextc%n*nextc%resp)/newn
                         currentCohort%npp         = (currentCohort%n*currentCohort%npp         + nextc%n*nextc%npp)/newn
                         currentCohort%gpp         = (currentCohort%n*currentCohort%gpp         + nextc%n*nextc%gpp)/newn
                         currentCohort%canopy_trim = (currentCohort%n*currentCohort%canopy_trim + nextc%n*nextc%canopy_trim)/newn
                         currentCohort%dmort       = (currentCohort%n*currentCohort%dmort       + nextc%n*nextc%dmort)/newn
                         currentCohort%fire_mort   = (currentCohort%n*currentCohort%fire_mort   + nextc%n*nextc%fire_mort)/newn
                         currentCohort%leaf_litter = (currentCohort%n*currentCohort%leaf_litter + nextc%n*nextc%leaf_litter)/newn

                         do i=1, nlevcan_ed     
                            if (currentCohort%year_net_uptake(i) == 999._r8 .or. nextc%year_net_uptake(i) == 999._r8) then
                               currentCohort%year_net_uptake(i) = min(nextc%year_net_uptake(i),currentCohort%year_net_uptake(i))
                            else
                               currentCohort%year_net_uptake(i) = (currentCohort%n*currentCohort%year_net_uptake(i) + &
                                  nextc%n*nextc%year_net_uptake(i))/newn                
                            endif
                         enddo

                         currentCohort%n = newn     
                         !remove fused cohort from the list
                         nextc%taller%shorter => nextnextc        
                         if (.not. associated(nextc%shorter)) then !this is the shortest cohort. 
                            currentPatch%shortest => nextc%taller
                         else
                            nextnextc%taller => nextc%taller
                         endif
                         if (associated(nextc)) then       
                            deallocate(nextc)            
                         endif
                      endif !canopy layer
                   endif !pft
                endif  !index no. 
             endif !diff   

             if (associated(nextc)) then             
                nextc => nextc%shorter    
             else
                nextc => nextnextc !if we have removed next
             endif
          enddo !end checking nextc cohort loop

          if (associated (currentCohort%shorter)) then
             currentCohort => currentCohort%shorter
          endif
       enddo !end currentCohort cohort loop

       !---------------------------------------------------------------------!
       ! Is the number of cohorts larger than the maximum?                   !
       !---------------------------------------------------------------------!   
       nocohorts = 0
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))
          nocohorts = nocohorts + 1
          currentCohort => currentCohort%shorter
       enddo

       if (nocohorts > maxcohorts) then
          iterate = 1
          dynamic_fusion_tolerance = dynamic_fusion_tolerance * 1.1_r8
          !write(iulog,*) 'maxcohorts exceeded',dynamic_fusion_tolerance
          !---------------------------------------------------------------------!
          ! Making profile tolerance larger means that more fusion will happen  !
          !---------------------------------------------------------------------!        
       else
          iterate = 0
       endif

    enddo !do while nocohorts>maxcohorts

    endif ! patch. 

    if (fusion_took_place == 1) then  ! if fusion(s) occured sort cohorts 
       call sort_cohorts(currentPatch)
    endif

  end subroutine fuse_cohorts

!-------------------------------------------------------------------------------------!

  subroutine sort_cohorts(patchptr)  
    ! ============================================================================
    !                 sort cohorts into the correct order   DO NOT CHANGE THIS IT WILL BREAK
    ! ============================================================================

    type(ed_patch_type) , intent(inout), target :: patchptr

    type(ed_patch_type) , pointer :: current_patch
    type(ed_cohort_type), pointer :: current_c, next_c
    type(ed_cohort_type), pointer :: shortestc, tallestc 
    type(ed_cohort_type), pointer :: storesmallcohort 
    type(ed_cohort_type), pointer :: storebigcohort   
    integer :: snull,tnull

    current_patch => patchptr
    tallestc  => NULL()
    shortestc => NULL()
    storebigcohort   => null()
    storesmallcohort => null()
    current_c => current_patch%tallest 

    do while (associated(current_c))  
       next_c => current_c%shorter
       tallestc  => storebigcohort 
       shortestc => storesmallcohort   
       if (associated(tallestc)) then
          tnull = 0
       else
          tnull = 1
          tallestc => current_c
       endif

       if (associated(shortestc)) then
          snull = 0
       else
          snull = 1
          shortestc => current_c
       endif

       call insert_cohort(current_c, tallestc, shortestc, tnull, snull, storebigcohort, storesmallcohort)

       current_patch%tallest  => storebigcohort 
       current_patch%shortest => storesmallcohort
       current_c => next_c

    enddo

  end subroutine sort_cohorts

  !-------------------------------------------------------------------------------------!
  subroutine insert_cohort(pcc, ptall, pshort, tnull, snull, storebigcohort, storesmallcohort)
    !
    ! !DESCRIPTION:
    ! Insert cohort into linked list                  
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_cohort_type) , intent(inout), target          :: pcc
    type(ed_cohort_type) , intent(inout), target          :: ptall
    type(ed_cohort_type) , intent(inout), target          :: pshort
    integer              , intent(in)                     :: tnull
    integer              , intent(in)                     :: snull
    type(ed_cohort_type) , intent(inout),pointer,optional :: storesmallcohort ! storage of the smallest cohort for insertion routine
    type(ed_cohort_type) , intent(inout),pointer,optional :: storebigcohort   ! storage of the largest cohort for insertion routine 
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type),  pointer :: currentPatch
    type(ed_cohort_type), pointer :: current
    type(ed_cohort_type), pointer :: tallptr, shortptr, icohort
    type(ed_cohort_type), pointer :: ptallest, pshortest 
    real(r8) :: tsp
    integer :: tallptrnull,exitloop
    !----------------------------------------------------------------------

    currentPatch => pcc%patchptr
    ptallest => ptall
    pshortest => pshort

    if (tnull == 1) then
       ptallest => null()
    endif
    if (snull == 1) then
       pshortest => null()
    endif

    icohort => pcc ! assign address to icohort local name  
    !place in the correct place in the linked list of heights 
    !begin by finding cohort that is just taller than the new cohort 
    tsp = icohort%dbh

    current => pshortest
    exitloop = 0
    !starting with shortest tree on the grid, find tree just  
    !taller than tree being considered and return its pointer 
    if (associated(current)) then
       do while (associated(current).and.exitloop == 0)
          if (current%dbh < tsp) then
             current => current%taller   
          else
             exitloop = 1 
          endif
       enddo
    endif

    if (associated(current)) then
       tallptr => current
       tallptrnull = 0
    else
       tallptr => null()
       tallptrnull = 1
    endif

    !new cohort is tallest 
    if (.not.associated(tallptr)) then  
       !new shorter cohort to the new cohort is the old tallest cohort 
       shortptr => ptallest

       !new cohort is tallest cohort and next taller remains null 
       ptallest => icohort
       if (present(storebigcohort)) then
          storebigcohort => icohort
       end if
       currentPatch%tallest => icohort 
       icohort%patchptr%tallest => icohort  
       !new cohort is not tallest 
    else
       !next shorter cohort to new cohort is the next shorter cohort 
       !to the cohort just taller than the new cohort 
       shortptr => tallptr%shorter

       !new cohort becomes the next shorter cohort to the cohort 
       !just taller than the new cohort 
       tallptr%shorter => icohort
    endif

    !new cohort is shortest 
    if (.not.associated(shortptr)) then
       !next shorter reamins null 
       !cohort is placed at the bottom of the list 
       pshortest => icohort
       if (present(storesmallcohort)) then
          storesmallcohort => icohort 
       end if
       currentPatch%shortest => icohort  
       icohort%patchptr%shortest => icohort 
    else
       !new cohort is not shortest and becomes next taller cohort 
       !to the cohort just below it as defined in the previous block 
       shortptr%taller => icohort
    endif

    ! assign taller and shorter links for the new cohort 
    icohort%taller => tallptr
    if (tallptrnull == 1) then 
       icohort%taller=> null()
    endif
    icohort%shorter => shortptr

  end subroutine insert_cohort

  !-------------------------------------------------------------------------------------!
  subroutine copy_cohort( currentCohort,copyc )
    !
    ! !DESCRIPTION:
    ! Copies all the variables in one cohort into another empty cohort                                    
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_cohort_type), intent(inout) , target ::  copyc         ! New cohort argument.
    type(ed_cohort_type), intent(in)    , target ::  currentCohort ! Old cohort argument.
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer ::  n,o           ! New and old cohort pointers
    !----------------------------------------------------------------------

    o => currentCohort
    n => copyc

    udata%cohort_number = udata%cohort_number + 1
    n%indexnumber = udata%cohort_number

    ! VEGETATION STRUCTURE
    n%pft             = o%pft
    n%n               = o%n                         
    n%dbh             = o%dbh                                        
    n%hite            = o%hite                                
    n%b               = o%b                            
    n%balive          = o%balive
    n%bdead           = o%bdead                          
    n%bstore          = o%bstore
    n%laimemory       = o%laimemory
    n%bsw             = o%bsw                  
    n%bl              = o%bl
    n%br              = o%br
    n%lai             = o%lai                         
    n%sai             = o%sai  
    n%gscan           = o%gscan
    n%leaf_cost       = o%leaf_cost
    n%canopy_layer    = o%canopy_layer
    n%nv              = o%nv
    n%status_coh      = o%status_coh
    n%canopy_trim     = o%canopy_trim
    n%status_coh      = o%status_coh               
    n%excl_weight     = o%excl_weight               
    n%prom_weight     = o%prom_weight               

    ! CARBON FLUXES 
    n%gpp             = o%gpp
    n%gpp_acc         = o%gpp_acc
    n%gpp_clm         = o%gpp_clm
    n%npp             = o%npp
    n%npp_clm         = o%npp_clm
    n%npp_acc         = o%npp_acc
    n%resp_clm        = o%resp_clm
    n%resp_acc        = o%resp_acc
    n%resp            = o%resp
    n%year_net_uptake = o%year_net_uptake
    n%ts_net_uptake   = o%ts_net_uptake

    !RESPIRATION
    n%rd              = o%rd
    n%resp_m          = o%resp_m
    n%resp_g          = o%resp_g
    n%livestem_mr     = o%livestem_mr
    n%livecroot_mr    = o%livecroot_mr
    n%froot_mr        = o%froot_mr
 
    ! NITROGEN POOLS      
    n%livestemn       = o%livestemn
    n%livecrootn      = o%livecrootn
    n%frootn          = o%frootn

    ! ALLOCATION
    n%md              = o%md
    n%leaf_md         = o%leaf_md
    n%root_md         = o%root_md
    n%carbon_balance  = o%carbon_balance
    n%dmort           = o%dmort
    n%seed_prod       = o%seed_prod
    n%treelai         = o%treelai
    n%treesai         = o%treesai
    n%leaf_litter     = o%leaf_litter
    n%c_area          = o%c_area
    n%woody_turnover  = o%woody_turnover

    ! VARIABLES NEEDED FOR INTEGRATION 
    n%dndt            = o%dndt
    n%dhdt            = o%dhdt
    n%ddbhdt          = o%ddbhdt
    n%dbalivedt       = o%dbalivedt
    n%dbdeaddt        = o%dbdeaddt
    n%dbstoredt       = o%dbstoredt
    n%storage_flux    = o%storage_flux

    ! FIRE 
    n%cfa             = o%cfa
    n%fire_mort       = o%fire_mort
    n%crownfire_mort  = o%crownfire_mort
    n%cambial_mort    = o%cambial_mort

    !Pointers
    n%taller          => NULL()     ! pointer to next tallest cohort     
    n%shorter         => NULL()     ! pointer to next shorter cohort     
    n%patchptr        => o%patchptr ! pointer to patch that cohort is in 
    n%siteptr         => o%siteptr  ! pointer to site that cohort is in  

  end subroutine copy_cohort

  !-------------------------------------------------------------------------------------!
  function count_cohorts( currentPatch ) result ( backcount )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    !
    ! !ARGUMENTS    
    type(ed_patch_type), intent(inout), target :: currentPatch      !new site
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer ::currentCohort   !new patch
    integer backcount
    !----------------------------------------------------------------------

    currentCohort => currentPatch%shortest

    currentPatch%countcohorts = 0
    do while (associated(currentCohort)) 
       currentPatch%countcohorts = currentPatch%countcohorts + 1 
       currentCohort => currentCohort%taller  
    enddo

    backcount = 0
    currentCohort => currentPatch%tallest
    do while (associated(currentCohort)) 
       backcount = backcount + 1
       currentCohort => currentCohort%shorter    
    enddo

    if (backcount /= currentPatch%countcohorts) then
       write(iulog,*) 'problem with linked list, not symmetrical' 
    endif

  end function count_cohorts

  !-------------------------------------------------------------------------------------!
  function countCohorts( bounds, ed_allsites_inst ) result ( totNumCohorts ) 
    !
    ! !DESCRIPTION:
    !  counts the total number of cohorts over all p levels (ed_patch_type) so we
    ! can allocate vectors, copy from LL -> vector and read/write restarts.
    !
    ! !USES:
    use decompMod, only : bounds_type
    !
    ! !ARGUMENTS    
    type(bounds_type)  , intent(in)            :: bounds 
    type(ed_site_type) , intent(inout), target :: ed_allsites_inst( bounds%begg: )
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type)  , pointer :: currentPatch
    type (ed_cohort_type) , pointer :: currentCohort
    integer :: g, totNumCohorts
    logical :: error
    !----------------------------------------------------------------------

    totNumCohorts = 0

    do g = bounds%begg,bounds%endg

       if (ed_allsites_inst(g)%istheresoil) then   

          currentPatch => ed_allsites_inst(g)%oldest_patch
          do while(associated(currentPatch))

             currentCohort => currentPatch%shortest
             do while(associated(currentCohort))        
                totNumCohorts = totNumCohorts + 1
                currentCohort => currentCohort%taller
             enddo !currentCohort
             currentPatch => currentPatch%younger
          end do

       end if
    end do

  end function countCohorts

end module EDCohortDynamicsMod
