module EDPatchDynamicsMod

  ! ============================================================================
  ! Controls formation, creation, fusing and termination of patch level processes. 
  ! ============================================================================

  use shr_kind_mod         , only : r8 => shr_kind_r8;
  use clm_varpar           , only : nclmax
  use clm_varctl           , only : iulog 

  use EcophysConType       , only : ecophyscon
  use EDEcophysContype     , only : EDecophyscon
  use EDtypesMod           , only : site, patch, cohort, udata, AREA, DBHMAX
  use EDtypesMod           , only : DG_SF, GRIDCELLEDSTATE, NCWD, N_DBH_BINS, NTOL, NUMPFT_ED
  use EDCohortDynamicsMod  , only : fuse_cohorts, terminate_cohorts, sort_cohorts, insert_cohort
  use EDCohortDynamicsMod  , only : zero_cohort, copy_cohort
  use EDGrowthFunctionsMod , only : mortality_rates

  implicit none
  save
  private

  public :: spawn_patches
  public :: create_patch
  public :: zero_patch
  public :: fuse_patches
  public :: fuse_2_patches
  public :: terminate_patches
  public :: patch_pft_size_profile
  public :: disturbance_rates
  public :: check_patch_area
  public :: set_patchno
  
  ! ============================================================================
  ! 10/30/09: Created by Rosie Fisher
  ! ============================================================================

contains

   subroutine disturbance_rates( site_in )
  ! ============================================================================
  ! Calculates the fire and mortality related disturbance rates for each patch,
  ! and then determines which is the larger at the patch scale (for now, there an only
  ! be one disturbance type for each timestep.  
  ! all disturbance rates here are per daily timestep. 
  ! ============================================================================

    use EDGrowthFunctionsMod, only : c_area

    implicit none

    type (site) , intent(inout) :: site_in

    type (patch) , pointer :: currentPatch
    type (cohort), pointer :: currentCohort

    !MORTALITY
    site_in%disturbance_mortality = 0.0_r8

    currentPatch => site_in%oldest_patch

    do while (associated(currentPatch))   

       currentCohort => currentPatch%shortest

       do while(associated(currentCohort))        
          ! Mortality for trees in the understorey.
          currentCohort%patchptr => currentPatch

          currentCohort%dmort  = mortality_rates(currentCohort)
          currentCohort%c_area = c_area(currentCohort)

          if(currentCohort%canopy_layer == 1)then

             currentPatch%disturbance_rates(1) = currentPatch%disturbance_rates(1) + &
                  min(1.0_r8,currentCohort%dmort)*udata%deltat*currentCohort%c_area/currentPatch%area

          endif

          currentCohort => currentCohort%taller

       enddo !currentCohort

       ! if fires occur at site 
       ! Fudge - fires can't burn the whole patch, as this causes /0 errors.
       ! This is accumulating the daily fires over the whole 30 day patch generation phase.  
       currentPatch%disturbance_rates(2) = min(0.99_r8,currentPatch%disturbance_rates(2) + currentPatch%frac_burnt)

       if(currentPatch%disturbance_rates(2) > 0.98_r8)then
          write(iulog,*) 'very high fire areas',currentPatch%disturbance_rates(2),currentPatch%frac_burnt
       endif

       !Only use larger of two natural disturbance modes WHY?
       if(currentPatch%disturbance_rates(2) > currentPatch%disturbance_rates(1))then  ! DISTURBANCE IS FIRE
          currentPatch%disturbance_rate = currentPatch%disturbance_rates(2)
       else  
          currentPatch%disturbance_rate = currentPatch%disturbance_rates(1)            ! DISTURBANCE IS MORTALITY
       endif

       site_in%disturbance_mortality = site_in%disturbance_mortality + &
            currentPatch%disturbance_rates(1)*currentPatch%area/area  
       currentPatch => currentPatch%younger

    enddo !patch loop 

    ! FIRE
    site_in%disturbance_fire = site_in%frac_burnt/AREA

    ! Use largest disturbance mode and ignore the other... This is necessary to 
    ! have a single type of disturbance and to calculate the survival rates etc... 
    if  (site_in%disturbance_fire > site_in%disturbance_mortality) then
       site_in%disturbance_rate =  site_in%disturbance_fire
       site_in%dist_type = 2
    else
       site_in%disturbance_rate = site_in%disturbance_mortality
       site_in%dist_type = 1
    endif

  end subroutine disturbance_rates

!-------------------------------------------------------------------------------!

  subroutine spawn_patches( sitetar )
    ! ============================================================================
    ! In this subroutine, the following happens
    ! 1) the total area disturbed is calculated
    ! 2) a new patch is created
    ! 3) properties are averaged
    ! 4) litter fluxes from fire and mortality are added 
    ! 5) For mortality, plants in existing patch canopy are killed. 
    ! 6) For mortality, Plants in new and existing understorey are killed
    ! 7) For fire, burned plants are killed, and unburned plants are added to new patch. 
    ! 8) New cohorts are added to new patch and sorted. 
    ! 9) New patch is added into linked list
    ! 10) Area checked, and patchno recalculated. 
    ! ============================================================================

    use EDParamsMod, only : ED_val_maxspread, &
         ED_val_understorey_death

    implicit none 

    type (site),  intent(inout), target :: sitetar

    type (site), pointer :: currentSite

    type (patch), pointer :: new_patch, currentPatch
    type (cohort),pointer :: currentCohort, nc

    real(r8) :: site_areadis        ! total area disturbed in m2 per site per day
    real(r8) :: patch_site_areadis  ! total area disturbed in m2 per patch per day
    real(r8) :: age                 ! notional age of this patch in years
    integer  :: tnull               ! is there a tallest cohort?
    integer  :: snull               ! is there a shortest cohort?
    real(r8) :: root_litter_local(numpft_ed) ! initial value of root litter. KgC/m2
    real(r8) :: leaf_litter_local(numpft_ed) ! initial value of leaf litter. KgC/m2
    real(r8) :: cwd_ag_local(ncwd)           ! initial value of above ground coarse woody debris. KgC/m2
    real(r8) :: cwd_bg_local(ncwd)           ! initial value of below ground coarse woody debris. KgC/m2
    real(r8) :: seed_bank_local(numpft_ed)   ! initial value of seed bank. KgC/m2
    real(r8) :: spread_local(nclmax)         ! initial value of canopy spread parameter.no units 

    currentSite => sitetar
    ! calculate area of disturbed land, in this timestep, by summing contributions from each existing patch. 
    currentPatch => currentSite%youngest_patch
    currentSite%cwd_ag_burned       = 0.0_r8
    currentSite%leaf_litter_burned  = 0.0_r8

    site_areadis = 0.0_r8
    do while(associated(currentPatch))

       !FIX(RF,032414) Does using the max(fire,mort) actually make sense here?
       site_areadis = site_areadis + currentPatch%area * min(1.0_r8,currentPatch%disturbance_rate) 
       currentPatch => currentPatch%older     

    enddo ! end loop over patches. sum area disturbed for all patches. 

    if (site_areadis > 0.0_r8) then  
       cwd_ag_local = 0.0_r8
       cwd_bg_local = 0.0_r8
       leaf_litter_local = 0.0_r8
       root_litter_local = 0.0_r8
       spread_local(1:nclmax) = ED_val_maxspread
       age = 0.0_r8
       seed_bank_local = 0.0_r8

       allocate(new_patch)

       call zero_patch(new_patch)

       call create_patch(currentSite,new_patch,age,site_areadis,spread_local,cwd_ag_local,cwd_bg_local,leaf_litter_local, &
            root_litter_local,seed_bank_local)

       new_patch%tallest  => null()
       new_patch%shortest => null()

       currentPatch => currentSite%oldest_patch
       ! loop round all the patches that contribute surviving indivduals and litter pools to the new patch.     
       do while(associated(currentPatch))   
          patch_site_areadis = currentPatch%area * currentPatch%disturbance_rate ! how much land is disturbed in this donor patch? 

          call average_patch_properties(currentPatch,new_patch,patch_site_areadis)               
          if (currentSite%disturbance_mortality > currentSite%disturbance_fire) then !mortality is dominant disturbance
             call mortality_litter_fluxes(currentPatch,new_patch,patch_site_areadis)  
          else
             call fire_litter_fluxes(currentPatch,new_patch,patch_site_areadis)  
          endif

          !INSERT SURVIVORS FROM DISTURBANCE INTO NEW PATCH 
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))       

             allocate(nc)             
             call zero_cohort(nc)

             ! nc is the new cohort that goes in the disturbed patch (new_patch)... currentCohort
             ! is the curent cohort that stays in the donor patch (currentPatch) 
             call copy_cohort(currentCohort,nc)

             !this is the case as the new patch probably doesn't have a closed canopy, and
             ! even if it does, that will be sorted out in canopy_structure. 
             nc%canopy_layer = 1 

             !mortality is dominant disturbance              
             if(currentPatch%disturbance_rates(1) > currentPatch%disturbance_rates(2))then 
                if(currentCohort%canopy_layer == 1)then         
                   ! keep the trees that didn't die
                   currentCohort%n = currentCohort%n * (1.0_r8 - min(1.0_r8,currentCohort%dmort * udata%deltat))
                   nc%n = 0.0_r8 ! kill all of the trees who caused the disturbance.         
                else 
                   if(ecophyscon%woody(currentCohort%pft) == 1)then

                      ! remaining of understory plants of those that are knocked over by the overstorey trees dying...  
                      nc%n = (1.0_r8 - ED_val_understorey_death) * currentCohort%n * patch_site_areadis/currentPatch%area 
                      ! understory trees that might potentially be knocked over in the disturbance. 
                      currentCohort%n = currentCohort%n * (1._r8 -  patch_site_areadis/currentPatch%area)
                      ! grass is not killed by mortality disturbance events. Just move it into the new patch area. 

                   else 

                      ! remaining of understory plants of those that are knocked over by the overstorey trees dying...  
                      nc%n = currentCohort%n * patch_site_areadis/currentPatch%area  
                      ! understory trees that might potentially be knocked over in the disturbance.                                
                      currentCohort%n = currentCohort%n * (1._r8 - patch_site_areadis/currentPatch%area)

                   endif
                endif
             else !fire

                ! loss of individual from fire in new patch.
                nc%n = currentCohort%n * patch_site_areadis/currentPatch%area * (1.0_r8 - currentCohort%fire_mort) 
                ! loss of individuals from source patch
                currentCohort%n = currentCohort%n * (1._r8 - patch_site_areadis/currentPatch%area) 

             endif

             if(nc%n > 0.0_r8)then   
                udata%storebigcohort   =>  new_patch%tallest
                udata%storesmallcohort =>  new_patch%shortest 
                if(associated(new_patch%tallest))then
                   tnull = 0
                else
                   tnull = 1
                   new_patch%tallest => nc
                   nc%taller => null()
                endif

                if(associated(new_patch%shortest))then
                   snull = 0
                else
                   snull = 1
                   new_patch%shortest => nc
                   nc%shorter => null()
                endif
                nc%patchptr => new_patch
                call insert_cohort(nc,new_patch%tallest,new_patch%shortest,tnull,snull)
                new_patch%tallest  => udata%storebigcohort 
                new_patch%shortest => udata%storesmallcohort   
             else
                deallocate(nc) !get rid of the new memory.
             endif

             currentCohort => currentCohort%taller      
          enddo ! currentCohort 
          call sort_cohorts(currentPatch) 

          !zero disturbance accumulators
          currentPatch%disturbance_rate  = 0._r8
          currentPatch%disturbance_rates = 0._r8

          !update area of donor patch
          currentPatch%area = currentPatch%area - patch_site_areadis

          !sort out the cohorts, since some of them may be so small as to need removing. 
          call fuse_cohorts(currentPatch)
          call terminate_cohorts(currentPatch)
          call sort_cohorts(currentPatch)

          currentPatch => currentPatch%younger

       enddo ! currentPatch patch loop. 

       !*************************/
       !**  INSERT NEW PATCH INTO LINKED LIST    
       !**********`***************/        
       currentPatch => currentSite%youngest_patch
       new_patch%older => currentPatch
       new_patch%younger => NULL()
       currentPatch%younger => new_patch
       currentSite%youngest_patch => new_patch

       call fuse_cohorts(new_patch)
       call terminate_cohorts(new_patch)
       call sort_cohorts(new_patch)

    endif !end new_patch area 

    call check_patch_area(currentSite)
    call set_patchno(currentSite)

  end subroutine spawn_patches

!-------------------------------------------------------------------------------!

  subroutine check_patch_area( currentSite )

  ! ============================================================================
  !  Check to see that total area is not exceeded.  
  ! ============================================================================

    use clm_varctl, only : iulog 

    implicit none      

    type(site),intent(in), pointer  :: currentSite 

    real(r8) :: areatot
    type(patch), pointer :: currentPatch 

    areatot = 0._r8
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
       areatot = areatot + currentPatch%area
       currentPatch => currentPatch%younger
       if (( areatot - area ) > 0._r8 ) then 
          write(iulog,*) 'trimming patch area - is too big' , areatot-area
          currentSite%oldest_patch%area = currentSite%oldest_patch%area - (areatot - area)
       endif
    enddo

  end subroutine check_patch_area

!-------------------------------------------------------------------------------!

  subroutine set_patchno( currentSite )

  ! ============================================================================
  !  Give patches an order number from the oldest to youngest. 
  ! ============================================================================

    implicit none      

    type(site),intent(in), pointer :: currentSite 

    type(patch), pointer :: currentPatch 

    integer patchno

    patchno = 1
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
       currentPatch%patchno = patchno
       patchno = patchno + 1
       currentPatch => currentPatch%younger
    enddo

  end subroutine set_patchno

!-------------------------------------------------------------------------------!

  subroutine average_patch_properties( currentPatch, newPatch, patch_site_areadis )

  ! ============================================================================
  ! Average together the state properties of all of the donor patches that
  ! make up the new patch. 
  ! ============================================================================

    implicit none      

    type(patch),intent(in),    pointer :: currentPatch
    type(patch),intent(inout), pointer :: newPatch 
    real(r8), intent(out)              :: patch_site_areadis   ! amount of land disturbed in this patch. m2

    integer  :: c,p ! counters for PFT and litter size class. 

    patch_site_areadis = currentPatch%area * currentPatch%disturbance_rate ! how much land is disturbed in this donor patch? 
 
    do p=1,numpft_ed
       newPatch%seed_bank(p) = newPatch%seed_bank(p) + currentPatch%seed_bank(p) * patch_site_areadis/newPatch%area
    enddo

    do c = 1,ncwd !move litter pool en mass into the new patch. 
       newPatch%cwd_ag(c) = newPatch%cwd_ag(c) + currentPatch%cwd_ag(c) * patch_site_areadis/newPatch%area
       newPatch%cwd_bg(c) = newPatch%cwd_bg(c) + currentPatch%cwd_bg(c) * patch_site_areadis/newPatch%area
    enddo

    do p = 1,numpft_ed !move litter pool en mass into the new patch
       newPatch%root_litter(p) = newPatch%root_litter(p) + currentPatch%root_litter(p) * patch_site_areadis/newPatch%area
       newPatch%leaf_litter(p) = newPatch%leaf_litter(p) + currentPatch%leaf_litter(p) * patch_site_areadis/newPatch%area
    enddo

    newPatch%spread = newPatch%spread + currentPatch%spread * patch_site_areadis/newPatch%area    

  end subroutine average_patch_properties

!-------------------------------------------------------------------------------!

  subroutine fire_litter_fluxes(cp_target,new_patch_target,patch_site_areadis)
  ! ============================================================================
  !  CWD pool burned by a fire. 
  !  Carbon going from burned trees into CWD pool
  !  Burn parts of trees that don't die in fire
  !  Burn live grasses and kill them. 
  ! ============================================================================

    use EDParamsMod,          only : ED_val_ag_biomass
    use SFParamsMod,          only : SF_VAL_CWD_FRAC
    use EDGrowthFunctionsMod, only : c_area

    implicit none

    type(patch), intent(inout), target :: cp_target
    type(patch), intent(inout), target :: new_patch_target
    real(r8),    intent(inout)         :: patch_site_areadis

    type(site), pointer   :: currentSite
    type(patch), pointer  :: currentPatch
    type(patch), pointer  :: new_patch
    type(cohort), pointer :: currentCohort

    real(r8) bcroot               ! amount of below ground coarse root per cohort  kgC. (goes into CWD_BG)
    real(r8) bstem                ! amount of above ground stem biomass per cohort  kgC.(goes into CWG_AG)
    real(r8) dead_tree_density    ! no trees killed by fire per m2
    reaL(r8) burned_litter        ! amount of each litter pool burned by fire.  kgC/m2/day
    real(r8) burned_leaves        ! amount of tissue consumed by fire for grass. KgC/individual/day

    integer c, p 
    !check that total area is not exceeded. 
    currentPatch => cp_target
    new_patch => new_patch_target

    if ( currentPatch%fire  ==  1 ) then !only do this if there was a fire in this actual patch. 
       patch_site_areadis = currentPatch%area * currentPatch%disturbance_rate ! how much land is disturbed in this donor patch? 
       currentSite => currentPatch%siteptr

       !************************************/ 
       !PART 1)  Burn the fractions of existing litter in the new patch that were consumed by the fire. 
       !************************************/ 
       do c = 1,ncwd
          burned_litter = new_patch%cwd_ag(c) * patch_site_areadis/new_patch%area * currentPatch%burnt_frac_litter(c+1) !kG/m2/day
          new_patch%cwd_ag(c) = new_patch%cwd_ag(c) - burned_litter
          currentSite%flux_out = currentSite%flux_out + burned_litter * new_patch%area !kG/site/day
       enddo

       do p = 1,numpft_ed
          burned_litter = new_patch%leaf_litter(p) * patch_site_areadis/new_patch%area * currentPatch%burnt_frac_litter(dg_sf)
          new_patch%leaf_litter(p) = new_patch%leaf_litter(p) - burned_litter
          currentSite%flux_out = currentSite%flux_out + burned_litter * new_patch%area !kG/site/dat
      enddo

       !************************************/     
       !PART 2) Put unburned parts of plants that died in the fire into the litter pool of new and old patches 
       ! This happens BEFORE the plant numbers have been updated. So we are working with the 
       ! pre-fire population of plants, which is the right way round. 
       !************************************/ 
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))
          p = currentCohort%pft
          if(ecophyscon%woody(p) == 1)then !DEAD (FROM FIRE) TREES
             !************************************/ 
             ! Number of trees that died because of the fire, per m2 of ground. 
             ! Divide their litter into the four litter streams, and spread evenly across ground surface. 
             !************************************/  
             ! stem biomass per tree
             bstem  = (currentCohort%bsw + currentCohort%bdead) * ED_val_ag_biomass           
             ! coarse root biomass per tree
             bcroot = (currentCohort%bsw + currentCohort%bdead) * (1.0_r8 - ED_val_ag_biomass) 
             ! density of dead trees per m2. 
             dead_tree_density  = (currentCohort%fire_mort * currentCohort%n*patch_site_areadis/currentPatch%area) / AREA  

             ! Unburned parts of dead tree pool. 
             ! Unburned leaves and roots    
             
             new_patch%leaf_litter(p) = new_patch%leaf_litter(p) + dead_tree_density * (currentCohort%bl) &
             * (1.0_r8-currentCohort%cfa)
             new_patch%root_litter(p) = new_patch%root_litter(p) + dead_tree_density * (currentCohort%br+currentCohort%bstore)
             currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + dead_tree_density * &
                  (currentCohort%bl) * (1.0_r8-currentCohort%cfa)
             currentPatch%root_litter(p) = currentPatch%root_litter(p) + dead_tree_density * &
                  (currentCohort%br+currentCohort%bstore)
      
             ! below ground coarse woody debris from burned trees
             do c = 1,ncwd
                new_patch%cwd_bg(c) = new_patch%cwd_bg(c) + dead_tree_density * SF_val_CWD_frac(c) * bcroot
                currentPatch%cwd_bg(c) = currentPatch%cwd_bg(c) + dead_tree_density * SF_val_CWD_frac(c) * bcroot
             enddo

             ! above ground coarse woody debris from unburned twigs and small branches
             do c = 1,2
                new_patch%cwd_ag(c) = new_patch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * bstem &
                * (1.0_r8-currentCohort%cfa)
                currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * &
                     bstem * (1.0_r8-currentCohort%cfa)
             enddo
             
             ! above ground coarse woody debris from large branches and stems: these do not burn in crown fires. 
             do c = 3,4
                new_patch%cwd_ag(c) = new_patch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * bstem
                currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * bstem
             enddo
             
             ! Burned parts of dead tree pool.  
             ! Burned twigs and small branches. 
             do c = 1,2

                currentSite%cwd_ag_burned(c) = currentSite%cwd_ag_burned(c) + dead_tree_density * &
                     SF_val_CWD_frac(c) * bstem * currentCohort%cfa
                currentSite%flux_out  = currentSite%flux_out + dead_tree_density * &
                     AREA * SF_val_CWD_frac(c) * bstem * currentCohort%cfa

             enddo
             
             !burned leaves. 
             do p = 1,numpft_ed                  

                currentSite%leaf_litter_burned(p) = currentSite%leaf_litter_burned(p) + &
                     dead_tree_density * currentCohort%bl * currentCohort%cfa
                currentSite%flux_out  = currentSite%flux_out + &
                     dead_tree_density * AREA * currentCohort%bl * currentCohort%cfa

             enddo

         endif

          currentCohort => currentCohort%taller

       enddo  ! currentCohort

       !************************************/     
       ! PART 3) Burn parts of trees that did *not* die in the fire.
       ! PART 4) Burn parts of grass that are consumed by the fire. 
       ! grasses are not killed directly by fire. They die by losing all of their leaves and starving. 
       !************************************/ 
       currentCohort => new_patch%shortest
       do while(associated(currentCohort))

          currentCohort%c_area = c_area(currentCohort) 
          if(ecophyscon%woody(currentCohort%pft) == 1)then
             burned_leaves = (currentCohort%bl+currentCohort%bsw) * currentCohort%cfa
          else
             burned_leaves = (currentCohort%bl+currentCohort%bsw) * currentPatch%burnt_frac_litter(6)
          endif
          if (burned_leaves > 0_r8) then

             currentCohort%balive = max(currentCohort%br,currentCohort%balive - burned_leaves)
             currentCohort%bl     = max(0.00001_r8,   currentCohort%bl - burned_leaves)
             !KgC/gridcell/day
             currentSite%flux_out = currentSite%flux_out + burned_leaves * currentCohort%n * &
                  patch_site_areadis/currentPatch%area * AREA 

          endif
          currentCohort%cfa = 0.0_r8        

          currentCohort => currentCohort%taller

       enddo

    endif !currentPatch%fire. 

  end subroutine fire_litter_fluxes

!-------------------------------------------------------------------------------!

  subroutine mortality_litter_fluxes(cp_target,new_patch_target,patch_site_areadis)
 
  ! ============================================================================
  !  Carbon going from ongoing mortality into CWD pools. 
  ! ============================================================================

    use EDParamsMod,  only : ED_val_ag_biomass, ED_val_understorey_death
    use SFParamsMod,  only : SF_val_cwd_frac

    implicit none      

    type(patch), intent(inout), target :: cp_target 
    type(patch), intent(inout), target :: new_patch_target
    real(r8),    intent(in)            :: patch_site_areadis

    real(r8) :: cwd_litter_density
    real(r8) :: litter_area ! area over which to distribute this litter. 

    type(cohort), pointer :: currentCohort
    type(patch), pointer :: currentPatch 

    type(patch), pointer :: new_patch 
    real(r8) :: understorey_dead  !Number of individual dead from the canopy layer /day
    real(r8) :: canopy_dead       !Number of individual dead from the understorey layer /day
    real(r8) :: np_mult           !Fraction of the new patch which came from the current patch (and so needs the same litter) 

    integer p,c

    currentPatch => cp_target
    new_patch => new_patch_target
    currentPatch%canopy_mortality_woody_litter    = 0.0_r8 ! mortality generated litter. KgC/m2/day
    currentPatch%canopy_mortality_leaf_litter(:)  = 0.0_r8
    currentPatch%canopy_mortality_root_litter(:)  = 0.0_r8

    currentCohort => currentPatch%shortest
    do while(associated(currentCohort))       
       p = currentCohort%pft
       if(currentPatch%disturbance_rates(1) > currentPatch%disturbance_rates(2))then !mortality is dominant disturbance 
          if(currentCohort%canopy_layer == 1)then         
             !currentCohort%dmort = mortality_rates(currentCohort) 
             !the disturbance calculations are done with the previous n, c_area and d_mort. So it's probably &
             !not right to recalcualte dmort here.
             canopy_dead = currentCohort%n * min(1.0_r8,currentCohort%dmort*udata%deltat)

             currentPatch%canopy_mortality_woody_litter   = currentPatch%canopy_mortality_woody_litter  + &
                  canopy_dead*(currentCohort%bdead+currentCohort%bsw)
             currentPatch%canopy_mortality_leaf_litter(p) = currentPatch%canopy_mortality_leaf_litter(p)+ &
                  canopy_dead*(currentCohort%bl)
             currentPatch%canopy_mortality_root_litter(p) = currentPatch%canopy_mortality_root_litter(p)+ &
                  canopy_dead*(currentCohort%br+currentCohort%bstore)

         else 
             if(ecophyscon%woody(currentCohort%pft) == 1)then

                understorey_dead = ED_val_understorey_death * currentCohort%n * (patch_site_areadis/currentPatch%area)  !kgC/site/day
                currentPatch%canopy_mortality_woody_litter  = currentPatch%canopy_mortality_woody_litter  + &
                     understorey_dead*(currentCohort%bdead+currentCohort%bsw)  
                currentPatch%canopy_mortality_leaf_litter(p)= currentPatch%canopy_mortality_leaf_litter(p)+ &
                     understorey_dead* currentCohort%bl 
                currentPatch%canopy_mortality_root_litter(p)= currentPatch%canopy_mortality_root_litter(p)+ &
                      understorey_dead*(currentCohort%br+currentCohort%bstore)

             ! FIX(SPM,040114) - clarify this comment
             ! grass is not killed by canopy mortality disturbance events.
             ! Just move it into the new patch area. 
             else 
                ! no-op
             endif
          endif
       endif

       currentCohort => currentCohort%taller      

    enddo !currentCohort         

    !************************************/
    !Evenly distribute the litter from the trees that died across the new and old patches
    !************************************/
    !************************************/
    !Evenly distribute the litter from the trees that died across the new and old patches
    !'litter' fluxes here are in KgC
    !************************************/
    litter_area = currentPatch%area 
    np_mult =  patch_site_areadis/new_patch%area
    ! This litter is distributed between the current and new patches, &
    ! not to any other patches. This is really the eventually area of the current patch &
    ! (currentPatch%area-patch_site_areadis) +patch_site_areadis...
    ! For the new patch, only some fraction of its land area (patch_areadis/np%area) is derived from the current patch
    ! so we need to multiply by patch_areadis/np%area
    do c = 1,ncwd
    
       cwd_litter_density = SF_val_CWD_frac(c) * currentPatch%canopy_mortality_woody_litter / litter_area
       
       new_patch%cwd_ag(c)    = new_patch%cwd_ag(c)    + ED_val_ag_biomass         * cwd_litter_density * np_mult
       currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + ED_val_ag_biomass         * cwd_litter_density
       new_patch%cwd_bg(c)    = new_patch%cwd_bg(c)    + (1._r8-ED_val_ag_biomass) * cwd_litter_density * np_mult 
       currentPatch%cwd_bg(c) = currentPatch%cwd_bg(c) + (1._r8-ED_val_ag_biomass) * cwd_litter_density 
       
    enddo 

    do p = 1,numpft_ed
    
       new_patch%leaf_litter(p) = new_patch%leaf_litter(p) + currentPatch%canopy_mortality_leaf_litter(p) / litter_area * np_mult
       new_patch%root_litter(p) = new_patch%root_litter(p) + currentPatch%canopy_mortality_root_litter(p) / litter_area * np_mult 
       currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + currentPatch%canopy_mortality_leaf_litter(p) / litter_area
       currentPatch%root_litter(p) = currentPatch%root_litter(p) + currentPatch%canopy_mortality_root_litter(p) / litter_area
       
    enddo

  end subroutine mortality_litter_fluxes

!-------------------------------------------------------------------------------!


  subroutine create_patch(currentSite, new_patch, age, areap, spread_local,cwd_ag_local,cwd_bg_local, &
       leaf_litter_local,root_litter_local,seed_bank_local)
  ! ============================================================================
  !   Set default values for creating a new patch
  ! ============================================================================

    use clm_varpar      , only : nlevgrnd 

    implicit none      

    type(site),  target, intent(inout) :: currentSite
    type(patch), target, intent(inout) :: new_patch

    real(r8), intent(in) :: age                 ! notional age of this patch in years
    real(r8), intent(in) :: areap               ! initial area of this patch in m2. 
    real(r8), intent(in) :: cwd_ag_local(:)     ! initial value of above ground coarse woody debris. KgC/m2
    real(r8), intent(in) :: cwd_bg_local(:)     ! initial value of below ground coarse woody debris. KgC/m2
    real(r8), intent(in) :: root_litter_local(:)! initial value of root litter. KgC/m2
    real(r8), intent(in) :: leaf_litter_local(:)! initial value of leaf litter. KgC/m2
    real(r8), intent(in) :: spread_local(:)     ! initial value of canopy spread parameter.no units 
    real(r8), intent(in) :: seed_bank_local(:)  ! initial value of seed bank. KgC/m2
    
    call zero_patch(new_patch) !The nan value in here is not working??

    new_patch%tallest => null()       ! pointer to patch's tallest cohort    
    new_patch%shortest => null()      ! pointer to patch's shortest cohort   
    new_patch%older => null()            ! pointer to next older patch   
    new_patch%younger => null()       ! pointer to next shorter patch      
    new_patch%siteptr => null()          ! pointer to the site that the patch is in

    ! assign known patch attributes 

    new_patch%siteptr    => currentSite 

    new_patch%age         = age   
    new_patch%area        = areap 
    new_patch%spread      = spread_local
    new_patch%cwd_ag      = cwd_ag_local
    new_patch%cwd_bg      = cwd_bg_local
    new_patch%leaf_litter = leaf_litter_local
    new_patch%root_litter = root_litter_local
    new_patch%seed_bank   = seed_bank_local
 
    !zeroing things because of the surfacealbedo problem... shouldnt really be necesary
    new_patch%cwd_ag_in(:) = 0._r8
    new_patch%cwd_bg_in(:) = 0._r8

    new_patch%f_sun = 0._r8
    new_patch%ed_laisun_z(:,:,:) = 0._r8 
    new_patch%ed_laisha_z(:,:,:) = 0._r8 
    new_patch%ed_parsun_z(:,:,:) = 0._r8 
    new_patch%ed_parsha_z(:,:,:) = 0._r8 
    new_patch%fabi = 0._r8
    new_patch%fabd = 0._r8
    new_patch%tr_soil_dir(:) = 1._r8
    new_patch%tr_soil_dif(:) = 1._r8
    new_patch%tr_soil_dir_dif(:) = 0._r8
    new_patch%fabd_sun_z(:,:,:) = 0._r8 
    new_patch%fabd_sha_z(:,:,:) = 0._r8 
    new_patch%fabi_sun_z(:,:,:) = 0._r8 
    new_patch%fabi_sha_z(:,:,:) = 0._r8  
    new_patch%frac_burnt = 0._r8  
    new_patch%total_tree_area = 0.0_r8  
    new_patch%NCL_p = 1

    allocate(new_patch%rootfr_ft(numpft_ed,nlevgrnd))
    allocate(new_patch%rootr_ft(numpft_ed,nlevgrnd)) 

  end subroutine create_patch

!-------------------------------------------------------------------------------!

  subroutine zero_patch(cp_p)
  ! ============================================================================
  !  Sets all the variables in the patch to nan or zero (this needs to be two seperate routines, one for nan & one for zero
  ! ============================================================================

    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)  

    implicit none   

    type(patch), target, intent(inout) :: cp_p

    type(patch), pointer :: currentPatch

    currentPatch  => cp_p  

    currentPatch%tallest => null()          
    currentPatch%shortest => null()         
    currentPatch%older => null()               
    currentPatch%younger => null()           
    currentPatch%siteptr => null()             

    currentPatch%patchno  = 999                            
    currentPatch%clm_pno  = 999                          

    currentPatch%age = nan                          
    currentPatch%area = nan                                           
    currentPatch%canopy_layer_lai(:) = nan               
    currentPatch%total_canopy_area = nan
    currentPatch%canopy_area = nan                                 
    currentPatch%bare_frac_area = nan                             

    currentPatch%tlai_profile(:,:,:) = nan 
    currentPatch%elai_profile(:,:,:) = nan 
    currentPatch%tsai_profile(:,:,:) = nan 
    currentPatch%esai_profile(:,:,:) = nan       
    currentPatch%canopy_area_profile(:,:,:) = nan       

    currentPatch%fabd_sun_z(:,:,:) = nan 
    currentPatch%fabd_sha_z(:,:,:) = nan 
    currentPatch%fabi_sun_z(:,:,:) = nan 
    currentPatch%fabi_sha_z(:,:,:) = nan  

    currentPatch%ed_laisun_z(:,:,:) = nan 
    currentPatch%ed_laisha_z(:,:,:) = nan 
    currentPatch%ed_parsun_z(:,:,:) = nan 
    currentPatch%ed_parsha_z(:,:,:) = nan 
    currentPatch%psn_z(:,:,:) = nan   

    currentPatch%f_sun(:,:,:) = nan
    currentPatch%tr_soil_dir(:) = nan    ! fraction of incoming direct  radiation that is transmitted to the soil as direct
    currentPatch%tr_soil_dif(:) = nan    ! fraction of incoming diffuse radiation that is transmitted to the soil as diffuse
    currentPatch%tr_soil_dir_dif(:) = nan! fraction of incoming direct  radiation that is transmitted to the soil as diffuse
    currentPatch%fab(:) = nan            ! fraction of incoming total   radiation that is absorbed by the canopy
    currentPatch%fabd(:) = nan           ! fraction of incoming direct  radiation that is absorbed by the canopy
    currentPatch%fabi(:) = nan           ! fraction of incoming diffuse radiation that is absorbed by the canopy

    currentPatch%present(:,:) = 999              ! is there any of this pft in this layer?
    currentPatch%nrad(:,:) = 999                 ! number of exposed leaf layers for each canopy layer and pft
    currentPatch%ncan(:,:) = 999                 ! number of total leaf layers for each canopy layer and pft
    currentPatch%lai = nan                       ! leaf area index of patch
    currentPatch%spread(:) = nan                 ! dynamic ratio of dbh to canopy area.
    currentPatch%pft_agb_profile(:,:) = nan    
    currentPatch%gpp = 0._r8 
    currentPatch%npp = 0._r8                
    currentPatch%seed_bank(:) = 0._r8                    
    currentPatch%dseed_dt(:) = 0._r8                    

    ! DISTURBANCE 
    currentPatch%disturbance_rates = 0._r8 
    currentPatch%disturbance_rate = 0._r8 

    ! LITTER
    currentPatch%cwd_ag(:) = 0.0_r8               ! above ground coarse woody debris gc/m2. 
    currentPatch%cwd_bg(:) = 0.0_r8               ! below ground coarse woody debris
    currentPatch%root_litter(:) = 0.0_r8
    currentPatch%leaf_litter(:) = 0.0_r8

    ! FIRE
    currentPatch%fuel_eff_moist = 0.0_r8    ! average fuel moisture content of the ground fuel 
                                            ! (incl. live grasses. omits 1000hr fuels)
    currentPatch%livegrass = 0.0_r8         ! total ag grass biomass in patch. 1=c3 grass, 2=c4 grass. gc/m2
    currentPatch%sum_fuel = 0.0_r8          ! total ground fuel related to ros (omits 1000hr fuels). gc/m2
    currentPatch%fuel_bulkd  = 0.0_r8       ! average fuel bulk density of the ground fuel 
                                            ! (incl. live grasses. omits 1000hr fuels). kgc/m3
    currentPatch%fuel_sav  = 0.0_r8         ! average surface area to volume ratio of the ground fuel 
                                            ! (incl. live grasses. omits 1000hr fuels).
    currentPatch%fuel_mef  = 0.0_r8         ! average moisture of extinction factor of the ground fuel
                                            ! (incl. live grasses. omits 1000hr fuels).
    currentPatch%ros_front  = 0.0_r8        ! average rate of forward spread of each fire in the patch. m/min.
    currentPatch%effect_wspeed = 0.0_r8     ! dailywind modified by fraction of relative grass and tree cover. m/min.
    currentPatch%tau_l = 0.0_r8             ! mins p&r(1986)
    currentPatch%fuel_frac(:) = 0.0_r8      ! fraction of each litter class in the sum_fuel 
                                            !- for purposes of calculating weighted averages. 
    currentPatch%tfc_ros = 0.0_r8           ! used in fi calc
    currentPatch%fi = 0._r8                 ! average fire intensity of flaming front during day.  
                                            ! backward ros plays no role. kj/m/s or kw/m.
    currentPatch%fire = 999                 ! sr decide_fire.1=fire hot enough to proceed. 0=stop everything- no fires today
    currentPatch%fd = 0.0_r8                ! fire duration (mins)
    currentPatch%ros_back = 0.0_r8          ! backward ros (m/min)
    currentPatch%ab = 0.0_r8                ! area burnt daily m2
    currentPatch%nf = 0.0_r8                ! number of fires initiated daily 
    currentPatch%sh  = 0.0_r8               ! average scorch height for the patch(m)
    currentPatch%frac_burnt = 0.0_r8        ! fraction burnt in each timestep. 
    currentPatch%burnt_frac_litter(:) = 0.0_r8 
    currentPatch%btran_ft(:) = 0.0_r8

    currentPatch%canopy_layer_lai(:) = 0.0_r8
    currentPatch%seeds_in(:)         = 0.0_r8
    currentPatch%seed_decay(:)       = 0.0_r8
    currentPatch%seed_germination(:) = 0.0_r8
    currentPatch%fab(:)              = 0.0_r8
    currentPatch%sabs_dir(:)         = 0.0_r8
    currentPatch%sabs_dif(:)         = 0.0_r8


  end subroutine zero_patch

!-------------------------------------------------------------------------------!

  subroutine fuse_patches( csite )
  ! ============================================================================
  !      Decide to fuse patches if their cohort structures are similar           
  ! ============================================================================

    use clm_varctl,   only : iulog    

    implicit none   

    type(site), target, intent(inout)  :: csite

    type(site), pointer :: currentSite
    type(patch), pointer :: currentPatch,tpp,tmpptr
    integer ft,z            !counters for pft and height class
    real(r8) :: norm        !normalized difference between biomass profiles
    real(r8) :: profiletol  !tolerance of patch fusion routine. Starts off high and is reduced if there are too many patches.
    integer  :: maxpatch    !maximum number of allowed patches. FIX-RF. These should be namelist variables. 
    integer  :: nopatches   !number of patches presently in gridcell
    integer  :: iterate     !switch of patch reduction iteration scheme. 1 to keep going, 0 to stop
    integer  :: fuse_flag   !do patches get fused (1) or not (0). 

    maxpatch = 4  

    currentSite => csite 

    profiletol = 0.6_r8 !start off with a very small profile tol, or a predefined parameter? 

    nopatches = 0
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       nopatches = nopatches +1
       currentPatch => currentPatch%older
    enddo
    !---------------------------------------------------------------------!
    !  We only really care about fusing patches if nopatches > 1           !
    !---------------------------------------------------------------------!
    iterate = 1

    !---------------------------------------------------------------------!
    !  Keep doing this until nopatches >= maxpatch                         !
    !---------------------------------------------------------------------!

    do while(iterate == 1)
       !---------------------------------------------------------------------!
       ! Calculate the biomass profile of each patch                         !
       !---------------------------------------------------------------------!  
       currentPatch => currentSite%youngest_patch
       do while(associated(currentPatch))
          call patch_pft_size_profile(currentPatch)
          currentPatch => currentPatch%older
       enddo

       !---------------------------------------------------------------------!
       ! Loop round current & target (currentPatch,tpp) patches to assess combinations !
       !---------------------------------------------------------------------!   
       currentPatch => currentSite%youngest_patch
       do while(associated(currentPatch))      
          tpp => currentSite%youngest_patch
          do while(associated(tpp))

             if(.not.associated(currentPatch))then
                write(iulog,*) 'ED: issue with currentPatch'
             endif

             if(associated(tpp).and.associated(currentPatch))then
                fuse_flag = 1 !the default is to fuse the patches
                if(currentPatch%patchno /= tpp%patchno) then   !these should be the same patch

                   !---------------------------------------------------------------------!
                   ! Calculate the difference criteria for each pft and dbh class        !
                   !---------------------------------------------------------------------!   
                   do ft = 1,numpft_ed        ! loop over pfts
                      do z = 1,n_dbh_bins      ! loop over hgt bins 
                         !is there biomass in this category?
                         if(currentPatch%pft_agb_profile(ft,z)  > 0.0_r8.or.tpp%pft_agb_profile(ft,z) > 0.0_r8)then 
                            norm = abs(currentPatch%pft_agb_profile(ft,z) - tpp%pft_agb_profile(ft,z))/(0.5_r8*&
                                 &(currentPatch%pft_agb_profile(ft,z) + tpp%pft_agb_profile(ft,z)))
                            !---------------------------------------------------------------------!
                            ! Look for differences in profile biomass, above the minimum biomass  !
                            !---------------------------------------------------------------------!

                            if(norm  > profiletol)then
                               !looking for differences between profile density.                   
                               if(currentPatch%pft_agb_profile(ft,z) > NTOL.or.tpp%pft_agb_profile(ft,z) > NTOL)then
                                  fuse_flag = 0 !do not fuse  - keep apart. 
                               endif
                            endif ! profile tol           
                         endif ! NTOL 
                      enddo !ht bins
                   enddo ! PFT

                   !---------------------------------------------------------------------!
                   ! Call the patch fusion routine if there is a meaningful difference   !
                   ! any of the pft x height categories                                  !
                   !---------------------------------------------------------------------!

                   if(fuse_flag  ==  1)then 
                      tmpptr => currentPatch%older       
                      call fuse_2_patches(currentPatch,tpp)                  
                      call fuse_cohorts(tpp)
                      call sort_cohorts(tpp)
                      currentPatch => tmpptr
                   else
                     ! write(iulog,*) 'patches not fused'
                   endif
                endif  !are both patches associated?        
             endif    !are these different patches?   
             tpp => tpp%older
          enddo !tpp loop

          if(associated(currentPatch))then 
             currentPatch => currentPatch%older 
          else
             currentPatch => null()
          endif !associated currentPatch

       enddo ! currentPatch loop

       !---------------------------------------------------------------------!
       ! Is the number of patches larger than the maximum?                   !
       !---------------------------------------------------------------------!   
       nopatches = 0
       currentPatch => currentSite%youngest_patch
       do while(associated(currentPatch))
          nopatches = nopatches +1
          currentPatch => currentPatch%older
       enddo

       if(nopatches > maxpatch)then
          iterate = 1
          profiletol = profiletol * 1.1_r8
          write(iulog,*) 'maxpatch exceeded, triggering patch fusion iteration.',profiletol,nopatches
          !---------------------------------------------------------------------!
          ! Making profile tolerance larger means that more fusion will happen  !
          !---------------------------------------------------------------------!        
       else
          iterate = 0
       endif

    enddo !do while nopatches>maxpatch
 
  end subroutine fuse_patches

!-------------------------------------------------------------------------------!

  subroutine fuse_2_patches(dp, rp)
  ! ============================================================================
  ! This function fuses the two patches specified in the argument.
  ! It fuses the first patch in the argument (the "donor") into the second
  ! patch in the argument (the "recipient"), and frees the memory 
  ! associated with the secnd patch
  ! ============================================================================

    implicit none

    type (patch),  intent(inout), pointer  :: dp ! Donor Patch
    type (patch),  intent(inout), pointer  :: rp ! Recipient Patch

    type (cohort), pointer :: currentCohort ! Current Cohort
    type (cohort), pointer :: nextc         ! Remembers next cohort in list 

    integer :: c,p !counters for pft and litter size class. 
    integer :: tnull,snull  ! are the tallest and shortest cohorts associated?

    !area weighted average of ages & litter & seed bank
     rp%age = (dp%age * dp%area + rp%age * rp%area)/(dp%area + rp%area)  

    do p = 1,numpft_ed
       rp%seed_bank(p)        = (rp%seed_bank(p)*rp%area + dp%seed_bank(p)*dp%area)/(rp%area + dp%area)
       rp%seeds_in(p)         = (rp%seeds_in(p)*rp%area + dp%seeds_in(p)*dp%area)/(rp%area + dp%area)
       rp%seed_decay(p)       = (rp%seed_decay(p)*rp%area + dp%seed_decay(p)*dp%area)/(rp%area + dp%area)
       rp%seed_germination(p) = (rp%seed_germination(p)*rp%area + dp%seed_germination(p)*dp%area)/(rp%area + dp%area)
    enddo

    do c = 1,ncwd
       rp%cwd_ag(c) = (dp%cwd_ag(c)*dp%area + rp%cwd_ag(c)*rp%area)/(dp%area + rp%area)
       rp%cwd_bg(c) = (dp%cwd_bg(c)*dp%area + rp%cwd_bg(c)*rp%area)/(dp%area + rp%area)
    enddo

    do p = 1,numpft_ed
       rp%leaf_litter(p) = (dp%leaf_litter(p)*dp%area + rp%leaf_litter(p)*rp%area)/(dp%area + rp%area)
       rp%root_litter(p) = (dp%root_litter(p)*dp%area + rp%root_litter(p)*rp%area)/(dp%area + rp%area)
    enddo
    
    rp%fuel_eff_moist       = (dp%fuel_eff_moist*dp%area + rp%fuel_eff_moist*rp%area)/(dp%area + rp%area)
    rp%livegrass            = (dp%livegrass*dp%area + rp%livegrass*rp%area)/(dp%area + rp%area)
    rp%sum_fuel             = (dp%sum_fuel*dp%area + rp%sum_fuel*rp%area)/(dp%area + rp%area)
    rp%fuel_bulkd           = (dp%fuel_bulkd*dp%area + rp%fuel_bulkd*rp%area)/(dp%area + rp%area)
    rp%fuel_sav             = (dp%fuel_sav*dp%area + rp%fuel_sav*rp%area)/(dp%area + rp%area)
    rp%fuel_mef             = (dp%fuel_mef*dp%area + rp%fuel_mef*rp%area)/(dp%area + rp%area)
    rp%ros_front            = (dp%ros_front*dp%area + rp%ros_front*rp%area)/(dp%area + rp%area)
    rp%effect_wspeed        = (dp%effect_wspeed*dp%area + rp%effect_wspeed*rp%area)/(dp%area + rp%area)
    rp%tau_l                = (dp%tau_l*dp%area + rp%tau_l*rp%area)/(dp%area + rp%area)
    rp%fuel_frac(:)         = (dp%fuel_frac(:)*dp%area + rp%fuel_frac(:)*rp%area)/(dp%area + rp%area)
    rp%tfc_ros              = (dp%tfc_ros*dp%area + rp%tfc_ros*rp%area)/(dp%area + rp%area)
    rp%fi                   = (dp%fi*dp%area + rp%fi*rp%area)/(dp%area + rp%area)
    rp%fd                   = (dp%fd*dp%area + rp%fd*rp%area)/(dp%area + rp%area)
    rp%ros_back             = (dp%ros_back*dp%area + rp%ros_back*rp%area)/(dp%area + rp%area)
    rp%ab                   = (dp%ab*dp%area + rp%ab*rp%area)/(dp%area + rp%area)
    rp%nf                   = (dp%nf*dp%area + rp%nf*rp%area)/(dp%area + rp%area)
    rp%sh                   = (dp%sh*dp%area + rp%sh*rp%area)/(dp%area + rp%area)
    rp%frac_burnt           = (dp%frac_burnt*dp%area + rp%frac_burnt*rp%area)/(dp%area + rp%area)
    rp%burnt_frac_litter(:) = (dp%burnt_frac_litter(:)*dp%area + rp%burnt_frac_litter(:)*rp%area)/(dp%area + rp%area)
    rp%btran_ft(:)          = (dp%btran_ft(:)*dp%area + rp%btran_ft(:)*rp%area)/(dp%area + rp%area)
    
    rp%area = rp%area + dp%area !THIS MUST COME AT THE END!

    !insert donor cohorts into recipient patch
    if(associated(dp%shortest))then

       currentCohort => dp%shortest
       if(associated(currentCohort)) then
          nextc => currentCohort%taller
       endif

       do while(associated(dp%shortest))

          udata%storebigcohort   => rp%tallest
          udata%storesmallcohort =>  rp%shortest

          if(associated(rp%tallest))then
             tnull = 0
          else
             tnull = 1
             rp%tallest => currentCohort
          endif

          if(associated(rp%shortest))then
             snull = 0
          else
             snull = 1
             rp%shortest => currentCohort
          endif

          call insert_cohort(currentCohort,rp%tallest,rp%shortest,tnull,snull)

          rp%tallest => udata%storebigcohort 
          rp%shortest => udata%storesmallcohort    
          currentCohort%patchptr => rp
          currentCohort => nextc

          dp%shortest => currentCohort

          if(associated(currentCohort)) then
             nextc => currentCohort%taller
          endif

       enddo !cohort
    endif !are there any cohorts?

    call patch_pft_size_profile(rp) ! Recalculate the patch size profile for the resulting patch

    ! FIX(SPM,032414) dangerous code here.  Passing in dp as a pointer allows the code below
    ! to effect the currentPatch that is the actual argument when in reality, dp should be 
    ! intent in only with these pointers being set on the actual argument
    ! outside of this routine (in fuse_patches).  basically this should be split
    ! into a copy, then change pointers, then delete.

    if(associated(dp%younger)) then 
       dp%younger%older => dp%older
    else 
       dp%siteptr%youngest_patch => dp%older !youngest
    endif
    if(associated(dp%older)) then 
       dp%older%younger => dp%younger
    else 
       dp%siteptr%oldest_patch => dp%younger  !oldest
    endif

    deallocate(dp)

  end subroutine fuse_2_patches

!-------------------------------------------------------------------------------!

  subroutine terminate_patches(cs_pnt)
  ! ============================================================================
  !  Terminate Patches if they  are too small                          
  ! ============================================================================

    use clm_varctl      ,only : iulog    

    implicit none

    type(site), target, intent(in) :: cs_pnt

    type(site),  pointer :: currentSite
    type(patch), pointer :: currentPatch

    real(r8) areatot ! variable for checking whether the total patch area is wrong. 
 
    currentSite => cs_pnt

    currentPatch => currentSite%oldest_patch
    
    !fuse patches if one of them is very small.... 
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch)) 
       if(currentPatch%area <= 0.001_r8)then
          if(associated(currentPatch%older).and.currentPatch%patchno /= currentSite%youngest_patch%patchno)then
            ! Do not force the fusion of the youngest patch to its neighbour. 
            ! This is only really meant for very old patches. 
             write(iulog,*) 'fusing patches because one is too small',currentPatch%area, currentPatch%lai, &
                  currentPatch%older%area,currentPatch%older%lai,currentPatch%seed_bank(1)
             call fuse_2_patches(currentPatch%older,currentPatch)
             deallocate(currentPatch%older)
             write(iulog,*) 'after fusion',currentPatch%area,currentPatch%seed_bank(1)
           endif
       endif

       currentPatch => currentPatch%older

    enddo

    !check area is not exceeded
    areatot = 0._r8
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
       areatot = areatot + currentPatch%area
       currentPatch => currentPatch%younger
       if((areatot-area) > 0.0000001_r8)then
          write(iulog,*) 'ED: areatot too large. end terminate', areatot,currentSite%clmgcell
       endif
    enddo

  end subroutine terminate_patches

!-------------------------------------------------------------------------------!

  subroutine patch_pft_size_profile(cp_pnt)
  ! ============================================================================
  !        Binned patch size profiles generated for patch fusion routine        
  ! ============================================================================
 
    implicit none  

    type(patch), target, intent(inout) :: cp_pnt

    type(patch), pointer  :: currentPatch
    type(cohort), pointer :: currentCohort

    real(r8) :: mind(N_DBH_BINS) ! Bottom of DBH bin 
    real(r8) :: maxd(N_DBH_BINS) ! Top of DBH bin
    real(r8) :: delta_dbh   ! Size of DBH bin
    integer  :: p    ! Counter for PFT 
    integer  :: j    ! Counter for DBH bins 

    currentPatch => cp_pnt

    delta_dbh = (DBHMAX/N_DBH_BINS)

    do p = 1,numpft_ed
       do j = 1,N_DBH_BINS
          currentPatch%pft_agb_profile(p,j) = 0.0_r8
       enddo
    enddo

    do j = 1,N_DBH_BINS   
        if (j == 1) then
           mind(j) = 0.0_r8
           maxd(j) = delta_dbh
        else 
           mind(j) = (j-1) * delta_dbh
           maxd(j) = (j)*delta_dbh
        endif
    enddo

    currentCohort => currentPatch%shortest
    do while(associated(currentCohort))    
       do j = 1,N_DBH_BINS   
          if((currentCohort%dbh  >  mind(j)) .AND. (currentCohort%dbh  <=  maxd(j)))then

             currentPatch%pft_agb_profile(currentCohort%pft,j) = currentPatch%pft_agb_profile(currentCohort%pft,j) + &
                  currentCohort%bdead*currentCohort%n/currentPatch%area

          endif
       enddo ! dbh bins

       ! Deal with largest dbh bin
       j = N_DBH_BINS-1
       if(currentCohort%dbh  >  j*delta_dbh)then

          currentPatch%pft_agb_profile(currentCohort%pft,j) = currentPatch%pft_agb_profile(currentCohort%pft,j) + &
               currentCohort%bdead*currentCohort%n/currentPatch%area

       endif !  

       currentCohort => currentCohort%taller

    enddo !currentCohort 
   
  end subroutine patch_pft_size_profile

!-------------------------------------------------------------------------------!

  function countPatches( bounds ) result ( totNumPatches ) 

  ! ============================================================================
  !        Loop over all Patches to count how many there are
  ! ============================================================================

    use decompMod           , only : bounds_type
    use abortutils          , only : endrun

    implicit none

    type(bounds_type), intent(in) :: bounds

    type (site) , pointer :: currentSite
    type (patch), pointer :: currentPatch
    integer :: g              ! gridcell
    integer :: totNumPatches  ! total number of patches.  

    totNumPatches = 0

    if (allocated(gridCellEdState)) then
       do g = bounds%begg,bounds%endg
          currentSite => gridCellEdState(g)%spnt
          if(currentSite%istheresoil == 1)then
             currentPatch => currentSite%oldest_patch
             do while(associated(currentPatch))
                totNumPatches = totNumPatches + 1
                currentPatch => currentPatch%younger
             enddo
          endif
       enddo
    else
      call endrun(' ED :: countPatches called, but gridCellEdState not yet allocated') 
    endif

   end function countPatches

  ! ============================================================================
end module EDPatchDynamicsMod
