
module EDCanopyStructureMod

  ! ============================================================================
  ! Code to determine whether the canopy is closed, and which plants are either in the understorey or overstorey
  ! This is obviosuly far too complicated for it's own good and needs re-writing.  
  ! ============================================================================

  use shr_kind_mod          , only : r8 => shr_kind_r8;
  use clm_varpar            , only : nclmax
  use clm_varctl            , only : iulog
  use EcophysConType        , only : ecophyscon

  use EDGrowthFunctionsMod  , only : c_area
  use EDCohortDynamicsMod   , only : copy_cohort, terminate_cohorts, fuse_cohorts
  use EDtypesMod            , only : site, patch, cohort, ncwd

  implicit none
  save
  private

  public :: canopy_structure
  public :: canopy_spread

  ! ============================================================================
  ! 10/30/09: Created by Rosie Fisher
  ! ============================================================================

contains

  subroutine canopy_structure( currentSite )
  ! ============================================================================
  ! This routine allocates the 'canopy_layer' attribute to each cohort
  ! All top leaves in the same canopy layer get the same light resources.
  ! The first canopy layer is the 'canopy' or 'overstorey'. The second is the 'understorey'.
  ! More than two layers is not permitted at the moment
  ! Seeds germinating into the 3rd or higher layers are automatically removed. 
  
  ! ------Perfect Plasticity-----
  ! The idea of these canopy layers derives originally from Purves et al. 2009
  ! Their concept is that, given enoughplasticity in canopy position, size, shape and depth
  ! all of the gound area will be filled perfectly by leaves, and additional leaves will have
  ! to exist in the understorey. 
  ! Purves et al. use the concept of 'Z*' to assume that the height required to attain a place in the
  ! canopy is spatially uniform. In this implementation, described in Fisher et al. (2010, New Phyt) we
  ! extent that concept to assume that position in the canopy has some random element, and that BOTH height
  ! and chance combine to determine whether trees get into the canopy. 
  ! Thus, when the canopy is closed and there is excess area, some of it must be demoted
  ! If we demote -all- the trees less than a given height, there is a massive advantage in being the cohort that is 
  ! the biggest when the canopy is closed. 
  ! In this implementation, the amount demoted, ('weight') is a function of the height weighted by the competitive exclusion
  ! parameter (ED_val_comp_excln). 
  
  ! Complexity in this routine results from a few things. 
  ! Firstly, the complication of the demotion amount sometimes being larger than the cohort area (for a very small, short cohort)
  ! Second, occasionaly, disturbance (specifically fire) can cause the canopy layer to become less than closed, 
  ! without changing the area of the patch. If this happens, then some of the plants in the lower layer need to be 'promoted' so 
  ! all of the routine has to happen in both the downwards and upwards directions. 
  
  
  ! The order of events here is therefore:
  ! (The entire subroutine has a single outer 'patch' loop. 
  ! Section 1: figure out the total area, and whether there are >1 canopy layers at all. 
  
  
  ! Sorts out cohorts into canopy and understorey layers...                              
  ! ============================================================================

    use clm_varpar,  only : nlevcan_ed

    use EDParamsMod, only : ED_val_comp_excln, ED_val_ag_biomass
    use SFParamsMod, only : SF_val_cwd_frac

    implicit none 

    type(site), intent(inout), pointer   :: currentSite

    type(patch), pointer  :: currentPatch
    type(cohort), pointer :: currentCohort,copyc

    integer :: i,j
    integer :: z     ! Current number of canopy layers. (1= canopy, 2 = understorey) 
    
    real(r8) :: checkarea
    real(r8) :: cc_loss
    real(r8) :: lossarea
    real(r8) :: newarea
    real(r8) :: arealayer(nlevcan_ed) ! Amount of plant area currently in each canopy layer
    real(r8) :: sumdiff(nlevcan_ed)   ! The total of the exclusion weights for all cohorts in layer z 
    real(r8) :: weight                ! The amount of the total lost area that comes from this cohort
    real(r8) :: sum_weights(nlevcan_ed)
    real(r8) :: new_total_area_check
    real(r8) :: missing_area, promarea,cc_gain,sumgain
    integer  :: promswitch,lower_cohort_switch
    integer  :: c
    real(r8) :: sumloss,excess_area
    integer  :: count_mi

    currentPatch => currentSite%oldest_patch
    
! Section 1: Check  total canopy area.    
    new_total_area_check = 0._r8
    do while (associated(currentPatch)) ! Patch loop    
       excess_area = 1.0_r8   
        
       ! Does any layer have excess area in it? Keep going until it does not...
       
       do while(excess_area > 0.000001_r8)
          
          ! Calculate the area currently in each canopy layer. 
          z = 1 
          arealayer = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))  
             currentCohort%c_area = c_area(currentCohort) ! Reassess cohort area. 
             arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer) + currentCohort%c_area
             z = max(z,currentCohort%canopy_layer) ! What is the current number of canopy layers? 
             currentCohort => currentCohort%shorter
          enddo
          
          ! Does the bottom layer have more than a full canopy? If so we need to make another layer.
          
          if(arealayer(z) > currentPatch%area)then  ! Do we have too much area in either layer?  
              !write(iulog,*) 'CANOPY CLOSURE', z
             z = z + 1
          endif

          currentPatch%NCL_p = min(nclmax,z)   ! Set current canopy layer occupancy indicator.  

          do i = 1,z ! Loop around the currently occupied canopy layers. 
             
             do while((arealayer(i)-currentPatch%area) > 0.000001_r8) 
             ! Is this layer currently over-occupied? 
             ! In that case, we need to work out which cohorts to demote. 
             
                sumloss = 0.0_r8
                new_total_area_check = 0.0_r8
                sumdiff(i) = 0.0_r8    
                currentCohort => currentPatch%tallest 
                do while (associated(currentCohort))
                   currentCohort%c_area = c_area(currentCohort)
                   if(arealayer(i) > currentPatch%area.and.currentCohort%canopy_layer == i)then
                      currentCohort%excl_weight = 1.0_r8/(currentCohort%dbh**ED_val_comp_excln)  
                      sumdiff(i) = sumdiff(i) + currentCohort%excl_weight
                   endif
                   currentCohort => currentCohort%shorter  
                enddo !currentCohort

                lossarea = arealayer(i) - currentPatch%area  !how much do we have to lose?
                sum_weights(i) = 0.0_r8
                currentCohort => currentPatch%tallest    !start from the tallest cohort
                 
                ! Correct the demoted cohorts for  
                do while (associated(currentCohort))
                   if(currentCohort%canopy_layer  ==  i) then
                      weight = currentCohort%excl_weight/sumdiff(i)     
                      currentCohort%excl_weight = min(currentCohort%c_area/lossarea, weight)
                      sum_weights(i) = sum_weights(i) + currentCohort%excl_weight
                   endif
                   currentCohort => currentCohort%shorter      
                enddo

                currentCohort => currentPatch%tallest
                do while (associated(currentCohort))      
                   if(currentCohort%canopy_layer == i)then !All the trees in this layer need to lose some area...
                      weight = currentCohort%excl_weight/sum_weights(i)
                      cc_loss = lossarea*weight !what this cohort has to lose. 
                      !-----------Split and copy boundary cohort-----------------!
                      if(cc_loss < currentCohort%c_area)then
                         allocate(copyc)

                         call copy_cohort(currentCohort,copyc) !makes an identical copy...
                         ! n.b this needs to happen BEFORE the cohort goes into the new layer, 
                         ! otherwise currentPatch%spread(i+1) will be higher and the area will change...!!! 
                         sumloss = sumloss + cc_loss 

                         newarea = currentCohort%c_area - cc_loss
                         copyc%n = currentCohort%n*newarea/currentCohort%c_area   !
                         currentCohort%n = currentCohort%n - (currentCohort%n*newarea/currentCohort%c_area) !     

                         copyc%canopy_layer = i !the taller cohort is the copy
                         currentCohort%canopy_layer = i + 1 !demote the current cohort to the understory.           
                         ! seperate cohorts. 
                         ! - 0.000000000001_r8 !needs to be a very small number to avoid 
                         ! causing non-linearity issues with c_area.  is this really required? 
                         currentCohort%dbh = currentCohort%dbh 
                         copyc%dbh = copyc%dbh !+ 0.000000000001_r8
                         !kill the ones which go into canopy layers that are not allowed... (default nclmax=2) 
                         if(i+1 > nclmax)then 
                           !put the litter from the terminated cohorts into the fragmenting pools
                          ! write(iulog,*) '3rd canopy layer'
                            do c=1,ncwd

                               currentPatch%CWD_AG(c)  = currentPatch%CWD_AG(c) + (currentCohort%bdead+currentCohort%bsw) * &
                                    ED_val_ag_biomass * &
                                    SF_val_CWD_frac(c)*currentCohort%n/currentPatch%area  
         
                               currentPatch%CWD_BG(c)  = currentPatch%CWD_BG(c) + (currentCohort%bdead+currentCohort%bsw) * &
                                    (1.0_r8-ED_val_ag_biomass) * &
                                    SF_val_CWD_frac(c)*currentCohort%n/currentPatch%area !litter flux per m2.

                            enddo

                                currentPatch%leaf_litter(currentCohort%pft)  = &
                                     currentPatch%leaf_litter(currentCohort%pft) + (currentCohort%bl)* &
                                          currentCohort%n/currentPatch%area ! leaf litter flux per m2.

                                currentPatch%root_litter(currentCohort%pft)  = &
                                     currentPatch%root_litter(currentCohort%pft) + &
                                     (currentCohort%br+currentCohort%bstore)*currentCohort%n/currentPatch%area
   									
                            currentCohort%n = 0.0_r8
                            currentCohort%c_area = 0._r8
                         else  
                            currentCohort%c_area = c_area(currentCohort)       
                         endif
                         copyc%c_area = c_area(copyc)
                         new_total_area_check = new_total_area_check+copyc%c_area

                         !----------- Insert copy into linked list ------------------------!                         
                         copyc%shorter => currentCohort
                         if(associated(currentCohort%taller))then
                            copyc%taller => currentCohort%taller
                            currentCohort%taller%shorter => copyc
                         else
                            currentPatch%tallest => copyc
                            copyc%taller => null()
                         endif
                         currentCohort%taller => copyc                  
                      else
                         currentCohort%canopy_layer = i + 1 !the whole cohort becomes demoted
                         sumloss = sumloss + currentCohort%c_area 

                         !kill the ones which go into canopy layers that are not allowed... (default nclmax=2) 
                         if(i+1 > nclmax)then  
                           !put the litter from the terminated cohorts into the fragmenting pools
                            do c=1,ncwd

                               currentPatch%CWD_AG(c)  = currentPatch%CWD_AG(c) + (currentCohort%bdead+currentCohort%bsw) * &
                                    ED_val_ag_biomass * &
                                    SF_val_CWD_frac(c)*currentCohort%n/currentPatch%area           
                               currentPatch%CWD_BG(c)  = currentPatch%CWD_BG(c) + (currentCohort%bdead+currentCohort%bsw) * &
                                    (1.0_r8-ED_val_ag_biomass) * &
                                    SF_val_CWD_frac(c)*currentCohort%n/currentPatch%area !litter flux per m2.

                            enddo

                                currentPatch%leaf_litter(currentCohort%pft)  = &
                                     currentPatch%leaf_litter(currentCohort%pft) + currentCohort%bl* &
                                          currentCohort%n/currentPatch%area ! leaf litter flux per m2.

                                currentPatch%root_litter(currentCohort%pft)  = &
                                     currentPatch%root_litter(currentCohort%pft) + &
                                     (currentCohort%br+currentCohort%bstore)*currentCohort%n/currentPatch%area
                            currentCohort%n = 0.0_r8
                            currentCohort%c_area = 0._r8

                         else  
                            currentCohort%c_area = c_area(currentCohort)       
                         endif

                         !write(iulog,*) 'demoting whole cohort', currentCohort%c_area,cc_loss, &
                              !currentCohort%canopy_layer,currentCohort%dbh

                      endif
                     ! call terminate_cohorts(currentPatch) 

                      !----------- End of cohort splitting ------------------------------!             
                   endif !canopy layer = i

                   currentCohort => currentCohort%shorter

                enddo !currentCohort 
                
                call terminate_cohorts(currentPatch)
                arealayer(i) = arealayer(i) - sumloss
                !Update arealayer for diff calculations of layer below. 
                arealayer(i + 1) = arealayer(i + 1) + sumloss 

             enddo !arealayer loop
             if(arealayer(i)-currentPatch%area > 0.00001_r8)then
                write(iulog,*) 'lossarea problem', lossarea,sumloss,z,currentPatch%patchno,currentPatch%clm_pno
             endif

          enddo !z  

          z = 1
          arealayer = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))  
             currentCohort%c_area = c_area(currentCohort)
             arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer) + currentCohort%c_area
             z = max(z,currentCohort%canopy_layer)
             currentCohort => currentCohort%shorter
          enddo

          !does the bottom layer have more than a full canopy? If so we need to make another layer.
          if(arealayer(z) > currentPatch%area)then
             z = z + 1
          endif
          excess_area = 0.0_r8
          do j=1,z
             if(arealayer(j) > currentPatch%area)then 
                excess_area = arealayer(j)-currentPatch%area
             endif
          enddo
          currentPatch%ncl_p = min(z,nclmax)

       enddo !is there still excess area in any layer?      

       call terminate_cohorts(currentPatch)
       call fuse_cohorts(currentPatch)
       call terminate_cohorts(currentPatch)

       ! ----------- Check cohort area ------------------------------!
       do i = 1,z
          checkarea = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))
             if(currentCohort%canopy_layer == i)then
                checkarea = checkarea + c_area(currentCohort)
             endif

             currentCohort=>currentCohort%shorter

          enddo

       enddo ! 


       ! ----------- Check whether the intended 'full' layers are actually filling all the space.
       ! If not, promote some fraction of cohorts upwards ------------------------------!     
       ! THIS SECTION MIGHT BE TRIGGERED BY A FIRE OR MORTALITY EVENT, FOLLOWED BY A PATCH FUSION, 
       ! SO THE TOP LAYER IS NO LONGER FULL...  

       promswitch = 0 

       missing_area=1.0_r8    
       count_mi = 0
       !does any layer have excess area in it? keep going until it does not...
       do while(missing_area > 0.000001_r8.and.z > 1) 
          count_mi = count_mi +1
          do i = 1,z-1 ! if z is greater than one, there is a possibility of too many plants in the understorey. 
             lower_cohort_switch = 1
             ! is the area of the layer less than the area of the patch, if it is supposed to be closed (z>1) 
             do while((arealayer(i)-currentPatch%area) < -0.000001_r8.and.lower_cohort_switch == 1) 

                if(arealayer(i+1) <= 0.000001_r8)then
                   currentCohort => currentPatch%tallest 
                   arealayer = 0._r8
                   do while (associated(currentCohort))            
                      if(currentCohort%canopy_layer == i+1)then !look at the cohorts in the canopy layer below... 
                         currentCohort%canopy_layer = i   
                         currentCohort%c_area = c_area(currentCohort)

                        ! write(iulog,*) 'promoting very small cohort', currentCohort%c_area,currentCohort%canopy_layer
                      endif
                      arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer)+currentCohort%c_area 
                      currentCohort => currentCohort%shorter   
                   enddo

                endif !promoting all of the small amount of area in the lower layers. 


                lower_cohort_switch = 0 
                sumgain = 0.0_r8
                sumdiff(i) = 0.0_r8    
                ! figure out with what weighting we need to promote cohorts.
                ! This is the opposite of the demotion weighting... 
                currentCohort => currentPatch%tallest 
                do while (associated(currentCohort))
                   currentCohort%c_area = c_area(currentCohort)
                   if(currentCohort%canopy_layer == i+1)then !look at the cohorts in the canopy layer below... 
                      currentCohort%prom_weight = currentCohort%dbh**ED_val_comp_excln   !as opposed to 1/(dbh^C_e) 
                      sumdiff(i) = sumdiff(i) + currentCohort%prom_weight
                   endif
                   currentCohort => currentCohort%shorter  
                enddo !currentCohort

                promarea =  currentPatch%area -arealayer(i) !how much do we need to gain?
                sum_weights(i) = 0.0_r8
                currentCohort => currentPatch%tallest    !start from the tallest cohort

                do while (associated(currentCohort))
                   if(currentCohort%canopy_layer  ==  i+1) then !still looking at the layer beneath. 
                      weight = currentCohort%prom_weight/sumdiff(i)
                      if(promarea > 0._r8)then    
                         currentCohort%prom_weight = min(currentCohort%c_area/promarea, weight)
                      else
                         currentCohort%prom_weight = 0._r8
                      endif
                      sum_weights(i) = sum_weights(i) + currentCohort%prom_weight
                   endif
                   currentCohort => currentCohort%shorter      
                enddo

                currentCohort => currentPatch%tallest
                do while (associated(currentCohort))      
                   if(currentCohort%canopy_layer == i+1)then !All the trees in this layer need to promote some area upwards... 
                      lower_cohort_switch = 1
                      weight = currentCohort%prom_weight/sum_weights(i)
                      cc_gain = promarea*weight !what this cohort has to promote. 
                      !-----------Split and copy boundary cohort-----------------!
                      if(cc_gain < currentCohort%c_area)then
                         allocate(copyc)

                         call copy_cohort(currentCohort,copyc) !makes an identical copy...
                         ! n.b this needs to happen BEFORE the cohort goes into the new layer, otherwise currentPatch
                         ! %spread(+1) will be higher and the area will change...!!!
                         sumgain = sumgain + cc_gain


                         newarea = currentCohort%c_area - cc_gain !new area of existing cohort
                         copyc%n = currentCohort%n*cc_gain/currentCohort%c_area   !number of individuals in promoted cohort. 
                         ! number of individuals in cohort remianing in understorey    
                         currentCohort%n = currentCohort%n - (currentCohort%n*cc_gain/currentCohort%c_area) 

                         currentCohort%canopy_layer = i+1  !keep current cohort in the understory.        
                         copyc%canopy_layer = i ! promote copy to the higher canopy layer. 

                         ! seperate cohorts. 
                         ! needs to be a very small number to avoid causing non-linearity issues with c_area. 
                         ! is this really required? 
                         currentCohort%dbh = currentCohort%dbh - 0.000000000001_r8 
                         copyc%dbh = copyc%dbh + 0.000000000001_r8

                         currentCohort%c_area = c_area(currentCohort)          
                         copyc%c_area = c_area(copyc)

                         !----------- Insert copy into linked list ------------------------!                         
                         copyc%shorter => currentCohort
                         if(associated(currentCohort%taller))then
                            copyc%taller => currentCohort%taller
                            currentCohort%taller%shorter => copyc
                         else
                            currentPatch%tallest => copyc
                            copyc%taller => null()
                         endif
                         currentCohort%taller => copyc                  
                      else
                         currentCohort%canopy_layer = i  !the whole cohort becomes promoted
                         sumgain = sumgain + currentCohort%c_area  !inserting deliberate mistake to see how far we make it... 
                         ! update area AFTER we sum up the losses. the cohort may shrink at this point,
                         ! if the upper canopy spread is smaller. this shold be dealt with by the 'excess area' loop.  
                         currentCohort%c_area = c_area(currentCohort) 

                         promswitch = 1

                        ! write(iulog,*) 'promoting whole cohort', currentCohort%c_area,cc_gain,currentCohort%canopy_layer, &
                              !currentCohort%pft,currentPatch%patchno

                      endif
                      !call terminate_cohorts(currentPatch) 
                      if(promswitch == 1)then
                        ! write(iulog,*) 'cohort loop',currentCohort%pft,currentCohort%indexnumber,currentPatch%patchno
                      endif
                      !----------- End of cohort splitting ------------------------------!             
                   else
                      if(promswitch == 1)then
                         ! write(iulog,*) 'cohort list',currentCohort%pft,currentCohort%indexnumber, &
                             ! currentCohort%canopy_layer,currentCohort%c_area
                      endif
                   endif

                   currentCohort => currentCohort%shorter
                enddo !currentCohort 
                arealayer(i) = arealayer(i) + sumgain
                arealayer(i + 1) = arealayer(i + 1) - sumgain !Update arealayer for diff calculations of layer below. 

                if(promswitch == 1)then
                  ! write(iulog,*) 'arealayer loop',arealayer(1:3),currentPatch%area,promarea,sumgain, &
                        !currentPatch%patchno,z,i,lower_cohort_switch
                endif
                if(promswitch == 1.and.associated(currentPatch%tallest))then
                   ! write(iulog,*) 'cohorts',currentCohort%pft,currentCohort%indexnumber,currentPatch%patchno, &
                        !currentCohort%c_area
                endif
             enddo !arealayer loop

             if(currentPatch%area-arealayer(i) < 0.000001_r8)then
                !write(iulog,*) 'gainarea problem',sumgain,arealayer(i),currentPatch%area,z, &
                     !currentPatch%patchno,currentPatch%clm_pno,currentPatch%area - arealayer(i),i,missing_area,count_mi
             endif
             if(promswitch == 1)then
               ! write(iulog,*) 'z loop',arealayer(1:3),currentPatch%patchno,z
             endif
          enddo !z  

          z = 1
          arealayer = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))  
             currentCohort%c_area = c_area(currentCohort)
             arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer) + currentCohort%c_area
             z = max(z,currentCohort%canopy_layer)
             currentCohort => currentCohort%shorter
          enddo

          missing_area = 0.0_r8
          do j=1,z-1
             if(arealayer(j) < currentPatch%area)then  !this is the amount of area that we still have spare in this layer. 
                missing_area = currentPatch%area - arealayer(j)
                if(missing_area <= 0.000001_r8.and.missing_area > 0._r8)then
                   missing_area = 0.0_r8
                  ! write(iulog,*) 'correcting MI',j,currentPatch%area - arealayer(j)
                endif
             endif
          enddo
          currentPatch%ncl_p = min(z,nclmax)
          if(promswitch == 1)then
            ! write(iulog,*) 'missingarea loop',arealayer(1:3),currentPatch%patchno,missing_area,z
          endif
       enddo !is there still not enough canopy area in any layer?         

       call terminate_cohorts(currentPatch)
       call fuse_cohorts(currentPatch)
       call terminate_cohorts(currentPatch)

       if(promswitch == 1)then
          !write(iulog,*) 'going into cohort check',currentPatch%clm_pno
       endif
       ! ----------- Check cohort area ------------------------------!
       do i = 1,z
          checkarea = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))
             if(currentCohort%canopy_layer == i)then
                checkarea = checkarea + c_area(currentCohort)
             endif

             currentCohort => currentCohort%shorter

          enddo

          if(((checkarea-currentPatch%area)) > 0.0001)then
             write(iulog,*) 'problem with canopy area', checkarea,currentPatch%area,checkarea-currentPatch%area,i,z,missing_area 
             currentCohort => currentPatch%tallest
             do while (associated(currentCohort))
             if(currentCohort%canopy_layer == i)then
                write(iulog,*) 'c_areas in top layer', c_area(currentCohort)
             endif
             currentCohort => currentCohort%shorter

          enddo
                  
          endif

          if ( i  >  1) then
             if ( (arealayer(i)  -  arealayer(i-1) )>1e-11 ) then
                write(iulog,*) 'smaller top layer than bottom layer ',arealayer(i),arealayer(i-1), &
                     currentPatch%area,currentPatch%spread(i-1:i)
             endif
          endif
       enddo ! 

       if(promswitch == 1)then 
         ! write(iulog,*) 'end patch loop',currentSite%clmgcell
       endif

       currentPatch => currentPatch%younger
    enddo !patch  

    if(promswitch == 1)then
      ! write(iulog,*) 'end  canopy structure',currentSite%clmgcell
    endif

  end subroutine canopy_structure

  ! ============================================================================
  !  Calculates the spatial spread of tree canopies based on canopy closure.                             
  ! ============================================================================
  subroutine canopy_spread( currentSite )

    use clm_varpar  , only : nlevcan_ed
    use EDParamsMod , only : ED_val_maxspread, ED_val_minspread 

    implicit none

    type (site), intent(inout),   pointer :: currentSite

    type (cohort), pointer :: currentCohort
    type (patch),  pointer :: currentPatch

    real(r8) :: arealayer(nlevcan_ed) ! Amount of canopy in each layer. 
    real(r8) :: inc                   ! Arbitrary daily incremental change in canopy area 
    integer z

    inc = 0.005_r8

    currentPatch => currentSite%oldest_patch

    do while (associated(currentPatch))

       !calculate canopy area in each canopy storey...
       arealayer = 0.0_r8   
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))
          currentCohort%c_area = c_area(currentCohort) 
          if(ecophyscon%woody(currentCohort%pft) == 1)then
             arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer) + currentCohort%c_area
          endif
          currentCohort => currentCohort%shorter
       enddo

       !If the canopy area is approaching closure, squash the tree canopies and make them taller and thinner
       do z = 1,nclmax  
         
          if(arealayer(z)/currentPatch%area > 0.9_r8)then
             currentPatch%spread(z) = currentPatch%spread(z) - inc
          else 
             currentPatch%spread(z) = currentPatch%spread(z) + inc 
          endif
          if(currentPatch%spread(z) >= ED_val_maxspread)then 
             currentPatch%spread(z) = ED_val_maxspread
          endif
          if(currentPatch%spread(z) <=  ED_val_minspread)then
             currentPatch%spread(z) = ED_val_minspread
          endif
        enddo !z
        !write(iulog,*) 'spread',currentPatch%spread(1:2)
        !currentPatch%spread(:) = ED_val_maxspread
        !FIX(RF,033114) spread is off
        !write(iulog,*) 'canopy_spread',currentPatch%area,currentPatch%spread(1:2)
        currentPatch => currentPatch%younger

    enddo !currentPatch

  end subroutine canopy_spread


  ! ============================================================================
end module EDCanopyStructureMod
