module ChecksBalancesMod

   use shr_kind_mod , only : r8 => shr_kind_r8
   use shr_const_mod, only: SHR_CONST_CDAY
   use EDtypesMod   , only : ed_site_type,ed_patch_type,ed_cohort_type
   use EDTypesMod   , only : AREA

   implicit none
   
   private
   public :: SummarizeNetFluxes
   public :: FATES_BGC_Carbon_Balancecheck
   public :: SiteCarbonStock

contains
   
   !------------------------------------------------------------------------
   
   subroutine SummarizeNetFluxes( nsites, sites, bc_in, is_beg_day )
      
      ! Summarize the combined production and decomposition fluxes into net fluxes
      ! This is done on the fast timestep, and to be called after both daily ED calls and 
      ! fast BGC calls.  Does not include summarization of fast-timestsp productivity calls 
      ! because these must be summarized prior to daily ED calls
      !
      ! Written by Charlie Koven, Feb 2016
      !
      ! !USES: 
      use FatesInterfaceMod , only : bc_in_type
     
      use EDtypesMod        , only : AREA
      !
      implicit none   
      !
      ! !ARGUMENTS    
      
      integer                                 , intent(in)    :: nsites
      type(ed_site_type)                      , intent(inout), target :: sites(nsites)
      type(bc_in_type)                        , intent(in)    :: bc_in(nsites)
      logical                                 , intent(in)    :: is_beg_day
      
      !
      ! !LOCAL VARIABLES:
      integer :: s
      
      type (ed_patch_type)  , pointer :: currentPatch
      type (ed_cohort_type) , pointer :: currentCohort
      real(r8) :: n_perm2     ! individuals per m2 of the whole column
      
      do s = 1,nsites
         
         sites(s)%fire_c_to_atm   = 0._r8    ! REMOVE THIS VARIABLE?
         sites(s)%ed_litter_stock = 0._r8
         sites(s)%cwd_stock       = 0._r8
         sites(s)%biomass_stock   = 0._r8
         
         ! map ed site-level fire fluxes to clm column fluxes
         sites(s)%fire_c_to_atm = sites(s)%total_burn_flux_to_atm / &
               ( AREA * SHR_CONST_CDAY * 1.e3_r8)
         
         currentPatch => sites(s)%oldest_patch
         do while(associated(currentPatch))
            
            ! map litter, CWD, and seed pools to column level
            sites(s)%cwd_stock = sites(s)%cwd_stock + &
                  (currentPatch%area / AREA) * &
                  (sum(currentPatch%cwd_ag)+sum(currentPatch%cwd_bg)) * 1.e3_r8
            
            sites(s)%ed_litter_stock = sites(s)%ed_litter_stock + &
                  (currentPatch%area / AREA) * &
                  (sum(currentPatch%leaf_litter)+sum(currentPatch%root_litter)) * 1.e3_r8
            
            currentCohort => currentPatch%tallest
            do while(associated(currentCohort))
               ! for quantities that are natively at column level or higher, 
               ! calculate plant density using whole area (for grid cell averages)
               n_perm2   = currentCohort%n/AREA                    
               
               ! map biomass pools to column level
               sites(s)%biomass_stock = sites(s)%biomass_stock + &
                     (currentCohort%bdead + currentCohort%balive + &
                     currentCohort%bstore) * n_perm2 * 1.e3_r8
               
               currentCohort => currentCohort%shorter
            enddo !currentCohort
            currentPatch => currentPatch%younger
         end do ! patch loop
         
         ! calculate NEP and NBP fluxes
         sites(s)%nep = sites(s)%npp - bc_in(s)%tot_het_resp
         sites(s)%nbp = sites(s)%npp - ( bc_in(s)%tot_het_resp + sites(s)%fire_c_to_atm )
         
         ! FATES stocks
         sites(s)%totfatesc = sites(s)%ed_litter_stock + sites(s)%cwd_stock + &
               (sum(sites(s)%seed_bank) * 1.e3_r8) + sites(s)%biomass_stock
         
         ! BGC stocks (used for error checking, totlitc should be zero?)
         sites(s)%totbgcc = bc_in(s)%tot_somc +  bc_in(s)%tot_litc
         
         ! Total Ecosystem Carbon Stocks
         sites(s)%totecosysc = sites(s)%totfatesc + sites(s)%totbgcc
         
      end do
      
      ! in FATES timesteps, because of offset between when ED and BGC reconcile the gain 
      ! and loss of litterfall carbon, (i.e. FATES reconciles it instantly, while BGC 
      ! reconciles it incrementally over the subsequent day) calculate the total 
      ! ED -> BGC flux and keep track of the last day's info for balance checking purposes
      if ( is_beg_day ) then
         !
         do s = 1,nsites
            ! order of operations in the next to lines is quite important ;)
            sites(s)%fates_to_bgc_last_ts = sites(s)%fates_to_bgc_this_ts
            sites(s)%fates_to_bgc_this_ts = 0._r8
            sites(s)%tot_seed_rain_flux   = 0._r8
            
            currentPatch => sites(s)%oldest_patch
            do while(associated(currentPatch))
               !
               sites(s)%fates_to_bgc_this_ts = sites(s)%fates_to_bgc_this_ts + &
                     (sum(currentPatch%CWD_AG_out) + sum(currentPatch%CWD_BG_out) + &
                     sum(currentPatch%seed_decay) + sum(currentPatch%leaf_litter_out) + &
                     sum(currentPatch%root_litter_out)) * &
                     ( currentPatch%area/AREA ) * 1.e3_r8 / ( 365.0_r8*SHR_CONST_CDAY )
               !
               sites(s)%tot_seed_rain_flux = sites(s)%tot_seed_rain_flux + &
                     sum(sites(s)%seed_rain_flux) * 1.e3_r8 / ( 365.0_r8*SHR_CONST_CDAY )
               !
               currentPatch => currentPatch%younger
            end do !currentPatch
         end do
      endif
      
      return
   end subroutine SummarizeNetFluxes
   
   ! ====================================================================================
   
   subroutine FATES_BGC_Carbon_Balancecheck(nsites, sites, bc_in, is_beg_day, dtime, nstep)
      
      ! Integrate in time the fluxes into and out of the ecosystem, and compare these 
      ! on a daily timestep to the chagne in carbon stocks of the ecosystem
      !
      ! Written by Charlie Koven, Feb 2016
      !
      ! !USES: 
      use FatesInterfaceMod , only : bc_in_type
      use EDtypesMod        , only : ed_site_type
      !
      implicit none   
      !
      ! !ARGUMENTS    
      integer                                 , intent(in)    :: nsites
      type(ed_site_type)                      , intent(inout), target :: sites(nsites)
      type(bc_in_type)                        , intent(in)    :: bc_in(nsites)
      logical                                 , intent(in)    :: is_beg_day
      real(r8)                                , intent(in)    :: dtime  ! time-step length (s)
      integer                                 , intent(in)    :: nstep  ! time-step index
      
      ! !LOCAL VARIABLES:
      real(r8) :: error_tolerance = 1.e-6_r8
      integer  :: s
      
      ! TODO: THIS INITIALIZATION SHOULD BE IN AN INITIALIZATION PART OF THE CODE
      ! COLD-START PORTION, NSTEP IS >1 FOR RESTARTS RIGHT? (RGK)
      if (nstep .le. 1) then
         ! when starting up the model, initialize the integrator variables
         do s = 1,nsites
            sites(s)%totecosysc_old       = sites(s)%totecosysc
            sites(s)%totfatesc_old        = sites(s)%totfatesc
            sites(s)%totbgcc_old          = sites(s)%totbgcc
            sites(s)%nep_timeintegrated   = 0._r8
            sites(s)%hr_timeintegrated    = 0._r8
            sites(s)%npp_timeintegrated   = 0._r8
            !
            ! also initialize the ed-BGC flux variables
            sites(s)%fates_to_bgc_this_ts = 0._r8
            sites(s)%fates_to_bgc_last_ts = 0._r8
            !
            sites(s)%cbal_err_fates = 0._r8
            sites(s)%cbal_err_bgc   = 0._r8
            sites(s)%cbal_err_tot = 0._r8
         end do
      endif
      
      do s = 1,nsites
         sites(s)%nep_timeintegrated = sites(s)%nep_timeintegrated + sites(s)%nep * dtime
         sites(s)%hr_timeintegrated  = sites(s)%hr_timeintegrated  + bc_in(s)%tot_het_resp * dtime
         sites(s)%npp_timeintegrated = sites(s)%npp_timeintegrated + sites(s)%npp * dtime
      end do
      
      ! If this is on the dynamics time-step, then we calculate the balance checks
      
      if ( is_beg_day ) then
         
         do s = 1,nsites
            
            ! NBP can only be updated when dynamics level information is available
            sites(s)%nbp_integrated  = sites(s)%nep_timeintegrated - &
                  sites(s)%fire_c_to_atm * SHR_CONST_CDAY + &
                  sites(s)%tot_seed_rain_flux * SHR_CONST_CDAY

         
            sites(s)%cbal_err_fates = sites(s)%totfatesc - & 
                  sites(s)%totfatesc_old - &
                  (sites(s)%npp_timeintegrated + &
                  sites(s)%tot_seed_rain_flux*SHR_CONST_CDAY - &
                  sites(s)%fates_to_bgc_this_ts*SHR_CONST_CDAY - &
                  sites(s)%fire_c_to_atm*SHR_CONST_CDAY)
            sites(s)%cbal_err_fates = sites(s)%cbal_err_fates / SHR_CONST_CDAY
            
            sites(s)%cbal_err_bgc = sites(s)%totbgcc - &
                  sites(s)%totbgcc_old - &
                  (sites(s)%fates_to_bgc_last_ts*SHR_CONST_CDAY - &
                  sites(s)%hr_timeintegrated)
            sites(s)%cbal_err_bgc = sites(s)%cbal_err_bgc / SHR_CONST_CDAY
            
            sites(s)%cbal_err_tot = sites(s)%totecosysc - &
                  sites(s)%totecosysc_old - &
                  (sites(s)%nbp_integrated + &
                  sites(s)%fates_to_bgc_last_ts*SHR_CONST_CDAY - &
                  sites(s)%fates_to_bgc_this_ts*SHR_CONST_CDAY)
            sites(s)%cbal_err_tot = sites(s)%cbal_err_tot / SHR_CONST_CDAY
            
            ! Send the current to the previous/last
            sites(s)%totecosysc_old = sites(s)%totecosysc
            sites(s)%totfatesc_old  = sites(s)%totfatesc
            sites(s)%totbgcc_old    = sites(s)%totbgcc
            sites(s)%nep_timeintegrated = 0._r8
            sites(s)%npp_timeintegrated = 0._r8
            sites(s)%hr_timeintegrated = 0._r8
            
         end do
         
      endif
      
      return
   end subroutine FATES_BGC_Carbon_Balancecheck
   
   ! ==============================================================================================

   subroutine SiteCarbonStock(currentSite,total_stock,biomass_stock,litter_stock,seed_stock)
     
     type(ed_site_type),intent(inout),target :: currentSite
     real(r8),intent(out)                    :: total_stock
     real(r8),intent(out)                    :: litter_stock
     real(r8),intent(out)                    :: biomass_stock
     real(r8),intent(out)                    :: seed_stock

     type(ed_patch_type), pointer :: currentPatch
     type(ed_cohort_type), pointer :: currentCohort
     
     litter_stock  = 0.0_r8
     biomass_stock = 0.0_r8
     seed_stock   =  sum(currentSite%seed_bank)*AREA

     currentPatch => currentSite%oldest_patch 
     do while(associated(currentPatch))
        litter_stock = litter_stock + currentPatch%area * &
              (sum(currentPatch%cwd_ag)      + &
               sum(currentPatch%cwd_bg)      + &
               sum(currentPatch%leaf_litter) + &
               sum(currentPatch%root_litter))

        currentCohort => currentPatch%tallest
        do while(associated(currentCohort))
           biomass_stock =  biomass_stock + (currentCohort%bdead + currentCohort%balive + &
                 currentCohort%bstore) * currentCohort%n
           currentCohort => currentCohort%shorter
        enddo !end cohort loop 
        currentPatch => currentPatch%younger
     enddo !end patch loop
     
     total_stock = biomass_stock + seed_stock + litter_stock

     return
  end subroutine SiteCarbonStock



   
end module ChecksBalancesMod
