module EDBtranMod
   
   !-------------------------------------------------------------------------------------
   ! Description:
   ! 
   ! ------------------------------------------------------------------------------------
   
   use EDPftvarcon       , only : EDPftvarcon_inst
   use FatesConstantsMod , only : tfrz => t_water_freeze_k_1atm 
   use FatesConstantsMod , only : itrue,ifalse
   use EDTypesMod        , only : ed_site_type,       &
                                  ed_patch_type,      &
                                  ed_cohort_type,     &
                                  maxpft
   use FatesInterfaceMod , only : hlm_numlevgrnd
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use FatesInterfaceMod , only : bc_in_type, &
                                  bc_out_type, &
                                  numpft
   use FatesInterfaceMod , only : hlm_use_planthydro
   use FatesGlobals      , only : fates_log

   !
   implicit none
   private
   
   public :: btran_ed
   public :: get_active_suction_layers
   
contains 
   
  ! ====================================================================================

  logical function check_layer_water(h2o_liq_vol, tempk)
    
    implicit none
    ! Arguments
    real(r8),intent(in) :: h2o_liq_vol
    real(r8),intent(in) :: tempk
    
    check_layer_water = .false.

    if ( h2o_liq_vol .gt. 0._r8 ) then
       if ( tempk .gt. tfrz-2._r8) then
          check_layer_water = .true.
       end if
    end if
    return
  end function check_layer_water

  ! =====================================================================================
  
  subroutine get_active_suction_layers(nsites, sites, bc_in, bc_out)
    
    ! Arguments
    
    integer,intent(in)                      :: nsites
    type(ed_site_type),intent(inout),target :: sites(nsites)
    type(bc_in_type),intent(in)             :: bc_in(nsites)
    type(bc_out_type),intent(inout)         :: bc_out(nsites)
    
    ! !LOCAL VARIABLES:
    integer  :: s                 ! site
    integer  :: j                 ! soil layer
    !------------------------------------------------------------------------------
    
      do s = 1,nsites
         if (bc_in(s)%filter_btran) then
            do j = 1,hlm_numlevgrnd
               bc_out(s)%active_suction_gl(j) = check_layer_water( bc_in(s)%h2o_liqvol_gl(j),bc_in(s)%tempk_gl(j) )
            end do
         else
            bc_out(s)%active_suction_gl(:) = .false.
         end if
      end do

  end subroutine get_active_suction_layers
  
  ! =====================================================================================

  subroutine btran_ed( nsites, sites, bc_in, bc_out)

    use FatesPlantHydraulicsMod, only : BTranForHLMDiagnosticsFromCohortHydr

      
      ! ---------------------------------------------------------------------------------
      ! Calculate the transpiration wetness function (BTRAN) and the root uptake
      ! distribution (ROOTR).
      ! Boundary conditions in: bc_in(s)%eff_porosity_gl(j)    unfrozen porosity
      !                         bc_in(s)%watsat_gl(j)          porosity
      !                         bc_in(s)%active_uptake_gl(j)   frozen/not frozen
      !                         bc_in(s)%smp_gl(j)             suction
      ! Boundary conditions out: bc_out(s)%rootr_pagl          root uptake distribution
      !                          bc_out(s)%btran_pa            wetness factor
      ! ---------------------------------------------------------------------------------
      
      ! Arguments
      
      integer,intent(in)                      :: nsites
      type(ed_site_type),intent(inout),target :: sites(nsites)
      type(bc_in_type),intent(in)             :: bc_in(nsites)
      type(bc_out_type),intent(inout)         :: bc_out(nsites)
      
      !
      ! !LOCAL VARIABLES:
      type(ed_patch_type),pointer             :: cpatch ! Current Patch Pointer
      type(ed_cohort_type),pointer            :: ccohort ! Current cohort pointer
      integer  :: s                 ! site
      integer  :: j                 ! soil layer
      integer  :: ifp               ! patch vector index for the site
      integer  :: ft                ! plant functional type index
      real(r8) :: smp_node          ! matrix potential
      real(r8) :: rresis            ! suction limitation to transpiration independent
                                    ! of root density
      real(r8) :: pftgs(maxpft)     ! pft weighted stomatal conductance s/m
      real(r8) :: temprootr              
      real(r8) :: balive_patch
      real(r8) :: sum_pftgs         ! sum of weighted conductances (for normalization)
      !------------------------------------------------------------------------------
      
      associate(                                 &
            smpsc     => EDPftvarcon_inst%smpsc          , &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
            smpso     => EDPftvarcon_inst%smpso            &  ! INTERF-TODO: THESE SHOULD BE FATES PARAMETERS
            )
        
        do s = 1,nsites

           bc_out(s)%rootr_pagl(:,:) = 0._r8

           ifp = 0
           cpatch => sites(s)%oldest_patch
           do while (associated(cpatch))                 
              ifp=ifp+1
              
              ! THIS SHOULD REALLY BE A COHORT LOOP ONCE WE HAVE rootfr_ft FOR COHORTS (RGK)
              
              do ft = 1,numpft
                 cpatch%btran_ft(ft) = 0.0_r8
                 do j = 1,hlm_numlevgrnd
                    
                    ! Calculations are only relevant where liquid water exists
                    ! see clm_fates%wrap_btran for calculation with CLM/ALM
                    
                    if ( check_layer_water(bc_in(s)%h2o_liqvol_gl(j),bc_in(s)%tempk_gl(j)) )  then
                       
                       smp_node = max(smpsc(ft), bc_in(s)%smp_gl(j))
                       
                       rresis  = min( (bc_in(s)%eff_porosity_gl(j)/bc_in(s)%watsat_gl(j))*               &
                            (smp_node - smpsc(ft)) / (smpso(ft) - smpsc(ft)), 1._r8)
                       
                       cpatch%rootr_ft(ft,j) = cpatch%rootfr_ft(ft,j)*rresis
                       
                       ! root water uptake is not linearly proportional to root density,
                       ! to allow proper deep root funciton. Replace with equations from SPA/Newman. FIX(RF,032414)
                       ! cpatch%rootr_ft(ft,j) = cpatch%rootfr_ft(ft,j)**0.3*rresis_ft(ft,j)/ &
                       ! sum(cpatch%rootfr_ft(ft,1:nlevgrnd)**0.3)
                       cpatch%btran_ft(ft) = cpatch%btran_ft(ft) + cpatch%rootr_ft(ft,j)
                       
                    else
                       cpatch%rootr_ft(ft,j) = 0._r8
                    end if
                    
                 end do !j
                 
                 ! Normalize root resistances to get layer contribution to ET
                 do j = 1,hlm_numlevgrnd    
                    if (cpatch%btran_ft(ft)  >  0.0_r8) then
                       cpatch%rootr_ft(ft,j) = cpatch%rootr_ft(ft,j)/cpatch%btran_ft(ft)
                    else
                       cpatch%rootr_ft(ft,j) = 0._r8
                    end if
                 end do
                 
              end do !PFT
              
              ! PFT-averaged point level root fraction for extraction purposese.
              ! This probably needs to be weighted by actual transpiration from each pft. FIX(RF,032414).
              pftgs(:) = 0._r8
              ccohort => cpatch%tallest
              do while(associated(ccohort))
                 pftgs(ccohort%pft) = pftgs(ccohort%pft) + ccohort%gscan * ccohort%n    
                 ccohort => ccohort%shorter
              enddo
              
              ! Process the boundary output, this is necessary for calculating the soil-moisture
              ! sink term across the different layers in driver/host.  Photosynthesis will
              ! pass the host a total transpiration for the patch.  This needs rootr to be
              ! distributed over the soil layers.
              sum_pftgs = sum(pftgs(1:numpft))

              do j = 1,hlm_numlevgrnd
                 bc_out(s)%rootr_pagl(ifp,j) = 0._r8
                 do ft = 1,numpft
                    if( sum_pftgs > 0._r8)then !prevent problem with the first timestep - might fail
                       !bit-retart test as a result? FIX(RF,032414)  
                       bc_out(s)%rootr_pagl(ifp,j) = bc_out(s)%rootr_pagl(ifp,j) + &
                            cpatch%rootr_ft(ft,j) * pftgs(ft)/sum_pftgs
                    else
                       bc_out(s)%rootr_pagl(ifp,j) = bc_out(s)%rootr_pagl(ifp,j) + &
                            cpatch%rootr_ft(ft,j) * 1./numpft
                    end if
                 enddo
              enddo
              
              ! Calculate the BTRAN that is passed back to the HLM
              ! used only for diagnostics. If plant hydraulics is turned off
              ! we are using the patchxpft level btran calculation
              
              if(hlm_use_planthydro.eq.ifalse) then
                 !weight patch level output BTRAN for the
                 bc_out(s)%btran_pa(ifp) = 0.0_r8
                 do ft = 1,numpft
                    if( sum_pftgs > 0._r8)then !prevent problem with the first timestep - might fail
                       !bit-retart test as a result? FIX(RF,032414)   
                       bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft)  * pftgs(ft)/sum_pftgs
                    else
                       bc_out(s)%btran_pa(ifp)   = bc_out(s)%btran_pa(ifp) + cpatch%btran_ft(ft) * 1./numpft
                    end if
                 enddo
              end if

              temprootr = sum(bc_out(s)%rootr_pagl(ifp,1:hlm_numlevgrnd))

              if(abs(1.0_r8-temprootr) > 1.0e-10_r8 .and. temprootr > 1.0e-10_r8)then
                 write(fates_log(),*) 'error with rootr in canopy fluxes',temprootr,sum_pftgs
                 do j = 1,hlm_numlevgrnd
                    bc_out(s)%rootr_pagl(ifp,j) = bc_out(s)%rootr_pagl(ifp,j)/temprootr
                 enddo
              end if
              
              cpatch => cpatch%younger
           end do
        
        end do
           
        if(hlm_use_planthydro.eq.itrue) then
           call BTranForHLMDiagnosticsFromCohortHydr(nsites,sites,bc_out)
        end if
        
      end associate
      
    end subroutine btran_ed


end module EDBtranMod
