module CNFireMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNFireMod
!
! !DESCRIPTION:
! Module holding routines fire mod
! nitrogen code.
!
! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_PI,SHR_CONST_TKFRZ
  use pft2colMod   , only: p2c
  use clm_varctl   , only: iulog, use_cn, use_cndv
  implicit none
  save
  private
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNFireArea
  public :: CNFireFluxes
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNFireArea
!
! !INTERFACE:
subroutine CNFireArea (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! Computes column-level area affected by fire in each timestep
! based on statistical fire model in Thonicke et al. 2001.
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size, get_nstep, get_days_per_year
   use clm_varpar      , only: max_pft_per_col
   use clm_varcon      , only: secspday
   use clm_varctl      , only: use_nofire
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! pft-level
   real(r8), pointer :: wtcol(:)        ! pft weight on the column
   integer , pointer :: ivt(:)          ! vegetation type for this pft
   real(r8), pointer :: woody(:)        ! binary flag for woody lifeform (1=woody, 0=not woody)
   ! column-level
   integer , pointer :: npfts(:)        ! number of pfts on the column
   integer , pointer :: pfti(:)         ! pft index array
   real(r8), pointer :: pwtgcell(:)     ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: wf(:)           ! soil water as frac. of whc for top 0.5 m
   real(r8), pointer :: t_grnd(:)       ! ground temperature (Kelvin)
   real(r8), pointer :: totlitc(:)      ! (gC/m2) total litter C (not including cwdc)
   real(r8), pointer :: cwdc(:)         ! (gC/m2) coarse woody debris C
   ! PET 5/20/08, test to increase fire area
   real(r8), pointer :: totvegc(:)    ! (gC/m2) total veg C (column-level mean)
   ! pointers for column averaging
!
! local pointers to implicit in/out scalars
!
   ! column-level
   real(r8), pointer :: me(:)               ! column-level moisture of extinction (proportion)
   real(r8), pointer :: fire_prob(:)        ! daily fire probability (0-1)
   real(r8), pointer :: mean_fire_prob(:)   ! e-folding mean of daily fire probability (0-1)
   real(r8), pointer :: fireseasonl(:)      ! annual fire season length (days, <= days/year)
   real(r8), pointer :: farea_burned(:)     ! fractional area burned in this timestep (proportion)
   real(r8), pointer :: ann_farea_burned(:) ! annual total fractional area burned (proportion)
!
! !OTHER LOCAL VARIABLES:
!   real(r8), parameter:: minfuel = 200.0_r8 ! dead fuel threshold to carry a fire (gC/m2)
! PET, 5/30/08: changed from 200 to 100 gC/m2, since the original paper didn't specify
! the units as carbon, I am assuming that they were in dry biomass, so carbon would be ~50%
   real(r8), parameter:: minfuel = 100.0_r8 ! dead fuel threshold to carry a fire (gC/m2)
   real(r8), parameter:: me_woody = 0.3_r8  ! moisture of extinction for woody PFTs (proportion)
   real(r8), parameter:: me_herb  = 0.2_r8  ! moisture of extinction for herbaceous PFTs (proportion)
   real(r8), parameter:: ef_time = 1.0_r8   ! e-folding time constant (years)
   integer :: fc,c,pi,p ! index variables
   real(r8):: dt        ! time step variable (s)
   real(r8):: fuelc     ! temporary column-level litter + cwd C (gC/m2)
   integer :: nef       ! number of e-folding timesteps
   real(r8):: ef_nsteps ! number of e-folding timesteps (real)
   integer :: nstep     ! current timestep number
   real(r8):: m         ! top-layer soil moisture (proportion)
   real(r8):: mep       ! pft-level moisture of extinction [proportion]
   real(r8):: s2        ! (mean_fire_prob - 1.0)
   real(r8):: dayspyr   ! days per year
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type members (pft-level)
   wtcol            => pft%wtcol
   ivt              => pft%itype
   pwtgcell         => pft%wtgcell  
   woody            => pftcon%woody

   ! assign local pointers to derived type members (column-level)
   npfts            => col%npfts
   pfti             => col%pfti
   wf               => cps%wf
   me               => cps%me
   fire_prob        => cps%fire_prob
   mean_fire_prob   => cps%mean_fire_prob
   fireseasonl      => cps%fireseasonl
   farea_burned     => cps%farea_burned
   ann_farea_burned => cps%ann_farea_burned
   t_grnd           => ces%t_grnd
   totlitc          => ccs%totlitc
   cwdc             => ccs%cwdc
   ! PET 5/20/08, test to increase fire area
   totvegc          => pcs_a%totvegc

   ! pft to column average for moisture of extinction
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      me(c) = 0._r8
   end do
   mep = me_woody
   do pi = 1,max_pft_per_col
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1
            if (pwtgcell(p)>0._r8) then
               if (woody(ivt(p)) == 1) then
                  mep = me_woody
               else
                  mep = me_herb
               end if
            end if
            me(c) = me(c) + mep*wtcol(p)
         end if
      end do
   end do

   ! Get model step size
   dt      = real( get_step_size(), r8 )

   ! Set the number of timesteps for e-folding.
   ! When the simulation has run fewer than this number of steps,
   ! re-scale the e-folding time to get a stable early estimate.
   nstep   = get_nstep()
   dayspyr = get_days_per_year()
   nef = (ef_time*dayspyr*secspday)/dt
   ef_nsteps = max(1,min(nstep,nef))
   
   ! test code, added 6/6/05, PET
   ! setting ef_nsteps to full count regardless of nstep, to see if this
   ! gets rid of transient in fire stats for initial run from spunup 
   ! initial conditions
   ef_nsteps = nef

   ! begin column loop to calculate fractional area affected by fire

   do fc = 1, num_soilc
      c = filter_soilc(fc)

      ! dead fuel C (total litter + CWD)
      fuelc = totlitc(c) + cwdc(c)
      ! PET 5/20/08, test to increase fire area
      ! PET, 5/30/08. going back to original treatment using dead fuel only
      ! fuelc = fuelc + totvegc(c)

      ! m is the fractional soil mositure in the top layer (taken here
      ! as the top 0.5 m)
      ! PET 5/30/08 - note that this has been changed in Hydrology to use top 5 cm.
      m = max(0._r8,wf(c))


      ! Calculate the probability of at least one fire in a day
      ! in the gridcell. minfuel is the limit for dead fuels below which
      ! fire is assumed unable to spread.

      if (t_grnd(c)>SHR_CONST_TKFRZ .and. fuelc>minfuel .and. me(c)>0._r8 .and. m<=me(c)) then
         fire_prob(c) = exp(-SHR_CONST_PI * (m/me(c))**2)
      else
         fire_prob(c) = 0._r8
      end if

      ! Use e-folding to keep a running mean of daily fire probability,
      ! which is then used to calculate annual fractional area burned.
      ! mean_fire_prob corresponds to the variable s from Thonicke.
      ! fireseasonl corresponds to the variable N from Thonicke.
      ! ann_farea_burned corresponds to the variable A from Thonicke.

      mean_fire_prob(c) = (mean_fire_prob(c)*(ef_nsteps-1._r8) + fire_prob(c))/ef_nsteps
      fireseasonl(c) = mean_fire_prob(c) * dayspyr
      s2 = mean_fire_prob(c)-1._r8
      ann_farea_burned(c) = mean_fire_prob(c)*exp(s2/(0.45_r8*(s2**3) + 2.83_r8*(s2**2) + 2.96_r8*s2 + 1.04_r8))

      ! Estimate the fractional area of the column affected by fire in this time step.
      ! Over a year this should sum to a value near the annual
      ! fractional area burned from equations above.

      if (fireseasonl(c) > 0._r8) then
         farea_burned(c) = (fire_prob(c)/fireseasonl(c)) * ann_farea_burned(c) * (dt/secspday)
      else
         farea_burned(c) = 0._r8
      end if

      if (use_nofire) then
         ! set the fire area 0 if NOFIRE flag is on
         farea_burned(c) = 0._r8
      end if

   end do  ! end of column loop

end subroutine CNFireArea
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNFireFluxes
!
! !INTERFACE:
subroutine CNFireFluxes (num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Fire effects routine for coupled carbon-nitrogen code (CN).
! Relies primarily on estimate of fractional area burned in this
! timestep, from CNFireArea().
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn()
!
! !REVISION HISTORY:
! 7/23/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   real(r8), pointer :: nind(:)         ! number of individuals (#/m2)
   integer , pointer :: ivt(:)          ! pft vegetation type
   real(r8), pointer :: woody(:)        ! binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: resist(:)       ! resistance to fire (no units)
   integer , pointer :: pcolumn(:)      ! pft's column index
   real(r8), pointer :: farea_burned(:) ! timestep fractional area burned (proportion)
   real(r8), pointer :: m_cwdc_to_fire(:)
   real(r8), pointer :: m_deadcrootc_to_cwdc_fire(:)
   real(r8), pointer :: m_deadstemc_to_cwdc_fire(:)
   real(r8), pointer :: m_litr1c_to_fire(:)             
   real(r8), pointer :: m_litr2c_to_fire(:)             
   real(r8), pointer :: m_litr3c_to_fire(:)             
   real(r8), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real(r8), pointer :: m_cwdn_to_fire(:)              
   real(r8), pointer :: m_deadcrootn_to_cwdn_fire(:)
   real(r8), pointer :: m_deadstemn_to_cwdn_fire(:)
   real(r8), pointer :: m_litr1n_to_fire(:)             
   real(r8), pointer :: m_litr2n_to_fire(:)             
   real(r8), pointer :: m_litr3n_to_fire(:)             
   real(r8), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real(r8), pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real(r8), pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real(r8), pointer :: litr3n(:)             ! (gN/m2) litter lignin N
   real(r8), pointer :: m_deadcrootc_storage_to_fire(:) 
   real(r8), pointer :: m_deadcrootc_to_fire(:)         
   real(r8), pointer :: m_deadcrootc_to_litter_fire(:)         
   real(r8), pointer :: m_deadcrootc_xfer_to_fire(:)
   real(r8), pointer :: m_deadstemc_storage_to_fire(:)  
   real(r8), pointer :: m_deadstemc_to_fire(:)
   real(r8), pointer :: m_deadstemc_to_litter_fire(:)
   real(r8), pointer :: m_deadstemc_to_litter(:)
   real(r8), pointer :: m_livestemc_to_litter(:)
   real(r8), pointer :: m_deadcrootc_to_litter(:)
   real(r8), pointer :: m_livecrootc_to_litter(:)
   real(r8), pointer :: m_deadstemc_xfer_to_fire(:) 
   real(r8), pointer :: m_frootc_storage_to_fire(:)     
   real(r8), pointer :: m_frootc_to_fire(:)             
   real(r8), pointer :: m_frootc_xfer_to_fire(:)    
   real(r8), pointer :: m_gresp_storage_to_fire(:)      
   real(r8), pointer :: m_gresp_xfer_to_fire(:)    
   real(r8), pointer :: m_leafc_storage_to_fire(:)      
   real(r8), pointer :: m_leafc_to_fire(:)             
   real(r8), pointer :: m_leafc_xfer_to_fire(:)     
   real(r8), pointer :: m_livecrootc_storage_to_fire(:) 
   real(r8), pointer :: m_livecrootc_to_fire(:)         
   real(r8), pointer :: m_livecrootc_xfer_to_fire(:)
   real(r8), pointer :: m_livestemc_storage_to_fire(:)  
   real(r8), pointer :: m_livestemc_to_fire(:)          
   real(r8), pointer :: m_livestemc_xfer_to_fire(:) 
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    !(gC/m2) dead coarse root C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafcmax(:)           ! (gC/m2) ann max leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    !(gC/m2) live coarse root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: m_deadcrootn_storage_to_fire(:) 
   real(r8), pointer :: m_deadcrootn_to_fire(:)         
   real(r8), pointer :: m_deadcrootn_to_litter_fire(:)         
   real(r8), pointer :: m_deadcrootn_xfer_to_fire(:)
   real(r8), pointer :: m_deadstemn_storage_to_fire(:)  
   real(r8), pointer :: m_deadstemn_to_fire(:)          
   real(r8), pointer :: m_deadstemn_to_litter_fire(:)          
   real(r8), pointer :: m_deadstemn_xfer_to_fire(:) 
   real(r8), pointer :: m_frootn_storage_to_fire(:)     
   real(r8), pointer :: m_frootn_to_fire(:)             
   real(r8), pointer :: m_frootn_xfer_to_fire(:)    
   real(r8), pointer :: m_leafn_storage_to_fire(:)      
   real(r8), pointer :: m_leafn_to_fire(:)              
   real(r8), pointer :: m_leafn_xfer_to_fire(:)     
   real(r8), pointer :: m_livecrootn_storage_to_fire(:) 
   real(r8), pointer :: m_livecrootn_to_fire(:)         
   real(r8), pointer :: m_livecrootn_xfer_to_fire(:)
   real(r8), pointer :: m_livestemn_storage_to_fire(:)  
   real(r8), pointer :: m_livestemn_to_fire(:)          
   real(r8), pointer :: m_livestemn_xfer_to_fire(:) 
   real(r8), pointer :: m_retransn_to_fire(:)           
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N 
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
!
! !OTHER LOCAL VARIABLES:
   !real(r8), parameter:: wcf = 0.2_r8 ! wood combustion fraction
   real(r8), parameter:: wcf = 0.4_r8 ! wood combustion fraction
   integer :: c,p                  ! indices
   integer :: fp,fc                ! filter indices
   real(r8):: f                    ! rate for fire effects (1/s)
   real(r8):: dt                   ! time step variable (s)
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers

    nind                           => pdgvs%nind
    ivt                            => pft%itype
    pcolumn                        => pft%column
    woody                          => pftcon%woody
    resist                         => pftcon%resist
    farea_burned                   => cps%farea_burned
    m_cwdc_to_fire                 => ccf%m_cwdc_to_fire
    m_deadcrootc_to_cwdc_fire      => ccf%m_deadcrootc_to_cwdc_fire
    m_deadstemc_to_cwdc_fire       => ccf%m_deadstemc_to_cwdc_fire
    m_litr1c_to_fire               => ccf%m_litr1c_to_fire
    m_litr2c_to_fire               => ccf%m_litr2c_to_fire
    m_litr3c_to_fire               => ccf%m_litr3c_to_fire
    cwdc                           => ccs%cwdc
    litr1c                         => ccs%litr1c
    litr2c                         => ccs%litr2c
    litr3c                         => ccs%litr3c
    m_cwdn_to_fire                 => cnf%m_cwdn_to_fire
    m_deadcrootn_to_cwdn_fire      => cnf%m_deadcrootn_to_cwdn_fire
    m_deadstemn_to_cwdn_fire       => cnf%m_deadstemn_to_cwdn_fire
    m_litr1n_to_fire               => cnf%m_litr1n_to_fire
    m_litr2n_to_fire               => cnf%m_litr2n_to_fire
    m_litr3n_to_fire               => cnf%m_litr3n_to_fire
    cwdn                           => cns%cwdn
    litr1n                         => cns%litr1n
    litr2n                         => cns%litr2n
    litr3n                         => cns%litr3n
    m_deadcrootc_storage_to_fire   => pcf%m_deadcrootc_storage_to_fire
    m_deadcrootc_to_fire           => pcf%m_deadcrootc_to_fire
    m_deadcrootc_to_litter_fire    => pcf%m_deadcrootc_to_litter_fire
    m_deadcrootc_xfer_to_fire      => pcf%m_deadcrootc_xfer_to_fire
    m_deadstemc_storage_to_fire    => pcf%m_deadstemc_storage_to_fire
    m_deadstemc_to_fire            => pcf%m_deadstemc_to_fire
    m_deadstemc_to_litter_fire     => pcf%m_deadstemc_to_litter_fire
    m_deadstemc_to_litter          => pcf%m_deadstemc_to_litter
    m_livestemc_to_litter          => pcf%m_livestemc_to_litter
    m_deadcrootc_to_litter         => pcf%m_deadcrootc_to_litter
    m_livecrootc_to_litter         => pcf%m_livecrootc_to_litter
    m_deadstemc_xfer_to_fire       => pcf%m_deadstemc_xfer_to_fire
    m_frootc_storage_to_fire       => pcf%m_frootc_storage_to_fire
    m_frootc_to_fire               => pcf%m_frootc_to_fire
    m_frootc_xfer_to_fire          => pcf%m_frootc_xfer_to_fire
    m_gresp_storage_to_fire        => pcf%m_gresp_storage_to_fire
    m_gresp_xfer_to_fire           => pcf%m_gresp_xfer_to_fire
    m_leafc_storage_to_fire        => pcf%m_leafc_storage_to_fire
    m_leafc_to_fire                => pcf%m_leafc_to_fire
    m_leafc_xfer_to_fire           => pcf%m_leafc_xfer_to_fire
    m_livecrootc_storage_to_fire   => pcf%m_livecrootc_storage_to_fire
    m_livecrootc_to_fire           => pcf%m_livecrootc_to_fire
    m_livecrootc_xfer_to_fire      => pcf%m_livecrootc_xfer_to_fire
    m_livestemc_storage_to_fire    => pcf%m_livestemc_storage_to_fire
    m_livestemc_to_fire            => pcf%m_livestemc_to_fire
    m_livestemc_xfer_to_fire       => pcf%m_livestemc_xfer_to_fire
    deadcrootc                     => pcs%deadcrootc
    deadcrootc_storage             => pcs%deadcrootc_storage
    deadcrootc_xfer                => pcs%deadcrootc_xfer
    deadstemc                      => pcs%deadstemc
    deadstemc_storage              => pcs%deadstemc_storage
    deadstemc_xfer                 => pcs%deadstemc_xfer
    frootc                         => pcs%frootc
    frootc_storage                 => pcs%frootc_storage
    frootc_xfer                    => pcs%frootc_xfer
    gresp_storage                  => pcs%gresp_storage
    gresp_xfer                     => pcs%gresp_xfer
    leafc                          => pcs%leafc
    leafcmax                       => pcs%leafcmax
    leafc_storage                  => pcs%leafc_storage
    leafc_xfer                     => pcs%leafc_xfer
    livecrootc                     => pcs%livecrootc
    livecrootc_storage             => pcs%livecrootc_storage
    livecrootc_xfer                => pcs%livecrootc_xfer
    livestemc                      => pcs%livestemc
    livestemc_storage              => pcs%livestemc_storage
    livestemc_xfer                 => pcs%livestemc_xfer
    m_deadcrootn_storage_to_fire   => pnf%m_deadcrootn_storage_to_fire
    m_deadcrootn_to_fire           => pnf%m_deadcrootn_to_fire
    m_deadcrootn_to_litter_fire    => pnf%m_deadcrootn_to_litter_fire
    m_deadcrootn_xfer_to_fire      => pnf%m_deadcrootn_xfer_to_fire
    m_deadstemn_storage_to_fire    => pnf%m_deadstemn_storage_to_fire
    m_deadstemn_to_fire            => pnf%m_deadstemn_to_fire
    m_deadstemn_to_litter_fire     => pnf%m_deadstemn_to_litter_fire
    m_deadstemn_xfer_to_fire       => pnf%m_deadstemn_xfer_to_fire
    m_frootn_storage_to_fire       => pnf%m_frootn_storage_to_fire
    m_frootn_to_fire               => pnf%m_frootn_to_fire
    m_frootn_xfer_to_fire          => pnf%m_frootn_xfer_to_fire
    m_leafn_storage_to_fire        => pnf%m_leafn_storage_to_fire
    m_leafn_to_fire                => pnf%m_leafn_to_fire
    m_leafn_xfer_to_fire           => pnf%m_leafn_xfer_to_fire
    m_livecrootn_storage_to_fire   => pnf%m_livecrootn_storage_to_fire
    m_livecrootn_to_fire           => pnf%m_livecrootn_to_fire
    m_livecrootn_xfer_to_fire      => pnf%m_livecrootn_xfer_to_fire
    m_livestemn_storage_to_fire    => pnf%m_livestemn_storage_to_fire
    m_livestemn_to_fire            => pnf%m_livestemn_to_fire
    m_livestemn_xfer_to_fire       => pnf%m_livestemn_xfer_to_fire
    m_retransn_to_fire             => pnf%m_retransn_to_fire
    deadcrootn                     => pns%deadcrootn
    deadcrootn_storage             => pns%deadcrootn_storage
    deadcrootn_xfer                => pns%deadcrootn_xfer
    deadstemn                      => pns%deadstemn
    deadstemn_storage              => pns%deadstemn_storage
    deadstemn_xfer                 => pns%deadstemn_xfer
    frootn                         => pns%frootn
    frootn_storage                 => pns%frootn_storage
    frootn_xfer                    => pns%frootn_xfer
    leafn                          => pns%leafn
    leafn_storage                  => pns%leafn_storage
    leafn_xfer                     => pns%leafn_xfer
    livecrootn                     => pns%livecrootn
    livecrootn_storage             => pns%livecrootn_storage
    livecrootn_xfer                => pns%livecrootn_xfer
    livestemn                      => pns%livestemn
    livestemn_storage              => pns%livestemn_storage
    livestemn_xfer                 => pns%livestemn_xfer
    retransn                       => pns%retransn


   ! Get model step size

   dt = real( get_step_size(), r8 )

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      ! get the column-level fractional area burned for this timestep
      ! and convert to a rate per second, then scale by the pft-level
      ! fire resistance
      f = (farea_burned(c) / dt) * (1._r8 - resist(ivt(p)))
      
      ! apply this rate to the pft state variables to get flux rates

      ! NOTE: the deadstem and deadcroot pools are only partly consumed
      ! by fire, and the remaining affected fraction goes to the column-level
      ! as litter (coarse woody debris). This is controlled by wcf, the woody
      ! combustion fraction.

      ! carbon fluxes
      m_leafc_to_fire(p)               =  leafc(p)              * f
      m_leafc_storage_to_fire(p)       =  leafc_storage(p)      * f
      m_leafc_xfer_to_fire(p)          =  leafc_xfer(p)         * f
      m_frootc_to_fire(p)              =  frootc(p)             * f
      m_frootc_storage_to_fire(p)      =  frootc_storage(p)     * f
      m_frootc_xfer_to_fire(p)         =  frootc_xfer(p)        * f
      m_livestemc_to_fire(p)           =  livestemc(p)          * f
      m_livestemc_storage_to_fire(p)   =  livestemc_storage(p)  * f
      m_livestemc_xfer_to_fire(p)      =  livestemc_xfer(p)     * f
      m_deadstemc_to_fire(p)           =  deadstemc(p)          * f*wcf
      m_deadstemc_to_litter_fire(p)    =  deadstemc(p)          * f*(1._r8 - wcf)
      m_deadstemc_storage_to_fire(p)   =  deadstemc_storage(p)  * f
      m_deadstemc_xfer_to_fire(p)      =  deadstemc_xfer(p)     * f
      m_livecrootc_to_fire(p)          =  livecrootc(p)         * f
      m_livecrootc_storage_to_fire(p)  =  livecrootc_storage(p) * f
      m_livecrootc_xfer_to_fire(p)     =  livecrootc_xfer(p)    * f
      m_deadcrootc_to_fire(p)          =  deadcrootc(p)         * f*wcf
      m_deadcrootc_to_litter_fire(p)   =  deadcrootc(p)         * f*(1._r8 - wcf)
      m_deadcrootc_storage_to_fire(p)  =  deadcrootc_storage(p) * f
      m_deadcrootc_xfer_to_fire(p)     =  deadcrootc_xfer(p)    * f
      m_gresp_storage_to_fire(p)       =  gresp_storage(p)      * f
      m_gresp_xfer_to_fire(p)          =  gresp_xfer(p)         * f

      ! nitrogen fluxes
      m_leafn_to_fire(p)               =  leafn(p)              * f
      m_leafn_storage_to_fire(p)       =  leafn_storage(p)      * f
      m_leafn_xfer_to_fire(p)          =  leafn_xfer(p)         * f
      m_frootn_to_fire(p)              =  frootn(p)             * f
      m_frootn_storage_to_fire(p)      =  frootn_storage(p)     * f
      m_frootn_xfer_to_fire(p)         =  frootn_xfer(p)        * f
      m_livestemn_to_fire(p)           =  livestemn(p)          * f
      m_livestemn_storage_to_fire(p)   =  livestemn_storage(p)  * f
      m_livestemn_xfer_to_fire(p)      =  livestemn_xfer(p)     * f
      m_deadstemn_to_fire(p)           =  deadstemn(p)          * f*wcf
      m_deadstemn_to_litter_fire(p)    =  deadstemn(p)          * f*(1._r8 - wcf)
      m_deadstemn_storage_to_fire(p)   =  deadstemn_storage(p)  * f
      m_deadstemn_xfer_to_fire(p)      =  deadstemn_xfer(p)     * f
      m_livecrootn_to_fire(p)          =  livecrootn(p)         * f
      m_livecrootn_storage_to_fire(p)  =  livecrootn_storage(p) * f
      m_livecrootn_xfer_to_fire(p)     =  livecrootn_xfer(p)    * f
      m_deadcrootn_to_fire(p)          =  deadcrootn(p)         * f*wcf
      m_deadcrootn_to_litter_fire(p)   =  deadcrootn(p)         * f*(1._r8 - wcf)
      m_deadcrootn_storage_to_fire(p)  =  deadcrootn_storage(p) * f
      m_deadcrootn_xfer_to_fire(p)     =  deadcrootn_xfer(p)    * f
      m_retransn_to_fire(p)            =  retransn(p)           * f

      if (use_cndv) then
         ! Carbon per individual (c) remains constant in gap mortality & fire
         ! but individuals are removed from the population P (#/m2 naturally
         ! vegetated area), so
         !
         ! c = Cnew*FPC/Pnew = Cold*FPC/Pold
         !
         ! where C = carbon/m2 pft area & FPC = pft area/naturally vegetated area.
         ! FPC does not change from mortality or fire. FPC changes from Light and
         ! Establishment at the end of the year. So...
         !
         ! Pnew = Pold * Cnew / Cold
         !
         ! where "new" refers to after mortality & fire, while "old" refers to
         ! before mortality & fire. For C I use total wood. (slevis)
         !
         ! nind calculation placed here for convenience; nind could be updated
         ! once per year instead if we saved Cold for that calculation;
         ! as is, nind slowly decreases through the year, while fpcgrid remains
         ! unchanged; this affects the htop calculation in CNVegStructUpdate
         
         if (woody(ivt(p)) == 1._r8) then
            if (livestemc(p)+deadstemc(p)+m_livestemc_to_litter(p)*dt+ &
                 m_deadstemc_to_litter(p)*dt > 0._r8) then
               nind(p) = nind(p) * (livestemc(p)  + deadstemc(p) +       &
                    livecrootc(p) + deadcrootc(p) - dt * &
                    (m_livestemc_to_fire(p)  +           &
                    m_livecrootc_to_fire(p) +           &
                    m_deadstemc_to_fire(p)  +           &
                    m_deadcrootc_to_fire(p) +           &
                    m_deadcrootc_to_litter_fire(p) +    &
                    m_deadstemc_to_litter_fire(p))) /   &
                    (livestemc(p)  + deadstemc(p) +       &
                    livecrootc(p) + deadcrootc(p) + dt * &
                    (m_livestemc_to_litter(p)  +         &
                    m_livecrootc_to_litter(p) +         &
                    m_deadcrootc_to_litter(p) +         &
                    m_deadstemc_to_litter(p)))
            else
               nind(p) = 0._r8
            end if
         end if
         
         ! annual dgvm calculations use lm_ind = leafcmax * fpcgrid / nind
         ! leafcmax is reset to 0 once per yr
         ! could calculate leafcmax in CSummary instead; if so, should remove
         ! subtraction of m_leafc_to_fire(p)*dt from the calculation (slevis)
         
         leafcmax(p) = max(leafc(p)-m_leafc_to_fire(p)*dt, leafcmax(p))
         if (ivt(p) == 0) leafcmax(p) = 0._r8
      end if

   end do  ! end of pfts loop

   ! send the fire affected but uncombusted woody fraction to the column-level cwd fluxes
   ! use p2c for weighted averaging from pft to column
   call p2c(num_soilc, filter_soilc, m_deadstemc_to_litter_fire, m_deadstemc_to_cwdc_fire)
   call p2c(num_soilc, filter_soilc, m_deadcrootc_to_litter_fire, m_deadcrootc_to_cwdc_fire)
   call p2c(num_soilc, filter_soilc, m_deadstemn_to_litter_fire, m_deadstemn_to_cwdn_fire)
   call p2c(num_soilc, filter_soilc, m_deadcrootn_to_litter_fire, m_deadcrootn_to_cwdn_fire)

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! get the column-level fractional area burned for this timestep
      ! and convert to a rate per second, then scale by the pft-level
      ! fire resistance

      f = farea_burned(c) / dt

      ! apply this rate to the column state variables to get flux rates

      ! NOTE: the coarse woody debris pools are only partly consumed
      ! by fire. This is controlled by wcf, the woody
      ! combustion fraction. For now using the same fraction for standing
      ! wood (deadstem and deadcroot pools) and woody litter (cwd pools).
      ! May be a good idea later to modify this to use different fractions
      ! for different woody pools, or make the combustion fraction a dynamic
      ! variable.

      ! carbon fluxes
      m_litr1c_to_fire(c) = litr1c(c) * f
      m_litr2c_to_fire(c) = litr2c(c) * f
      m_litr3c_to_fire(c) = litr3c(c) * f
      m_cwdc_to_fire(c)   = cwdc(c)   * f*wcf

      ! nitrogen fluxes
      m_litr1n_to_fire(c) = litr1n(c) * f
      m_litr2n_to_fire(c) = litr2n(c) * f
      m_litr3n_to_fire(c) = litr3n(c) * f
      m_cwdn_to_fire(c)   = cwdn(c)   * f*wcf

   end do  ! end of column loop

end subroutine CNFireFluxes
!-----------------------------------------------------------------------

end module CNFireMod
