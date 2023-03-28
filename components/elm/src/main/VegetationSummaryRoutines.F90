module VegetationSummaryRoutinesMod

   use elm_varpar
   !use VegetationDataType , only : veg_cf, veg_cs, veg_pf, veg_ps
   use VegetationDataType , only : vegetation_carbon_flux, vegetation_carbon_state
   use VegetationDataType , only : vegetation_nitrogen_flux, vegetation_nitrogen_state
   use VegetationDataType , only : vegetation_phosphorus_flux, vegetation_phosphorus_state
   use VegetationType , only : veg_pp
   use ColumnType     , only : col_pp
   use ColumnDataType , only : column_carbon_flux, column_carbon_state
   use ColumnDataType , only : column_nitrogen_flux, column_nitrogen_state
   use ColumnDataType , only : column_phosphorus_flux, column_phosphorus_state
   use CNDecompCascadeConType

   use elm_varpar
   use elm_varctl
   use elm_varcon

   public :: veg_cf_summary_acc
   public :: veg_cs_summary_acc
   public :: veg_nf_summary_acc
   public :: veg_ns_summary_acc
   public :: veg_pf_summary_acc
   public :: veg_ps_summary_acc
   public :: summary_p2c
   public :: veg_cf_setvalues_acc
   public :: veg_pf_setvalues_acc
   public :: veg_nf_setvalues_acc

   public :: summary_veg_flux_p2c, summary_veg_state_p2c
   public :: veg_cf_summary_rr

contains
   !------------------------------------------------------------
  subroutine veg_cf_summary_rr(this, num_soilp, filter_soilp, num_soilc, filter_soilc, col_cf_input)
   !
   ! !DESCRIPTION:
   ! summarize root respiration
   !
   ! !USES:
   use subgridAveMod, only : p2c_1d_filter_parallel
   !
   ! !ARGUMENTS:
   type(vegetation_carbon_flux) :: this
   integer, intent(in) :: num_soilp
   integer, intent(in) :: filter_soilp(:)
   integer, intent(in) :: num_soilc
   integer, intent(in) :: filter_soilc(:)
   type(column_carbon_flux), intent(inout) :: col_cf_input
   !
   ! !LOCAL VARIABLES
   integer :: fp, p
   !------------------------------------------------------------
   !$acc parallel loop independent gang vector private(p) default(present)
   do fp = 1,num_soilp
     p = filter_soilp(fp)
     ! root respiration (RR)
     this%rr(p) = &
     this%froot_mr(p) + &
     this%cpool_froot_gr(p) + &
     this%cpool_livecroot_gr(p) + &
     this%cpool_deadcroot_gr(p) + &
     this%transfer_froot_gr(p) + &
     this%transfer_livecroot_gr(p) + &
     this%transfer_deadcroot_gr(p) + &
     this%cpool_froot_storage_gr(p) + &
     this%cpool_livecroot_storage_gr(p) + &
     this%cpool_deadcroot_storage_gr(p)
   enddo
   call p2c_1d_filter_parallel(num_soilc, filter_soilc, &
           this%rr, col_cf_input%rr)

 end subroutine veg_cf_summary_rr


   subroutine veg_cf_summary_for_ch4_acc( this,num_soilp,filter_soilp )
     !
     ! !DESCRIPTION:
     ! summarize vegetation-level fluxes for methane calculation
     !
     ! !USES:
     !
     ! !ARGUMENTS:
     type(vegetation_carbon_flux) :: this
     integer , intent(in) :: num_soilp 
     integer , intent(in) :: filter_soilp(:)
     !!LOCAL VARIABLES:
     integer :: fp, p 

    !$acc parallel loop independent gang vector default(present) 
    do fp = 1, num_soilp 
      p = filter_soilp(fp)
        !------------------------------------------------------------
        ! aboveground NPP: leaf, live stem, dead stem (AGNPP)
        ! This is supposed to correspond as closely as possible to
        ! field measurements of AGNPP, so it ignores the storage pools
        ! and only treats the fluxes into displayed pools.

        this%agnpp(p) = &
             this%cpool_to_leafc(p)                  + &
             this%leafc_xfer_to_leafc(p)             + &
             this%cpool_to_livestemc(p)              + &
             this%livestemc_xfer_to_livestemc(p)     + &
             this%cpool_to_deadstemc(p)              + &
             this%deadstemc_xfer_to_deadstemc(p)

        if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
           this%agnpp(p) =                    &
                this%agnpp(p)               + &
                this%cpool_to_grainc(p)     + &
                this%grainc_xfer_to_grainc(p)
        endif

        ! belowground NPP: fine root, live coarse root, dead coarse root (BGNPP)
        ! This is supposed to correspond as closely as possible to
        ! field measurements of BGNPP, so it ignores the storage pools
        ! and only treats the fluxes into displayed pools.
        this%bgnpp(p) = &
             this%cpool_to_frootc(p)                   + &
             this%frootc_xfer_to_frootc(p)             + &
             this%cpool_to_livecrootc(p)               + &
             this%livecrootc_xfer_to_livecrootc(p)     + &
             this%cpool_to_deadcrootc(p)               + &
             this%deadcrootc_xfer_to_deadcrootc(p)

        this%agwdnpp(p) = &
             this%cpool_to_livestemc(p)              + &
             this%livestemc_xfer_to_livestemc(p)     + &
             this%cpool_to_deadstemc(p)              + &
             this%deadstemc_xfer_to_deadstemc(p)
   end do 

  end subroutine veg_cf_summary_for_ch4_acc

  !-----------------------------------------------------------------------
  subroutine veg_cf_summary_acc(this, num_soilp,filter_soilp, isotope)
    !
    ! !DESCRIPTION:
    ! patch-level carbon flux summary calculations
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(vegetation_carbon_flux)                 :: this
    integer , intent(in) :: num_soilp 
    integer , intent(in) :: filter_soilp(:)
    integer ,intent(in)  :: isotope
    !!LOCAL VARIABLES:
    integer :: fp, p 

    !$acc data copyin(isotope)

   !$acc parallel loop independent gang vector default(present) 
   do fp = 1, num_soilp 
      p = filter_soilp(fp)
      ! calculate pft-level summary carbon fluxes and states

       ! gross primary production (GPP)
       this%gpp(p) = &
            this%psnsun_to_cpool(p) + &
            this%psnshade_to_cpool(p)

       ! maintenance respiration (MR)
       this%leaf_mr(p)      = this%leaf_curmr(p)      + this%leaf_xsmr(p)
       this%froot_mr(p)     = this%froot_curmr(p)     + this%froot_xsmr(p)
       this%livestem_mr(p)  = this%livestem_curmr(p)  + this%livestem_xsmr(p)
       this%livecroot_mr(p) = this%livecroot_curmr(p) + this%livecroot_xsmr(p)

       this%mr(p)  = &
            this%leaf_mr(p)     + &
            this%froot_mr(p)    + &
            this%livestem_mr(p) + &
            this%livecroot_mr(p)

       ! growth respiration (GR)
       ! current GR is respired this time step for new growth displayed in this timestep
       this%current_gr(p) = &
            this%cpool_leaf_gr(p)      + &
            this%cpool_froot_gr(p)     + &
            this%cpool_livestem_gr(p)  + &
            this%cpool_deadstem_gr(p)  + &
            this%cpool_livecroot_gr(p) + &
            this%cpool_deadcroot_gr(p)

       ! transfer GR is respired this time step for transfer growth displayed in this timestep
       this%transfer_gr(p) = &
            this%transfer_leaf_gr(p)      + &
            this%transfer_froot_gr(p)     + &
            this%transfer_livestem_gr(p)  + &
            this%transfer_deadstem_gr(p)  + &
            this%transfer_livecroot_gr(p) + &
            this%transfer_deadcroot_gr(p)

       ! storage GR is respired this time step for growth sent to storage for later display
       this%storage_gr(p) = &
            this%cpool_leaf_storage_gr(p)      + &
            this%cpool_froot_storage_gr(p)     + &
            this%cpool_livestem_storage_gr(p)  + &
            this%cpool_deadstem_storage_gr(p)  + &
            this%cpool_livecroot_storage_gr(p) + &
            this%cpool_deadcroot_storage_gr(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%mr(p) = &
               this%mr(p) + &
               this%grain_mr(p)

          this%current_gr(p) = &
               this%current_gr(p) + &
               this%cpool_grain_gr(p)

          this%transfer_gr(p) = &
               this%transfer_gr(p) + &
               this%transfer_grain_gr(p)

          this%storage_gr(p) = &
               this%storage_gr(p) + &
               this%cpool_grain_storage_gr(p)
       end if

       ! GR is the sum of current + transfer + storage GR
       this%gr(p) = &
            this%current_gr(p)  + &
            this%transfer_gr(p) + &
            this%storage_gr(p)

       ! autotrophic respiration (AR)
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%ar(p) = &
               this%mr(p) + &
               this%gr(p) + &
               this%xr(p) + &
               this%xsmrpool_to_atm(p) ! xsmr... is -ve (slevis)
          if (nu_com .ne. 'RD' ) then
             this%ar(p) = this%ar(p) + &
                  this%xsmrpool_turnover(p)
          end if
       else
          this%ar(p) = &
               this%mr(p) + &
               this%gr(p) + &
               this%xr(p)
           if (nu_com .ne. 'RD' ) then
              this%ar(p) = this%ar(p) + &
                   this%xsmrpool_turnover(p)
           end if
       end if

       ! net primary production (NPP)
       this%npp(p) = &
            this%gpp(p) - &
            this%ar(p)

       ! update the annual NPP accumulator, for use in allocation code
       if (isotope == 2) then
          this%tempsum_npp(p) = &
               this%tempsum_npp(p) + &
               this%npp(p)
       end if

       ! litterfall (LITFALL)

       this%litfall(p) = &
            this%leafc_to_litter(p)                     + &
            this%frootc_to_litter(p)                    + &
            this%m_leafc_to_litter(p)                   + &
            this%m_leafc_storage_to_litter(p)           + &
            this%m_leafc_xfer_to_litter(p)              + &
            this%m_frootc_to_litter(p)                  + &
            this%m_frootc_storage_to_litter(p)          + &
            this%m_frootc_xfer_to_litter(p)             + &
            this%m_livestemc_to_litter(p)               + &
            this%m_livestemc_storage_to_litter(p)       + &
            this%m_livestemc_xfer_to_litter(p)          + &
            this%m_deadstemc_to_litter(p)               + &
            this%m_deadstemc_storage_to_litter(p)       + &
            this%m_deadstemc_xfer_to_litter(p)          + &
            this%m_livecrootc_to_litter(p)              + &
            this%m_livecrootc_storage_to_litter(p)      + &
            this%m_livecrootc_xfer_to_litter(p)         + &
            this%m_deadcrootc_to_litter(p)              + &
            this%m_deadcrootc_storage_to_litter(p)      + &
            this%m_deadcrootc_xfer_to_litter(p)         + &
            this%m_gresp_storage_to_litter(p)           + &
            this%m_gresp_xfer_to_litter(p)              + &
            this%m_leafc_to_litter_fire(p)              + &
            this%m_leafc_storage_to_litter_fire(p)      + &
            this%m_leafc_xfer_to_litter_fire(p)         + &
            this%m_livestemc_to_litter_fire(p)          + &
            this%m_livestemc_storage_to_litter_fire(p)  + &
            this%m_livestemc_xfer_to_litter_fire(p)     + &
            this%m_deadstemc_to_litter_fire(p)          + &
            this%m_deadstemc_storage_to_litter_fire(p)  + &
            this%m_deadstemc_xfer_to_litter_fire(p)     + &
            this%m_frootc_to_litter_fire(p)             + &
            this%m_frootc_storage_to_litter_fire(p)     + &
            this%m_frootc_xfer_to_litter_fire(p)        + &
            this%m_livecrootc_to_litter_fire(p)         + &
            this%m_livecrootc_storage_to_litter_fire(p) + &
            this%m_livecrootc_xfer_to_litter_fire(p)    + &
            this%m_deadcrootc_to_litter_fire(p)         + &
            this%m_deadcrootc_storage_to_litter_fire(p) + &
            this%m_deadcrootc_xfer_to_litter_fire(p)    + &
            this%m_gresp_storage_to_litter_fire(p)      + &
            this%m_gresp_xfer_to_litter_fire(p)         + &
            this%hrv_leafc_to_litter(p)                 + &
            this%hrv_leafc_storage_to_litter(p)         + &
            this%hrv_leafc_xfer_to_litter(p)            + &
            this%hrv_frootc_to_litter(p)                + &
            this%hrv_frootc_storage_to_litter(p)        + &
            this%hrv_frootc_xfer_to_litter(p)           + &
            this%hrv_livestemc_to_litter(p)             + &
            this%hrv_livestemc_storage_to_litter(p)     + &
            this%hrv_livestemc_xfer_to_litter(p)        + &
            this%hrv_deadstemc_storage_to_litter(p)     + &
            this%hrv_deadstemc_xfer_to_litter(p)        + &
            this%hrv_livecrootc_to_litter(p)            + &
            this%hrv_livecrootc_storage_to_litter(p)    + &
            this%hrv_livecrootc_xfer_to_litter(p)       + &
            this%hrv_deadcrootc_to_litter(p)            + &
            this%hrv_deadcrootc_storage_to_litter(p)    + &
            this%hrv_deadcrootc_xfer_to_litter(p)       + &
            this%hrv_gresp_storage_to_litter(p)         + &
            this%hrv_gresp_xfer_to_litter(p)            + &
            this%hrv_cpool_to_litter(p)

       ! patch-level fire losses (VEGFIRE)
       this%vegfire(p) = 0._r8

       ! patch-level wood harvest
       this%wood_harvestc(p) = &
            this%hrv_deadstemc_to_prod10c(p) + &
            this%hrv_deadstemc_to_prod100c(p)
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%wood_harvestc(p) = &
               this%wood_harvestc(p) + &
               this%hrv_cropc_to_prod1c(p)
       end if

       ! patch-level carbon losses to fire changed by F. Li and S. Levis
       this%fire_closs(p) = &
            this%m_leafc_to_fire(p)                + &
            this%m_leafc_storage_to_fire(p)        + &
            this%m_leafc_xfer_to_fire(p)           + &
            this%m_frootc_to_fire(p)               + &
            this%m_frootc_storage_to_fire(p)       + &
            this%m_frootc_xfer_to_fire(p)          + &
            this%m_livestemc_to_fire(p)            + &
            this%m_livestemc_storage_to_fire(p)    + &
            this%m_livestemc_xfer_to_fire(p)       + &
            this%m_deadstemc_to_fire(p)            + &
            this%m_deadstemc_storage_to_fire(p)    + &
            this%m_deadstemc_xfer_to_fire(p)       + &
            this%m_livecrootc_to_fire(p)           + &
            this%m_livecrootc_storage_to_fire(p)   + &
            this%m_livecrootc_xfer_to_fire(p)      + &
            this%m_deadcrootc_to_fire(p)           + &
            this%m_deadcrootc_storage_to_fire(p)   + &
            this%m_deadcrootc_xfer_to_fire(p)      + &
            this%m_gresp_storage_to_fire(p)        + &
            this%m_gresp_xfer_to_fire(p)           + &
            this%m_cpool_to_fire(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%litfall(p) =                  &
               this%litfall(p)             + &
               this%livestemc_to_litter(p) + &
               this%grainc_to_food(p)
       end if

       ! new summary variables for CLAMP

       ! (FROOTC_ALLOC) - fine root C allocation
       this%frootc_alloc(p) = &
            this%frootc_xfer_to_frootc(p)    + &
            this%cpool_to_frootc(p)

       ! (FROOTC_LOSS) - fine root C loss changed by F. Li and S. Levis
       this%frootc_loss(p) = &
            this%m_frootc_to_litter(p)       + &
            this%m_frootc_to_fire(p)         + &
            this%m_frootc_to_litter_fire(p)  + &
            this%hrv_frootc_to_litter(p)     + &
            this%frootc_to_litter(p)

       ! (LEAFC_ALLOC) - leaf C allocation
       this%leafc_alloc(p) = &
            this%leafc_xfer_to_leafc(p)    + &
            this%cpool_to_leafc(p)

       ! (LEAFC_LOSS) - leaf C loss changed by F. Li and S. Levis
       this%leafc_loss(p) = &
            this%m_leafc_to_litter(p)      + &
            this%m_leafc_to_fire(p)        + &
            this%m_leafc_to_litter_fire(p) + &
            this%hrv_leafc_to_litter(p)    + &
            this%leafc_to_litter(p)

       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%leafc_loss(p) = &
               this%leafc_loss(p) + &
               this%hrv_leafc_to_prod1c(p)
       end if


       ! (WOODC_ALLOC) - wood C allocation
       this%woodc_alloc(p) = &
            this%livestemc_xfer_to_livestemc(p)   + &
            this%deadstemc_xfer_to_deadstemc(p)   + &
            this%livecrootc_xfer_to_livecrootc(p) + &
            this%deadcrootc_xfer_to_deadcrootc(p) + &
            this%cpool_to_livestemc(p)            + &
            this%cpool_to_deadstemc(p)            + &
            this%cpool_to_livecrootc(p)           + &
            this%cpool_to_deadcrootc(p)

       ! (WOODC_LOSS) - wood C loss
       this%woodc_loss(p) = &
            this%m_livestemc_to_litter(p)            + &
            this%m_deadstemc_to_litter(p)            + &
            this%m_livecrootc_to_litter(p)           + &
            this%m_deadcrootc_to_litter(p)           + &
            this%m_livestemc_to_fire(p)              + &
            this%m_deadstemc_to_fire(p)              + &
            this%m_livecrootc_to_fire(p)             + &
            this%m_deadcrootc_to_fire(p)             + &
            this%hrv_livestemc_to_litter(p)          + &
            this%hrv_livestemc_storage_to_litter(p)  + &
            this%hrv_livestemc_xfer_to_litter(p)     + &
            this%hrv_deadstemc_to_prod10c(p)         + &
            this%hrv_deadstemc_to_prod100c(p)        + &
            this%hrv_deadstemc_storage_to_litter(p)  + &
            this%hrv_deadstemc_xfer_to_litter(p)     + &
            this%hrv_livecrootc_to_litter(p)         + &
            this%hrv_livecrootc_storage_to_litter(p) + &
            this%hrv_livecrootc_xfer_to_litter(p)    + &
            this%hrv_deadcrootc_to_litter(p)         + &
            this%hrv_deadcrootc_storage_to_litter(p) + &
            this%hrv_deadcrootc_xfer_to_litter(p)
       ! putting the harvested crop stem and grain in the wood loss bdrewniak
       if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%woodc_loss(p) = &
               this%woodc_loss(p) + &
               this%hrv_grainc_to_prod1c(p) + &
               this%hrv_livestemc_to_prod1c(p)
       end if
   end do 
   !$acc end data 
  end subroutine veg_cf_summary_acc

  !-----------------------------------------------------------------------
  subroutine veg_cs_summary_acc(this, num_soilp, filter_soilp)
    !
    ! !DESCRIPTION:
    ! Vegetation-level carbon state summary calculations
    ! !ARGUMENTS:
    type(vegetation_carbon_state)                :: this
    integer , intent(in) :: num_soilp 
    integer , intent(in) :: filter_soilp(:)

    !!LOCAL VARIABLES:
    integer :: fp, p 

   !$acc parallel loop independent gang vector default(present) 
   do fp = 1, num_soilp 
      p = filter_soilp(fp)

      ! displayed vegetation carbon, excluding storage and cpool (DISPVEGC)
      this%dispvegc(p) =        &
           this%leafc(p)      + &
           this%frootc(p)     + &
           this%livestemc(p)  + &
           this%deadstemc(p)  + &
           this%livecrootc(p) + &
           this%deadcrootc(p)

      ! stored vegetation carbon, excluding cpool (STORVEGC)
      this%storvegc(p) =                &
           this%cpool(p)              + &
           this%leafc_storage(p)      + &
           this%frootc_storage(p)     + &
           this%livestemc_storage(p)  + &
           this%deadstemc_storage(p)  + &
           this%livecrootc_storage(p) + &
           this%deadcrootc_storage(p) + &
           this%leafc_xfer(p)         + &
           this%frootc_xfer(p)        + &
           this%livestemc_xfer(p)     + &
           this%deadstemc_xfer(p)     + &
           this%livecrootc_xfer(p)    + &
           this%deadcrootc_xfer(p)    + &
           this%gresp_storage(p)      + &
           this%gresp_xfer(p)

      if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%storvegc(p) =            &
               this%storvegc(p)       + &
               this%grainc_storage(p) + &
               this%grainc_xfer(p)

          this%dispvegc(p) =            &
               this%dispvegc(p)       + &
               this%grainc(p)
      end if

      ! total vegetation carbon, excluding cpool (TOTVEGC)
      this%totvegc(p) = &
           this%dispvegc(p) + &
           this%storvegc(p)

      ! total pft-level carbon, including xsmrpool, ctrunc
      this%totpftc(p) = &
           this%totvegc(p) + &
           this%xsmrpool(p) + &
           this%ctrunc(p)

      ! (WOODC) - wood C
      this%woodc(p) = &
           this%deadstemc(p)    + &
           this%livestemc(p)    + &
           this%deadcrootc(p)   + &
           this%livecrootc(p)

      this%totvegc_abg(p) = &
               this%leafc(p)              + &
               this%leafc_storage(p)      + &
               this%leafc_xfer(p)         + &
               this%livestemc(p)          + &
               this%livestemc_storage(p)  + &
               this%livestemc_xfer(p)     + &
               this%deadstemc(p)          + &
               this%deadstemc_storage(p)  + &
               this%deadstemc_xfer(p)
   end do 

  end subroutine veg_cs_summary_acc

subroutine summary_veg_state_p2c ( numfc, filterc,veg_cs, col_cs, &
                  veg_ps,col_ps, veg_ns, col_ns)
  !
  ! !DESCRIPTION:
  ! perform pft to column averaging for single level pft arrays
  !
  ! !ARGUMENTS:
  integer , intent(in)  :: numfc
  integer , intent(in)  :: filterc(numfc)
  type(vegetation_carbon_state), intent(inout)   :: veg_cs
  type(column_carbon_state)    , intent(inout)   :: col_cs  ! column-level state for p2c
  type(vegetation_phosphorus_state), intent(inout)::veg_ps
  type(column_phosphorus_state), intent(inout)   :: col_ps
  type(vegetation_nitrogen_state), intent(inout) :: veg_ns
  type(column_nitrogen_state)  , intent(inout)   :: col_ns

  ! !LOCAL VARIABLES:
  integer :: fc,c,p  ! indices
  real(r8) :: sum1, sum2, sum3, sum4
  real(r8) :: psum1, psum2, psum3
  real(r8) :: nsum1, nsum2, nsum3, nsum4


  !-----------------------------------------------------------------------
  !$acc enter data create(sum1,sum2,sum3,sum4,psum1,psum2,psum3,&
   !$acc nsum1, nsum2, nsum3, nsum4)
  !$acc parallel loop gang worker default(present) private(c)
  do fc = 1,numfc
     c = filterc(fc)
     sum1= 0._r8
     sum2= 0._r8
     sum3= 0._r8
     sum4= 0._r8
     psum1 = 0._r8; psum2 = 0._r8; psum3 = 0._r8;
     nsum1 = 0._r8; nsum2 = 0._r8; nsum3 = 0._r8; nsum4 = 0._r8
     !$acc loop vector reduction(+:sum1,sum2,sum3,sum4,psum1,psum2,psum3,nsum1, nsum2, nsum3, nsum4)
     do p = col_pp%pfti(c), col_pp%pftf(c)
        if (veg_pp%active(p)) then
            sum1 = sum1 + veg_cs%totpftc(p)           * veg_pp%wtcol(p)
            sum2 = sum2 + veg_cs%totvegc(p)           * veg_pp%wtcol(p)
            sum3 = sum3 + veg_cs%totvegc_abg(p)       * veg_pp%wtcol(p)
            sum4 = sum4 + veg_cs%cropseedc_deficit(p) * veg_pp%wtcol(p)
            ! phosphorus state:
            psum1 = psum1 + veg_ps%totvegp(p) * veg_pp%wtcol(p)
            psum2 = psum2 + veg_ps%totpftp(p) * veg_pp%wtcol(p)
            psum3 = psum3 + veg_ps%cropseedp_deficit(p) * veg_pp%wtcol(p)
            !Nitrogen State:
            nsum1 = nsum1 + veg_ns%totvegn(p) * veg_pp%wtcol(p)
            nsum2 = nsum2 + veg_ns%totpftn(p) * veg_pp%wtcol(p)
            nsum3 = nsum3 + veg_ns%cropseedn_deficit(p) * veg_pp%wtcol(p)
            nsum4 = nsum4 + veg_ns%plant_n_buffer(p) * veg_pp%wtcol(p)
         end if
     end do
     col_cs%totpftc(c)          = sum1
     col_cs%totvegc(c)          = sum2
     col_cs%totvegc_abg(c)      = sum3
     col_cs%cropseedc_deficit(c)= sum4
     ! phosphorus state:
     col_ps%totpftp(c)          = psum1
     col_ps%totvegp(c)          = psum2
     col_ps%cropseedp_deficit(c)= psum3
     ! nitrogen state:
     col_ns%totpftn(c)           = nsum1
     col_ns%totvegn(c)           = nsum2
     col_ns%cropseedn_deficit(c) = nsum3
     col_ns%plant_n_buffer(c)    = nsum4

  end do

end subroutine summary_veg_state_p2c

subroutine summary_veg_flux_p2c(numfc, filterc, veg_cf, col_cf,&
                     veg_pf,col_pf, veg_nf, col_nf)

   integer , intent(in)  :: numfc
   integer , intent(in)  :: filterc(numfc)
   type(vegetation_carbon_flux), intent(inout)   :: veg_cf
   type (column_carbon_flux), intent(inout) :: col_cf          ! column-level state for p2c
   type(vegetation_phosphorus_flux), intent(inout) :: veg_pf
   type(column_phosphorus_flux)  , intent(inout)   :: col_pf
   type(vegetation_nitrogen_flux), intent(inout)   :: veg_nf
   type(column_nitrogen_flux)    , intent(inout)   :: col_nf

   ! !LOCAL VARIABLES:
   integer :: fc,c,p  ! indices
   real(r8) :: sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8
   real(r8) :: psum1, psum2, nsum1, nsum2


   !$acc enter data create(sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,&
   !$acc                   psum1,psum2)

   !$acc parallel loop gang worker default(present) private(c)
   do fc = 1,numfc
     c = filterc(fc)
     sum1= 0._r8; sum2= 0._r8; sum3= 0._r8
     sum4= 0._r8; sum5= 0._r8; sum6= 0._r8
     sum7= 0._r8; sum8= 0._r8

     psum1 = 0._r8; psum2 = 0._r8;
     !$acc loop vector reduction(+:sum1,sum2,sum3,sum4 ,&
     !$acc  sum5, sum6, sum7, sum8, psum1, psum2)
     do p = col_pp%pfti(c), col_pp%pftf(c)
         if (veg_pp%active(p)) then
             ! carbon flux :
             sum1 = sum1 + veg_cf%gpp(p)     * veg_pp%wtcol(p)
             sum2 = sum2 + veg_cf%ar (p)     * veg_pp%wtcol(p)
             sum3 = sum3 + veg_cf%npp(p)     * veg_pp%wtcol(p)
             sum4 = sum4 + veg_cf%vegfire(p) * veg_pp%wtcol(p)
             sum5 = sum5 + veg_cf%wood_harvestc(p) * veg_pp%wtcol(p)
             sum6 = sum6 + veg_cf%fire_closs(p) * veg_pp%wtcol(p)
             sum7 = sum7 + veg_cf%litfall(p) * veg_pp%wtcol(p)
             sum8 = sum8 + veg_cf%hrv_xsmrpool_to_atm(p) * veg_pp%wtcol(p)
             ! phophorus flux :
             psum1 = psum1 + veg_pf%fire_ploss(p)   * veg_pp%wtcol(p)
             psum2 = psum2 + veg_pf%wood_harvestp(p)* veg_pp%wtcol(p)
             ! nitrogen flux :
             nsum1 = nsum1 + veg_nf%fire_nloss(p)   * veg_pp%wtcol(p)
             nsum2 = nsum2 + veg_nf%wood_harvestn(p)* veg_pp%wtcol(p)
          end if
     end do
     col_cf%gpp(c)  = sum1; col_cf%ar (c)  = sum2
     col_cf%npp(c)  = sum3; col_cf%vegfire(c) = sum4
     col_cf%wood_harvestc(c)  = sum5
     col_cf%fire_closs_p2c(c) = sum6
     col_cf%litfall(c) = sum7; col_cf%hrv_xsmrpool_to_atm(c) = sum8
     ! phosphorus flux :
     col_pf%fire_ploss(c)    = psum1;
     col_pf%wood_harvestp(c) = psum2;
     ! nitrogen flux :
     col_nf%fire_nloss(c) = nsum1;
     col_nf%wood_harvestn(c) = nsum2;
   end do

end subroutine summary_veg_flux_p2c

!-----------------------------------------------------------------------
subroutine veg_nf_summary_acc(this, num_soilp, filter_soilp)
  !
  ! !DESCRIPTION:
  ! Vegetation-level nitrogen flux summary calculations
  !
  ! !ARGUMENTS:
  type(vegetation_nitrogen_flux)             :: this
  integer , intent(in) :: num_soilp
  integer , intent(in) :: filter_soilp(:) 
  !
  integer :: fp, p 

  !$acc parallel loop independent gang vector default(present) 
  do fp = 1, num_soilp 
      p = filter_soilp(fp)
     ! total N deployment (from sminn and retranslocated N pool) (NDEPLOY)
     this%ndeploy(p) = &
          this%sminn_to_npool(p) + &
          this%retransn_to_npool(p)

     ! pft-level wood harvest
     this%wood_harvestn(p) = &
          this%hrv_deadstemn_to_prod10n(p) + &
          this%hrv_deadstemn_to_prod100n(p)
     if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%wood_harvestn(p) = &
          this%wood_harvestn(p) + &
          this%hrv_cropn_to_prod1n(p)
     end if

     ! total pft-level fire N losses
     this%fire_nloss(p) = &
          this%m_leafn_to_fire(p)               + &
          this%m_leafn_storage_to_fire(p)       + &
          this%m_leafn_xfer_to_fire(p)          + &
          this%m_frootn_to_fire(p)              + &
          this%m_frootn_storage_to_fire(p)      + &
          this%m_frootn_xfer_to_fire(p)         + &
          this%m_livestemn_to_fire(p)           + &
          this%m_livestemn_storage_to_fire(p)   + &
          this%m_livestemn_xfer_to_fire(p)      + &
          this%m_deadstemn_to_fire(p)           + &
          this%m_deadstemn_storage_to_fire(p)   + &
          this%m_deadstemn_xfer_to_fire(p)      + &
          this%m_livecrootn_to_fire(p)          + &
          this%m_livecrootn_storage_to_fire(p)  + &
          this%m_livecrootn_xfer_to_fire(p)     + &
          this%m_deadcrootn_to_fire(p)          + &
          this%m_deadcrootn_storage_to_fire(p)  + &
          this%m_deadcrootn_xfer_to_fire(p)     + &
          this%m_retransn_to_fire(p)            + &
          this%m_npool_to_fire(p)

    this%gap_nloss_litter(p) = &
         this%m_leafn_to_litter(p)              + &
         this%m_leafn_storage_to_litter(p)      + &
         this%m_leafn_xfer_to_litter(p)         + &
         this%m_frootn_to_litter(p)             + &
         this%m_frootn_storage_to_litter(p)     + &
         this%m_frootn_xfer_to_litter(p)        + &
         this%m_livestemn_to_litter(p)          + &
         this%m_livestemn_storage_to_litter(p)  + &
         this%m_livestemn_xfer_to_litter(p)     + &
         this%m_deadstemn_to_litter(p)          + &
         this%m_deadstemn_storage_to_litter(p)  + &
         this%m_deadstemn_xfer_to_litter(p)     + &
         this%m_livecrootn_to_litter(p)         + &
         this%m_livecrootn_storage_to_litter(p) + &
         this%m_livecrootn_xfer_to_litter(p)    + &
         this%m_deadcrootn_to_litter(p)         + &
         this%m_deadcrootn_storage_to_litter(p) + &
         this%m_deadcrootn_xfer_to_litter(p)    + &
         this%m_retransn_to_litter(p)           + &
         this%m_npool_to_litter(p)

    this%fire_nloss_litter(p) = &
         this%m_deadstemn_to_litter_fire(p)     + &
         this%m_deadcrootn_to_litter_fire(p)    + &
         this%m_retransn_to_litter_fire(p)      + &
         this%m_npool_to_litter_fire(p)         + &
         this%m_leafn_to_litter_fire(p)         + &
         this%m_frootn_to_litter_fire(p)        + &
         this%m_livestemn_to_litter_fire(p)     + &
         this%m_livecrootn_to_litter_fire(p)    + &
         this%m_leafn_storage_to_litter_fire(p) + &
         this%m_frootn_storage_to_litter_fire(p)       + &
         this%m_livestemn_storage_to_litter_fire(p)    + &
         this%m_deadstemn_storage_to_litter_fire(p)    + &
         this%m_livecrootn_storage_to_litter_fire(p)   + &
         this%m_deadcrootn_storage_to_litter_fire(p)   + &
         this%m_leafn_xfer_to_litter_fire(p)           + &
         this%m_frootn_xfer_to_litter_fire(p)          + &
         this%m_livestemn_xfer_to_litter_fire(p)       + &
         this%m_deadstemn_xfer_to_litter_fire(p)       + &
         this%m_livecrootn_xfer_to_litter_fire(p)      + &
         this%m_deadcrootn_xfer_to_litter_fire(p)

    this%hrv_nloss_litter(p) = &
         this%hrv_retransn_to_litter(p)          + &
         this%hrv_npool_to_litter(p)             + &
         this%hrv_leafn_to_litter(p)             + &
         this%hrv_leafn_storage_to_litter(p)     + &
         this%hrv_leafn_xfer_to_litter(p)        + &
         this%hrv_frootn_to_litter(p)            + &
         this%hrv_frootn_storage_to_litter(p)    + &
         this%hrv_frootn_xfer_to_litter(p)       + &
         this%hrv_livestemn_to_litter(p)         + &
         this%hrv_livestemn_storage_to_litter(p) + &
         this%hrv_livestemn_xfer_to_litter(p)    + &
         this%hrv_deadstemn_storage_to_litter(p) + &
         this%hrv_deadstemn_xfer_to_litter(p)    + &
         this%hrv_livecrootn_to_litter(p)        + &
         this%hrv_livecrootn_storage_to_litter(p)+ &
         this%hrv_livecrootn_xfer_to_litter(p)   + &
         this%hrv_deadcrootn_to_litter(p)        + &
         this%hrv_deadcrootn_storage_to_litter(p)+ &
         this%hrv_deadcrootn_xfer_to_litter(p)
    if (crop_prog) then
       this%sen_nloss_litter(p) = &
           this%livestemn_to_litter(p)            + &
           this%leafn_to_litter(p)                + &
           this%frootn_to_litter(p)
    else
       this%sen_nloss_litter(p) = &
           this%leafn_to_litter(p)                + &
           this%frootn_to_litter(p)
    end if
   end do 

  ! call p2c(bounds, num_soilc, filter_soilc, &
  !      fire_nloss_patch(bounds%begp:bounds%endp)    , &
  !      fire_nloss_col(bounds%begc:bounds%endc))
  !
  ! call p2c(bounds, num_soilc, filter_soilc, &
  !      wood_harvestn_patch(bounds%begp:bounds%endp) , &
  !      wood_harvestn_col(bounds%begc:bounds%endc))

end subroutine veg_nf_summary_acc

!-----------------------------------------------------------------------
subroutine veg_ns_summary_acc(this,num_soilp, filter_soilp)
  !
  ! !DESCRIPTION:
  ! Vegetation-level nitrogen state summary calculations
  !
  ! !ARGUMENTS:
  type(vegetation_nitrogen_state)            :: this
  integer , intent(in) :: num_soilp
  integer , intent(in) :: filter_soilp(:) 
  ! !LOCAL VARIABLES:
  !
  integer :: fp, p 

  !$acc parallel loop independent gang vector default(present) 
  do fp = 1, num_soilp 
      p = filter_soilp(fp)

     ! displayed vegetation nitrogen, excluding storage (DISPVEGN)
     this%dispvegn(p) = &
          this%leafn(p)      + &
          this%frootn(p)     + &
          this%livestemn(p)  + &
          this%deadstemn(p)  + &
          this%livecrootn(p) + &
          this%deadcrootn(p)

    ! stored vegetation nitrogen, including retranslocated N pool (STORVEGN)
    this%storvegn(p) = &
         this%leafn_storage(p)      + &
         this%frootn_storage(p)     + &
         this%livestemn_storage(p)  + &
         this%deadstemn_storage(p)  + &
         this%livecrootn_storage(p) + &
         this%deadcrootn_storage(p) + &
         this%leafn_xfer(p)         + &
         this%frootn_xfer(p)        + &
         this%livestemn_xfer(p)     + &
         this%deadstemn_xfer(p)     + &
         this%livecrootn_xfer(p)    + &
         this%deadcrootn_xfer(p)    + &
         this%npool(p)              + &
         this%retransn(p)

    if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
       this%dispvegn(p) = &
            this%dispvegn(p) + &
            this%grainn(p)

       this%storvegn(p) = &
            this%storvegn(p) + &
            this%grainn_storage(p)     + &
            this%grainn_xfer(p)
    end if

    ! total vegetation nitrogen (TOTVEGN)
    this%totvegn(p) = &
         this%dispvegn(p) + &
         this%storvegn(p)

    ! total pft-level carbon (add pft_ntrunc)
    this%totpftn(p) = &
         this%totvegn(p) + &
         this%ntrunc(p)
   end do 

 ! call p2c(bounds, num_soilc, filter_soilc, &
 !      plant_n_buffer_patch(bounds%begp:bounds%endp)  , &
 !      plant_n_buffer_col(bounds%begc:bounds%endc))
 !
 ! call p2c(bounds, num_soilc, filter_soilc, &
 !      totvegn_patch(bounds%begp:bounds%endp) , &
 !      totvegn_col(bounds%begc:bounds%endc))
 !
 ! call p2c(bounds, num_soilc, filter_soilc, &
 !      totpftn_patch(bounds%begp:bounds%endp) , &
 !      totpftn_col(bounds%begc:bounds%endc))
 !
 ! call p2c(bounds, num_soilc, filter_soilc, &
 !      cropseedn_deficit_patch(bounds%begp:bounds%endp) , &
 !      cropseedn_deficit_col(bounds%begc:bounds%endc))

end subroutine veg_ns_summary_acc

!-----------------------------------------------------------------------
subroutine veg_ps_summary_acc(this, num_soilp, filter_soilp)
  ! !ARGUMENTS:
  type (vegetation_phosphorus_state) :: this
  integer , intent(in) :: num_soilp
  integer , intent(in) :: filter_soilp(:) 
  ! !LOCAL VARIABLES:
  !
  integer :: fp, p 
  !$acc parallel loop independent gang vector default(present) 
  do fp = 1, num_soilp 
      p = filter_soilp(fp)
     ! displayed vegetation phosphorus, excluding storage (DISPVEGN)
     this%dispvegp(p) = &
          this%leafp(p)      + &
          this%frootp(p)     + &
          this%livestemp(p)  + &
          this%deadstemp(p)  + &
          this%livecrootp(p) + &
          this%deadcrootp(p)

    ! stored vegetation phosphorus, including retranslocated N pool (STORVEGN)
    this%storvegp(p) = &
         this%leafp_storage(p)      + &
         this%frootp_storage(p)     + &
         this%livestemp_storage(p)  + &
         this%deadstemp_storage(p)  + &
         this%livecrootp_storage(p) + &
         this%deadcrootp_storage(p) + &
         this%leafp_xfer(p)         + &
         this%frootp_xfer(p)        + &
         this%livestemp_xfer(p)     + &
         this%deadstemp_xfer(p)     + &
         this%livecrootp_xfer(p)    + &
         this%deadcrootp_xfer(p)    + &
         this%ppool(p)              + &
         this%retransp(p)

    if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
      this%dispvegp(p) = &
           this%dispvegp(p) + &
           this%grainp(p)

      this%storvegp(p) = &
           this%storvegp(p) + &
           this%grainp_storage(p)     + &
           this%grainp_xfer(p)
    end if

    ! total vegetation phosphorus (TOTVEGN)
    this%totvegp(p) = this%dispvegp(p) + this%storvegp(p)

    ! total pft-level carbon (add pft_ntrunc)
    this%totpftp(p) = this%totvegp(p) + this%ptrunc(p)

   end do 

 ! call p2c(bounds, num_soilc, filter_soilc, &
 !     totvegp_patch(bounds%begp:bounds%endp)  , &
 !     totvegp_col(bounds%begc:bounds%endc) )
 !
 ! call p2c(bounds, num_soilc, filter_soilc, &
 !     totpftp_patch(bounds%begp:bounds%endp) , &
 !     totpftp_col(bounds%begc:bounds%endc) )
 !
 ! call p2c(bounds, num_soilc, filter_soilc, &
 !     cropseedp_deficit_patch(bounds%begp:bounds%endp) , &
 !     cropseedp_deficit_col(bounds%begc:bounds%endc) )

end subroutine veg_ps_summary_acc



!-----------------------------------------------------------------------
subroutine veg_pf_summary_acc(this,num_soilp, filter_soilp)
  !
  ! !ARGUMENTS:
  type (vegetation_phosphorus_flux) :: this
  integer , intent(in) :: num_soilp
  integer , intent(in) :: filter_soilp(:) 
  ! !LOCAL VARIABLES:
  !
  integer :: fp, p
   
  !$acc parallel loop independent gang vector default(present) 
  do fp = 1, num_soilp 
      p = filter_soilp(fp)
      ! total P deployment (from sminn and retranslocated P pool) (PDEPLOY)
      this%pdeploy(p) = &
          this%sminp_to_ppool(p) + &
          this%retransp_to_ppool(p)

     ! pft-level wood harvest
     this%wood_harvestp(p) = &
          this%hrv_deadstemp_to_prod10p(p) + &
          this%hrv_deadstemp_to_prod100p(p)
     if ( crop_prog .and. veg_pp%itype(p) >= npcropmin )then
          this%wood_harvestp(p) = &
          this%wood_harvestp(p) + &
          this%hrv_cropp_to_prod1p(p)
     end if

     ! total pft-level fire P losses
     this%fire_ploss(p) = &
          this%m_leafp_to_fire(p)               + &
          this%m_leafp_storage_to_fire(p)       + &
          this%m_leafp_xfer_to_fire(p)          + &
          this%m_frootp_to_fire(p)              + &
          this%m_frootp_storage_to_fire(p)      + &
          this%m_frootp_xfer_to_fire(p)         + &
          this%m_livestemp_to_fire(p)           + &
          this%m_livestemp_storage_to_fire(p)   + &
          this%m_livestemp_xfer_to_fire(p)      + &
          this%m_deadstemp_to_fire(p)           + &
          this%m_deadstemp_storage_to_fire(p)   + &
          this%m_deadstemp_xfer_to_fire(p)      + &
          this%m_livecrootp_to_fire(p)          + &
          this%m_livecrootp_storage_to_fire(p)  + &
          this%m_livecrootp_xfer_to_fire(p)     + &
          this%m_deadcrootp_to_fire(p)          + &
          this%m_deadcrootp_storage_to_fire(p)  + &
          this%m_deadcrootp_xfer_to_fire(p)     + &
          this%m_retransp_to_fire(p)            + &
          this%m_ppool_to_fire(p)

    this%gap_ploss_litter(p) = &
         this%m_leafp_to_litter(p)              + &
         this%m_leafp_storage_to_litter(p)      + &
         this%m_leafp_xfer_to_litter(p)         + &
         this%m_frootp_to_litter(p)             + &
         this%m_frootp_storage_to_litter(p)     + &
         this%m_frootp_xfer_to_litter(p)        + &
         this%m_livestemp_to_litter(p)          + &
         this%m_livestemp_storage_to_litter(p)  + &
         this%m_livestemp_xfer_to_litter(p)     + &
         this%m_deadstemp_to_litter(p)          + &
         this%m_deadstemp_storage_to_litter(p)  + &
         this%m_deadstemp_xfer_to_litter(p)     + &
         this%m_livecrootp_to_litter(p)         + &
         this%m_livecrootp_storage_to_litter(p) + &
         this%m_livecrootp_xfer_to_litter(p)    + &
         this%m_deadcrootp_to_litter(p)         + &
         this%m_deadcrootp_storage_to_litter(p) + &
         this%m_deadcrootp_xfer_to_litter(p)    + &
         this%m_retransp_to_litter(p)           + &
         this%m_ppool_to_litter(p)

    this%fire_ploss_litter(p) = &
         this%m_deadstemp_to_litter_fire(p)     + &
         this%m_deadcrootp_to_litter_fire(p)    + &
         this%m_retransp_to_litter_fire(p)      + &
         this%m_ppool_to_litter_fire(p)         + &
         this%m_leafp_to_litter_fire(p)         + &
         this%m_frootp_to_litter_fire(p)        + &
         this%m_livestemp_to_litter_fire(p)     + &
         this%m_livecrootp_to_litter_fire(p)    + &
         this%m_leafp_storage_to_litter_fire(p) + &
         this%m_frootp_storage_to_litter_fire(p)       + &
         this%m_livestemp_storage_to_litter_fire(p)    + &
         this%m_deadstemp_storage_to_litter_fire(p)    + &
         this%m_livecrootp_storage_to_litter_fire(p)   + &
         this%m_deadcrootp_storage_to_litter_fire(p)   + &
         this%m_leafp_xfer_to_litter_fire(p)           + &
         this%m_frootp_xfer_to_litter_fire(p)          + &
         this%m_livestemp_xfer_to_litter_fire(p)       + &
         this%m_deadstemp_xfer_to_litter_fire(p)       + &
         this%m_livecrootp_xfer_to_litter_fire(p)      + &
         this%m_deadcrootp_xfer_to_litter_fire(p)

    this%hrv_ploss_litter(p) = &
         this%hrv_retransp_to_litter(p)         + &
         this%hrv_ppool_to_litter(p)            + &
         this%hrv_leafp_to_litter(p)            + &
         this%hrv_leafp_storage_to_litter(p)    + &
         this%hrv_leafp_xfer_to_litter(p)       + &
         this%hrv_frootp_to_litter(p)           + &
         this%hrv_frootp_storage_to_litter(p)   + &
         this%hrv_frootp_xfer_to_litter(p)      + &
         this%hrv_livestemp_to_litter(p)        + &
         this%hrv_livestemp_storage_to_litter(p)+ &
         this%hrv_livestemp_xfer_to_litter(p)   + &
         this%hrv_deadstemp_storage_to_litter(p)+ &
         this%hrv_deadstemp_xfer_to_litter(p)   + &
         this%hrv_livecrootp_to_litter(p)       + &
         this%hrv_livecrootp_storage_to_litter(p)+ &
         this%hrv_livecrootp_xfer_to_litter(p)  + &
         this%hrv_deadcrootp_to_litter(p)       + &
         this%hrv_deadcrootp_storage_to_litter(p)+ &
         this%hrv_deadcrootp_xfer_to_litter(p)

    if (crop_prog) then
       this%sen_ploss_litter(p) = &
           this%livestemp_to_litter(p)            + &
           this%leafp_to_litter(p)                + &
           this%frootp_to_litter(p)
    else
       this%sen_ploss_litter(p) = &
           this%leafp_to_litter(p)                + &
           this%frootp_to_litter(p)
    end if
   end do 
  ! call p2c(bounds, num_soilc, filter_soilc, &
  !      fire_ploss_patch(bounds%begp:bounds%endp)     , &
  !      fire_ploss_col(bounds%begc:bounds%endc) )
  !
  ! call p2c(bounds, num_soilc, filter_soilc, &
  !      wood_harvestp_patch(bounds%begp:bounds%endp) , &
  !      wood_harvestp_col(bounds%begc:bounds%endc))

end subroutine veg_pf_summary_acc


!-----------------------------------------------------------------------
subroutine veg_cf_setvalues_acc ( this, num_soilp,filter_soilp, value_patch)
  !
  ! !DESCRIPTION:
  ! Set vegetation-level carbon fluxes
  ! !ARGUMENTS:
  type (vegetation_carbon_flux) :: this
  integer , intent(in) :: num_soilp 
  integer , intent(in) :: filter_soilp(:) 
  real(r8), intent(in) :: value_patch

  integer :: fp, i
  !$acc data copyin(value_patch)
  if(.not.use_fates) then
   !$acc parallel loop independent gang vector default(present) 
   do fp = 1, num_soilp 
      i = filter_soilp(fp) 
      this%m_leafc_to_litter(i)                   = value_patch
      this%m_frootc_to_litter(i)                  = value_patch
      this%m_leafc_storage_to_litter(i)           = value_patch
      this%m_frootc_storage_to_litter(i)          = value_patch
      this%m_livestemc_storage_to_litter(i)       = value_patch
      this%m_deadstemc_storage_to_litter(i)       = value_patch
      this%m_livecrootc_storage_to_litter(i)      = value_patch
      this%m_deadcrootc_storage_to_litter(i)      = value_patch
      this%m_leafc_xfer_to_litter(i)              = value_patch
      this%m_frootc_xfer_to_litter(i)             = value_patch
      this%m_livestemc_xfer_to_litter(i)          = value_patch
      this%m_deadstemc_xfer_to_litter(i)          = value_patch
      this%m_livecrootc_xfer_to_litter(i)         = value_patch
      this%m_deadcrootc_xfer_to_litter(i)         = value_patch
      this%m_livestemc_to_litter(i)               = value_patch
      this%m_deadstemc_to_litter(i)               = value_patch
      this%m_livecrootc_to_litter(i)              = value_patch
      this%m_deadcrootc_to_litter(i)              = value_patch
      this%m_gresp_storage_to_litter(i)           = value_patch
      this%m_gresp_xfer_to_litter(i)              = value_patch
      this%m_cpool_to_litter(i)                   = value_patch
      this%hrv_leafc_to_litter(i)                 = value_patch
      this%hrv_leafc_storage_to_litter(i)         = value_patch
      this%hrv_leafc_xfer_to_litter(i)            = value_patch
      this%hrv_frootc_to_litter(i)                = value_patch
      this%hrv_frootc_storage_to_litter(i)        = value_patch
      this%hrv_frootc_xfer_to_litter(i)           = value_patch
      this%hrv_livestemc_to_litter(i)             = value_patch
      this%hrv_livestemc_storage_to_litter(i)     = value_patch
      this%hrv_livestemc_xfer_to_litter(i)        = value_patch
      this%hrv_deadstemc_to_prod10c(i)            = value_patch
      this%hrv_deadstemc_to_prod100c(i)           = value_patch
      this%hrv_leafc_to_prod1c(i)                 = value_patch
      this%hrv_livestemc_to_prod1c(i)             = value_patch
      this%hrv_grainc_to_prod1c(i)                = value_patch
      this%hrv_cropc_to_prod1c(i)                 = value_patch
      this%hrv_deadstemc_storage_to_litter(i)     = value_patch
      this%hrv_deadstemc_xfer_to_litter(i)        = value_patch
      this%hrv_livecrootc_to_litter(i)            = value_patch
      this%hrv_livecrootc_storage_to_litter(i)    = value_patch
      this%hrv_livecrootc_xfer_to_litter(i)       = value_patch
      this%hrv_deadcrootc_to_litter(i)            = value_patch
      this%hrv_deadcrootc_storage_to_litter(i)    = value_patch
      this%hrv_deadcrootc_xfer_to_litter(i)       = value_patch
      this%hrv_gresp_storage_to_litter(i)         = value_patch
      this%hrv_gresp_xfer_to_litter(i)            = value_patch
      this%hrv_xsmrpool_to_atm(i)                 = value_patch
      this%hrv_cpool_to_litter(i)                 = value_patch

      this%m_leafc_to_fire(i)                     = value_patch
      this%m_leafc_storage_to_fire(i)             = value_patch
      this%m_leafc_xfer_to_fire(i)                = value_patch
      this%m_livestemc_to_fire(i)                 = value_patch
      this%m_livestemc_storage_to_fire(i)         = value_patch
      this%m_livestemc_xfer_to_fire(i)            = value_patch
      this%m_deadstemc_to_fire(i)                 = value_patch
      this%m_deadstemc_storage_to_fire(i)         = value_patch
      this%m_deadstemc_xfer_to_fire(i)            = value_patch
      this%m_frootc_to_fire(i)                    = value_patch
      this%m_frootc_storage_to_fire(i)            = value_patch
      this%m_frootc_xfer_to_fire(i)               = value_patch
      this%m_livecrootc_to_fire(i)                = value_patch
      this%m_livecrootc_storage_to_fire(i)        = value_patch
      this%m_livecrootc_xfer_to_fire(i)           = value_patch
      this%m_deadcrootc_to_fire(i)                = value_patch
      this%m_deadcrootc_storage_to_fire(i)        = value_patch
      this%m_deadcrootc_xfer_to_fire(i)           = value_patch
      this%m_gresp_storage_to_fire(i)             = value_patch
      this%m_gresp_xfer_to_fire(i)                = value_patch
      this%m_cpool_to_fire(i)                     = value_patch

      this%m_leafc_to_litter_fire(i)              = value_patch
      this%m_leafc_storage_to_litter_fire(i)      = value_patch
      this%m_leafc_xfer_to_litter_fire(i)         = value_patch
      this%m_livestemc_to_litter_fire(i)          = value_patch
      this%m_livestemc_storage_to_litter_fire(i)  = value_patch
      this%m_livestemc_xfer_to_litter_fire(i)     = value_patch
      this%m_livestemc_to_deadstemc_fire(i)       = value_patch
      this%m_deadstemc_to_litter_fire(i)          = value_patch
      this%m_deadstemc_storage_to_litter_fire(i)  = value_patch
      this%m_deadstemc_xfer_to_litter_fire(i)     = value_patch
      this%m_frootc_to_litter_fire(i)             = value_patch
      this%m_frootc_storage_to_litter_fire(i)     = value_patch
      this%m_frootc_xfer_to_litter_fire(i)        = value_patch
      this%m_livecrootc_to_litter_fire(i)         = value_patch
      this%m_livecrootc_storage_to_litter_fire(i) = value_patch
      this%m_livecrootc_xfer_to_litter_fire(i)    = value_patch
      this%m_livecrootc_to_deadcrootc_fire(i)     = value_patch
      this%m_deadcrootc_to_litter_fire(i)         = value_patch
      this%m_deadcrootc_storage_to_litter_fire(i) = value_patch
      this%m_deadcrootc_xfer_to_litter_fire(i)    = value_patch
      this%m_gresp_storage_to_litter_fire(i)      = value_patch
      this%m_gresp_xfer_to_litter_fire(i)         = value_patch
      this%m_cpool_to_litter_fire(i)              = value_patch
      this%leafc_xfer_to_leafc(i)                 = value_patch
      this%frootc_xfer_to_frootc(i)               = value_patch
      this%livestemc_xfer_to_livestemc(i)         = value_patch
      this%deadstemc_xfer_to_deadstemc(i)         = value_patch
      this%livecrootc_xfer_to_livecrootc(i)       = value_patch
      this%deadcrootc_xfer_to_deadcrootc(i)       = value_patch
      this%leafc_to_litter(i)                     = value_patch
      this%frootc_to_litter(i)                    = value_patch
      this%leaf_mr(i)                             = value_patch
      this%froot_mr(i)                            = value_patch
      this%livestem_mr(i)                         = value_patch
      this%livecroot_mr(i)                        = value_patch
      this%grain_mr(i)                            = value_patch
      this%leaf_curmr(i)                          = value_patch
      this%froot_curmr(i)                         = value_patch
      this%livestem_curmr(i)                      = value_patch
      this%livecroot_curmr(i)                     = value_patch
      this%grain_curmr(i)                         = value_patch
      this%leaf_xsmr(i)                           = value_patch
      this%froot_xsmr(i)                          = value_patch
      this%livestem_xsmr(i)                       = value_patch
      this%livecroot_xsmr(i)                      = value_patch
      this%grain_xsmr(i)                          = value_patch
      this%xr(i)                                  = value_patch
      this%psnsun_to_cpool(i)                     = value_patch
      this%psnshade_to_cpool(i)                   = value_patch
      this%cpool_to_xsmrpool(i)                   = value_patch
      this%cpool_to_leafc(i)                      = value_patch
      this%cpool_to_leafc_storage(i)              = value_patch
      this%cpool_to_frootc(i)                     = value_patch
      this%cpool_to_frootc_storage(i)             = value_patch
      this%cpool_to_livestemc(i)                  = value_patch
      this%cpool_to_livestemc_storage(i)          = value_patch
      this%cpool_to_deadstemc(i)                  = value_patch
      this%cpool_to_deadstemc_storage(i)          = value_patch
      this%cpool_to_livecrootc(i)                 = value_patch
      this%cpool_to_livecrootc_storage(i)         = value_patch
      this%cpool_to_deadcrootc(i)                 = value_patch
      this%cpool_to_deadcrootc_storage(i)         = value_patch
      this%cpool_to_gresp_storage(i)              = value_patch
      this%cpool_leaf_gr(i)                       = value_patch
      this%cpool_leaf_storage_gr(i)               = value_patch
      this%transfer_leaf_gr(i)                    = value_patch
      this%cpool_froot_gr(i)                      = value_patch
      this%cpool_froot_storage_gr(i)              = value_patch
      this%transfer_froot_gr(i)                   = value_patch
      this%cpool_livestem_gr(i)                   = value_patch
      this%cpool_livestem_storage_gr(i)           = value_patch
      this%transfer_livestem_gr(i)                = value_patch
      this%cpool_deadstem_gr(i)                   = value_patch
      this%cpool_deadstem_storage_gr(i)           = value_patch
      this%transfer_deadstem_gr(i)                = value_patch
      this%cpool_livecroot_gr(i)                  = value_patch
      this%cpool_livecroot_storage_gr(i)          = value_patch
      this%transfer_livecroot_gr(i)               = value_patch
      this%cpool_deadcroot_gr(i)                  = value_patch
      this%cpool_deadcroot_storage_gr(i)          = value_patch
      this%transfer_deadcroot_gr(i)               = value_patch
      this%leafc_storage_to_xfer(i)               = value_patch
      this%frootc_storage_to_xfer(i)              = value_patch
      this%livestemc_storage_to_xfer(i)           = value_patch
      this%deadstemc_storage_to_xfer(i)           = value_patch
      this%livecrootc_storage_to_xfer(i)          = value_patch
      this%deadcrootc_storage_to_xfer(i)          = value_patch
      this%gresp_storage_to_xfer(i)               = value_patch
      this%livestemc_to_deadstemc(i)              = value_patch
      this%livecrootc_to_deadcrootc(i)            = value_patch
      this%gpp(i)                                 = value_patch
      this%gpp_before_downreg(i)                  = value_patch
      this%mr(i)                                  = value_patch
      this%current_gr(i)                          = value_patch
      this%transfer_gr(i)                         = value_patch
      this%storage_gr(i)                          = value_patch
      this%gr(i)                                  = value_patch
      this%ar(i)                                  = value_patch
      this%rr(i)                                  = value_patch
      this%npp(i)                                 = value_patch
      this%agnpp(i)                               = value_patch
      this%bgnpp(i)                               = value_patch
      this%agwdnpp(i)                             = value_patch
      this%litfall(i)                             = value_patch
      this%vegfire(i)                             = value_patch
      this%wood_harvestc(i)                       = value_patch
      this%cinputs(i)                             = value_patch
      this%coutputs(i)                            = value_patch
      this%fire_closs(i)                          = value_patch
      this%frootc_alloc(i)                        = value_patch
      this%frootc_loss(i)                         = value_patch
      this%leafc_alloc(i)                         = value_patch
      this%leafc_loss(i)                          = value_patch
      this%woodc_alloc(i)                         = value_patch
      this%woodc_loss(i)                          = value_patch
      this%xsmrpool_turnover(i)                   = value_patch
   end do 
  end if !(.not.use_fates)

  if ( crop_prog )then
      !$acc parallel loop independent gang vector default(present) 
      do fp = 1, num_soilp 
         i = filter_soilp(fp) 
        this%xsmrpool_to_atm(i)         = value_patch
        this%livestemc_to_litter(i)     = value_patch
        this%grainc_to_food(i)          = value_patch
        this%grainc_xfer_to_grainc(i)   = value_patch
        this%cpool_to_grainc(i)         = value_patch
        this%cpool_to_grainc_storage(i) = value_patch
        this%cpool_grain_gr(i)          = value_patch
        this%cpool_grain_storage_gr(i)  = value_patch
        this%transfer_grain_gr(i)       = value_patch
        this%grainc_storage_to_xfer(i)  = value_patch
        this%crop_seedc_to_leaf(i)      = value_patch
      end do 
  end if

  !$acc end data

end subroutine veg_cf_setvalues_acc

!-----------------------------------------------------------------------
subroutine veg_nf_setvalues_acc( this, num_soilp,filter_soilp, value_patch)
  !
  ! !DESCRIPTION:
  ! Set vegetation-level nitrogen fluxes
  ! !ARGUMENTS:
  type (vegetation_nitrogen_flux) :: this
  integer , intent(in) :: num_soilp
  integer , intent(in) :: filter_soilp(:)
  real(r8), intent(in) :: value_patch
  !
  integer :: fp, i
  !$acc data copyin(value_patch)
  !------------------------------------------------------------------------
  !$acc parallel loop independent gang vector default(present) 
  do fp = 1, num_soilp 
      i = filter_soilp(fp) 
     this%m_leafn_to_litter(i)                   = value_patch
     this%m_frootn_to_litter(i)                  = value_patch
     this%m_leafn_storage_to_litter(i)           = value_patch
     this%m_frootn_storage_to_litter(i)          = value_patch
     this%m_livestemn_storage_to_litter(i)       = value_patch
     this%m_deadstemn_storage_to_litter(i)       = value_patch
     this%m_livecrootn_storage_to_litter(i)      = value_patch
     this%m_deadcrootn_storage_to_litter(i)      = value_patch
     this%m_leafn_xfer_to_litter(i)              = value_patch
     this%m_frootn_xfer_to_litter(i)             = value_patch
     this%m_livestemn_xfer_to_litter(i)          = value_patch
     this%m_deadstemn_xfer_to_litter(i)          = value_patch
     this%m_livecrootn_xfer_to_litter(i)         = value_patch
     this%m_deadcrootn_xfer_to_litter(i)         = value_patch
     this%m_livestemn_to_litter(i)               = value_patch
     this%m_deadstemn_to_litter(i)               = value_patch
     this%m_livecrootn_to_litter(i)              = value_patch
     this%m_deadcrootn_to_litter(i)              = value_patch
     this%m_retransn_to_litter(i)                = value_patch
     this%m_npool_to_litter(i)                   = value_patch
     this%hrv_leafn_to_litter(i)                 = value_patch
     this%hrv_frootn_to_litter(i)                = value_patch
     this%hrv_leafn_storage_to_litter(i)         = value_patch
     this%hrv_frootn_storage_to_litter(i)        = value_patch
     this%hrv_livestemn_storage_to_litter(i)     = value_patch
     this%hrv_deadstemn_storage_to_litter(i)     = value_patch
     this%hrv_livecrootn_storage_to_litter(i)    = value_patch
     this%hrv_deadcrootn_storage_to_litter(i)    = value_patch
     this%hrv_leafn_xfer_to_litter(i)            = value_patch
     this%hrv_frootn_xfer_to_litter(i)           = value_patch
     this%hrv_livestemn_xfer_to_litter(i)        = value_patch
     this%hrv_deadstemn_xfer_to_litter(i)        = value_patch
     this%hrv_livecrootn_xfer_to_litter(i)       = value_patch
     this%hrv_deadcrootn_xfer_to_litter(i)       = value_patch
     this%hrv_livestemn_to_litter(i)             = value_patch
     this%hrv_deadstemn_to_prod10n(i)            = value_patch
     this%hrv_deadstemn_to_prod100n(i)           = value_patch

     this%hrv_leafn_to_prod1n(i)                 = value_patch
     this%hrv_livestemn_to_prod1n(i)             = value_patch
     this%hrv_grainn_to_prod1n(i)                = value_patch
     this%hrv_cropn_to_prod1n(i)                 = value_patch
     this%hrv_livecrootn_to_litter(i)            = value_patch
     this%hrv_deadcrootn_to_litter(i)            = value_patch
     this%hrv_retransn_to_litter(i)              = value_patch
     this%hrv_npool_to_litter(i)                 = value_patch

     this%m_leafn_to_fire(i)                     = value_patch
     this%m_leafn_storage_to_fire(i)             = value_patch
     this%m_leafn_xfer_to_fire(i)                = value_patch
     this%m_livestemn_to_fire(i)                 = value_patch
     this%m_livestemn_storage_to_fire(i)         = value_patch
     this%m_livestemn_xfer_to_fire(i)            = value_patch
     this%m_deadstemn_to_fire(i)                 = value_patch
     this%m_deadstemn_storage_to_fire(i)         = value_patch
     this%m_deadstemn_xfer_to_fire(i)            = value_patch
     this%m_frootn_to_fire(i)                    = value_patch
     this%m_frootn_storage_to_fire(i)            = value_patch
     this%m_frootn_xfer_to_fire(i)               = value_patch
     this%m_livecrootn_to_fire(i)                = value_patch
     this%m_livecrootn_storage_to_fire(i)        = value_patch
     this%m_livecrootn_xfer_to_fire(i)           = value_patch
     this%m_deadcrootn_to_fire(i)                = value_patch
     this%m_deadcrootn_storage_to_fire(i)        = value_patch
     this%m_deadcrootn_xfer_to_fire(i)           = value_patch
     this%m_retransn_to_fire(i)                  = value_patch
     this%m_npool_to_fire(i)                     = value_patch

     this%m_leafn_to_litter_fire(i)              = value_patch
     this%m_leafn_storage_to_litter_fire(i)      = value_patch
     this%m_leafn_xfer_to_litter_fire(i)         = value_patch
     this%m_livestemn_to_litter_fire(i)          = value_patch
     this%m_livestemn_storage_to_litter_fire(i)  = value_patch
     this%m_livestemn_xfer_to_litter_fire(i)     = value_patch
     this%m_livestemn_to_deadstemn_fire(i)       = value_patch
     this%m_deadstemn_to_litter_fire(i)          = value_patch
     this%m_deadstemn_storage_to_litter_fire(i)  = value_patch
     this%m_deadstemn_xfer_to_litter_fire(i)     = value_patch
     this%m_frootn_to_litter_fire(i)             = value_patch
     this%m_frootn_storage_to_litter_fire(i)     = value_patch
     this%m_frootn_xfer_to_litter_fire(i)        = value_patch
     this%m_livecrootn_to_litter_fire(i)         = value_patch
     this%m_livecrootn_storage_to_litter_fire(i) = value_patch
     this%m_livecrootn_xfer_to_litter_fire(i)    = value_patch
     this%m_livecrootn_to_deadcrootn_fire(i)     = value_patch
     this%m_deadcrootn_to_litter_fire(i)         = value_patch
     this%m_deadcrootn_storage_to_litter_fire(i) = value_patch
     this%m_deadcrootn_xfer_to_litter_fire(i)    = value_patch
     this%m_retransn_to_litter_fire(i)           = value_patch
     this%m_npool_to_litter_fire(i)              = value_patch

     this%leafn_xfer_to_leafn(i)                 = value_patch
     this%frootn_xfer_to_frootn(i)               = value_patch
     this%livestemn_xfer_to_livestemn(i)         = value_patch
     this%deadstemn_xfer_to_deadstemn(i)         = value_patch
     this%livecrootn_xfer_to_livecrootn(i)       = value_patch
     this%deadcrootn_xfer_to_deadcrootn(i)       = value_patch
     this%leafn_to_litter(i)                     = value_patch
     this%leafn_to_retransn(i)                   = value_patch
     this%frootn_to_litter(i)                    = value_patch
     this%retransn_to_npool(i)                   = value_patch
     this%sminn_to_npool(i)                      = value_patch
     this%npool_to_leafn(i)                      = value_patch
     this%npool_to_leafn_storage(i)              = value_patch
     this%npool_to_frootn(i)                     = value_patch
     this%npool_to_frootn_storage(i)             = value_patch
     this%npool_to_livestemn(i)                  = value_patch
     this%npool_to_livestemn_storage(i)          = value_patch
     this%npool_to_deadstemn(i)                  = value_patch
     this%npool_to_deadstemn_storage(i)          = value_patch
     this%npool_to_livecrootn(i)                 = value_patch
     this%npool_to_livecrootn_storage(i)         = value_patch
     this%npool_to_deadcrootn(i)                 = value_patch
     this%npool_to_deadcrootn_storage(i)         = value_patch
     this%leafn_storage_to_xfer(i)               = value_patch
     this%frootn_storage_to_xfer(i)              = value_patch
     this%livestemn_storage_to_xfer(i)           = value_patch
     this%deadstemn_storage_to_xfer(i)           = value_patch
     this%livecrootn_storage_to_xfer(i)          = value_patch
     this%deadcrootn_storage_to_xfer(i)          = value_patch
     this%livestemn_to_deadstemn(i)              = value_patch
     this%livestemn_to_retransn(i)               = value_patch
     this%livecrootn_to_deadcrootn(i)            = value_patch
     this%livecrootn_to_retransn(i)              = value_patch
     this%ndeploy(i)                             = value_patch
     this%wood_harvestn(i)                       = value_patch
     this%fire_nloss(i)                          = value_patch
     this%nfix_to_plantn(i)                      = value_patch
     this%gap_nloss_litter(i)                    = value_patch
     this%fire_nloss_litter(i)                   = value_patch
     this%hrv_nloss_litter(i)                    = value_patch
     this%sen_nloss_litter(i)                    = value_patch
     this%crop_seedn_to_leaf(i)                  = value_patch
     this%livestemn_to_litter(i)                 = value_patch

      if ( crop_prog )then
        this%grainn_to_food(i)                   = value_patch
        this%grainn_xfer_to_grainn(i)            = value_patch
        this%npool_to_grainn(i)                  = value_patch
        this%npool_to_grainn_storage(i)          = value_patch
        this%grainn_storage_to_xfer(i)           = value_patch
        this%soyfixn(i)                          = value_patch
        this%frootn_to_retransn(i)               = value_patch
      end if
   end do 
   !$acc end data 

end subroutine veg_nf_setvalues_acc

!-----------------------------------------------------------------------
subroutine veg_pf_setvalues_acc ( this, num_soilp, filter_soilp, value_patch)
  !
  ! !DESCRIPTION:
  ! Set phosphorus flux variables
  ! !ARGUMENTS:
  type (vegetation_phosphorus_flux) :: this
  integer, intent(in) :: num_soilp 
  integer, intent(in) :: filter_soilp(:)  
  real(r8), intent(in) :: value_patch

  integer :: fp, i
  !$acc data copyin(value_patch)

  !------------------------------------------------------------------------
  !$acc parallel loop independent gang vector default(present) 
  do fp = 1, num_soilp 
      i = filter_soilp(fp) 
     this%m_leafp_to_litter(i)                   = value_patch
     this%m_frootp_to_litter(i)                  = value_patch
     this%m_leafp_storage_to_litter(i)           = value_patch
     this%m_frootp_storage_to_litter(i)          = value_patch
     this%m_livestemp_storage_to_litter(i)       = value_patch
     this%m_deadstemp_storage_to_litter(i)       = value_patch
     this%m_livecrootp_storage_to_litter(i)      = value_patch
     this%m_deadcrootp_storage_to_litter(i)      = value_patch
     this%m_leafp_xfer_to_litter(i)              = value_patch
     this%m_frootp_xfer_to_litter(i)             = value_patch
     this%m_livestemp_xfer_to_litter(i)          = value_patch
     this%m_deadstemp_xfer_to_litter(i)          = value_patch
     this%m_livecrootp_xfer_to_litter(i)         = value_patch
     this%m_deadcrootp_xfer_to_litter(i)         = value_patch
     this%m_livestemp_to_litter(i)               = value_patch
     this%m_deadstemp_to_litter(i)               = value_patch
     this%m_livecrootp_to_litter(i)              = value_patch
     this%m_deadcrootp_to_litter(i)              = value_patch
     this%m_retransp_to_litter(i)                = value_patch
     this%m_ppool_to_litter(i)                   = value_patch
     this%hrv_leafp_to_litter(i)                 = value_patch
     this%hrv_frootp_to_litter(i)                = value_patch
     this%hrv_leafp_storage_to_litter(i)         = value_patch
     this%hrv_frootp_storage_to_litter(i)        = value_patch
     this%hrv_livestemp_storage_to_litter(i)     = value_patch
     this%hrv_deadstemp_storage_to_litter(i)     = value_patch
     this%hrv_livecrootp_storage_to_litter(i)    = value_patch
     this%hrv_deadcrootp_storage_to_litter(i)    = value_patch
     this%hrv_leafp_xfer_to_litter(i)            = value_patch
     this%hrv_frootp_xfer_to_litter(i)           = value_patch
     this%hrv_livestemp_xfer_to_litter(i)        = value_patch
     this%hrv_deadstemp_xfer_to_litter(i)        = value_patch
     this%hrv_livecrootp_xfer_to_litter(i)       = value_patch
     this%hrv_deadcrootp_xfer_to_litter(i)       = value_patch
     this%hrv_livestemp_to_litter(i)             = value_patch
     this%hrv_deadstemp_to_prod10p(i)            = value_patch
     this%hrv_deadstemp_to_prod100p(i)           = value_patch
     this%hrv_leafp_to_prod1p(i)                 = value_patch
     this%hrv_livestemp_to_prod1p(i)             = value_patch
     this%hrv_grainp_to_prod1p(i)                = value_patch
     this%hrv_cropp_to_prod1p(i)                 = value_patch
     this%hrv_livecrootp_to_litter(i)            = value_patch
     this%hrv_deadcrootp_to_litter(i)            = value_patch
     this%hrv_retransp_to_litter(i)              = value_patch
     this%hrv_ppool_to_litter(i)                 = value_patch

     this%m_leafp_to_fire(i)                     = value_patch
     this%m_leafp_storage_to_fire(i)             = value_patch
     this%m_leafp_xfer_to_fire(i)                = value_patch
     this%m_livestemp_to_fire(i)                 = value_patch
     this%m_livestemp_storage_to_fire(i)         = value_patch
     this%m_livestemp_xfer_to_fire(i)            = value_patch
     this%m_deadstemp_to_fire(i)                 = value_patch
     this%m_deadstemp_storage_to_fire(i)         = value_patch
     this%m_deadstemp_xfer_to_fire(i)            = value_patch
     this%m_frootp_to_fire(i)                    = value_patch
     this%m_frootp_storage_to_fire(i)            = value_patch
     this%m_frootp_xfer_to_fire(i)               = value_patch
     this%m_livecrootp_to_fire(i)                = value_patch
     this%m_livecrootp_storage_to_fire(i)        = value_patch
     this%m_livecrootp_xfer_to_fire(i)           = value_patch
     this%m_deadcrootp_to_fire(i)                = value_patch
     this%m_deadcrootp_storage_to_fire(i)        = value_patch
     this%m_deadcrootp_xfer_to_fire(i)           = value_patch
     this%m_retransp_to_fire(i)                  = value_patch
     this%m_ppool_to_fire(i)                     = value_patch

     this%m_leafp_to_litter_fire(i)              = value_patch
     this%m_leafp_storage_to_litter_fire(i)      = value_patch
     this%m_leafp_xfer_to_litter_fire(i)         = value_patch
     this%m_livestemp_to_litter_fire(i)          = value_patch
     this%m_livestemp_storage_to_litter_fire(i)  = value_patch
     this%m_livestemp_xfer_to_litter_fire(i)     = value_patch
     this%m_livestemp_to_deadstemp_fire(i)       = value_patch
     this%m_deadstemp_to_litter_fire(i)          = value_patch
     this%m_deadstemp_storage_to_litter_fire(i)  = value_patch
     this%m_deadstemp_xfer_to_litter_fire(i)     = value_patch
     this%m_frootp_to_litter_fire(i)             = value_patch
     this%m_frootp_storage_to_litter_fire(i)     = value_patch
     this%m_frootp_xfer_to_litter_fire(i)        = value_patch
     this%m_livecrootp_to_litter_fire(i)         = value_patch
     this%m_livecrootp_storage_to_litter_fire(i) = value_patch
     this%m_livecrootp_xfer_to_litter_fire(i)    = value_patch
     this%m_livecrootp_to_deadcrootp_fire(i)     = value_patch
     this%m_deadcrootp_to_litter_fire(i)         = value_patch
     this%m_deadcrootp_storage_to_litter_fire(i) = value_patch
     this%m_deadcrootp_xfer_to_litter_fire(i)    = value_patch
     this%m_retransp_to_litter_fire(i)           = value_patch
     this%m_ppool_to_litter_fire(i)              = value_patch

     this%leafp_xfer_to_leafp(i)                 = value_patch
     this%frootp_xfer_to_frootp(i)               = value_patch
     this%livestemp_xfer_to_livestemp(i)         = value_patch
     this%deadstemp_xfer_to_deadstemp(i)         = value_patch
     this%livecrootp_xfer_to_livecrootp(i)       = value_patch
     this%deadcrootp_xfer_to_deadcrootp(i)       = value_patch
     this%leafp_to_litter(i)                     = value_patch
     this%leafp_to_retransp(i)                   = value_patch
     this%frootp_to_litter(i)                    = value_patch
     this%retransp_to_ppool(i)                   = value_patch
     this%sminp_to_ppool(i)                      = value_patch
     this%ppool_to_leafp(i)                      = value_patch
     this%ppool_to_leafp_storage(i)              = value_patch
     this%ppool_to_frootp(i)                     = value_patch
     this%ppool_to_frootp_storage(i)             = value_patch
     this%ppool_to_livestemp(i)                  = value_patch
     this%ppool_to_livestemp_storage(i)          = value_patch
     this%ppool_to_deadstemp(i)                  = value_patch
     this%ppool_to_deadstemp_storage(i)          = value_patch
     this%ppool_to_livecrootp(i)                 = value_patch
     this%ppool_to_livecrootp_storage(i)         = value_patch
     this%ppool_to_deadcrootp(i)                 = value_patch
     this%ppool_to_deadcrootp_storage(i)         = value_patch
     this%leafp_storage_to_xfer(i)               = value_patch
     this%frootp_storage_to_xfer(i)              = value_patch
     this%livestemp_storage_to_xfer(i)           = value_patch
     this%deadstemp_storage_to_xfer(i)           = value_patch
     this%livecrootp_storage_to_xfer(i)          = value_patch
     this%deadcrootp_storage_to_xfer(i)          = value_patch
     this%livestemp_to_deadstemp(i)              = value_patch
     this%livestemp_to_retransp(i)               = value_patch
     this%livecrootp_to_deadcrootp(i)            = value_patch
     this%livecrootp_to_retransp(i)              = value_patch
     this%pdeploy(i)                             = value_patch
     this%wood_harvestp(i)                       = value_patch
     this%fire_ploss(i)                          = value_patch
     this%biochem_pmin_to_plant(i)               = value_patch
     this%gap_ploss_litter(i)                    = value_patch
     this%fire_ploss_litter(i)                   = value_patch
     this%hrv_ploss_litter(i)                    = value_patch
     this%sen_ploss_litter(i)                    = value_patch
     this%livestemp_to_litter(i)                 = value_patch

      if ( crop_prog )then
         this%grainp_to_food(i)                   = value_patch
         this%grainp_xfer_to_grainp(i)            = value_patch
         this%ppool_to_grainp(i)                  = value_patch
         this%ppool_to_grainp_storage(i)          = value_patch
         this%grainp_storage_to_xfer(i)           = value_patch
         this%frootp_to_retransp(i)               = value_patch
         this%crop_seedp_to_leaf(i)               = value_patch
      end if
end do 
!$acc end data
end subroutine veg_pf_setvalues_acc

end module VegetationSummaryRoutinesMod
