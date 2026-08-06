module CNFireEmissionsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Gathers carbon emissions from fire sources to be sent to CAM-Chem via
  ! the coupler .... 
  ! Created by F. Vitt, and revised by F. Li
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use GridcellType, only : grc_pp                
  use ColumnType,   only : col_pp                
  use decompMod,    only : bounds_type
  use shr_fire_emis_mod,  only : shr_fire_emis_comps_n, shr_fire_emis_comp_t, shr_fire_emis_linkedlist
  use shr_fire_emis_mod,  only : shr_fire_emis_mechcomps_n, shr_fire_emis_mechcomps
  use spmdMod,            only : masterproc
  use elm_varctl,         only : iulog
  use elm_varcon,         only : spval
  use VegetationType,     only : veg_pp                
  use VegetationDataType, only : veg_cs
  use VegetationDataType, only : veg_cf 
  use ColumnDataType,     only : col_cf                
  use ColumnDataType,     only : col_cs                
  !
  implicit none
  save
  private 
  !
  ! !PUBLIC MEMBER FUNCTIONS:
   public :: CNFireEmisUpdate
  !
  ! !PRIVATE TYPES:
  type, private :: emis_t
     real(r8), pointer :: emis(:)
  end type emis_t
  !
  ! !PUBLIC TYPES:
  type, public :: fireemis_type
     real(r8),     pointer, public  :: fireflx_patch(:,:) ! carbon flux from fire sources (kg/m2/sec)
     real(r8),     pointer, public  :: ztop_patch(:)      ! height of the smoke plume (meters)
     type(emis_t), pointer, private :: comp(:)            ! fire emissions component (corresponds to emis factors table input file)
     type(emis_t), pointer, private :: mech(:)            ! cam-chem mechism species emissions
     type(emis_t),          private :: totfire            ! sum of all species emissions
     type(emis_t),          private :: fire_crop          ! sum of all species emissions
     type(emis_t),          private :: fire_ncrop         ! sum of all species emissions
     type(emis_t),          private :: fire_deforest      ! sum of all species emissions
     type(emis_t),          private :: fire_brl_forest        ! sum of all species emissions
     type(emis_t),          private :: fire_tmp_forest        ! sum of all species emissions
     type(emis_t),          private :: fire_trp_forest        ! sum of all species emissions
     type(emis_t),          private :: fire_tmp_shrub         ! sum of all species emissions
     type(emis_t),          private :: fire_tmp_grass         ! sum of all species emissions
     type(emis_t),          private :: fire_savanna         ! sum of all species emissions
     type(emis_t),          private :: cpool_other         ! sum of all species emissions
   contains
     procedure, public  :: Init
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
  end type fireemis_type
  !------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)

    use shr_fire_emis_mod,  only : shr_fire_emis_factors_file
    use FireEmisFactorsMod, only : fire_emis_factors_init, fire_emis_factors_get
    use elm_varpar,         only : numpft

    implicit none

    ! args
    class(fireemis_type) :: this
    type(bounds_type), intent(in) :: bounds

    ! local vars
    integer :: nmech, nemis
    real(r8) :: factors(numpft)
    real(r8) :: molec_wght
    type(shr_fire_emis_comp_t), pointer :: emis_cmp
    
    if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'fire_emis_nl settings:'
      write(iulog,*) '  shr_fire_emis_mechcomps_n  = ', shr_fire_emis_mechcomps_n 
      write(iulog,*) '  shr_fire_emis_factors_file  = ', shr_fire_emis_factors_file
!      write(iulog,*) '  fire_emis_coeff  = ', emis_cmp%coeff
!      write(iulog,*) '  shr_fire_emis_linkedlist  = ', shr_fire_emis_linkedlist
      write(iulog,*) ' '
    endif
    
    if ( shr_fire_emis_mechcomps_n < 1) return

    call fire_emis_factors_init( shr_fire_emis_factors_file )

    emis_cmp => shr_fire_emis_linkedlist
    do while(associated(emis_cmp))
       allocate(emis_cmp%emis_factors(numpft))
       call fire_emis_factors_get( trim(emis_cmp%name), factors, molec_wght )
       emis_cmp%emis_factors = factors*1.e-3_r8 ! convert g/kg dry fuel to kg/kg
       emis_cmp%molec_weight = molec_wght
       emis_cmp => emis_cmp%next_emiscomp
    enddo

    call this%InitAllocate(bounds) 
    call this%InitHistory(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! Allocate memory for module datatypes
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use elm_varcon      , only : spval

    ! !ARGUMENTS:
    class(fireemis_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp, i
    !---------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp


    allocate(this%totfire%emis(begp:endp)); this%totfire%emis(:) = nan
    allocate(this%fire_crop%emis(begp:endp)); this%fire_crop%emis(:) = nan
    allocate(this%fire_ncrop%emis(begp:endp)); this%fire_ncrop%emis(:) = nan
    allocate(this%fire_deforest%emis(begp:endp)); this%fire_deforest%emis(:) = nan

    allocate(this%fire_brl_forest%emis(begp:endp)); this%fire_brl_forest%emis(:) = nan
    allocate(this%fire_tmp_forest%emis(begp:endp)); this%fire_tmp_forest%emis(:) = nan
    allocate(this%fire_trp_forest%emis(begp:endp)); this%fire_trp_forest%emis(:) = nan
    allocate(this%fire_tmp_shrub%emis(begp:endp));  this%fire_tmp_shrub%emis(:) = nan
    allocate(this%fire_tmp_grass%emis(begp:endp));  this%fire_tmp_grass%emis(:) = nan
    allocate(this%fire_savanna%emis(begp:endp));    this%fire_savanna%emis(:) = nan
    allocate(this%cpool_other%emis(begp:endp));    this%cpool_other%emis(:) = nan

    if (shr_fire_emis_mechcomps_n>0) then
       allocate(this%fireflx_patch(begp:endp,1:shr_fire_emis_mechcomps_n)); this%fireflx_patch(:,:) = nan
       allocate(this%ztop_patch(begp:endp)); this%ztop_patch(:) = nan
 
       allocate(this%mech(shr_fire_emis_mechcomps_n))
       do i = 1, shr_fire_emis_mechcomps_n
          allocate(this%mech(i)%emis(begp:endp)); this%mech(i)%emis(:) = nan
       enddo
    endif

    if (shr_fire_emis_comps_n>0) then
       allocate(this%comp(shr_fire_emis_comps_n))
       do i = 1, shr_fire_emis_comps_n
          allocate(this%comp(i)%emis(begp:endp)); this%comp(i)%emis(:) = nan
       enddo
    endif

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    use elm_varcon  , only : spval
    use histFileMod , only : hist_addfld1d

    ! !ARGUMENTS:
    class(fireemis_type) :: this
    type(bounds_type), intent(in) :: bounds  

    ! !LOCAL VARIABLES
    integer :: begp, endp
    integer :: imech, icomp
    type(shr_fire_emis_comp_t), pointer :: emis_cmp

    begp = bounds%begp; endp = bounds%endp
 
   if (shr_fire_emis_mechcomps_n>0) then

       emis_cmp => shr_fire_emis_linkedlist

       ! loop over fire components
       emis_cmp_loop: do while(associated(emis_cmp))

          icomp = emis_cmp%index

          this%comp(icomp)%emis(begp:endp) = spval
          call hist_addfld1d (fname='FireComp_'//trim(emis_cmp%name), units='kg/m2/sec', &
               avgflag='A', long_name='fire emissions flux of '//trim(emis_cmp%name), &
               ptr_patch=this%comp(icomp)%emis)

          emis_cmp => emis_cmp%next_emiscomp

       enddo emis_cmp_loop


       ! loop over atm chem mechanism species
       do imech = 1,shr_fire_emis_mechcomps_n


          this%mech(imech)%emis(begp:endp) = spval
	  call hist_addfld1d (fname='FireMech_'//trim(shr_fire_emis_mechcomps(imech)%name), units='kg/m2/sec', &
               avgflag='A', long_name='fire emissions flux of '//trim(shr_fire_emis_mechcomps(imech)%name), &
               ptr_patch=this%mech(imech)%emis)

       enddo

       this%totfire%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_TOT', units='gC/m2/sec', &
            avgflag='A', long_name='Total fire emissions flux ', &
            ptr_patch=this%totfire%emis)

!LXu@01/2026+++
       this%fire_crop%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_Crop', units='gC/m2/sec', &
            avgflag='A', long_name='Fire emissions flux from agriculture burnings ', &
            ptr_patch=this%fire_crop%emis)

       this%fire_ncrop%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_NCrop', units='gC/m2/sec', &
            avgflag='A', long_name='Fire emissions flux from vegetation fires ', &
            ptr_patch=this%fire_ncrop%emis)

       this%fire_deforest%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_Deforest', units='gC/m2/sec', &
            avgflag='A', long_name='Fire emissions flux from deforestation fires ', &
            ptr_patch=this%fire_deforest%emis)

       this%fire_brl_forest%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_Brl_Forest', units='gC/m2/sec', &
            avgflag='A', long_name='Fire emissions flux from boreal forest fires ', &
            ptr_patch=this%fire_brl_forest%emis)

       this%fire_tmp_forest%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_Tmp_Forest', units='gC/m2/sec', &
            avgflag='A', long_name='Fire emissions flux from temperate forest fires ', &
            ptr_patch=this%fire_tmp_forest%emis)

       this%fire_trp_forest%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_Trp_Forest', units='gC/m2/sec', &
            avgflag='A', long_name='Fire emissions flux from tropical forest fires ', &
            ptr_patch=this%fire_trp_forest%emis)

       this%fire_tmp_shrub%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_Tmp_Shrub', units='gC/m2/sec', &
            avgflag='A', long_name='Fire emissions flux from temperate shrub fires ', &
            ptr_patch=this%fire_tmp_shrub%emis)

       this%fire_tmp_grass%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_Tmp_Grass', units='gC/m2/sec', &
            avgflag='A', long_name='Fire emissions flux from temperate grass fires ', &
            ptr_patch=this%fire_tmp_grass%emis)

       this%fire_savanna%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_Savanna', units='gC/m2/sec', &
            avgflag='A', long_name='Fire emissions flux from savanna fires ', &
            ptr_patch=this%fire_savanna%emis)

       this%cpool_other%emis(begp:endp) = spval
       call hist_addfld1d (fname='Cpool_other', units='gC/m2', &
            avgflag='A', long_name='Carbon in other pools except leaf, livestem and deadstem ', &
            ptr_patch=this%cpool_other%emis)
!LXu@01/2026---

       this%ztop_patch(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_ZTOP', units='m', &
            avgflag='A', long_name='Top of vertical fire emissions distribution ', &
            ptr_patch=this%ztop_patch)
    endif

 
  end subroutine InitHistory

  !-----------------------------------------------------------------------
!LXu@02/20+++++
!  subroutine CNFireEmisUpdate(bounds, num_soilp, filter_soilp, cnveg_cf_vars, cnveg_cs_inst, fireemis_inst )
!  subroutine CNFireEmisUpdate(bounds, num_soilp, filter_soilp, cnveg_cf_vars, cnveg_cs_vars, fireemis_vars )
  subroutine CNFireEmisUpdate(bounds, num_soilp, filter_soilp, cnstate_vars, fireemis_vars )
!  subroutine CNFireEmisUpdate(bounds, num_soilp, filter_soilp, col_cf, veg_cf, veg_cs, fireemis_vars )
!LXu@02/20-----

!    use CNVegcarbonfluxType,  only : cnveg_carbonflux_type
!    use CNVegCarbonStateType, only : cnveg_carbonstate_type 
!    use VegetationDataType     , only : vegetation_carbon_flux
!    use VegetationDataType     , only : vegetation_carbon_state      

!    use CNCarbonFluxType       , only : carbonflux_type
!    use CNCarbonStateType      , only : carbonstate_type
    use CNStateType,           only : cnstate_type
    use elm_varpar,            only : ndecomp_pools, nlevdecomp
    use elm_varcon,            only : dzsoi_decomp
    use dynSubgridControlMod , only : run_has_transient_landcover
    use pftvarcon,             only : nbrdlf_evr_shrub, nbrdlf_dcd_tmp_shrub, nbrdlf_dcd_brl_shrub
    use pftvarcon,             only : nc3_arctic_grass, nc3_nonarctic_grass, nc4_grass, nc3crop
    use pftvarcon,             only : ndllf_evr_tmp_tree, ndllf_evr_brl_tree
    use pftvarcon,             only : ndllf_dcd_brl_tree, nbrdlf_evr_tmp_tree
    use pftvarcon,             only : nbrdlf_dcd_tmp_tree, nbrdlf_dcd_brl_tree
    use pftvarcon,             only : nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree

    !ARGUMENTS:
    type(bounds_type),           intent(in)     :: bounds                  
    integer,                     intent(in)     :: num_soilp       ! number of soil pfts in filter
    integer,                     intent(in)     :: filter_soilp(:) ! filter for soil pfts
!    integer,                     intent(in)    :: filter_soilp(num_solip) ! filter for soil pfts
!    type(carbonflux_type),intent(in)     :: cnveg_cf_vars
!    type(carbonstate_type),intent(in)    :: cnveg_cs_vars 
!    type(carbonflux_type),       intent(in)    :: cnveg_cf_vars
    type(cnstate_type),          intent(in)    :: cnstate_vars 
    type(fireemis_type),         intent(inout) :: fireemis_vars

    !LOCAL VARIABLES:
    real(r8) :: fire_flux
    real(r8) :: fire_flux_lf 
    real(r8) :: fire_flux_lf1 
    real(r8) :: fire_flux_crop, fire_flux_ncrop, fire_flux_deforest 
    real(r8) :: fire_flux_brl_forest,  fire_flux_tmp_forest,fire_flux_trp_forest
    real(r8) :: fire_flux_tmp_shrub, fire_flux_tmp_grass, fire_flux_savanna, flux_cpool_other
    type(shr_fire_emis_comp_t), pointer :: emis_cmp
    real(r8) :: emis_flux(shr_fire_emis_comps_n)
    integer  :: fp,p,g,c                ! indices
    real(r8) :: epsilon                 ! emission factor [ug m-2 h-1]
    integer  :: i, ii, icomp, imech, n_emis_comps, l, j
    logical           :: transient_landcover  ! whether this run has any prescribed transient landcover
    integer  :: kyr, kmo, kda, mcsec, nstep

!LXu@05/20+++++
    real(r8)              :: dmr                    ! the ratio of DM/Carbon_flux from fire emissions [kg/kg]
    real(r8), parameter   :: eqas_latS = -10.0_r8   ! Latitude for Equatorial Asia peat fires
    real(r8), parameter   :: eqas_latN =   8.0_r8   ! Latitude for Equatorial Asia peat fires
    real(r8), parameter   :: eqas_lonL =  95.0_r8   ! Longitude for Equatorial Asia peat fires
    real(r8), parameter   :: eqas_lonR = 160.0_r8   ! Longitude for Equatorial Asia peat fires
!LXu@05/20 derived from GFED4 emission factors
!    real(r8), parameter   :: ef_peat(6) = (/0.10_r8, 14.2_r8, 0.1075_r8, 4.3_r8, 3.88_r8, 335.4_r8/)   ! BC,OC,SO4,SO2,SOAG,CO g/kg(DM)
!    real(r8), parameter   :: ef_peat(8) = (/0.10_r8, 14.2_r8, 0.1075_r8, 4.3_r8, 3.88_r8, 335.4_r8, 0.0104_r8, 0.0337_r8/)   ! BC,OC,SO4,SO2,SOAG,CO,PO4_a4,PO4_a3 g/kg(DM)
!    real(r8), parameter   :: dmr_peat   = 0.57_r8   ! C/DM ratio
!    real(r8), parameter   :: dm_ratio(16) = &     
!        				 (/ 0.50_r8, 0.49_r8, 0.49_r8, 0.50_r8, 0.50_r8, 0.50_r8, 0.50_r8, 0.49_r8,  &
!        				    0.50_r8, 0.50_r8, 0.49_r8, 0.49_r8, 0.49_r8, 0.49_r8, 0.44_r8, 0.44_r8 /)
!LXu@03/2025 derived from GFED5 emission factors , unit is g/kg(DM)
! BC,OC,SO4,SO2,SOAG (using OC)
! CO, C2H4,C2H6, C3H8, CH2O,  
! CH3CHO, CH3COCH3, ISOP, C10H16, NO
    real(r8), parameter   :: ef_peat(15) = (/  0.02_r8, 13.17_r8,  2.06_r8,  2.06_r8, 13.17_r8, &
                                             225.00_r8,  1.53_r8,  2.52_r8,  4.06_r8,  1.35_r8, &
					       1.81_r8, 0.926_r8, 0.646_r8, 0.241_r8,  0.93_r8/)   
    real(r8), parameter   :: dmr_peat   = 0.53_r8   ! C/DM ratio
    real(r8), parameter   :: dm_ratio(16) = &     
        				 (/ 0.48_r8, 0.49_r8, 0.49_r8, 0.50_r8, 0.48_r8, 0.50_r8, 0.48_r8, 0.49_r8,  &
        				    0.48_r8, 0.48_r8, 0.49_r8, 0.48_r8, 0.48_r8, 0.48_r8, 0.42_r8, 0.42_r8 /)
!LXu@05/20-----

    if ( shr_fire_emis_mechcomps_n < 1) return

    associate( & 
         cropf_col                   => cnstate_vars%cropf_col            , & ! Input:  [real(r8) (:)     ]  cropland fraction in veg column
         baf_peatf                   => cnstate_vars%baf_peatf_col        , & ! Output: [real(r8) (:)     ]  burned area fraction for peatland (/sec)
         trotr1_col                  => cnstate_vars%trotr1_col           , & ! Output: [real(r8) (:)     ]  pft weight of BET on the gridcell (0-1)
         trotr2_col                  => cnstate_vars%trotr2_col           , & ! Output: [real(r8) (:)     ]  pft weight of BDT on the gridcell (0-1)
         dtrotr_col                  => cnstate_vars%dtrotr_col           , & ! Input:  [real(r8) (:)     ]  ann. decreased frac. coverage of BET+BDT (0-1) on GC
         totvegc_col                 => col_cs%totvegc                    , & ! Output: [real(r8) (:)     ]  totvegc at column level
         totvegc                     => veg_cs%totvegc                    , & ! Input:  [real(r8) (:)     ]  (gC/m2) total vegetation carbon, excluding cpool
         leafc_storage                       =>    veg_cs%leafc_storage       , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C storage
         leafc_xfer                          =>    veg_cs%leafc_xfer          , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C transfer
         livestemc_storage                   =>    veg_cs%livestemc_storage   , & ! Input:  [real(r8) (:)     ]  (gC/m2) live stem C storage
         livestemc_xfer                      =>    veg_cs%livestemc_xfer      , & ! Input:  [real(r8) (:)     ]  (gC/m2) live stem C transfer
         deadstemc_storage                   =>    veg_cs%deadstemc_storage   , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead stem C storage
         deadstemc_xfer                      =>    veg_cs%deadstemc_xfer      , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead stem C transfer
         frootc                              =>    veg_cs%frootc              , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C
         frootc_storage                      =>    veg_cs%frootc_storage      , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C storage
         frootc_xfer                         =>    veg_cs%frootc_xfer         , & ! Input:  [real(r8) (:)     ]  (gC/m2) fine root C transfer
         livecrootc                          =>    veg_cs%livecrootc          , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C
         livecrootc_storage                  =>    veg_cs%livecrootc_storage  , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C storage
         livecrootc_xfer                     =>    veg_cs%livecrootc_xfer     , & ! Input:  [real(r8) (:)     ]  (gC/m2) live coarse root C transfer
         deadcrootc                          =>    veg_cs%deadcrootc          , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C
         deadcrootc_storage                  =>    veg_cs%deadcrootc_storage  , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C storage
         deadcrootc_xfer                     =>    veg_cs%deadcrootc_xfer     , & ! Input:  [real(r8) (:)     ]  (gC/m2) dead coarse root C transfer
         gresp_storage                       =>    veg_cs%gresp_storage       , & ! Input:  [real(r8) (:)     ]  (gC/m2) growth respiration storage
         gresp_xfer                          =>    veg_cs%gresp_xfer          , & ! Input:  [real(r8) (:)     ]  (gC/m2) growth respiration transfer
         m_decomp_cpools_to_fire_vr  => col_cf%m_decomp_cpools_to_fire_vr , & ! Output: [real(r8) (:,:,:) ]  (gC/m3/s) VR decomp. C fire loss
         somc_fire                   => col_cf%somc_fire                  , & ! Input:  [real(r8) (:)     ]  (gC/m2/s) fire C emissions due to peat burning
         fire_emis                   => fireemis_vars%fireflx_patch       , &
         totfire                     => fireemis_vars%totfire             , &
         fire_crop                   => fireemis_vars%fire_crop           , &
         fire_ncrop                  => fireemis_vars%fire_ncrop          , &
         fire_deforest               => fireemis_vars%fire_deforest       , &
         fire_brl_forest             => fireemis_vars%fire_brl_forest    , &
         fire_tmp_forest             => fireemis_vars%fire_tmp_forest    , &
         fire_trp_forest             => fireemis_vars%fire_trp_forest    , &
         fire_tmp_shrub              => fireemis_vars%fire_tmp_shrub          , &
         fire_tmp_grass              => fireemis_vars%fire_tmp_grass          , &
         fire_savanna                => fireemis_vars%fire_savanna          , &
         cpool_other                 => fireemis_vars%cpool_other          , &
         mech                        => fireemis_vars%mech                , &
         comp                        => fireemis_vars%comp                , &
         ztop                        => fireemis_vars%ztop_patch            &
         )

     transient_landcover = run_has_transient_landcover()

      ! initialize to zero ...
      fire_emis(bounds%begp:bounds%endp,:) = 0._r8
      totfire%emis(bounds%begp:bounds%endp) =  0._r8
      ztop(bounds%begp:bounds%endp) =  0._r8

      fire_crop%emis(bounds%begp:bounds%endp) =  0._r8
      fire_ncrop%emis(bounds%begp:bounds%endp) =  0._r8
      fire_deforest%emis(bounds%begp:bounds%endp) =  0._r8

      fire_brl_forest%emis(bounds%begp:bounds%endp) =  0._r8
      fire_tmp_forest%emis(bounds%begp:bounds%endp) =  0._r8
      fire_trp_forest%emis(bounds%begp:bounds%endp) =  0._r8
      fire_tmp_shrub%emis(bounds%begp:bounds%endp) =  0._r8
      fire_tmp_grass%emis(bounds%begp:bounds%endp) =  0._r8
      fire_savanna%emis(bounds%begp:bounds%endp) =  0._r8
      cpool_other%emis(bounds%begp:bounds%endp) =  0._r8

      do i = 1, shr_fire_emis_mechcomps_n
         mech(i)%emis(bounds%begp:bounds%endp) =  0._r8
      enddo

      do i = 1, shr_fire_emis_comps_n
         comp(i)%emis(bounds%begp:bounds%endp) =  0._r8
      enddo

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         g = veg_pp%gridcell(p)
         c = veg_pp%column(p)

         ! initialize EF
         epsilon=0._r8
         emis_flux(:) = 0._r8
	 dmr = 0.5_r8


         ! calculate fire emissions for non-bare ground PFTs
         if (veg_pp%itype(p) > 0)then

          ! calculate fire emissions for non-bare ground PFTs
            if(totvegc_col(c) > 0._r8)then
               fire_flux_lf1=0._r8 
               do l = 1, ndecomp_pools
                  do j = 1, nlevdecomp
                     fire_flux_lf1 = fire_flux_lf1 + &
                          m_decomp_cpools_to_fire_vr(c,j,l)*dzsoi_decomp(j) !fire_flux_lf1 is column-level variable
                  end do
               end do
!               fire_flux_lf = fire_flux_lf1*veg_cs%totvegc_patch(p)/col_cs%totvegc_col(c)
               fire_flux_lf = fire_flux_lf1*totvegc(p)/totvegc_col(c)
            else
               fire_flux_lf=0._r8 
            end if
             

!           fire_flux =  veg_cf%m_leafc_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc
           if ( veg_pp%itype(p) < nc3crop .and. cropf_col(c) < 1.0_r8) then
                fire_flux_ncrop =  fire_flux_lf &
                 + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                 + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                 + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                 + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                 + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                 + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                 + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                 + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                 + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                 + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                 + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                 + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                 + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                 + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                 + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                 + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                 + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                 + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                 + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                 + veg_cf%m_gresp_xfer_to_fire                (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 

               flux_cpool_other =  leafc_storage(p) & ! (gC/m2/s) fire C emissions from leafc_storage
                 + leafc_xfer(p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                 + livestemc_storage(p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                 + livestemc_xfer(p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                 + deadstemc_storage(p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                 + deadstemc_xfer(p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                 + frootc(p) & ! (gC/m2/s) fire C emissions from frootc
                 + frootc_storage(p) & ! (gC/m2/s) fire C emissions from frootc_storage
                 + frootc_xfer(p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                 + livecrootc(p) & ! (gC/m2/s) fire C emissions from livecrootc
                 + livecrootc_storage(p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                 + livecrootc_xfer(p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                 + deadcrootc(p) & ! (gC/m2/s) fire C emissions from deadcrootc
                 + deadcrootc_storage(p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                 + deadcrootc_xfer(p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                 + gresp_storage(p) & ! (gC/m2/s) fire C emissions from gresp_storage
                 + gresp_xfer(p)   ! (gC/m2/s) fire C emissions from gresp_xfer 
	                              
           ! for diagnostics
!           if( .not. (kmo == 1 .and. kda == 1 .and. mcsec == 0) )then
!              if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 .and. dtrotr_col(c) > 0._r8 .and. &
!                   lfc(c) > 0._r8 .and. fbac1(c) == 0._r8) then
!                 lfc2(c) = max(0._r8, min(lfc(c), (farea_burned(c)-baf_crop(c) - &
!                      baf_peatf(c))/2.0*dt))/(dtrotr_col(c)*dayspyr*secspday/dt)/dt
!                 lfc(c)  = lfc(c) - max(0._r8, min(lfc(c), (farea_burned(c)-baf_crop(c) - &
!                      baf_peatf(c))*dt/2.0_r8))
!              end if
!           end if
               if (transient_landcover) then    !true when landuse change data is used
        	  if( trotr1_col(c)+trotr2_col(c) > 0.6_r8 .and. dtrotr_col(c) > 0._r8)then
        	   fire_flux_deforest = fire_flux_lf &
                    + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                    + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                    + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                    + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                    + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                    + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                    + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                    + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                    + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                    + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                    + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                    + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                    + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                    + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                    + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                    + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                    + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                    + veg_cf%m_gresp_xfer_to_fire                (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 
        	  else
		   fire_flux_deforest = 0._r8
                  end if
	       end if
! boreal forest
               if (veg_pp%itype(p) == ndllf_evr_brl_tree  &
	      .or. veg_pp%itype(p) == ndllf_dcd_brl_tree  &
	      .or. veg_pp%itype(p) == nbrdlf_dcd_brl_tree &    
	      .or. veg_pp%itype(p) == nbrdlf_dcd_brl_shrub  &
	      .or. veg_pp%itype(p) == nc3_arctic_grass ) then
                  fire_flux_brl_forest = fire_flux_lf &
                    + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                    + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                    + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                    + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                    + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                    + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                    + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                    + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                    + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                    + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                    + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                    + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                    + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                    + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                    + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                    + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                    + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                    + veg_cf%m_gresp_xfer_to_fire                (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 
               else
		   fire_flux_brl_forest = 0._r8
	       end if
! temperate forest
               if (veg_pp%itype(p) == ndllf_evr_tmp_tree   &
	      .or. veg_pp%itype(p) == nbrdlf_evr_tmp_tree  &
	      .or. veg_pp%itype(p) == nbrdlf_dcd_tmp_tree ) then    !forest
        	   fire_flux_tmp_forest = fire_flux_lf &
                    + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                    + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                    + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                    + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                    + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                    + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                    + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                    + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                    + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                    + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                    + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                    + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                    + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                    + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                    + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                    + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                    + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                    + veg_cf%m_gresp_xfer_to_fire                (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 
               else
		   fire_flux_tmp_forest = 0._r8
	       end if
! tropical forest
               if (veg_pp%itype(p) == nbrdlf_evr_trp_tree  ) then    
        	   fire_flux_trp_forest = fire_flux_lf &
                    + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                    + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                    + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                    + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                    + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                    + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                    + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                    + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                    + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                    + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                    + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                    + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                    + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                    + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                    + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                    + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                    + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                    + veg_cf%m_gresp_xfer_to_fire                (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 
               else
		   fire_flux_trp_forest = 0._r8
	       end if
!temperate shrub
               if (veg_pp%itype(p) == nbrdlf_evr_shrub .or. veg_pp%itype(p) == nbrdlf_dcd_tmp_shrub) then    
        	   fire_flux_tmp_shrub = fire_flux_lf &
                    + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                    + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                    + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                    + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                    + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                    + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                    + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                    + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                    + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                    + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                    + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                    + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                    + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                    + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                    + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                    + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                    + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                    + veg_cf%m_gresp_xfer_to_fire                (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 
               else
		   fire_flux_tmp_shrub = 0._r8
	       end if
!savanna
               if (veg_pp%itype(p) == nc4_grass .or. veg_pp%itype(p) == nbrdlf_dcd_trp_tree) then    
        	   fire_flux_savanna = fire_flux_lf &
                    + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                    + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                    + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                    + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                    + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                    + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                    + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                    + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                    + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                    + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                    + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                    + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                    + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                    + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                    + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                    + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                    + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                    + veg_cf%m_gresp_xfer_to_fire                (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 
               else
		   fire_flux_savanna = 0._r8
	       end if
!temperate grass
               if (veg_pp%itype(p) == nc3_arctic_grass ) then    !forest
                  if(grc_pp%latdeg(g) < -23._r8 .or. grc_pp%latdeg(g) > 23._r8 )then
        	   fire_flux_tmp_grass = fire_flux_lf &
                    + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                    + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                    + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                    + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                    + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                    + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                    + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                    + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                    + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                    + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                    + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                    + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                    + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                    + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                    + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                    + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                    + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                    + veg_cf%m_gresp_xfer_to_fire                (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 
                  else
        	   fire_flux_savanna = fire_flux_savanna + fire_flux_lf &
                    + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                    + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                    + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                    + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                    + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                    + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                    + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                    + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                    + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                    + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                    + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                    + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                    + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                    + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                    + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                    + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                    + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                    + veg_cf%m_gresp_xfer_to_fire                (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 
		  end if
	       else
		   fire_flux_tmp_grass = 0._r8
	       end if
           else
                fire_flux_crop = fire_flux_lf  &
                    + veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
                    + veg_cf%m_leafc_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from leafc_storage
                    + veg_cf%m_leafc_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from leafc_xfer
                    + veg_cf%m_livestemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from livestemc
                    + veg_cf%m_livestemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from livestemc_storage
                    + veg_cf%m_livestemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from livestemc_xfer
                    + veg_cf%m_deadstemc_to_fire                 (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_deadstemc_storage_to_fire         (p) & ! (gC/m2/s) fire C emissions from deadstemc_storage
                    + veg_cf%m_deadstemc_xfer_to_fire            (p) & ! (gC/m2/s) fire C emissions from deadstemc_xfer
                    + veg_cf%m_frootc_to_fire                    (p) & ! (gC/m2/s) fire C emissions from frootc
                    + veg_cf%m_frootc_storage_to_fire            (p) & ! (gC/m2/s) fire C emissions from frootc_storage
                    + veg_cf%m_frootc_xfer_to_fire               (p) & ! (gC/m2/s) fire C emissions from frootc_xfer
                    + veg_cf%m_livecrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from livecrootc
                    + veg_cf%m_livecrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from livecrootc_storage 
                    + veg_cf%m_livecrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from livecrootc_xfer
                    + veg_cf%m_deadcrootc_to_fire                (p) & ! (gC/m2/s) fire C emissions from deadcrootc
                    + veg_cf%m_deadcrootc_storage_to_fire        (p) & ! (gC/m2/s) fire C emissions from deadcrootc_storage
                    + veg_cf%m_deadcrootc_xfer_to_fire           (p) & ! (gC/m2/s) fire C emissions from deadcrootc_xfer
                    + veg_cf%m_gresp_storage_to_fire             (p) & ! (gC/m2/s) fire C emissions from gresp_storage
                    + veg_cf%m_gresp_xfer_to_fire                (p) 
	   end if


!LXu@01/26, diagnose fire c flux 
           fire_deforest%emis(p) = fire_flux_deforest  !  gC/m2/sec
           fire_ncrop%emis(p) = fire_flux_ncrop        !  gC/m2/sec
           fire_crop%emis(p) = fire_flux_crop          !  gC/m2/sec

           fire_brl_forest%emis(p) = fire_flux_brl_forest      !  gC/m2/sec
           fire_tmp_forest%emis(p) = fire_flux_tmp_forest      !  gC/m2/sec
           fire_trp_forest%emis(p) = fire_flux_trp_forest      !  gC/m2/sec
           fire_tmp_shrub%emis(p) = fire_flux_tmp_shrub        !  gC/m2/sec
           fire_tmp_grass%emis(p) = fire_flux_tmp_grass        !  gC/m2/sec
           fire_savanna%emis(p) = fire_flux_savanna       !  gC/m2/sec
           cpool_other%emis(p) = flux_cpool_other        !  gC/m2/sec

!              if(baf_peatf(c) > 0._r8)then
!                 fire_flux_peat =  somc_fire(c) *totvegc(p)/totvegc_col(c)  ! (gC/m2/s) fire C emissions from peat
!              else
!                 fire_flux_peat =  0._r8
!	      end if
!            fire_peat%emis(p) = 0._r8!  gC/m2/sec
!LXu@01/26, add the peat fire flux---
	    fire_flux = fire_flux_ncrop + fire_flux_crop  

            ! for diagnostics
            totfire%emis(p) = fire_flux !  gC/m2/sec


            ! loop over fire components
            emis_cmp => shr_fire_emis_linkedlist
            emis_cmp_loop: do while(associated(emis_cmp))

               icomp = emis_cmp%index
!               epsilon = emis_cmp%emis_factors(veg_pp%itype(p))
!	       dmr = dm_ratio(veg_pp%itype(p))

!               comp(icomp)%emis(p) = epsilon * fire_flux* 1.e-3_r8/0.5_r8  ! (to convert gC/m2/sec to kg species/m2/sec)
!               comp(icomp)%emis(p) = epsilon * (fire_flux/dmr) * 1.e-3_r8 ! EF * DM (to convert gC/m2/sec to kg species/m2/sec) 
!               emis_flux(icomp) = emis_cmp%coeff*comp(icomp)%emis(p)

	       ! updated the fire emission in the Indonesia region for tropical peat-land fires
	       ! use tropical peat fire emissio nfactors
	       !  ef_peat(15) = (/  0.02_r8, 13.17_r8,  2.06_r8,  2.06_r8, 13.17_r8, &
	       !                  225.00_r8,  1.53_r8,  2.52_r8,  4.06_r8,  1.35_r8, &
	       !	   1.81_r8, 0.926_r8, 0.646_r8, 0.241_r8,  0.93_r8/)   
	       !  dmr_peat   = 0.53_r8   ! C/DM ratio
!               if (      grc_pp%latdeg(g) > eqas_latS .and. grc_pp%latdeg(g) < eqas_latN  &
!	           .and. grc_pp%londeg(g) > eqas_lonL .and. grc_pp%londeg(g) < eqas_lonR   ) then
!	       
!        	  epsilon = ef_peat(icomp)*1.e-3_r8 ! convert g/kg dry fuel to kg/kg(DM)
!		  dmr = dmr_peat
!        	  comp(icomp)%emis(p) = epsilon * (fire_flux/dmr) * 1.e-3_r8 ! EF * DM (to convert gC/m2/sec to kg species/m2/sec) 
!	       
!	       else

        	  epsilon = emis_cmp%emis_factors(veg_pp%itype(p))
		  dmr = dm_ratio(veg_pp%itype(p))
        	  comp(icomp)%emis(p) = epsilon * (fire_flux/dmr) * 1.e-3_r8 ! EF * DM (to convert gC/m2/sec to kg species/m2/sec) 
	       
!	       end if

               emis_flux(icomp) = emis_cmp%coeff*comp(icomp)%emis(p)
!       if (masterproc) then
!           write(iulog,*) '  fire_emis_coeff  = ', emis_cmp%coeff
!       end if
! No fire emission from ELM (harded coded)
!               emis_flux(icomp) = emis_cmp%coeff*comp(icomp)%emis(p)*0.0_r8
 
               emis_cmp => emis_cmp%next_emiscomp

            enddo emis_cmp_loop

            ! sum up the emissions compontent fluxes for the fluxes of chem mechanism compounds 
            do imech = 1,shr_fire_emis_mechcomps_n
               n_emis_comps = shr_fire_emis_mechcomps(imech)%n_emis_comps
               do icomp = 1,n_emis_comps ! loop over number of emission components that make up the nth mechanism compoud
                  ii = shr_fire_emis_mechcomps(imech)%emis_comps(icomp)%ptr%index
                  fire_emis(p,imech) = fire_emis(p,imech) + emis_flux(ii)
                  mech(imech)%emis(p) = fire_emis(p,imech)
               enddo
            enddo

            ztop(p) = vert_dist_top( veg_pp%itype(p) )

	 end if ! ivt(1:15 only)         if (veg_pp%itype(p) > 0)then

      enddo ! fp 
    end associate

  end subroutine CNFireEmisUpdate

! Private methods
!-----------------------------------------------------------------------
!ztop compiled from Val Martin et al ACP 2010, Tosca et al. JGR  2011 and Jian et al., ACP 2013
!st ztop updated based on Val Martin pers. communication Jan2015 
!-----------------------------------------------------------------------
!   not_vegetated    500 m                      
!PFT1: needleleaf_evergreen_temperate_tree     4000 m
!2: needleleaf_evergreen_boreal_tree    4000 m
!3: needleleaf_deciduous_boreal_tree    3000 m    
!4: broadleaf_evergreen_tropical_tree     2500 m  
!5: broadleaf_evergreen_temperate_tree   3000 m   
!6: broadleaf_deciduous_tropical_tree     2500 m  
!7: broadleaf_deciduous_temperate_tree  3000 m    
!8: broadleaf_deciduous_boreal_tree      3000 m   
!9: broadleaf_evergreen_shrub   2000 m            
!10: broadleaf_deciduous_temperate_shrub  2000 m  
!11: broadleaf_deciduous_boreal_shrub    2000 m    
!12: c3_arctic_grass   1000 m                      
!13: c3_non-arctic_grass  1000 m              
!14: c4_grass   1000 m                             
!15: c3_crop      1000 m
!(and all new crops: 1000m)

  function vert_dist_top( veg_type ) result(ztop)
    use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
    use pftvarcon       , only : noveg, ndllf_evr_tmp_tree, ndllf_evr_brl_tree
    use pftvarcon       , only : ndllf_dcd_brl_tree, nbrdlf_evr_tmp_tree
    use pftvarcon       , only : nbrdlf_dcd_tmp_tree, nbrdlf_dcd_brl_tree
    use pftvarcon       , only : nbrdlf_evr_trp_tree, nbrdlf_dcd_trp_tree
    use pftvarcon       , only : nbrdlf_evr_shrub, nbrdlf_dcd_brl_shrub
    use pftvarcon       , only : nc3_arctic_grass, nc3_nonarctic_grass
    use pftvarcon       , only : nc3crop, nc3irrig
    use pftvarcon       , only : npcropmin, npcropmax
    implicit none
    integer, intent(in) :: veg_type

    real(r8) :: ztop

    ! Bare soil, won't be used
    if (      veg_type == noveg ) then
       ztop = nan
    ! temperate and boreal evergreen needleleaf trees
    else if ( veg_type == ndllf_evr_tmp_tree  .or.  veg_type == ndllf_evr_brl_tree   ) then
       ztop = 4.e3_r8 ! m
    ! temperate and boreal trees
    else if ( veg_type == ndllf_dcd_brl_tree  .or.  veg_type == nbrdlf_evr_tmp_tree .or. &
              veg_type == nbrdlf_dcd_tmp_tree .or.  veg_type == nbrdlf_dcd_brl_tree  ) then
       ztop = 3.e3_r8 ! m
    ! tropical broadleaf trees (evergreen and decidious)
    else if ( veg_type == nbrdlf_evr_trp_tree .or.  veg_type == nbrdlf_dcd_trp_tree  ) then
       ztop = 2.5e3_r8 ! m
    ! shrubs
    else if ( veg_type >= nbrdlf_evr_shrub    .and. veg_type <= nbrdlf_dcd_brl_shrub ) then
       ztop = 2.e3_r8 ! m
    ! grasses
    else if ( veg_type >= nc3_arctic_grass    .and. veg_type <= nc3_nonarctic_grass  ) then
       ztop = 1.e3_r8 ! m
    ! generic unmanaged crops
    else if ( veg_type == nc3crop             .or.  veg_type <= nc3irrig             ) then
       ztop = 1.e3_r8 ! m
    ! Prognostic crops
    else if ( veg_type >= npcropmin           .and. veg_type <= npcropmax            ) then
       ztop = 1.e3_r8 ! m
    else
       call endrun('ERROR:: undefined veg_type' )
    end if

  end function vert_dist_top

end module CNFireEmissionsMod

