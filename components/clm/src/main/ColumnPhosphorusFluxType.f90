module ColumnPhosphorusFluxType

! ----------------------------------------------------------------------------
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod            , only : errMsg => shr_log_errMsg
  use clm_varpar             , only : ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar             , only : nlevdecomp_full, nlevdecomp, crop_prog
  use clm_varcon             , only : spval, ispval, dzsoi_decomp
  use decompMod              , only : bounds_type
  use clm_varctl             , only : use_nitrif_denitrif, use_vertsoilc
  use CNDecompCascadeConType , only : decomp_cascade_con
  use abortutils             , only : endrun
  use LandunitType           , only : lun                
  use ColumnType             , only : col                
  use PatchType              , only : pft
  !! bgc interface & pflotran:
  use clm_varctl             , only : use_bgc_interface, use_pflotran, pf_cmode, pf_hmode, use_vertsoilc
  ! 
  ! !PUBLIC TYPES:
  implicit none
  save
  private
	
! -----------------------------------------------------------------------------
	type, public :: soilcol_phosphorus_flux

     real(r8), pointer :: m_p_to_litr_met_fire_col                  (:,:)   ! col P from leaf, froot, xfer and storage P to litter labile P by fire (gP/m3/s) 
     real(r8), pointer :: m_p_to_litr_cel_fire_col                  (:,:)   ! col P from leaf, froot, xfer and storage P to litter cellulose P by fire (gP/m3/s) 
     real(r8), pointer :: m_p_to_litr_lig_fire_col                  (:,:)   ! col P from leaf, froot, xfer and storage P to litter lignin P by fire (gP/m3/s) 
     real(r8), pointer :: harvest_p_to_litr_met_p_col               (:,:)   ! col P fluxes associated with harvest to litter metabolic pool (gP/m3/s)
     real(r8), pointer :: harvest_p_to_litr_cel_p_col               (:,:)   ! col P fluxes associated with harvest to litter cellulose pool (gP/m3/s)
     real(r8), pointer :: harvest_p_to_litr_lig_p_col               (:,:)   ! col P fluxes associated with harvest to litter lignin pool (gP/m3/s)
     real(r8), pointer :: harvest_p_to_cwdp_col                     (:,:)   ! col P fluxes associated with harvest to CWD pool (gP/m3/s)

     ! crop harvest
     real(r8), pointer :: hrv_cropp_to_prod1p_col                   (:)     ! crop P harvest mortality to 1-yr product pool (gP/m2/s)

     ! fire P fluxes 
     real(r8), pointer :: m_decomp_ppools_to_fire_vr_col            (:,:,:) ! col vertically-resolved decomposing P fire loss (gP/m3/s)
     real(r8), pointer :: m_decomp_ppools_to_fire_col               (:,:)   ! col vertically-integrated (diagnostic) decomposing P fire loss (gP/m2/s)

     real(r8), pointer :: fire_ploss_col                            (:)     ! col total column-level fire P loss (gP/m2/s)

     real(r8), pointer :: fire_ploss_p2c_col                        (:)     ! col patch2col column-level fire P loss (gP/m2/s) (p2c)
     real(r8), pointer :: fire_mortality_p_to_cwdp_col              (:,:)   ! col P fluxes associated with fire mortality to CWD pool (gP/m3/s)

     ! summary (diagnostic) flux variables, not involved in mass balance

     real(r8), pointer :: wood_harvestp_col                         (:)     ! col total P losses to wood product pools (gP/m2/s) (p2c)

     ! deposition fluxes
     real(r8), pointer :: pdep_to_sminp_col                         (:)     ! col atmospheric P deposition to soil mineral P (gP/m2/s)
     real(r8), pointer :: fert_p_to_sminp_col                            (:)     ! col fertilizer P to soil mineral P (gP/m2/s)

     ! phenology: litterfall and crop fluxes
     real(r8), pointer :: phenology_p_to_litr_met_p_col             (:,:)   ! col P fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gP/m3/s)
     real(r8), pointer :: phenology_p_to_litr_cel_p_col             (:,:)   ! col P fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gP/m3/s)
     real(r8), pointer :: phenology_p_to_litr_lig_p_col             (:,:)   ! col P fluxes associated with phenology (litterfall and crop) to litter lignin pool (gP/m3/s)

     ! gap mortality fluxes
     real(r8), pointer :: gap_mortality_p_to_litr_met_p_col         (:,:)   ! col P fluxes associated with gap mortality to litter metabolic pool (gP/m3/s)
     real(r8), pointer :: gap_mortality_p_to_litr_cel_p_col         (:,:)   ! col P fluxes associated with gap mortality to litter cellulose pool (gP/m3/s)
     real(r8), pointer :: gap_mortality_p_to_litr_lig_p_col         (:,:)   ! col P fluxes associated with gap mortality to litter lignin pool (gP/m3/s)
     real(r8), pointer :: gap_mortality_p_to_cwdp_col               (:,:)   ! col P fluxes associated with gap mortality to CWD pool (gP/m3/s)

     ! decomposition fluxes
     real(r8), pointer :: decomp_cascade_ptransfer_vr_col           (:,:,:) ! col vert-res transfer of P from donor to receiver pool along decomp. cascade (gP/m3/s)
     real(r8), pointer :: decomp_cascade_ptransfer_col              (:,:)   ! col vert-int (diagnostic) transfer of P from donor to receiver pool along decomp. cascade (gP/m2/s)
     real(r8), pointer :: decomp_cascade_sminp_flux_vr_col          (:,:,:) ! col vert-res mineral P flux for transition along decomposition cascade (gP/m3/s)
     real(r8), pointer :: decomp_cascade_sminp_flux_col             (:,:)   ! col vert-int (diagnostic) mineral P flux for transition along decomposition cascade (gP/m2/s)
                                                                            ! Used to update concentrations concurrently with vertical transport
     ! vertically-resolved immobilization fluxes
     real(r8), pointer :: potential_immob_p_vr_col                    (:,:)   ! col vertically-resolved potential P immobilization (gP/m3/s) at each level
     real(r8), pointer :: potential_immob_p_col                       (:)     ! col vert-int (diagnostic) potential P immobilization (gP/m2/s)
     real(r8), pointer :: actual_immob_p_vr_col                       (:,:)   ! col vertically-resolved actual P immobilization (gP/m3/s) at each level
     real(r8), pointer :: actual_immob_p_col                          (:)     ! col vert-int (diagnostic) actual P immobilization (gP/m2/s)
     real(r8), pointer :: sminp_to_plant_vr_col                     (:,:)   ! col vertically-resolved plant uptake of soil mineral P (gP/m3/s)
     real(r8), pointer :: sminp_to_plant_col                        (:)     ! col vert-int (diagnostic) plant uptake of soil mineral P (gP/m2/s)
     real(r8), pointer :: supplement_to_sminp_vr_col                (:,:)   ! col vertically-resolved supplemental P supply (gP/m3/s)
     real(r8), pointer :: supplement_to_sminp_col                   (:)     ! col vert-int (diagnostic) supplemental P supply (gP/m2/s)
     real(r8), pointer :: gross_pmin_vr_col                         (:,:)   ! col vertically-resolved gross rate of P mineralization (gP/m3/s)
     real(r8), pointer :: gross_pmin_col                            (:)     ! col vert-int (diagnostic) gross rate of P mineralization (gP/m2/s)
     real(r8), pointer :: net_pmin_vr_col                           (:,:)   ! col vertically-resolved net rate of P mineralization (gP/m3/s)
     real(r8), pointer :: net_pmin_col                              (:)     ! col vert-int (diagnostic) net rate of P mineralization (gP/m2/s)

     real(r8), pointer :: biochem_pmin_ppools_vr_col                (:,:,:) ! col vertically-resolved biochemical P mineralization for each soi pool (gP/m3/s)
     real(r8), pointer :: biochem_pmin_vr_col                       (:,:)   ! col vertically-resolved total biochemical P mineralization (gP/m3/s)
     real(r8), pointer :: biochem_pmin_col                          (:)   ! col vert-int (diagnostic) total biochemical P mineralization (gP/m3/s)

     ! new variables for phosphorus code
     ! inorganic P transformation fluxes
     real(r8), pointer :: primp_to_labilep_vr_col                     (:,:)   ! col (gP/m3/s) flux of P from primary mineral to labile 
     real(r8), pointer :: primp_to_labilep_col                        (:)     ! col (gP/m3/s) flux of P from primary mineral to labile 
     real(r8), pointer :: labilep_to_secondp_vr_col                   (:,:)   ! col (gP/m3/s) flux of labile P to secondary mineral P 
     real(r8), pointer :: labilep_to_secondp_col                      (:)     ! col (gP/m3/s) flux of labile P to secondary mineral P 
     real(r8), pointer :: secondp_to_labilep_vr_col                   (:,:)   ! col (gP/m3/s) flux of the desorption of secondary mineral P to labile P
     real(r8), pointer :: secondp_to_labilep_col                      (:)     ! col (gP/m3/s) flux of the desorption of secondary mineral P to labile P
     real(r8), pointer :: secondp_to_occlp_vr_col                     (:,:)   ! col (gP/m3/s) flux of the occlusion of secondary P to occluded P
     real(r8), pointer :: secondp_to_occlp_col                        (:)     ! col (gP/m3/s) flux of the occlusion of secondary P to occluded P

     ! leaching fluxes
     real(r8), pointer :: sminp_leached_vr_col                      (:,:)   ! col vertically-resolved soil mineral P pool loss to leaching (gP/m3/s)
     real(r8), pointer :: sminp_leached_col                         (:)     ! col soil mineral P pool loss to leaching (gP/m2/s)

     ! dynamic landcover fluxes
     real(r8), pointer :: dwt_seedp_to_leaf_col                     (:)     ! col (gP/m2/s) seed source to PFT-level
     real(r8), pointer :: dwt_seedp_to_deadstem_col                 (:)     ! col (gP/m2/s) seed source to PFT-level
     real(r8), pointer :: dwt_conv_pflux_col                        (:)     ! col (gP/m2/s) conversion P flux (immediate loss to atm)
     real(r8), pointer :: dwt_prod10p_gain_col                      (:)     ! col (gP/m2/s) addition to 10-yr wood product pool
     real(r8), pointer :: dwt_prod100p_gain_col                     (:)     ! col (gP/m2/s) addition to 100-yr wood product pool
     real(r8), pointer :: dwt_frootp_to_litr_met_p_col              (:,:)   ! col (gP/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootp_to_litr_cel_p_col              (:,:)   ! col (gP/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_frootp_to_litr_lig_p_col              (:,:)   ! col (gP/m3/s) fine root to litter due to landcover change
     real(r8), pointer :: dwt_livecrootp_to_cwdp_col                (:,:)   ! col (gP/m3/s) live coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_deadcrootp_to_cwdp_col                (:,:)   ! col (gP/m3/s) dead coarse root to CWD due to landcover change
     real(r8), pointer :: dwt_ploss_col                             (:)     ! col (gP/m2/s) total phosphorus loss from product pools and conversion

     ! wood product pool loss fluxes
     real(r8), pointer :: prod1p_loss_col                           (:)     ! col (gP/m2/s) decomposition loss from 1-yr crop product pool
     real(r8), pointer :: prod10p_loss_col                          (:)     ! col (gP/m2/s) decomposition loss from 10-yr wood product pool
     real(r8), pointer :: prod100p_loss_col                         (:)     ! col (gP/m2/s) decomposition loss from 100-yr wood product pool
     real(r8), pointer :: product_ploss_col                         (:)     ! col (gP/m2/s) total wood product phosphorus loss

     ! summary (diagnostic) flux variables, not involved in mass balance
     real(r8), pointer :: pinputs_col                               (:)     ! col column-level P inputs (gP/m2/s)
     real(r8), pointer :: poutputs_col                              (:)     ! col column-level P outputs (gP/m2/s)
     real(r8), pointer :: som_p_leached_col                         (:)     ! col total SOM P loss from vertical transport (gP/m^2/s)
     real(r8), pointer :: decomp_ppools_leached_col                 (:,:)   ! col P loss from vertical transport from each decomposing P pool (gP/m^2/s)
     real(r8), pointer :: decomp_ppools_transport_tendency_col      (:,:,:) ! col P tendency due to vertical transport in decomposing P pools (gP/m^3/s)

     ! all n pools involved in decomposition
     real(r8), pointer :: decomp_ppools_sourcesink_col              (:,:,:) ! col (gP/m3) change in decomposing P pools

     ! clm_bgc_interface & pflotran
     !------------------------------------------------------------------------
     real(r8), pointer :: plant_pdemand_col                         (:)     ! col P flux required to support initial GPP (gN/m2/s)
     real(r8), pointer :: plant_pdemand_vr_col                      (:,:)   ! col vertically-resolved P flux required to support initial GPP (gP/m3/s)
     ! for PF-bgc mass-balance error checking
     real(r8), pointer :: externalp_to_decomp_ppools_col            (:,:,:) ! col net N fluxes associated with litter/som-adding/removal to decomp pools (gP/m3/s)
                                                                            ! (sum of all external P additions and removals, excluding decomposition/hr).
     real(r8), pointer :: externalp_to_decomp_delta_col             (:)     ! col summarized net N i/o changes associated with litter/som-adding/removal to decomp pools  btw time-step (gP/m2)
     real(r8), pointer :: sminp_net_transport_vr_col                (:,:)   ! col net sminp transport associated with runoff/leaching (gP/m3/s)
     real(r8), pointer :: sminp_net_transport_delta_col             (:)     ! col summarized net change of column-level sminp leaching bwtn time-step (for balance checking) (gP/m2)
     !------------------------------------------------------------------------
     real(r8), pointer :: adsorb_to_labilep_col                     (:)
     real(r8), pointer :: desorb_to_solutionp_col                   (:)

  contains
      procedure, public :: Init => init_col_pf
      procedure, public :: InitAllocate => initallocate_col_pf
      procedure, public :: InitHistory => inithistory_col_pf
      procedure, public :: InitCold => initcold_col_pf
      procedure, public :: Clean => clean_col_pf
  end type soilcol_phosphorus_flux


  !------------------------------------------------------------------------
  subroutine init_col_pf(this, bounds)

    class(phosphorusflux_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate (bounds)
    call this%InitHistory (bounds)
    call this%InitCold (bounds)

  end subroutine init_col_pf

  !------------------------------------------------------------------------


  !------------------------------------------------------------------------
  subroutine inithistory_col_pf(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : nlevsno, nlevgrnd, crop_prog 
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp
    !
    ! !ARGUMENTS:
    class(phosphorusflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer        :: k,l
    integer        :: begp, endp
    integer        :: begc, endc
    character(10)  :: active
    character(24)  :: fieldname
    character(100) :: longname
    character(8)   :: vr_suffix
    character(1)   :: aa 
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) 
    ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------   
    begc = bounds%begc; endc= bounds%endc

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    !-------------------------------
    ! P flux variables - native to column
    !-------------------------------

    this%pdep_to_sminp_col(begc:endc) = spval
    call hist_addfld1d (fname='PDEP_TO_SMINP', units='gP/m^2/s', &
         avgflag='A', long_name='atmospheric P deposition to soil mineral P', &
         ptr_col=this%pdep_to_sminp_col)


    do k = 1, ndecomp_pools
       if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
          this%m_decomp_ppools_to_fire_col(begc:endc,k) = spval
          data1dptr => this%m_decomp_ppools_to_fire_col(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'P_TO_FIRE'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' P fire loss'
          call hist_addfld1d (fname=fieldname, units='gP/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')

          if ( nlevdecomp_full > 1 ) then
             this%m_decomp_ppools_to_fire_vr_col(begc:endc,:,k) = spval
             data2dptr => this%m_decomp_ppools_to_fire_vr_col(:,:,k)
             fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'P_TO_FIRE'//trim(vr_suffix)
             longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' P fire loss'
             call hist_addfld_decomp (fname=fieldname, units='gP/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       endif

       !! biochemical P mineralization for each soil pool at each soil level for
       !! each column
       if ( k >= 5 )then

          if ( nlevdecomp_full > 1 ) then
             this%biochem_pmin_ppools_vr_col(begc:endc,:,k) = spval
             data2dptr => this%biochem_pmin_ppools_vr_col(:,:,k)
             write(aa,'(i1)') k 
             fieldname = 'BIOCHEM_PMIN_PPOOL'//aa//trim(vr_suffix)
             longname  = 'Biochemical mineralization of ppool'//aa
             call hist_addfld_decomp (fname=fieldname, units='gP/m^2/s',type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif
       endif
    end do

    do l = 1, ndecomp_cascade_transitions
       ! vertically integrated fluxes
       !-- mineralization/immobilization fluxes (none from CWD)
       if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
          this%decomp_cascade_sminp_flux_col(begc:endc,l) = spval
          data1dptr => this%decomp_cascade_sminp_flux_col(:,l)
          if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
             fieldname = 'SMINP_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'P_'//&
                  trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))
             longname =  'mineral P flux for decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                  'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
          else
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                  //'P_TO_SMINP'
             longname =  'mineral P flux for decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
          endif
          call hist_addfld1d (fname=fieldname, units='gP/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       end if

       !-- transfer fluxes (none from terminal pool, if present)
       if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
          this%decomp_cascade_ptransfer_col(begc:endc,l) = spval
          data1dptr => this%decomp_cascade_ptransfer_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'P_TO_'//&
               trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'P'
          longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
               ' P to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' N'
          call hist_addfld1d (fname=fieldname, units='gP/m^2',  &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
       end if

       ! vertically resolved fluxes
       if ( nlevdecomp_full > 1 ) then
          !-- mineralization/immobilization fluxes (none from CWD)
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
             this%decomp_cascade_sminp_flux_vr_col(begc:endc,:,l) = spval
             data2dptr => this%decomp_cascade_sminp_flux_vr_col(:,:,l)
             if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
                fieldname = 'SMINP_TO_'&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'P_'//&
                     trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_donor_pool(l)))//trim(vr_suffix)
                longname =  'mineral P flux for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//&
                     'to '//trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))
             else
                fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                     //'P_TO_SMINP'//trim(vr_suffix)
                longname =  'mineral P flux for decomp. of '&
                     //trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))
             endif
             call hist_addfld_decomp (fname=fieldname, units='gP/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif

          !-- transfer fluxes (none from terminal pool, if present)
          if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
             this%decomp_cascade_ptransfer_vr_col(begc:endc,:,l) = spval
             data2dptr => this%decomp_cascade_ptransfer_vr_col(:,:,l)
             fieldname = trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'P_TO_'//&
                  trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                  //'P'//trim(vr_suffix)
             longname =  'decomp. of '&
                  //trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                  ' P to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' P'
             call hist_addfld_decomp (fname=fieldname, units='gP/m^3',  type2d='levdcmp', &
                  avgflag='A', long_name=longname, &
                  ptr_col=data2dptr, default='inactive')
          endif

       endif
    end do

    this%som_p_leached_col(begc:endc) = spval
    call hist_addfld1d (fname='SOM_P_LEACHED', units='gP/m^2/s', &
         avgflag='A', long_name='total flux of P from SOM pools due to leaching', &
         ptr_col=this%som_p_leached_col, default='inactive')

    do k = 1, ndecomp_pools
       if ( .not. decomp_cascade_con%is_cwd(k) ) then
          this%decomp_ppools_leached_col(begc:endc,k) = spval
          data1dptr => this%decomp_ppools_leached_col(:,k)
          fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'P_TO_LEACHING'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' P leaching loss'
          call hist_addfld1d (fname=fieldname, units='gP/m^2/s', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr, default='inactive')

          this%decomp_ppools_transport_tendency_col(begc:endc,:,k) = spval
          data2dptr => this%decomp_ppools_transport_tendency_col(:,:,k)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(k))//'P_TNDNCY_VERT_TRANSPORT'
          longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' P tendency due to vertical transport'
          call hist_addfld_decomp (fname=fieldname, units='gP/m^3/s',  type2d='levdcmp', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       end if
    end do

    this%primp_to_labilep_col(begc:endc) = spval
    call hist_addfld1d (fname='PRIMP_TO_LABILEP', units='gP/m^2/s',   &
         avgflag='A', long_name='PRIMARY MINERAL P TO LABILE P', &
         ptr_col=this%primp_to_labilep_col)

    if ( nlevdecomp_full > 1 ) then  
       this%primp_to_labilep_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='PRIMP_TO_LABILEP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='PRIMARY MINERAL P TO LABILE P', &
            ptr_col=this%primp_to_labilep_vr_col, default='inactive')
    endif

    this%labilep_to_secondp_col(begc:endc) = spval
    call hist_addfld1d (fname='LABILEP_TO_SECONDP', units='gP/m^2/s',   &
         avgflag='A', long_name='LABILE P TO SECONDARY MINERAL P', &
         ptr_col=this%labilep_to_secondp_col)

    if ( nlevdecomp_full > 1 ) then  
       this%labilep_to_secondp_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='LABILEP_TO_SECONDP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='LABILE P TO SECONDARY MINERAL P', &
            ptr_col=this%labilep_to_secondp_vr_col, default='inactive')
    endif


    this%secondp_to_labilep_col(begc:endc) = spval
    call hist_addfld1d (fname='SECONDP_TO_LABILEP', units='gP/m^2/s',   &
         avgflag='A', long_name='SECONDARY MINERAL P TO LABILE P', &
         ptr_col=this%secondp_to_labilep_col)

    if ( nlevdecomp_full > 1 ) then  
       this%secondp_to_labilep_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SECONDP_TO_LABILEP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='SECONDARY MINERAL P TO LABILE P', &
            ptr_col=this%secondp_to_labilep_vr_col, default='inactive')
    endif

    this%secondp_to_occlp_col(begc:endc) = spval
    call hist_addfld1d (fname='SECONDP_TO_OCCLP', units='gP/m^2/s',   &
         avgflag='A', long_name='SECONDARY MINERAL P TO OCCLUDED P', &
         ptr_col=this%secondp_to_occlp_col)

    if ( nlevdecomp_full > 1 ) then  
       this%secondp_to_occlp_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SECONDP_TO_OCCLP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='SECONDARY MINERAL P TO OCCLUDED P', &
            ptr_col=this%secondp_to_occlp_vr_col, default='inactive')
    endif

    this%sminp_leached_col(begc:endc) = spval
    call hist_addfld1d (fname='SMINP_LEACHED', units='gP/m^2/s',   &
         avgflag='A', long_name='soil mineral P pool loss to leaching', &
         ptr_col=this%sminp_leached_col)

    if ( nlevdecomp_full > 1 ) then  
       this%sminp_leached_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINP_LEACHED'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='soil mineral P pool loss to leaching', &
            ptr_col=this%sminp_leached_vr_col, default='inactive')
    endif


    if ( nlevdecomp_full > 1 ) then
       this%potential_immob_p_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='POTENTIAL_IMMOB_P'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='potential P immobilization', &
            ptr_col=this%potential_immob_p_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%actual_immob_p_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='ACTUAL_IMMOB_P'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='actual P immobilization', &
            ptr_col=this%actual_immob_p_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%sminp_to_plant_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SMINP_TO_PLANT'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='plant uptake of soil mineral P', &
            ptr_col=this%sminp_to_plant_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%supplement_to_sminp_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='SUPPLEMENT_TO_SMINP'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='supplemental P supply', &
            ptr_col=this%supplement_to_sminp_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%gross_pmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='GROSS_PMIN'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='gross rate of P mineralization', &
            ptr_col=this%gross_pmin_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%net_pmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='NET_PMIN'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='net rate of P mineralization', &
            ptr_col=this%net_pmin_vr_col, default='inactive')
    end if

    if ( nlevdecomp_full > 1 ) then
       this%biochem_pmin_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='BIOCHEM_PMIN'//trim(vr_suffix), units='gP/m^3/s',  type2d='levdcmp', &
            avgflag='A', long_name='biochemical rate of P mineralization', &
            ptr_col=this%biochem_pmin_vr_col, default='inactive')
    end if
         
    this%potential_immob_p_col(begc:endc) = spval
    call hist_addfld1d (fname='POTENTIAL_IMMOB_P', units='gP/m^2/s', &
         avgflag='A', long_name='potential P immobilization', &
         ptr_col=this%potential_immob_p_col)

    this%actual_immob_p_col(begc:endc) = spval
    call hist_addfld1d (fname='ACTUAL_IMMOB_P', units='gP/m^2/s', &
         avgflag='A', long_name='actual P immobilization', &
         ptr_col=this%actual_immob_p_col)

    this%sminp_to_plant_col(begc:endc) = spval
    call hist_addfld1d (fname='SMINP_TO_PLANT', units='gP/m^2/s', &
         avgflag='A', long_name='plant uptake of soil mineral P', &
         ptr_col=this%sminp_to_plant_col)

    this%supplement_to_sminp_col(begc:endc) = spval
    call hist_addfld1d (fname='SUPPLEMENT_TO_SMINP', units='gP/m^2/s', &
         avgflag='A', long_name='supplemental P supply', &
         ptr_col=this%supplement_to_sminp_col)

    this%gross_pmin_col(begc:endc) = spval
    call hist_addfld1d (fname='GROSS_PMIN', units='gP/m^2/s', &
         avgflag='A', long_name='gross rate of P mineralization', &
         ptr_col=this%gross_pmin_col)

    this%net_pmin_col(begc:endc) = spval
    call hist_addfld1d (fname='NET_PMIN', units='gP/m^2/s', &
         avgflag='A', long_name='net rate of P mineralization', &
         ptr_col=this%net_pmin_col)

    this%biochem_pmin_col(begc:endc) = spval
    call hist_addfld1d (fname='BIOCHEM_PMIN', units='gP/m^2/s', &
         avgflag='A', long_name='biochemical rate of P mineralization', &
         ptr_col=this%biochem_pmin_col)

    this%fire_ploss_col(begc:endc) = spval
    call hist_addfld1d (fname='COL_FIRE_PLOSS', units='gP/m^2/s', &
         avgflag='A', long_name='total column-level fire P loss', &
         ptr_col=this%fire_ploss_col, default='inactive')

    this%dwt_seedp_to_leaf_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_LEAF', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level leaf', &
         ptr_col=this%dwt_seedp_to_leaf_col, default='inactive')

    this%dwt_seedp_to_deadstem_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_SEEDP_TO_DEADSTEM', units='gP/m^2/s', &
         avgflag='A', long_name='seed source to PFT-level deadstem', &
         ptr_col=this%dwt_seedp_to_deadstem_col, default='inactive')

    this%dwt_conv_pflux_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_CONV_PFLUX', units='gP/m^2/s', &
         avgflag='A', long_name='conversion P flux (immediate loss to atm)', &
         ptr_col=this%dwt_conv_pflux_col, default='inactive')

    this%dwt_prod10p_gain_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_PROD10P_GAIN', units='gP/m^2/s', &
         avgflag='A', long_name='addition to 10-yr wood product pool', &
         ptr_col=this%dwt_prod10p_gain_col, default='inactive')

    this%prod10p_loss_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD10P_LOSS', units='gP/m^2/s', &
         avgflag='A', long_name='loss from 10-yr wood product pool', &
         ptr_col=this%prod10p_loss_col, default='inactive')

    this%dwt_prod100p_gain_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_PROD100P_GAIN', units='gP/m^2/s', &
         avgflag='A', long_name='addition to 100-yr wood product pool', &
         ptr_col=this%dwt_prod100p_gain_col, default='inactive')

    this%prod100p_loss_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD100P_LOSS', units='gP/m^2/s', &
         avgflag='A', long_name='loss from 100-yr wood product pool', &
         ptr_col=this%prod100p_loss_col, default='inactive')

    this%prod1p_loss_col(begc:endc) = spval
    call hist_addfld1d (fname='PROD1P_LOSS', units='gP/m^2/s', &
         avgflag='A', long_name='loss from 1-yr crop product pool', &
         ptr_col=this%prod1p_loss_col)

    this%product_ploss_col(begc:endc) = spval
    call hist_addfld1d (fname='PRODUCT_PLOSS', units='gP/m^2/s', &
         avgflag='A', long_name='total P loss from wood product pools', &
         ptr_col=this%product_ploss_col, default='inactive')

    this%dwt_frootp_to_litr_met_p_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_FROOTP_TO_LITR_MET_P', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=this%dwt_frootp_to_litr_met_p_col, default='inactive')

    this%dwt_frootp_to_litr_cel_p_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_FROOTP_TO_LITR_CEL_P', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=this%dwt_frootp_to_litr_cel_p_col, default='inactive')

    this%dwt_frootp_to_litr_lig_p_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_FROOTP_TO_LITR_LIG_P', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='fine root to litter due to landcover change', &
         ptr_col=this%dwt_frootp_to_litr_lig_p_col, default='inactive')

    this%dwt_livecrootp_to_cwdp_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_LIVECROOTP_TO_CWDP', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='live coarse root to CWD due to landcover change', &
         ptr_col=this%dwt_livecrootp_to_cwdp_col, default='inactive')

    this%dwt_deadcrootp_to_cwdp_col(begc:endc,:) = spval
    call hist_addfld_decomp (fname='DWT_DEADCROOTP_TO_CWDP', units='gP/m^2/s',  type2d='levdcmp', &
         avgflag='A', long_name='dead coarse root to CWD due to landcover change', &
         ptr_col=this%dwt_deadcrootp_to_cwdp_col, default='inactive')

    this%dwt_ploss_col(begc:endc) = spval
    call hist_addfld1d (fname='DWT_PLOSS', units='gP/m^2/s', &
         avgflag='A', long_name='total phosphorus loss from landcover conversion', &
         ptr_col=this%dwt_ploss_col, default='inactive')

    if (crop_prog) then
       this%fert_p_to_sminp_col(begc:endc) = spval
       call hist_addfld1d (fname='FERT_TO_LABILEP', units='gP/m^2/s', &
            avgflag='A', long_name='fertilizer to soil mineral P', &
            ptr_col=this%fert_p_to_sminp_col)
    end if

    !! bgc interface
    this%plant_pdemand_col(begc:endc) = spval
    call hist_addfld1d (fname='PLANT_PDEMAND_COL', units='gN/m^2/s', &
        avgflag='A', long_name='P flux required to support initial GPP', &
        ptr_col=this%plant_pdemand_col)

    this%adsorb_to_labilep_col(begc:endc) = spval
    call hist_addfld1d (fname='ADSORBTION_P', units='gP/m^2/s', &
         avgflag='A', long_name='adsorb P flux', &
         ptr_col=this%adsorb_to_labilep_col, default='active')
         
    this%desorb_to_solutionp_col(begc:endc) = spval
    call hist_addfld1d (fname='DESORPTION_P', units='gP/m^2/s', &
         avgflag='A', long_name='desorp P flux', &
         ptr_col=this%desorb_to_solutionp_col, default='active')
         
  end subroutine inithistory_col_pf

  !-----------------------------------------------------------------------



  subroutine initallocate_col_pf(this, bounds)
    class(soilcol_phosphorus_flux) :: this
 		
 		type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begp,endp
    integer           :: begc,endc
    !------------------------------------------------------------------------
   
    begc = bounds%begc; endc = bounds%endc   

    allocate(this%pdep_to_sminp_col             (begc:endc))    ; this%pdep_to_sminp_col     (:) = nan
    allocate(this%fert_p_to_sminp_col             (begc:endc))    ; this%fert_p_to_sminp_col     (:) = nan
    allocate(this%hrv_deadstemp_to_prod10p_col  (begc:endc))    ; this%hrv_deadstemp_to_prod10p_col  (:) = nan
    allocate(this%hrv_deadstemp_to_prod100p_col (begc:endc))    ; this%hrv_deadstemp_to_prod100p_col (:) = nan
    allocate(this%hrv_cropp_to_prod1p_col       (begc:endc))    ; this%hrv_cropp_to_prod1p_col       (:) = nan
    allocate(this%sminp_to_plant_col            (begc:endc))    ; this%sminp_to_plant_col     (:) = nan
    allocate(this%potential_immob_p_col           (begc:endc))    ; this%potential_immob_p_col           (:) = nan
    allocate(this%actual_immob_p_col              (begc:endc))    ; this%actual_immob_p_col              (:) = nan
    allocate(this%gross_pmin_col                (begc:endc))    ; this%gross_pmin_col                (:) = nan
    allocate(this%net_pmin_col                  (begc:endc))    ; this%net_pmin_col                  (:) = nan
    allocate(this%supplement_to_sminp_col       (begc:endc))    ; this%supplement_to_sminp_col       (:) = nan
    allocate(this%prod1p_loss_col               (begc:endc))    ; this%prod1p_loss_col              (:) = nan
    allocate(this%prod10p_loss_col              (begc:endc))    ; this%prod10p_loss_col              (:) = nan
    allocate(this%prod100p_loss_col             (begc:endc))    ; this%prod100p_loss_col     (:) = nan
    allocate(this%product_ploss_col             (begc:endc))    ; this%product_ploss_col     (:) = nan
    allocate(this%pinputs_col                   (begc:endc))    ; this%pinputs_col                   (:) = nan
    allocate(this%poutputs_col                  (begc:endc))    ; this%poutputs_col                  (:) = nan
    allocate(this%fire_ploss_col                (begc:endc))    ; this%fire_ploss_col                (:) = nan
    allocate(this%fire_ploss_p2c_col            (begc:endc))    ; this%fire_ploss_p2c_col            (:) = nan
    allocate(this%som_p_leached_col             (begc:endc))    ; this%som_p_leached_col     (:) = nan

    allocate(this%m_p_to_litr_met_fire_col   (begc:endc,1:nlevdecomp_full)) ; this%m_p_to_litr_met_fire_col   (:,:) = nan
    allocate(this%m_p_to_litr_cel_fire_col   (begc:endc,1:nlevdecomp_full)) ; this%m_p_to_litr_cel_fire_col   (:,:) = nan
    allocate(this%m_p_to_litr_lig_fire_col   (begc:endc,1:nlevdecomp_full)) ; this%m_p_to_litr_lig_fire_col   (:,:) = nan
    allocate(this%potential_immob_p_vr_col   (begc:endc,1:nlevdecomp_full)) ; this%potential_immob_p_vr_col     (:,:) = nan
    allocate(this%actual_immob_p_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%actual_immob_p_vr_col        (:,:) = nan
    allocate(this%sminp_to_plant_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%sminp_to_plant_vr_col      (:,:) = nan
    allocate(this%supplement_to_sminp_vr_col (begc:endc,1:nlevdecomp_full)) ; this%supplement_to_sminp_vr_col (:,:) = nan
    allocate(this%gross_pmin_vr_col          (begc:endc,1:nlevdecomp_full)) ; this%gross_pmin_vr_col          (:,:) = nan
    allocate(this%net_pmin_vr_col            (begc:endc,1:nlevdecomp_full)) ; this%net_pmin_vr_col            (:,:) = nan

    allocate(this%biochem_pmin_ppools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    allocate(this%biochem_pmin_vr_col       (begc:endc,1:nlevdecomp_full))
    allocate(this%biochem_pmin_col          (begc:endc))

    this%biochem_pmin_ppools_vr_col  (:,:,:) = nan
    this%biochem_pmin_vr_col         (:,:) = nan
    this%biochem_pmin_col            (:) = nan

    allocate(this%dwt_seedp_to_leaf_col      (begc:endc))                   ; this%dwt_seedp_to_leaf_col      (:)   = nan
    allocate(this%dwt_seedp_to_deadstem_col  (begc:endc))                   ; this%dwt_seedp_to_deadstem_col  (:)   = nan
    allocate(this%dwt_conv_pflux_col         (begc:endc))                   ; this%dwt_conv_pflux_col         (:)   = nan
    allocate(this%dwt_prod10p_gain_col       (begc:endc))                   ; this%dwt_prod10p_gain_col       (:)   = nan
    allocate(this%dwt_prod100p_gain_col      (begc:endc))                   ; this%dwt_prod100p_gain_col      (:)   = nan
    allocate(this%dwt_ploss_col              (begc:endc))                   ; this%dwt_ploss_col              (:)   = nan
    allocate(this%wood_harvestp_col          (begc:endc))                   ; this%wood_harvestp_col          (:)   = nan

    allocate(this%dwt_frootp_to_litr_met_p_col(begc:endc,1:nlevdecomp_full)) ; this%dwt_frootp_to_litr_met_p_col     (:,:) = nan
    allocate(this%dwt_frootp_to_litr_cel_p_col(begc:endc,1:nlevdecomp_full)) ; this%dwt_frootp_to_litr_cel_p_col     (:,:) = nan
    allocate(this%dwt_frootp_to_litr_lig_p_col(begc:endc,1:nlevdecomp_full)) ; this%dwt_frootp_to_litr_lig_p_col     (:,:) = nan
    allocate(this%dwt_livecrootp_to_cwdp_col  (begc:endc,1:nlevdecomp_full)) ; this%dwt_livecrootp_to_cwdp_col       (:,:) = nan
    allocate(this%dwt_deadcrootp_to_cwdp_col  (begc:endc,1:nlevdecomp_full)) ; this%dwt_deadcrootp_to_cwdp_col       (:,:) = nan


    allocate(this%decomp_cascade_ptransfer_vr_col   (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%decomp_cascade_sminp_flux_vr_col  (begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions ))
    allocate(this%m_decomp_ppools_to_fire_vr_col    (begc:endc,1:nlevdecomp_full,1:ndecomp_pools               ))
    allocate(this%decomp_cascade_ptransfer_col      (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%decomp_cascade_sminp_flux_col     (begc:endc,1:ndecomp_cascade_transitions                   ))
    allocate(this%m_decomp_ppools_to_fire_col       (begc:endc,1:ndecomp_pools                                 ))

    this%decomp_cascade_ptransfer_vr_col  (:,:,:) = nan
    this%decomp_cascade_sminp_flux_vr_col (:,:,:) = nan
    this%m_decomp_ppools_to_fire_vr_col   (:,:,:) = nan
    this%m_decomp_ppools_to_fire_col      (:,:)   = nan
    this%decomp_cascade_ptransfer_col     (:,:)   = nan
    this%decomp_cascade_sminp_flux_col    (:,:)   = nan

    allocate(this%phenology_p_to_litr_met_p_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%phenology_p_to_litr_cel_p_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%phenology_p_to_litr_lig_p_col     (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_p_to_litr_met_p_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_p_to_litr_cel_p_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_p_to_litr_lig_p_col (begc:endc, 1:nlevdecomp_full))
    allocate(this%gap_mortality_p_to_cwdp_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%fire_mortality_p_to_cwdp_col      (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_p_to_litr_met_p_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_p_to_litr_cel_p_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_p_to_litr_lig_p_col       (begc:endc, 1:nlevdecomp_full))
    allocate(this%harvest_p_to_cwdp_col             (begc:endc, 1:nlevdecomp_full))

    this%phenology_p_to_litr_met_p_col     (:,:) = nan
    this%phenology_p_to_litr_cel_p_col     (:,:) = nan
    this%phenology_p_to_litr_lig_p_col     (:,:) = nan
    this%gap_mortality_p_to_litr_met_p_col (:,:) = nan
    this%gap_mortality_p_to_litr_cel_p_col (:,:) = nan
    this%gap_mortality_p_to_litr_lig_p_col (:,:) = nan
    this%gap_mortality_p_to_cwdp_col       (:,:) = nan
    this%fire_mortality_p_to_cwdp_col      (:,:) = nan
    this%harvest_p_to_litr_met_p_col       (:,:) = nan
    this%harvest_p_to_litr_cel_p_col       (:,:) = nan
    this%harvest_p_to_litr_lig_p_col       (:,:) = nan
    this%harvest_p_to_cwdp_col             (:,:) = nan


    
    allocate(this%primp_to_labilep_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%primp_to_labilep_col                    (begc:endc                                                 ))
    allocate(this%labilep_to_secondp_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%labilep_to_secondp_col                    (begc:endc                                                 ))
    allocate(this%secondp_to_labilep_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%secondp_to_labilep_col                    (begc:endc                                                 ))
    allocate(this%secondp_to_occlp_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%secondp_to_occlp_col                    (begc:endc                                                 ))

    this%primp_to_labilep_vr_col                 (:,:)   = nan
    this%primp_to_labilep_col                    (:)     = nan
    this%labilep_to_secondp_vr_col                 (:,:)   = nan
    this%labilep_to_secondp_col                    (:)     = nan
    this%secondp_to_labilep_vr_col                 (:,:)   = nan
    this%secondp_to_labilep_col                    (:)     = nan
    this%secondp_to_occlp_vr_col                 (:,:)   = nan
    this%secondp_to_occlp_col                    (:)     = nan

    allocate(this%sminp_leached_vr_col                 (begc:endc,1:nlevdecomp_full                               ))
    allocate(this%sminp_leached_col                    (begc:endc                                                 ))
    allocate(this%decomp_ppools_leached_col            (begc:endc,1:ndecomp_pools                                 ))
    allocate(this%decomp_ppools_transport_tendency_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools               ))

    this%sminp_leached_vr_col                 (:,:)   = nan
    this%sminp_leached_col                    (:)     = nan
    this%decomp_ppools_leached_col            (:,:)   = nan
    this%decomp_ppools_transport_tendency_col (:,:,:) = nan  

    allocate(this%decomp_ppools_sourcesink_col (begc:endc,1:nlevdecomp_full,1:ndecomp_pools))
    this%decomp_ppools_sourcesink_col (:,:,:) = nan

    allocate(this%adsorb_to_labilep_col       (begc:endc                   )) ; this%adsorb_to_labilep_col       (:)   = nan
    allocate(this%desorb_to_solutionp_col     (begc:endc                   )) ; this%desorb_to_solutionp_col     (:)   = nan

    allocate(this%smin_p_to_plant_col         (begc:endc                   )) ; this%smin_p_to_plant_col         (:)   = nan

    ! clm_bgc_interface & pflotran
    !------------------------------------------------------------------------
    allocate(this%plant_pdemand_col                 (begc:endc))                                    ; this%plant_pdemand_col                 (:)     = nan
    allocate(this%plant_pdemand_vr_col              (begc:endc,1:nlevdecomp_full))                  ; this%plant_pdemand_vr_col (:,:) = nan
    allocate(this%externalp_to_decomp_ppools_col    (begc:endc, 1:nlevdecomp_full, 1:ndecomp_pools)); this%externalp_to_decomp_ppools_col    (:,:,:) = spval
    allocate(this%externalp_to_decomp_delta_col     (begc:endc))                                    ; this%externalp_to_decomp_delta_col     (:)     = spval
    allocate(this%sminp_net_transport_vr_col        (begc:endc, 1:nlevdecomp_full))                 ; this%sminp_net_transport_vr_col        (:,:)   = spval
    allocate(this%sminp_net_transport_delta_col     (begc:endc))                                    ; this%sminp_net_transport_delta_col     (:)     = spval

  end subroutine initallocate_col_pf  

end ColumnPhosphorusFluxType