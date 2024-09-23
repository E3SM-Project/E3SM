module CNFireEmissionsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Gathers carbon emissions from fire sources to be sent to CAM-Chem via
  ! the coupler .... 
  ! Created by F. Vitt, and revised by F. Li
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
!  use PatchType,    only : patch                
  use GridcellType, only : grc_pp                
  use ColumnType,   only : col_pp                
  use VegetationType, only : veg_pp                
  use decompMod,    only : bounds_type
  use shr_fire_emis_mod,  only : shr_fire_emis_comps_n, shr_fire_emis_comp_t, shr_fire_emis_linkedlist
  use shr_fire_emis_mod,  only : shr_fire_emis_mechcomps_n, shr_fire_emis_mechcomps
  use spmdMod,            only : masterproc
  use elm_varctl,         only : iulog
  use elm_varcon        , only : spval
  use VegetationDataType     , only : veg_cs, veg_cf 
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
!      write(iulog,*) '  shr_fire_emis_factors_file  = ', shr_fire_emis_factors_file
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

!          call hist_addfld1d (fname='FireComp_'//trim(emis_cmp%name), units='kg/m2/sec', &
!               avgflag='A', long_name='fire emissions flux of '//trim(emis_cmp%name), &
!               ptr_patch=this%comp(icomp)%emis, default='inactive')
          this%comp(icomp)%emis(begp:endp) = spval
          call hist_addfld1d (fname='FireComp_'//trim(emis_cmp%name), units='kg/m2/sec', &
               avgflag='A', long_name='fire emissions flux of '//trim(emis_cmp%name), &
               ptr_patch=this%comp(icomp)%emis)

          emis_cmp => emis_cmp%next_emiscomp

       enddo emis_cmp_loop


       ! loop over atm chem mechanism species
       do imech = 1,shr_fire_emis_mechcomps_n

!          call hist_addfld1d (fname='FireMech_'//trim(shr_fire_emis_mechcomps(imech)%name), units='kg/m2/sec', &
!               avgflag='A', long_name='fire emissions flux of '//trim(shr_fire_emis_mechcomps(imech)%name), &
!               ptr_patch=this%mech(imech)%emis, default='inactive')

          this%mech(imech)%emis(begp:endp) = spval
	  call hist_addfld1d (fname='FireMech_'//trim(shr_fire_emis_mechcomps(imech)%name), units='kg/m2/sec', &
               avgflag='A', long_name='fire emissions flux of '//trim(shr_fire_emis_mechcomps(imech)%name), &
               ptr_patch=this%mech(imech)%emis)

       enddo

!       call hist_addfld1d (fname='FireEmis_TOT', units='gC/m2/sec', &
!            avgflag='A', long_name='Total fire emissions flux ', &
!            ptr_patch=this%totfire%emis, default='inactive')
       this%totfire%emis(begp:endp) = spval
       call hist_addfld1d (fname='FireEmis_TOT', units='gC/m2/sec', &
            avgflag='A', long_name='Total fire emissions flux ', &
            ptr_patch=this%totfire%emis)

!       call hist_addfld1d (fname='FireEmis_ZTOP', units='m', &
!            avgflag='A', long_name='Top of vertical fire emissions distribution ', &
!            ptr_patch=this%ztop_patch, default='inactive')

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
  subroutine CNFireEmisUpdate(bounds, num_soilp, filter_soilp, cnveg_cs_vars, fireemis_vars )
!LXu@02/20-----

!    use CNVegcarbonfluxType,  only : cnveg_carbonflux_type
!    use CNVegCarbonStateType, only : cnveg_carbonstate_type 
!    use VegetationDataType     , only : vegetation_carbon_flux
!    use VegetationDataType     , only : vegetation_carbon_state      
!    use CNCarbonFluxType       , only : carbonflux_type
!    use CNCarbonStateType      , only : carbonstate_type
    use CNStateType,          only : cnstate_type
    use elm_varpar,           only : ndecomp_pools, nlevdecomp
    use elm_varcon,           only : dzsoi_decomp

    !ARGUMENTS:
    type(bounds_type),           intent(in)     :: bounds                  
    integer,                     intent(in)     :: num_soilp       ! number of soil pfts in filter
    integer,                     intent(in)     :: filter_soilp(:) ! filter for soil pfts
!    integer,                     intent(in)    :: filter_soilp(num_solip) ! filter for soil pfts
!    type(carbonflux_type),intent(in)     :: cnveg_cf_vars
!    type(carbonstate_type),intent(in)    :: cnveg_cs_vars 
    type(cnstate_type),          intent(in)    :: cnveg_cs_vars 
    type(fireemis_type),         intent(inout) :: fireemis_vars

    !LOCAL VARIABLES:
    real(r8) :: fire_flux
    real(r8) :: fire_flux_lf 
    real(r8) :: fire_flux_lf1 
    type(shr_fire_emis_comp_t), pointer :: emis_cmp
    real(r8) :: emis_flux(shr_fire_emis_comps_n)
    integer  :: fp,p,g,c                ! indices
    real(r8) :: epsilon                 ! emission factor [ug m-2 h-1]
    integer  :: i, ii, icomp, imech, n_emis_comps, l, j

!LXu@05/20+++++
    real(r8)              :: dmr                    ! the ratio of DM/Carbon_flux from fire emissions [kg/kg]
    real(r8), parameter   :: eqas_latS = -10.0_r8   ! Latitude for Equatorial Asia peat fires
    real(r8), parameter   :: eqas_latN =   8.0_r8   ! Latitude for Equatorial Asia peat fires
    real(r8), parameter   :: eqas_lonL =  95.0_r8   ! Longitude for Equatorial Asia peat fires
    real(r8), parameter   :: eqas_lonR = 160.0_r8   ! Longitude for Equatorial Asia peat fires
!    real(r8), parameter   :: ef_peat(6) = (/0.10_r8, 14.2_r8, 0.1075_r8, 4.3_r8, 3.88_r8, 335.4_r8/)   ! BC,OC,SO4,SO2,SOAG,CO g/kg(DM)
    real(r8), parameter   :: ef_peat(8) = (/0.10_r8, 14.2_r8, 0.1075_r8, 4.3_r8, 3.88_r8, 335.4_r8, 0.0104_r8, 0.0337_r8/)   ! BC,OC,SO4,SO2,SOAG,CO,PO4_a4,PO4_a3 g/kg(DM)
    real(r8), parameter   :: dmr_peat   = 0.57_r8   ! C/DM ratio
    real(r8), parameter   :: dm_ratio(16) = &     
        				 (/ 0.50_r8, 0.49_r8, 0.49_r8, 0.50_r8, 0.50_r8, 0.50_r8, 0.50_r8, 0.49_r8,  &
        				    0.50_r8, 0.50_r8, 0.49_r8, 0.49_r8, 0.49_r8, 0.49_r8, 0.44_r8, 0.44_r8 /)
!LXu@05/20-----

    if ( shr_fire_emis_mechcomps_n < 1) return

    associate( & 
         fire_emis => fireemis_vars%fireflx_patch, &
         totfire   => fireemis_vars%totfire, &
         mech      => fireemis_vars%mech, &
         comp      => fireemis_vars%comp, &
         ztop      => fireemis_vars%ztop_patch &
         )

      ! initialize to zero ...
      fire_emis(bounds%begp:bounds%endp,:) = 0._r8
      totfire%emis(bounds%begp:bounds%endp) =  0._r8
      ztop(bounds%begp:bounds%endp) =  0._r8

      do i = 1, shr_fire_emis_mechcomps_n
         mech(i)%emis(bounds%begp:bounds%endp) =  0._r8
      enddo

      do i = 1, shr_fire_emis_comps_n
         comp(i)%emis(bounds%begp:bounds%endp) =  0._r8
      enddo

      ! Begin loop over points
      !_______________________________________________________________________________
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
!oringional in CESM2
! vertically-resolved decomposing C fire loss   
!            if(cnveg_cs_vars%totvegc_col(c) > 0._r8)then
!               fire_flux_lf1=0._r8 
!               do l = 1, ndecomp_pools
!                  do j = 1, nlevdecomp
!                     fire_flux_lf1 = fire_flux_lf1 + &
!                          cnveg_cf_vars%m_decomp_cpools_to_fire_vr_col(c,j,l)*dzsoi_decomp(j)
!                  enddo
!               end do
!               fire_flux_lf = fire_flux_lf1*cnveg_cs_vars%totvegc_patch(p)/cnveg_cs_vars%totvegc_col(c)
!            else
!               fire_flux_lf=0._r8 
!            end if
!               fire_flux_lf=0._r8 
!            fire_flux = 
!                 + cnveg_cf_vars%m_leafc_to_fire_patch                     (p) & ! (gC/m2/s) fire C emissions from leafc

!             write(iulog,*) cnveg_cf_vars%m_leafc_to_fire_patch(p),      &
!			    cnveg_cf_vars%m_livestemc_to_fire_patch(p),  &
!			    cnveg_cf_vars%m_deadstemc_to_fire_patch(p),  &
!			    cnveg_cf_vars%m_frootc_to_fire_patch(p),     &
!			    cnveg_cf_vars%m_livecrootc_to_fire_patch(p), &
!			    cnveg_cf_vars%m_deadcrootc_to_fire_patch(p), &
!			    cnveg_cf_vars%m_cpool_to_fire_patch(p)
            fire_flux =  &
                   veg_cf%m_leafc_to_fire                     (p) & ! (gC/m2/s) fire C emissions from leafc
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
                 + veg_cf%m_gresp_xfer_to_fire                (p) & ! (gC/m2/s) fire C emissions from gresp_xfer 
                 + veg_cf%m_cpool_to_fire                     (p)   ! (gC/m2/s) fire C emissions from gresp_xfer 

!            fire_flux = fire_flux_lf &
!                 + cnveg_cf_vars%m_leafc_to_fire_patch                     (p) !& 
!                 + cnveg_cf_vars%m_leafc_xfer_to_fire_patch                (p) & 
!                 + cnveg_cf_vars%m_livestemc_to_fire_patch                 (p) & 
!                 + cnveg_cf_vars%m_livestemc_storage_to_fire_patch         (p) & 
!                 + cnveg_cf_vars%m_livestemc_xfer_to_fire_patch            (p) & 
!                 + cnveg_cf_vars%m_deadstemc_to_fire_patch                 (p) & 
!                 + cnveg_cf_vars%m_deadstemc_storage_to_fire_patch         (p) & 
!                 + cnveg_cf_vars%m_deadstemc_xfer_to_fire_patch            (p) & 
!                 + cnveg_cf_vars%m_frootc_to_fire_patch                    (p) & 
!                 + cnveg_cf_vars%m_frootc_storage_to_fire_patch            (p) & 
!                 + cnveg_cf_vars%m_frootc_xfer_to_fire_patch               (p) & 
!                 + cnveg_cf_vars%m_livecrootc_to_fire_patch                (p) & 
!                 + cnveg_cf_vars%m_livecrootc_storage_to_fire_patch        (p) & 
!                 + cnveg_cf_vars%m_livecrootc_xfer_to_fire_patch           (p) & 
!                 + cnveg_cf_vars%m_deadcrootc_to_fire_patch                (p) & 
!                 + cnveg_cf_vars%m_deadcrootc_storage_to_fire_patch        (p) & 
!                 + cnveg_cf_vars%m_deadcrootc_xfer_to_fire_patch           (p) & 
!                 + cnveg_cf_vars%m_gresp_storage_to_fire_patch             (p) & 
!                 + cnveg_cf_vars%m_gresp_xfer_to_fire_patch                (p) & 
!                 + cnveg_cf_vars%m_cpool_to_fire_patch                     (p)    
!            fire_flux = fire_flux_lf  &
!                 + cnveg_cf_vars%m_leafc_to_fire_patch                     (p) &
!                 + cnveg_cf_vars%m_cpool_to_fire_patch                     (p)    
!            fire_flux = 0._r8
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

	       !updated the fire emission in the Indonesia region
               if (      grc_pp%latdeg(g) > eqas_latS .and. grc_pp%latdeg(g) < eqas_latN  &
	           .and. grc_pp%londeg(g) > eqas_lonL .and. grc_pp%londeg(g) < eqas_lonR   ) then
	       
        	  epsilon = ef_peat(icomp)*1.e-3_r8 ! convert g/kg dry fuel to kg/kg(DM)
		  dmr = dmr_peat
        	  comp(icomp)%emis(p) = epsilon * (fire_flux/dmr) * 1.e-3_r8 ! EF * DM (to convert gC/m2/sec to kg species/m2/sec) 
	       
	       else

        	  epsilon = emis_cmp%emis_factors(veg_pp%itype(p))
		  dmr = dm_ratio(veg_pp%itype(p))
        	  comp(icomp)%emis(p) = epsilon * (fire_flux/dmr) * 1.e-3_r8 ! EF * DM (to convert gC/m2/sec to kg species/m2/sec) 
	       
	       end if

               emis_flux(icomp) = emis_cmp%coeff*comp(icomp)%emis(p)
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

         end if ! ivt(1:15 only)

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

