module SurfaceRadiationMod

  !------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar fluxes absorbed by vegetation and ground surface
  !
  ! !USES:
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use clm_varctl        , only : use_snicar_frc, use_ed
  use decompMod         , only : bounds_type
  use clm_varcon        , only : namec
  use atm2lndType       , only : atm2lnd_type
  use WaterstateType    , only : waterstate_type
  use CanopyStateType   , only : canopystate_type
  use SurfaceAlbedoType , only : surfalb_type
  use SolarAbsorbedType , only : solarabs_type
  use GridcellType      , only : grc                
  use LandunitType      , only : lun                
  use ColumnType        , only : col                
  use PatchType         , only : pft                
  use EDVecPatchtype    , only : EDpft
  use EDtypesMod

  !
  ! !PRIVATE TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SurfaceRadiation         ! Solar fluxes absorbed by veg and ground surface
  !
  ! !PRIVATE DATA:
  type, public :: surfrad_type
     real(r8), pointer, private  :: sfc_frc_aer_patch     (:) ! patch surface forcing of snow with all aerosols (pft) [W/m2]
     real(r8), pointer, private  :: sfc_frc_bc_patch      (:) ! patch surface forcing of snow with BC (pft) [W/m2]
     real(r8), pointer, private  :: sfc_frc_oc_patch      (:) ! patch surface forcing of snow with OC (pft) [W/m2]
     real(r8), pointer, private  :: sfc_frc_dst_patch     (:) ! patch surface forcing of snow with dust (pft) [W/m2]
     real(r8), pointer, private  :: sfc_frc_aer_sno_patch (:) ! patch surface forcing of snow with all aerosols, averaged only when snow is present (pft) [W/m2]
     real(r8), pointer, private  :: sfc_frc_bc_sno_patch  (:) ! patch surface forcing of snow with BC, averaged only when snow is present (pft) [W/m2]
     real(r8), pointer, private  :: sfc_frc_oc_sno_patch  (:) ! patch surface forcing of snow with OC, averaged only when snow is present (pft) [W/m2]
     real(r8), pointer, private  :: sfc_frc_dst_sno_patch (:) ! patch surface forcing of snow with dust, averaged only when snow is present (pft) [W/m2]

     real(r8), pointer, private  :: parveg_ln_patch       (:) ! patch  absorbed par by vegetation at local noon (W/m**2)

     real(r8), pointer, private  :: fsr_sno_vd_patch      (:) ! patch reflected direct beam vis solar radiation from snow (W/m**2)
     real(r8), pointer, private  :: fsr_sno_nd_patch      (:) ! patch reflected direct beam NIR solar radiation from snow (W/m**2)
     real(r8), pointer, private  :: fsr_sno_vi_patch      (:) ! patch reflected diffuse vis solar radiation from snow (W/m**2)
     real(r8), pointer, private  :: fsr_sno_ni_patch      (:) ! patch reflected diffuse NIR solar radiation from snow (W/m**2)

     real(r8), pointer, private  :: fsr_vis_d_patch       (:) ! patch reflected direct beam vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsr_vis_i_patch       (:) ! patch reflected diffuse vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsr_vis_d_ln_patch    (:) ! patch reflected direct beam vis solar radiation at local noon (W/m**2)

     real(r8), pointer, private  :: fsds_sno_vd_patch     (:) ! patch incident visible, direct radiation on snow  (for history files)  [W/m2]
     real(r8), pointer, private  :: fsds_sno_nd_patch     (:) ! patch incident near-IR, direct radiation on snow  (for history files)  [W/m2]
     real(r8), pointer, private  :: fsds_sno_vi_patch     (:) ! patch incident visible, diffuse radiation on snow (for history files) [W/m2]
     real(r8), pointer, private  :: fsds_sno_ni_patch     (:) ! patch incident near-IR, diffuse radiation on snow (for history files) [W/m2]

     real(r8), pointer, private  :: fsds_vis_d_patch      (:) ! patch incident direct beam vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsds_vis_i_patch      (:) ! patch incident diffuse vis solar radiation (W/m**2)
     real(r8), pointer, private  :: fsds_vis_d_ln_patch   (:) ! patch incident direct beam vis solar radiation at local noon (W/m**2)
     real(r8), pointer, private  :: fsds_vis_i_ln_patch   (:) ! patch incident diffuse beam vis solar radiation at local noon (W/m**2)

   contains

     procedure, public  :: Init
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory  
     procedure, private :: InitCold     

  end type surfrad_type
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%sfc_frc_aer_patch     (begp:endp))              ; this%sfc_frc_aer_patch     (:)   = nan
    allocate(this%sfc_frc_bc_patch      (begp:endp))              ; this%sfc_frc_bc_patch      (:)   = nan
    allocate(this%sfc_frc_oc_patch      (begp:endp))              ; this%sfc_frc_oc_patch      (:)   = nan
    allocate(this%sfc_frc_dst_patch     (begp:endp))              ; this%sfc_frc_dst_patch     (:)   = nan
    allocate(this%sfc_frc_aer_sno_patch (begp:endp))              ; this%sfc_frc_aer_sno_patch (:)   = nan
    allocate(this%sfc_frc_bc_sno_patch  (begp:endp))              ; this%sfc_frc_bc_sno_patch  (:)   = nan
    allocate(this%sfc_frc_oc_sno_patch  (begp:endp))              ; this%sfc_frc_oc_sno_patch  (:)   = nan
    allocate(this%sfc_frc_dst_sno_patch (begp:endp))              ; this%sfc_frc_dst_sno_patch (:)   = nan

    allocate(this%parveg_ln_patch       (begp:endp))              ; this%parveg_ln_patch       (:)   = nan

    allocate(this%fsr_vis_d_patch       (begp:endp))              ; this%fsr_vis_d_patch       (:)   = nan
    allocate(this%fsr_vis_d_ln_patch    (begp:endp))              ; this%fsr_vis_d_ln_patch    (:)   = nan
    allocate(this%fsr_vis_i_patch       (begp:endp))              ; this%fsr_vis_i_patch       (:)   = nan
    allocate(this%fsr_sno_vd_patch      (begp:endp))              ; this%fsr_sno_vd_patch      (:)   = nan
    allocate(this%fsr_sno_nd_patch      (begp:endp))              ; this%fsr_sno_nd_patch      (:)   = nan
    allocate(this%fsr_sno_vi_patch      (begp:endp))              ; this%fsr_sno_vi_patch      (:)   = nan
    allocate(this%fsr_sno_ni_patch      (begp:endp))              ; this%fsr_sno_ni_patch      (:)   = nan

    allocate(this%fsds_vis_d_patch      (begp:endp))              ; this%fsds_vis_d_patch      (:)   = nan
    allocate(this%fsds_vis_i_patch      (begp:endp))              ; this%fsds_vis_i_patch      (:)   = nan
    allocate(this%fsds_vis_d_ln_patch   (begp:endp))              ; this%fsds_vis_d_ln_patch   (:)   = nan
    allocate(this%fsds_vis_i_ln_patch   (begp:endp))              ; this%fsds_vis_i_ln_patch   (:)   = nan
    allocate(this%fsds_sno_vd_patch     (begp:endp))              ; this%fsds_sno_vd_patch     (:)   = nan
    allocate(this%fsds_sno_nd_patch     (begp:endp))              ; this%fsds_sno_nd_patch     (:)   = nan
    allocate(this%fsds_sno_vi_patch     (begp:endp))              ; this%fsds_sno_vi_patch     (:)   = nan
    allocate(this%fsds_sno_ni_patch     (begp:endp))              ; this%fsds_sno_ni_patch     (:)   = nan

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! History fields initialization
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    use clm_varcon    , only : spval
    use histFileMod   , only : hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc

    if (use_snicar_frc) then
       this%sfc_frc_aer_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOAERFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of all aerosols in snow (land) ', &
            ptr_patch=this%sfc_frc_aer_patch, set_urb=spval)

       this%sfc_frc_aer_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOAERFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of all aerosols in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_aer_sno_patch, set_urb=spval)

       this%sfc_frc_bc_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOBCFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of BC in snow (land) ', &
            ptr_patch=this%sfc_frc_bc_patch, set_urb=spval)

       this%sfc_frc_bc_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOBCFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of BC in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_bc_sno_patch, set_urb=spval)

       this%sfc_frc_oc_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOOCFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of OC in snow (land) ', &
            ptr_patch=this%sfc_frc_oc_patch, set_urb=spval)

       this%sfc_frc_oc_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNOOCFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of OC in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_oc_sno_patch, set_urb=spval)

       this%sfc_frc_dst_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNODSTFRCL', units='W/m^2', &
            avgflag='A', long_name='surface forcing of dust in snow (land) ', &
            ptr_patch=this%sfc_frc_dst_patch, set_urb=spval)

       this%sfc_frc_dst_sno_patch(begp:endp) = spval
       call hist_addfld1d (fname='SNODSTFRC2L', units='W/m^2', &
            avgflag='A', long_name='surface forcing of dust in snow, averaged only when snow is present (land)', &
            ptr_patch=this%sfc_frc_dst_sno_patch, set_urb=spval)
    end if

    this%fsds_vis_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation', &
         ptr_patch=this%fsds_vis_d_patch)

    this%fsds_vis_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation', &
         ptr_patch=this%fsds_vis_i_patch)

    this%fsr_vis_d_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation', &
         ptr_patch=this%fsr_vis_d_patch, c2l_scale_type='urbanf')

    this%fsr_vis_i_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation', &
         ptr_patch=this%fsr_vis_i_patch, c2l_scale_type='urbanf')

    this%fsds_vis_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation at local noon', &
         ptr_patch=this%fsds_vis_d_ln_patch)

    this%fsds_vis_i_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSDSVILN', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation at local noon', &
         ptr_patch=this%fsds_vis_i_ln_patch)

    this%parveg_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='PARVEGLN', units='W/m^2',  &
         avgflag='A', long_name='absorbed par by vegetation at local noon', &
         ptr_patch=this%parveg_ln_patch)

    this%fsr_vis_d_ln_patch(begp:endp) = spval
    call hist_addfld1d (fname='FSRVDLN', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation at local noon', &
         ptr_patch=this%fsr_vis_d_ln_patch, c2l_scale_type='urbanf')

    this%fsds_sno_vd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_vd_patch, default='inactive')

    this%fsds_sno_nd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSND', units='W/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_nd_patch, default='inactive')

    this%fsds_sno_vi_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_vi_patch, default='inactive')

    this%fsds_sno_ni_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSDSNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation on snow', &
         ptr_patch=this%fsds_sno_ni_patch, default='inactive')

    this%fsr_sno_vd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRVD', units='W/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_vd_patch, default='inactive')

    this%fsr_sno_nd_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRND', units='W/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_nd_patch, default='inactive')

    this%fsr_sno_vi_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRVI', units='W/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_vi_patch, default='inactive')

    this%fsr_sno_ni_patch(begp:endp) = spval
    call hist_addfld1d (fname='SNOFSRNI', units='W/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation from snow', &
         ptr_patch=this%fsr_sno_ni_patch, default='inactive')

  end subroutine InitHistory

  !------------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(surfrad_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: p,l
    !-----------------------------------------------------------------------

    ! nothing for now

  end subroutine InitCold

  !------------------------------------------------------------------------------
  subroutine SurfaceRadiation(bounds, num_nourbanp, filter_nourbanp, &
       num_urbanp, filter_urbanp, num_urbanc, filter_urbanc, &
       atm2lnd_vars, waterstate_vars, canopystate_vars, surfalb_vars, &
       solarabs_vars, surfrad_vars)
     !
     ! !DESCRIPTION: 
     ! Solar fluxes absorbed by vegetation and ground surface
     ! Note possible problem when land is on different grid than atmosphere.
     ! Land may have sun above the horizon (coszen > 0) but atmosphere may
     ! have sun below the horizon (forc_solad = 0 and forc_solai = 0). This is okay
     ! because all fluxes (absorbed, reflected, transmitted) are multiplied
     ! by the incoming flux and all will equal zero.
     ! Atmosphere may have sun above horizon (forc_solad > 0 and forc_solai > 0) but
     ! land may have sun below horizon. This is okay because fabd, fabi,
     ! ftdd, ftid, and ftii all equal zero so that sabv=sabg=fsa=0. Also,
     ! albd and albi equal one so that fsr=forc_solad+forc_solai. In other words, all
     ! the radiation is reflected. NDVI should equal zero in this case.
     ! However, the way the code is currently implemented this is only true
     ! if (forc_solad+forc_solai)|vis = (forc_solad+forc_solai)|nir.
     ! Output variables are parsun,parsha,sabv,sabg,fsa,fsr,ndvi
     !
     ! !USES:
     use clm_varpar       , only : numrad, nlevsno
     use clm_varcon       , only : spval, degpsec, isecspday
     use landunit_varcon  , only : istsoil, istcrop 
     use clm_varctl       , only : subgridflag, use_snicar_frc, iulog
     use clm_time_manager , only : get_curr_date, get_step_size
     use SnowSnicarMod    , only : DO_SNO_OC
     use abortutils       , only : endrun
     !
     ! !ARGUMENTS:
     type(bounds_type)      , intent(in)    :: bounds             
     integer                , intent(in)    :: num_nourbanp       ! number of patches in non-urban points in pft filter
     integer                , intent(in)    :: filter_nourbanp(:) ! patch filter for non-urban points
     integer                , intent(in)    :: num_urbanp         ! number of patches in non-urban points in pft filter
     integer                , intent(in)    :: filter_urbanp(:)   ! patch filter for non-urban points
     integer                , intent(in)    :: num_urbanc         ! number of urban columns in clump
     integer                , intent(in)    :: filter_urbanc(:)   ! urban column filter
     type(atm2lnd_type)     , intent(in)    :: atm2lnd_vars
     type(waterstate_type)  , intent(in)    :: waterstate_vars
     type(surfalb_type)     , intent(in)    :: surfalb_vars
     type(canopystate_type) , intent(inout) :: canopystate_vars
     type(solarabs_type)    , intent(inout) :: solarabs_vars
     type(surfrad_type)     , intent(inout) :: surfrad_vars
     !
     ! !LOCAL VARIABLES:
     integer , parameter :: nband = numrad           ! number of solar radiation waveband classes
     real(r8), parameter :: mpe = 1.e-06_r8          ! prevents overflow for division by zero
     integer  :: fp                                  ! non-urban filter pft index
     integer  :: p                                   ! patch index
     integer  :: c                                   ! column index
     integer  :: l                                   ! landunit index
     integer  :: g                                   ! grid cell index
     integer  :: ib                                  ! waveband number (1=vis, 2=nir)
     integer  :: iv                                  ! canopy layer
     real(r8) :: absrad                              ! absorbed solar radiation (W/m**2)
     integer  :: i                                   ! layer index [idx]
     real(r8) :: rnir                                ! reflected solar radiation [nir] (W/m**2)
     real(r8) :: rvis                                ! reflected solar radiation [vis] (W/m**2)
     real(r8) :: trd(bounds%begp:bounds%endp,numrad) ! transmitted solar radiation: direct (W/m**2)
     real(r8) :: tri(bounds%begp:bounds%endp,numrad) ! transmitted solar radiation: diffuse (W/m**2)
     real(r8) :: cad(bounds%begp:bounds%endp,numrad) ! direct beam absorbed by canopy (W/m**2)
     real(r8) :: cai(bounds%begp:bounds%endp,numrad) ! diffuse radiation absorbed by canopy (W/m**2)
     integer  :: local_secp1                         ! seconds into current date in local time
     real(r8) :: dtime                               ! land model time step (sec)
     integer  :: year,month,day,secs                 ! calendar info for current time step
     real(r8) :: sabg_snl_sum                        ! temporary, absorbed energy in all active snow layers [W/m2]
     real(r8) :: absrad_pur                          ! temp: absorbed solar radiation by pure snow [W/m2]
     real(r8) :: absrad_bc                           ! temp: absorbed solar radiation without BC [W/m2]
     real(r8) :: absrad_oc                           ! temp: absorbed solar radiation without OC [W/m2]
     real(r8) :: absrad_dst                          ! temp: absorbed solar radiation without dust [W/m2]
     real(r8) :: sabg_pur(bounds%begp:bounds%endp)   ! solar radiation absorbed by ground with pure snow [W/m2]
     real(r8) :: sabg_bc(bounds%begp:bounds%endp)    ! solar radiation absorbed by ground without BC [W/m2]
     real(r8) :: sabg_oc(bounds%begp:bounds%endp)    ! solar radiation absorbed by ground without OC [W/m2]
     real(r8) :: sabg_dst(bounds%begp:bounds%endp)   ! solar radiation absorbed by ground without dust [W/m2]
     real(r8) :: parveg(bounds%begp:bounds%endp)     ! absorbed par by vegetation (W/m**2)
     !
     integer, parameter :: noonsec   = isecspday / 2 ! seconds at local noon
     !
     !ED specific variables 
     real(r8)               :: errsol(bounds%begp:bounds%endp) ! solar radiation error Wm-2
     real(r8)               :: sunlai                          ! intermediate for calculating canopy fsun
     real(r8)               :: shalai                          ! intermediate for calculating canopy fsha
     integer                :: CL                              ! Canopy Layer index
     integer                :: FT                              ! clm patch index
     real                   :: gaib, rib                       ! for debugging
     type (site) ,  pointer :: currentSite 
     type (patch),  pointer :: currentPatch                    ! Import fapar matrix for each patch from ED data structure.
     !------------------------------------------------------------------------------

     associate(                                                     & 
          snl             =>    col%snl                           , & ! Input:  [integer  (:)   ] negative number of snow layers [nbr]     

          forc_solad      =>    atm2lnd_vars%forc_solad_grc       , & ! Input:  [real(r8) (:,:) ] direct beam radiation (W/m**2)        
          forc_solai      =>    atm2lnd_vars%forc_solai_grc       , & ! Input:  [real(r8) (:,:) ] diffuse radiation (W/m**2)            

          snow_depth      =>    waterstate_vars%snow_depth_col    , & ! Input:  [real(r8) (:)   ] snow height (m)                         
          frac_sno        =>    waterstate_vars%frac_sno_col      , & ! Input:  [real(r8) (:)   ] fraction of ground covered by snow (0 to 1)
          
          nrad            =>    surfalb_vars%nrad_patch           , & ! Input:  [integer  (:)   ] number of canopy layers, above snow for radiative transfer
          fsun_z          =>    surfalb_vars%fsun_z_patch         , & ! Input:  [real(r8) (:,:) ] sunlit fraction of canopy layer       
          tlai_z          =>    surfalb_vars%tlai_z_patch         , & ! Input:  [real(r8) (:,:) ] tlai increment for canopy layer       
          tsai_z          =>    surfalb_vars%tsai_z_patch         , & ! Input:  [real(r8) (:,:) ] tsai increment for canopy layer       
          coszen          =>    surfalb_vars%coszen_col           , & ! Input:  [real(r8) (:)   ] column cosine of solar zenith angle            
          albgrd          =>    surfalb_vars%albgrd_col           , & ! Input:  [real(r8) (:,:) ] ground albedo (direct)                
          albgri          =>    surfalb_vars%albgri_col           , & ! Input:  [real(r8) (:,:) ] ground albedo (diffuse)               
          albsod          =>    surfalb_vars%albsod_col           , & ! Input:  [real(r8) (:,:) ] direct-beam soil albedo (col,bnd) [frc]
          albgrd_oc       =>    surfalb_vars%albgrd_oc_col        , & ! Input:  [real(r8) (:,:) ] ground albedo without OC (direct) (col,bnd)
          albgri_oc       =>    surfalb_vars%albgri_oc_col        , & ! Input:  [real(r8) (:,:) ] ground albedo without OC (diffuse) (col,bnd)
          albgrd_dst      =>    surfalb_vars%albgrd_dst_col       , & ! Input:  [real(r8) (:,:) ] ground albedo without dust (direct) (col,bnd)
          albgri_dst      =>    surfalb_vars%albgri_dst_col       , & ! Input:  [real(r8) (:,:) ] ground albedo without dust (diffuse) (col,bnd)
          albsnd_hst      =>    surfalb_vars%albsnd_hst_col       , & ! Input:  [real(r8) (:,:) ] snow albedo, direct, for history files (col,bnd) [frc]
          albsni_hst      =>    surfalb_vars%albsni_hst_col       , & ! Input:  [real(r8) (:,:) ] snow ground albedo, diffuse, for history files (col,bnd
          flx_absdv       =>    surfalb_vars%flx_absdv_col        , & ! Input:  [real(r8) (:,:) ] direct flux absorption factor (col,lyr): VIS [frc] 
          flx_absdn       =>    surfalb_vars%flx_absdn_col        , & ! Input:  [real(r8) (:,:) ] direct flux absorption factor (col,lyr): NIR [frc]
          flx_absiv       =>    surfalb_vars%flx_absiv_col        , & ! Input:  [real(r8) (:,:) ] diffuse flux absorption factor (col,lyr): VIS [frc]
          flx_absin       =>    surfalb_vars%flx_absin_col        , & ! Input:  [real(r8) (:,:) ] diffuse flux absorption factor (col,lyr): NIR [frc]
          albsoi          =>    surfalb_vars%albsoi_col           , & ! Input:  [real(r8) (:,:) ] diffuse soil albedo (col,bnd) [frc] 
          albd            =>    surfalb_vars%albd_patch           , & ! Input:  [real(r8) (:,:) ] surface albedo (direct)               
          albi            =>    surfalb_vars%albi_patch           , & ! Input:  [real(r8) (:,:) ] surface albedo (diffuse)              
          fabd            =>    surfalb_vars%fabd_patch           , & ! Input:  [real(r8) (:,:) ] flux absorbed by canopy per unit direct flux
          fabd_sun        =>    surfalb_vars%fabd_sun_patch       , & ! Input:  [real(r8) (:,:) ] flux absorbed by sunlit canopy per unit direct flux
          fabd_sha        =>    surfalb_vars%fabd_sha_patch       , & ! Input:  [real(r8) (:,:) ] flux absorbed by shaded canopy per unit direct flux
          fabi            =>    surfalb_vars%fabi_patch           , & ! Input:  [real(r8) (:,:) ] flux absorbed by canopy per unit diffuse flux
          fabi_sun        =>    surfalb_vars%fabi_sun_patch       , & ! Input:  [real(r8) (:,:) ] flux absorbed by sunlit canopy per unit diffuse flux
          fabi_sha        =>    surfalb_vars%fabi_sha_patch       , & ! Input:  [real(r8) (:,:) ] flux absorbed by shaded canopy per unit diffuse flux
          ftdd            =>    surfalb_vars%ftdd_patch           , & ! Input:  [real(r8) (:,:) ] down direct flux below canopy per unit direct flux
          ftid            =>    surfalb_vars%ftid_patch           , & ! Input:  [real(r8) (:,:) ] down diffuse flux below canopy per unit direct flux
          ftii            =>    surfalb_vars%ftii_patch           , & ! Input:  [real(r8) (:,:) ] down diffuse flux below canopy per unit diffuse flux
          fabd_sun_z      =>    surfalb_vars%fabd_sun_z_patch     , & ! Input:  [real(r8) (:,:) ] absorbed sunlit leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabd_sha_z      =>    surfalb_vars%fabd_sha_z_patch     , & ! Input:  [real(r8) (:,:) ] absorbed shaded leaf direct  PAR (per unit lai+sai) for each canopy layer
          fabi_sun_z      =>    surfalb_vars%fabi_sun_z_patch     , & ! Input:  [real(r8) (:,:) ] absorbed sunlit leaf diffuse PAR (per unit lai+sai) for each canopy layer
          fabi_sha_z      =>    surfalb_vars%fabi_sha_z_patch     , & ! Input:  [real(r8) (:,:) ] absorbed shaded leaf diffuse PAR (per unit lai+sai) for each canopy layer
          albgrd_pur      =>    surfalb_vars%albgrd_pur_col       , & ! Input:  [real(r8) (:,:) ] pure snow ground albedo (direct)      
          albgri_pur      =>    surfalb_vars%albgri_pur_col       , & ! Input:  [real(r8) (:,:) ] pure snow ground albedo (diffuse)     
          albgrd_bc       =>    surfalb_vars%albgrd_bc_col        , & ! Input:  [real(r8) (:,:) ] ground albedo without BC (direct) (col,bnd)
          albgri_bc       =>    surfalb_vars%albgri_bc_col        , & ! Input:  [real(r8) (:,:) ] ground albedo without BC (diffuse) (col,bnd)
          tlai            =>    canopystate_vars%tlai_patch       , & ! Input:  [real(r8) (:)   ] one-sided leaf area index
          elai            =>    canopystate_vars%elai_patch       , & ! Input:  [real(r8) (:)   ] one-sided leaf area index with burying by snow
          esai            =>    canopystate_vars%esai_patch       , & ! Input:  [real(r8) (:)   ] one-sided stem area index with burying by snow
          laisun          =>    canopystate_vars%laisun_patch     , & ! Output: [real(r8) (:)   ] sunlit leaf area                        
          laisha          =>    canopystate_vars%laisha_patch     , & ! Output: [real(r8) (:)   ] shaded leaf area                        
          laisun_z        =>    canopystate_vars%laisun_z_patch   , & ! Output: [real(r8) (:,:) ] sunlit leaf area for canopy layer     
          laisha_z        =>    canopystate_vars%laisha_z_patch   , & ! Output: [real(r8) (:,:) ] shaded leaf area for canopy layer     
          fsun            =>    canopystate_vars%fsun_patch       , & ! Output: [real(r8) (:)   ] sunlit fraction of canopy               

          fsa             =>    solarabs_vars%fsa_patch           , & ! Output: [real(r8) (:)   ] solar radiation absorbed (total) (W/m**2)
          fsr             =>    solarabs_vars%fsr_patch           , & ! Output: [real(r8) (:)   ] solar radiation reflected (W/m**2)      
          sabv            =>    solarabs_vars%sabv_patch          , & ! Output: [real(r8) (:)   ] solar radiation absorbed by vegetation (W/m**2)
          sabg            =>    solarabs_vars%sabg_patch          , & ! Output: [real(r8) (:)   ] solar radiation absorbed by ground (W/m**2)
          sabg_pen        =>    solarabs_vars%sabg_pen_patch      , & ! Output: [real(r8) (:)   ] solar (rural) radiation penetrating top soisno layer (W/m**2)
          sabg_soil       =>    solarabs_vars%sabg_soil_patch     , & ! Output: [real(r8) (:)   ] solar radiation absorbed by soil (W/m**2)
          sabg_snow       =>    solarabs_vars%sabg_snow_patch     , & ! Output: [real(r8) (:)   ] solar radiation absorbed by snow (W/m**2)
          sabg_lyr        =>    solarabs_vars%sabg_lyr_patch      , & ! Output: [real(r8) (:,:) ] absorbed radiative flux (pft,lyr) [W/m2]
          parsun_z        =>    solarabs_vars%parsun_z_patch      , & ! Output: [real(r8) (:,:) ] absorbed PAR for sunlit leaves in canopy layer
          parsha_z        =>    solarabs_vars%parsha_z_patch      , & ! Output: [real(r8) (:,:) ] absorbed PAR for shaded leaves in canopy layer
          fsr_nir_d       =>    solarabs_vars%fsr_nir_d_patch     , & ! Output: [real(r8) (:)   ] reflected direct beam nir solar radiation (W/m**2)
          fsr_nir_i       =>    solarabs_vars%fsr_nir_i_patch     , & ! Output: [real(r8) (:)   ] reflected diffuse nir solar radiation (W/m**2)
          fsr_nir_d_ln    =>    solarabs_vars%fsr_nir_d_ln_patch  , & ! Output: [real(r8) (:)   ] reflected direct beam nir solar rad at local noon (W/m**2)
          fsds_nir_d      =>    solarabs_vars%fsds_nir_d_patch    , & ! Output: [real(r8) (:)   ] incident direct beam nir solar radiation (W/m**2)
          fsds_nir_d_ln   =>    solarabs_vars%fsds_nir_d_ln_patch , & ! Output: [real(r8) (:)   ] incident direct beam nir solar rad at local noon (W/m**2)
          fsds_nir_i      =>    solarabs_vars%fsds_nir_i_patch    , & ! Output: [real(r8) (:)   ] incident diffuse nir solar radiation (W/m**2)
          fsa_r           =>    solarabs_vars%fsa_r_patch         , & ! Output: [real(r8) (:)   ] rural solar radiation absorbed (total) (W/m**2)
          sub_surf_abs_SW =>    solarabs_vars%sub_surf_abs_SW_col , & ! Output: [real(r8) (:)   ] percent of solar radiation absorbed below first snow layer (W/M**2)

          parveg_ln       =>    surfrad_vars%parveg_ln_patch      , & ! Output: [real(r8) (:)   ] absorbed par by vegetation at local noon (W/m**2)
          fsr_vis_d       =>    surfrad_vars%fsr_vis_d_patch      , & ! Output: [real(r8) (:)   ] reflected direct beam vis solar radiation (W/m**2)
          fsr_vis_i       =>    surfrad_vars%fsr_vis_i_patch      , & ! Output: [real(r8) (:)   ] reflected diffuse vis solar radiation (W/m**2)
          fsds_vis_i_ln   =>    surfrad_vars%fsds_vis_i_ln_patch  , & ! Output: [real(r8) (:)   ] incident diffuse beam vis solar rad at local noon (W/m**2)
          fsr_vis_d_ln    =>    surfrad_vars%fsr_vis_d_ln_patch   , & ! Output: [real(r8) (:)   ] reflected direct beam vis solar rad at local noon (W/m**2)
          fsds_vis_d      =>    surfrad_vars%fsds_vis_d_patch     , & ! Output: [real(r8) (:)   ] incident direct beam vis solar radiation (W/m**2)
          fsds_vis_i      =>    surfrad_vars%fsds_vis_i_patch     , & ! Output: [real(r8) (:)   ] incident diffuse vis solar radiation (W/m**2)
          fsds_vis_d_ln   =>    surfrad_vars%fsds_vis_d_ln_patch  , & ! Output: [real(r8) (:)   ] incident direct beam vis solar rad at local noon (W/m**2)
          sfc_frc_aer     =>    surfrad_vars%sfc_frc_aer_patch    , & ! Output: [real(r8) (:)   ] surface forcing of snow with all aerosols (pft) [W/m2]
          sfc_frc_aer_sno =>    surfrad_vars%sfc_frc_aer_sno_patch, & ! Output: [real(r8) (:)   ] surface forcing of snow with all aerosols, averaged only when snow is present (pft) [W/m2]
          sfc_frc_bc      =>    surfrad_vars%sfc_frc_bc_patch     , & ! Output: [real(r8) (:)   ] surface forcing of snow with BC (pft) [W/m2]
          sfc_frc_bc_sno  =>    surfrad_vars%sfc_frc_bc_sno_patch , & ! Output: [real(r8) (:)   ] surface forcing of snow with BC, averaged only when snow is present (pft) [W/m2]
          sfc_frc_oc      =>    surfrad_vars%sfc_frc_oc_patch     , & ! Output: [real(r8) (:)   ] surface forcing of snow with OC (pft) [W/m2]
          sfc_frc_oc_sno  =>    surfrad_vars%sfc_frc_oc_sno_patch , & ! Output: [real(r8) (:)   ] surface forcing of snow with OC, averaged only when snow is present (pft) [W/m2]
          sfc_frc_dst     =>    surfrad_vars%sfc_frc_dst_patch    , & ! Output: [real(r8) (:)   ] surface forcing of snow with dust (pft) [W/m2]
          sfc_frc_dst_sno =>    surfrad_vars%sfc_frc_dst_sno_patch, & ! Output: [real(r8) (:)   ] surface forcing of snow with dust, averaged only when snow is present (pft) [W/m2]
          fsr_sno_vd      =>    surfrad_vars%fsr_sno_vd_patch     , & ! Output: [real(r8) (:)   ] reflected visible, direct radiation from snow (for history files) (pft) [W/m2]
          fsr_sno_nd      =>    surfrad_vars%fsr_sno_nd_patch     , & ! Output: [real(r8) (:)   ] reflected near-IR, direct radiation from snow (for history files) (pft) [W/m2]
          fsr_sno_vi      =>    surfrad_vars%fsr_sno_vi_patch     , & ! Output: [real(r8) (:)   ] reflected visible, diffuse radiation from snow (for history files) (pft) [W/m2]
          fsr_sno_ni      =>    surfrad_vars%fsr_sno_ni_patch     , & ! Output: [real(r8) (:)   ] reflected near-IR, diffuse radiation from snow (for history files) (pft) [W/m2]
          fsds_sno_vd     =>    surfrad_vars%fsds_sno_vd_patch    , & ! Output: [real(r8) (:)   ] incident visible, direct radiation on snow (for history files) (pft) [W/m2]
          fsds_sno_nd     =>    surfrad_vars%fsds_sno_nd_patch    , & ! Output: [real(r8) (:)   ] incident near-IR, direct radiation on snow (for history files) (pft) [W/m2]
          fsds_sno_vi     =>    surfrad_vars%fsds_sno_vi_patch    , & ! Output: [real(r8) (:)   ] incident visible, diffuse radiation on snow (for history files) (pft) [W/m2]
          fsds_sno_ni     =>    surfrad_vars%fsds_sno_ni_patch      & ! Output: [real(r8) (:)   ] incident near-IR, diffuse radiation on snow (for history files) (pft) [W/m2]
          )

       ! Determine seconds off current time step
     
       dtime = get_step_size()
       call get_curr_date (year, month, day, secs)

       ! Initialize fluxes

       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          l = pft%landunit(p)
          g = pft%gridcell(p)

          sabg_soil(p)  = 0._r8
          sabg_snow(p)  = 0._r8
          sabg(p)       = 0._r8
          sabv(p)       = 0._r8
          fsa(p)        = 0._r8
          if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
             fsa_r(p) = 0._r8
          end if
          sabg_lyr(p,:) = 0._r8
          sabg_pur(p)   = 0._r8
          sabg_bc(p)    = 0._r8
          sabg_oc(p)    = 0._r8
          sabg_dst(p)   = 0._r8


          if( use_ed )then

             if ( EDpft%ED_patch(p) == 1 )then !#1
                currentSite => gridCellEdState(g)%spnt      
                currentPatch => gridCellEdState(g)%spnt%oldest_patch    
                do while(p /= currentPatch%clm_pno)
                   currentPatch => currentPatch%younger
                enddo
                currentPatch%ed_parsun_z(:,:,:) = 0._r8
                currentPatch%ed_parsha_z(:,:,:) = 0._r8
                currentPatch%ed_laisun_z(:,:,:) = 0._r8     
                currentPatch%ed_laisha_z(:,:,:) = 0._r8
                fsun(p) = 0._r8
             endif

          else ! not use_ed

             do iv = 1, nrad(p)
                parsun_z(p,iv) = 0._r8
                parsha_z(p,iv) = 0._r8
                laisun_z(p,iv) = 0._r8
                laisha_z(p,iv) = 0._r8
             end do

          end if ! end of if-use_ed

       end do

       ! Loop over patches to calculate laisun_z and laisha_z for each layer.
       ! Derive canopy laisun, laisha, and fsun from layer sums.
       ! If sun/shade big leaf code, nrad=1 and fsun_z(p,1) and tlai_z(p,1) from
       ! SurfaceAlbedo is canopy integrated so that layer value equals canopy value.

       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          g = pft%gridcell(p)

          if( use_ed )then

             ! currentPatch%f_sun is calculated in the surface_albedo routine...
             if(EDpft%ED_patch(p).eq.1)then
                fsun(p) = 0._r8
                sunlai = 0._r8
                shalai = 0._r8
                currentSite => gridCellEdState(g)%spnt      
                currentPatch => gridCellEdState(g)%spnt%oldest_patch    
                do while(p /= currentPatch%clm_pno)
                   currentPatch => currentPatch%younger
                enddo
                do CL = 1, currentPatch%NCL_p
                   do FT = 1,numpft_ed
                      do iv = 1, currentPatch%nrad(CL,ft) !NORMAL CASE. 
                         ! FIX(SPM,040114) ** Should this be elai or tlai? Surely we only do radiation for elai? 
                         currentPatch%ed_laisun_z(CL,ft,iv) = currentPatch%elai_profile(CL,ft,iv) * &
                              currentPatch%f_sun(CL,ft,iv)
                         currentPatch%ed_laisha_z(CL,ft,iv) = currentPatch%elai_profile(CL,ft,iv) * &
                              (1._r8 - currentPatch%f_sun(CL,ft,iv))
                      end do
                      sunlai = sunlai + sum(currentPatch%ed_laisun_z(CL,ft,1: currentPatch%nrad(CL,ft)))
                      shalai = shalai + sum(currentPatch%ed_laisha_z(CL,ft,1: currentPatch%nrad(CL,ft)))
                      !needed for the VOC emissions, etc. 
                   end do
                end do
                if(sunlai+shalai > 0._r8)then
                   fsun(p) = sunlai / (sunlai+shalai) 
                else
                   fsun(p) = 0._r8
                endif
                if(fsun(p) > 1._r8)then
                   write(iulog,*) 'too much leaf area in profile', fsun(p),currentPatch%lai,sunlai,shalai
                endif

             else ! not ed patch

               fsun(p) = 0.0_r8

             end if !ED_patch   

          else ! use_ed false.  revert to normal multi-layer canopy.

             laisun(p) = 0._r8
             laisha(p) = 0._r8
             do iv = 1, nrad(p)
                laisun_z(p,iv) = tlai_z(p,iv) * fsun_z(p,iv)
                laisha_z(p,iv) = tlai_z(p,iv) * (1._r8 - fsun_z(p,iv))
                laisun(p) = laisun(p) + laisun_z(p,iv) 
                laisha(p) = laisha(p) + laisha_z(p,iv) 
             end do
             if (elai(p) > 0._r8) then
                fsun(p) = laisun(p) / elai(p)
             else
                fsun(p) = 0._r8
             end if

          end if ! end of if-use_ed  

       end do ! end of fp = 1,num_nourbanp loop

       do ib = 1, numrad
          do fp = 1,num_urbanp
             p = filter_urbanp(fp)
             if (ib == 1) then
                fsun(p) = 0._r8
             end if
          end do
       end do

       ! Loop over nband wavebands
       do ib = 1, nband
          do fp = 1,num_nourbanp
             p = filter_nourbanp(fp)
             c = pft%column(p)
             l = pft%landunit(p)
             g = pft%gridcell(p)

             ! Absorbed by canopy

             cad(p,ib) = forc_solad(g,ib)*fabd(p,ib)
             cai(p,ib) = forc_solai(g,ib)*fabi(p,ib)
             sabv(p) = sabv(p) + cad(p,ib) + cai(p,ib)
             fsa(p)  = fsa(p)  + cad(p,ib) + cai(p,ib)
             if (ib == 1) then
                parveg(p) = cad(p,ib) + cai(p,ib)
             end if
             if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + cad(p,ib) + cai(p,ib)
             end if

             ! Absorbed PAR profile through canopy
             ! If sun/shade big leaf code, nrad=1 and fluxes from SurfaceAlbedo
             ! are canopy integrated so that layer values equal big leaf values.

             if (ib == 1) then

                if ( use_ed ) then   
                   if (EDpft%ED_patch(p).eq.1 )then

                      currentSite => gridCellEdState(g)%spnt      
                      currentPatch => gridCellEdState(g)%spnt%oldest_patch    
                      do while(p /= currentPatch%clm_pno)
                         currentPatch => currentPatch%younger
                      enddo
                      do CL = 1, currentPatch%NCL_p
                         do FT = 1,numpft_ed
                            do iv = 1, currentPatch%nrad(CL,ft)
                               currentPatch%ed_parsun_z(CL,ft,iv) = forc_solad(g,ib)*currentPatch%fabd_sun_z(CL,ft,iv) + &
                                    forc_solai(g,ib)*currentPatch%fabi_sun_z(CL,ft,iv) 
                               currentPatch%ed_parsha_z(CL,ft,iv) = forc_solad(g,ib)*currentPatch%fabd_sha_z(CL,ft,iv) + &
                                    forc_solai(g,ib)*currentPatch%fabi_sha_z(CL,ft,iv)          
                            end do !iv
                         end do !FT
                      end do !CL
                   end if ! ED_patch check

                else ! not use_ed

                   do iv = 1, nrad(p)
                      parsun_z(p,iv) = forc_solad(g,ib)*fabd_sun_z(p,iv) + forc_solai(g,ib)*fabi_sun_z(p,iv)
                      parsha_z(p,iv) = forc_solad(g,ib)*fabd_sha_z(p,iv) + forc_solai(g,ib)*fabi_sha_z(p,iv)
                   end do

                end if  ! end of if-use_ed 

             end if   ! end of if ib is 1

             ! Transmitted = solar fluxes incident on ground

             trd(p,ib) = forc_solad(g,ib)*ftdd(p,ib)
             tri(p,ib) = forc_solad(g,ib)*ftid(p,ib) + forc_solai(g,ib)*ftii(p,ib)
             ! Solar radiation absorbed by ground surface
             ! calculate absorbed solar by soil/snow separately
             absrad  = trd(p,ib)*(1._r8-albsod(c,ib)) + tri(p,ib)*(1._r8-albsoi(c,ib))
             sabg_soil(p) = sabg_soil(p) + absrad
             absrad  = trd(p,ib)*(1._r8-albsnd_hst(c,ib)) + tri(p,ib)*(1._r8-albsni_hst(c,ib))
             sabg_snow(p) = sabg_snow(p) + absrad
             absrad  = trd(p,ib)*(1._r8-albgrd(c,ib)) + tri(p,ib)*(1._r8-albgri(c,ib))
             sabg(p) = sabg(p) + absrad
             fsa(p)  = fsa(p)  + absrad
             if (lun%itype(l)==istsoil .or. lun%itype(l)==istcrop) then
                fsa_r(p)  = fsa_r(p)  + absrad
             end if
             if (snl(c) == 0) then
                sabg_snow(p) = sabg(p)
                sabg_soil(p) = sabg(p)
             endif
             ! if no subgrid fluxes, make sure to set both components equal to weighted average
             if (subgridflag == 0) then 
                sabg_snow(p) = sabg(p)
                sabg_soil(p) = sabg(p)
             endif

             if (use_snicar_frc) then
                ! Solar radiation absorbed by ground surface without BC
                absrad_bc = trd(p,ib)*(1._r8-albgrd_bc(c,ib)) + tri(p,ib)*(1._r8-albgri_bc(c,ib))
                sabg_bc(p) = sabg_bc(p) + absrad_bc

                ! Solar radiation absorbed by ground surface without OC
                absrad_oc = trd(p,ib)*(1._r8-albgrd_oc(c,ib)) + tri(p,ib)*(1._r8-albgri_oc(c,ib))
                sabg_oc(p) = sabg_oc(p) + absrad_oc

                ! Solar radiation absorbed by ground surface without dust
                absrad_dst = trd(p,ib)*(1._r8-albgrd_dst(c,ib)) + tri(p,ib)*(1._r8-albgri_dst(c,ib))
                sabg_dst(p) = sabg_dst(p) + absrad_dst

                ! Solar radiation absorbed by ground surface without any aerosols
                absrad_pur = trd(p,ib)*(1._r8-albgrd_pur(c,ib)) + tri(p,ib)*(1._r8-albgri_pur(c,ib))
                sabg_pur(p) = sabg_pur(p) + absrad_pur
             end if

          end do ! end of pft loop
       end do ! end nbands loop   

       !   compute absorbed flux in each snow layer and top soil layer,
       !   based on flux factors computed in the radiative transfer portion of SNICAR.

       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          c = pft%column(p)
          l = pft%landunit(p)
          sabg_snl_sum = 0._r8

          sub_surf_abs_SW(c) = 0._r8

          ! CASE1: No snow layers: all energy is absorbed in top soil layer
          if (snl(c) == 0) then
             sabg_lyr(p,:) = 0._r8
             sabg_lyr(p,1) = sabg(p)
             sabg_snl_sum  = sabg_lyr(p,1)

             ! CASE 2: Snow layers present: absorbed radiation is scaled according to 
             ! flux factors computed by SNICAR
          else
             do i = -nlevsno+1,1,1
                sabg_lyr(p,i) = flx_absdv(c,i)*trd(p,1) + flx_absdn(c,i)*trd(p,2) + &
                     flx_absiv(c,i)*tri(p,1) + flx_absin(c,i)*tri(p,2)
                ! summed radiation in active snow layers:
                if (i >= snl(c)+1) then
                   sabg_snl_sum = sabg_snl_sum + sabg_lyr(p,i)
                endif
                if (i > snl(c)+1) then ! if snow layer is below surface snow layer
                   !accumulate subsurface flux as a diagnostic for history file
                   sub_surf_abs_SW(c) = sub_surf_abs_SW(c) + sabg_lyr(p,i)
                endif
             enddo

             ! Divide absorbed by total, to get % absorbed in subsurface
             if (sabg_snl_sum /= 0._r8) then
                sub_surf_abs_SW(c) = sub_surf_abs_SW(c)/sabg_snl_sum
             else
                sub_surf_abs_SW(c) = 0._r8
             endif

             ! Error handling: The situation below can occur when solar radiation is 
             ! NOT computed every timestep.
             ! When the number of snow layers has changed in between computations of the 
             ! absorbed solar energy in each layer, we must redistribute the absorbed energy
             ! to avoid physically unrealistic conditions. The assumptions made below are 
             ! somewhat arbitrary, but this situation does not arise very frequently. 
             ! This error handling is implemented to accomodate any value of the
             ! radiation frequency.
             ! change condition to match sabg_snow isntead of sabg
             if (abs(sabg_snl_sum-sabg_snow(p)) > 0.00001_r8) then
                if (snl(c) == 0) then
                   sabg_lyr(p,-4:0) = 0._r8
                   sabg_lyr(p,1) = sabg(p)
                elseif (snl(c) == -1) then
                   sabg_lyr(p,-4:-1) = 0._r8
                   sabg_lyr(p,0) = sabg_snow(p)*0.6_r8
                   sabg_lyr(p,1) = sabg_snow(p)*0.4_r8
                else
                   sabg_lyr(p,:) = 0._r8
                   sabg_lyr(p,snl(c)+1) = sabg_snow(p)*0.75_r8
                   sabg_lyr(p,snl(c)+2) = sabg_snow(p)*0.25_r8
                endif
             endif

             ! If shallow snow depth, all solar radiation absorbed in top or top two snow layers
             ! to prevent unrealistic timestep soil warming 
             if (subgridflag == 0) then 
                if (snow_depth(c) < 0.10_r8) then
                   if (snl(c) == 0) then
                      sabg_lyr(p,-4:0) = 0._r8
                      sabg_lyr(p,1) = sabg(p)
                   elseif (snl(c) == -1) then
                      sabg_lyr(p,-4:-1) = 0._r8
                      sabg_lyr(p,0) = sabg(p)
                      sabg_lyr(p,1) = 0._r8
                   else
                      sabg_lyr(p,:) = 0._r8
                      sabg_lyr(p,snl(c)+1) = sabg(p)*0.75_r8
                      sabg_lyr(p,snl(c)+2) = sabg(p)*0.25_r8
                   endif
                endif
             endif
          endif

          ! This situation should not happen:
          if (abs(sum(sabg_lyr(p,:))-sabg_snow(p)) > 0.00001_r8) then
             write(iulog,*)"SNICAR ERROR: Absorbed ground radiation not equal to summed snow layer radiation"
             write(iulog,*)"Diff        = ",sum(sabg_lyr(p,:))-sabg_snow(p)
             write(iulog,*)"sabg_snow(p)= ",sabg_snow(p)
             write(iulog,*)"sabg_sum(p) = ",sum(sabg_lyr(p,:))
             write(iulog,*)"snl(c)      = ",snl(c)
             write(iulog,*)"flx_absdv1  = ",trd(p,1)*(1.-albgrd(c,1))
             write(iulog,*)"flx_absdv2  = ",sum(flx_absdv(c,:))*trd(p,1)
             write(iulog,*)"flx_absiv1  = ",tri(p,1)*(1.-albgri(c,1))
             write(iulog,*)"flx_absiv2  = ",sum(flx_absiv(c,:))*tri(p,1)
             write(iulog,*)"flx_absdn1  = ",trd(p,2)*(1.-albgrd(c,2))
             write(iulog,*)"flx_absdn2  = ",sum(flx_absdn(c,:))*trd(p,2)
             write(iulog,*)"flx_absin1  = ",tri(p,2)*(1.-albgri(c,2))
             write(iulog,*)"flx_absin2  = ",sum(flx_absin(c,:))*tri(p,2)
             write(iulog,*)"albgrd_nir  = ",albgrd(c,2)
             write(iulog,*)"coszen      = ",coszen(c)
             call endrun(decomp_index=c, clmlevel=namec, msg=errmsg(__FILE__, __LINE__))
          endif

          ! Diagnostic: shortwave penetrating ground (e.g. top layer)
          if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
             sabg_pen(p) = sabg(p) - sabg_lyr(p, snl(c)+1)
          end if

          if (use_snicar_frc) then

             ! BC aerosol forcing (pft-level):
             sfc_frc_bc(p) = sabg(p) - sabg_bc(p)

             ! OC aerosol forcing (pft-level):
             if (DO_SNO_OC) then
                sfc_frc_oc(p) = sabg(p) - sabg_oc(p)
             else
                sfc_frc_oc(p) = 0._r8
             endif

             ! dust aerosol forcing (pft-level):
             sfc_frc_dst(p) = sabg(p) - sabg_dst(p)

             ! all-aerosol forcing (pft-level):
             sfc_frc_aer(p) = sabg(p) - sabg_pur(p)        

             ! forcings averaged only over snow:
             if (frac_sno(c) > 0._r8) then
                sfc_frc_bc_sno(p)  = sfc_frc_bc(p)/frac_sno(c)
                sfc_frc_oc_sno(p)  = sfc_frc_oc(p)/frac_sno(c)
                sfc_frc_dst_sno(p) = sfc_frc_dst(p)/frac_sno(c)
                sfc_frc_aer_sno(p) = sfc_frc_aer(p)/frac_sno(c)
             else
                sfc_frc_bc_sno(p)  = spval
                sfc_frc_oc_sno(p)  = spval
                sfc_frc_dst_sno(p) = spval
                sfc_frc_aer_sno(p) = spval
             endif
          end if
       enddo

       ! Radiation diagnostics
       
       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          g = pft%gridcell(p)

          ! NDVI and reflected solar radiation

          rvis = albd(p,1)*forc_solad(g,1) + albi(p,1)*forc_solai(g,1)
          rnir = albd(p,2)*forc_solad(g,2) + albi(p,2)*forc_solai(g,2)
          fsr(p) = rvis + rnir

          fsds_vis_d(p) = forc_solad(g,1)
          fsds_nir_d(p) = forc_solad(g,2)
          fsds_vis_i(p) = forc_solai(g,1)
          fsds_nir_i(p) = forc_solai(g,2)
          fsr_vis_d(p)  = albd(p,1)*forc_solad(g,1)
          fsr_nir_d(p)  = albd(p,2)*forc_solad(g,2)
          fsr_vis_i(p)  = albi(p,1)*forc_solai(g,1)
          fsr_nir_i(p)  = albi(p,2)*forc_solai(g,2)

          local_secp1 = secs + nint((grc%londeg(g)/degpsec)/dtime)*dtime
          local_secp1 = mod(local_secp1,isecspday)
          if (local_secp1 == isecspday/2) then
             fsds_vis_d_ln(p) = forc_solad(g,1)
             fsds_nir_d_ln(p) = forc_solad(g,2)
             fsr_vis_d_ln(p) = albd(p,1)*forc_solad(g,1)
             fsr_nir_d_ln(p) = albd(p,2)*forc_solad(g,2)
             fsds_vis_i_ln(p) = forc_solai(g,1)
             parveg_ln(p)     = parveg(p)
          else
             fsds_vis_d_ln(p) = spval
             fsds_nir_d_ln(p) = spval
             fsr_vis_d_ln(p) = spval
             fsr_nir_d_ln(p) = spval
             fsds_vis_i_ln(p) = spval
             parveg_ln(p)     = spval
          end if

          ! diagnostic variables (downwelling and absorbed radiation partitioning) for history files
          ! (OPTIONAL)
          c = pft%column(p)
          if (snl(c) < 0) then
             fsds_sno_vd(p) = forc_solad(g,1)
             fsds_sno_nd(p) = forc_solad(g,2)
             fsds_sno_vi(p) = forc_solai(g,1)
             fsds_sno_ni(p) = forc_solai(g,2)

             fsr_sno_vd(p) = fsds_vis_d(p)*albsnd_hst(c,1)
             fsr_sno_nd(p) = fsds_nir_d(p)*albsnd_hst(c,2)
             fsr_sno_vi(p) = fsds_vis_i(p)*albsni_hst(c,1)
             fsr_sno_ni(p) = fsds_nir_i(p)*albsni_hst(c,2)
          else
             fsds_sno_vd(p) = spval
             fsds_sno_nd(p) = spval
             fsds_sno_vi(p) = spval
             fsds_sno_ni(p) = spval

             fsr_sno_vd(p) = spval
             fsr_sno_nd(p) = spval
             fsr_sno_vi(p) = spval
             fsr_sno_ni(p) = spval
          endif
       end do

       do fp = 1,num_urbanp
          p = filter_urbanp(fp)
          g = pft%gridcell(p)

          local_secp1 = secs + nint((grc%londeg(g)/degpsec)/dtime)*dtime
          local_secp1 = mod(local_secp1,isecspday)

        if(elai(p)==0.0_r8.and.fabd(p,1)>0._r8)then
           ! FIX(SPM, 051314) - is this necessary ?  puts lots of info in
           ! lnd.log
           write(iulog,*) 'absorption without LAI',elai(p),tlai(p),fabd(p,1),p
        endif
          ! Solar incident 

          fsds_vis_d(p) = forc_solad(g,1)
          fsds_nir_d(p) = forc_solad(g,2)    
          fsds_vis_i(p) = forc_solai(g,1)
          fsds_nir_i(p) = forc_solai(g,2)

          ! Determine local noon incident solar
          if (local_secp1 == noonsec) then
             fsds_vis_d_ln(p) = forc_solad(g,1)
             fsds_nir_d_ln(p) = forc_solad(g,2)
             fsds_vis_i_ln(p) = forc_solai(g,1)
             parveg_ln(p)     = 0._r8
          else
             fsds_vis_d_ln(p) = spval 
             fsds_nir_d_ln(p) = spval 
             fsds_vis_i_ln(p) = spval
             parveg_ln(p)     = spval
          endif

          ! Solar reflected 
          ! per unit ground area (roof, road) and per unit wall area (sunwall, shadewall)

          fsr_vis_d(p) = albd(p,1) * forc_solad(g,1)
          fsr_nir_d(p) = albd(p,2) * forc_solad(g,2)
          fsr_vis_i(p) = albi(p,1) * forc_solai(g,1)
          fsr_nir_i(p) = albi(p,2) * forc_solai(g,2)

          ! Determine local noon reflected solar
          if (local_secp1 == noonsec) then
             fsr_vis_d_ln(p) = fsr_vis_d(p)
             fsr_nir_d_ln(p) = fsr_nir_d(p)
          else
             fsr_vis_d_ln(p) = spval 
             fsr_nir_d_ln(p) = spval 
          endif
          fsr(p) = fsr_vis_d(p) + fsr_nir_d(p) + fsr_vis_i(p) + fsr_nir_i(p)  
       end do

       do fp = 1,num_nourbanp
          p = filter_nourbanp(fp)
          g = pft%gridcell(p)
          if (use_ed) then
             errsol(p) = (fsa(p) + fsr(p)  - (forc_solad(g,1) + forc_solad(g,2) + forc_solai(g,1) + forc_solai(g,2)))
             if(abs(errsol(p)) > 0.1_r8)then
                g = pft%gridcell(p)
                write(iulog,*) 'sol error in surf rad',p,g, errsol(p),EDpft%ed_patch(p)
             endif
          end if
       end do

     end associate

   end subroutine SurfaceRadiation
  
end module SurfaceRadiationMod
